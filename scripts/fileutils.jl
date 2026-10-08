using CSV, DataFrames
using RCall

# ─────────────────────────────────────────────────────────────────────────────
# Shared helpers for finding & loading per-taxon files
#
# File naming conventions used throughout:
#   raw occurrences      <taxon>-<date>.csv                 (ex. Genus_species-2026_10_06.csv)
#   cleaned occurrences  <taxon>-<date>_cleaned.csv
#   missing coordinates  <taxon>-<date>_no_coords.csv
#   georef merged        <taxon>_georef_merged.csv
# where <taxon> is the taxon name with spaces replaced by underscores.
# ─────────────────────────────────────────────────────────────────────────────

"""
    synonym_names(syns) -> Vector{String}

Get the names out of a taxon's synonym list, which can either be a vector of
names or a vector of `(name, authority)` tuples.
"""
synonym_names(syns) = String[syn isa AbstractString ? syn : first(syn) for syn in syns]

"""
    taxon_stem(fname::String) -> String

Strip suffixes from a filename (like .csv, _georef_merged, or date strings)
to get the bare taxon stem name.
"""
function taxon_stem(fname::String)::String
    base = replace(fname, r"\.(csv|xlsx)$" => "")
    base = replace(base, r"_georef_merged$" => "")
    base = replace(base, r"-\d{4}[-_]\d{2}[-_]\d{2}(_cleaned|_no_coords)?$" => "")
    return base
end

"""
    latest_raw_file(filename, raw_dir) -> (path, date_suffix)

Return the most recently pulled raw CSV for a taxon in `raw_dir` (files named
`<filename>-<date>.csv`), plus its date as a `yyyy_mm_dd` string.
The date in the filename can be separated by `_` or `-`; ties on the date are
broken by the file's modification time. Returns `("", "")` if nothing is found.
"""
function latest_raw_file(filename::String, raw_dir::String)
    pat = Regex("^\\Q$(filename)\\E-(\\d{4})[-_](\\d{2})[-_](\\d{2})\\.csv\$")
    candidates = filter(f -> occursin(pat, f), readdir(raw_dir))
    isempty(candidates) && return ("", "")

    filedate(f) = join(match(pat, f).captures, "_")
    latest = last(sort(candidates, by = f -> (filedate(f), mtime(joinpath(raw_dir, f)))))
    return (joinpath(raw_dir, latest), filedate(latest))
end

"""
    resolve_clean_file(filename, clean_dir) -> String

Return the best available cleaned CSV path for a taxon in `clean_dir`:
prefers `*_georef_merged.csv`, falls back to the most recent `*_cleaned.csv`.
Returns an empty string if nothing is found.
"""
function resolve_clean_file(filename::String, clean_dir::String)
    stem = taxon_stem(filename)

    # match on the whole taxon stem, so a species doesn't pick up its infraspecific taxa's files
    taxon_files = filter(f -> endswith(f, ".csv") && taxon_stem(f) == stem, readdir(clean_dir))

    georef_candidates = filter(f -> endswith(f, "_georef_merged.csv"), taxon_files)
    cleaned_candidates = filter(f -> occursin("cleaned", f), taxon_files)
    if !isempty(georef_candidates)
        return joinpath(clean_dir, first(georef_candidates))
    elseif !isempty(cleaned_candidates)
        return joinpath(clean_dir, last(sort(cleaned_candidates)))
    else
        return ""
    end
end

"""
    load_taxa_data(traits_path; species_ranges_path=nothing) -> (taxa_traits, taxa_to_nativerange_dict)

Load a taxa CSV (needs `scientificName` & `botanicalCountries` columns) and build
a scientificName → botanical country codes lookup.

If `species_ranges_path` is given (a CSV with the same two columns), those taxa
are appended to `taxa_traits` — ex. species-level taxa that were pulled for taxa
sampled at the infraspecific level.
"""
function load_taxa_data(traits_path::String; species_ranges_path::Union{String,Nothing}=nothing)
    taxa_traits = DataFrame(CSV.File(traits_path))
    if !isnothing(species_ranges_path)
        taxa_traits = vcat(taxa_traits, DataFrame(CSV.File(species_ranges_path)); cols=:union)
    end
    # botanicalCountries is written out as a string (ex. ["CLC", "CLS"]), so parse out the codes
    taxa_to_nativerange_dict = Dict(
        row.scientificName => String[m.match for m in eachmatch(r"[A-Z]{3}", coalesce(row.botanicalCountries, ""))]
        for row in eachrow(taxa_traits)
    )
    return taxa_traits, taxa_to_nativerange_dict
end

"""
    load_georef_df(taxa_filename, georef_dir) -> DataFrame or nothing

Look for a georeferenced points file (.csv or .xlsx) for a taxon, filter to
viable rows, and return a DataFrame with Float64 lat/lon columns — or `nothing`
if no usable file is found.
"""
function load_georef_df(taxa_filename::String, georef_dir::String)
    isdir(georef_dir) || return nothing
    stem = taxon_stem(taxa_filename)

    for ext in [".csv", ".xlsx"]
        candidates = filter(
            f -> endswith(f, ext) && taxon_stem(f) == stem,
            readdir(georef_dir)
        )
        isempty(candidates) && continue

        georef_path = joinpath(georef_dir, first(candidates))
        println("  Found georeferenced file: $georef_path")

        try
            df = if ext == ".csv"
                DataFrame(CSV.File(georef_path))
            else
                @rput georef_path
                R"georef_r <- readxl::read_excel(georef_path)"
                @rget georef_r
                georef_r
            end

            viable_col = findfirst(c -> uppercase(string(c)) == "VIABLE", names(df))
            if isnothing(viable_col)
                println("  WARNING: No 'Viable' column found in $georef_path — skipping georef points")
                return nothing
            end

            rename!(df, names(df)[viable_col] => :Viable)
            df[!, :Viable] = map(v -> uppercase(strip(string(v))), df.Viable)
            filter!(row -> row.Viable == "TRUE", df)

            df[!, :latitude] = map(v -> ismissing(v) ? missing : tryparse(Float64, string(v)), df.latitude)
            df[!, :longitude] = map(v -> ismissing(v) ? missing : tryparse(Float64, string(v)), df.longitude)
            filter!(row -> !ismissing(row.latitude) && !ismissing(row.longitude), df)

            println("  Georeferenced viable points: $(nrow(df))")
            return nrow(df) > 0 ? df : nothing

        catch e
            println("  WARNING: Could not read georef file $georef_path: $e")
            return nothing
        end
    end
    return nothing
end

"""
    load_pt_occs_df(clean_dir) -> DataFrame

Stack all cleaned occurrence CSVs into a single DataFrame, selecting the best
available file per taxon (`_georef_merged.csv` preferred over `_cleaned.csv`,
via `resolve_clean_file`).  Two extra columns are prepended:

    taxon       – bare taxon name (spaces → underscores, no date/suffix)
    source_file – basename of the file that was loaded

All other columns are passed through as-is.
"""
function load_pt_occs_df(clean_dir::String)
    all_files = filter(f -> endswith(f, ".csv"), readdir(clean_dir))
    taxa = sort(unique(taxon_stem.(all_files)))

    chunks = DataFrame[]
    for t in taxa
        fpath = resolve_clean_file(t, clean_dir)
        isempty(fpath) && continue

        df = CSV.read(fpath, DataFrame; missingstring=["", "NA"])
        insertcols!(df, 1,
            :taxon => t,
            :source_file => basename(fpath)
        )
        push!(chunks, df)
    end

    return vcat(chunks...; cols=:union)
end
