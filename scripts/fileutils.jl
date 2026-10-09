using CSV, DataFrames
using RCall

# ─────────────────────────────────────────────────────────────────────────────
# Shared helpers for finding & loading per-taxon files
#
# File naming conventions used throughout:
#   raw occurrences      <taxon>-<date>.csv                 (ex. Genus_species-2026_10_06.csv)
#   cleaned occurrences  <taxon>-<date>_cleaned.csv
#   missing coordinates  <taxon>-<date>_no_coords.csv
#   final occurrences    <taxon>.csv                        (cleaned + georeferenced + filtered)
# where <taxon> is the taxon name with spaces replaced by underscores.
# ─────────────────────────────────────────────────────────────────────────────

"""
    synonym_names(syns) -> Vector{String}

Get the (unique) names out of a taxon's synonym list, which can either be a
vector of names or a vector of `(name, authority)` tuples.
"""
synonym_names(syns) = unique(String[syn isa AbstractString ? syn : first(syn) for syn in syns])

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
    resolve_clean_file(filename, dir) -> String

Return the occurrence CSV for a taxon in `dir`: `<taxon>.csv`, or the most
recent `<taxon>-<date>_cleaned.csv` if the files are dated.
Returns an empty string if nothing is found.
"""
function resolve_clean_file(filename::String, dir::String)
    stem = taxon_stem(filename)

    # match on the whole taxon stem, so a species doesn't pick up its infraspecific taxa's files
    candidates = filter(readdir(dir)) do f
        endswith(f, ".csv") && !endswith(f, "_georef_merged.csv") && !endswith(f, "_no_coords.csv") && taxon_stem(f) == stem
    end
    return isempty(candidates) ? "" : joinpath(dir, last(sort(candidates)))
end

"""
    latest_no_coords_file(filename, no_coords_dir) -> String

Return the most recent `<taxon>-<date>_no_coords.csv` for a taxon in
`no_coords_dir`, or an empty string if nothing is found.
"""
function latest_no_coords_file(filename::String, no_coords_dir::String)
    stem = taxon_stem(filename)
    candidates = filter(f -> endswith(f, "_no_coords.csv") && taxon_stem(f) == stem, readdir(no_coords_dir))
    return isempty(candidates) ? "" : joinpath(no_coords_dir, last(sort(candidates)))
end

"""
    taxon_files(dir) -> Vector of taxon stem => path

One occurrence CSV per taxon in `dir` (via `resolve_clean_file`), sorted by taxon.
"""
function taxon_files(dir::String)
    stems = sort(unique(taxon_stem.(filter(f -> endswith(f, ".csv"), readdir(dir)))))
    return filter(p -> !isempty(last(p)), [s => resolve_clean_file(s, dir) for s in stems])
end

"""
    read_occs(path) -> DataFrame

Read an occurrence CSV, with `""` & `"NA"` as missing and the `ID` column as a
String (GBIF IDs are numbers & iDigBio IDs are UUIDs, so a file can have either or both).
"""
read_occs(path::String) = CSV.read(path, DataFrame; missingstring=["", "NA"], types=Dict(:ID => String))

"""
    normalize_id(id) -> String (or missing)

Get an occurrence ID as a plain string. GBIF IDs that have been through a
spreadsheet can come back as floats or as text in E-notation
(ex. "1.26039201E9"), which wouldn't match the IDs in the occurrence files.
"""
function normalize_id(id)
    ismissing(id) && return missing
    id isa Integer && return string(id)
    id isa AbstractFloat && return string(round(Int, id))
    s = strip(string(id))
    if occursin(r"^\d+\.\d+$", s) || occursin(r"^\d+(\.\d+)?[eE]\+?\d+$", s)
        return string(round(Int, parse(Float64, s)))
    end
    return String(s)
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
    read_georef_file(path) -> DataFrame

Read a georeferenced points file (.csv or .xlsx) with every row, viable or not.
`ID` is a plain string (see `normalize_id`), and the viable & notes columns are
named `VIABLE` & `NOTES` whatever their capitalisation in the file; `VIABLE` is
"TRUE", "FALSE" or "" (not filled in).
"""
function read_georef_file(path::String)
    df = if endswith(path, ".xlsx")
        georef_path = path
        @rput georef_path
        R"georef_r <- readxl::read_excel(georef_path)"
        @rget georef_r
        georef_r
    else
        DataFrame(CSV.File(path))
    end

    for col in ("VIABLE", "NOTES")
        i = findfirst(c -> uppercase(strip(c)) == col, names(df))
        isnothing(i) || names(df)[i] == col || rename!(df, names(df)[i] => col)
    end
    if "VIABLE" in names(df)
        df[!, :VIABLE] = String[ismissing(v) ? "" : uppercase(strip(string(v))) for v in df.VIABLE]
    end
    "ID" in names(df) && (df[!, :ID] = normalize_id.(df.ID))
    return df
end

"""
    load_georef_df(taxa_filename, georef_dir) -> DataFrame or nothing

Look for a georeferenced points file (.csv or .xlsx) for a taxon, filter to
viable rows, and return a DataFrame with Float64 lat/lon columns & `ID` as a
plain string (see `read_georef_file`) — or `nothing` if no usable file is found.
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
            df = read_georef_file(georef_path)

            if !("VIABLE" in names(df))
                println("  WARNING: No 'Viable' column found in $georef_path — skipping georef points")
                return nothing
            end
            filter!(row -> row.VIABLE == "TRUE", df)

            df[!, :latitude] = map(v -> ismissing(v) ? missing : tryparse(Float64, string(v)), df.latitude)
            df[!, :longitude] = map(v -> ismissing(v) ? missing : tryparse(Float64, string(v)), df.longitude)
            # tryparse gives `nothing` for coordinates that aren't numbers (ex. "NA"), so drop those too
            filter!(row -> row.latitude isa Number && row.longitude isa Number, df)
            df[!, :latitude] = Float64.(df.latitude)
            df[!, :longitude] = Float64.(df.longitude)

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
    load_pt_occs_df(dir) -> DataFrame

Stack the occurrence CSVs in `dir` into a single DataFrame, one file per taxon
(via `taxon_files`).  Two extra columns are prepended:

    taxon       – bare taxon name (spaces → underscores, no date/suffix)
    source_file – basename of the file that was loaded

All other columns are passed through as-is.
"""
function load_pt_occs_df(dir::String)
    chunks = DataFrame[]
    for (t, fpath) in taxon_files(dir)
        df = read_occs(fpath)
        insertcols!(df, 1,
            :taxon => t,
            :source_file => basename(fpath)
        )
        push!(chunks, df)
    end

    return vcat(chunks...; cols=:union)
end
