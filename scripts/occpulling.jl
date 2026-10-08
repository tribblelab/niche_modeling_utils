using CSV, DataFrames
using RCall
using Dates

# needs fileutils.jl (synonym_names, latest_raw_file)

R"""
library(gatoRs)
"""

"""
    pull_occurrences(synonymdict, raw_dir; kwargs...) -> (successful, no_data, failed)

Download GBIF / iDigBio occurrences (via `gatoRs::gators_download`) for every
taxon in `synonymdict`, writing one raw CSV per taxon to
`raw_dir/<taxon>-<date_suffix>.csv`.

Returns the taxa that were pulled (or already had a file), the taxa with no
occurrence records, and the taxa that errored.

# Arguments
- `synonymdict`: taxon name => its synonym list (a vector of names, or of
  `(name, authority)` tuples). All names in the list are searched.
- `raw_dir`: Directory to write the raw CSVs to. Created if absent.
- `date_suffix`: Date string for the filenames (default: today, as `yyyy_mm_dd`).
- `taxa`: Which taxa in `synonymdict` to pull (default: all, sorted).
- `overwrite`: If false (default), skip taxa whose file already exists.
- `sleep_after`: Seconds to wait between downloads (default: 3).

# Example
```julia
synonymdict = Dict("Calceolaria alba" => ["Calceolaria alba", "Fagelia alba"])
pulled = pull_occurrences(synonymdict, "data/occurrence_data/pt_occs_raw")
pulled.failed
```
"""
function pull_occurrences(synonymdict::AbstractDict, raw_dir::String;
    date_suffix::String=Dates.format(Dates.today(), "yyyy_mm_dd"),
    taxa=sort(collect(keys(synonymdict))),
    overwrite::Bool=false,
    sleep_after::Real=3)

    mkpath(raw_dir)

    successful = String[]
    no_data = String[]
    failed = String[]

    for (idx, taxon) in enumerate(taxa)
        synlist = synonym_names(synonymdict[taxon])
        filepath = joinpath(raw_dir, "$(replace(taxon, " " => "_"))-$(date_suffix).csv")
        println("Pulling $idx/$(length(taxa)): $taxon ($(length(synlist)) names)")

        if !overwrite && isfile(filepath)
            println("  Skipping - file already exists")
            push!(successful, taxon)
            continue
        end

        @rput synlist filepath
        try
            R"gators_download(synonyms.list = synlist, write.file = TRUE, filename = filepath)"
            println("  Success!")
            push!(successful, taxon)
            sleep(sleep_after)
        catch e
            if occursin("No records found", string(e))
                println("  No occurrence data found for this taxon")
                push!(no_data, taxon)
                sleep(1)
            else
                println("  ERROR: $e")
                push!(failed, taxon)
                sleep(10)
            end
        end
    end

    isempty(no_data) || println("\nTaxa with no occurrence records: $(join(no_data, ", "))")
    isempty(failed) || println("\nFailed taxa (errors): $(join(failed, ", "))")
    isempty(no_data) && isempty(failed) && println("\nAll pulls succeeded (or were skipped).")

    return (successful=successful, no_data=no_data, failed=failed)
end

"""
    append_occurrences(extra_df, raw_dir; taxon_col=:accepted_name) -> Vector{String}

Append extra occurrence records (ex. hand-picked iNaturalist observations) to
the most recently pulled raw CSV of each taxon in `raw_dir`. Records are grouped
by `taxon_col`; duplicated rows are dropped, so this is safe to re-run.

`extra_df` should use the same column names as the raw files. Returns the paths
that were written.
"""
function append_occurrences(extra_df::AbstractDataFrame, raw_dir::String; taxon_col::Symbol=:accepted_name)
    written = String[]

    for subdf in groupby(extra_df, taxon_col)
        taxon = subdf[1, taxon_col]

        # Find the most recently pulled raw file for this taxon
        fpath, _ = latest_raw_file(String(replace(taxon, " " => "_")), raw_dir)
        if isempty(fpath)
            @warn "No raw file found for $taxon — skipping"
            continue
        end

        # Read existing file and append
        existing_df = DataFrame(CSV.File(fpath; missingstring=["", "NA"]))
        combined_df = vcat(existing_df, DataFrame(subdf); cols=:union)
        unique!(combined_df)

        CSV.write(fpath, combined_df)
        println("Appended $(nrow(subdf)) records to $(fpath)")
        push!(written, fpath)
    end

    return written
end
