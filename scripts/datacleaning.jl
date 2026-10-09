using CSV, DataFrames
using RCall

# needs fileutils.jl (synonym_names, latest_raw_file)

R"""
library(gatoRs)
library(ggplot2)
library(sf)
library(ggspatial)
library(dplyr) #needs to be loaded in last, so it defaults to correct `filter` fxn
"""

"""
    find_homonyms_to_exclude(synonymdict, homonyms) -> Dict(taxon => ["name authority", ...])

For each taxon in `synonymdict`, collect the homonyms that belong to a
*different* taxon, so they can be filtered out of that taxon's occurrences
(see `clean_occurrences`). A homonym is a name that is used for more than one
taxon, under different authorities.

Example: "Calceolaria nitida" has the synonym "Calceolaria glandulosa Benth.",
so "Calceolaria glandulosa Poepp. ex Benth." (a different, accepted taxon) is
returned as a name to exclude for "Calceolaria nitida".

# Arguments
- `synonymdict`: taxon name => vector of `(name, authority)` tuples
  (authority can be `missing`).
- `homonyms`: DataFrame with one row per homonym & the columns
  `scientificName` (the shared name), `authority`, `accepted` (is this homonym
  an accepted name?) and `validName` (the accepted name it's a synonym of,
  if `accepted` is false).

A homonym is kept (not excluded) for a taxon if the taxon uses that name under
the same authority, or if the homonym is a name for that taxon — including,
for a species, any of its infraspecific taxa.
"""
function find_homonyms_to_exclude(synonymdict::AbstractDict, homonyms::AbstractDataFrame)
    homonyms_to_exclude = Dict{String, Vector{String}}()

    for (taxon, syns) in synonymdict
        # name -> authorities under which this taxon uses that name
        syn_auths = Dict{String, Vector{String}}()
        for (syn_name, syn_auth) in syns
            auths = get!(syn_auths, syn_name, String[])
            ismissing(syn_auth) || push!(auths, syn_auth)
        end

        for row in eachrow(homonyms)
            haskey(syn_auths, row.scientificName) || continue
            (ismissing(row.authority) || row.authority == "") && continue

            homonym_taxon = row.accepted ? row.scientificName : coalesce(row.validName, "")
            same_authority = row.authority in syn_auths[row.scientificName]
            same_taxon = homonym_taxon == taxon || startswith(homonym_taxon, taxon * " ")
            (same_authority || same_taxon) && continue

            push!(get!(homonyms_to_exclude, taxon, String[]), "$(row.scientificName) $(row.authority)")
        end
    end

    return homonyms_to_exclude
end

"""
    clean_occurrences(synonymdict, raw_dir, clean_dir, no_coords_dir; kwargs...) -> DataFrame

Batch-clean the most recently pulled raw occurrence CSV (see `latest_raw_file`)
of every taxon in `synonymdict`, using `gatoRs`. For each taxon:

1. standardize `basisOfRecord` & drop the record types in `basis_to_remove`
2. `gatoRs::taxa_clean`: keep records matching the taxon's synonym list
3. `gatoRs::remove_duplicates`: drop duplicates between GBIF & iDigBio
4. drop the homonyms in `homonyms_to_exclude` (exact match on name + authority)
5. split off the records with no coordinates (but some locality info), for georeferencing
6. `gatoRs::basic_locality_clean` on the records with coordinates
7. drop points outside the bounding box / in `countries_to_remove`

Writes `clean_dir/<taxon>-<date>_cleaned.csv` and, if there are any,
`no_coords_dir/<taxon>-<date>_no_coords.csv`, where `<date>` is the date of the
raw file. Returns a DataFrame with the number of rows left after each step.

# Arguments
- `synonymdict`: taxon name => its synonym list (a vector of names, or of
  `(name, authority)` tuples). The first name is used as the accepted name.
- `raw_dir`, `clean_dir`, `no_coords_dir`: Input & output directories.
- `homonyms_to_exclude`: taxon => full scientific names (with authority) to
  drop, ex. from `find_homonyms_to_exclude` (default: none).
- `basis_to_remove`: `basisOfRecord` values to drop.
- `lat_min`, `lat_max`, `lon_min`, `lon_max`: Bounding box (default: unbounded).
- `countries_to_remove`: Values of the `country` column to drop (default: none).
- `drop_missing_country`: Also drop records with no `country` (default: false).
- `digits`: Coordinates are rounded to this many decimal places (default: 2).
"""
function clean_occurrences(synonymdict::AbstractDict, raw_dir::String, clean_dir::String, no_coords_dir::String;
    homonyms_to_exclude::AbstractDict=Dict{String, Vector{String}}(),
    basis_to_remove::Vector{String}=["MATERIAL_SAMPLE", "MATERIAL_CITATION", "LIVING_SPECIMEN", "FOSSIL_SPECIMEN"],
    lat_min::Real=-Inf, lat_max::Real=Inf,
    lon_min::Real=-Inf, lon_max::Real=Inf,
    countries_to_remove::Vector{String}=String[],
    drop_missing_country::Bool=false,
    digits::Int=2)

    filtering_stats = DataFrame(
        taxa = String[],
        raw = Int[],
        after_taxa_clean = Int[],
        after_remove_duplicates = Int[],
        removed_homonyms = Int[],
        rows_missing_coords = Int[],
        rows_with_coords = Int[],
        after_locality_filt = Int[],
        after_remove_flagged = Int[]
    )

    mkpath(clean_dir)
    mkpath(no_coords_dir)

    @rput basis_to_remove lat_min lat_max lon_min lon_max countries_to_remove drop_missing_country digits

    for taxa in sort(collect(keys(synonymdict)))
        synlist = synonym_names(synonymdict[taxa])
        filename = replace(taxa, " " => "_")
        # use the most recently pulled raw file for this taxon
        filepath, raw_date = latest_raw_file(String(filename), raw_dir)
        if isempty(filepath)
            @warn "No raw file found for $taxa — skipping"
            continue
        end

        # homonyms (belonging to other taxa) to filter out
        exclude_names = get(homonyms_to_exclude, taxa, String[])
        outfile = joinpath(clean_dir, "$(filename)-$(raw_date)_cleaned.csv")
        outfile_na = joinpath(no_coords_dir, "$(filename)-$(raw_date)_no_coords.csv")

        @rput filepath synlist exclude_names outfile outfile_na

        R"""
        df <- read.csv(filepath)
        nrow_raw <- nrow(df)

        # standardize and filter basis of record
        df$basisOfRecord[df$basisOfRecord %in% c("preservedspecimen", "PreservedSpecimen", "Preserved Specimen")] <- "PRESERVED_SPECIMEN"
        df <- df[!df$basisOfRecord %in% basis_to_remove, ]

        df <- taxa_clean(df,
                         synonyms.list = synlist,
                         taxa.filter = "fuzzy",
                         accepted.name = synlist[1])
        nrow_taxa_clean <- nrow(df)

        # remove duplicates between GBIF and iDigBio
        df <- remove_duplicates(df,
                                remove.unparseable = TRUE)
        nrow_remove_dup <- nrow(df)

        # filter out homonyms that belong to a different taxon (matching on the full name + authority)
        df_before <- nrow(df)
        df <- df[!(trimws(df$scientificName) %in% exclude_names), ]
        nrow_removed_homonyms <- df_before - nrow(df)

        # separate rows with NA lat/long for potential georeferencing
        df_na <- df[is.na(df$latitude) | is.na(df$longitude), ]
        # filter out rows where absolutely no locality info
        df_na <- df_na[df_na$locality != "locality:  NA, occurrenceRemarks: NA, verbatimLocality: NA", ]
        df_clean <- df[!is.na(df$latitude) & !is.na(df$longitude), ]

        nrow_na <- nrow(df_na)
        nrow_with_coords <- nrow(df_clean)

        df_clean <- basic_locality_clean(df_clean,
                                    remove.zero = TRUE,
                                    precision = TRUE,
                                    digits = digits,
                                    remove.skewed = TRUE)
        nrow_locality <- nrow(df_clean)

        # remove points outside the lat/long bounding box
        df_clean <- df_clean %>%
             dplyr::filter(latitude  >= lat_min, latitude  <= lat_max,
                           longitude >= lon_min, longitude <= lon_max)
        # remove points from unwanted countries
        df_clean <- df_clean %>% dplyr::filter(!(country %in% countries_to_remove))
        if (drop_missing_country) {
            df_clean <- df_clean %>% dplyr::filter(!is.na(country))
        }
        nrow_flagged <- nrow(df_clean)

        write.csv(df_clean, outfile, row.names = FALSE)
        # save rows with NA coordinates (if any exist)
        if (nrow_na > 0) write.csv(df_na, outfile_na, row.names = FALSE)

        rm(df, df_clean, df_na)
        """

        @rget nrow_raw nrow_taxa_clean nrow_remove_dup nrow_removed_homonyms nrow_na nrow_with_coords nrow_locality nrow_flagged

        push!(filtering_stats, (taxa, nrow_raw, nrow_taxa_clean, nrow_remove_dup,
                                nrow_removed_homonyms,
                                nrow_na, nrow_with_coords, nrow_locality, nrow_flagged))

        println("  $taxa: Raw=$nrow_raw → Taxa clean=$nrow_taxa_clean → Remove dup=$nrow_remove_dup")
        if nrow_removed_homonyms > 0
            println("    ⚠️  Removed $nrow_removed_homonyms rows of homonyms: $(join(exclude_names, "; "))")
        end
        println("    Split: $nrow_with_coords with coords, $nrow_na missing coords")
    end

    return filtering_stats
end

"""
    merge_georef(clean_df, no_coords_df, georef_df) -> DataFrame

Add georeferenced records to a taxon's cleaned occurrences. `georef_df` is the
viable georeferenced points (from `load_georef_df`: `ID`, `latitude`,
`longitude` & optionally a notes column); the matching rows of `no_coords_df`
(the records that were pulled without coordinates) are given those coordinates
and appended to `clean_df`, with any notes in a `georef_notes` column.
"""
function merge_georef(clean_df::AbstractDataFrame, no_coords_df::AbstractDataFrame, georef_df::AbstractDataFrame)
    notes_col = findfirst(c -> uppercase(c) == "NOTES", names(georef_df))
    georef_lookup = Dict(
        string(row.ID) => (lat=row.latitude, lon=row.longitude, notes=isnothing(notes_col) ? missing : row[notes_col])
        for row in eachrow(georef_df)
    )

    georef_rows = filter(row -> haskey(georef_lookup, string(row.ID)), no_coords_df)
    georef_rows[!, :latitude] = [georef_lookup[string(id)].lat for id in georef_rows.ID]
    georef_rows[!, :longitude] = [georef_lookup[string(id)].lon for id in georef_rows.ID]
    georef_rows[!, :georef_notes] = Union{Missing,String}[ismissing(georef_lookup[string(id)].notes) ? missing : string(georef_lookup[string(id)].notes) for id in georef_rows.ID]

    merged = copy(clean_df)
    merged[!, :georef_notes] = Vector{Union{Missing,String}}(missing, nrow(merged))
    return vcat(merged, georef_rows; cols=:union)
end

"""
    carry_over_georef(geolocate_df, georef_df; cols) -> (DataFrame, n_filled)

Fill in the rows of a GEOLocate-format file (from `prepare_geolocate_files`)
that have already been georeferenced somewhere else — ex. a species-level file,
from the files of its infraspecific taxa. Rows are matched on `ID`, and only
rows of `georef_df` (from `read_georef_file`) with `VIABLE` filled in (TRUE or
FALSE) are used. `cols` are the columns copied over; `NOTES` & `VIABLE` are
added to `geolocate_df` if it doesn't have them.
"""
function carry_over_georef(geolocate_df::AbstractDataFrame, georef_df::AbstractDataFrame;
    cols=["latitude", "longitude", "correction status", "precision", "error polygon", "multiple results", "NOTES", "VIABLE"])

    done = Dict{String,DataFrameRow}()
    for row in eachrow(georef_df)
        (ismissing(row.ID) || !(row.VIABLE in ("TRUE", "FALSE"))) && continue
        get!(done, string(row.ID), row)
    end

    filled = copy(geolocate_df)
    cols = filter(c -> c in names(georef_df), cols)
    for c in cols
        filled[!, c] = c in names(filled) ?
            Union{Missing,String}[ismissing(v) ? missing : string(v) for v in filled[!, c]] :
            Vector{Union{Missing,String}}(missing, nrow(filled))
    end

    n_filled = 0
    for (i, id) in enumerate(filled.ID)
        row = get(done, string(id), nothing)
        isnothing(row) && continue
        for c in cols
            filled[i, c] = ismissing(row[c]) ? missing : string(row[c])
        end
        n_filled += 1
    end
    return filled, n_filled
end

# ─────────────────────────────────────────────────────────────────────────────
# Filters. Each `*_mask` function returns a vector that's `true` for the rows
# to keep; `remove_rows` applies a mask & keeps a record of what was removed.
# ─────────────────────────────────────────────────────────────────────────────

"""
    coords_mask(df; lat_min, lat_max, lon_min, lon_max) -> Vector{Bool}

`true` for rows inside the lat/lon bounding box. All bounds are optional —
specify only the axes you want to constrain. Rows with missing coordinates are
`false` when a bound is supplied for that axis.
"""
function coords_mask(df::AbstractDataFrame;
    lat_min=nothing, lat_max=nothing,
    lon_min=nothing, lon_max=nothing)

    if all(isnothing, (lat_min, lat_max, lon_min, lon_max))
        error("At least one of lat_min, lat_max, lon_min, lon_max must be specified")
    end

    function keep(row)
        lat = row.latitude
        lon = row.longitude
        if !isnothing(lat_min) || !isnothing(lat_max)
            (ismissing(lat) || !isa(lat, Number)) && return false
            !isnothing(lat_min) && lat < lat_min && return false
            !isnothing(lat_max) && lat > lat_max && return false
        end
        if !isnothing(lon_min) || !isnothing(lon_max)
            (ismissing(lon) || !isa(lon, Number)) && return false
            !isnothing(lon_min) && lon < lon_min && return false
            !isnothing(lon_max) && lon > lon_max && return false
        end
        return true
    end

    return Bool[keep(row) for row in eachrow(df)]
end

"""
    countries_mask(df, countries; country_col=:country) -> Vector{Bool}

`false` for rows from any of `countries` (or any other regional category in
`country_col`). Rows with no country are kept.
"""
function countries_mask(df::AbstractDataFrame, countries::Vector{String}; country_col::Symbol=:country)
    string(country_col) in names(df) || error("Column `$country_col` not found in dataframe")
    return Bool[ismissing(x) || !(x in countries) for x in df[!, country_col]]
end

"""
    names_mask(df, names_to_remove) -> Vector{Bool}

`false` for rows whose `scientificName` is one of `names_to_remove`
(matching on the full string, so including the authority if there is one).
"""
function names_mask(df::AbstractDataFrame, names_to_remove::Vector{String})
    remove_set = Set(names_to_remove)
    return Bool[ismissing(x) || !(x in remove_set) for x in df.scientificName]
end

"""
    remove_rows(df, keep, reason, removed) -> DataFrame

Drop the rows of `df` where `keep` is false. The dropped rows are pushed onto
`removed` (a vector of DataFrames) with a `removed_reason` column, so there's a
record of every point that was filtered out & why.
"""
function remove_rows(df::AbstractDataFrame, keep::AbstractVector{Bool}, reason::AbstractString, removed::Vector{DataFrame})
    all(keep) && return df
    gone = df[.!keep, :]
    gone[!, :removed_reason] .= String(reason)
    push!(removed, gone)
    return df[keep, :]
end

"""
    plot_kept_removed(kept_df, removed_df; title="")

Map the points that were kept (green) & removed (red ×) by a filter.
"""
function plot_kept_removed(kept_df::AbstractDataFrame, removed_df::AbstractDataFrame; title="")
    plot_title = String(title)
    kept_lat = Vector{Union{Missing,Float64}}(kept_df.latitude)
    kept_lon = Vector{Union{Missing,Float64}}(kept_df.longitude)
    removed_lat = Vector{Union{Missing,Float64}}(removed_df.latitude)
    removed_lon = Vector{Union{Missing,Float64}}(removed_df.longitude)
    @rput plot_title kept_lat kept_lon removed_lat removed_lon
    R"""
    kept_pts    <- data.frame(latitude  = as.numeric(kept_lat),
                              longitude = as.numeric(kept_lon))
    removed_pts <- data.frame(latitude  = as.numeric(removed_lat),
                              longitude = as.numeric(removed_lon))
    all_lat <- c(kept_pts$latitude,  removed_pts$latitude)
    all_lon <- c(kept_pts$longitude, removed_pts$longitude)
    world    <- annotation_borders(database="world", colour="gray80", fill="gray80")
    borders  <- annotation_borders(database="world", colour="gray40", fill=NA, size=0.5)
    p <- ggplot() + world + borders +
        geom_point(data=kept_pts,
                   aes(x=longitude, y=latitude),
                   color="darkgreen", size=2, alpha=0.7) +
        { if (nrow(removed_pts) > 0)
              geom_point(data=removed_pts,
                         aes(x=longitude, y=latitude),
                         color="red", shape=4, size=3.5, stroke=1.2, alpha=0.9)
          else NULL } +
        coord_sf(xlim=c(min(all_lon, na.rm=TRUE)-2, max(all_lon, na.rm=TRUE)+2),
                 ylim=c(min(all_lat, na.rm=TRUE)-2, max(all_lat, na.rm=TRUE)+2)) +
        labs(title=paste0(plot_title, "  —  green: kept (", nrow(kept_pts),
                          ")   red ×: removed (", nrow(removed_pts), ")"),
             x="Longitude", y="Latitude") +
        theme_minimal() +
        theme(plot.title=element_text(size=10))
    print(p)
    """
    return nothing
end

"""
    filter_coords(input_path::String, output_path::String; kwargs...)

Remove rows outside the specified lat/lon bounding box (see `coords_mask`).
Always plots a preview map (green = kept, red × = removed).
Set `save=true` to write the filtered result to `output_path`.

# Arguments
- `input_path`: Path to the input CSV file.
- `output_path`: Path to save the filtered CSV file.
- `taxon_label`: Label for the plot title (default: "Taxon").
- `lat_min`, `lat_max`, `lon_min`, `lon_max`: Coordinate bounds.
- `save`: Boolean to indicate whether to save the output (default: false).
"""
function filter_coords(input_path::String, output_path::String;
    taxon_label="Taxon",
    lat_min=nothing, lat_max=nothing,
    lon_min=nothing, lon_max=nothing,
    save::Bool=false)

    df = DataFrame(CSV.File(input_path))
    keep_mask = coords_mask(df; lat_min=lat_min, lat_max=lat_max, lon_min=lon_min, lon_max=lon_max)
    kept_df = df[keep_mask, :]
    removed_df = df[.!keep_mask, :]

    println("  Preview: $(nrow(removed_df)) point(s) will be removed, $(nrow(kept_df)) will be kept")
    plot_kept_removed(kept_df, removed_df; title=taxon_label)

    if !save
        println(" Preview only — re-run with save=true to apply.")
        return nothing
    end

    CSV.write(output_path, kept_df)
    println(" Written: $output_path")
end

"""
    filter_countries(input_path::String, output_path::String, countries::Vector{String}; country_col::Symbol=:country, save::Bool=true)

Remove rows from a dataset that match specific countries (or any other regional
category; see `countries_mask`). Rows with no country are kept.

# Arguments
- `input_path`: Path to the input CSV file.
- `output_path`: Path to save the filtered CSV file.
- `countries`: A vector of country names to filter out.
- `country_col`: The symbol of the column containing country data (default: `:country`).
- `save`: Boolean to indicate whether to save the output (default: true).
"""
function filter_countries(input_path::String, output_path::String, countries::Vector{String}; country_col::Symbol=:country, save::Bool=true)
    df = DataFrame(CSV.File(input_path))
    keep_mask = countries_mask(df, countries; country_col=country_col)

    if all(keep_mask)
        @warn "No input country provided in `countries` list matches countries existing in dataframe"
    end
    df = df[keep_mask, :]

    if save
        CSV.write(output_path, df)
        println(" Written: $output_path")
    end

    return df
end

"""
    prepare_geolocate_files(input_dir::String, output_dir::String)
 
Read every CSV in `input_dir`, reformat columns to GEOLocate batch input
format, and write the result to `output_dir` (same filename).
 
Expected input columns: locality, country, stateProvince, county, latitude,
longitude, ID, scientificName, basisOfRecord.
 
Output columns (GEOLocate format): "locality string", country, state, county,
latitude, longitude, "correction status", precision, "error polygon",
"multiple results", ID, name, basis.
"""
function prepare_geolocate_files(
    input_dir::String,
    output_dir::String
)
    isdir(output_dir) || mkpath(output_dir)

    files = DataFrames.filter(f -> endswith(f, ".csv"), readdir(input_dir))

    if isempty(files)
        println("No CSV files found in $input_dir")
        return 0
    end

    converted_count = 0

    for (idx, file) in enumerate(files)
        println("Processing $idx/$(length(files)): $file")

        input_file = joinpath(input_dir, file)
        output_file = joinpath(output_dir, file)

        @rput input_file output_file

        R"""
        rawdf_GeoRef <- read.csv(input_file)

        if (nrow(rawdf_GeoRef) == 0) {
            cat("  Skipping - empty file\n")
            next
        }

        rawdf_GeoRef <- rawdf_GeoRef %>%
            dplyr::select("locality string" = locality,
                          country,
                          state = stateProvince,
                          county,
                          latitude,
                          longitude,
                          ID,
                          name = scientificName,
                          basis = basisOfRecord)

        rawdf_GeoRef$'correction status' <- ""
        rawdf_GeoRef$precision           <- ""
        rawdf_GeoRef$'error polygon'     <- ""
        rawdf_GeoRef$'multiple results'  <- ""

        rawdf_GeoRef2 <- rawdf_GeoRef[, c("locality string", "country",
                                           "state", "county", "latitude",
                                           "longitude", "correction status",
                                           "precision", "error polygon",
                                           "multiple results", "ID",
                                           "name", "basis")]

        write.csv(rawdf_GeoRef2, output_file, row.names = FALSE)
        """

        converted_count += 1
    end

    println("\nDone. Converted $converted_count files → $output_dir")
    return converted_count
end

"""
    filter_scientific_names(input_path::String, output_path::String, names_to_remove::Vector{String}; save::Bool=true)

Remove rows whose `scientificName` matches any value in `names_to_remove`.
Prints a count of removed rows.

# Example
```julia
filter_scientific_names("clean/Calceolaria_undulata-2026_10_06_cleaned.csv", "filtered/Calceolaria_undulata-2026_10_06_cleaned.csv",
    ["Calceolaria foliosa", "Calceolaria foliosa Meyen ex Walp. & Schauer"])
```
"""
function filter_scientific_names(input_path::String, output_path::String, names_to_remove::Vector{String}; save::Bool=true)
    df = DataFrame(CSV.File(input_path; missingstring=["", "NA"]))
    keep_mask = names_mask(df, names_to_remove)

    n_removed = count(!, keep_mask)
    if n_removed == 0
        println("  No rows matched the given scientificName values — nothing removed.")
        return
    end

    kept_df = df[keep_mask, :]
    names_str = join(names_to_remove, ", ")
    println("  Removing $n_removed row(s) with scientificName in: $names_str")

    if save
        CSV.write(output_path, kept_df)
        println("  Written: $output_path")
    else
        println("  Preview only — re-run with save=true to apply.")
    end
end
