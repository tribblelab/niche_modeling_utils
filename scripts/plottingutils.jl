using CSV, DataFrames
using RCall

# needs fileutils.jl (load_taxa_data, load_georef_df, latest_raw_file, resolve_clean_file)

R"""
library(ggplot2)
library(sf)
library(ggspatial)
library(gridExtra)
"""

"""
    load_bot_regions(shapefile_path)

Load TDWG botanical country polygons into R (once, globally).
"""
function load_bot_regions(shapefile_path="data/bot_country_shapefiles/level3.shp")
    @rput shapefile_path
    R"bot_regions = st_read(shapefile_path, quiet = TRUE)"
end

"""
    push_species_to_r(taxa, raw_file, clean_file, country_codes, georef_df)

Transfer all data needed by the R plotting block to R's environment.
"""
function push_species_to_r(taxa, raw_file, clean_file, country_codes, georef_df, only_clean=false)
    @rput taxa raw_file clean_file country_codes only_clean
    has_georef = !isnothing(georef_df) && nrow(georef_df) > 0
    @rput has_georef
    if has_georef
        georef_ids = string.(georef_df.ID)
        @rput georef_ids
    end
end

# Core R plotting block (rendered into a grid.arrange object in R)

const RPLOT_BLOCK = """
    raw_df  <- read.csv(raw_file)
    clean_df <- read.csv(clean_file)
    raw_df  <- raw_df[!is.na(raw_df\$latitude) & !is.na(raw_df\$longitude), ]

    if (nrow(raw_df) == 0 || nrow(clean_df) == 0) stop("empty data")

    raw_lon_min <- min(raw_df\$longitude, na.rm=TRUE) - 3
    raw_lon_max <- max(raw_df\$longitude, na.rm=TRUE) + 3
    raw_lat_min <- min(raw_df\$latitude,  na.rm=TRUE) - 3
    raw_lat_max <- max(raw_df\$latitude,  na.rm=TRUE) + 3

    all_clean_lon <- clean_df\$longitude
    all_clean_lat <- clean_df\$latitude
    clean_lon_min <- min(all_clean_lon, na.rm=TRUE) - 2
    clean_lon_max <- max(all_clean_lon, na.rm=TRUE) + 2
    clean_lat_min <- min(all_clean_lat, na.rm=TRUE) - 2
    clean_lat_max <- max(all_clean_lat, na.rm=TRUE) + 2

    world    <- annotation_borders(database="world", colour="gray80", fill="gray80")
    countries <- annotation_borders(database="world", colour="gray40", fill=NA, size=0.5)

    bot_underlay <- bot_regions[bot_regions\$LEVEL3_COD %in% country_codes, ]
    underlay_layer <- if (nrow(bot_underlay) > 0) {
        geom_sf(data=bot_underlay, fill="khaki1", color="goldenrod4",
                alpha=0.25, linewidth=0.3)
    } else { NULL }

    p1 <- ggplot() +
        world + countries + underlay_layer +
        geom_point(data=raw_df, aes(x=longitude, y=latitude),
                   color="blue", size=1.5, alpha=0.6) +
        coord_sf(xlim=c(raw_lon_min, raw_lon_max),
                 ylim=c(raw_lat_min, raw_lat_max)) +
        labs(title=paste0("Raw (n=", nrow(raw_df), ")"),
             x="Longitude", y="Latitude") +
        theme_minimal() +
        theme(plot.title=element_text(size=10))

    clean_title <- paste0("Cleaned (n=", nrow(clean_df), ")")

    clean_df\$plot_category <- "PRESERVED_SPECIMEN"
    if ("basisOfRecord" %in% names(clean_df)) {
        clean_df\$plot_category[clean_df\$basisOfRecord == "HUMAN_OBSERVATION"] <- "HUMAN_OBSERVATION"
    }
    if (has_georef) {
        is_georefed <- as.character(clean_df\$ID) %in% as.character(georef_ids)
        clean_df\$plot_category[is_georefed & clean_df\$plot_category == "PRESERVED_SPECIMEN"] <- "PRESERVED_SPECIMEN (georeferenced)"
    }
    
    plot_colors <- c("HUMAN_OBSERVATION" = "blue",
                     "PRESERVED_SPECIMEN" = "darkgreen",
                     "PRESERVED_SPECIMEN (georeferenced)" = "darkorange")

    p2 <- ggplot() +
        world + countries + underlay_layer +
        geom_point(data=clean_df, aes(x=longitude, y=latitude, color=plot_category),
                   size=1.5, alpha=0.8) +
        scale_color_manual(values=plot_colors) +
        coord_sf(xlim=c(clean_lon_min, clean_lon_max),
                 ylim=c(clean_lat_min, clean_lat_max)) +
        labs(title=clean_title, x="Longitude", y="Latitude", color="Record Type") +
        annotation_scale(location="bl") +
        annotation_north_arrow(location="tl",
                               height=unit(0.8,"cm"), width=unit(0.8,"cm")) +
        theme_minimal() +
        theme(plot.title=element_text(size=10),
              legend.position="bottom",
              legend.title=element_text(size=9),
              legend.text=element_text(size=8))

    if (only_clean) {
        combined <- grid.arrange(p2, ncol=1, top=taxa)
    } else {
        combined <- grid.arrange(p1, p2, ncol=2, top=taxa)
    }
    print(combined)
"""

"""
    plot_species(taxa;
        traits_path,
        raw_dir,
        clean_dir,
        georef_dir,
        filtered_dir   = nothing,
        shapefile_path = "data/bot_country_shapefiles/level3.shp",
        append_date    = true,
        species_ranges_path = nothing)

Display a before/after occurrence map for a single `taxa` string in the R
graphics viewer (no PDF output).

If `append_date` is true, the raw file is the most recently pulled
`<taxon>-<date>.csv` in `raw_dir`; otherwise it's `<taxon>.csv`.
`species_ranges_path` is passed on to `load_taxa_data`.
If `filtered_dir` is given & has a file for the taxon, those points are plotted
instead of the ones in `clean_dir`.

# Example
```julia
plot_species("Calceolaria alba")
```
"""
function plot_species(taxa::String;
    traits_path::String,
    raw_dir::String,
    clean_dir::String,
    georef_dir::String,
    filtered_dir::Union{String,Nothing}=nothing,
    shapefile_path::String="data/bot_country_shapefiles/level3.shp",
    append_date::Bool=true,
    species_ranges_path::Union{String,Nothing}=nothing,
    only_clean=false,
    min_points=0
)
    _, taxa_to_nativerange_dict = load_taxa_data(traits_path; species_ranges_path=species_ranges_path)
    load_bot_regions(shapefile_path)

    filename = replace(taxa, " " => "_")
    raw_file = append_date ? first(latest_raw_file(filename, raw_dir)) : joinpath(raw_dir, "$(filename).csv")
    clean_file = resolve_clean_file(filename, clean_dir)

    isfile(raw_file) || error("Raw file not found: $raw_file")
    isempty(clean_file) && error("No cleaned file found for: $taxa")

    # plot the filtered version of the cleaned points, if there is one
    filtered_file = isnothing(filtered_dir) ? "" : resolve_clean_file(filename, filtered_dir)
    if !isempty(filtered_file)
        clean_file = filtered_file
    end

    raw_df = filter(
        row -> !ismissing(row.latitude) && !ismissing(row.longitude),
        DataFrame(CSV.File(raw_file))
    )
    clean_df = DataFrame(CSV.File(clean_file))

    nrow(raw_df) > 0 || error("No valid coordinates in raw file for $taxa")
    nrow(clean_df) > 0 || error("No rows in cleaned file for $taxa")
    nrow(clean_df) >= min_points || error("Only $(nrow(clean_df)) points in cleaned file, which is < min_points=$min_points")

    georef_df = load_georef_df(filename, georef_dir)
    country_codes = get(taxa_to_nativerange_dict, taxa, String[])

    push_species_to_r(taxa, raw_file, clean_file, country_codes, georef_df, only_clean)
    reval(RPLOT_BLOCK)

    println("Plot displayed for: $taxa")
    return nothing
end

"""
    plot_all_occurrence_maps(;
        traits_path,
        raw_dir,
        clean_dir,
        georef_dir,
        filtered_dir   = nothing,
        pdf_file,
        shapefile_path = "data/bot_country_shapefiles/level3.shp",
        append_date    = true,
        species_ranges_path = nothing)

Iterate over every taxon in `traits_path` (plus those in `species_ranges_path`,
if given), build a before/after occurrence map for each, and write all pages to
`pdf_file`.  Returns the number of taxa successfully plotted.
If `filtered_dir` is given, a taxon's file there is plotted instead of the one in `clean_dir`.
"""
function plot_all_occurrence_maps(;
    traits_path::String,
    raw_dir::String,
    clean_dir::String,
    georef_dir::String,
    filtered_dir::Union{String,Nothing}=nothing,
    shapefile_path::String="data/bot_country_shapefiles/level3.shp",
    pdf_file::String,
    append_date::Bool=true,
    species_ranges_path::Union{String,Nothing}=nothing,
    only_clean=false,
    min_points=15
)
    taxa_traits, taxa_to_nativerange_dict = load_taxa_data(traits_path; species_ranges_path=species_ranges_path)
    load_bot_regions(shapefile_path)

    @rput pdf_file
    R"pdf(pdf_file, width=14, height=7)"

    plotted_count = 0

    for (idx, taxa) in enumerate(taxa_traits.scientificName)
        println("Processing $idx/$(length(taxa_traits.scientificName)): $taxa")

        filename = replace(taxa, " " => "_")
        raw_file = append_date ? first(latest_raw_file(filename, raw_dir)) : joinpath(raw_dir, "$(filename).csv")
        clean_file = resolve_clean_file(filename, clean_dir)

        if !isfile(raw_file) || isempty(clean_file)
            println("  Skipping - files not found")
            continue
        end

        # plot the filtered version of the cleaned points, if there is one
        filtered_file = isnothing(filtered_dir) ? "" : resolve_clean_file(filename, filtered_dir)
        if !isempty(filtered_file)
            clean_file = filtered_file
        end

        raw_df = filter(
            row -> !ismissing(row.latitude) && !ismissing(row.longitude),
            DataFrame(CSV.File(raw_file))
        )
        clean_df = DataFrame(CSV.File(clean_file))

        if nrow(raw_df) == 0 || nrow(clean_df) == 0
            println("  Skipping - no valid coordinates")
            continue
        end

        if nrow(clean_df) < min_points
            println("  Skipping - fewer than $min_points points ($(nrow(clean_df)))")
            continue
        end

        georef_df = load_georef_df(filename, georef_dir)
        country_codes = get(taxa_to_nativerange_dict, taxa, String[])

        push_species_to_r(taxa, raw_file, clean_file, country_codes, georef_df, only_clean)
        reval(RPLOT_BLOCK)
        plotted_count += 1
    end

    R"dev.off()"
    println("\nDone. Plotted $plotted_count taxa → $pdf_file")
    return plotted_count
end