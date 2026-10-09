# niche_modeling_utils

## to instantiate Julia pkgs for the first time:

to instantiate / load packages from `\niche_modeling_utils` directory:

`julia --project`

```julia
using Pkg
# install the packages listed in the environment
Pkg.instantiate()
```

## R package dependencies + install

```R
packages <- c("gatoRs", "ggplot2", "sf", "ggspatial", "gridExtra", "CoordinateCleaner", "readxl", "dplyr")
new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)
```

## what's in `scripts/`

there's no package (yet), so `include` the files you need. `fileutils.jl` goes first, since the others use it:

```julia
include("niche_modeling_utils/scripts/fileutils.jl")
include("niche_modeling_utils/scripts/occpulling.jl")
include("niche_modeling_utils/scripts/datacleaning.jl")
include("niche_modeling_utils/scripts/plottingutils.jl")
```

| file | functions |
|---|---|
| `checklistutils.jl` | `build_checklist`: build a WCVP/POWO name checklist for a family + genus; then `checklist_synonyms`, `checklist_native_range` & `checklist_homonyms` to get each taxon's synonyms, native range & homonyms out of it |
| `occpulling.jl` | `pull_occurrences` (GBIF / iDigBio, via gatoRs), `append_occurrences` |
| `fileutils.jl` | find & load per-taxon files: `latest_raw_file`, `resolve_clean_file`, `latest_no_coords_file`, `taxon_files`, `taxon_stem`, `read_occs`, `normalize_id`, `load_taxa_data`, `read_georef_file`, `load_georef_df`, `load_pt_occs_df` |
| `datacleaning.jl` | `find_homonyms_to_exclude`, `clean_occurrences` (batch gatoRs cleaning), `prepare_geolocate_files`, `merge_georef`, `carry_over_georef`; filters as masks (`coords_mask`, `countries_mask`, `names_mask`, applied & logged with `remove_rows`) or file to file (`filter_coords`, `filter_countries`, `filter_scientific_names`); `plot_kept_removed` |
| `plottingutils.jl` | `plot_species`, `plot_all_occurrence_maps` |

the functions expect one file per taxon, named with the taxon (spaces → underscores) & the date it was pulled:

- raw occurrences: `<taxon>-<yyyy_mm_dd>.csv`
- cleaned occurrences: `<taxon>-<yyyy_mm_dd>_cleaned.csv`
- occurrences with no coordinates: `<taxon>-<yyyy_mm_dd>_no_coords.csv`
- final occurrences (cleaned + georeferenced + filtered): `<taxon>.csv`

a typical run:

```julia
# taxon => its synonyms (names, or (name, authority) tuples)
pull_occurrences(synonymdict, "data/pt_occs_raw")
stats = clean_occurrences(synonymdict, "data/pt_occs_raw", "data/pt_occs_clean", "data/pt_occs_to_georeference")
prepare_geolocate_files("data/pt_occs_to_georeference", "data/pt_occs_to_georeference_geolocate_format")
```
