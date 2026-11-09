## ============================================================
##  ##  Shared data preparation for richness and model analyses
##
##  Produces:
##    df_full      — full cleaned dataset with grid cells
##    df           — temporally replicated cells, event_id filtered, deduplicated
##    shared_cells — character vector of cell IDs sampled in ≥ 2 years
##    keep_species — character vector of species meeting modelling
##                   thresholds (used to align richness analysis)
##    sp_cell_year — aggregated effort/response table for models
##    GRID_M       — grid cell size in metres
##
## ============================================================

setwd("~/Library/CloudStorage/OneDrive-MWLR/Marion Donald Rutherford Postdoc - Methods, Data and Analysis/Data_working/output/data")

library(tidyverse)
library(sf)

## ── 1. Load  ────────────────────────────────────────

GRID_M <- 50000L
df_full <- read.csv("final_df_for_modelv2.csv")
## Confirm only the seven intended families remain before continuing.
stopifnot(all(unique(df_full$family) %in% FAMILIES))
cat("Families present after filtering:\n")
print(sort(unique(df_full$family)))
cat("Year range:", range(df_full$year), "\n\n")

## ── 2. NZTM grid cells ───────────────────────────────────────

xy             <- st_coordinates(
  st_transform(
    st_as_sf(df_full, coords = c("longitude", "latitude"), crs = 4326),
    2193))
df_full$cell_e <- floor(xy[, "X"] / GRID_M) * GRID_M
df_full$cell_n <- floor(xy[, "Y"] / GRID_M) * GRID_M
df_full$cell   <- paste(df_full$cell_e, df_full$cell_n, sep = "_")

## ── 3. Temporally replicated cells ───────────────────────────
##  Restrict to cells sampled in ≥ 2 distinct years.

shared_cells <- df_full |>
  group_by(cell) |>
  dplyr::summarise(n_years = n_distinct(year), .groups = "drop") |>
  filter(n_years >= 2) |>
  pull(cell)
cat("Temporally replicated cells:", length(shared_cells), "\n")

## ── 4. Clean df ──────────────────────────────────────────────
##
##  Filter to temporally replicated cells, drop year-only records
##  (no event_id), and deduplicate to one record per species per
##  collection event.

df <- df_full |>
  filter(cell %in% shared_cells) |>
  filter(!is.na(event_id)) |>
  distinct(event_id, species, cell,
           cell_e, cell_n, latitude, longitude, year,
           .keep_all = TRUE) 

cat("Records after dedup:", nrow(df), "\n")
cat("Years:  ", n_distinct(df$year),    "\n")
cat("Species:", n_distinct(df$species), "\n")
cat("Cells:  ", n_distinct(df$cell),    "\n\n")

## ── 5. Effort and aggregated response ────────────────────────
##
##  sp_cell_year is the unit of analysis for the GLMER models:
##  species × cell × year, with total effort (n records that
##  year in that cell) used to compute n_other.

effort_df <- df |>
  group_by(cell, year) |>
  summarise(effort = n(), .groups = "drop")

sp_cell_year <- df |>
  group_by(species, cell, year) |>     
  summarise(n_sp = n(), .groups = "drop") |>
  left_join(effort_df, by = c("cell", "year")) |>
  mutate(n_other = effort - n_sp,
         obs_id  = factor(seq_len(n())))

stopifnot(all(sp_cell_year$n_sp <= sp_cell_year$effort))

## ── 6. Species filtering ─────────────────────────────────────
##  Thresholds:
##    ≥ 30 records across ≥ 20 distinct years, ≥ 30-year span,
##      first record no later than 1948
##    ≥ 30 records total, across ≥ 15 cells
##    ≥ 30 rows in sp_cell_year (ensures enough cell × year
##      combinations to fit random effects)

sp_period_counts <- df |>
  group_by(species) |>
  summarise(
    n         = n(),
    n_years   = n_distinct(year),
    year_span = max(year) - min(year),
    year_min  = min(year),
    .groups   = "drop"
  ) |>
  filter(n >= 30, n_years >= 20, year_span >= 30, year_min <= 1948)

sp_spatial_counts <- df |>
  group_by(species) |>
  summarise(n_records = n(), n_cells = n_distinct(cell), .groups = "drop") |>
  filter(n_records >= 30, n_cells >= 15)

sp_row_counts <- sp_cell_year |>
  group_by(species) |>
  summarise(n_rows = n(), .groups = "drop") |>
  filter(n_rows >= 30)

keep_species <- Reduce(intersect, list(
  sp_period_counts$species,
  sp_spatial_counts$species,
  sp_row_counts$species
)) |> unique()

cat("Species retained:", length(keep_species), "\n")
cat("Species dropped: ", n_distinct(df$species) - length(keep_species), "\n")

cat("\n── data prep complete ───────────────────────────────────\n")
cat("Objects available: df_full, df, shared_cells,\n")
cat("                   keep_species, sp_cell_year, GRID_M\n")

## ── 7. Species summary (diagnostic / supplementary table) ────

species_summary <- sp_period_counts |>
  dplyr::select(species, n, n_years, year_span) |>
  # include all species, not just those passing the filter
  right_join(
    df |> distinct(species),
    by = "species"
  ) |>
  left_join(
    sp_spatial_counts |> dplyr::select(species, n_records, n_cells),
    by = "species"
  ) |>
  left_join(
    sp_row_counts |> dplyr::select(species, n_rows),
    by = "species"
  ) |>
  mutate(
    pass_temporal  = n >= 30 & n_years >= 20 & year_span >= 30,
    pass_spatial   = n_records >= 30 & n_cells >= 15,
    pass_rows      = n_rows >= 30,
    retained       = species %in% keep_species
  ) |>
  arrange(desc(retained), desc(n))

species_summary <- species_summary |>
  mutate(
    pass_temporal = coalesce(pass_temporal, FALSE),
    pass_spatial  = coalesce(pass_spatial,  FALSE),
    pass_rows     = coalesce(pass_rows,     FALSE)
  )

species_summary |>
  filter(!retained) |>
  summarise(
    fails_temporal_only = sum(!pass_temporal &  pass_spatial &  pass_rows),
    fails_spatial_only  = sum( pass_temporal & !pass_spatial &  pass_rows),
    fails_rows_only     = sum( pass_temporal &  pass_spatial & !pass_rows),
    fails_multiple      = sum((!pass_temporal + !pass_spatial + !pass_rows) > 1),
    fails_all           = sum(!pass_temporal & !pass_spatial & !pass_rows)
  ) |>
  print()

print(species_summary, n = Inf)

species_summary |>
  filter(retained) |>
  dplyr::select(species, n, n_years, year_span, n_cells, n_rows) |>
  arrange(n) |>
  print()


## 1. Collated dataset (before cell filter, dedup, event filter)
df_full |> summarise(records = n(), species = n_distinct(species),
                     year_min = min(year), year_max = max(year))

## 2. Bee vs fly split, analytical dataset
BEE_FAM <- c("Colletidae", "Apidae", "Halictidae", "Megachilidae")
FLY_FAM <- c("Syrphidae", "Calliphoridae", "Pollendiae")
df <- df |> mutate(group = if_else(family %in% BEE_FAM, "bee", "fly"))
df |> group_by(group) |> summarise(records = n(), species = n_distinct(species))
df |> filter(species %in% keep_species) |> distinct(species, group) |> count(group)

