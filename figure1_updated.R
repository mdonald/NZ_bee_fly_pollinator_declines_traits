## ============================================================
##  figure1_sampling.R
##  Figure 1: sampling coverage and temporal distribution of bee and
##  fly records, rebuilt from the full collated dataset of the five
##  focal families, BEFORE deduplication, the two-year cell rule and
##  species filtering (df_full, not df).
##
##  Run after sourcing data_prep_continuous.R (needs df_full and GRID_M).
##  Needs internet the first time (geodata downloads the NZ outline
##  to climate/, which the climate script has probably already done).
## ============================================================

setwd("~/Library/CloudStorage/OneDrive-MWLR/Marion Donald Rutherford Postdoc - Methods, Data and Analysis/Data_working/output/model")

source("data_prep_continuous.R")
detach("package:conflicted", unload = TRUE)

library(dplyr)
library(ggplot2)
library(sf)
library(patchwork)
library(viridis)

FAMILIES <- c("Syrphidae", "Calliphoridae", "Colletidae", "Apidae", "Halictidae")
BEE_FAM  <- c("Colletidae", "Apidae", "Halictidae")

## ── 1. Plot data: all records of the five focal families ─────
##  df_full = Apis mellifera, pre-1920 records and records without
##  coordinates or year removed. It keeps records without a month and
##  cells sampled in only one year. If it still contains other
##  families, they are dropped here.
plot_df <- df_full |>
  mutate(group = factor(if_else(family %in% BEE_FAM, "Bee", "Fly"),
                        levels = c("Bee", "Fly")))

cat("Records:", nrow(plot_df),
    "| Cells:", n_distinct(plot_df$cell),
    "| Species:", n_distinct(plot_df$species), "\n")

## ── Panel A: cells coloured by year span, from the same dataset as B/C ──

cell_span <- plot_df |>
  group_by(cell, cell_e, cell_n) |>
  summarise(year_span = max(year) - min(year), .groups = "drop")

cat("Cells:", nrow(cell_span), "\n")

cell_polys <- lapply(seq_len(nrow(cell_span)), function(i) {
  e <- cell_span$cell_e[i]
  n <- cell_span$cell_n[i]
  st_polygon(list(cbind(c(e, e + GRID_M, e + GRID_M, e,          e),
                        c(n, n,          n + GRID_M, n + GRID_M, n))))
})
cells_sf <- st_sf(cell_span, geometry = st_sfc(cell_polys, crs = 2193))

nz_land <- geodata::gadm(country = "NZ", level = 0, path = "climate/") |>
  st_as_sf() |>
  st_transform(2193) |>
  st_make_valid() |>
  st_geometry() |>
  st_union()

cells_land <- st_intersection(cells_sf, nz_land)

p_a <- ggplot() +
  geom_sf(data = st_sf(geometry = nz_land),
          fill = "grey95", colour = "grey60", linewidth = 0.2) +
  geom_sf(data = cells_land, aes(fill = year_span), colour = NA) +
  scale_fill_viridis_c(option = "mako", direction = -1, breaks = c(0, 25, 50, 75, 100),
                       name = "Year span\n(years)") +
  scale_x_continuous(breaks = seq(166, 178, 2)) +
  scale_y_continuous(breaks = seq(-46, -34, 2)) +
  coord_sf(xlim = c(166.3, 178.8), ylim = c(-47.5, -34),
           default_crs = st_crs(4326), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  theme_classic(base_size = 14) +
  theme(legend.position = "top",
        panel.border = element_rect(colour = "grey60", fill = NA))

## ── Panels B and C: records and species per year ──────────────

grp_cols <- c(Bee = "grey35", Fly = "grey67")   # Bee stacks on top of Fly
x_scale  <- scale_x_continuous(breaks = seq(1920, 2020, 20),
                               expand = expansion(mult = c(0.01, 0.02)))

p_b <- plot_df |>
  count(year, group) |>
  ggplot(aes(year, n, fill = group)) +
  geom_col(width = 1) +
  scale_fill_manual(values = grp_cols, name = NULL) +
  x_scale +
  scale_y_continuous(labels = scales::comma,
                     expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Year", y = "Number of records") +
  theme_classic(base_size = 14) +
  theme(legend.position = "top")

p_c <- plot_df |>
  group_by(year, group) |>
  summarise(n = n_distinct(species), .groups = "drop") |>
  ggplot(aes(year, n, fill = group)) +
  geom_col(width = 1) +
  scale_fill_manual(values = grp_cols, name = NULL) +
  x_scale +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Year", y = "Number of species") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none")

## ── Combine and save ────────────────────────────────────────

fig1 <- (p_a + (p_b / p_c) +
           plot_layout(widths = c(1, 1.3)) +
           plot_annotation(tag_levels = "A")) &
  theme(plot.tag = element_text(size = 22))

ggsave("fig1_sampling.png", fig1, width = 11, height = 6, dpi = 300)
ggsave("fig1_sampling.pdf", fig1, width = 11, height = 6)

## ── Numbers to check against the text ──────────────────────────
c(records = nrow(plot_df), species = n_distinct(plot_df$species),
  cells = n_distinct(plot_df$cell),
  cells_one_year = sum(cell_span$year_span == 0))
plot_df |> count(year, sort = TRUE) |> head(5)   # busiest years by records
