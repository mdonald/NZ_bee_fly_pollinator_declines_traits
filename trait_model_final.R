library(brms)
library(dplyr)
library(tidybayes)
library(ggplot2)
library(tidyverse)

setwd("~/Library/CloudStorage/OneDrive-MWLR/Marion Donald Rutherford Postdoc - Methods, Data and Analysis/Data_working/output/data/")

sp_results_cont <- read.csv("../data/species_trends_continuous4.csv")
traits          <- read.csv("traits_df.csv")
source("../model/final_scripts/for_publication/data_prep_continuous.R")
detach("package:conflicted", unload = TRUE)

library(dplyr)


range_size <- df |>
  dplyr::group_by(species) |>
  summarise(
    n_cells_occupied = n_distinct(cell),
    .groups          = "drop"
  )


sp_results_traits_complete <- sp_results_cont |>
  left_join(traits, by = "species") |>
  dplyr::select(-itd_citation) %>% 
  # dplyr::select(
  #   species, trend_pp_decade,
  #   taxonomic_group, native_introduced,
  #   mean_ITD_mm, mean_summer_temp, mean_annual_prec,
  #   temp_breadth_sd, prec_breadth_sd,
  #   n_plant_families, fp_duration
  # ) |>
  drop_na() |>
  mutate(
    mean_ITD_mm_z        = scale(mean_ITD_mm)[, 1],
    mean_summer_temp_z   = scale(mean_summer_temp)[, 1],
    mean_annual_temp_z   = scale(mean_annual_temp)[, 1],
    mean_annual_prec_z   = scale(mean_annual_prec)[, 1],
    temp_breadth_z       = scale(temp_breadth)[, 1],
    prec_breadth_z       = scale(prec_breadth)[, 1],
    n_plant_families_z   = scale(n_plant_families)[, 1],
    fp_duration_z        = scale(fp_duration)[, 1],
    temp_breadth_sd_z    = scale(temp_breadth_sd)[, 1],
    prec_breadth_sd_z    = scale(prec_breadth_sd)[, 1],
    temp_breadth_iqr_z   = scale(temp_breadth_iqr)[, 1],  # new
    prec_breadth_iqr_z   = scale(prec_breadth_iqr)[, 1]   # new
  )

sp_joined <- sp_results_cont |>
  left_join(traits, by = "species") |>
  dplyr::select(-itd_citation)

# species dropped by drop_na()
dropped <- sp_joined |>
  anti_join(sp_results_traits_complete, by = "species")



dropped_sp <- setdiff(sp_results_cont$species, sp_results_traits_complete$species)

sp_results_cont |>
  left_join(traits, by = "species") |>
  filter(species %in% dropped_sp) |>
  select(species, where(\(x) anyNA(x)))

# which traits
sp_joined |>
  summarise(across(everything(), ~ sum(is.na(.x)))) |>
  tidyr::pivot_longer(everything(), names_to = "column", values_to = "n_na") |>
  filter(n_na > 0) |>
  arrange(desc(n_na))

# Join to traits
sp_results_traits_complete <- sp_results_traits_complete |>
  left_join(range_size, by = "species") |>
  mutate(
    n_cells_z = scale(n_cells_occupied)[, 1]
  ) 

# Check distribution
sp_results_traits_complete |>
  dplyr::select(species, n_cells_occupied, taxonomic_group) |>
  arrange(desc(n_cells_occupied)) |>
  print()


cat("Species in trait analysis:", nrow(sp_results_traits_complete), "\n")

sp_results_traits_complete |>
  dplyr::select(mean_summer_temp_z, mean_annual_prec_z,
                temp_breadth_sd_z, prec_breadth_sd_z,
                n_plant_families_z, mean_ITD_mm_z,
                fp_duration_z, n_cells_z) |>
  cor(use = "pairwise.complete.obs") |>
  round(2)

## ── 2. Propagate stage-1 uncertainty ─────────────────────────
##
##  Extract posterior SD of the trend on the probability scale
##  for each species directly from the stage-1 model draws.
##  This is used as the known measurement error in | se(se_pp).

cont_models  <- readRDS("cont_models_v4.rds")
keep_species <- names(cont_models)[!sapply(cont_models, is.null)]

se_pp_exact <- bind_rows(lapply(keep_species, function(sp) {
  m <- cont_models[[sp]]
  if (is.null(m)) return(tibble(species = sp, se_pp_exact = NA_real_))
  draws       <- as_draws_df(m)
  trend_draws <- plogis(draws[["b_Intercept"]] + draws[["b_year_c"]]) -
    plogis(draws[["b_Intercept"]])
  tibble(species = sp, se_pp_exact = sd(trend_draws))
}))

sp_results_traits_complete <- sp_results_traits_complete |>
  left_join(se_pp_exact, by = "species") |>
  mutate(se_pp = se_pp_exact)

cat("Species with se_pp:", sum(!is.na(sp_results_traits_complete$se_pp)), "\n")
cat("Missing se_pp:     ", sum( is.na(sp_results_traits_complete$se_pp)), "\n")

## ── 3. Predictor collinearity check ──────────────────────────

sp_results_traits_complete |>
  dplyr::select(mean_annual_temp_z, mean_annual_prec_z,
                temp_breadth_sd_z, prec_breadth_sd_z,
                n_plant_families_z, mean_ITD_mm_z, fp_duration_z) |>
  cor(use = "pairwise.complete.obs") |>
  round(2)

## ── 4. Model settings ────────────────────────────────────────

PRIORS <- c(
  prior(normal(0, 0.03), class = "b"),
  prior(normal(0, 0.02), class = "Intercept"),
  prior(exponential(20),  class = "sigma")
)

PRIORS_NULL <- c(
  prior(normal(0, 0.02), class = "Intercept"),
  prior(exponential(20),  class = "sigma")
)

CHAINS <- 4; ITER <- 10000; WARMUP <- 2000; CORES <- 4; SEED <- 42

fit_trait_model <- function(formula, data = sp_results_traits_complete) {
  brm(
    formula = formula,
    data    = data,
    family  = gaussian(),
    prior   = PRIORS,
    chains  = CHAINS, iter = ITER, warmup = WARMUP,
    cores   = CORES,  seed = SEED,
    control = list(adapt_delta = 0.95),
    silent  = 2, refresh = 0
  )
}

## ── 5. Models ────────────────────────────────────────────────

m_final <- fit_trait_model(
  trend_pp_decade | se(se_pp, sigma = TRUE) ~
    native_introduced +
    mean_annual_temp_z +
    n_plant_families_z +
    mean_ITD_mm_z  +
    fp_duration_z +
    n_cells_z
)



library(tidybayes)

draws_df <- m_final |>
  gather_draws(
    b_native_introducednative,
    b_mean_annual_temp_z,
    b_n_plant_families_z,
    b_mean_ITD_mm_z,
    b_fp_duration_z,
    b_n_cells_z
  ) |>
  mutate(
    .variable = dplyr::recode(.variable,
                              "b_native_introducednative" = "Origin status (native vs introduced)",
                              "b_mean_annual_temp_z"      = "Mean Annual Temp",
                              "b_n_plant_families_z"      = "Diet breadth (n plant families)",
                              "b_mean_ITD_mm_z"           = "Body size (ITD)",
                              "b_fp_duration_z"           = "Flight period duration",
                              "b_n_cells_z"               = "Range size (n grid cells)"
    )
  )



library(brms)
library(dplyr)
library(tidybayes)
library(ggplot2)
library(brms)
library(dplyr)
library(tidybayes)
library(ggplot2)

## ── 1. Draws for all six predictors, decade-scale and 100-year extrapolation ──

draws_100yr <- m_final |>
  gather_draws(
    b_mean_ITD_mm_z,
    b_n_cells_z,
    b_mean_annual_temp_z,
    b_native_introducednative,
    b_n_plant_families_z,
    b_fp_duration_z
  ) |>
  mutate(
    .variable = dplyr::recode(.variable,
                              "b_mean_ITD_mm_z"           = "Body size (ITD)",
                              "b_n_cells_z"               = "Range size (# grid cells)",
                              "b_mean_annual_temp_z"      = "Habitat mean annual temperature",
                              "b_native_introducednative" = "Origin status (native vs introduced)",
                              "b_n_plant_families_z"      = "Diet breadth (# plant families)",
                              "b_fp_duration_z"           = "Adult activity period (# days)"
    ),
    value_100yr = .value * 10   # linear scaling only — not a model estimate
  )

## ── 2. Labels: pd + 100-year estimate and CI only ──────────────────────

label_100yr <- draws_100yr |>
  group_by(.variable) |>
  summarise(
    mean_decade = mean(.value),   # kept only for ordering, not shown in label
    mean_100yr  = mean(value_100yr),
    lo_100yr    = quantile(value_100yr, 0.025),
    hi_100yr    = quantile(value_100yr, 0.975),
    pd          = max(mean(.value > 0), mean(.value < 0)),
    .groups = "drop"
  ) |>
  mutate(
    label = sprintf("%+.1f%%, (95%% CI: %.1f, %.1f), pd = %.0f%%",
                    
                    mean_100yr * 100, lo_100yr * 100, hi_100yr * 100,  pd * 100)
  )

## ── 3. Consistent ordering by |decade-scale effect| ────────────────────

effect_order_100 <- label_100yr |>
  arrange(desc(mean_decade)) |>   # most positive first = bottom, most negative = top
  pull(.variable) |>
  as.character()

draws_100yr <- draws_100yr |> mutate(.variable = factor(.variable, levels = effect_order_100))
label_100yr <- label_100yr |> mutate(.variable = factor(.variable, levels = effect_order_100))

## ── 4. Figure ────────────────────────────────────────────────────────

ggplot(draws_100yr, aes(x = value_100yr, y = .variable)) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.8, linetype = "dashed") +
  stat_slab(
    aes(fill = after_stat(x > 0)),
    slab_alpha = 0.8,
    normalize  = "groups",
    height     = 0.6
  ) +
  stat_pointinterval(
    .width = c(0.70, 0.95),
    point_interval = mean_qi,
    colour = "black",
    point_size = 2.5,
    interval_size = 1.5
  ) +
  geom_text(
    data = label_100yr,
    aes(x = Inf, y = .variable, label = label),
    inherit.aes = FALSE,
    hjust = 1.02, vjust = -1.8, size = 3.2, colour = "grey30"
  ) +
  scale_fill_manual(values = c("TRUE" = "#B896B8", "FALSE" = "#5C3A5C"), guide = "none")+
  scale_x_continuous(
    labels = scales::percent_format(accuracy = 1),
    name   = "Effect on relative occurrence over 100 years",
    expand = expansion(mult = c(0.05, 0.3))
  ) +
  labs(
    y = NULL)+ 
  theme_classic(base_size = 16) +
  theme(
    panel.grid   = element_blank(),
    axis.line.y  = element_blank(),
    axis.ticks.y = element_blank(),
    plot.caption = element_text(hjust = 0, size = 9, colour = "grey40", face = "italic")
  )

ggsave("fig_trait_100yr_extrapolated.png", width = 11, height = 5, dpi = 300)



effect_table <- draws_100yr |>
  mutate(.variable = as.character(.variable)) |>
  group_by(.variable) |>
  summarise(
    mean_decade = mean(.value),
    mean_100yr  = mean(value_100yr),
    lo_100yr    = quantile(value_100yr, 0.025),
    hi_100yr    = quantile(value_100yr, 0.975),
    pd          = max(mean(.value > 0), mean(.value < 0)),
    .groups = "drop"
  ) |>
  arrange(desc(abs(mean_decade))) |>
  transmute(
    Predictor                              = .variable,
    `Posterior mean (decade-scale)`        = signif(mean_decade, 3),
    `Posterior mean (%, per decade)`       = round(mean_decade * 100, 2),
    `100-year mean (%)`                    = round(mean_100yr * 100, 2),
    `100-year 95% CI lower (%)`            = round(lo_100yr * 100, 2),
    `100-year 95% CI upper (%)`            = round(hi_100yr * 100, 2),
    `Probability of direction (pd)`        = round(pd, 3)
  )

print(effect_table, width = Inf)
write_csv(effect_table, "effect_table.csv")
