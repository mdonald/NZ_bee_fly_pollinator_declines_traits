## ============================================================
##  Bayesian continuous-year trend model
##
##  year_c = (year - 1975) / 10   (decades, centred on 1975)
##
##  Key coefficient:
##    year_c   log-odds change per decade
##
##  The intercept = log-odds of detection at year 1975.
##  Predicted trend on the probability scale is the marginal
##  effect of +1 decade evaluated at 1975 (p(1985) - p(1975)).
## ============================================================

library(brms)
library(dplyr)
library(tidyr)
library(ggplot2)

setwd("~/Library/CloudStorage/OneDrive-MWLR/Marion Donald Rutherford Postdoc - Methods, Data and Analysis/Data_working/output/model/final_scripts/for_publication/")

source("data_prep_continuous.R")
detach("package:conflicted", unload = TRUE)


library(conflicted)
conflicts_prefer(brms::ar)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::lag)
conflicts_prefer(dplyr::summarise)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::mutate)
conflicts_prefer(dplyr::rename)
conflicts_prefer(dplyr::arrange)
conflicts_prefer(base::intersect)
conflicts_prefer(stats::chisq.test)
conflicts_prefer(stats::fisher.test)

## ── 0. brms options ──────────────────────────────────────────

CHAINS  <- 4L
ITER    <- 10000L   # 2000 → 4000: doubles post-warmup draws to 8000
WARMUP  <- 2000L   # 1000 → 2000: more adaptation for complex geometry
CORES   <- min(4L, parallel::detectCores() - 1L)
SEED    <- 42L

## ── 1. Add year_c to model data ──────────────────────────────
##
##  Centre on 1975: intercept = baseline detection at 1975.
##  Divide by 10: slope = log-odds change per decade
model_data <- sp_cell_year |>
  filter(species %in% keep_species) |>
  mutate(across(c(species, cell), factor))

model_data_cont <- model_data |>
  mutate(year_c = (year - 1975) / 10) %>% 
  mutate(obs_id = factor(seq_len(n())))




nrow(model_data_cont)                    
n_distinct(model_data_cont$cell)        
n_distinct(model_data_cont$species)     
length(shared_cells)        


## ── 2. Priors ─────────────────────────────────────────────────
##  Normal(0, 1) for year_c. ##
##  Normal(0, 2.5) for the intercept (baseline log-odds at 1975).
##  Exponential(1) for random-effect SDs.

priors_cont <- c(
  prior(normal(0, 2.5), class = "Intercept"),
  prior(normal(0, 1),   class = "b"),
  prior(exponential(1), class = "sd")
)

## ── 3. Fit function ───────────────────────────────────────────

fit_one_species_cont <- function(sp, data,
                                 priors = priors_cont,
                                 chains = CHAINS,
                                 iter   = ITER,
                                 warmup = WARMUP,
                                 cores  = CORES,
                                 seed   = SEED) {

  d <- dplyr::filter(data, species == sp)

  null_out <- function(status)
    tibble(
      species         = sp,       status          = status,
      coef_decade     = NA_real_, se_decade       = NA_real_,
      ci_lo           = NA_real_, ci_hi           = NA_real_,
      ci_lo_50        = NA_real_, ci_hi_50        = NA_real_,
      p_direction     = NA_real_, p_1975          = NA_real_,
      trend_pp_decade = NA_real_,
      n_years         = NA_integer_, year_min     = NA_integer_,
      year_max        = NA_integer_, n_cells      = NA_integer_,
      rhat_max        = NA_real_,    direction    = NA_character_
    )

  if (n_distinct(d$cell) < 2) return(null_out("too_few_cells"))
  if (n_distinct(d$year) < 3) return(null_out("too_few_years"))


  formula_brms <- bf(n_sp | trials(effort) ~ year_c + (1 | cell) + (1 | obs_id))
  
  m <- tryCatch(
    suppressMessages(brm(
      formula = formula_brms,
      data    = d,
      family  = binomial(link = "logit"),
      prior   = priors,
      chains  = chains,
      iter    = iter,
      warmup  = warmup,
      cores   = cores,
      seed    = seed,
      control = list(adapt_delta = 0.99, max_treedepth = 15),
      silent  = 2,
      refresh = 0
    )),
    error = function(e) {
      message(sprintf("  brms failed for %s: %s", sp, conditionMessage(e)))
      NULL
    }
  )

  if (is.null(m)) return(null_out("failed"))

  ## Convergence
  rhat_vals <- brms::rhat(m)
  rhat_max  <- max(rhat_vals, na.rm = TRUE)
  if (rhat_max > 1.05)
    warning(sprintf("R-hat > 1.01 for %s (max = %.3f)",
                    sp, rhat_max))

  ## Draw extraction
  draws     <- as_draws_df(m)
  col_int   <- grep("^b_Intercept", names(draws), value = TRUE)[1]
  col_trend <- grep("^b_year_c",    names(draws), value = TRUE)[1]

  int_draws   <- draws[[col_int]]
  trend_draws <- draws[[col_trend]]

  ## Posterior summary of slope (log-odds / decade)
  coef_mean <- mean(trend_draws)
  coef_sd   <- sd(trend_draws)
  ci_95     <- quantile(trend_draws, c(0.025, 0.975))
  ci_50     <- quantile(trend_draws, c(0.25,  0.75))

  ## P(direction): probability the slope is in the direction of the mean
  p_direction <- if (coef_mean >= 0) mean(trend_draws > 0) else mean(trend_draws < 0)

  ## Baseline detection at 1975 (year_c = 0, population level)
  p_1975_draws <- plogis(int_draws)
  p_1975_mean  <- mean(p_1975_draws)

  ## Marginal trend on probability scale: p(1985) - p(1975)
  ##  i.e. the change in detection probability for +1 decade at 1975
  trend_pp_draws <- plogis(int_draws + trend_draws) - p_1975_draws
  trend_pp_mean  <- mean(trend_pp_draws)
  trend_pp_ci95  <- quantile(trend_pp_draws, c(0.025, 0.975))
  trend_pp_ci50  <- quantile(trend_pp_draws, c(0.25,  0.75))

  direction <- dplyr::case_when(
    p_direction > 0.95 & coef_mean < 0 ~ "Declining",
    p_direction > 0.95 & coef_mean > 0 ~ "Increasing",
    TRUE                               ~ "No change"
  )

  list(
    model   = m,
    results = tibble(
      species         = sp,
      rhat_max        = rhat_max,
      n_years         = n_distinct(d$year),
      year_min        = min(d$year),
      year_max        = max(d$year),
      n_cells         = n_distinct(d$cell),
      ## Log-odds slope
      coef_decade     = coef_mean,
      se_decade       = coef_sd,
      ## Marginal trend on probability scale (per decade at 1975)
      trend_pp_decade = trend_pp_mean,
      ci_lo           = trend_pp_ci95[1],
      ci_hi           = trend_pp_ci95[2],
      ci_lo_50        = trend_pp_ci50[1],
      ci_hi_50        = trend_pp_ci50[2],
      ## Baseline
      p_1975          = p_1975_mean,
      ## Evidence
      p_direction     = p_direction,
      direction       = direction
    )
  )
}

## ── 4. Run all species ────────────────────────────────────────

cat("\nFitting continuous-year Bayesian models (",
    length(keep_species), "species)...\n")

cont_list <- lapply(keep_species, function(sp) {
  cat(sprintf("  Fitting: %s\n", sp))
  fit_one_species_cont(sp, data = model_data_cont, priors = priors_cont)
})

cont_models     <- setNames(lapply(cont_list, `[[`, "model"), keep_species)
sp_results_cont <- bind_rows(lapply(cont_list, function(x) {
  if (is.data.frame(x)) x else x$results
}))

saveRDS(cont_models, "cont_models_v4.rds")
readRDS("../data/cont_models_v4.rds")
## ── 5. Convergence report ─────────────────────────────────────

cat("\n── Convergence summary ──────────────────────────────────\n")
sp_results_cont |>
  dplyr::filter(!is.na(rhat_max)) |>
  summarise(
    n_ok       = sum(rhat_max < 1.01),
    n_poor     = sum(rhat_max >= 1.01)
  ) |>
  print()

poor <- dplyr::filter(sp_results_cont, !is.na(rhat_max), rhat_max > 1.01)
if (nrow(poor) > 0) {
  cat("\nSpecies with R-hat > 1.05 (rerun with iter = 4000):\n")
  print(dplyr::select(poor, species, rhat_max, status))
}

sp_results_cont |>
  dplyr::filter(!is.na(rhat_max), rhat_max > 1.01) |>
  dplyr::select(species, rhat_max, n_years, year_min, year_max, n_cells) |>
  arrange(desc(rhat_max)) |>
  print()

## ── 6. Results summary ───────────────────────────────────────

cat("\n── Trend summary (P(direction) > 0.95) ─────────────────\n")
print(count(sp_results_cont, direction))

write.csv(sp_results_cont, "species_trends_continuous4.csv", row.names = FALSE)

sp_results_cont <- read.csv("species_trends_continuous4.csv")

## ── 4b. Save posterior predicted trend lines ──────────────────
##
##  Saves population-level predictions for all species so that
##  trend line plots can be reproduced without the model objects.

year_grid_all <- seq(
  min(model_data_cont$year),
  max(model_data_cont$year),
  by = 1
)

pred_lines_all <- bind_rows(lapply(keep_species, function(sp) {
  m <- cont_models[[sp]]
  if (is.null(m)) return(NULL)
  
  ## Clip to species' own data range
  sp_years <- dplyr::filter(model_data_cont, species == sp)
  yg       <- seq(min(sp_years$year), max(sp_years$year), by = 1)
  
  newdata <- tibble(
    year_c = (yg - 1975) / 10,
    effort = 1L
  )
  
  ep <- posterior_epred(m, newdata = newdata, re_formula = NA)
  
  tibble(
    species  = sp,
    year     = yg,
    pred_med = apply(ep, 2, median),
    pred_lo  = apply(ep, 2, quantile, 0.025),
    pred_hi  = apply(ep, 2, quantile, 0.975)
  )
})) |>
  left_join(dplyr::select(sp_results_cont, species, direction,
                          p_direction),
            by = "species")

write.csv(pred_lines_all, "predicted_trend_lines4.csv", row.names = FALSE)
cat("  predicted_trend_lines.csv\n")


sp_results_cont |>
  filter(!is.na(p_1975)) |>
  arrange(p_1975) |>
  dplyr::select(species, p_1975, coef_decade, trend_pp_decade) |>
  tibble::as_tibble() |>
  print(n = Inf)
ggplot(sp_results_cont |> filter(!is.na(p_1975)),
       aes(x = p_1975)) +
  geom_histogram(bins = 20, fill = "steelblue", colour = "white") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Baseline detection probability at 1975", y = "Number of species") +
  theme_minimal()

## ── 7. Dotplot ───────────────────────────────────────────────
##
##  x-axis: marginal trend per decade on the probability scale.
##  Directly comparable to the binary model's change_pp dotplot,
##  but now the unit is percentage points per decade rather than
##  total change pre→post.
sp_results_cont |>
  dplyr::filter(!is.na(trend_pp_decade)) |>
  mutate(
    species    = reorder(species, trend_pp_decade),
    evidence   = case_when(
      ci_lo > 0 | ci_hi < 0       ~ "strong",    # 95% CI excludes zero
      ci_lo_50 > 0 | ci_hi_50 < 0 ~ "moderate",  # 50% CI excludes zero
      TRUE                         ~ "weak"        # 50% CI includes zero
    ),
    colour_var = case_when(
      evidence == "weak"   ~ "None",
      trend_pp_decade < 0  ~ "Negative",
      TRUE                 ~ "Positive"
    )
  ) |>
  #filter(species != "Leioproctus pekanui") %>% 
  ggplot(aes(x = trend_pp_decade, y = species,
             colour = colour_var)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = ci_lo,    xmax = ci_hi),
                 height = 0, linewidth = 0.5) +
  geom_errorbarh(aes(xmin = ci_lo_50, xmax = ci_hi_50),
                 height = 0, linewidth = 2) +
  geom_point(size = 2.5, colour = "white") +
  geom_point(size = 1.8) +
  scale_colour_manual(
    values = c("Negative" = "#cc0000", "Positive" = "#0055cc", "None" = "grey60"),
    guide  = "none"
  ) +
  # scale_alpha_manual(
  #   values = c("strong" = 1.0, "moderate" = 0.6, "weak" = .8),
  #   guide  = "none"
  # ) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  labs(
    x        = "Change in relative occurrence per decade \n (percentage points)",
    y        = NULL,
   # title    = "NZ pollinator occurrence: continuous year trend (Bayesian)",
    #subtitle = "Posterior mean | thick = 50% CI, thin = 95% CI | grey = 50% CI includes zero | faded = 95% CI includes zero"
  ) +
  theme_test(base_size = 18) +
  theme(axis.text.y = element_text(face = "italic"))

ggsave("continuous_trend_dotplot_v2.png",
       width  = 8,
       height = max(6, nrow(dplyr::filter(sp_results_cont,
                                          !is.na(trend_pp_decade))) * 0.28),
       dpi    = 300)




## ── Predicted trends: all species ────────────────────────────
year_grid <- seq(min(model_data_cont$year),
                 max(model_data_cont$year),
                 by = 1)

## All species with evidence classification
all_species_cont <- sp_results_cont |>
  dplyr::filter(!is.na(trend_pp_decade)) |>
  mutate(
    evidence   = case_when(
      ci_lo > 0 | ci_hi < 0       ~ "strong",
      ci_lo_50 > 0 | ci_hi_50 < 0 ~ "moderate",
      TRUE                         ~ "weak"
    ),
    colour_var = case_when(
      evidence == "weak"  ~ "None",
      trend_pp_decade < 0 ~ "Negative",
      TRUE                ~ "Positive"
    )
  ) 


all_pred_lines <- bind_rows(lapply(keep_species, function(sp) {
  m <- cont_models[[sp]]
  if (is.null(m)) return(NULL)
  
  newdata <- tibble(
    year_c = (year_grid - 1975) / 10,
    effort = 1L
  )
  
  ep <- posterior_epred(m, newdata = newdata, re_formula = NA)
  
  tibble(
    species  = sp,
    year     = year_grid,
    pred_med = apply(ep, 2, median)
  )
})) |>
  left_join(dplyr::select(all_species_cont, species, colour_var),
            by = "species") 

ggplot(all_pred_lines,
       aes(x = year, y = pred_med,
           group = species, colour = colour_var)) +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(
    values = c("Negative" = "#cc0000", "Positive" = "#0055cc", "None" = "grey60"),
    guide  = "none"
  ) +
  lims(x = c(1920, 2026))+
  scale_y_continuous(
    limits = c(0, .75),
    labels = scales::percent_format(accuracy = 1)
  )+

  labs(
    x        = "Year",
    y        = "Predicted relative occurrence (%)",
   # title    = "NZ pollinator trends",
  #  subtitle = "Population-level posterior median. Colour = direction | Faded = moderate evidence | Grey = weak evidence."
  ) +
  theme_classic(base_size = 18) +
  theme(legend.position = "top")

ggsave("all_species_trends_v4.png", width = 10, height = 6, dpi = 300)



