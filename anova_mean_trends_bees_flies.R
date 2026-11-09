# suggested by the reviewer to statistically compare the mean trends of bee and fly declines
# 6 Nov 2026
# Marion Donald

library(dplyr)


sp_results_cont <- read.csv("../data/species_trends_continuous4.csv")

bee_genera <- c("Lasioglossum", "Leioproctus", "Bombus", "Hylaeus")
fly_genera <- c("Lucilia", "Calliphora", "Melanostoma", "Melangyna",
                "Helophilus", "Eristalis", "Xenocalliphora", "Orthoprosopa")

d <- sp_results_cont |>
  mutate(
    genus = sub(" .*", "", species),
    taxon = case_when(
      genus %in% bee_genera ~ "bee",
      genus %in% fly_genera ~ "fly"
    )
  )


# Sample sizes and group means
d |>
  group_by(taxon) |>
  summarise(n = n(), mean_trend = mean(trend_pp_decade), sd = sd(trend_pp_decade))

# One-way ANOVA on the species-level trend estimates
fit <- aov(trend_pp_decade ~ taxon, data = d)
summary(fit)

# Assumption checks and a non-parametric robustness check
shapiro.test(residuals(fit))
wilcox.test(trend_pp_decade ~ taxon, data = d)


