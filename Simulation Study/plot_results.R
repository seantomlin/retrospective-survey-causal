# plot_results.R
#
# Plots bias, MSE, and 95% CI coverage across gamma for each estimator,
# using the summary written by simulate.R. Does not rerun the simulation.
#
# Output: figures/Figure2_simulation_results.png   Figure 2 (delta = 1)

library(dplyr)
library(tidyr)
library(ggplot2)

summary_df <- readRDS("results/summary.rds")

estimator_labels <- c(
  naive   = "Naive",
  naive2  = "Naive2",
  wca_k5  = "WCA (K=5)",
  wca_k10 = "WCA (K=10)",
  wca_k15 = "WCA (K=15)",
  wca_k20 = "WCA (K=20)",
  oracle  = "Oracle"
)
estimator_shapes <- c(
  naive = 4, naive2 = 17, wca_k5 = 16, wca_k10 = 15,
  wca_k15 = 18, wca_k20 = 8, oracle = 3
)

# One panel per delta x metric, each with its own y-axis. Dashed reference
# lines mark zero bias and nominal 95% coverage.
plot_performance <- function(summary_df,
                             estimators = levels(summary_df$estimator),
                             deltas     = c(0, 1)) {
  plot_df <- summary_df %>%
    filter(estimator %in% estimators, delta %in% deltas) %>%
    select(gamma, delta, estimator, Bias = bias, MSE = mse, Coverage = coverage) %>%
    pivot_longer(c(Bias, MSE, Coverage), names_to = "metric") %>%
    mutate(
      metric    = factor(metric, levels = c("Bias", "MSE", "Coverage")),
      delta     = factor(paste0("delta == ", delta)),  # plotmath: shown as the Greek delta
      estimator = droplevels(factor(estimator, levels = names(estimator_labels)))
    )

  ref_lines <- tibble(
    metric = factor(c("Bias", "Coverage"), levels = levels(plot_df$metric)),
    yint   = c(0, 0.95)
  )

  ggplot(plot_df, aes(gamma, value, color = estimator, shape = estimator)) +
    geom_hline(data = ref_lines, aes(yintercept = yint),
               linetype = "dashed", color = "gray40") +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.2) +
    scale_color_discrete(labels = estimator_labels) +
    scale_shape_manual(values = estimator_shapes, labels = estimator_labels) +
    facet_wrap(delta ~ metric, nrow = length(deltas), scales = "free_y",
               labeller = label_parsed) +
    labs(x = expression(gamma), y = NULL, color = NULL, shape = NULL) +
    theme_bw(base_size = 13) +
    theme(
      legend.position  = "bottom",
      strip.text       = element_text(face = "bold"),
      panel.spacing    = unit(1.1, "lines")
    )
}

dir.create("figures", showWarnings = FALSE)

# [Not in the paper] All estimators, both effect sizes; uncomment to reproduce
# p_all <- plot_performance(summary_df)
# ggsave("figures/performance_all.png", p_all,
#        width = 9, height = 6, units = "in", dpi = 300)

# Figure 2: delta = 1 only
p_delta1 <- plot_performance(summary_df, deltas = 1)
ggsave("figures/Figure2_simulation_results.png", p_delta1,
       width = 6.5, height = 3.2, units = "in", dpi = 600)
