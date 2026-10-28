################################################################################
# NYTS 2021 data application: survey bootstrap for the primary PATE
#
# Paper:    Characterizing Retrospective Constraints of Survey Designs for
#           Estimating Population-level Causal Effects. Tomlin, Andridge, and
#           Lu. Under review.
#
# Purpose:  The standard errors from nyts2021_analysis.R use Taylor
#           linearization, which accounts for the NYTS multistage design but
#           treats the WCA factor and the propensity score as fixed. This
#           script re-estimates both in each replicate of a design-based
#           survey bootstrap, for the primary specification (pooled WCA, full
#           propensity score model), and reports a percentile 95% CI
#           (reported in the text of Section 6 of the paper).
#
#           Bootstrap: Rao and Wu (1988) rescaling bootstrap of PSUs within
#           strata, via survey::as.svrepdesign(type = "subbootstrap").
#           For each replicate b = 1, ..., B:
#             1. replicate weight w^(b) carries the design-based variability
#             2. recompute the pooled WCA factor c_hat^(b)(W, Z)
#             3. refit the propensity score e_hat^(b)(X) using w^(b) * c_hat^(b)
#             4. re-estimate the PATE from the IPW-weighted regression of Y on Z
#
# Input:    Objects created by nyts2021_analysis.R (sourced below if they are
#           not already in the session): neverSmkr_design, ps_covariates,
#           pate_results.
#
# Output:   output/tables/Section6_bootstrap_pate.csv   Section 6 text: bootstrap vs
#                                                       Taylor SE and 95% CI
#           (The replicate-level estimates and a histogram of the bootstrap
#           distribution are not in the paper; that code is commented out.)
#
# Usage:    Run with the working directory set to Data Application/, after
#           placing the NYTS files there (see README.md). B = 1000 replicates
#           take about 15 minutes.
################################################################################

library(here)
library(survey)
library(dplyr)

## Run the main analysis first if its objects are not already loaded
if (!exists("pate_results") || !exists("neverSmkr_design") || !exists("ps_covariates")) {
  source(here("nyts2021_analysis.R"))
}

table_dir  <- here("output", "tables")
dir.create(table_dir,  recursive = TRUE, showWarnings = FALSE)


## =============================================================================
## 1. Replicate-weight design
## =============================================================================

B <- 1000
set.seed(20260502)

rep_design <- as.svrepdesign(
  design     = neverSmkr_design,
  type       = "subbootstrap",
  replicates = B
)

rep_wt_matrix <- weights(rep_design, type = "analysis")   # n x B
df_rep        <- rep_design$variables

cat(sprintf("Replicate-weight design: %d records x %d replicates\n",
            nrow(rep_wt_matrix), ncol(rep_wt_matrix)))


## =============================================================================
## 2. Bootstrap loop
## =============================================================================

ps_formula_full <- reformulate(ps_covariates, response = "Z")

pate_boot <- rep(NA_real_, B)
fail_flag <- logical(B)

pb <- txtProgressBar(min = 0, max = B, style = 3)

for (b in seq_len(B)) {

  ## (a) Replicate weights; units in PSUs not drawn have weight zero
  w_b  <- rep_wt_matrix[, b]
  in_b <- w_b > 0

  if (sum(in_b) < 100) {
    fail_flag[b] <- TRUE
    setTxtProgressBar(pb, b)
    next
  }

  df_b <- df_rep[in_b, ]
  df_b$w_b <- w_b[in_b]

  ## (b) Pooled WCA factor within the replicate (cells: Location_bin x Z)
  cell_b <- df_b %>%
    group_by(Location_bin, Z) %>%
    summarise(n_cell = n(), wt_cell = sum(w_b), .groups = "drop") %>%
    mutate(
      prop_sample = n_cell  / sum(n_cell),
      prop_wtd    = wt_cell / sum(wt_cell),
      c_hat       = prop_sample / prop_wtd
    )

  if (any(!is.finite(cell_b$c_hat))) {
    fail_flag[b] <- TRUE
    setTxtProgressBar(pb, b)
    next
  }

  ## (c) Apply WCA and normalise to preserve the replicate's weighted total
  df_b <- df_b %>%
    left_join(cell_b %>% select(Location_bin, Z, c_hat),
              by = c("Location_bin", "Z")) %>%
    mutate(w_b_adj = w_b * c_hat)
  df_b$w_b_adj_norm <- df_b$w_b_adj * sum(df_b$w_b) / sum(df_b$w_b_adj)

  ## (d) Propensity score on the replicate
  rep_design_b <- svydesign(id = ~psu, strata = ~stratum, weights = ~w_b_adj_norm,
                            nest = TRUE, data = df_b)

  ps_fit_b <- tryCatch(
    svyglm(ps_formula_full, design = rep_design_b,
           family = quasibinomial(link = "logit")),
    error = function(e) NULL
  )

  if (is.null(ps_fit_b)) {
    fail_flag[b] <- TRUE
    setTxtProgressBar(pb, b)
    next
  }

  ## (e) IPW, with mild trimming of extreme propensity scores
  ps_b <- pmin(pmax(fitted(ps_fit_b), 0.001), 0.999)
  df_b$w_pate_b <- df_b$w_b_adj_norm * ifelse(df_b$Z == 1, 1 / ps_b, 1 / (1 - ps_b))

  ## (f) PATE: survey-weighted regression of Y on Z
  pate_design_b <- svydesign(id = ~psu, strata = ~stratum, weights = ~w_pate_b,
                             nest = TRUE, data = df_b)

  pate_fit_b <- tryCatch(
    svyglm(Y ~ Z, design = pate_design_b, family = gaussian()),
    error = function(e) NULL
  )

  if (is.null(pate_fit_b)) {
    fail_flag[b] <- TRUE
    setTxtProgressBar(pb, b)
    next
  }

  pate_boot[b] <- unname(coef(pate_fit_b)["Z"])
  setTxtProgressBar(pb, b)
}

close(pb)


## =============================================================================
## 3. Bootstrap summary and comparison with Taylor linearization
## =============================================================================

pate_boot_ok <- pate_boot[!fail_flag]
n_success    <- length(pate_boot_ok)

primary      <- pate_results$pooled_full    # point estimate and Taylor SE/CI
boot_se      <- sd(pate_boot_ok)
boot_ci      <- unname(quantile(pate_boot_ok, probs = c(0.025, 0.975)))

boot_summary <- data.frame(
  Method      = c("Taylor linearization (weights fixed)",
                  "Survey bootstrap (WCA and PS re-estimated)"),
  Estimate    = primary$estimate,
  SE          = c(primary$se, boot_se),
  CI_lower    = c(primary$ci[1], boot_ci[1]),
  CI_upper    = c(primary$ci[2], boot_ci[2]),
  CI_type     = c("Wald", "Percentile"),
  Replicates  = c(NA, n_success),
  Failed      = c(NA, sum(fail_flag))
) %>%
  mutate(across(c(Estimate, SE, CI_lower, CI_upper), ~ round(.x, 4)))

cat("\n========== Survey bootstrap, primary PATE (pooled WCA, full PS) ==========\n")
print(boot_summary)
cat(sprintf("Bootstrap mean %.4f | median %.4f | SE ratio (bootstrap / Taylor) %.3f\n",
            mean(pate_boot_ok), median(pate_boot_ok), boot_se / primary$se))

write.csv(boot_summary, file.path(table_dir, "Section6_bootstrap_pate.csv"),
          row.names = FALSE)

## [Not in the paper] Replicate-level estimates; uncomment to save them
# write.csv(data.frame(replicate = seq_len(B), pate = pate_boot, failed = fail_flag),
#           file.path(table_dir, "bootstrap_replicates.csv"), row.names = FALSE)


## [Not in the paper] Section 4 plots the bootstrap distribution, which is
## not reported in the manuscript. Uncomment (and load ggplot2) to reproduce it.

# ## =============================================================================
# ## 4. Figure: bootstrap distribution
# ## =============================================================================

# p_boot <- ggplot(data.frame(pate = pate_boot_ok), aes(x = pate)) +
#   geom_histogram(bins = 40, fill = "steelblue4", color = "white", alpha = 0.85) +
#   geom_vline(xintercept = primary$estimate, color = "darkred",
#              linewidth = 0.8, linetype = "dashed") +
#   geom_vline(xintercept = boot_ci, color = "grey40",
#              linewidth = 0.5, linetype = "dotted") +
#   labs(
#     title    = "Survey bootstrap distribution of the PATE",
#     subtitle = sprintf("Rao-Wu rescaling bootstrap, %d replicates (pooled WCA, full PS model)",
#                        n_success),
#     x        = "PATE (risk difference)",
#     y        = "Frequency",
#     caption  = sprintf("Dashed: point estimate (%.3f). Dotted: 95%% percentile CI (%.3f, %.3f).",
#                        primary$estimate, boot_ci[1], boot_ci[2])
#   ) +
#   theme_bw(base_size = 12) +
#   theme(panel.grid.minor = element_blank())

# ggsave(here("output", "figures", "bootstrap_distribution.pdf"), p_boot,
#        width = 8, height = 5)
