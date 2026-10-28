# simulate.R
#
# Paper: Characterizing Retrospective Constraints of Survey Designs for
#        Estimating Population-level Causal Effects. Tomlin, Andridge, and Lu.
#        Under review. Section 5 (simulation study).
#
# Simulation study: estimating a population average treatment effect (PATE)
# from a sample whose selection depends on treatment.
#
# Each replicate draws a finite population with covariates X, treatment Z,
# and potential outcomes Y(0), Y(1). Units enter the sample with a
# probability that depends on X and, through gamma, on Z. All estimators are
# inverse-probability-weighted differences in means that combine a selection
# probability P(S = 1 | .) with a propensity score P(Z = 1 | X). They differ
# only in how the selection probability is modeled:
#
#   naive    S ~ 1        ignores selection
#   naive2   S ~ X        ignores the dependence of selection on Z
#   wca_kK   S ~ X, plus a weight-class adjustment (WCA) using the observed
#            treatment share within K classes of P(S = 1 | X)
#   oracle   S ~ X + Z    correctly specified selection model
#
# Outputs (written to results/):
#   sim_results.rds                 one row per replicate x scenario
#   summary.rds                     bias, variance, MSE, and coverage by scenario
#                                   (read by plot_results.R for Figure 2)
#   TableS1_simulation_summary.csv  Table S1: the same summary as a table
#
# Replicates run in parallel. Results are reproducible for a given seed and
# number of workers (n_cores); changing n_cores changes the random streams.
# Set n_sim to 50 for a quick check.

library(MASS)      # mvrnorm()
library(survey)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(parallel)

# Settings ----------------------------------------------------------------

seed       <- 2025
n_cores    <- max(1, detectCores() - 1)
n_pop      <- 50000                  # finite population size
n_sim      <- 1000                   # replicates per scenario
gamma_list <- c(0, 0.5, 1.0, 1.5)    # log-odds effect of Z on selection
delta_list <- c(0, 1)                # PATE, in SD units of Y(0)
wca_k_list <- c(5, 15)               # number of weight classes for WCA

# Treatment model: about 40% treated
treat_intercept <- -0.42
treat_slopes    <- c(x1 = 0.2, x2 = 0.3)

# Selection model: about 10,000 sampled units when gamma = 0
select_intercept <- -1.6
select_slopes    <- c(x1 = 0.5, x2 = -0.5, x3 = 0.5)

estimators <- c("naive", "naive2", paste0("wca_k", wca_k_list), "oracle")

# Data generation ---------------------------------------------------------

generate_population <- function(n_pop, delta, gamma) {
  x <- MASS::mvrnorm(n_pop, mu = rep(0, 3), Sigma = diag(3))
  colnames(x) <- paste0("x", 1:3)
  pop <- as.data.frame(x)

  # Potential outcomes. The effect is delta SDs of Y(0) on average and
  # varies with x2, so the realized PATE differs slightly across replicates.
  pop$y_0 <- with(pop, x1 + 0.5 * x2 + rnorm(n_pop))
  pop$y_1 <- with(pop, x1 + x2 + rnorm(n_pop)) + delta * sd(pop$y_0)

  pop$z <- rbinom(n_pop, 1, plogis(
    treat_intercept + treat_slopes["x1"] * pop$x1 + treat_slopes["x2"] * pop$x2
  ))

  lin_s <- select_intercept +
    drop(as.matrix(pop[, names(select_slopes)]) %*% select_slopes) +
    gamma * pop$z
  pop$s <- rbinom(n_pop, 1, plogis(lin_s))

  pop
}

# Weights -----------------------------------------------------------------

# Inverse of P(selected) x P(observed treatment arm)
ipw_weights <- function(p_select, ps, z) {
  ifelse(z == 1, 1 / (p_select * ps), 1 / (p_select * (1 - ps)))
}

# WCA: group sampled units into K classes of P(S = 1 | X), then also divide
# by the observed share of treated (or control) units within each class.
wca_weights <- function(p_select, ps, z, K) {
  weight_class <- dplyr::ntile(p_select, K)
  p_treated    <- ave(z, weight_class)   # treated share within each class
  ifelse(z == 1,
         1 / (p_select * p_treated * ps),
         1 / (p_select * (1 - p_treated) * (1 - ps)))
}

# Estimator ---------------------------------------------------------------

# Weighted (Hajek) difference in means between arms, with a 95% CI.
# The SE treats the estimated weights as fixed and the two arm means as
# independent.
ipw_estimate <- function(obs, w, truth) {
  obs$w  <- w
  design <- svydesign(ids = ~1, weights = ~w, data = obs)
  m1     <- svymean(~y_obs, subset(design, z == 1))
  m0     <- svymean(~y_obs, subset(design, z == 0))

  est <- coef(m1)[[1]] - coef(m0)[[1]]
  se  <- sqrt(vcov(m1)[1, 1] + vcov(m0)[1, 1])

  c(est = est, cover = abs(est - truth) <= 1.96 * se)
}

# One replicate -----------------------------------------------------------

run_replicate <- function(n_pop, delta, gamma) {
  pop   <- generate_population(n_pop, delta, gamma)
  truth <- mean(pop$y_1 - pop$y_0)   # realized finite-population PATE

  # Propensity score and selection models are fit on the full population
  pop$ps       <- fitted(glm(z ~ x1 + x2 + x3, data = pop, family = binomial()))
  pop$p_naive  <- mean(pop$s)
  pop$p_naive2 <- fitted(glm(s ~ x1 + x2 + x3, data = pop, family = binomial()))
  pop$p_oracle <- fitted(glm(s ~ x1 + x2 + x3 + z, data = pop, family = binomial()))

  obs <- pop %>%
    filter(s == 1) %>%
    mutate(y_obs = if_else(z == 1, y_1, y_0))

  out <- list(
    naive  = ipw_estimate(obs, ipw_weights(obs$p_naive,  obs$ps, obs$z), truth),
    naive2 = ipw_estimate(obs, ipw_weights(obs$p_naive2, obs$ps, obs$z), truth)
  )
  for (K in wca_k_list) {
    out[[paste0("wca_k", K)]] <-
      ipw_estimate(obs, wca_weights(obs$p_naive2, obs$ps, obs$z, K), truth)
  }
  out$oracle <- ipw_estimate(obs, ipw_weights(obs$p_oracle, obs$ps, obs$z), truth)

  c(unlist(out), pate_actual = truth, n_sample = nrow(obs))
}

# Parallel setup ----------------------------------------------------------

# PSOCK clusters work on Windows, macOS, and Linux. Each worker gets its own
# L'Ecuyer-CMRG random number stream from the seed.
cl <- makeCluster(n_cores)
invisible(clusterEvalQ(cl, {
  library(survey)
  library(dplyr)
}))
clusterExport(cl, c(
  "n_pop", "treat_intercept", "treat_slopes", "select_intercept",
  "select_slopes", "wca_k_list", "generate_population", "ipw_weights",
  "wca_weights", "ipw_estimate", "run_replicate"
))
clusterSetRNGStream(cl, seed)

# Run all scenarios -------------------------------------------------------

scenarios <- expand_grid(gamma = gamma_list, delta = delta_list)

sim_results <- pmap_dfr(scenarios, function(gamma, delta) {
  reps <- parSapply(cl, seq_len(n_sim),
                    function(i, gamma, delta) run_replicate(n_pop, delta, gamma),
                    gamma = gamma, delta = delta)
  as_tibble(t(reps)) %>%
    mutate(gamma = gamma, delta = delta, rep = row_number(), .before = 1)
})

stopCluster(cl)

# Summarize ---------------------------------------------------------------

summary_df <- sim_results %>%
  pivot_longer(
    cols      = -c(gamma, delta, rep, pate_actual, n_sample),
    names_to  = c("estimator", ".value"),
    names_sep = "\\."
  ) %>%
  group_by(gamma, delta, estimator) %>%
  summarise(
    bias      = mean(est - pate_actual),
    variance  = var(est),
    mse       = mean((est - pate_actual)^2),
    coverage  = mean(cover),
    mcse_cov  = sqrt(coverage * (1 - coverage) / n()),  # Monte Carlo SE
    .groups   = "drop"
  ) %>%
  mutate(estimator = factor(estimator, levels = estimators)) %>%
  arrange(gamma, delta, estimator)

dir.create("results", showWarnings = FALSE)
saveRDS(sim_results, "results/sim_results.rds")
saveRDS(summary_df,  "results/summary.rds")

# Table S1: bias, variance, MSE, and coverage for every scenario and estimator
write.csv(summary_df %>% select(-mcse_cov),
          "results/TableS1_simulation_summary.csv", row.names = FALSE)

print(summary_df, n = Inf)
