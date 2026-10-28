# Weight-class adjustment for treatment-dependent sample selection

Code for the simulation study (Section 5) in *Characterizing Retrospective Constraints of Survey Designs for Estimating Population-level Causal Effects* (Tomlin, Andridge, and Lu, under review).

## Question
When a sample's selection probabilities depend on treatment status, standard inverse-probability weighting that models selection only on covariates gives biased estimates of the population average treatment effect (PATE). This study compares how well several weighting strategies recover the PATE as the dependence of selection on treatment (gamma) increases:

| Estimator | Selection model |
|---|---|
| Naive | S ~ 1 (ignores selection) |
| Naive2 | S ~ X (ignores dependence on treatment) |
| WCA (K) | S ~ X, plus a weight-class adjustment using the treated share within K classes of P(S = 1 \| X) |
| Oracle | S ~ X + Z (correctly specified) |

WCA is run with K = 5 and K = 15 weight classes. Each scenario is a combination of gamma (0, 0.5, 1, 1.5) and true effect size delta (0 or 1 SD of the control outcome), with 1,000 replicates each. Estimators are compared on bias, variance, MSE, and 95% CI coverage.

## Files

Output file names start with the number of the table or figure they produce in the manuscript (Table S1 is in the Supplementary Material).

```
Simulation Study/
├── README.md
├── simulate.R             generates the populations, computes all estimators, writes results/
├── plot_results.R         reads results/summary.rds, writes figures/
├── results/               created when simulate.R runs
│   ├── sim_results.rds                  one row per replicate and scenario
│   ├── summary.rds                      bias, variance, MSE, and coverage by scenario
│   └── TableS1_simulation_summary.csv   Table S1: full simulation results
└── figures/               created when plot_results.R runs
    └── Figure2_simulation_results.png   Figure 2: bias, MSE, and coverage for delta = 1
```

A version of the figure with both effect sizes is not in the manuscript; its code is commented out in `plot_results.R`.

## Running

Set the R working directory to this folder and run

```r
source("simulate.R")
source("plot_results.R")
```

Replicates run in parallel across all but one CPU core, using a PSOCK cluster that works on Windows, macOS, and Linux. For a quick check, set `n_sim <- 50` at the top of `simulate.R`.

Each worker gets its own random number seed, so results are reproducible for a given seed and number of workers. To reproduce results exactly on a different machine, set `n_cores` to the same value used originally.

Requires R and the packages MASS, survey, dplyr, tidyr, purrr, tibble, parallel, and ggplot2.
