# Characterizing Retrospective Constraints of Survey Designs for Estimating Population-level Causal Effects

Reproducibility code for

> Tomlin, S., Andridge, R., and Lu, B. Characterizing Retrospective Constraints of Survey Designs for Estimating Population-level Causal Effects. Under review.

The paper studies post hoc estimation of the population average treatment effect (PATE) from complex survey data when the released design weights do not account for exposure-dependent selection. It proposes a weight class adjustment (WCA) to the released weights and evaluates it in a simulation study and an application to the 2021 National Youth Tobacco Survey (NYTS).

## Repository structure

Output file names start with the number of the table or figure they produce in the manuscript (Table S1, Figure S3, and so on refer to the Supplementary Material).

```
.
├── README.md
├── LICENSE
├── Simulation Study/                          Section 5
│   ├── README.md
│   ├── simulate.R                             runs the simulation, writes results/
│   ├── plot_results.R                         reads results/summary.rds, writes figures/
│   ├── results/TableS1_simulation_summary.csv     Table S1   (created when run)
│   └── figures/Figure2_simulation_results.png     Figure 2   (created when run)
└── Data Application/                          Section 6
    ├── README.md
    ├── .here                                  marks this folder as the root for here()
    ├── nyts2021_analysis.R                    main analysis
    ├── nyts2021_survey_bootstrap.R            survey bootstrap for the primary PATE
    └── output/
        ├── tables/
        │   ├── Table1_pate_estimates.csv              Table 1
        │   ├── Section6_bootstrap_pate.csv            Section 6 (text)
        │   ├── TableS2_sample_characteristics.tex     Table S2
        │   ├── TableS4_ps_coefficients.tex            Table S4
        │   └── TableS5_balance.tex                    Table S5
        ├── figures/
        │   ├── Figure3_love_plot.pdf                  Figure 3
        │   └── FigureS3_wca_prevalence.pdf            Figure S3
        └── sessionInfo_nyts.txt                       R session information
```

Code for additional diagnostics that are not reported in the manuscript is kept in the scripts but commented out, with a note explaining how to restore it.

## Map of the manuscript to the code

| Manuscript | Script | Output file |
|---|---|---|
| Figure 1 (structural scenarios) | none (conceptual diagram) | |
| Figure 2 (simulation results, δ = 1) | `Simulation Study/plot_results.R` | `figures/Figure2_simulation_results.png` |
| Table 1 (PATE under seven specifications) | `Data Application/nyts2021_analysis.R` | `output/tables/Table1_pate_estimates.csv` |
| Figure 3 (covariate balance) | `Data Application/nyts2021_analysis.R` | `output/figures/Figure3_love_plot.pdf` |
| Section 6 text (survey bootstrap SE and CI) | `Data Application/nyts2021_survey_bootstrap.R` | `output/tables/Section6_bootstrap_pate.csv` |
| Table S1 (full simulation results) | `Simulation Study/simulate.R` | `results/TableS1_simulation_summary.csv` |
| Figure S1 (NYTS cohort selection) | none (diagram) | |
| Table S2 (sample characteristics) | `Data Application/nyts2021_analysis.R` | `output/tables/TableS2_sample_characteristics.tex` |
| Figure S2 (NYTS causal model) | none (diagram) | |
| Table S3 (cell-level WCA factors) | `Data Application/nyts2021_analysis.R` | printed to the console (Section 7 of the script) |
| Figure S3 (prevalence before and after WCA) | `Data Application/nyts2021_analysis.R` | `output/figures/FigureS3_wca_prevalence.pdf` |
| Table S4 (propensity score model) | `Data Application/nyts2021_analysis.R` | `output/tables/TableS4_ps_coefficients.tex` |
| Table S5 (covariate balance) | `Data Application/nyts2021_analysis.R` | `output/tables/TableS5_balance.tex` |

## Requirements

R (version 4.3.1 was used for the paper) and the following CRAN packages.

| Component | Packages |
|---|---|
| Simulation Study | `MASS`, `survey`, `dplyr`, `tidyr`, `purrr`, `tibble`, `parallel`, `ggplot2` |
| Data Application | `here`, `haven`, `dplyr`, `tidyr`, `tibble`, `survey`, `tableone`, `cobalt`, `ggplot2`, `knitr`, `kableExtra`, `scales` |

Install them all with

```r
install.packages(c("MASS", "survey", "dplyr", "tidyr", "purrr", "tibble",
                   "ggplot2", "here", "haven", "tableone", "cobalt",
                   "knitr", "kableExtra", "scales"))
```

The exact package versions used for the data application are listed in `Data Application/output/sessionInfo_nyts.txt`.

## Data

The simulation study needs no external data.

The NYTS 2021 public-use data are freely available from the CDC Office on Smoking and Health but are not redistributed in this repository. Download the 2021 SAS dataset (`nyts2021.sas7bdat`) and format catalogue (`formats.sas7bcat`), along with the codebook and methodology report, from
<https://www.cdc.gov/tobacco/about-data/surveys/historical-nyts-data-and-documentation.html>,
and place the two data files in `Data Application/`. On the first run the analysis script caches the data as `nyts2021.rds` in the same folder.

## Reproducing the results

**Simulation study.** Set the R working directory to `Simulation Study/` and run

```r
source("simulate.R")       # 4 values of gamma x 2 values of delta x 1,000 replicates
source("plot_results.R")
```

Replicates run in parallel on all but one CPU core. Results are exactly reproducible for a given seed and number of workers (`n_cores`); see `Simulation Study/README.md`. For a quick check, set `n_sim <- 50` at the top of `simulate.R`.

**Data application.** After placing the NYTS files as described above, set the R working directory to `Data Application/` and run

```r
source("nyts2021_analysis.R")          # a few minutes
source("nyts2021_survey_bootstrap.R")  # B = 1,000 replicates; about 15 minutes
```

The bootstrap script reuses the objects created by the main analysis, and runs it first if they are not already loaded.

## Citation

If you use this code, please cite the paper above.

## License

Released under the MIT License; see `LICENSE`.

## Contact

Sean Tomlin. Questions and issues are welcome through the GitHub issue tracker.
