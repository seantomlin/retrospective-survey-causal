# NYTS 2021 data application

Code for the data application (Section 6) in *Characterizing Retrospective Constraints of Survey Designs for Estimating Population-level Causal Effects* (Tomlin, Andridge, and Lu, under review).
It estimates the population average treatment effect (PATE) of prior e-cigarette use
on future intention to smoke among US middle and high school students who have never
smoked a cigarette, with and without the weight class adjustment (WCA).

## Files

Output file names start with the number of the table or figure they produce in the manuscript (S numbers refer to the Supplementary Material).

```
Data Application/
├── README.md
├── .here                                   marks this folder as the root for here()
├── nyts2021_analysis.R                     main analysis
├── nyts2021_survey_bootstrap.R             survey bootstrap for the primary PATE
└── output/
    ├── tables/
    │   ├── Table1_pate_estimates.csv              Table 1: PATE under seven specifications
    │   ├── Section6_bootstrap_pate.csv            Section 6 (text): bootstrap SE and 95% CI
    │   ├── TableS2_sample_characteristics.tex     Table S2: sample characteristics by prior e-cig use
    │   ├── TableS4_ps_coefficients.tex            Table S4: propensity score model
    │   └── TableS5_balance.tex                    Table S5: covariate balance before and after IPW
    ├── figures/
    │   ├── Figure3_love_plot.pdf                  Figure 3: covariate balance
    │   └── FigureS3_wca_prevalence.pdf            Figure S3: P(Z = 1 | W) before and after WCA
    └── sessionInfo_nyts.txt                       R session information
```

The cell-level WCA factors in Table S3 are printed to the console by `nyts2021_analysis.R` (Section 7 of the script).

Code for additional diagnostics that are not in the manuscript (selection diagnostics, a forest plot of the Table 1 estimates, and the bootstrap replicates and histogram) is kept in the scripts but commented out, with a note explaining how to restore it.

## Data

The 2021 National Youth Tobacco Survey is public, but it is not redistributed here.

1. Download the 2021 SAS dataset and format catalogue from the CDC Office on Smoking and Health:
   <https://www.cdc.gov/tobacco/about-data/surveys/historical-nyts-data-and-documentation.html>
   (2021 data and documentation). The codebook and methodology report are available from the same page.
2. Place the two files in this folder (`Data Application/`), using these names:

   ```
   nyts2021.sas7bdat
   formats.sas7bcat
   ```

On the first run the script caches the data as `nyts2021.rds` in this folder.

## Running

Set the R working directory to this folder and run

```r
source("nyts2021_analysis.R")          # a few minutes
source("nyts2021_survey_bootstrap.R")  # B = 1,000 replicates; about 15 minutes
```

The bootstrap script reuses the objects created by the main analysis, and runs it first if they are not already loaded. Its random seed is fixed, so it reproduces the bootstrap results in the paper.

Required packages: `here`, `haven`, `dplyr`, `tidyr`, `tibble`, `survey`, `tableone`,
`cobalt`, `ggplot2`, `knitr`, `kableExtra`, `scales`. The versions used for the paper
(R 4.3.1) are listed in `output/sessionInfo_nyts.txt`.

## Variable definitions

| Paper | Construction | NYTS 2021 items |
|---|---|---|
| Population | Never cigarette smokers | `QN35 == 2` |
| Z (exposure) | Ever used an e-cigarette | `Q6 == "01"` |
| Y (outcome) | Any intention to smoke (0 only if both items are code `004`) | `Q49`, `Q50` |
| W (undisclosed design variable) | Survey location: on school campus vs at home/other place | `Location` |
| X (confounders) | Demographics, tobacco and e-cigarette ad exposure, secondhand smoke and vapour exposure, e-cigarette use near school, household tobacco use, home language, difficulty concentrating | `Q1`, `Q2`, `Q4A`–`Q4E`, `Q5A`–`Q5E`, `Q128`–`Q131`, `Q139`–`Q142`, `Q144`–`Q147`, `Q150`, `Q151`, `Q152B`–`Q152F`, `Q153A`–`Q153J`, `Q154`, `Q166` |

Records with the NYTS missing codes `.N` or `.Z` on any analysis variable are excluded
(complete-case analysis). Design weights are then redistributed within design strata
so the analytic sample reproduces each stratum's full-sample weight total.
