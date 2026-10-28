################################################################################
# NYTS 2021 data application
#
# Paper:    Characterizing Retrospective Constraints of Survey Designs for
#           Estimating Population-level Causal Effects. Tomlin, Andridge, and
#           Lu. Under review.
#
# Purpose:  Estimates the population average treatment effect (PATE) of prior
#           e-cigarette use (Z) on future intention to smoke (Y) among US
#           middle and high school students who have never smoked a cigarette,
#           using the 2021 National Youth Tobacco Survey (NYTS). Compares:
#             (a) a crude survey-weighted estimate,
#             (b) survey-weighted IPW without selection correction,
#             (c) IPW combined with the weight class adjustment (WCA) for
#                 exposure-dependent selection on survey location (W),
#                 with WCA cells pooled and within design strata,
#           each under a full and a simple (pre-exposure demographics only)
#           propensity score specification.
#
# Input:    nyts2021.sas7bdat and formats.sas7bcat, placed in this folder
#           (download instructions in README.md; the data are not
#           redistributed in this repository)
#
# Output (file names give the manuscript table or figure number):
#           output/tables/Table1_pate_estimates.csv          Table 1   PATE under seven specifications
#           output/figures/Figure3_love_plot.pdf             Figure 3  Covariate balance
#           output/tables/TableS2_sample_characteristics.tex Table S2  Sample characteristics by Z
#           output/figures/FigureS3_wca_prevalence.pdf       Figure S3 P(Z = 1 | W) before/after WCA
#           output/tables/TableS4_ps_coefficients.tex        Table S4  Propensity score model
#           output/tables/TableS5_balance.tex                Table S5  Covariate balance before/after IPW
#           output/sessionInfo_nyts.txt                                R session information
#           The cell-level WCA factors in Table S3 are printed to the console
#           (Section 7).
#
# Usage:    Set the working directory to this folder (the .here file marks it
#           as the project root for here::here()). Runtime: a few minutes.
#
# Sections: 1 Setup | 2 Data import | 3 Analytic sample | 4 Weight
#           redistribution | 5 Labels and Table S2 | 6 Selection diagnostics
#           (not in the paper; commented out) | 7 Pooled WCA (Table S3,
#           Figure S3) | 8 Stratum WCA | 9 PATE estimation (Table 1) |
#           10 Balance (Figure 3, Table S5) | 11 Forest plot (not in the
#           paper; commented out) | 12 Table S4 | 13 Session info
################################################################################


## =============================================================================
## 1. Setup
## =============================================================================

library(here)
library(haven)
library(dplyr)
library(tidyr)
library(tibble)
library(survey)
library(tableone)
library(cobalt)
library(ggplot2)
library(knitr)
library(kableExtra)

options(survey.lonely.psu = "adjust")
options(survey.adjust.domain.lonely = TRUE)

data_dir   <- here()
table_dir  <- here("output", "tables")
figure_dir <- here("output", "figures")
dir.create(table_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

# Figures are saved as vector PDFs. cairo_pdf handles the non-ASCII
# characters (e.g. "ĉ") used in some axis labels.
save_figure <- function(plot, file, width = 10, height = 4.5) {
  ggsave(file.path(figure_dir, file), plot = plot,
         width = width, height = height, device = cairo_pdf)
}


## =============================================================================
## 2. Data import
##    The SAS file is read once (applying the CDC format catalogue) and cached
##    as .rds for faster reruns.
## =============================================================================

rds_file <- file.path(data_dir, "nyts2021.rds")

if (file.exists(rds_file)) {
  nyts2021_raw <- readRDS(rds_file)
} else {
  nyts2021_raw <- read_sas(
    file.path(data_dir, "nyts2021.sas7bdat"),
    catalog_file = file.path(data_dir, "formats.sas7bcat")
  )
  saveRDS(nyts2021_raw, rds_file)
}


## =============================================================================
## 3. Analytic sample and variable construction
##    NYTS uses the special missing codes ".N" and ".Z"; records with these
##    codes on any analysis variable are excluded.
## =============================================================================

nyts_design <- svydesign(
  id      = ~psu,
  strata  = ~stratum,
  weights = ~finwgt,
  nest    = TRUE,
  data    = nyts2021_raw
)

## --- Never cigarette smokers (QN35 == 2) with non-missing Location ----------
neverSmkr_design <- subset(
  nyts_design,
  QN35 == 2 &
    !is.na(Location) &
    !(Location %in% c(".N", ".Z", ""))
)

## --- W: survey location, binary ----------------------------------------------
neverSmkr_design <- update(
  neverSmkr_design,
  Location = factor(
    Location,
    levels = c("At home (virtual learning)",
               "In a school building/classroom",
               "Some other place")
  ),
  Location_bin = factor(
    case_when(
      Location == "In a school building/classroom" ~ "On School Campus",
      Location %in% c("At home (virtual learning)", "Some other place") ~ "At Home/Other Place",
      TRUE ~ NA_character_
    ),
    levels = c("On School Campus", "At Home/Other Place")
  )
)
neverSmkr_design <- subset(neverSmkr_design, !is.na(Location_bin))

## --- Y: future intention to smoke (Q49, Q50) ----------------------------------
## Y = 0 if the response to both Q49 and Q50 is code "004", 1 otherwise.
neverSmkr_design <- subset(
  neverSmkr_design,
  !(Q49 %in% c(".N", ".Z")) &
    !(Q50 %in% c(".N", ".Z"))
)
neverSmkr_design <- update(
  neverSmkr_design,
  Y = 1 - (Q49 == "004") * (Q50 == "004")
)

## --- Z: ever used an e-cigarette (Q6 == "01") ---------------------------------
neverSmkr_design <- subset(neverSmkr_design, !(Q6 %in% c(".N", ".Z")))
neverSmkr_design <- update(neverSmkr_design, Z = as.numeric(Q6 == "01"))

## --- Age (Q1) and sex (Q2) ------------------------------------------------------
neverSmkr_design <- subset(neverSmkr_design, Q1 != ".N")
neverSmkr_design <- update(
  neverSmkr_design,
  Age = factor(
    case_when(
      Q1 %in% c("01", "02", "03") ~ "11 or younger",
      Q1 == "04" ~ "12",
      Q1 == "05" ~ "13",
      Q1 == "06" ~ "14",
      Q1 == "07" ~ "15",
      Q1 == "08" ~ "16",
      Q1 == "09" ~ "17",
      Q1 %in% c("10", "11") ~ "18 or older",
      TRUE ~ NA_character_
    ),
    levels = c("11 or younger", "12", "13", "14", "15", "16", "17", "18 or older")
  )
)

neverSmkr_design <- subset(neverSmkr_design, Q2 != ".N")
neverSmkr_design <- update(
  neverSmkr_design,
  Sex = factor(case_when(Q2 == "01" ~ "Male",
                         Q2 == "02" ~ "Female"))
)

## --- Hispanic ethnicity (Q4A-Q4E; Q4A == "1" is not Hispanic) -------------------
neverSmkr_design <- subset(
  neverSmkr_design,
  !(Q4A %in% c(".N", ".Z")) &
    !(Q4B %in% c(".N", ".Z")) &
    !(Q4C %in% c(".N", ".Z")) &
    !(Q4D %in% c(".N", ".Z")) &
    !(Q4E %in% c(".N", ".Z"))
)
neverSmkr_design <- update(
  neverSmkr_design,
  Hispanic = factor(ifelse(Q4A == "1", "No", "Yes"))
)

## --- Race (Q5A-Q5E) ---------------------------------------------------------------
## More than one race selected -> "Other"; no item answered -> "Not answered".
neverSmkr_design <- subset(
  neverSmkr_design,
  !(Q5A %in% c(".Z")) &
    !(Q5B %in% c(".Z")) &
    !(Q5C %in% c(".Z")) &
    !(Q5D %in% c(".Z")) &
    !(Q5E %in% c(".Z"))
)
neverSmkr_design <- update(
  neverSmkr_design,
  race_count      = (Q5A == "1") + (Q5B == "1") + (Q5C == "1") + (Q5D == "1") + (Q5E == "1"),
  raceWhite       = Q5E == "1",
  raceBlack       = Q5C == "1",
  raceOther       = (Q5A == "1" | Q5B == "1" | Q5D == "1"),
  raceNotAnswered = Q5A == ".N" & Q5B == ".N" & Q5C == ".N" & Q5D == ".N" & Q5E == ".N",
  Race = factor(
    ifelse(raceNotAnswered, "Not answered",
           ifelse(race_count > 1, "Other",
                  ifelse(raceWhite, "White",
                         ifelse(raceBlack, "Black",
                                ifelse(raceOther, "Other", NA_character_))))),
    levels = c("White", "Black", "Other", "Not answered")
  )
)

## --- E-cigarette use near school (any of Q152B-Q152F) --------------------------------
neverSmkr_design <- subset(
  neverSmkr_design,
  !(Q152B %in% c(".N", ".Z")) &
    !(Q152C %in% c(".N", ".Z")) &
    !(Q152D %in% c(".N", ".Z")) &
    !(Q152E %in% c(".N", ".Z")) &
    !(Q152F %in% c(".N", ".Z"))
)
neverSmkr_design <- update(
  neverSmkr_design,
  eCig_at_School = factor(ifelse(
    Q152B == "1" | Q152C == "1" | Q152D == "1" | Q152E == "1" | Q152F == "1",
    "Yes", "No"
  ))
)

## --- Household tobacco use (any of Q153A-Q153J) -------------------------------------
neverSmkr_design <- subset(
  neverSmkr_design,
  !(Q153A %in% c(".N", ".Z")) &
    !(Q153B %in% c(".N", ".Z")) &
    !(Q153C %in% c(".N", ".Z")) &
    !(Q153D %in% c(".N", ".Z")) &
    !(Q153E %in% c(".N", ".Z")) &
    !(Q153F %in% c(".N", ".Z")) &
    !(Q153G %in% c(".N", ".Z")) &
    !(Q153H %in% c(".N", ".Z")) &
    !(Q153I %in% c(".N", ".Z")) &
    !(Q153J %in% c(".N", ".Z"))
)
neverSmkr_design <- update(
  neverSmkr_design,
  smoke_in_house = factor(ifelse(
    Q153A == "1" | Q153B == "1" | Q153C == "1" | Q153D == "1" | Q153E == "1" |
      Q153F == "1" | Q153G == "1" | Q153H == "1" | Q153I == "1" | Q153J == "1",
    "Yes", "No"
  ))
)

## --- Single-item covariates: exclude missing codes ----------------------------------
cols_single_2021 <- c(
  "Q139", "Q140", "Q141", "Q142",   # tobacco ads
  "Q128", "Q129", "Q130", "Q131",   # e-cigarette ads
  "Q144", "Q145",                   # secondhand smoke: home, vehicle (past 7 days)
  "Q146", "Q147",                   # public smoke exposure: indoor, outdoor
  "Q150", "Q151",                   # public vapour exposure: indoor, outdoor
  "Q154", "Q166",                   # non-English at home, difficulty concentrating
  "Location"
)
neverSmkr_design <- subset(
  neverSmkr_design,
  Reduce(`&`, lapply(cols_single_2021, function(v) !get(v) %in% c(".N", ".Z")))
)

## --- Derived and recoded covariates --------------------------------------------------
## Response codes 1-8 map to 0-7 days.
neverSmkr_design <- update(
  neverSmkr_design,
  Days_Q144 = as.numeric(Q144) - 1,
  Days_Q145 = as.numeric(Q145) - 1
)
neverSmkr_design <- update(
  neverSmkr_design,
  Q154 = factor(case_when(Q154 == "01" ~ "Yes", Q154 == "02" ~ "No")),
  Q166 = factor(case_when(Q166 == "01" ~ "Yes", Q166 == "02" ~ "No"))
)
neverSmkr_design <- update(
  neverSmkr_design,
  Race = relevel(Race, ref = "White"),
  Sex  = relevel(Sex,  ref = "Male")
)

analytic_df <- model.frame(neverSmkr_design)

cat(sprintf("\nFull NYTS 2021 sample:  %d\nAnalytic sample:        %d\n",
            nrow(nyts2021_raw), nrow(analytic_df)))


## =============================================================================
## 4. Redistribution of design weights to the analytic sample
##    Within each design stratum, weights of retained records are scaled so
##    they sum to the stratum's full-sample weight total. The survey design is
##    then rebuilt on the analytic sample with the redistributed weights.
## =============================================================================

nyts2021_raw <- nyts2021_raw %>%
  mutate(in_analytic = StudentLoginID %in% analytic_df$StudentLoginID)

stratum_redist <- nyts2021_raw %>%
  group_by(stratum) %>%
  summarise(
    wt_total_full = sum(finwgt, na.rm = TRUE),
    wt_total_kept = sum(finwgt[in_analytic], na.rm = TRUE),
    n_full        = n(),
    n_kept        = sum(in_analytic),
    .groups       = "drop"
  ) %>%
  mutate(redist_factor = if_else(wt_total_kept > 0,
                                 wt_total_full / wt_total_kept,
                                 NA_real_))

cat("\n--- Weight redistribution diagnostics ---\n")
cat(sprintf("Number of strata:               %d\n", nrow(stratum_redist)))
cat(sprintf("Strata with no retained records: %d\n", sum(is.na(stratum_redist$redist_factor))))
cat(sprintf("Redistribution factor: mean %.4f, median %.4f, range [%.3f, %.3f]\n",
            mean(stratum_redist$redist_factor, na.rm = TRUE),
            median(stratum_redist$redist_factor, na.rm = TRUE),
            min(stratum_redist$redist_factor, na.rm = TRUE),
            max(stratum_redist$redist_factor, na.rm = TRUE)))

## finwgt is replaced by the redistributed weight; the original is kept as finwgt_orig
analytic_df <- analytic_df %>%
  mutate(finwgt_orig = finwgt) %>%
  left_join(stratum_redist %>% select(stratum, redist_factor), by = "stratum") %>%
  mutate(finwgt = finwgt_orig * redist_factor) %>%
  select(-redist_factor)

cat(sprintf("Weighted N: full NYTS %.0f | analytic, original %.0f | analytic, redistributed %.0f\n",
            sum(nyts2021_raw$finwgt, na.rm = TRUE),
            sum(analytic_df$finwgt_orig),
            sum(analytic_df$finwgt)))

neverSmkr_design <- svydesign(
  id      = ~psu,
  strata  = ~stratum,
  weights = ~finwgt,
  nest    = TRUE,
  data    = analytic_df
)


## =============================================================================
## 5. Labels and Table S2 (sample characteristics)
## =============================================================================

## --- Covariates -----------------------------------------------------------------
## Location_bin (W) enters through the WCA, not the propensity score model.
ps_covariates <- c(
  "Age", "Sex", "Hispanic", "Race",
  "Q139", "Q140", "Q141", "Q142",
  "Q128", "Q129", "Q130", "Q131",
  "Days_Q144", "Days_Q145",
  "Q146", "Q147", "Q150", "Q151",
  "eCig_at_School", "smoke_in_house",
  "Q154", "Q166"
)
table1_vars <- c(ps_covariates, "Location_bin")

## --- Display labels for variables and response codes --------------------------------
var_labels <- c(
  Age            = "Age",
  Sex            = "Sex",
  Hispanic       = "Hispanic ethnicity",
  Race           = "Race",
  Q139           = "Tobacco ads exposure: Internet",
  Q140           = "Tobacco ads exposure: Newspapers/Magazines",
  Q141           = "Tobacco ads exposure: Retail stores",
  Q142           = "Tobacco ads exposure: TV/Streaming/Movies",
  Q128           = "E-cig ads exposure: Internet",
  Q129           = "E-cig ads exposure: Newspapers/Magazines",
  Q130           = "E-cig ads exposure: Retail stores",
  Q131           = "E-cig ads exposure: TV/Streaming/Movies",
  Days_Q144      = "Days someone smoked tobacco in home (past 7 days)",
  Days_Q145      = "Days rode in vehicle with someone smoking tobacco (past 7 days)",
  Q146           = "Days exposed to tobacco smoke indoors (past 30 days)",
  Q147           = "Days exposed to tobacco smoke outdoors (past 30 days)",
  Q150           = "Days exposed to e-cig vapor indoors (past 30 days)",
  Q151           = "Days exposed to e-cig vapor outdoors (past 30 days)",
  eCig_at_School = "E-cigarette use near school",
  smoke_in_house = "Household smoking behavior",
  Q154           = "Language other than English at home",
  Q166           = "Difficulty concentrating",
  Location_bin   = "Location during survey"
)

make_freq_scale <- function(label_001) {
  c("001" = label_001, "002" = "Never", "003" = "Rarely",
    "004" = "Sometimes", "005" = "Most of the time", "006" = "Always")
}
internet_levels <- make_freq_scale("I do not use the Internet")
print_levels    <- make_freq_scale("I do not read newspapers or magazines")
retail_levels   <- make_freq_scale("Never go to a convenience store/etc.")
tv_levels       <- make_freq_scale("Do not watch TV/streaming or movies")
days_exposure_levels <- c(
  "001" = "0 days", "002" = "1 or 2 days", "003" = "3 to 5 days",
  "004" = "6 to 9 days", "005" = "10 to 19 days", "006" = "20 to 29 days",
  "007" = "All 30 days"
)
codebook <- list(
  Q139 = internet_levels, Q140 = print_levels, Q141 = retail_levels, Q142 = tv_levels,
  Q128 = internet_levels, Q129 = print_levels, Q130 = retail_levels, Q131 = tv_levels,
  Q146 = days_exposure_levels, Q147 = days_exposure_levels,
  Q150 = days_exposure_levels, Q151 = days_exposure_levels
)

## --- Helpers to translate model terms into display labels ----------------------------
## parse_term() splits terms such as "Q139002" (svyglm), "Q139_002" (cobalt)
## or "eCig_at_SchoolYes" into a base variable and a level, matching the
## longest variable name first.
term_bases <- table1_vars[order(-nchar(table1_vars))]

parse_term <- function(term) {
  base <- term_bases[startsWith(term, term_bases)][1]
  if (is.na(base)) return(c(base = term, level = ""))
  c(base = base, level = sub("^_", "", substring(term, nchar(base) + 1)))
}
label_var <- function(base) {
  if (base %in% names(var_labels)) unname(var_labels[base]) else base
}
label_level <- function(base, level) {
  if (base %in% names(codebook) && level %in% names(codebook[[base]])) {
    unname(codebook[[base]][level])
  } else {
    level
  }
}

## --- Table 1: unweighted sample characteristics by Z ----------------------------------
table1_2021 <- CreateTableOne(
  vars   = table1_vars,
  strata = "Z",
  data   = analytic_df,
  test   = FALSE
)

tab1_mat <- print(
  table1_2021,
  showAllLevels = TRUE,
  smd           = FALSE,
  quote         = FALSE,
  noSpaces      = TRUE,
  printToggle   = FALSE
)

## Replace variable names and response codes with display labels
current_base <- ""
for (i in seq_len(nrow(tab1_mat))) {
  rn <- rownames(tab1_mat)[i]
  if (rn != "") {
    indicator    <- regmatches(rn, regexpr(" \\(.*\\)$", rn))
    current_base <- sub(" \\(.*\\)$", "", rn)
    rownames(tab1_mat)[i] <- paste0(label_var(current_base),
                                    if (length(indicator)) indicator else "")
  }
  tab1_mat[i, "level"] <- label_level(current_base, tab1_mat[i, "level"])
}

tab1_latex <- kable(
  tab1_mat,
  format    = "latex",
  booktabs  = TRUE,
  longtable = TRUE,
  escape    = TRUE,
  caption   = "Sample characteristics by prior e-cigarette use ($Z$) among never-smokers in the 2021 NYTS. Cell values are $n$ (\\%) for categorical variables and mean (SD) for continuous variables.",
  label     = "tbl1-nyts2021"
) %>%
  kable_styling(latex_options = c("striped", "hold_position", "repeat_header"),
                font_size = 8) %>%
  column_spec(1, width = "4.5cm") %>%
  column_spec(2, width = "5cm")

writeLines(tab1_latex, file.path(table_dir, "TableS2_sample_characteristics.tex"))


## [Not in the paper] Section 6 produces additional selection diagnostics
## that are not reported in the manuscript. It is commented out; uncomment
## to reproduce them. Nothing later in the script depends on it.

# ## =============================================================================
# ## 6. Diagnostics: Z-dependent selection within the survey design
# ##    If selection depends on Z beyond the design, the distribution of
# ##    selection probabilities (1 / finwgt) differs by Z within design cells.
# ## =============================================================================

# z_labels <- c("Z = 0 (no prior e-cig use)", "Z = 1 (prior e-cig use)")

# diag_df <- neverSmkr_design$variables %>%
#   mutate(
#     sel_prob   = 1 / finwgt,
#     log_finwgt = log(finwgt),
#     Z_label    = factor(Z, levels = c(0, 1), labels = z_labels)
#   )

# cat("\nStratum x PSU cells:", n_distinct(paste(diag_df$stratum, diag_df$psu)),
#     "\nStrata:", n_distinct(diag_df$stratum), "\n")

# print(
#   diag_df %>%
#     group_by(Z_label) %>%
#     summarise(n           = n(),
#               mean_prob   = mean(sel_prob),
#               median_prob = median(sel_prob),
#               q25         = quantile(sel_prob, 0.25),
#               q75         = quantile(sel_prob, 0.75),
#               .groups     = "drop")
# )

# ## Mean log-weight difference (Z = 1 minus Z = 0) within stratum x PSU cells
# ## that contain both exposure groups
# psu_strata_wide <- diag_df %>%
#   group_by(stratum, psu, Z) %>%
#   summarise(mean_log_wt = mean(log_finwgt), .groups = "drop") %>%
#   pivot_wider(id_cols = c(stratum, psu), names_from = Z,
#               values_from = mean_log_wt, names_prefix = "mean_log_wt_Z") %>%
#   filter(!is.na(mean_log_wt_Z0), !is.na(mean_log_wt_Z1)) %>%
#   mutate(diff_log_wt = mean_log_wt_Z1 - mean_log_wt_Z0)

# ## Figure: density of selection probabilities by Z
# p_sel_density <- ggplot(diag_df, aes(x = sel_prob, fill = Z_label)) +
#   geom_density(alpha = 0.5) +
#   scale_x_log10() +
#   scale_fill_manual(values = c("steelblue4", "darkorange3")) +
#   labs(
#     title    = "Selection probability distribution by exposure status",
#     subtitle = "Overall — log scale; differences indicate Z-dependent selection",
#     x        = "Selection probability  (1 / finwgt)  [log scale]",
#     y        = "Density",
#     fill     = NULL
#   ) +
#   theme_bw() +
#   theme(legend.position = "bottom", panel.grid.minor = element_blank())
# save_figure(p_sel_density, "nyts_selection_density.pdf")

# ## Figure: within-cell log-weight differences
# mean_diff <- mean(psu_strata_wide$diff_log_wt)
# p_sel_withincell <- ggplot(psu_strata_wide, aes(x = diff_log_wt)) +
#   geom_histogram(bins = 40, fill = "steelblue4", color = "white", alpha = 0.85) +
#   geom_vline(xintercept = 0, color = "black", linewidth = 0.6) +
#   geom_vline(xintercept = mean_diff, color = "darkorange3",
#              linewidth = 0.7, linetype = "dashed") +
#   annotate("text", x = mean_diff, y = Inf,
#            label = sprintf("mean = %.3f", mean_diff),
#            color = "darkorange3", hjust = -0.1, vjust = 1.5, size = 3.5) +
#   labs(
#     title    = "Within-cell log-weight difference (Z = 1 − Z = 0)",
#     subtitle = "Each bar = one PSU × stratum cell; nonzero mean = Z-dependent selection",
#     x        = "Mean log(finwgt | Z=1) − Mean log(finwgt | Z=0)",
#     y        = "Number of PSU × stratum cells"
#   ) +
#   theme_bw() +
#   theme(panel.grid.minor = element_blank())
# save_figure(p_sel_withincell, "nyts_selection_withincell.pdf")

# ## Figure: log-weights by Z in 12 randomly chosen strata (both Z groups, n >= 50)
# set.seed(2021)
# sample_strata <- diag_df %>%
#   group_by(stratum) %>%
#   filter(n_distinct(Z) == 2, n() >= 50) %>%
#   ungroup() %>%
#   pull(stratum) %>%
#   unique()
# sample_strata <- sample(sample_strata, min(12, length(sample_strata)))

# p_sel_strata <- diag_df %>%
#   filter(stratum %in% sample_strata) %>%
#   ggplot(aes(x = factor(stratum), y = log_finwgt, fill = Z_label)) +
#   geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
#   scale_fill_manual(values = c("steelblue4", "darkorange3")) +
#   labs(
#     title    = "log(finwgt) by Z within sampled strata",
#     subtitle = "Systematic vertical separation indicates Z-dependent selection within design strata",
#     x        = "Stratum (sampled)",
#     y        = "log(finwgt)",
#     fill     = NULL
#   ) +
#   theme_bw() +
#   theme(legend.position = "bottom",
#         axis.text.x = element_text(angle = 45, hjust = 1),
#         panel.grid.minor = element_blank())
# save_figure(p_sel_strata, "nyts_selection_strata.pdf")


## =============================================================================
## 7. Weight class adjustment (WCA), pooled cells (Table S3 values, Figure S3)
##    Cells are Location_bin x Z. The correction factor is
##      ĉ = (sample share of cell) / (finwgt-implied share of cell),
##    and adjusted weights are normalised to preserve the weighted total.
## =============================================================================

cell_adjust <- analytic_df %>%
  group_by(Location_bin, Z) %>%
  summarise(n_cell = n(), wt_cell = sum(finwgt), .groups = "drop") %>%
  mutate(
    n_total     = sum(n_cell),
    wt_total    = sum(wt_cell),
    prop_sample = n_cell  / n_total,
    prop_wtd    = wt_cell / wt_total,
    c_hat       = prop_sample / prop_wtd
  )
print(cell_adjust)   # values reported in Table S3

analytic_df <- analytic_df %>%
  left_join(cell_adjust %>% select(Location_bin, Z, c_hat),
            by = c("Location_bin", "Z")) %>%
  mutate(finwgt_adj = finwgt * c_hat)

scale_factor <- sum(analytic_df$finwgt) / sum(analytic_df$finwgt_adj)
analytic_df  <- analytic_df %>% mutate(finwgt_adj_norm = finwgt_adj * scale_factor)

cat(sprintf("\nWeighted N: original %.0f | adjusted %.0f | normalised %.0f\n",
            sum(analytic_df$finwgt), sum(analytic_df$finwgt_adj),
            sum(analytic_df$finwgt_adj_norm)))

## Sample restrictions were applied in Section 3; confirm before rebuilding
stopifnot(all(analytic_df$QN35 == 2), !anyNA(analytic_df$Location))

neverSmkr_design_adj <- svydesign(
  id      = ~psu,
  strata  = ~stratum,
  weights = ~finwgt_adj_norm,
  nest    = TRUE,
  data    = analytic_df
)

## --- Figure: P(Z = 1 | Location) before and after WCA --------------------------------
spec_before <- "Before WCA: finwgt"
spec_after  <- "After WCA: finwgt × ĉ(Location, Z)"

prev_combined <- bind_rows(
  svyby(~Z, ~Location_bin, design = neverSmkr_design, FUN = svymean, na.rm = TRUE) %>%
    as.data.frame() %>% mutate(Specification = spec_before),
  svyby(~Z, ~Location_bin, design = neverSmkr_design_adj, FUN = svymean, na.rm = TRUE) %>%
    as.data.frame() %>% mutate(Specification = spec_after)
) %>%
  rename(prop = Z) %>%
  mutate(
    ci_lo = prop - 1.96 * se,
    ci_hi = prop + 1.96 * se,
    Specification = factor(Specification, levels = c(spec_before, spec_after))
  )
print(prev_combined)

## Prevalence ratio, At Home/Other Place relative to On School Campus
print(
  prev_combined %>%
    group_by(Specification) %>%
    summarise(p_campus     = prop[Location_bin == "On School Campus"],
              p_home_other = prop[Location_bin == "At Home/Other Place"],
              ratio_home_other_to_campus = p_home_other / p_campus,
              .groups = "drop")
)

p_prev <- ggplot(prev_combined, aes(x = Specification, y = prop, fill = Location_bin)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.75),
           width = 0.65, alpha = 0.9) +
  geom_text(aes(label = sprintf("%.1f%%", prop * 100), y = prop + 0.012),
            position = position_dodge(width = 0.75),
            size = 3.6, fontface = "bold") +
  scale_fill_manual(values = c("steelblue4", "steelblue1"), name = "Survey setting") +
  scale_y_continuous(labels = scales::percent,
                     limits = c(0, max(prev_combined$ci_hi) * 1.15),
                     expand = c(0, 0)) +
  labs(title = "Conditional prevalence of prior e-cig use by survey setting",
       x = NULL, y = "P(Z = 1 | Location)") +
  theme_bw(base_size = 12) +
  theme(
    plot.title         = element_text(face = "bold", size = 13),
    axis.text.x        = element_text(face = "bold", size = 11),
    axis.title.y       = element_text(face = "bold"),
    legend.position    = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank()
  )
save_figure(p_prev, "FigureS3_wca_prevalence.pdf", width = 8, height = 5)


## =============================================================================
## 8. Sensitivity analysis: WCA within design strata
##    Cells are stratum x Location_bin x Z. Cells with fewer than 30 records
##    (or a non-finite factor) fall back to the pooled ĉ from Section 7.
## =============================================================================

pooled_adjust <- cell_adjust %>% select(Location_bin, Z, c_hat_pool = c_hat)

cell_adjust_stratum <- analytic_df %>%
  group_by(stratum, Location_bin, Z) %>%
  summarise(n_cell_stratum = n(), wt_cell_stratum = sum(finwgt), .groups = "drop") %>%
  group_by(stratum) %>%
  mutate(
    prop_sample_stratum = n_cell_stratum  / sum(n_cell_stratum),
    prop_wtd_stratum    = wt_cell_stratum / sum(wt_cell_stratum),
    c_hat_stratum_raw   = prop_sample_stratum / prop_wtd_stratum
  ) %>%
  ungroup() %>%
  mutate(cell_stable = n_cell_stratum >= 30 & is.finite(c_hat_stratum_raw)) %>%
  left_join(pooled_adjust, by = c("Location_bin", "Z")) %>%
  mutate(c_hat_stratum = if_else(cell_stable, c_hat_stratum_raw, c_hat_pool))

cat("\n--- Stratum-cell diagnostics ---\n")
cat(sprintf("Stratum cells: %d | stable (n >= 30): %d (%.1f%%)\n",
            nrow(cell_adjust_stratum), sum(cell_adjust_stratum$cell_stable),
            100 * mean(cell_adjust_stratum$cell_stable)))
cat(sprintf("Range of ĉ in stable cells: [%.3f, %.3f]\n",
            min(cell_adjust_stratum$c_hat_stratum_raw[cell_adjust_stratum$cell_stable], na.rm = TRUE),
            max(cell_adjust_stratum$c_hat_stratum_raw[cell_adjust_stratum$cell_stable], na.rm = TRUE)))

analytic_df <- analytic_df %>%
  left_join(cell_adjust_stratum %>% select(stratum, Location_bin, Z, c_hat_stratum),
            by = c("stratum", "Location_bin", "Z")) %>%
  mutate(finwgt_adj_stratum = finwgt * c_hat_stratum)

scale_factor_stratum <- sum(analytic_df$finwgt) / sum(analytic_df$finwgt_adj_stratum)
analytic_df <- analytic_df %>%
  mutate(finwgt_adj_norm_stratum = finwgt_adj_stratum * scale_factor_stratum)

neverSmkr_design_adj_stratum <- svydesign(
  id      = ~psu,
  strata  = ~stratum,
  weights = ~finwgt_adj_norm_stratum,
  nest    = TRUE,
  data    = analytic_df
)

cat("\n--- P(Z = 1 | Location) under stratum WCA ---\n")
print(svyby(~Z, ~Location_bin, design = neverSmkr_design_adj_stratum,
            FUN = svymean, na.rm = TRUE))

cat("\n--- Weight distributions ---\nOriginal finwgt:\n")
print(summary(analytic_df$finwgt))
cat("Pooled WCA weight:\n");  print(summary(analytic_df$finwgt_adj_norm))
cat("Stratum WCA weight:\n"); print(summary(analytic_df$finwgt_adj_norm_stratum))


## =============================================================================
## 9. PATE estimation (Table 1)
##    For each design: fit a survey-weighted propensity score model, form ATE
##    inverse probability weights, multiply by the design (or WCA) weight, and
##    regress Y on Z with an identity link. The coefficient on Z is the PATE
##    on the risk difference scale.
## =============================================================================

ps_formula_full   <- reformulate(ps_covariates, response = "Z")
## Simple specification: pre-exposure demographics and home language only
## (sensitivity check for covariates possibly measured after exposure)
ps_formula_simple <- Z ~ Age + Sex + Hispanic + Race + Q154

estimate_pate <- function(design, ps_formula) {
  ps_model <- svyglm(ps_formula, design = design,
                     family = quasibinomial(link = "logit"))
  dat <- design$variables
  ps  <- fitted(ps_model)
  stopifnot(length(ps) == nrow(dat))  # no records dropped for missing covariates
  
  dat$ipw_ate <- ifelse(dat$Z == 1, 1 / ps, 1 / (1 - ps))
  dat$w_pate  <- weights(design) * dat$ipw_ate
  
  pate_design <- svydesign(id = ~psu, strata = ~stratum, weights = ~w_pate,
                           nest = TRUE, data = dat)
  fit <- svyglm(Y ~ Z, design = pate_design, family = gaussian())
  
  list(ps_model = ps_model,
       data     = dat,
       estimate = unname(coef(fit)["Z"]),
       se       = unname(SE(fit)["Z"]),
       ci       = unname(confint(fit)["Z", ]),
       p_value  = summary(fit)$coefficients["Z", "Pr(>|t|)"])
}

## --- Crude: design weights only, no propensity score ----------------------------------
crude_fit <- svyglm(Y ~ Z, design = neverSmkr_design, family = gaussian())
crude <- list(estimate = unname(coef(crude_fit)["Z"]),
              se       = unname(SE(crude_fit)["Z"]),
              ci       = unname(confint(crude_fit)["Z", ]),
              p_value  = summary(crude_fit)$coefficients["Z", "Pr(>|t|)"])

## --- Propensity score specifications ---------------------------------------------------
pate_results <- list(
  crude          = crude,
  noWCA_simple   = estimate_pate(neverSmkr_design,             ps_formula_simple),
  noWCA_full     = estimate_pate(neverSmkr_design,             ps_formula_full),
  pooled_simple  = estimate_pate(neverSmkr_design_adj,         ps_formula_simple),
  pooled_full    = estimate_pate(neverSmkr_design_adj,         ps_formula_full),
  stratum_simple = estimate_pate(neverSmkr_design_adj_stratum, ps_formula_simple),
  stratum_full   = estimate_pate(neverSmkr_design_adj_stratum, ps_formula_full)
)

spec_labels <- c(
  crude          = "Marginal (crude): finwgt only, no PS",
  noWCA_simple   = "Simple PS-adjusted: finwgt × IPW(X)",
  noWCA_full     = "PS-adjusted: finwgt × IPW(X)",
  pooled_simple  = "Simple PS-adjusted: finwgt × ĉ(W,Z) × IPW(X)",
  pooled_full    = "PS-adjusted: finwgt × ĉ(W,Z) × IPW(X)",
  stratum_simple = "Simple PS-adjusted: finwgt × ĉ_stratum(W,Z) × IPW(X)",
  stratum_full   = "PS-adjusted: finwgt × ĉ_stratum(W,Z) × IPW(X)"
)

seven_way <- data.frame(
  Specification = unname(spec_labels[names(pate_results)]),
  Estimate      = sapply(pate_results, `[[`, "estimate"),
  SE            = sapply(pate_results, `[[`, "se"),
  CI_lower      = sapply(pate_results, function(r) r$ci[1]),
  CI_upper      = sapply(pate_results, function(r) r$ci[2]),
  p_value       = sapply(pate_results, `[[`, "p_value"),
  row.names     = NULL
) %>%
  mutate(across(c(Estimate, SE, CI_lower, CI_upper), ~ round(.x, 4)),
         p_value = signif(p_value, 3))

cat("\n========== PATE (risk difference) under seven specifications ==========\n")
print(seven_way)

write.csv(seven_way, file.path(table_dir, "Table1_pate_estimates.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

## Change from adding the pooled WCA (full PS specification)
delta_est <- pate_results$pooled_full$estimate - pate_results$noWCA_full$estimate
cat(sprintf("\nChange in PATE, pooled WCA vs no WCA: %+.4f (%+.2f%%)\n",
            delta_est, 100 * delta_est / pate_results$noWCA_full$estimate))


## =============================================================================
## 10. Covariate balance, pooled WCA and full PS specification (Figure 3, Table S5)
## =============================================================================

pooled_full <- pate_results$pooled_full

## Note: with un = TRUE, cobalt's "unadjusted" column is the unweighted sample.
bal_pop_ipw <- bal.tab(
  x        = pooled_full$data[, ps_covariates],
  treat    = pooled_full$data$Z,
  weights  = pooled_full$data$w_pate,
  estimand = "ATE",
  un       = TRUE
)
print(bal_pop_ipw)

## Display labels for balance terms
bal_terms  <- rownames(bal_pop_ipw$Balance)
bal_parsed <- t(sapply(bal_terms, parse_term))
bal_labels <- data.frame(
  old      = bal_terms,
  variable = sapply(bal_parsed[, "base"], label_var),
  level    = mapply(label_level, bal_parsed[, "base"], bal_parsed[, "level"]),
  stringsAsFactors = FALSE,
  row.names = NULL
)
bal_labels$new <- ifelse(bal_labels$level == "", bal_labels$variable,
                         paste0(bal_labels$variable, ": ", bal_labels$level))

## --- Figure: love plot ----------------------------------------------------------------
## v_map labels each balance term with its variable label only (no response
## level), e.g. "Q139_002" -> "Tobacco ads exposure: Internet".
v_map <- data.frame(
  old = bal_labels$old,
  new = bal_labels$variable,
  stringsAsFactors = FALSE
)
rownames(v_map) <- v_map$old

love_plot <- love.plot(
  bal_pop_ipw,
  title         = "Population-level Covariate Balance",
  abs           = FALSE,
  var.names     = v_map,
  stars         = "std",
  drop.distance = TRUE,
  wrap          = 1,
  thresholds    = c(m = .1),
  limits        = c(-0.35, 0.35),
  alpha         = 0.5,
  size          = 3,
  shapes        = c("circle", "triangle"),
  colors        = c("grey40", "steelblue4"),
  sample.names  = c("Before IPW", "After IPW")
) +
  xlab("Standardized Mean Difference") +
  theme_minimal() +
  theme(panel.spacing    = unit(2, "lines"),
        strip.text       = element_text(size = 12, face = "bold"),
        panel.grid.minor = element_blank(),
        legend.position  = "bottom")
save_figure(love_plot, "Figure3_love_plot.pdf")

## --- Table: balance before and after IPW -------------------------------------------------
bal_table <- data.frame(
  Variable     = bal_labels$variable,
  Level        = bal_labels$level,
  `Before IPW` = sprintf("%.3f", bal_pop_ipw$Balance$Diff.Un),
  `After IPW`  = sprintf("%.3f", bal_pop_ipw$Balance$Diff.Adj),
  check.names  = FALSE
) %>%
  group_by(Variable) %>%
  mutate(Variable = ifelse(row_number() == 1, Variable, "")) %>%
  ungroup()

bal_latex <- kable(
  bal_table,
  format    = "latex",
  booktabs  = TRUE,
  longtable = TRUE,
  escape    = TRUE,
  align     = c("l", "l", "r", "r"),
  caption   = "Population-level covariate balance before and after IPW. Values are standardized mean differences. Conventional balance threshold is $|\\text{SMD}| < 0.1$.",
  label     = "balance-pop-ipw"
) %>%
  kable_styling(latex_options = c("striped", "hold_position", "repeat_header"),
                font_size = 9) %>%
  column_spec(1, width = "5cm") %>%
  column_spec(2, width = "5.5cm")

writeLines(bal_latex, file.path(table_dir, "TableS5_balance.tex"))


## [Not in the paper] Section 11 produces a forest plot of the Table 1
## estimates that is not reported in the manuscript. It is commented out;
## uncomment to reproduce it. Nothing later in the script depends on it.

# ## =============================================================================
# ## 11. Figure: PATE forest plot (percentage points)
# ## =============================================================================

# forest_df <- data.frame(
#   Specification = c("Marginal (crude)",
#                     "Simple PS-adjusted (no WCA)",
#                     "PS-adjusted (no WCA)",
#                     "Simple PS-adjusted + pooled WCA",
#                     "PS-adjusted + pooled WCA",
#                     "Simple PS-adjusted + stratum WCA",
#                     "PS-adjusted + stratum WCA"),
#   Estimate = seven_way$Estimate * 100,
#   CI_lower = seven_way$CI_lower * 100,
#   CI_upper = seven_way$CI_upper * 100
# )
# forest_df$AdjType <- factor(
#   case_when(
#     grepl("Marginal", forest_df$Specification)            ~ "Marginal",
#     grepl("^Simple PS-adjusted", forest_df$Specification) ~ "Simple PS-adjusted",
#     TRUE                                                  ~ "PS-adjusted"
#   ),
#   levels = c("Marginal", "Simple PS-adjusted", "PS-adjusted")
# )
# forest_df$Specification <- factor(forest_df$Specification,
#                                   levels = rev(forest_df$Specification))

# forest_plot <- ggplot(forest_df, aes(x = Estimate, y = Specification, color = AdjType)) +
#   geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.4) +
#   geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.22, linewidth = 0.7) +
#   geom_point(size = 4, shape = 18) +
#   geom_text(aes(label = sprintf("%.2f (%.2f, %.2f)", Estimate, CI_lower, CI_upper),
#                 x = 33),
#             hjust = 0, size = 3.4, color = "black", family = "sans") +
#   scale_color_manual(values = c("Marginal"           = "grey45",
#                                 "Simple PS-adjusted" = "#a3c4d9",
#                                 "PS-adjusted"        = "steelblue4"),
#                      name = NULL) +
#   scale_x_continuous(breaks = seq(0, 30, 5), limits = c(-2, 47), expand = c(0, 0)) +
#   labs(
#     title    = "PATE of Prior E-cigarette Use on Future Intention to Smoke",
#     subtitle = "NYTS 2021 never-smokers; risk difference scale; specification comparison",
#     x        = "Risk difference (percentage points)",
#     y        = NULL
#   ) +
#   theme_bw(base_size = 12) +
#   theme(
#     plot.title         = element_text(face = "bold", size = 13),
#     plot.subtitle      = element_text(face = "italic", color = "grey30", size = 11),
#     axis.title.x       = element_text(face = "bold"),
#     axis.text.y        = element_text(size = 11),
#     legend.position    = "bottom",
#     panel.grid.major.y = element_blank(),
#     panel.grid.minor   = element_blank()
#   )
# save_figure(forest_plot, "nyts_pate_forest.pdf")


## =============================================================================
## 12. Table S4: propensity score model
##     (pooled WCA design, full specification; odds ratios with Wald CIs on
##     the design degrees of freedom)
## =============================================================================

ps_svy_model <- pooled_full$ps_model
ps_ci        <- confint(ps_svy_model, ddf = degf(neverSmkr_design_adj))
ps_parsed    <- t(sapply(names(coef(ps_svy_model)), parse_term))

ps_table <- data.frame(
  Variable = sapply(ps_parsed[, "base"], label_var),
  Level    = mapply(label_level, ps_parsed[, "base"], ps_parsed[, "level"]),
  OR_CI    = sprintf("%.2f (%.2f, %.2f)",
                     exp(coef(ps_svy_model)), exp(ps_ci[, 1]), exp(ps_ci[, 2])),
  row.names = NULL
) %>%
  group_by(Variable) %>%
  mutate(Variable = ifelse(row_number() == 1, Variable, "")) %>%
  ungroup()

ps_latex <- kable(
  ps_table,
  format    = "latex",
  booktabs  = TRUE,
  longtable = TRUE,
  escape    = TRUE,
  align     = c("l", "l", "r"),
  col.names = c("Variable", "Level", "OR (95% CI)"),
  caption   = "Propensity score model coefficients for prior e-cigarette use ($Z$) among never-smokers in the 2021 NYTS. Estimates are odds ratios with 95\\% confidence intervals, computed from survey-weighted logistic regression using the design degrees of freedom.",
  label     = "ps-coefs"
) %>%
  kable_styling(latex_options = c("striped", "hold_position", "repeat_header"),
                font_size = 9) %>%
  column_spec(1, width = "5cm") %>%
  column_spec(2, width = "5.5cm")

writeLines(ps_latex, file.path(table_dir, "TableS4_ps_coefficients.tex"))


## =============================================================================
## 13. Session info
## =============================================================================

writeLines(capture.output(sessionInfo()),
           here("output", "sessionInfo_nyts.txt"))

