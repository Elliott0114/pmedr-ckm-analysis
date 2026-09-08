# Cardiovascular-kidney-metabolic syndrome in United States adults

R code for the final national analyses and exploratory geographic projections in
*Recognition and Control Gaps in Advanced Cardiovascular-Kidney-Metabolic Syndrome
Among United States Adults*.

The analysis uses the National Health and Nutrition Examination Survey (NHANES),
August 2021–August 2023, and the 2023 Behavioral Risk Factor Surveillance System
(BRFSS). The scripts read the public survey files, derive clinical indicators and
survey weights, fit the national and state models, and write numerical results.

## Data and environment

Download the survey files listed in [data/README.md](data/README.md). File hashes
are recorded in [data_sources.csv](data_sources.csv). Data are stored locally and
are excluded from Git.

Use R 4.5.3 and restore the package versions in `renv.lock`:

```r
install.packages("renv")  # Once, if renv is unavailable.
renv::restore(prompt = FALSE)
```

The lock file contains the analysis packages and their dependencies. No package
updates are needed to run this version. Run all commands below from the repository
root. An output directory must be empty or absent; use a new directory for a rerun.

## Run the analyses

```sh
Rscript scripts/01_prepare_national.R data/nhanes outputs/prepared
Rscript scripts/02_national_analysis.R outputs/prepared/national_data.rds data/nhanes outputs/national
Rscript scripts/03_fit_state_models.R outputs/prepared/national_data.rds outputs/models
Rscript scripts/04_validate_state_models.R outputs/models outputs/validation
Rscript scripts/05_estimate_states.R outputs/models data/brfss/LLCP2023.XPT outputs/states
```

| Script | Analyses and outputs |
| --- | --- |
| `01_prepare_national.R` | Survey joins, clinical definitions, selection weights, PREVENT risk and final staging; sample flow and analysis data |
| `02_national_analysis.R` | National proportions, fixed hypertension–diabetes subgroup, kidney definitions, fasting and risk-input sensitivities, adjusted stage comparisons and socioeconomic models |
| `03_fit_state_models.R` | Primary and sensitivity models on the same training sample; coefficients, model dimensions and exclusions |
| `04_validate_state_models.R` | Twenty repetitions of five-fold primary-sampling-unit cross-validation, 200 bootstrap repetitions, and calibration with group counts |
| `05_estimate_states.R` | Reported cardiovascular disease and modeled stage-3 contributions, age standardization, missing-predictor weighting, model sensitivity and survey overlap |

## Interpretation

National estimates use the phlebotomy weight. Fasting comparisons use the same
eligible fasting sample and `WTSAF2YR` for both definitions. PREVENT calculations
are restricted to ages 30–79; boundary adjustments affect risk inputs only.
Uncomputable risk within this age range does not assign a low stage by default.

Descriptive proportions use Korn–Graubard beta intervals with the domain degrees
of freedom and retain reliability flags. Boundary proportions have no estimated
interval. Stage comparisons and prevalence ratios retain their specified model
interval methods, including the design-degree-of-freedom approximation for the
socioeconomic models. Model size and covariance-rank diagnostics accompany these
estimates; the full socioeconomic model does not support an omnibus joint test.

The geographic outputs are exploratory projections. `stage3` and `stage3_std`
denote modeled stage-3 contributions among all adults with known common
cardiovascular-disease status, not directly measured state prevalence. The state
intervals combine survey sampling variance and coefficient variance. They do not
cover all transport, measurement or model-specification uncertainty. The finite
degrees-of-freedom calculation is an interval sensitivity approximation.

This repository covers the final clinical and geographic analyses. Manuscript
formatting, artwork assembly, historical alternatives and correspondence are
outside its scope. Additional external state-context source data are documented
in the manuscript's source-data workbook.
