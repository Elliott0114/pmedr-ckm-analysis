source("R/setup.R", encoding = "UTF-8")

args <- commandArgs(trailingOnly = TRUE)

stopifnot(length(args) == 3L)

model_objects <- readRDS(file.path(args[1], "models.rds"))

train <- readRDS(file.path(args[1], "training_data.rds"))

out <- new_output_directory(args[3])

write_table <- function(x, n) readr::write_csv(x, file.path(out, n), na = "NA")

raw <- haven::read_xpt(args[2], col_select = c("_STATE", "_LLCPWT", "_STSTR", "_PSU", "_AGE80",
    "SEXVAR", "_IMPRACE", "_BMI5CAT", "_EDUCAG", "BPHIGH6", "DIABETE4", "CHCKDNY2", "CVDINFR4",
    "CVDCRHD4", "CVDSTRK3", "CHOLCHK3", "TOLDHI3", "_SMOKER3", "_HLTHPL1", "PERSDOC3", "GENHLTH"))

raw <- haven::zap_labels(raw)

raw <- raw[!is.na(raw$`_AGE80`) & raw$`_AGE80` >= 20, ]

br <- prepare_brfss_bridge(raw)

br <- br[br$state_abbr %in% c(state.abb, "DC") & is.finite(br$brfss_weight) & br$brfss_weight >
    0, ]

unknown <- summarise(group_by(br, state_abbr, state_name), eligible_n = n(), unknown_cvd_n = sum(is.na(sr_cvd)),
    known_cvd_n = sum(!is.na(sr_cvd)), unknown_cvd_weighted_pct = 100 * sum(brfss_weight * is.na(sr_cvd))/sum(brfss_weight),
    .groups = "drop")

write_table(unknown, "unknown_cvd_by_state.csv")

br <- br[!is.na(br$sr_cvd), ]

br$row_id <- seq_len(nrow(br))

br$prevent_age_eligible <- as.integer(br$age_years >= 30 & br$age_years <= 79)

br$age_std <- cut(br$age_years, c(20, 45, 65, 80, Inf), right = FALSE, labels = c("20-44", "45-64",
    "65-79", "80+"))

ref <- tapply(br$brfss_weight, br$age_std, sum)

ref <- ref/sum(ref)

predictors <- all.vars(delete.response(terms(state_formulas()$lean_core)))

no <- br[br$sr_cvd == 0, ]

response <- fit_response_model(no, predictors)

no <- response$data

cc <- no[no$analysis_complete == 1, ]

stopifnot(nrow(br) == 413672L, nrow(cc) == 343715L, length(unique(br$state_abbr)) == 49L)

write_table(state_response_diagnostics(no), "state_response_diagnostics.csv")

write_table(data.frame(brfss_n = nrow(br), n_no_common_cvd = nrow(no), n_complete = nrow(cc), weighted_complete = sum(no$brfss_weight *
    no$analysis_complete)/sum(no$brfss_weight), n_jurisdictions = length(unique(br$state_abbr)),
    p_lower = response$lower, p_upper = response$upper), "sample_flow.csv")

write_table(data.frame(age_group = names(ref), standard_weight = as.numeric(ref)), "age_standardization_reference.csv")

rows <- list()

for (model_id in names(model_objects)) for (method in c("CC", "IPW")) {
    message("State estimates: ", model_id, " / ", method)
    rows[[length(rows) + 1L]] <- estimate_states(model_objects[[model_id]], method, model_id, br,
        cc, ref)$estimates
}

allstates <- bind_rows(rows)

stopifnot(nrow(allstates) == 1176L, all(is.finite(allstates$se)))

primary <- filter(allstates, model == "lean_core", method == "IPW")

write_table(allstates, "state_component_estimates.csv")

write_table(primary, "primary_state_estimates.csv")

ss <- summarise(group_by(primary, estimand), min = min(estimate), max = max(estimate), between_state_sd = sd(estimate),
    median_se = median(se), median_ci_width = median(high - low), .groups = "drop")

write_table(ss, "geographic_summary.csv")

wide <- tidyr::pivot_wider(select(primary, state_abbr, estimand, estimate), names_from = estimand,
    values_from = estimate)

write_table(data.frame(comparison = c("stage3 vs CVD crude", "stage3 vs CVD age-standardized", "crude vs standardized stage3"),
    spearman = c(cor(wide$stage3, wide$cvd, method = "spearman"), cor(wide$stage3_std, wide$cvd_std,
        method = "spearman"), cor(wide$stage3, wide$stage3_std, method = "spearman"))), "age_standardization_correlations.csv")

comparison <- summarise(group_by(left_join(filter(allstates, method == "IPW", estimand %in% c("stage3",
    "stage3_std")), select(primary, state_abbr, estimand, primary = estimate), by = c("state_abbr",
    "estimand")), model, estimand), spearman = cor(estimate, primary, method = "spearman"), mean_absolute_difference_pp = 100 *
    mean(abs(estimate - primary)), max_absolute_difference_pp = 100 * max(abs(estimate - primary)),
    .groups = "drop")

write_table(comparison, "model_dependence.csv")

ccipw <- mutate(tidyr::pivot_wider(select(filter(allstates, model == "lean_core"), state_abbr, estimand,
    method, estimate), names_from = method, values_from = estimate), ipw_minus_cc_pp = 100 * (IPW -
    CC))

write_table(ccipw, "normalized_cc_ipw_comparison.csv")

overlap <- transport_overlap(train, cc)

write_table(overlap$smd, "overlap_smd.csv")

write_table(overlap$metrics, "overlap_metrics.csv")

write_table(overlap$state_support, "state_support.csv")

write_table(finite_df_intervals(primary, model_objects$lean_core), "state_finite_df_sensitivity.csv")

capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
