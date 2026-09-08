national_selection_weights <- function(base) {
    base$complete_ckm <- base$ckm_primary_complete %in% TRUE
    base$complete_num <- as.integer(base$complete_ckm)
    base$sex_selection <- factor(ifelse(base$RIAGENDR == 1, "Male", "Female"))
    base$race_selection <- addNA(factor(base$race_ethnicity), ifany = TRUE)
    levels(base$race_selection)[is.na(levels(base$race_selection))] <- "Missing"
    base$education_selection <- addNA(factor(base$education_group), ifany = TRUE)
    levels(base$education_selection)[is.na(levels(base$education_selection))] <- "Missing"
    base$pir_missing <- as.integer(is.na(base$INDFMPIR))
    base$pir_selection <- ifelse(is.na(base$INDFMPIR), stats::median(base$INDFMPIR, na.rm = TRUE),
        base$INDFMPIR)
    factor3 <- function(x) factor(ifelse(is.na(x), "Missing", ifelse(x, "Yes", "No")), levels = c("No",
        "Yes", "Missing"))
    base$insurance_selection <- factor(ifelse(is.na(base$no_health_insurance), "Missing", ifelse(base$no_health_insurance,
        "Uninsured", "Insured")))
    base$hypertension_selection <- factor3(base$hypertension_dx)
    base$diabetes_selection <- factor3(base$diabetes_dx)
    base$kidney_selection <- factor3(base$kidney_condition_dx)
    base$cvd_selection <- factor3(base$stage4_clinical_cvd)
    base$smoking_selection <- factor3(base$current_smoking)
    base$weight_scaled <- base$WTPH2YR_CKM/mean(base$WTPH2YR_CKM)
    propensity_formula <- complete_num ~ splines::ns(RIDAGEYR, df = 3) + sex_selection + race_selection +
        education_selection + pir_selection + pir_missing + insurance_selection + hypertension_selection +
        diabetes_selection + kidney_selection + cvd_selection + smoking_selection
    propensity_model <- suppressWarnings(stats::glm(propensity_formula, data = base, family = quasibinomial(),
        weights = weight_scaled, na.action = stats::na.exclude))
    base$selection_probability_raw <- as.numeric(stats::predict(propensity_model, newdata = base,
        type = "response"))
    complete_prob <- base$selection_probability_raw[base$complete_ckm]
    limits <- stats::quantile(complete_prob, c(0.01, 0.99), na.rm = TRUE)
    lower <- max(0.05, limits[[1]])
    upper <- min(0.995, limits[[2]])
    base$selection_probability_truncated <- pmin(pmax(base$selection_probability_raw, lower), upper)
    base$WTPH2YR_CKM_IPW_RAW <- base$WTPH2YR_CKM/pmax(base$selection_probability_raw, 0.01)
    base$WTPH2YR_CKM_IPW_TRUNC <- base$WTPH2YR_CKM/base$selection_probability_truncated
    base
}
