clinical_indicators <- function(d) {
    d$bp_gap_130 <- ifelse(d$hypertension_dx_or_med %in% TRUE, as.integer(d$mean_sbp >= 130 | d$mean_dbp >=
        80), NA_integer_)
    d$bp_gap_140 <- ifelse(d$hypertension_dx_or_med %in% TRUE, as.integer(d$mean_sbp >= 140 | d$mean_dbp >=
        90), NA_integer_)
    d$a1c_gap_7 <- ifelse(d$diabetes_dx %in% TRUE, as.integer(d$LBXGH >= 7), NA_integer_)
    d$a1c_gap_8 <- ifelse(d$diabetes_dx %in% TRUE, as.integer(d$LBXGH >= 8), NA_integer_)
    d$kidney_gap <- ifelse(d$ckd %in% TRUE & !is.na(d$kidney_condition_dx), as.integer(d$kidney_condition_dx %in%
        FALSE), NA_integer_)
    d$bp_applicable <- d$hypertension_dx_or_med %in% TRUE
    d$glycemia_applicable <- d$diabetes_dx %in% TRUE
    d$kidney_applicable <- d$ckd %in% TRUE
    applicable <- cbind(d$bp_applicable, d$glycemia_applicable, d$kidney_applicable)
    d$applicable_domains <- rowSums(applicable)
    primary_gaps <- cbind(d$bp_gap_130, d$a1c_gap_7, d$kidney_gap)
    primary_observed <- rowSums(!is.na(primary_gaps) & applicable)
    d$primary_all_observed <- primary_observed == d$applicable_domains
    primary_count <- rowSums(primary_gaps == 1L, na.rm = TRUE)
    primary_count[!d$primary_all_observed] <- NA_integer_
    d$primary_gap_ge1 <- as.integer(primary_count >= 1L)
    d$primary_gap_ge2 <- as.integer(primary_count >= 2L)
    strict_gaps <- cbind(d$bp_gap_140, d$a1c_gap_8, d$kidney_gap)
    strict_observed <- rowSums(!is.na(strict_gaps) & applicable)
    d$strict_all_observed <- strict_observed == d$applicable_domains
    strict_count <- rowSums(strict_gaps == 1L, na.rm = TRUE)
    strict_count[!d$strict_all_observed] <- NA_integer_
    d$strict_gap_ge1 <- as.integer(strict_count >= 1L)
    d$strict_gap_ge2 <- as.integer(strict_count >= 2L)
    d$stage <- factor(d$ckm_stage_observed, levels = 2:4)
    d$age10 <- (d$RIDAGEYR - 60)/10
    d$sex_adjust <- factor(d$RIAGENDR)
    d$race_adjust <- factor(d$race_ethnicity)
    d$fixed_bp_a1c <- d$bp_applicable & d$glycemia_applicable
    d$joint_primary <- ifelse(d$fixed_bp_a1c, as.integer(d$bp_gap_130 == 1 & d$a1c_gap_7 == 1),
        NA)
    d$joint_strict <- ifelse(d$fixed_bp_a1c, as.integer(d$bp_gap_140 == 1 & d$a1c_gap_8 == 1), NA)
    d$app_ge1 <- d$applicable_domains >= 1
    d$app_ge2 <- d$applicable_domains >= 2
    d$source3 <- ifelse(d$ckm_stage_observed != 3, NA, ifelse(isTRUE_vector(d$stage3_prevent_high_risk),
        ifelse(isTRUE_vector(d$ckd_very_high_risk_proxy), "Both", "PREVENT only"), "Kidney only"))
    d$domain_pattern <- paste0(ifelse(d$bp_applicable, "BP", ""), ifelse(d$glycemia_applicable,
        "+A1c", ""), ifelse(d$kidney_applicable, "+Kidney", ""))
    d
}
