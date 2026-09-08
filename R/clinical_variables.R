derive_ckm_variables <- function(data) {
    data$mean_sbp <- mean_available(data, c("BPXOSY1", "BPXOSY2", "BPXOSY3", "BPXSY1", "BPXSY2",
        "BPXSY3"))
    data$mean_dbp <- mean_available(data, c("BPXODI1", "BPXODI2", "BPXODI3", "BPXDI1", "BPXDI2",
        "BPXDI3"))
    data$egfr <- egfr_ckdepi_2021(col_or_na(data, "LBXSCR"), col_or_na(data, "RIDAGEYR"), col_or_na(data,
        "RIAGENDR"))
    data$uacr <- uacr_mg_g(data)
    data$sex_label <- factor(ifelse(data$RIAGENDR == 1, "Male", ifelse(data$RIAGENDR == 2, "Female",
        NA_character_)), levels = c("Male", "Female"))
    data$prevent_sex <- ifelse(data$RIAGENDR == 1, "male", ifelse(data$RIAGENDR == 2, "female",
        NA_character_))
    data$prevent_age <- data$RIDAGEYR
    data$race_ethnicity <- factor(dplyr::case_when(data$RIDRETH3 == 1 ~ "Mexican American", data$RIDRETH3 ==
        2 ~ "Other Hispanic", data$RIDRETH3 == 3 ~ "Non-Hispanic White", data$RIDRETH3 == 4 ~ "Non-Hispanic Black",
        data$RIDRETH3 == 6 ~ "Non-Hispanic Asian", data$RIDRETH3 == 7 ~ "Other/Multi-racial", TRUE ~
            NA_character_))
    data$age_group <- cut(data$RIDAGEYR, breaks = c(20, 40, 60, 80, Inf), right = FALSE, labels = c("20-39",
        "40-59", "60-79", "80+"))
    data$education_group <- factor(dplyr::case_when(data$DMDEDUC2 %in% c(1, 2) ~ "<High school",
        data$DMDEDUC2 == 3 ~ "High school/GED", data$DMDEDUC2 == 4 ~ "Some college/AA", data$DMDEDUC2 ==
            5 ~ "College graduate+", TRUE ~ NA_character_), levels = c("<High school", "High school/GED",
        "Some college/AA", "College graduate+"))
    data$pir_group <- cut(data$INDFMPIR, breaks = c(-Inf, 1.3, 3.5, Inf), labels = c("<1.30", "1.30-3.49",
        ">=3.50"), right = FALSE)
    data$hypertension_dx <- yes_no(col_or_na(data, "BPQ020"))
    data$antihypertensive_use <- dplyr::case_when(col_or_na(data, "BPQ150") == 1 ~ TRUE, col_or_na(data,
        "BPQ150") == 2 ~ FALSE, col_or_na(data, "BPQ020") == 2 ~ FALSE, TRUE ~ NA)
    data$measured_hypertension <- data$mean_sbp >= 130 | data$mean_dbp >= 80
    data$hypertension <- data$measured_hypertension | isTRUE_vector(data$hypertension_dx) | isTRUE_vector(data$antihypertensive_use)
    data$diabetes_dx <- col_or_na(data, "DIQ010") == 1
    data$borderline_diabetes_dx <- col_or_na(data, "DIQ010") == 3
    data$diabetes_medication <- dplyr::case_when(col_or_na(data, "DIQ050") == 1 | col_or_na(data,
        "DIQ070") == 1 ~ TRUE, col_or_na(data, "DIQ010") == 2 ~ FALSE, col_or_na(data, "DIQ050") ==
        2 & col_or_na(data, "DIQ070") == 2 ~ FALSE, TRUE ~ NA)
    data$lab_diabetes_a1c <- data$LBXGH >= 6.5
    data$diabetes <- isTRUE_vector(data$diabetes_dx) | isTRUE_vector(data$diabetes_medication) |
        isTRUE_vector(data$lab_diabetes_a1c)
    data$prediabetes <- !isTRUE_vector(data$diabetes) & ((data$LBXGH >= 5.7 & data$LBXGH < 6.5) |
        col_or_na(data, "DIQ160") == 1 | isTRUE_vector(data$borderline_diabetes_dx))
    data$lipid_lowering_med_proxy <- yes_no(col_or_na(data, "BPQ101D"))
    data$high_cholesterol_dx <- yes_no(col_or_na(data, "BPQ080"))
    data$low_hdl <- (data$RIAGENDR == 1 & data$LBDHDD < 40) | (data$RIAGENDR == 2 & data$LBDHDD <
        50)
    data$high_total_cholesterol <- data$LBXTC >= 240
    data$dyslipidemia_cholesterol_proxy <- isTRUE_vector(data$low_hdl) | isTRUE_vector(data$high_total_cholesterol) |
        isTRUE_vector(data$high_cholesterol_dx) | isTRUE_vector(data$lipid_lowering_med_proxy)
    data$ckd <- data$egfr < 60 | data$uacr >= 30
    data$ckd_very_high_risk_proxy <- data$egfr < 30 | (data$egfr >= 30 & data$egfr < 45 & data$uacr >=
        30) | (data$egfr >= 45 & data$egfr < 60 & data$uacr >= 300)
    data$kidney_condition_dx <- yes_no(col_or_na(data, "KIQ022"))
    data$stage4_clinical_cvd <- row_any_yes(data, c("MCQ160B", "MCQ160C", "MCQ160D", "MCQ160E",
        "MCQ160F"))
    data$general_obesity <- data$BMXBMI >= 30
    data$overweight_or_obesity <- data$BMXBMI >= 25
    data$abdominal_obesity <- (data$RIAGENDR == 1 & data$BMXWAIST > 102) | (data$RIAGENDR == 2 &
        data$BMXWAIST > 88)
    data$adiposity_stage1 <- isTRUE_vector(data$overweight_or_obesity) | isTRUE_vector(data$abdominal_obesity)
    metabolic_syndrome_components <- cbind(abdominal_obesity = isTRUE_vector(data$abdominal_obesity),
        hypertension = isTRUE_vector(data$hypertension), dysglycemia = isTRUE_vector(data$prediabetes) |
            isTRUE_vector(data$diabetes), low_hdl = isTRUE_vector(data$low_hdl))
    data$metabolic_syndrome_no_tg_proxy <- rowSums(metabolic_syndrome_components, na.rm = TRUE) >=
        3
    data$current_smoking <- dplyr::case_when(col_or_na(data, "SMQ020") == 2 ~ FALSE, col_or_na(data,
        "SMQ020") == 1 & col_or_na(data, "SMQ040") %in% c(1, 2) ~ TRUE, col_or_na(data, "SMQ020") ==
        1 & col_or_na(data, "SMQ040") == 3 ~ FALSE, TRUE ~ NA)
    data$stage2_metabolic_risk_or_ckd <- isTRUE_vector(data$hypertension) | isTRUE_vector(data$diabetes) |
        isTRUE_vector(data$ckd) | isTRUE_vector(data$metabolic_syndrome_no_tg_proxy)
    data$stage1_adiposity_or_prediabetes <- isTRUE_vector(data$adiposity_stage1) | isTRUE_vector(data$prediabetes)
    data$ckm_primary_complete <- stats::complete.cases(data[, c("RIDAGEYR", "RIAGENDR", "BMXBMI",
        "BMXWAIST", "mean_sbp", "mean_dbp", "LBXGH", "LBXTC", "LBDHDD", "LBXSCR", "egfr", "uacr",
        "hypertension_dx", "antihypertensive_use", "lipid_lowering_med_proxy", "stage4_clinical_cvd"),
        drop = FALSE])
    data$undiagnosed_hypertension <- isTRUE_vector(data$measured_hypertension) & data$hypertension_dx ==
        FALSE
    data$hypertension_dx_or_med <- isTRUE_vector(data$hypertension_dx) | isTRUE_vector(data$antihypertensive_use)
    data$bp_uncontrolled_among_diagnosed <- isTRUE_vector(data$hypertension_dx_or_med) & (data$mean_sbp >=
        130 | data$mean_dbp >= 80)
    data$undiagnosed_diabetes <- isTRUE_vector(data$lab_diabetes_a1c) & data$diabetes_dx == FALSE
    data$a1c_uncontrolled_among_diagnosed <- isTRUE_vector(data$diabetes_dx) & data$LBXGH >= 7
    data$unrecognized_ckd <- isTRUE_vector(data$ckd) & data$kidney_condition_dx == FALSE
    data$no_health_insurance <- col_or_na(data, "HIQ011") == 2
    data$no_usual_care_place <- col_or_na(data, "HUQ030") == 2
    data$low_or_very_low_food_security <- col_or_na(data, "FSDAD") %in% c(3, 4)
    data
}
