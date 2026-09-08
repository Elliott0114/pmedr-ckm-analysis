state_fips_lookup <- function() {
    data.frame(state_fips = c(1, 2, 4, 5, 6, 8, 9, 10, 11, 12, 13, 15, 16, 17, 18, 19, 20, 21, 22,
        23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 44, 45,
        46, 47, 48, 49, 50, 51, 53, 54, 55, 56, 66, 72, 78), state_abbr = c("AL", "AK", "AZ", "AR",
        "CA", "CO", "CT", "DE", "DC", "FL", "GA", "HI", "ID", "IL", "IN", "IA", "KS", "KY", "LA",
        "ME", "MD", "MA", "MI", "MN", "MS", "MO", "MT", "NE", "NV", "NH", "NJ", "NM", "NY", "NC",
        "ND", "OH", "OK", "OR", "PA", "RI", "SC", "SD", "TN", "TX", "UT", "VT", "VA", "WA", "WV",
        "WI", "WY", "GU", "PR", "VI"), state_name = c("Alabama", "Alaska", "Arizona", "Arkansas",
        "California", "Colorado", "Connecticut", "Delaware", "District of Columbia", "Florida",
        "Georgia", "Hawaii", "Idaho", "Illinois", "Indiana", "Iowa", "Kansas", "Kentucky", "Louisiana",
        "Maine", "Maryland", "Massachusetts", "Michigan", "Minnesota", "Mississippi", "Missouri",
        "Montana", "Nebraska", "Nevada", "New Hampshire", "New Jersey", "New Mexico", "New York",
        "North Carolina", "North Dakota", "Ohio", "Oklahoma", "Oregon", "Pennsylvania", "Rhode Island",
        "South Carolina", "South Dakota", "Tennessee", "Texas", "Utah", "Vermont", "Virginia", "Washington",
        "West Virginia", "Wisconsin", "Wyoming", "Guam", "Puerto Rico", "US Virgin Islands"), stringsAsFactors = FALSE)
}

fixed_levels <- function() {
    list(sex = c("Female", "Male"), race = c("Non-Hispanic White", "Non-Hispanic Black", "Non-Hispanic Asian",
        "Hispanic", "Other/Multi-racial"), bmi = c("Normal/underweight", "Overweight", "Obesity"),
        education = c("College graduate+", "Some college/AA", "High school/GED", "<High school"))
}

binary_yes_no <- function(x, yes = 1, no = 2) {
    out <- rep(NA_integer_, length(x))
    out[x %in% yes] <- 1L
    out[x %in% no] <- 0L
    out
}

brfss_high_cholesterol_dx <- function(chol_check, told_high) {
    out <- rep(NA_integer_, length(told_high))
    out[told_high == 1] <- 1L
    out[told_high == 2 | chol_check == 1] <- 0L
    out
}

prepare_nhanes_brfss_bridge <- function(ckm) {
    levels <- fixed_levels()
    out <- data.frame(SEQN = ckm$SEQN, advanced_ckm = ckm$advanced_ckm_stage_3_4_num, weight = ckm$WTPH2YR_CKM,
        strata = ckm$SDMVSTRA, psu = ckm$SDMVPSU, age_years = pmin(ckm$RIDAGEYR, 80), sex = factor(ifelse(ckm$RIAGENDR ==
            1, "Male", "Female"), levels = levels$sex), race = factor(dplyr::case_when(ckm$RIDRETH3 ==
            3 ~ "Non-Hispanic White", ckm$RIDRETH3 == 4 ~ "Non-Hispanic Black", ckm$RIDRETH3 ==
            6 ~ "Non-Hispanic Asian", ckm$RIDRETH3 %in% c(1, 2) ~ "Hispanic", ckm$RIDRETH3 == 7 ~
            "Other/Multi-racial", TRUE ~ NA_character_), levels = levels$race), bmi_group = factor(dplyr::case_when(ckm$BMXBMI <
            25 ~ "Normal/underweight", ckm$BMXBMI >= 25 & ckm$BMXBMI < 30 ~ "Overweight", ckm$BMXBMI >=
            30 ~ "Obesity", TRUE ~ NA_character_), levels = levels$bmi), education = factor(dplyr::case_when(ckm$DMDEDUC2 ==
            5 ~ "College graduate+", ckm$DMDEDUC2 == 4 ~ "Some college/AA", ckm$DMDEDUC2 == 3 ~
            "High school/GED", ckm$DMDEDUC2 %in% c(1, 2) ~ "<High school", TRUE ~ NA_character_),
            levels = levels$education), sr_hypertension = as.integer(ckm$hypertension_dx), sr_diabetes = as.integer(ckm$diabetes_dx),
        sr_prediabetes = as.integer(!ckm$diabetes_dx & (ckm$borderline_diabetes_dx | ckm$DIQ160 ==
            1)), sr_ckd = as.integer(ckm$kidney_condition_dx), sr_cvd = as.integer(row_any_yes(ckm,
            c("MCQ160C", "MCQ160E", "MCQ160F"))), sr_high_chol = as.integer(ckm$high_cholesterol_dx),
        current_smoking = as.integer(ckm$current_smoking), no_health_insurance = as.integer(ckm$no_health_insurance),
        no_usual_care = as.integer(ckm$no_usual_care_place), general_health_fair_poor = as.integer(ckm$HUQ010 %in%
            c(4, 5)))
    out$age_years_sq <- out$age_years^2
    out
}

prepare_brfss_bridge <- function(brfss) {
    levels <- fixed_levels()
    lookup <- state_fips_lookup()
    out <- data.frame(state_fips = brfss$`_STATE`, brfss_weight = brfss$`_LLCPWT`, brfss_strata = brfss$`_STSTR`,
        brfss_psu = brfss$`_PSU`, age_years = pmin(brfss$`_AGE80`, 80), sex = factor(dplyr::case_when(brfss$SEXVAR ==
            1 ~ "Male", brfss$SEXVAR == 2 ~ "Female", TRUE ~ NA_character_), levels = levels$sex),
        race = factor(dplyr::case_when(brfss$`_IMPRACE` == 1 ~ "Non-Hispanic White", brfss$`_IMPRACE` ==
            2 ~ "Non-Hispanic Black", brfss$`_IMPRACE` == 3 ~ "Non-Hispanic Asian", brfss$`_IMPRACE` ==
            5 ~ "Hispanic", brfss$`_IMPRACE` %in% c(4, 6) ~ "Other/Multi-racial", TRUE ~ NA_character_),
            levels = levels$race), bmi_group = factor(dplyr::case_when(brfss$`_BMI5CAT` %in% c(1,
            2) ~ "Normal/underweight", brfss$`_BMI5CAT` == 3 ~ "Overweight", brfss$`_BMI5CAT` ==
            4 ~ "Obesity", TRUE ~ NA_character_), levels = levels$bmi), education = factor(dplyr::case_when(brfss$`_EDUCAG` ==
            4 ~ "College graduate+", brfss$`_EDUCAG` == 3 ~ "Some college/AA", brfss$`_EDUCAG` ==
            2 ~ "High school/GED", brfss$`_EDUCAG` == 1 ~ "<High school", TRUE ~ NA_character_),
            levels = levels$education), sr_hypertension = binary_yes_no(brfss$BPHIGH6, yes = 1,
            no = c(2, 3, 4)), sr_diabetes = binary_yes_no(brfss$DIABETE4, yes = 1, no = c(2, 3,
            4)), sr_prediabetes = binary_yes_no(brfss$DIABETE4, yes = 4, no = c(1, 2, 3)), sr_ckd = binary_yes_no(brfss$CHCKDNY2,
            yes = 1, no = 2), sr_cvd = as.integer(row_any_yes(brfss, c("CVDINFR4", "CVDCRHD4", "CVDSTRK3"))),
        sr_high_chol = brfss_high_cholesterol_dx(brfss$CHOLCHK3, brfss$TOLDHI3), current_smoking = binary_yes_no(brfss$`_SMOKER3`,
            yes = c(1, 2), no = c(3, 4)), no_health_insurance = binary_yes_no(brfss$`_HLTHPL1`,
            yes = 2, no = 1), no_usual_care = binary_yes_no(brfss$PERSDOC3, yes = 3, no = c(1, 2)),
        general_health_fair_poor = binary_yes_no(brfss$GENHLTH, yes = c(4, 5), no = c(1, 2, 3)))
    out$age_years_sq <- out$age_years^2
    out <- dplyr::left_join(out, lookup, by = "state_fips")
    out
}
