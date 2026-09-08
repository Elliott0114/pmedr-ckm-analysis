state_formulas <- function() {
    list(lean_core = outcome ~ splines::ns(age_years, df = 2) + prevent_age_eligible + sex + race +
        sr_hypertension + sr_diabetes + sr_ckd + current_smoking, lean_without_reported_ckd = outcome ~
        splines::ns(age_years, df = 2) + prevent_age_eligible + sex + race + sr_hypertension + sr_diabetes +
            current_smoking, demographic_only = outcome ~ splines::ns(age_years, df = 2) + prevent_age_eligible +
        sex + race)
}
