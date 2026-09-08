col_or_na <- function(data, name) {
    if (name %in% names(data)) {
        return(data[[name]])
    }
    rep(NA_real_, nrow(data))
}

mean_available <- function(data, vars) {
    vars <- intersect(vars, names(data))
    if (length(vars) == 0L) {
        return(rep(NA_real_, nrow(data)))
    }
    values <- data[, vars, drop = FALSE]
    out <- rowMeans(values, na.rm = TRUE)
    out[rowSums(!is.na(values)) == 0L] <- NA_real_
    out
}

yes_no <- function(x) {
    out <- rep(NA, length(x))
    out[x == 1] <- TRUE
    out[x == 2] <- FALSE
    out
}

row_any_yes <- function(data, vars) {
    vars <- intersect(vars, names(data))
    if (length(vars) == 0L) {
        return(rep(NA, nrow(data)))
    }
    mat <- as.data.frame(lapply(data[, vars, drop = FALSE], yes_no))
    out <- rep(FALSE, nrow(data))
    has_unknown <- rep(FALSE, nrow(data))
    for (v in names(mat)) {
        out <- out | isTRUE_vector(mat[[v]])
        has_unknown <- has_unknown | is.na(mat[[v]])
    }
    out[!out & has_unknown] <- NA
    out
}

isTRUE_vector <- function(x) {
    !is.na(x) & x
}

egfr_ckdepi_2021 <- function(scr_mg_dl, age_years, sex_code) {
    female <- sex_code == 2
    kappa <- ifelse(female, 0.7, 0.9)
    alpha <- ifelse(female, -0.241, -0.302)
    ratio <- scr_mg_dl/kappa
    egfr <- 142 * pmin(ratio, 1, na.rm = FALSE)^alpha * pmax(ratio, 1, na.rm = FALSE)^(-1.2) * 0.9938^age_years *
        ifelse(female, 1.012, 1)
    egfr[is.na(scr_mg_dl) | is.na(age_years) | is.na(sex_code)] <- NA_real_
    egfr
}

uacr_mg_g <- function(data) {
    direct <- col_or_na(data, "URDACT")
    albumin_mg_l <- col_or_na(data, "URXUMA")
    creatinine_mg_dl <- col_or_na(data, "URXUCR")
    calculated <- albumin_mg_l/creatinine_mg_dl * 100
    ifelse(!is.na(direct), direct, calculated)
}
