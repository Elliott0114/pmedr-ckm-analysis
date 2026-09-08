# Proportions use Korn-Graubard beta intervals with the domain degrees-of-freedom
# adjustment. Clinical values stay on their measured scale; bounds apply only
# to PREVENT risk inputs. The adapter retains the preventr 0.12.0 equations.
# Stage 3-4 contrasts use the fitted model residual degrees of freedom.

options(stringsAsFactors = FALSE, survey.lonely.psu = "adjust")

is_yes <- function(x) x %in% TRUE

design_of <- function(d, w = "WTPH2YR_CKM") survey::svydesign(ids = ~SDMVPSU, strata = ~SDMVSTRA,
    weights = as.formula(paste0("~", w)), data = d, nest = TRUE)

kg_interval <- function(p, v, n, df) {
    if (!is.finite(p) || p <= 0 || p >= 1 || !is.finite(v) || v <= 0 || n <= 1 || df <= 0)
        return(c(low = NA, high = NA, neff = NA, neff_adj = NA))
    raw <- p * (1 - p)/v
    ne <- min(n, raw)
    nadj <- min(n, raw * (qt(0.975, n - 1)/qt(0.975, df))^2)
    x <- nadj * p
    c(low = qbeta(0.025, x, nadj - x + 1), high = qbeta(0.975, x + 1, nadj - x), neff = ne, neff_adj = nadj)
}

prop_row <- function(design, keep, value, family, variant = "Primary", method = "CC", group = "",
    endpoint = "") {
    eligible <- keep %in% TRUE
    good <- eligible & !is.na(value)
    dd <- design[good, ]
    n <- sum(good)
    events <- sum(value[good])
    df <- survey::degf(dd)
    if (n == 0)
        return(data.frame(family, variant, method, group, endpoint, eligible_n = sum(eligible),
            n = 0, events = 0, estimate = NA, se = NA, low = NA, high = NA, wald_low = NA, wald_high = NA,
            design_df = df, neff = NA, neff_adj = NA, ci_width = NA, relative_width = NA, reliability = "No observed denominator"))
    dd$variables$.audit_flag <- as.numeric(value[good])
    a <- survey::svymean(~.audit_flag, dd)
    p <- unname(coef(a)[1])
    se <- unname(survey::SE(a)[1])
    ci <- kg_interval(p, se^2, n, df)
    width <- unname(ci["high"] - ci["low"])
    rw <- width/p
    reason <- character()
    if (p == 0 || p == 1)
        reason <- c(reason, "Boundary proportion; interval not estimable")
    if (n < 30)
        reason <- c(reason, "N<30")
    if (is.finite(ci["neff"]) && ci["neff"] < 30)
        reason <- c(reason, "Effective N<30")
    if (is.finite(width) && width >= 0.3)
        reason <- c(reason, "CI width>=30pp")
    if (is.finite(width) && width > 0.05 && width < 0.3 && rw > 1.3)
        reason <- c(reason, "Relative CI width>130%")
    if (p > 0.5 && p < 1 && is.finite(width) && width > 0.05 && width < 0.3 && width/(1 - p) > 1.3)
        reason <- c(reason, "Complement relative CI width>130%")
    if (df < 8)
        reason <- c(reason, "Design df<8; interpretation requires caution")
    if (!is.finite(ci["low"]) && p > 0 && p < 1)
        reason <- c(reason, "Interval not estimable")
    if (!length(reason))
        reason <- "Meets numerical presentation criteria"
    data.frame(family, variant, method, group, endpoint, eligible_n = sum(eligible), n, events,
        estimate = p, se, low = ci["low"], high = ci["high"], wald_low = p - qnorm(0.975) * se,
        wald_high = p + qnorm(0.975) * se, design_df = df, neff = ci["neff"], neff_adj = ci["neff_adj"],
        ci_width = width, relative_width = rw, reliability = paste(reason, collapse = "; "), row.names = NULL)
}

aha_env <- new.env(parent = asNamespace("preventr"))

aha_env$is_valid_sbp <- function(sbp, quiet = TRUE) {
    if (is.numeric(sbp) && length(sbp) == 1L && is.finite(sbp) && sbp >= 90 && sbp <= 200)
        TRUE
    else preventr:::is_valid_sbp(sbp, quiet)
}

aha_env$is_valid_hba1c <- function(hba1c, allow_empty = FALSE, quiet = TRUE) {
    if (is.numeric(hba1c) && length(hba1c) == 1L && is.finite(hba1c) && hba1c >= 3 && hba1c <= 15)
        TRUE
    else preventr:::is_valid_hba1c(hba1c, allow_empty, quiet)
}

aha_env$select_model <- getFromNamespace("select_model", "preventr")

environment(aha_env$select_model) <- aha_env

aha_estimate <- preventr::estimate_risk

environment(aha_estimate) <- aha_env

risk_one <- function(d, i, optional = TRUE, bounded = TRUE) {
    adj <- function(x, l, h) if (bounded)
        pmin(pmax(x, l), h)
    else x
    a <- list(age = d$RIDAGEYR[i], sex = d$prevent_sex[i], sbp = adj(d$mean_sbp[i], 90, 200), bp_tx = d$antihypertensive_use[i],
        total_c = adj(d$LBXTC[i], 130, 320), hdl_c = adj(d$LBDHDD[i], 20, 100), statin = d$lipid_lowering_med_proxy[i],
        dm = d$diabetes[i], smoking = d$current_smoking[i], egfr = adj(d$egfr[i], 15, 140), bmi = adj(d$BMXBMI[i],
            18.5, 39.9), time = "10yr", quiet = TRUE, collapse = TRUE)
    if (optional) {
        a$hba1c <- d$LBXGH[i]
        a$uacr <- d$uacr[i]
    }
    r <- do.call(aha_estimate, a)
    list(risk = r$total_cvd[1], model = as.character(r$model[1]), problem = as.character(r$input_problems[1]))
}

risks <- function(d, optional = TRUE, bounded = TRUE) {
    ix <- which(d$RIDAGEYR >= 30 & d$RIDAGEYR <= 79 & d$stage4_clinical_cvd %in% FALSE)
    ans <- data.frame(SEQN = d$SEQN, risk = NA_real_, model = NA_character_, problem = NA_character_)
    for (k in seq_along(ix)) {
        i <- ix[k]
        z <- risk_one(d, i, optional, bounded)
        ans$risk[i] <- z$risk
        ans$model[i] <- z$model
        ans$problem[i] <- z$problem
        if (k%%1000 == 0)
            message("Risk calculation ", k, "/", length(ix), "; optional=", optional, "; bounded=",
                bounded)
    }
    ans
}

assign_stage <- function(d, risk, stage1 = d$stage1_adiposity_or_prediabetes, stage2 = d$stage2_metabolic_risk_or_ckd,
    kidney_only = FALSE) {
    st <- rep(0L, nrow(d))
    st[is_yes(stage1)] <- 1L
    st[is_yes(stage2)] <- 2L
    st[is_yes(d$ckd_very_high_risk_proxy) | (!kidney_only & is_yes(risk >= 0.2))] <- 3L
    if (!kidney_only) {
        unresolved <- d$RIDAGEYR >= 30 & d$RIDAGEYR <= 79 & d$stage4_clinical_cvd %in% FALSE & d$ckd_very_high_risk_proxy %in%
            FALSE & !is.finite(risk)
        st[unresolved] <- NA_integer_
    }
    st[is_yes(d$stage4_clinical_cvd)] <- 4L
    st
}

stage_rows <- function(d, st, variant = "Primary", w = "WTPH2YR_CKM", method = "CC") {
    des <- design_of(d, w)
    dplyr::bind_rows(lapply(c(as.character(0:4), "3-4", "Unresolved"), function(k) {
        val <- if (k == "3-4")
            st %in% 3:4
        else if (k == "Unresolved")
            is.na(st)
        else st %in% as.integer(k)
        prop_row(des, rep(TRUE, nrow(d)), val, "Stage distribution", variant, method, k, "stage")
    }))
}

specs <- data.frame(endpoint = c("bp_gap_130", "bp_gap_140", "a1c_gap_7", "a1c_gap_8", "kidney_gap",
    "primary_gap_ge1", "primary_gap_ge2", "strict_gap_ge1", "strict_gap_ge2", "joint_primary", "joint_strict"),
    den = c("bp_applicable", "bp_applicable", "glycemia_applicable", "glycemia_applicable", "kidney_applicable",
        "app_ge1", "app_ge2", "app_ge1", "app_ge2", "fixed_bp_a1c", "fixed_bp_a1c"))

fixed_prediction_matrix <- function(model, data) {
    tt <- delete.response(terms(model))
    mf <- model.frame(tt, data = data, xlev = model$xlevels, na.action = na.pass)
    model.matrix(tt, mf, contrasts.arg = model$contrasts)[, names(coef(model)), drop = FALSE]
}

adjusted_difference <- function(model, data, weights, stagevar = "stage") {
    a <- b <- data
    a[[stagevar]] <- factor("3", levels = levels(data[[stagevar]]))
    b[[stagevar]] <- factor("4", levels = levels(data[[stagevar]]))
    xa <- fixed_prediction_matrix(model, a)
    xb <- fixed_prediction_matrix(model, b)
    ba <- coef(model)
    pa <- plogis(xa %*% ba)[, 1]
    pb <- plogis(xb %*% ba)[, 1]
    w <- weights/sum(weights)
    grad <- colSums(xa * (w * pa * (1 - pa))) - colSums(xb * (w * pb * (1 - pb)))
    se <- sqrt(as.numeric(t(grad) %*% vcov(model) %*% grad))
    est <- sum(w * (pa - pb))
    crit <- qt(0.975, model$df.residual)
    c(pd34 = est, pd34_se = se, pd34_low = est - crit * se, pd34_high = est + crit * se, p3 = sum(w *
        pa), p4 = sum(w * pb))
}
