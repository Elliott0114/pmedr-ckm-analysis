source("R/setup.R", encoding = "UTF-8")

args <- commandArgs(trailingOnly = TRUE)

stopifnot(length(args) == 3L)

d <- readRDS(args[1])

raw_dir <- normalizePath(args[2], winslash = "/", mustWork = TRUE)

out <- new_output_directory(args[3])

write_table <- function(d, n) readr::write_csv(d, file.path(out, n), na = "NA")

st0 <- d$ckm_stage_observed

stopifnot(nrow(d) == 5006L, !anyDuplicated(d$SEQN), sum(st0 %in% 3) == 294L, sum(st0 %in% 4) ==
    601L, sum(is.na(st0)) == 4L)

methods <- c(CC = "WTPH2YR_CKM", IPW = "WTPH2YR_CKM_IPW_TRUNC")

master <- list()

add <- function(x) master[[length(master) + 1L]] <<- x

message("National proportions and reliability")

for (method in names(methods)) {
    w <- methods[[method]]
    des <- design_of(d, w)
    add(stage_rows(d, st0, w = w, method = method))
    for (j in seq_len(nrow(specs))) for (st in c("2", "3", "4", "3-4")) {
        keep <- if (st == "3-4")
            st0 %in% 3:4
        else st0 %in% as.integer(st)
        add(prop_row(des, keep & is_yes(d[[specs$den[j]]]), d[[specs$endpoint[j]]], "Markers", "Primary",
            method, st, specs$endpoint[j]))
    }
    den <- list(kidney = d$ckd, bp = d$hypertension_dx_or_med, a1c = d$diabetes_dx, lipid = d$high_total_cholesterol |
        d$high_cholesterol_dx, food = rep(TRUE, nrow(d)), care = rep(TRUE, nrow(d)), insurance = rep(TRUE,
        nrow(d)), bp140 = d$hypertension_dx_or_med, a1c8 = d$diabetes_dx)
    val <- list(kidney = d$kidney_condition_dx == FALSE, bp = d$mean_sbp >= 130 | d$mean_dbp >=
        80, a1c = d$LBXGH >= 7, lipid = d$lipid_lowering_med_proxy == FALSE, food = ifelse(d$FSDAD %in%
        1:4, d$FSDAD %in% 3:4, NA), care = ifelse(d$HUQ030 %in% 1:3, d$HUQ030 == 2, NA), insurance = ifelse(d$HIQ011 %in%
        1:2, d$HIQ011 == 2, NA), bp140 = d$mean_sbp >= 140 | d$mean_dbp >= 90, a1c8 = d$LBXGH >=
        8)
    for (k in names(val)) add(prop_row(des, st0 %in% 3:4 & is_yes(den[[k]]), val[[k]], "Original domains",
        "Primary", method, "3-4", k))
    for (age in levels(d$age_group)) for (sex in levels(d$sex_label)) add(prop_row(des, d$age_group %in%
        age & d$sex_label %in% sex, st0 %in% 3:4, "Age-sex", "Primary", method, paste(age, sex,
        sep = ":"), "advanced"))
    kd <- list(`eGFR <60 or UACR >=30` = d$egfr < 60 | d$uacr >= 30, `eGFR <60` = d$egfr < 60, `UACR >=30` = d$uacr >=
        30, `eGFR <45 or UACR >=300` = d$egfr < 45 | d$uacr >= 300, `eGFR <45` = d$egfr < 45, `UACR >=300` = d$uacr >=
        300, `eGFR <60 and UACR >=30` = d$egfr < 60 & d$uacr >= 30)
    for (k in names(kd)) add(prop_row(des, st0 %in% 3:4 & is_yes(kd[[k]]), d$kidney_condition_dx ==
        FALSE, "Kidney severity", "Primary", method, "3-4", k))
    for (st in 2:4) for (pat in unique(d$domain_pattern)) add(prop_row(des, st0 %in% st, d$domain_pattern ==
        pat, "Domain composition", "Primary", method, as.character(st), ifelse(pat == "", "None",
        pat)))
    for (src in c("PREVENT only", "Kidney only", "Both")) {
        flag <- d$source3 %in% src
        add(prop_row(des, st0 %in% 3, flag, "Stage-3 pathways", "Primary", method, src, "source_fraction"))
        for (j in 1:5) add(prop_row(des, flag & is_yes(d[[specs$den[j]]]), d[[specs$endpoint[j]]],
            "Stage-3 pathway markers", "Primary", method, src, specs$endpoint[j]))
    }
}

message("PREVENT base-only and raw-range sensitivities")

base_risk <- risks(d, optional = FALSE)

write_table(base_risk, "prevent_base_records.csv")

raw_risk <- risks(d, optional = TRUE, bounded = FALSE)

write_table(raw_risk, "prevent_raw_range_records.csv")

base_stage <- assign_stage(d, base_risk$risk)

raw_stage <- assign_stage(d, raw_risk$risk)

kid_stage <- assign_stage(d, rep(NA_real_, nrow(d)), kidney_only = TRUE)

write_table(as.data.frame(table(primary = st0, base_only = base_stage, useNA = "always")), "stage_transition_base.csv")

for (method in names(methods)) {
    w <- methods[[method]]
    des <- design_of(d, w)
    for (route in c("Base PREVENT", "Raw ranges", "Kidney-only stage 3")) {
        st <- switch(route, `Base PREVENT` = base_stage, `Raw ranges` = raw_stage, `Kidney-only stage 3` = kid_stage)
        add(stage_rows(d, st, route, w, method))
    }
    for (j in c(1, 2, 3, 4, 5, 10, 11)) for (st in c("3", "4", "3-4")) {
        keep <- if (st == "3-4")
            base_stage %in% 3:4
        else base_stage %in% as.integer(st)
        add(prop_row(des, keep & is_yes(d[[specs$den[j]]]), d[[specs$endpoint[j]]], "Markers", "Base PREVENT",
            method, st, specs$endpoint[j]))
    }
}

stopifnot(sum(is.na(raw_stage)) == 574)

message("Fasting sensitivity with matched sample and fasting weights")

tg <- haven::zap_labels(haven::read_xpt(file.path(raw_dir, "TRIGLY_L.xpt")))

gl <- haven::zap_labels(haven::read_xpt(file.path(raw_dir, "GLU_L.xpt")))

stopifnot(!anyDuplicated(tg$SEQN), !anyDuplicated(gl$SEQN))

fd <- d

fd$WTSAF2YR <- tg$WTSAF2YR[match(d$SEQN, tg$SEQN)]

fd$LBXTLG <- tg$LBXTLG[match(d$SEQN, tg$SEQN)]

fd$LBXGLU <- gl$LBXGLU[match(d$SEQN, gl$SEQN)]

keep <- is.finite(fd$WTSAF2YR) & fd$WTSAF2YR > 0 & complete.cases(fd[, c("LBXTLG", "LBXGLU", "SDMVSTRA",
    "SDMVPSU")])

fd <- fd[keep, ]

fd$diabetes <- is_yes(fd$diabetes) | fd$LBXGLU >= 126

fd$prediabetes <- !fd$diabetes & (is_yes(fd$prediabetes) | (fd$LBXGLU >= 100 & fd$LBXGLU < 126))

ms <- rowSums(cbind(is_yes(fd$abdominal_obesity), is_yes(fd$hypertension), fd$diabetes | fd$prediabetes,
    is_yes(fd$low_hdl), fd$LBXTLG >= 150)) >= 3

f1 <- is_yes(fd$adiposity_stage1) | fd$prediabetes

f2 <- is_yes(fd$hypertension) | fd$diabetes | is_yes(fd$ckd) | fd$LBXTLG >= 150 | ms

frisk <- risks(fd)

fst <- assign_stage(fd, frisk$risk, f1, f2)

add(stage_rows(fd, fd$ckm_stage_observed, "Primary in fasting sample", "WTSAF2YR", "Fasting weight"))

add(stage_rows(fd, fst, "TG/glucose augmented", "WTSAF2YR", "Fasting weight"))

write_table(data.frame(SEQN = fd$SEQN, primary_stage = fd$ckm_stage_observed, fasting_stage = fst,
    primary_diabetes = d$diabetes[match(fd$SEQN, d$SEQN)], fasting_diabetes = fd$diabetes, primary_risk = fd$prevent_total_cvd_10yr,
    fasting_risk = frisk$risk, WTSAF2YR = fd$WTSAF2YR), "fasting_classification.csv")

write_table(data.frame(n_primary = nrow(d), n_fasting = nrow(fd), new_diabetes = sum(fd$diabetes &
    !is_yes(d$diabetes[match(fd$SEQN, d$SEQN)])), stage_changes = sum(fst != fd$ckm_stage_observed,
    na.rm = TRUE), unresolved = sum(is.na(fst))), "fasting_flow.csv")

message("Restricted-support stage comparisons")

support <- list()

contrasts <- list()

for (method in names(methods)) for (j in seq_len(nrow(specs))) {
    endpoint <- specs$endpoint[j]
    eligible <- st0 %in% 3:4 & is_yes(d[[specs$den[j]]]) & !is.na(d[[endpoint]])
    ageok <- eligible & d$RIDAGEYR >= 45 & d$RIDAGEYR <= 79
    cell <- paste(d$RIAGENDR, d$race_ethnicity, sep = ":")
    both <- intersect(unique(cell[ageok & st0 %in% 3]), unique(cell[ageok & st0 %in% 4]))
    keep <- ageok & cell %in% both
    w <- methods[[method]]
    for (st in 3:4) {
        z <- eligible & st0 %in% st
        k <- keep & st0 %in% st
        support[[length(support) + 1]] <- data.frame(method, endpoint, stage = st, eligible_n = sum(z),
            retained_n = sum(k), excluded_n = sum(z) - sum(k), weighted_retained = sum(d[[w]][k])/sum(d[[w]][z]),
            retained_cells = paste(both, collapse = ";"))
    }
    dd <- design_of(d, w)[keep, ]
    dd$variables$stage <- droplevels(factor(st0[keep], levels = 3:4))
    dd$variables$age10 <- (dd$variables$RIDAGEYR - 60)/10
    dd$variables$sex_adjust <- droplevels(factor(dd$variables$RIAGENDR))
    dd$variables$race_adjust <- droplevels(factor(dd$variables$race_ethnicity))
    n <- sum(keep)
    ev <- sum(d[[endpoint]][keep])
    model <- tryCatch(survey::svyglm(as.formula(paste0(endpoint, " ~ stage+age10+sex_adjust+race_adjust")),
        dd, family = quasibinomial()), error = function(e) e)
    fail <- ""
    vals <- rep(NA_real_, 6)
    names(vals) <- c("pd34", "pd34_se", "pd34_low", "pd34_high", "p3", "p4")
    slopes <- rdf <- rank <- NA_real_
    if (inherits(model, "error"))
        fail <- conditionMessage(model)
    else {
        slopes <- length(coef(model)) - 1
        rdf <- model$df.residual
        rank <- qr(vcov(model))$rank
        if (!isTRUE(model$converged))
            fail <- "Nonconvergence"
        if (any(!is.finite(coef(model))) || any(!is.finite(vcov(model))))
            fail <- "Non-estimable coefficient/variance"
        if (rank < length(coef(model)))
            fail <- "Covariance rank deficient"
        if (rdf <= 0)
            fail <- "Residual design df<=0"
        if (fail == "")
            vals <- adjusted_difference(model, dd$variables, weights(dd))
    }
    contrasts[[length(contrasts) + 1]] <- data.frame(method, endpoint, n, events = ev, slopes, design_df = survey::degf(dd),
        residual_df = rdf, covariance_rank = rank, reportable = fail == "", reason = fail, as.list(vals),
        check.names = FALSE)
    for (st in 3:4) add(prop_row(design_of(d, w), keep & st0 %in% st, d[[endpoint]], "Restricted support",
        "Ages 45-79; common sex-race cells", method, as.character(st), endpoint))
}

write_table(bind_rows(support), "support_flow.csv")

write_table(bind_rows(contrasts), "support_contrasts.csv")

message("Equity models on the same complete-case sample")

ed <- d

ed$advanced_ckm_stage_3_4_num <- ifelse(is.na(st0), NA_integer_, as.integer(st0 %in% 3:4))

ed$sex_label <- factor(ed$sex_label, levels = c("Female", "Male"))

ed$race_ethnicity <- factor(ed$race_ethnicity, levels = c("Non-Hispanic White", "Non-Hispanic Black",
    "Mexican American", "Other Hispanic", "Non-Hispanic Asian", "Other/Multi-racial"))

ed$education_group <- factor(ed$education_group, levels = c("College graduate+", "Some college/AA",
    "High school/GED", "<High school"))

ed$pir_group_model <- factor(ed$pir_group, levels = c(">=3.50", "1.30-3.49", "<1.30"))

ed$low_or_very_low_food_security_num <- ifelse(ed$FSDAD %in% 1:4, as.integer(ed$FSDAD %in% 3:4),
    NA_integer_)

ed$no_health_insurance_num <- ifelse(ed$HIQ011 %in% 1:2, as.integer(ed$HIQ011 == 2), NA_integer_)

vars <- c("advanced_ckm_stage_3_4_num", "age_group", "sex_label", "race_ethnicity", "education_group",
    "pir_group_model", "no_health_insurance_num", "low_or_very_low_food_security_num")

eqdes <- design_of(ed)[complete.cases(ed[, vars]), ]

stopifnot(nrow(eqdes$variables) == 4382, sum(eqdes$variables$advanced_ckm_stage_3_4_num) == 770)

forms <- list(Full = advanced_ckm_stage_3_4_num ~ age_group + sex_label + race_ethnicity + education_group +
    pir_group_model + no_health_insurance_num + low_or_very_low_food_security_num)

for (x in c("education_group", "pir_group_model", "no_health_insurance_num", "low_or_very_low_food_security_num")) forms[[paste0("Demographic_",
    x)]] <- as.formula(paste0("advanced_ckm_stage_3_4_num ~ splines::ns(RIDAGEYR,df=2)+sex_label+race_ethnicity+",
    x))

eq <- list()

eqstr <- list()

for (nm in names(forms)) {
    f <- survey::svyglm(forms[[nm]], eqdes, family = quasipoisson(link = "log"))
    b <- coef(f)
    se <- sqrt(diag(vcov(f)))
    df <- survey::degf(eqdes)
    crit <- qt(0.975, df)
    stopifnot(isTRUE(f$converged), all(is.finite(b)), all(is.finite(se)))
    eq[[nm]] <- data.frame(model = nm, term = names(b), pr = exp(b), low = exp(b - crit * se), high = exp(b +
        crit * se), se_log = se, design_df = df, ci_method = "Individual-coefficient t(design df) approximation")
    eqstr[[nm]] <- data.frame(model = nm, n = nrow(eqdes$variables), events = sum(eqdes$variables$advanced_ckm_stage_3_4_num),
        slopes = length(b) - 1, design_df = df, residual_df = f$df.residual, covariance_rank = qr(vcov(f))$rank,
        converged = f$converged)
}

eq <- bind_rows(eq)

write_table(eq, "equity_models.csv")

write_table(bind_rows(eqstr), "equity_model_structure.csv")

ans <- dplyr::bind_rows(master)

stopifnot(all(ans$low[is.finite(ans$low)] >= 0), all(ans$high[is.finite(ans$high)] <= 1))

joint <- filter(ans, family == "Markers", variant == "Primary", method == "CC", group == "3-4",
    endpoint == "joint_primary")

stopifnot(joint$n == 277, joint$events == 91)

write_table(ans, "national_results.csv")

write_table(stage_contrasts(d), "stage_contrasts.csv")

capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
