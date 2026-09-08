source("R/setup.R", encoding = "UTF-8")

args <- commandArgs(trailingOnly = TRUE)

stopifnot(length(args) == 2L)

d <- readRDS(args[1])

out <- new_output_directory(args[2])

d$advanced_ckm_stage_3_4_num <- as.integer(d$ckm_stage_observed %in% 3:4)

nh <- prepare_nhanes_brfss_bridge(d)

nh$outcome <- as.integer(d$ckm_stage_observed == 3)

nh$known_cvd <- d$stage4_clinical_cvd

nh$prevent_age_eligible <- as.integer(nh$age_years >= 30 & nh$age_years <= 79)

formulas <- state_formulas()

predictors <- all.vars(delete.response(terms(formulas$lean_core)))

keep <- nh$known_cvd %in% FALSE & complete.cases(nh[, c("outcome", "weight", "strata", "psu", predictors)]) &
    nh$weight > 0

train <- nh[keep, ]

stopifnot(nrow(train) == 4395L, sum(train$outcome) == 293L)

models <- lapply(formulas, function(f) fit_survey_logistic(train, f))

stopifnot(all(vapply(models, function(m) isTRUE(m$converged), logical(1))))

saveRDS(train, file.path(out, "training_data.rds"))

saveRDS(models, file.path(out, "models.rds"))

structure <- bind_rows(lapply(names(models), function(n) model_structure_row(n, models[[n]], train,
    formulas[[n]])))

coefficients <- bind_rows(lapply(names(models), function(n) data.frame(model = n, term = names(coef(models[[n]])),
    coefficient = as.numeric(coef(models[[n]])), se = sqrt(diag(vcov(models[[n]]))))))

readr::write_csv(structure, file.path(out, "model_structure.csv"))

readr::write_csv(coefficients, file.path(out, "model_coefficients.csv"))

readr::write_csv(verify_fixed_prediction(models$lean_core, train), file.path(out, "prediction_basis_verification.csv"))

excluded <- nh[!nh$known_cvd %in% TRUE & !keep, ]

needed <- c("outcome", "weight", "strata", "psu", predictors)

excluded$reason <- apply(is.na(excluded[, needed]), 1, function(a) paste(needed[a], collapse = "; "))

reasons <- summarise(group_by(excluded, reason), n = n(), stage3 = sum(outcome == 1, na.rm = TRUE),
    .groups = "drop")

readr::write_csv(reasons, file.path(out, "training_exclusions.csv"))

capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
