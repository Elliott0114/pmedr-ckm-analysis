# Cross-validation and bootstrap resample primary sampling units within strata.

clip_probability <- function(p) {
    pmin(pmax(as.numeric(p), 1e-06), 1 - 1e-06)
}

logit_probability <- function(p) {
    p <- clip_probability(p)
    log(p/(1 - p))
}

weighted_auc <- function(y, p, w) {
    keep <- stats::complete.cases(y, p, w) & is.finite(w) & w > 0
    y <- as.integer(y[keep])
    p <- as.numeric(p[keep])
    w <- as.numeric(w[keep])
    if (length(unique(y)) != 2L)
        return(NA_real_)
    grouped <- dplyr::mutate(dplyr::arrange(dplyr::summarise(dplyr::group_by(data.frame(y = y, p = p,
        w = w), .data$p), case_w = sum(.data$w[.data$y == 1L]), control_w = sum(.data$w[.data$y ==
        0L]), .groups = "drop"), .data$p), control_before = dplyr::lag(cumsum(.data$control_w),
        default = 0))
    sum(grouped$case_w * (grouped$control_before + 0.5 * grouped$control_w))/(sum(grouped$case_w) *
        sum(grouped$control_w))
}

weighted_quantile_group <- function(x, w, groups = 10L) {
    keep <- is.finite(x) & is.finite(w) & w > 0
    out <- rep(NA_integer_, length(x))
    ord <- order(x[keep])
    cumulative_weight <- cumsum(w[keep][ord])/sum(w[keep][ord])
    out[which(keep)[ord]] <- pmin(groups, pmax(1L, ceiling(cumulative_weight * groups)))
    out
}

make_design <- function(data, weight = "weight", strata = "strata", psu = "psu") {
    survey::svydesign(ids = stats::as.formula(paste0("~", psu)), strata = stats::as.formula(paste0("~",
        strata)), weights = stats::as.formula(paste0("~", weight)), data = data, nest = TRUE)
}

metric_set <- function(data, prediction) {
    prediction <- clip_probability(prediction)
    y <- as.integer(data$outcome)
    w <- as.numeric(data$weight)
    link <- logit_probability(prediction)
    calibration_data <- data.frame(outcome = y, link = link, weight = w, strata = data$strata, psu = data$psu)
    calibration_design <- make_design(calibration_data)
    slope_fit <- suppressWarnings(survey::svyglm(outcome ~ link, calibration_design, family = quasibinomial()))
    intercept_fit <- suppressWarnings(survey::svyglm(outcome ~ 1 + offset(link), calibration_design,
        family = quasibinomial()))
    decile <- weighted_quantile_group(prediction, w)
    calibration <- dplyr::summarise(dplyr::group_by(data.frame(y = y, p = prediction, w = w, decile = decile),
        .data$decile), observed = stats::weighted.mean(.data$y, .data$w), predicted = stats::weighted.mean(.data$p,
        .data$w), total_weight = sum(.data$w), .groups = "drop")
    prevalence <- stats::weighted.mean(y, w)
    brier <- stats::weighted.mean((y - prediction)^2, w)
    null_brier <- stats::weighted.mean((y - prevalence)^2, w)
    c(auc = weighted_auc(y, prediction, w), brier = brier, brier_skill_vs_null = 1 - brier/null_brier,
        calibration_intercept = unname(stats::coef(intercept_fit)[1]), calibration_slope = unname(stats::coef(slope_fit)[2]),
        expected_calibration_error = stats::weighted.mean(abs(calibration$observed - calibration$predicted),
            calibration$total_weight), observed_prevalence = prevalence, predicted_prevalence = stats::weighted.mean(prediction,
            w))
}

bootstrap_optimism <- function(data, formula, apparent, reps, seed) {
    strata_values <- unique(data$strata)
    rows <- vector("list", reps)
    for (replicate_id in seq_len(reps)) {
        set.seed(seed + replicate_id)
        sampled <- lapply(strata_values, function(stratum_value) {
            stratum_data <- data[data$strata == stratum_value, , drop = FALSE]
            psus <- unique(stratum_data$psu)
            selected <- sample(psus, length(psus), replace = TRUE)
            dplyr::bind_rows(lapply(seq_along(selected), function(draw_position) {
                unit <- stratum_data[stratum_data$psu == selected[draw_position], , drop = FALSE]
                unit$psu <- paste0(as.character(stratum_value), "_boot_", draw_position)
                unit
            }))
        })
        boot_data <- dplyr::bind_rows(sampled)
        fit <- tryCatch(fit_survey_logistic(boot_data, formula), error = function(e) NULL)
        if (!is.null(fit)) {
            p_train <- tryCatch(as.numeric(stats::predict(fit, newdata = boot_data, type = "response")),
                error = function(e) rep(NA_real_, nrow(boot_data)))
            p_test <- tryCatch(as.numeric(stats::predict(fit, newdata = data, type = "response")),
                error = function(e) rep(NA_real_, nrow(data)))
            if (all(is.finite(p_train)) && all(is.finite(p_test))) {
                train_metric <- metric_set(boot_data, p_train)
                test_metric <- metric_set(data, p_test)
                rows[[replicate_id]] <- data.frame(replicate = replicate_id, metric = names(train_metric),
                  bootstrap_apparent = as.numeric(train_metric), original_test = as.numeric(test_metric),
                  optimism = as.numeric(train_metric - test_metric), corrected = as.numeric(apparent -
                    (train_metric - test_metric)), stringsAsFactors = FALSE)
            }
        }
        if (replicate_id%%10L == 0L)
            message("Lean PSU bootstrap: ", replicate_id, "/", reps)
    }
    dplyr::bind_rows(rows)
}

repeated_psu_cv_predictions <- function(data, formulas, repeats, folds, seed) {
    cluster_map <- unique(data.frame(key = interaction(data$strata, data$psu, drop = TRUE), strata = data$strata,
        psu = data$psu))
    prediction_rows <- vector("list", repeats * length(formulas))
    metric_rows <- vector("list", repeats * length(formulas))
    assignment_rows <- vector("list", repeats)
    row_id <- 0L
    for (repeat_id in seq_len(repeats)) {
        set.seed(seed + repeat_id)
        cluster_map$fold <- NA_integer_
        for (stratum_value in unique(cluster_map$strata)) {
            index <- which(cluster_map$strata == stratum_value)
            cluster_map$fold[index] <- sample(seq_len(folds), length(index), replace = length(index) >
                folds)
        }
        row_key <- interaction(data$strata, data$psu, drop = TRUE)
        fold_id <- cluster_map$fold[match(as.character(row_key), as.character(cluster_map$key))]
        assignment_rows[[repeat_id]] <- data.frame(repeat_id = repeat_id, cluster_map, stringsAsFactors = FALSE)
        for (model_name in names(formulas)) {
            prediction <- rep(NA_real_, nrow(data))
            for (fold in seq_len(folds)) {
                train <- data[fold_id != fold, , drop = FALSE]
                test <- data[fold_id == fold, , drop = FALSE]
                if (nrow(test) == 0L)
                  next
                fit <- tryCatch(fit_survey_logistic(train, formulas[[model_name]]), error = function(e) NULL)
                if (!is.null(fit)) {
                  prediction[fold_id == fold] <- tryCatch(as.numeric(stats::predict(fit, newdata = test,
                    type = "response")), error = function(e) rep(NA_real_, nrow(test)))
                }
            }
            row_id <- row_id + 1L
            prediction_rows[[row_id]] <- data.frame(repeat_id = repeat_id, model = model_name, SEQN = data$SEQN,
                fold = fold_id, outcome = data$outcome, weight = data$weight, prediction = prediction,
                stringsAsFactors = FALSE)
            if (all(is.finite(prediction))) {
                metrics <- metric_set(data, prediction)
                metric_rows[[row_id]] <- data.frame(repeat_id = repeat_id, model = model_name, metric = names(metrics),
                  estimate = as.numeric(metrics), stringsAsFactors = FALSE)
            }
        }
        if (repeat_id%%5L == 0L)
            message("PSU cross-validation: ", repeat_id, "/", repeats)
    }
    list(predictions = dplyr::bind_rows(prediction_rows), metrics = dplyr::bind_rows(metric_rows),
        assignments = dplyr::bind_rows(assignment_rows))
}

model_structure_row <- function(model_name, model, data, formula) {
    coefficients <- stats::coef(model)
    parameter_count <- sum(!is.na(coefficients))
    slopes <- parameter_count - 1L
    design_df <- survey::degf(make_design(data))
    events <- sum(data$outcome == 1L)
    data.frame(model = model_name, formula = paste(deparse(formula), collapse = " "), n = nrow(data),
        events = events, prevalence_weighted = stats::weighted.mean(data$outcome, data$weight),
        coefficients_including_intercept = parameter_count, slopes = slopes, events_per_slope = events/slopes,
        design_df = design_df, residual_design_df = design_df - slopes, covariance_rank = qr(stats::vcov(model))$rank,
        converged = isTRUE(model$converged), stringsAsFactors = FALSE)
}

subgroup_long <- function(data) {
    data$age_calibration_group <- cut(data$age_years, breaks = c(20, 30, 45, 65, 80, Inf), right = FALSE,
        labels = c("20-29", "30-44", "45-64", "65-79", "80+"))
    list(age = as.character(data$age_calibration_group), sex = as.character(data$sex), race_ethnicity = as.character(data$race),
        reported_hypertension = ifelse(data$sr_hypertension == 1L, "Yes", "No"), reported_diabetes = ifelse(data$sr_diabetes ==
            1L, "Yes", "No"), reported_kidney_disease = ifelse(data$sr_ckd == 1L, "Yes", "No"),
        current_smoking = ifelse(data$current_smoking == 1L, "Yes", "No"))
}

summarize_subgroup_calibration <- function(data, cv_predictions, model_name) {
    subgroup_values <- subgroup_long(data)
    prediction_data <- dplyr::mutate(dplyr::filter(cv_predictions, .data$model == model_name), row_id = match(.data$SEQN,
        data$SEQN))
    rows <- list()
    row_id <- 0L
    for (subgroup_name in names(subgroup_values)) {
        group_vector <- subgroup_values[[subgroup_name]]
        for (level in unique(group_vector[!is.na(group_vector)])) {
            indices <- which(group_vector == level)
            observed <- stats::weighted.mean(data$outcome[indices], data$weight[indices])
            for (repeat_id in sort(unique(prediction_data$repeat_id))) {
                current <- prediction_data[prediction_data$repeat_id == repeat_id, , drop = FALSE]
                selected <- match(indices, current$row_id)
                selected <- selected[!is.na(selected)]
                predicted <- stats::weighted.mean(current$prediction[selected], current$weight[selected],
                  na.rm = TRUE)
                brier <- stats::weighted.mean((current$outcome[selected] - current$prediction[selected])^2,
                  current$weight[selected], na.rm = TRUE)
                row_id <- row_id + 1L
                rows[[row_id]] <- data.frame(subgroup = subgroup_name, level = level, repeat_id = repeat_id,
                  n = length(indices), events = sum(data$outcome[indices] == 1L), observed = observed,
                  predicted = predicted, calibration_difference = predicted - observed, observed_to_expected = ifelse(predicted >
                    0, observed/predicted, NA_real_), brier = brier, stringsAsFactors = FALSE)
            }
        }
    }
    draws <- dplyr::bind_rows(rows)
    summary <- dplyr::summarise(dplyr::group_by(draws, .data$subgroup, .data$level), n = dplyr::first(.data$n),
        events = dplyr::first(.data$events), observed = dplyr::first(.data$observed), predicted_low_across_repeats = stats::quantile(.data$predicted,
            0.025), predicted_high_across_repeats = stats::quantile(.data$predicted, 0.975), predicted = mean(.data$predicted),
        calibration_difference = mean(.data$calibration_difference), observed_to_expected = mean(.data$observed_to_expected),
        brier = mean(.data$brier), stable_cell = .data$n >= 100L & .data$events >= 10L, .groups = "drop")
    list(draws = draws, summary = summary)
}
