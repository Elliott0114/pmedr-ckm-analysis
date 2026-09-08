# The modeled stage-3 contribution is (1 - reported CVD) times conditional
# predicted stage-3 probability. Survey sampling and coefficient variance
# are combined using a first-order approximation.

weighted_mean_safe <- function(x, w) {
    keep <- !is.na(x) & is.finite(w) & w > 0
    if (!any(keep))
        return(NA_real_)
    stats::weighted.mean(x[keep], w[keep])
}

weighted_var_safe <- function(x, w) {
    keep <- !is.na(x) & is.finite(w) & w > 0
    if (sum(keep) < 2L)
        return(NA_real_)
    x <- as.numeric(x[keep])
    w <- as.numeric(w[keep])
    mu <- sum(w * x)/sum(w)
    sum(w * (x - mu)^2)/sum(w)
}

kish_ess <- function(w) {
    w <- w[is.finite(w) & w > 0]
    sum(w)^2/sum(w^2)
}

fit_response_model <- function(brfss, lean_predictors) {
    brfss$analysis_complete <- as.integer(stats::complete.cases(brfss[, lean_predictors, drop = FALSE]))
    brfss$sex_response <- addNA(factor(brfss$sex), ifany = TRUE)
    brfss$race_response <- addNA(factor(brfss$race), ifany = TRUE)
    brfss$state_response <- factor(brfss$state_abbr)
    brfss$response_weight_scaled <- brfss$brfss_weight/mean(brfss$brfss_weight)
    response_model <- suppressWarnings(stats::glm(analysis_complete ~ splines::ns(age_years, df = 3) +
        sex_response + race_response + state_response, data = brfss, family = quasibinomial(), weights = response_weight_scaled,
        na.action = stats::na.exclude))
    brfss$response_probability_raw <- as.numeric(stats::predict(response_model, newdata = brfss,
        type = "response"))
    complete_probability <- brfss$response_probability_raw[brfss$analysis_complete == 1L]
    probability_limits <- stats::quantile(complete_probability, c(0.01, 0.99), na.rm = TRUE)
    lower <- max(0.05, probability_limits[[1]])
    upper <- min(0.995, probability_limits[[2]])
    brfss$response_probability_truncated <- pmin(pmax(brfss$response_probability_raw, lower), upper)
    list(data = brfss, model = response_model, lower = lower, upper = upper)
}

state_response_diagnostics <- function(response_data) {
    dplyr::bind_rows(lapply(split(response_data, response_data$state_abbr), function(d) {
        complete <- d$analysis_complete == 1L
        raw_weight <- d$brfss_weight[complete]/pmax(d$response_probability_raw[complete], 0.01)
        trunc_weight <- d$brfss_weight[complete]/d$response_probability_truncated[complete]
        data.frame(state_abbr = as.character(d$state_abbr[1]), state_name = as.character(d$state_name[1]),
            n_no_cvd = nrow(d), n_complete = sum(complete), weighted_complete_fraction = weighted_mean_safe(d$analysis_complete,
                d$brfss_weight), response_probability_min_complete = min(d$response_probability_raw[complete]),
            response_probability_p01_complete = unname(stats::quantile(d$response_probability_raw[complete],
                0.01)), response_probability_median_complete = stats::median(d$response_probability_raw[complete]),
            response_probability_p99_complete = unname(stats::quantile(d$response_probability_raw[complete],
                0.99)), response_probability_max_complete = max(d$response_probability_raw[complete]),
            original_weight_ess_complete = kish_ess(d$brfss_weight[complete]), raw_ipw_ess = kish_ess(raw_weight),
            truncated_ipw_ess = kish_ess(trunc_weight), truncated_ipw_max_share = max(trunc_weight)/sum(trunc_weight),
            stringsAsFactors = FALSE)
    }))
}

smd_row <- function(variable, level, nhanes_value, brfss_value, type = c("binary", "continuous"),
    nhanes_sd = NA, brfss_sd = NA) {
    type <- match.arg(type)
    if (type == "continuous") {
        denominator <- sqrt((nhanes_sd^2 + brfss_sd^2)/2)
    }
    else {
        denominator <- sqrt((nhanes_value * (1 - nhanes_value) + brfss_value * (1 - brfss_value))/2)
    }
    data.frame(variable = variable, level = level, nhanes = nhanes_value, brfss = brfss_value, standardized_mean_difference = (brfss_value -
        nhanes_value)/denominator, abs_smd = abs((brfss_value - nhanes_value)/denominator), stringsAsFactors = FALSE)
}

transport_overlap <- function(nhanes, brfss) {
    rows <- list()
    age_n <- weighted_mean_safe(nhanes$age_years, nhanes$weight)
    age_b <- weighted_mean_safe(brfss$age_years, brfss$brfss_weight)
    rows[[1]] <- smd_row("Age", "continuous", age_n, age_b, "continuous", sqrt(weighted_var_safe(nhanes$age_years,
        nhanes$weight)), sqrt(weighted_var_safe(brfss$age_years, brfss$brfss_weight)))
    variables <- c("prevent_age_eligible", "sex", "race", "sr_hypertension", "sr_diabetes", "sr_ckd",
        "current_smoking")
    row_id <- 1L
    for (variable in variables) {
        levels <- if (is.factor(nhanes[[variable]]) || is.factor(brfss[[variable]])) {
            union(as.character(unique(nhanes[[variable]])), as.character(unique(brfss[[variable]])))
        }
        else {
            sort(unique(c(nhanes[[variable]], brfss[[variable]])))
        }
        levels <- levels[!is.na(levels)]
        for (level in levels) {
            p_n <- weighted_mean_safe(as.integer(as.character(nhanes[[variable]]) == as.character(level)),
                nhanes$weight)
            p_b <- weighted_mean_safe(as.integer(as.character(brfss[[variable]]) == as.character(level)),
                brfss$brfss_weight)
            row_id <- row_id + 1L
            rows[[row_id]] <- smd_row(variable, as.character(level), p_n, p_b, "binary")
        }
    }
    smd <- dplyr::arrange(dplyr::bind_rows(rows), dplyr::desc(.data$abs_smd))
    nhanes_source <- dplyr::transmute(nhanes, source = 0L, age_years, prevent_age_eligible, sex,
        race, sr_hypertension, sr_diabetes, sr_ckd, current_smoking, source_weight = weight/sum(weight))
    brfss_source <- dplyr::transmute(brfss, source = 1L, age_years, prevent_age_eligible, sex, race,
        sr_hypertension, sr_diabetes, sr_ckd, current_smoking, source_weight = brfss_weight/sum(brfss_weight))
    combined <- dplyr::bind_rows(nhanes_source, brfss_source)
    combined$source_weight <- combined$source_weight * nrow(combined)/2
    source_model <- stats::glm(source ~ splines::ns(age_years, df = 2) + prevent_age_eligible +
        sex + race + sr_hypertension + sr_diabetes + sr_ckd + current_smoking, data = combined,
        family = quasibinomial(), weights = source_weight)
    combined$source_probability <- clip_probability(stats::predict(source_model, type = "response"))
    source_auc <- weighted_auc(combined$source, combined$source_probability, combined$source_weight)
    breaks <- seq(0, 1, length.out = 101)
    combined$bin <- cut(combined$source_probability, breaks = breaks, include.lowest = TRUE, labels = FALSE)
    masses <- tidyr::pivot_wider(tidyr::complete(dplyr::ungroup(dplyr::mutate(dplyr::group_by(dplyr::summarise(dplyr::group_by(combined,
        .data$source, .data$bin), mass = sum(.data$source_weight), .groups = "drop"), .data$source),
        mass = .data$mass/sum(.data$mass))), source, bin = 1:100, fill = list(mass = 0)), names_from = source,
        values_from = mass, names_prefix = "source_")
    overlap_coefficient <- sum(pmin(masses$source_0, masses$source_1))
    pattern_vars <- c("sex", "race", "sr_hypertension", "sr_diabetes", "sr_ckd", "current_smoking")
    make_pattern <- function(data) do.call(paste, c(lapply(data[, pattern_vars, drop = FALSE], as.character),
        sep = "|"))
    nhanes$pattern <- make_pattern(nhanes)
    brfss$pattern <- make_pattern(brfss)
    support <- dplyr::summarise(dplyr::group_by(nhanes, .data$pattern), age_min = min(.data$age_years),
        age_max = max(.data$age_years), nhanes_n = dplyr::n(), .groups = "drop")
    brfss <- dplyr::mutate(dplyr::left_join(brfss, support, by = "pattern"), exact_pattern_supported = !is.na(.data$nhanes_n),
        pattern_and_age_supported = .data$exact_pattern_supported & .data$age_years >= .data$age_min &
            .data$age_years <= .data$age_max)
    state_support <- dplyr::summarise(dplyr::group_by(brfss, .data$state_abbr, .data$state_name),
        n = dplyr::n(), weighted_exact_pattern_support = stats::weighted.mean(.data$exact_pattern_supported,
            .data$brfss_weight), weighted_pattern_age_support = stats::weighted.mean(.data$pattern_and_age_supported,
            .data$brfss_weight), .groups = "drop")
    metrics <- data.frame(metric = c("Source classifier weighted AUC", "Propensity histogram overlap coefficient",
        "Maximum absolute weighted SMD", "Median state exact-pattern support", "Minimum state exact-pattern support",
        "Median state pattern-and-age support", "Minimum state pattern-and-age support"), value = c(source_auc,
        overlap_coefficient, max(smd$abs_smd, na.rm = TRUE), stats::median(state_support$weighted_exact_pattern_support),
        min(state_support$weighted_exact_pattern_support), stats::median(state_support$weighted_pattern_age_support),
        min(state_support$weighted_pattern_age_support)), stringsAsFactors = FALSE)
    list(smd = smd, metrics = metrics, state_support = state_support, source_model = source_model)
}

model_matrix_newdata <- function(model, newdata) {
    matrix <- stats::model.matrix(stats::delete.response(stats::terms(model)), data = newdata, contrasts.arg = model$contrasts,
        xlev = model$xlevels)
    matrix[, names(stats::coef(model)), drop = FALSE]
}

estimate_states <- function(model, method, model_name, br, cc, ref) {
    p <- as.numeric(predict(model, newdata = cc, type = "response"))
    X <- model_matrix_newdata(model, cc)
    pi <- if (method == "CC")
        rep(1, nrow(cc))
    else cc$response_probability_truncated
    full <- br
    full$r <- 0
    full$s <- 0
    full$r[cc$row_id] <- 1/pi
    full$s[cc$row_id] <- p/pi
    full$one <- 1
    full$c <- full$sr_cvd
    ans <- list()
    gradients <- list()
    for (state in sort(unique(full$state_abbr))) {
        sf <- full[full$state_abbr == state, ]
        ix <- which(cc$state_abbr == state)
        sf$psu_key <- interaction(sf$brfss_strata, sf$brfss_psu, drop = TRUE)
        vars <- character()
        for (g in 0:4) {
            mask <- if (g == 0)
                rep(1, nrow(sf))
            else as.integer(sf$age_std == names(ref)[g])
            for (v in c("one", "c", "r", "s")) {
                nm <- paste0("g", g, "_", v)
                sf[[nm]] <- mask * sf[[v]]
                vars <- c(vars, nm)
            }
        }
        des <- survey::svydesign(ids = ~psu_key, strata = ~brfss_strata, weights = ~brfss_weight,
            data = sf, nest = TRUE)
        tot <- survey::svytotal(reformulate(vars), des)
        tv <- coef(tot)
        V <- vcov(tot)
        G <- matrix(0, 4, length(tv), dimnames = list(c("cvd", "stage3", "cvd_std", "stage3_std"),
            names(tv)))
        GB <- matrix(0, 4, ncol(X), dimnames = list(rownames(G), colnames(X)))
        vals <- setNames(rep(0, 4), rownames(G))
        for (g in 0:4) {
            nms <- paste0("g", g, "_", c("one", "c", "r", "s"))
            t <- unname(tv[nms])
            T <- t[1]
            C <- t[2]
            R <- t[3]
            S <- t[4]
            stopifnot(T > 0, R > 0)
            a <- 1 - C/T
            q <- S/R
            cv <- C/T
            st3 <- a * q
            dc <- c(-C/T^2, 1/T, 0, 0)
            ds <- c(C/T^2 * q, -q/T, -a * S/R^2, a/R)
            rows <- if (g == 0)
                1:2
            else 3:4
            alpha <- if (g == 0)
                1
            else ref[g]
            vals[rows] <- vals[rows] + alpha * c(cv, st3)
            G[rows[1], nms] <- G[rows[1], nms] + alpha * dc
            G[rows[2], nms] <- G[rows[2], nms] + alpha * ds
            ii <- if (g == 0)
                ix
            else ix[cc$age_std[ix] == names(ref)[g]]
            gb <- a * colSums(X[ii, , drop = FALSE] * (cc$brfss_weight[ii]/pi[ii] * p[ii] * (1 -
                p[ii])))/R
            GB[rows[2], ] <- GB[rows[2], ] + alpha * gb
        }
        sv <- G %*% V %*% t(G)
        mv <- GB %*% vcov(model) %*% t(GB)
        joint <- sv + mv
        independent_point <- (1 - weighted.mean(sf$sr_cvd, sf$brfss_weight)) * weighted.mean(p[ix],
            cc$brfss_weight[ix]/pi[ix])
        if (abs(vals["stage3"] - independent_point) > 1e-10)
            stop("State ratio identity failed: ", state, " / ", method, "; totals=", vals["stage3"],
                "; direct=", independent_point)
        se <- sqrt(diag(joint))
        for (k in seq_along(vals)) ans[[length(ans) + 1]] <- data.frame(state_abbr = state, state_name = sf$state_name[1],
            model = model_name, method, estimand = names(vals)[k], estimate = vals[k], sampling_se = sqrt(sv[k,
                k]), model_se = sqrt(mv[k, k]), se = se[k], low = vals[k] - qnorm(0.975) * se[k],
            high = vals[k] + qnorm(0.975) * se[k], n = nrow(sf), design_df = survey::degf(des))
        gradients[[state]] <- GB
    }
    list(estimates = bind_rows(ans), gradients = gradients, predictions = p)
}
