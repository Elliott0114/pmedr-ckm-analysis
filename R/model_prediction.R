# Use the training spline knots in every prediction call.

predict.pmedr_svyglm <- function(object, newdata = NULL, type = c("link", "response"), ...) {
    type <- match.arg(type)
    if (is.null(newdata))
        newdata <- object$survey.design$variables
    tt <- stats::delete.response(stats::terms(object))
    mf <- stats::model.frame(tt, data = newdata, xlev = object$xlevels, na.action = stats::na.pass)
    mm <- stats::model.matrix(tt, mf, contrasts.arg = object$contrasts)
    beta <- stats::coef(object)
    stopifnot(all(names(beta) %in% colnames(mm)), all(is.finite(beta)))
    eta <- as.numeric(mm[, names(beta), drop = FALSE] %*% beta)
    if (type == "response")
        object$family$linkinv(eta)
    else eta
}

fit_survey_logistic <- function(data, formula) {
    fit <- survey::svyglm(formula, make_design(data), family = quasibinomial())
    class(fit) <- c("pmedr_svyglm", class(fit))
    fit
}

verify_fixed_prediction <- function(model, newdata) {
    p <- as.numeric(stats::predict(model, newdata = newdata, type = "response"))
    halves <- split(seq_len(nrow(newdata)), rep(1:2, length.out = nrow(newdata)))
    partitioned <- numeric(nrow(newdata))
    for (ix in halves) partitioned[ix] <- stats::predict(model, newdata = newdata[ix, ], type = "response")
    mm <- model_matrix_newdata(model, newdata)
    explicit <- as.numeric(plogis(mm %*% coef(model)))
    stopifnot(max(abs(p - partitioned)) < 1e-12, max(abs(p - explicit)) < 1e-12)
    data.frame(max_partition_difference = max(abs(p - partitioned)), max_matrix_difference = max(abs(p -
        explicit)))
}
