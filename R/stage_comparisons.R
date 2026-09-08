stage_contrasts <- function(d) {
    methods <- c(CC = "WTPH2YR_CKM", IPW = "WTPH2YR_CKM_IPW_TRUNC")
    contrasts <- list()
    for (method in names(methods)) {
        design <- design_of(d, methods[[method]])
        for (j in seq_len(nrow(specs))) {
            y <- specs$endpoint[j]
            den <- specs$den[j]
            keep <- d$ckm_stage_observed %in% 2:4 & d[[den]] & !is.na(d[[y]])
            keep[is.na(keep)] <- FALSE
            dd <- design[keep, ]
            f <- as.formula(paste0(y, " ~ stage + age10 + sex_adjust + race_adjust"))
            logfit <- survey::svyglm(f, dd, family = quasibinomial())
            prfit <- survey::svyglm(f, dd, family = quasipoisson(link = "log"))
            adjust <- survey::svyglm(as.formula(paste0(y, " ~ age10 + sex_adjust + race_adjust")),
                dd, family = quasibinomial())
            margins <- survey::svypredmeans(adjust, ~stage, predictat = factor(2:4, levels = 2:4))
            delta <- c(`2` = 0, `3` = 1, `4` = -1)
            pd <- survey::svycontrast(margins, delta)
            rdf <- logfit$df.residual
            crit <- qt(0.975, df = rdf)
            beta <- coef(prfit)
            vc <- vcov(prfit)
            z <- setNames(rep(0, length(beta)), names(beta))
            z["stage3"] <- 1
            z["stage4"] <- -1
            lp <- sum(z * beta)
            lpse <- sqrt(as.numeric(t(z) %*% vc %*% z))
            ev <- sum(dd$variables[[y]])
            non <- nrow(dd$variables) - ev
            slopes <- length(beta) - 1
            stable <- rdf >= 1 && min(ev, non) >= 10 * slopes && qr(vc)$rank == length(beta) &&
                isTRUE(logfit$converged) && isTRUE(prfit$converged)
            contrasts[[length(contrasts) + 1]] <- data.frame(method, outcome = y, n = nrow(dd$variables),
                events = ev, nonevents = non, slopes, residual_df = rdf, adjusted_p2 = coef(margins)["2"],
                adjusted_p3 = coef(margins)["3"], adjusted_p4 = coef(margins)["4"], pd34 = coef(pd)[1],
                pd34_se = survey::SE(pd)[1], pd34_low = coef(pd)[1] - crit * survey::SE(pd)[1],
                pd34_high = coef(pd)[1] + crit * survey::SE(pd)[1], pr34 = exp(lp), pr34_low = exp(lp -
                  qt(0.975, df = prfit$df.residual) * lpse), pr34_high = exp(lp + qt(0.975, df = prfit$df.residual) *
                  lpse), stable, check.names = FALSE)
        }
    }
    dplyr::bind_rows(contrasts)
}
