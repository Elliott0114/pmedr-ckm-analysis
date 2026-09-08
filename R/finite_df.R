finite_df_intervals <- function(a, model) {
    vs <- a$sampling_se^2
    vm <- a$model_se^2
    a$nhanes_design_df <- survey::degf(model$survey.design)
    a$nhanes_residual_df <- model$df.residual
    for (nu in c(a$nhanes_design_df[1], a$nhanes_residual_df[1])) {
        df <- (vs + vm)^2/(vs^2/a$design_df + vm^2/nu)
        crit <- qt(0.975, df)
        a[[paste0("effective_df_", nu)]] <- df
        a[[paste0("critical_", nu)]] <- crit
        a[[paste0("low_ws_", nu)]] <- a$estimate - crit * a$se
        a[[paste0("high_ws_", nu)]] <- a$estimate + crit * a$se
    }
    a$model_variance_share <- vm/(vs + vm)
    a
}
