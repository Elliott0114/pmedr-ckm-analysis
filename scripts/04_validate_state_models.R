source("R/setup.R", encoding = "UTF-8")

args <- commandArgs(trailingOnly = TRUE)

stopifnot(length(args) == 2L)

train <- readRDS(file.path(args[1], "training_data.rds"))

models <- readRDS(file.path(args[1], "models.rds"))

formulas <- state_formulas()

out <- new_output_directory(args[2])

write_table <- function(x, n) readr::write_csv(x, file.path(out, n), na = "NA")

cv <- repeated_psu_cv_predictions(train, formulas, 20, 5, 20260907)

stopifnot(nrow(cv$metrics) == 20 * length(models) * 8, all(is.finite(cv$predictions$prediction)))

write_table(cv$metrics, "cv_metrics_draws.csv")

write_table(cv$assignments, "cv_psu_assignments.csv")

cv_summary <- summarise(group_by(cv$metrics, model, metric), partition_low = quantile(estimate,
    0.025), partition_high = quantile(estimate, 0.975), estimate = mean(estimate), .groups = "drop")

write_table(cv_summary, "cv_metrics.csv")

app <- metric_set(train, predict(models$lean_core, newdata = train, type = "response"))

boot <- bootstrap_optimism(train, formulas$lean_core, app, 200, 20260907)

stopifnot(length(unique(boot$replicate)) >= 190)

write_table(boot, "bootstrap_optimism_draws.csv")

boot_summary <- summarise(group_by(boot, metric), estimate = mean(corrected), mean_optimism = mean(optimism),
    successful_replicates = n(), .groups = "drop")

write_table(boot_summary, "bootstrap_optimism_summary.csv")

sub <- summarize_subgroup_calibration(train, cv$predictions, "lean_core")

write_table(sub$summary, "subgroup_calibration.csv")

oof <- summarise(group_by(filter(cv$predictions, model == "lean_core"), SEQN), pred = mean(prediction),
    .groups = "drop")

cal <- train

cal$pred <- oof$pred[match(train$SEQN, oof$SEQN)]

cal$bin <- weighted_quantile_group(cal$pred, cal$weight)

bins <- as.data.frame(survey::svyby(~outcome + pred, ~bin, make_design(cal), survey::svymean, vartype = c("se",
    "ci")))

counts <- summarise(group_by(cal, bin), n = n(), events = sum(outcome), prediction_min = min(pred),
    prediction_max = max(pred), .groups = "drop")

bins <- left_join(bins, counts, by = "bin")

bins$interval_note <- ifelse(bins$events == 0 | bins$events == bins$n, "Boundary observed proportion: interval not estimable",
    "Normal-Wald uncertainty in observed proportion only")

bins$ci_l.outcome[bins$events == 0 | bins$events == bins$n] <- NA_real_

bins$ci_u.outcome[bins$events == 0 | bins$events == bins$n] <- NA_real_

write_table(bins, "calibration_with_counts.csv")

capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
