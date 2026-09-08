source("R/setup.R", encoding = "UTF-8")

args <- commandArgs(trailingOnly = TRUE)

stopifnot(length(args) == 2L)

raw_dir <- normalizePath(args[1], winslash = "/", mustWork = TRUE)

out <- new_output_directory(args[2])

tables <- c("DEMO_L", "BMX_L", "BPXO_L", "BPQ_L", "DIQ_L", "GHB_L", "TCHOL_L", "HDL_L", "BIOPRO_L",
    "ALB_CR_L", "MCQ_L", "KIQ_U_L", "SMQ_L", "HIQ_L", "HUQ_L", "FSQ_L")

data <- setNames(lapply(tables, function(n) {
    z <- haven::zap_labels(haven::read_xpt(file.path(raw_dir, paste0(n, ".xpt"))))
    stopifnot(!anyDuplicated(z$SEQN))
    z
}), tables)

stopifnot("WTPH2YR" %in% names(data$GHB_L))

names(data$GHB_L)[names(data$GHB_L) == "WTPH2YR"] <- "WTPH2YR_CKM"

for (n in c("TCHOL_L", "HDL_L", "BIOPRO_L")) data[[n]]$WTPH2YR <- NULL

merged <- Reduce(function(x, y) dplyr::full_join(x, y, by = "SEQN"), data)

adult <- merged[!is.na(merged$RIDAGEYR) & merged$RIDAGEYR >= 20, ]

base <- filter(derive_ckm_variables(adult), is.finite(WTPH2YR_CKM), WTPH2YR_CKM > 0, !is.na(SDMVSTRA),
    !is.na(SDMVPSU))

base <- national_selection_weights(base)

d <- base[base$complete_ckm, ]

risk <- risks(d)

d$prevent_total_cvd_10yr <- risk$risk

d$prevent_input_problem <- risk$problem

d$stage3_prevent_high_risk <- d$prevent_total_cvd_10yr >= 0.2

d$ckm_stage_observed <- assign_stage(d, risk$risk)

d <- clinical_indicators(d)

stopifnot(nrow(d) == 5006L, sum(d$ckm_stage_observed %in% 3) == 294L, sum(d$ckm_stage_observed %in%
    4) == 601L, sum(is.na(d$ckm_stage_observed)) == 4L)

saveRDS(d, file.path(out, "national_data.rds"))

readr::write_csv(data.frame(step = c("Survey records", "Adults aged 20 years or older", "Positive phlebotomy weight and survey design",
    "Complete clinical inputs", "Stage unresolved"), n = c(nrow(merged), nrow(adult), nrow(base),
    nrow(d), sum(is.na(d$ckm_stage_observed)))), file.path(out, "national_sample_flow.csv"))

capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
