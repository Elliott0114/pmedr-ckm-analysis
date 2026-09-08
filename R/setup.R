options(stringsAsFactors = FALSE, survey.lonely.psu = "adjust")

required <- c("dplyr", "haven", "readr", "survey", "preventr", "tidyr")

missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing)) stop("Restore renv.lock before running: ", paste(missing, collapse = ", "))

suppressPackageStartupMessages(library(dplyr))

for (path in sort(list.files("R", pattern = "[.]R$", full.names = TRUE))) {
    if (basename(path) != "setup.R")
        source(path, encoding = "UTF-8")
}
