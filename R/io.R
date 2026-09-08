new_output_directory <- function(path) {
    if (dir.exists(path) && length(list.files(path, all.files = TRUE, no.. = TRUE))) {
        stop("Output directory is not empty: ", path)
    }
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    normalizePath(path, winslash = "/", mustWork = TRUE)
}
