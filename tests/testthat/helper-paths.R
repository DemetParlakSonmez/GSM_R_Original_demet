#' @title Test data path helper
#' @description Resolves paths to test data files. Works both when the package
#'   is installed (via system.file) and when sourced directly during development.
#' @keywords internal
test_data_path <- function(filename) {
  # Try system.file first (installed package)
  path <- system.file("extdata", filename, package = "gsmR")
  if (nchar(path) > 0 && file.exists(path)) return(path)

  # Fall back to relative paths from package root
  candidates <- c(
    file.path("inst", "extdata", filename),
    file.path("..", "..", "inst", "extdata", filename),
    file.path("../../inst/extdata", filename)
  )
  for (cand in candidates) {
    if (file.exists(cand)) return(normalizePath(cand))
  }

  stop("Test data file not found: ", filename,
       ". Run tests from the package root directory.")
}
