#' Launch the gsmR Shiny Application
#'
#' Opens an interactive web interface for the G-S-M workflow. The app provides
#' a GUI for loading data, collapsing replicates, pre-scoring groups, and
#' scoring groups with machine learning classifiers.
#'
#' @param ... Additional arguments passed to [shiny::runApp()].
#' @return This function does not return; it runs the Shiny app until the
#'   browser window is closed.
#' @export
#'
#' @examples
#' if (interactive()) {
#'   run_app()
#' }
run_app <- function(...) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("Package 'shiny' is required. Install with: install.packages('shiny')")
  }
  if (!requireNamespace("bslib", quietly = TRUE)) {
    stop("Package 'bslib' is required. Install with: install.packages('bslib')")
  }
  if (!requireNamespace("DT", quietly = TRUE)) {
    stop("Package 'DT' is required. Install with: install.packages('DT')")
  }

  app_dir <- system.file("shiny", package = "gsmR")
  if (app_dir == "") {
    # Dev mode fallback
    app_dir <- file.path(getwd(), "inst", "shiny")
  }
  if (!dir.exists(app_dir)) {
    stop("Cannot find the Shiny app directory. ",
         "Make sure gsmR is installed or run from the package root.")
  }

  shiny::runApp(app_dir, ...)
}
