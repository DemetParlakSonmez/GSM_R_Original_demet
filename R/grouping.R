#' @title Grouping-method dispatch
#' @description A single entry point for the G (Grouping) step that selects the
#'   grouping method and returns groups in the common representation used
#'   throughout gsmR: a **named list of character vectors**
#'   (`group name -> gene symbols`). This is the dispatch mechanism new grouping
#'   methods register into. The interactive app, and later the batch and CLI
#'   paths, all call this so the method logic lives in one place instead of being
#'   duplicated.
#' @name grouping
NULL

#' Compute a grouping by the selected method
#'
#' Dispatches to a grouping method and returns a named list of gene vectors,
#' identical in shape to what [read_groups()] / [read_gmt()] return, so the
#' Scoring and Modeling steps consume the result unchanged.
#'
#' @param method Grouping method. One of:
#'   * `"file"` — read a predefined grouping from a CSV/GMT file (the existing
#'     behavior). Requires `file`; `format` selects the parser.
#'   * `"bridge"` — derive groups data-drivenly with [bridge_grouping()].
#'     Requires `X` and `y`; the prior graph comes from `prior_groups`,
#'     `prior_file`, or neither (`B = 0`).
#' @param file Path to a grouping file (method `"file"`, or as the BRIDGE prior
#'   source when `prior_source = "file"`).
#' @param format Parser for `"file"`: `"auto"` (default, by extension), `"csv"`,
#'   or `"gmt"`.
#' @param X Genes-by-samples expression matrix (method `"bridge"`).
#' @param y Class-label vector/factor of length `ncol(X)` (method `"bridge"`).
#' @param prior_source For BRIDGE, where the prior graph `B` comes from:
#'   `"none"` (`B = 0`), `"groups"` (derive from a grouping named list supplied
#'   in `prior_groups`), or `"file"` (read an edge-list/grouping file at
#'   `prior_file`). Default `"none"`.
#' @param prior_groups A named grouping list to use as the BRIDGE prior when
#'   `prior_source = "groups"`.
#' @param prior_file Path to an edge-list CSV (`gene_i, gene_j, weight`) or a
#'   2-column grouping CSV to use as the BRIDGE prior when
#'   `prior_source = "file"`.
#' @param bridge_args A named list of extra arguments forwarded to
#'   [bridge_grouping()] (e.g. `lambda`, `beta`, `resolution`, `tau`,
#'   `data_type`).
#' @return A named list of character vectors (`group name -> gene symbols`).
#' @examples
#' \donttest{
#' # File-based (unchanged behavior):
#' f <- system.file("extdata", "test_groups.csv", package = "gsmR")
#' if (nzchar(f)) g <- compute_grouping("file", file = f)
#' }
#' @export
compute_grouping <- function(method = c("file", "bridge"),
                             file = NULL, format = "auto",
                             X = NULL, y = NULL,
                             prior_source = c("none", "groups", "file"),
                             prior_groups = NULL, prior_file = NULL,
                             bridge_args = list()) {
  method <- match.arg(method)

  if (method == "file") {
    if (is.null(file)) stop("method = \"file\" requires `file`.", call. = FALSE)
    fmt <- format
    if (fmt == "auto") {
      ext <- tolower(tools::file_ext(file))
      fmt <- if (ext == "gmt") "gmt" else "csv"
    }
    return(if (fmt == "gmt") read_gmt(file) else read_groups(file))
  }

  # method == "bridge"
  if (is.null(X) || is.null(y)) {
    stop("method = \"bridge\" requires `X` (genes x samples) and `y`.",
         call. = FALSE)
  }
  prior_source <- match.arg(prior_source)
  B <- switch(prior_source,
    none = NULL,
    groups = {
      if (is.null(prior_groups)) {
        stop("prior_source = \"groups\" requires `prior_groups` (a named ",
             "grouping list).", call. = FALSE)
      }
      prior_groups
    },
    file = {
      if (is.null(prior_file)) {
        stop("prior_source = \"file\" requires `prior_file`.", call. = FALSE)
      }
      bridge_read_prior_file(prior_file)
    }
  )

  do.call(bridge_grouping, c(list(X = X, y = y, B = B), bridge_args))
}

#' Read a BRIDGE prior file (edge list or 2-column grouping CSV)
#'
#' A 3-or-more-column file is treated as an edge list (`gene_i, gene_j, weight`)
#' and returned as-is for [bridge_grouping()] to convert. A 2-column file is
#' treated as a grouping file and read with [read_groups()] into a named list.
#'
#' @param file Path to the prior file.
#' @return Either a `data.frame` edge list or a named grouping list.
#' @keywords internal
bridge_read_prior_file <- function(file) {
  stopifnot(file.exists(file))
  peek <- data.table::fread(file, nrows = 1, header = TRUE, check.names = FALSE)
  if (ncol(peek) >= 3) {
    as.data.frame(
      data.table::fread(file, header = TRUE, check.names = FALSE),
      stringsAsFactors = FALSE
    )
  } else {
    read_groups(file)
  }
}
