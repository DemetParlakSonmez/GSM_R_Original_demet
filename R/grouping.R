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

#' Derive the groups to use for one cross-validation fold
#'
#' The leakage-free entry point for cross-validation. For a fixed (file-based)
#' grouping it simply returns the dataset's pre-computed groups. For a
#' data-driven grouping carrying a `grouping_spec` (e.g. BRIDGE), it
#' **recomputes the groups from the training samples only**, so the class labels
#' are never seen for the held-out test samples of the fold.
#'
#' @param dataset A `gsm_dataset`. If `dataset$grouping_spec` is `NULL` or its
#'   method is not data-driven, `dataset$groups` is returned unchanged.
#' @param train_idx Integer vector of training-sample column indices for this
#'   fold. Only these columns of the expression matrix are used to derive groups.
#' @param expr Optional expression matrix to use instead of `dataset$expr`
#'   (e.g. a leakage-free train-fitted z-scored matrix). Defaults to
#'   `dataset$expr`.
#' @param fold_seed Optional integer seed so groups are reproducible per fold and
#'   vary across folds. Defaults to the spec's seed (or `1`).
#' @param with_graph If `TRUE`, also attach BRIDGE's differential network for
#'   this fold as a `"graph"` attribute (an edge list) on the returned groups.
#'   The `"gnn"` classifier convolves over it; every other backend ignores it.
#'   Defaults to `FALSE`, which keeps the per-fold path fast.
#' @return A named list of character vectors (`group name -> gene symbols`).
#' @export
grouping_for_fold <- function(dataset, train_idx, expr = NULL, fold_seed = NULL,
                              with_graph = FALSE) {
  spec <- dataset$grouping_spec
  if (is.null(spec) || is.null(spec$method) || spec$method != "bridge") {
    return(dataset$groups)  # fixed / file-based grouping: unchanged
  }
  X <- if (is.null(expr)) dataset$expr else expr
  Xtr <- X[, train_idx, drop = FALSE]
  ytr <- dataset$labels[train_idx]

  args <- spec$args %||% list()
  args$X <- Xtr
  args$y <- ytr
  args$B <- spec$prior
  if (!is.null(spec$prefilter))        args$prefilter        <- spec$prefilter
  if (!is.null(spec$precision_method)) args$precision_method <- spec$precision_method
  if (!is.null(spec$max_genes))        args$max_genes        <- spec$max_genes
  if (!is.null(fold_seed))             args$seed             <- fold_seed
  if (is.null(args$verbose))           args$verbose          <- FALSE
  # `details = TRUE` is what makes bridge_grouping() attach the "graph" attribute
  # (the blended differential network). Only pay for it when a caller wants it.
  # The display caps are lifted here: a GNN convolving over a truncated graph
  # would see isolated genes with self-loops only.
  if (isTRUE(with_graph)) {
    args$details         <- TRUE
    args$graph_max_nodes <- Inf
    args$graph_max_edges <- Inf
  }

  do.call(bridge_grouping, args)
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
