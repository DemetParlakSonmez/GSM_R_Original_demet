#' @title GSM Dataset construction and validation
#' @description Construct a validated GSM dataset from an expression matrix,
#'   class labels, and grouping information. The dataset bundles all inputs
#'   and intersects gene IDs so that groups only reference genes present in
#'   the expression matrix.
#' @name dataset
NULL

#' Construct a GSM dataset
#'
#' Validates and bundles an expression matrix, class labels, and gene groups
#' into a single S3 object of class `gsm_dataset`.
#'
#' @param expr Numeric matrix with genes as rows and samples as columns.
#' @param labels Named factor of class labels. Names must match `colnames(expr)`.
#' @param groups Named list of character vectors (group name -> gene symbols).
#' @param min_group_size Minimum number of genes in a group after intersection.
#'   Groups below this threshold are dropped. Default `5`.
#' @param max_group_size Maximum number of genes in a group. Groups above this
#'   threshold are dropped. Default `500`.
#' @return An S3 object of class `gsm_dataset`.
#' @export
gsm_dataset <- function(expr, labels, groups,
                        min_group_size = 5, max_group_size = 500) {
  # --- Validate inputs ---
  if (!is.matrix(expr) || !is.numeric(expr)) {
    stop("`expr` must be a numeric matrix (genes x samples).")
  }
  if (is.null(rownames(expr))) {
    stop("`expr` must have row names (gene identifiers).")
  }
  if (is.null(colnames(expr))) {
    stop("`expr` must have column names (sample identifiers).")
  }
  if (!is.factor(labels)) {
    labels <- factor(labels)
  }

  # Check sample matching
  check_sample_match(expr, labels)

  if (!is.list(groups) || is.null(names(groups))) {
    stop("`groups` must be a named list of character vectors.")
  }

  # --- Intersect genes ---
  available_genes <- rownames(expr)
  filtered_groups <- intersect_genes(available_genes, groups)

  # Count genes lost
  n_original_genes <- length(unique(unlist(groups)))
  n_available_genes <- length(unique(unlist(filtered_groups)))
  n_lost <- n_original_genes - n_available_genes
  if (n_lost > 0) {
    message(sprintf(
      "Note: %d/%d genes in groups are not in the expression matrix and were removed.",
      n_lost, n_original_genes
    ))
  }

  # --- Filter by size ---
  filtered_groups <- filter_groups_by_size(
    filtered_groups,
    min_size = min_group_size,
    max_size = max_group_size
  )

  if (length(filtered_groups) == 0) {
    stop("No groups remain after filtering. Check gene ID matching between ",
         "expression matrix and grouping file.")
  }

  # --- Build object ---
  obj <- list(
    expr   = expr,
    labels = labels,
    groups = filtered_groups,
    params = list(
      min_group_size = min_group_size,
      max_group_size = max_group_size,
      n_samples      = ncol(expr),
      n_genes        = nrow(expr),
      n_groups       = length(filtered_groups),
      n_classes      = length(levels(labels))
    )
  )
  class(obj) <- "gsm_dataset"
  obj
}


#' Print a GSM dataset
#' @param x A `gsm_dataset` object.
#' @param ... Ignored.
#' @export
print.gsm_dataset <- function(x, ...) {
  cat("GSM Dataset\n")
  cat("  Samples:", x$params$n_samples, "\n")
  cat("  Genes:  ", x$params$n_genes, "\n")
  cat("  Groups: ", x$params$n_groups, "\n")
  cat("  Classes:", paste(levels(x$labels), collapse = ", "), "\n")
  cat("  Class distribution:", paste(
    paste0(levels(x$labels), "=", table(x$labels)),
    collapse = ", "
  ), "\n")
  group_sizes <- vapply(x$groups, length, integer(1))
  cat("  Group sizes: min=", min(group_sizes),
      ", median=", stats::median(group_sizes),
      ", max=", max(group_sizes), "\n")
  invisible(x)
}


#' Summarise a GSM dataset
#' @param object A `gsm_dataset` object.
#' @param ... Ignored.
#' @export
summary.gsm_dataset <- function(object, ...) {
  group_sizes <- vapply(object$groups, length, integer(1))
  list(
    n_samples  = object$params$n_samples,
    n_genes    = object$params$n_genes,
    n_groups   = object$params$n_groups,
    classes    = levels(object$labels),
    class_dist = as.list(table(object$labels)),
    group_size_summary = summary(group_sizes)
  )
}


#' Validate a GSM dataset object
#'
#' Checks that the object has the correct structure and consistent dimensions.
#'
#' @param x An object to validate.
#' @return `TRUE` invisibly if valid. Stops with an error otherwise.
#' @export
validate_gsm_dataset <- function(x) {
  if (!inherits(x, "gsm_dataset")) {
    stop("Object is not a gsm_dataset.")
  }
  if (!is.matrix(x$expr) || !is.numeric(x$expr)) {
    stop("$expr is not a numeric matrix.")
  }
  if (!is.factor(x$labels)) {
    stop("$labels is not a factor.")
  }
  if (ncol(x$expr) != length(x$labels)) {
    stop("Number of samples in $expr (", ncol(x$expr),
         ") does not match length of $labels (", length(x$labels), ").")
  }
  if (!all(colnames(x$expr) == names(x$labels))) {
    stop("Sample names in $expr do not match names in $labels.")
  }
  if (!is.list(x$groups) || length(x$groups) == 0) {
    stop("$groups must be a non-empty named list.")
  }
  # Check all group genes are in the expression matrix
  all_group_genes <- unique(unlist(x$groups))
  missing <- setdiff(all_group_genes, rownames(x$expr))
  if (length(missing) > 0) {
    stop(length(missing), " genes in groups are not in expression matrix: ",
         paste(head(missing, 5), collapse = ", "), "...")
  }
  invisible(TRUE)
}
