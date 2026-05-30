#' @title Utility functions for GSM
#' @description Small, testable helper functions for gene filtering, group
#'   management, normalization, and data splitting.
#' @name utils
NULL

#' Filter groups by size
#'
#' Remove groups with fewer than `min_size` or more than `max_size` genes.
#'
#' @param groups Named list of character vectors.
#' @param min_size Minimum group size (inclusive). Default `5`.
#' @param max_size Maximum group size (inclusive). Default `500`.
#' @return Filtered named list.
#' @export
filter_groups_by_size <- function(groups, min_size = 5, max_size = 500) {
  sizes <- vapply(groups, length, integer(1))
  keep <- sizes >= min_size & sizes <= max_size
  n_dropped <- sum(!keep)
  if (n_dropped > 0) {
    message(sprintf("Dropped %d groups outside size range [%d, %d]. %d groups remain.",
                    n_dropped, min_size, max_size, sum(keep)))
  }
  groups[keep]
}


#' Intersect gene IDs between expression data and groups
#'
#' For each group, keep only genes that are present in the expression data.
#' Groups that become empty are removed.
#'
#' @param available_genes Character vector of available gene IDs (e.g.
#'   `rownames(expr)`).
#' @param groups Named list of character vectors.
#' @return Filtered named list with only available genes.
#' @export
intersect_genes <- function(available_genes, groups) {
  result <- lapply(groups, \(g) intersect(g, available_genes))
  # Remove empty groups
  sizes <- vapply(result, length, integer(1))
  result[sizes > 0]
}


#' Check that sample IDs match between expression matrix and labels
#'
#' @param expr Numeric matrix with sample IDs as column names.
#' @param labels Named factor of class labels.
#' @return `TRUE` invisibly if matching. Stops with an error otherwise.
check_sample_match <- function(expr, labels) {
  expr_samples <- colnames(expr)
  label_samples <- names(labels)

  if (is.null(expr_samples)) {
    stop("Expression matrix has no column names (sample IDs).")
  }
  if (is.null(label_samples)) {
    stop("Labels have no names (sample IDs).")
  }
  if (length(expr_samples) != length(label_samples)) {
    stop("Number of samples in expression matrix (", length(expr_samples),
         ") does not match number of labels (", length(label_samples), ").")
  }
  if (!all(expr_samples == label_samples)) {
    stop("Sample IDs in expression matrix do not match label names.")
  }
  invisible(TRUE)
}


#' Normalize expression data
#'
#' Apply normalization to each gene (row) across samples. This is an abstract
#' function that dispatches to the specified method.
#'
#' @param expr Numeric matrix (genes x samples).
#' @param method Normalization method. Currently `"zscore"` (default) or
#'   `"none"`.
#' @return Normalized numeric matrix with the same dimensions.
#' @export
normalize_expression <- function(expr, method = "zscore") {
  fn <- switch(method,
    zscore = normalize_zscore,
    none   = identity,
    stop("Unknown normalization method: ", method)
  )
  fn(expr)
}


#' Z-score normalization (per gene)
#' @param expr Numeric matrix (genes x samples).
#' @return Z-score normalized matrix.
#' @keywords internal
normalize_zscore <- function(expr) {
  row_means <- rowMeans(expr, na.rm = TRUE)
  row_sds <- apply(expr, 1, stats::sd, na.rm = TRUE)

  # Avoid division by zero for constant genes
  row_sds[row_sds == 0] <- 1

  (expr - row_means) / row_sds
}


#' Split samples into training and test sets
#'
#' Creates a random train/test split while preserving class proportions
#' (stratified sampling).
#'
#' @param labels Factor of class labels.
#' @param ratio Proportion of samples to use for training. Default `0.7`.
#' @param seed Random seed for reproducibility. Default `NULL` (no seed set).
#' @return A list with `$train` and `$test` integer index vectors.
#' @export
split_train_test <- function(labels, ratio = 0.7, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  n <- length(labels)
  train_idx <- integer(0)

  # Stratified split
  for (lev in levels(labels)) {
    idx <- which(labels == lev)
    n_train <- max(1, round(length(idx) * ratio))
    sampled <- sample(idx, n_train)
    train_idx <- c(train_idx, sampled)
  }

  train_idx <- sort(train_idx)
  test_idx <- setdiff(seq_len(n), train_idx)

  list(train = train_idx, test = test_idx)
}
