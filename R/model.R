#' @title Final modeling using top-ranked groups
#' @description Abstract interface for building a final predictive model
#'   from the combined genes of the top-ranked groups.
#' @name model
NULL

#' Build a final model from top-ranked groups
#'
#' An abstract function that combines genes from the top-k groups into a
#' single feature set, trains a classifier, and evaluates on the test set.
#'
#' @param expr Numeric matrix (genes x samples).
#' @param labels Factor of class labels.
#' @param top_groups Named list of character vectors (top-k groups).
#' @param train_idx Integer vector of training sample indices.
#' @param method Classifier method. Default `"rf"` (Random Forest).
#' @param ... Additional arguments passed to the classifier.
#' @return A list with:
#'   \describe{
#'     \item{metrics}{Named list of classification metrics.}
#'     \item{predictions}{Factor of predicted labels for test samples.}
#'     \item{actual}{Factor of actual labels for test samples.}
#'     \item{features}{Character vector of gene names used.}
#'     \item{n_groups}{Number of groups used.}
#'   }
#' @export
build_model <- function(expr, labels, top_groups, train_idx,
                        method = "rf", ...) {
  stopifnot(is.matrix(expr) && is.numeric(expr))
  stopifnot(is.factor(labels))
  stopifnot(is.list(top_groups) && length(top_groups) > 0)

  fn <- switch(method,
    rf = build_model_rf,
    stop("Unknown modeling method: '", method, "'. Available: 'rf'.")
  )
  fn(expr, labels, top_groups, train_idx, ...)
}


#' Build a Random Forest model from top groups
#'
#' @param expr Numeric matrix (genes x samples).
#' @param labels Factor of class labels.
#' @param top_groups Named list of character vectors.
#' @param train_idx Integer vector of training sample indices.
#' @param ... Additional arguments passed to [randomForest::randomForest()].
#' @return A list with `metrics`, `predictions`, `actual`, `features`,
#'   and `n_groups`.
#' @keywords internal
build_model_rf <- function(expr, labels, top_groups, train_idx, ...) {
  if (!requireNamespace("randomForest", quietly = TRUE)) {
    stop("Package 'randomForest' is required for method 'rf'. ",
         "Install it with: install.packages('randomForest')")
  }

  test_idx <- setdiff(seq_len(ncol(expr)), train_idx)
  if (length(test_idx) == 0) {
    stop("No test samples.")
  }

  # Combine unique genes from all top groups
  all_genes <- unique(unlist(top_groups))
  # Intersect with available genes
  all_genes <- intersect(all_genes, rownames(expr))

  if (length(all_genes) == 0) {
    return(list(
      metrics     = calc_all_metrics(labels[test_idx], labels[test_idx]),
      predictions = factor(rep(NA, length(test_idx)), levels = levels(labels)),
      actual      = labels[test_idx],
      features    = character(0),
      n_groups    = length(top_groups)
    ))
  }

  # Build data frames
  train_data <- as.data.frame(t(expr[all_genes, train_idx, drop = FALSE]))
  test_data  <- as.data.frame(t(expr[all_genes, test_idx, drop = FALSE]))

  train_labels <- labels[train_idx]
  test_labels  <- labels[test_idx]

  # Determine positive class (first level or "pos" if present)
  positive_class <- if ("pos" %in% levels(labels)) "pos" else levels(labels)[1]

  # Train
  model <- randomForest::randomForest(
    x = train_data,
    y = train_labels,
    ...
  )

  # Predict
  predictions <- stats::predict(model, newdata = test_data)

  # Metrics
  metrics <- calc_all_metrics(test_labels, predictions, positive_class)

  list(
    metrics     = metrics,
    predictions = predictions,
    actual      = test_labels,
    features    = all_genes,
    n_groups    = length(top_groups)
  )
}
