#' @title Group scoring using machine learning classifiers
#' @description Abstract interface for scoring individual gene groups by
#'   training a classifier on the group's genes and evaluating accuracy.
#' @name score
NULL

#' Score a single gene group
#'
#' An abstract function that trains a classifier on the genes of a single
#' group using training data and evaluates on test data. Dispatches to a
#' specific classifier method.
#'
#' @param expr Numeric matrix (genes x samples).
#' @param labels Factor of class labels.
#' @param genes Character vector of gene names for this group.
#' @param train_idx Integer vector of training sample indices.
#' @param method Classifier method. Default `"rf"` (Random Forest).
#' @param ... Additional arguments passed to the classifier.
#' @return A list with:
#'   \describe{
#'     \item{accuracy}{Classification accuracy on the test set.}
#'     \item{predictions}{Factor of predicted labels for test samples.}
#'     \item{actual}{Factor of actual labels for test samples.}
#'   }
#' @export
score_group <- function(expr, labels, genes, train_idx,
                        method = "rf", ...) {
  stopifnot(is.matrix(expr) && is.numeric(expr))
  stopifnot(is.factor(labels))
  stopifnot(is.character(genes) && length(genes) > 0)

  # Filter to genes present in expression matrix
  genes <- intersect(genes, rownames(expr))
  if (length(genes) == 0) {
    return(list(accuracy = NA_real_, predictions = NULL, actual = NULL))
  }

  fn <- switch(method,
    rf  = score_group_rf,
    svm = score_group_svm,
    knn = score_group_knn,
    stop("Unknown scoring method: '", method, "'. Available: 'rf', 'svm', 'knn'.")
  )
  fn(expr, labels, genes, train_idx, ...)
}


#' Score a group using Random Forest
#'
#' Trains a Random Forest classifier on the group's genes from training
#' samples and evaluates accuracy on test samples.
#'
#' @param expr Numeric matrix (genes x samples).
#' @param labels Factor of class labels.
#' @param genes Character vector of gene names.
#' @param train_idx Integer vector of training sample indices.
#' @param ... Additional arguments passed to [randomForest::randomForest()].
#' @return A list with `accuracy`, `predictions`, and `actual`.
#' @keywords internal
score_group_rf <- function(expr, labels, genes, train_idx, ...) {
  if (!requireNamespace("randomForest", quietly = TRUE)) {
    stop("Package 'randomForest' is required for method 'rf'. ",
         "Install it with: install.packages('randomForest')")
  }

  test_idx <- setdiff(seq_len(ncol(expr)), train_idx)
  if (length(test_idx) == 0) {
    stop("No test samples. Ensure train_idx does not include all samples.")
  }

  # Build data frames (samples x genes)
  train_data <- as.data.frame(t(expr[genes, train_idx, drop = FALSE]))
  test_data  <- as.data.frame(t(expr[genes, test_idx, drop = FALSE]))

  train_labels <- labels[train_idx]
  test_labels  <- labels[test_idx]

  # Train Random Forest
  model <- randomForest::randomForest(
    x = train_data,
    y = train_labels,
    ...
  )

  # Predict on test set
  predictions <- stats::predict(model, newdata = test_data)

  # Calculate accuracy
  accuracy <- calc_accuracy(test_labels, predictions)

  list(
    accuracy    = accuracy,
    predictions = predictions,
    actual      = test_labels
  )
}


#' Score a group using SVM
#'
#' @param expr Numeric matrix (genes x samples).
#' @param labels Factor of class labels.
#' @param genes Character vector of gene names.
#' @param train_idx Integer vector of training sample indices.
#' @param ... Additional arguments passed to [e1071::svm()].
#' @return A list with `accuracy`, `predictions`, and `actual`.
#' @keywords internal
score_group_svm <- function(expr, labels, genes, train_idx, ...) {
  if (!requireNamespace("e1071", quietly = TRUE)) {
    stop("Package 'e1071' is required for method 'svm'. ",
         "Install it with: install.packages('e1071')")
  }

  test_idx <- setdiff(seq_len(ncol(expr)), train_idx)
  if (length(test_idx) == 0) {
    stop("No test samples.")
  }

  train_data <- as.data.frame(t(expr[genes, train_idx, drop = FALSE]))
  test_data  <- as.data.frame(t(expr[genes, test_idx, drop = FALSE]))

  train_labels <- labels[train_idx]
  test_labels  <- labels[test_idx]

  # Train SVM
  model <- e1071::svm(
    x = train_data,
    y = train_labels,
    ...
  )

  # Predict on test set
  predictions <- stats::predict(model, newdata = test_data)

  # Calculate accuracy
  accuracy <- calc_accuracy(test_labels, predictions)

  list(
    accuracy    = accuracy,
    predictions = predictions,
    actual      = test_labels
  )
}


#' Score a group using KNN
#'
#' @param expr Numeric matrix (genes x samples).
#' @param labels Factor of class labels.
#' @param genes Character vector of gene names.
#' @param train_idx Integer vector of training sample indices.
#' @param ... Additional arguments passed to [class::knn()].
#' @return A list with `accuracy`, `predictions`, and `actual`.
#' @keywords internal
score_group_knn <- function(expr, labels, genes, train_idx, ...) {
  if (!requireNamespace("class", quietly = TRUE)) {
    stop("Package 'class' is required for method 'knn'.")
  }

  test_idx <- setdiff(seq_len(ncol(expr)), train_idx)
  if (length(test_idx) == 0) {
    stop("No test samples.")
  }

  train_data <- t(expr[genes, train_idx, drop = FALSE])
  test_data  <- t(expr[genes, test_idx, drop = FALSE])

  train_labels <- labels[train_idx]
  test_labels  <- labels[test_idx]

  # Predict using KNN (class::knn)
  predictions <- class::knn(
    train = train_data,
    test  = test_data,
    cl    = train_labels,
    ...
  )

  # Calculate accuracy
  accuracy <- calc_accuracy(test_labels, predictions)

  list(
    accuracy    = accuracy,
    predictions = predictions,
    actual      = test_labels
  )
}

