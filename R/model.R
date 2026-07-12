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
#'     \item{importance}{Named numeric of per-gene RandomForest importance
#'       (MeanDecreaseGini) for `method = "rf"`; `NULL` for SVM/KNN.}
#'   }
#' @export
build_model <- function(expr, labels, top_groups, train_idx,
                        method = "rf", ...) {
  stopifnot(is.matrix(expr) && is.numeric(expr))
  stopifnot(is.factor(labels))
  stopifnot(is.list(top_groups) && length(top_groups) > 0)

  fn <- switch(method,
    rf  = build_model_rf,
    svm = build_model_svm,
    knn = build_model_knn,
    xgb = build_model_xgb,
    mlp = build_model_mlp,
    gnn = build_model_gnn,
    stop("Unknown modeling method: '", method, "'. Available: ",
         "'rf', 'svm', 'knn', 'xgb', 'mlp', 'gnn'.")
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
      n_groups    = length(top_groups),
      importance  = NULL
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

  # Predict (class + positive-class probability for AUC)
  predictions <- stats::predict(model, newdata = test_data)
  prob_pos <- tryCatch({
    pm <- stats::predict(model, newdata = test_data, type = "prob")
    if (positive_class %in% colnames(pm)) pm[, positive_class] else NULL
  }, error = function(e) NULL)

  # Metrics
  metrics <- calc_all_metrics(test_labels, predictions, positive_class, prob = prob_pos)

  # Per-gene Random Forest importance (MeanDecreaseGini), named by gene.
  importance <- tryCatch({
    imp <- model$importance
    stats::setNames(as.numeric(imp[, ncol(imp)]), rownames(imp))
  }, error = function(e) NULL)

  list(
    metrics     = metrics,
    predictions = predictions,
    actual      = test_labels,
    features    = all_genes,
    n_groups    = length(top_groups),
    importance  = importance
  )
}


#' Build an SVM model from top groups
#'
#' @param expr Numeric matrix (genes x samples).
#' @param labels Factor of class labels.
#' @param top_groups Named list of character vectors.
#' @param train_idx Integer vector of training sample indices.
#' @param ... Additional arguments passed to [e1071::svm()].
#' @return A list with `metrics`, `predictions`, `actual`, `features`,
#'   and `n_groups`.
#' @keywords internal
build_model_svm <- function(expr, labels, top_groups, train_idx, ...) {
  if (!requireNamespace("e1071", quietly = TRUE)) {
    stop("Package 'e1071' is required for method 'svm'. ",
         "Install it with: install.packages('e1071')")
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
      n_groups    = length(top_groups),
      importance  = NULL
    ))
  }

  # Build data frames
  train_data <- as.data.frame(t(expr[all_genes, train_idx, drop = FALSE]))
  test_data  <- as.data.frame(t(expr[all_genes, test_idx, drop = FALSE]))

  train_labels <- labels[train_idx]
  test_labels  <- labels[test_idx]

  # Determine positive class
  positive_class <- if ("pos" %in% levels(labels)) "pos" else levels(labels)[1]

  # Train SVM (probability=TRUE enables positive-class probabilities for AUC)
  model <- e1071::svm(
    x = train_data,
    y = train_labels,
    probability = TRUE,
    ...
  )

  # Predict (class + positive-class probability)
  predictions <- stats::predict(model, newdata = test_data, probability = TRUE)
  prob_pos <- tryCatch({
    pm <- attr(predictions, "probabilities")
    if (!is.null(pm) && positive_class %in% colnames(pm)) pm[, positive_class] else NULL
  }, error = function(e) NULL)

  # Metrics
  metrics <- calc_all_metrics(test_labels, predictions, positive_class, prob = prob_pos)

  # Per-gene importance is well-defined only for the LINEAR kernel, where the
  # decision function is w·x + b and |w_g| measures gene g's contribution.
  # e1071 encodes the kernel as an integer (0 = linear). w = t(coefs) %*% SV.
  importance <- tryCatch({
    if (identical(as.integer(model$kernel), 0L)) {
      w <- as.numeric(t(model$coefs) %*% model$SV)
      stats::setNames(abs(w), colnames(model$SV))[all_genes]
    } else NULL
  }, error = function(e) NULL)

  list(
    metrics     = metrics,
    predictions = predictions,
    actual      = test_labels,
    features    = all_genes,
    n_groups    = length(top_groups),
    importance  = importance
  )
}


#' Build a KNN model from top groups
#'
#' @param expr Numeric matrix (genes x samples).
#' @param labels Factor of class labels.
#' @param top_groups Named list of character vectors.
#' @param train_idx Integer vector of training sample indices.
#' @param ... Additional arguments passed to [class::knn()].
#' @return A list with `metrics`, `predictions`, `actual`, `features`,
#'   and `n_groups`.
#' @keywords internal
build_model_knn <- function(expr, labels, top_groups, train_idx, ...) {
  if (!requireNamespace("class", quietly = TRUE)) {
    stop("Package 'class' is required for method 'knn'.")
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
      n_groups    = length(top_groups),
      importance  = NULL
    ))
  }

  # Build data frames
  train_data <- t(expr[all_genes, train_idx, drop = FALSE])
  test_data  <- t(expr[all_genes, test_idx, drop = FALSE])

  train_labels <- labels[train_idx]
  test_labels  <- labels[test_idx]

  # Determine positive class
  positive_class <- if ("pos" %in% levels(labels)) "pos" else levels(labels)[1]

  # Predict using KNN (prob=TRUE returns the winning-class vote proportion)
  predictions <- class::knn(
    train = train_data,
    test  = test_data,
    cl    = train_labels,
    prob  = TRUE,
    ...
  )
  prob_pos <- tryCatch({
    pr <- attr(predictions, "prob")
    if (is.null(pr)) NULL else ifelse(predictions == positive_class, pr, 1 - pr)
  }, error = function(e) NULL)

  # Metrics
  metrics <- calc_all_metrics(test_labels, predictions, positive_class, prob = prob_pos)

  list(
    metrics     = metrics,
    predictions = predictions,
    actual      = test_labels,
    features    = all_genes,
    n_groups    = length(top_groups),
    importance  = NULL
  )
}


#' Build an XGBoost model from top groups
#'
#' Combines the genes of the top-k groups and trains a binary gradient-boosted
#' tree classifier. Returns per-gene importance (XGBoost `Gain`), enabling the
#' consensus biomarker ranking in [gsm_biomarkers()]. Two-class problems only.
#'
#' @param expr Numeric matrix (genes x samples).
#' @param labels Factor of class labels.
#' @param top_groups Named list of character vectors.
#' @param train_idx Integer vector of training sample indices.
#' @param ... Additional arguments forwarded to [xgboost::xgboost()].
#' @return A list with `metrics`, `predictions`, `actual`, `features`,
#'   `n_groups`, and `importance` (per-gene Gain, named by gene).
#' @keywords internal
build_model_xgb <- function(expr, labels, top_groups, train_idx, ...) {
  test_idx <- setdiff(seq_len(ncol(expr)), train_idx)
  if (length(test_idx) == 0) stop("No test samples.")

  all_genes <- intersect(unique(unlist(top_groups)), rownames(expr))
  if (length(all_genes) == 0) {
    return(list(
      metrics     = calc_all_metrics(labels[test_idx], labels[test_idx]),
      predictions = factor(rep(NA, length(test_idx)), levels = levels(labels)),
      actual      = labels[test_idx],
      features    = character(0),
      n_groups    = length(top_groups),
      importance  = NULL
    ))
  }

  train_data <- t(expr[all_genes, train_idx, drop = FALSE])
  test_data  <- t(expr[all_genes, test_idx, drop = FALSE])
  test_labels <- labels[test_idx]
  positive_class <- if ("pos" %in% levels(labels)) "pos" else levels(labels)[1]

  res <- .xgb_fit(train_data, labels[train_idx], test_data, positive_class, ...)
  predictions <- factor(ifelse(res$prob > 0.5, res$pos, res$neg),
                        levels = levels(labels))
  metrics <- calc_all_metrics(test_labels, predictions, positive_class, prob = res$prob)

  # Per-gene importance (Gain), always named by gene. Genes the model never
  # split on are absent from xgb.importance(); leave them NA (not 0) so
  # gsm_biomarkers() treats them as not selected by this fold.
  importance <- tryCatch({
    im <- xgboost::xgb.importance(model = res$fit)
    gain <- stats::setNames(rep(NA_real_, length(all_genes)), all_genes)
    hit <- intersect(im$Feature, all_genes)
    gain[hit] <- im$Gain[match(hit, im$Feature)]
    gain
  }, error = function(e) NULL)

  list(
    metrics     = metrics,
    predictions = predictions,
    actual      = test_labels,
    features    = all_genes,
    n_groups    = length(top_groups),
    importance  = importance
  )
}

