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
    xgb = score_group_xgb,
    mlp = score_group_mlp,
    gnn = score_group_gnn,
    stop("Unknown scoring method: '", method, "'. Available: ",
         "'rf', 'svm', 'knn', 'xgb', 'mlp', 'gnn'.")
  )
  fn(expr, labels, genes, train_idx, ...)
}


#' Score a single group using TRAINING data only (mean F1 via k-fold CV)
#'
#' Computes a ranking score for a group as the mean F1 across stratified
#' k-fold cross-validation on the training partition only (no test leakage),
#' matching the G-S-M Scoring phase. This is the criterion used to rank groups
#' for the accumulative modeling step.
#'
#' @param expr Numeric matrix (genes x samples).
#' @param labels Factor of class labels.
#' @param genes Character vector of gene names for this group.
#' @param train_idx Integer vector of training sample indices.
#' @param method Classifier method (`"rf"`, `"svm"`, `"knn"`).
#' @param cv Number of stratified CV folds (default 3, per the paper).
#' @param ... Additional arguments passed to the classifier.
#' @return A single numeric mean-F1 score (NA if no usable genes).
#' @export
score_group_train <- function(expr, labels, genes, train_idx,
                              method = "rf", cv = 3, ...) {
  stopifnot(is.matrix(expr) && is.numeric(expr))
  stopifnot(is.factor(labels))
  genes <- intersect(genes, rownames(expr))
  if (length(genes) == 0) return(NA_real_)

  x <- t(expr[genes, train_idx, drop = FALSE])    # samples x genes
  y <- droplevels(labels[train_idx])
  if (nlevels(y) < 2) return(NA_real_)
  pos <- if ("pos" %in% levels(y)) "pos" else levels(y)[1]

  folds <- stratified_folds(y, k = cv)
  f1s <- vapply(seq_along(folds), function(i) {
    va <- folds[[i]]; tr <- setdiff(seq_along(y), va)
    if (length(unique(y[tr])) < 2 || length(va) == 0) return(NA_real_)
    pred <- .cv_predict(method, x[tr, , drop = FALSE], y[tr],
                        x[va, , drop = FALSE], ...)
    if (is.null(pred)) return(NA_real_)
    calc_f_measure(y[va], pred, positive_class = pos)
  }, numeric(1))
  mean(f1s, na.rm = TRUE)
}

# Train on (x_tr, y_tr), predict x_va, for one CV fold. Returns a factor or NULL.
.cv_predict <- function(method, x_tr, y_tr, x_va, ...) {
  if (method == "rf") {
    if (!requireNamespace("randomForest", quietly = TRUE))
      stop("Package 'randomForest' is required for method 'rf'.")
    fit <- randomForest::randomForest(x = as.data.frame(x_tr), y = y_tr, ...)
    stats::predict(fit, newdata = as.data.frame(x_va))
  } else if (method == "svm") {
    if (!requireNamespace("e1071", quietly = TRUE))
      stop("Package 'e1071' is required for method 'svm'.")
    fit <- e1071::svm(x = as.data.frame(x_tr), y = y_tr, ...)
    stats::predict(fit, newdata = as.data.frame(x_va))
  } else if (method == "knn") {
    if (!requireNamespace("class", quietly = TRUE))
      stop("Package 'class' is required for method 'knn'.")
    kk <- list(...)$k %||% 3
    class::knn(train = x_tr, test = x_va, cl = y_tr, k = kk)
  } else if (method == "xgb") {
    res <- .xgb_fit(x_tr, y_tr, x_va, ...)
    factor(ifelse(res$prob > 0.5, res$pos, res$neg), levels = levels(droplevels(y_tr)))
  } else if (method == "mlp") {
    res <- .mlp_fit(x_tr, y_tr, x_va, ...)
    factor(ifelse(res$prob > 0.5, res$pos, res$neg), levels = levels(droplevels(y_tr)))
  } else if (method == "gnn") {
    dots <- list(...)
    # x_tr is samples x genes and holds only this inner CV fold's training rows,
    # so t(x_tr) is a leakage-free source for the correlation fallback graph.
    adj <- gsm_gene_adjacency(colnames(x_tr), graph = dots$graph,
                              expr = t(x_tr), train_idx = NULL,
                              cor_threshold = dots$cor_threshold %||% 0.3)
    dots$graph <- NULL; dots$cor_threshold <- NULL
    res <- do.call(.gnn_fit, c(list(x_tr, y_tr, x_va, adj), dots))
    factor(ifelse(res$prob > 0.5, res$pos, res$neg), levels = levels(droplevels(y_tr)))
  } else {
    stop("Unknown scoring method: '", method, "'.")
  }
}

# Fit a binary XGBoost classifier and predict positive-class probabilities.
# Shared by score_group_xgb(), .cv_predict(), and build_model_xgb(). Returns a
# list(fit, prob = P(positive) on x_te, pos, neg, levels). Two-class only.
.xgb_fit <- function(x_tr, y_tr, x_te, positive_class = NULL, ...) {
  if (!requireNamespace("xgboost", quietly = TRUE)) {
    stop("Package 'xgboost' is required for method 'xgb'. ",
         "Install it with: install.packages('xgboost')")
  }
  y_tr <- droplevels(as.factor(y_tr))
  lev <- levels(y_tr)
  if (length(lev) != 2)
    stop("method 'xgb' currently supports two-class problems only.")
  if (is.null(positive_class) || !positive_class %in% lev)
    positive_class <- if ("pos" %in% lev) "pos" else lev[length(lev)]
  neg <- setdiff(lev, positive_class)
  label <- as.integer(y_tr == positive_class)        # 1 = positive class

  dots <- list(...)
  nrounds <- dots$nrounds %||% 50L
  params <- dots$params %||% list()
  params$objective <- params$objective %||% "binary:logistic"
  if (is.null(params$nthread)) params$nthread <- 1   # avoid nested parallelism in MCCV

  # Use the stable low-level xgb.train/xgb.DMatrix API (the high-level
  # xgboost() signature changed in 2.x). The DMatrix keeps the matrix column
  # names, so xgb.importance() reports genes rather than f0/f1/...
  dtrain <- xgboost::xgb.DMatrix(data = as.matrix(x_tr), label = label)
  fit <- xgboost::xgb.train(params = params, data = dtrain,
                            nrounds = nrounds, verbose = 0)
  prob <- as.numeric(stats::predict(fit, as.matrix(x_te)))
  list(fit = fit, prob = prob, pos = positive_class, neg = neg, levels = lev)
}

#' Stratified k-fold assignment (returns a list of validation-index vectors)
#' @param y Factor of labels.
#' @param k Number of folds.
#' @return List of length k; each element holds the row indices of that fold.
#' @export
stratified_folds <- function(y, k = 3) {
  k <- max(2, min(k, min(table(y))))              # k <= smallest class size
  fold_id <- integer(length(y))
  for (lev in levels(y)) {
    idx <- which(y == lev)
    fold_id[idx] <- (sample(seq_along(idx)) - 1) %% k + 1
  }
  split(seq_along(y), fold_id)
}

# Local null-coalesce (base R %||% only from 4.4.0).
`%||%` <- function(a, b) if (is.null(a)) b else a


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

  positive_class <- if ("pos" %in% levels(labels)) "pos" else levels(labels)[1]

  # Train Random Forest
  model <- randomForest::randomForest(
    x = train_data,
    y = train_labels,
    ...
  )

  # Predict on test set (class + positive-class probability for AUC)
  predictions <- stats::predict(model, newdata = test_data)
  prob_pos <- tryCatch({
    pm <- stats::predict(model, newdata = test_data, type = "prob")
    if (positive_class %in% colnames(pm)) pm[, positive_class] else NULL
  }, error = function(e) NULL)

  list(
    accuracy    = calc_accuracy(test_labels, predictions),
    auc         = if (is.null(prob_pos)) NA_real_
                  else calc_auc(test_labels, prob_pos, positive_class),
    predictions = predictions,
    prob        = prob_pos,
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

  positive_class <- if ("pos" %in% levels(labels)) "pos" else levels(labels)[1]

  # Train SVM (probability=TRUE enables positive-class probabilities for AUC)
  model <- e1071::svm(
    x = train_data,
    y = train_labels,
    probability = TRUE,
    ...
  )

  # Predict on test set (class + positive-class probability)
  predictions <- stats::predict(model, newdata = test_data, probability = TRUE)
  prob_pos <- tryCatch({
    pm <- attr(predictions, "probabilities")
    if (!is.null(pm) && positive_class %in% colnames(pm)) pm[, positive_class] else NULL
  }, error = function(e) NULL)

  list(
    accuracy    = calc_accuracy(test_labels, predictions),
    auc         = if (is.null(prob_pos)) NA_real_
                  else calc_auc(test_labels, prob_pos, positive_class),
    predictions = predictions,
    prob        = prob_pos,
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

  list(
    accuracy    = calc_accuracy(test_labels, predictions),
    auc         = if (is.null(prob_pos)) NA_real_
                  else calc_auc(test_labels, prob_pos, positive_class),
    predictions = predictions,
    prob        = prob_pos,
    actual      = test_labels
  )
}


#' Score a group using XGBoost (gradient-boosted trees)
#'
#' Trains a binary XGBoost classifier on the group's genes from the training
#' samples and evaluates on the test samples. Two-class problems only.
#'
#' @param expr Numeric matrix (genes x samples).
#' @param labels Factor of class labels.
#' @param genes Character vector of gene names.
#' @param train_idx Integer vector of training sample indices.
#' @param ... Additional arguments forwarded to [xgboost::xgboost()] (e.g.
#'   `nrounds`, or a `params` list with `max_depth`, `eta`, ...).
#' @return A list with `accuracy`, `auc`, `predictions`, `prob`, and `actual`.
#' @keywords internal
score_group_xgb <- function(expr, labels, genes, train_idx, ...) {
  test_idx <- setdiff(seq_len(ncol(expr)), train_idx)
  if (length(test_idx) == 0) stop("No test samples.")

  train_data <- t(expr[genes, train_idx, drop = FALSE])
  test_data  <- t(expr[genes, test_idx, drop = FALSE])
  test_labels <- labels[test_idx]
  positive_class <- if ("pos" %in% levels(labels)) "pos" else levels(labels)[1]

  res <- .xgb_fit(train_data, labels[train_idx], test_data, positive_class, ...)
  predictions <- factor(ifelse(res$prob > 0.5, res$pos, res$neg),
                        levels = levels(labels))

  list(
    accuracy    = calc_accuracy(test_labels, predictions),
    auc         = calc_auc(test_labels, res$prob, positive_class),
    predictions = predictions,
    prob        = res$prob,
    actual      = test_labels
  )
}

