#' @title Classification metrics
#' @description Functions to compute standard classification metrics matching
#'   the KNIME workflow output columns: Accuracy, Sensitivity, Specificity,
#'   Precision, F-measure, and Cohen's kappa.
#' @name metrics
NULL

#' Calculate classification accuracy
#'
#' @param actual Factor or vector of actual class labels.
#' @param predicted Factor or vector of predicted class labels.
#' @return Numeric accuracy (proportion correct).
#' @export
calc_accuracy <- function(actual, predicted) {
  stopifnot(length(actual) == length(predicted))
  sum(actual == predicted) / length(actual)
}


#' Calculate sensitivity (recall / true positive rate)
#'
#' @param actual Factor or vector of actual class labels.
#' @param predicted Factor or vector of predicted class labels.
#' @param positive_class The label of the positive class.
#' @return Numeric sensitivity.
#' @export
calc_sensitivity <- function(actual, predicted, positive_class = NULL) {
  if (is.null(positive_class)) positive_class <- levels(factor(actual))[1]
  tp <- sum(actual == positive_class & predicted == positive_class)
  fn <- sum(actual == positive_class & predicted != positive_class)
  denom <- tp + fn
  if (denom == 0) return(NA_real_)
  tp / denom
}


#' Calculate specificity (true negative rate)
#'
#' @param actual Factor or vector of actual class labels.
#' @param predicted Factor or vector of predicted class labels.
#' @param positive_class The label of the positive class.
#' @return Numeric specificity.
#' @export
calc_specificity <- function(actual, predicted, positive_class = NULL) {
  if (is.null(positive_class)) positive_class <- levels(factor(actual))[1]
  tn <- sum(actual != positive_class & predicted != positive_class)
  fp <- sum(actual != positive_class & predicted == positive_class)
  denom <- tn + fp
  if (denom == 0) return(NA_real_)
  tn / denom
}


#' Calculate precision (positive predictive value)
#'
#' @param actual Factor or vector of actual class labels.
#' @param predicted Factor or vector of predicted class labels.
#' @param positive_class The label of the positive class.
#' @return Numeric precision.
#' @export
calc_precision <- function(actual, predicted, positive_class = NULL) {
  if (is.null(positive_class)) positive_class <- levels(factor(actual))[1]
  tp <- sum(actual == positive_class & predicted == positive_class)
  fp <- sum(actual != positive_class & predicted == positive_class)
  denom <- tp + fp
  if (denom == 0) return(NA_real_)
  tp / denom
}


#' Calculate F-measure (F1 score)
#'
#' @param actual Factor or vector of actual class labels.
#' @param predicted Factor or vector of predicted class labels.
#' @param positive_class The label of the positive class.
#' @return Numeric F1 score.
#' @export
calc_f_measure <- function(actual, predicted, positive_class = NULL) {
  p <- calc_precision(actual, predicted, positive_class)
  r <- calc_sensitivity(actual, predicted, positive_class)
  if (is.na(p) || is.na(r) || (p + r) == 0) return(NA_real_)
  2 * p * r / (p + r)
}


#' Calculate Cohen's kappa
#'
#' Measures inter-rater agreement adjusted for chance.
#'
#' @param actual Factor or vector of actual class labels.
#' @param predicted Factor or vector of predicted class labels.
#' @return Numeric Cohen's kappa.
#' @export
calc_cohens_kappa <- function(actual, predicted) {
  stopifnot(length(actual) == length(predicted))
  n <- length(actual)
  if (n == 0) return(NA_real_)

  # Confusion counts
  lvls <- union(levels(factor(actual)), levels(factor(predicted)))
  actual <- factor(actual, levels = lvls)
  predicted <- factor(predicted, levels = lvls)
  cm <- table(actual, predicted)

  # Observed agreement
  po <- sum(diag(cm)) / n

  # Expected agreement by chance
  row_sums <- rowSums(cm)
  col_sums <- colSums(cm)
  pe <- sum(row_sums * col_sums) / (n * n)

  if (pe == 1) return(1.0)
  (po - pe) / (1 - pe)
}


#' Calculate all classification metrics
#'
#' Convenience function returning all metrics as a named list. Column names
#' match the KNIME workflow output.
#'
#' @param actual Factor or vector of actual class labels.
#' @param predicted Factor or vector of predicted class labels.
#' @param positive_class The label of the positive class.
#' @param prob Optional numeric vector of predicted positive-class probabilities
#'   (same length/order as `actual`); enables AUC-ROC.
#' @return Named list with: Accuracy, Sensitivity, Specificity, Precision,
#'   F_measure, Cohens_kappa, AUC.
#' @export
calc_all_metrics <- function(actual, predicted, positive_class = NULL,
                             prob = NULL) {
  if (is.null(positive_class)) positive_class <- levels(factor(actual))[1]
  list(
    Accuracy     = calc_accuracy(actual, predicted),
    Sensitivity  = calc_sensitivity(actual, predicted, positive_class),
    Specificity  = calc_specificity(actual, predicted, positive_class),
    Precision    = calc_precision(actual, predicted, positive_class),
    F_measure    = calc_f_measure(actual, predicted, positive_class),
    Cohens_kappa = calc_cohens_kappa(actual, predicted),
    AUC          = if (is.null(prob)) NA_real_
                   else calc_auc(actual, prob, positive_class)
  )
}

#' Area under the ROC curve (rank-based / Mann-Whitney)
#'
#' @param actual Factor or vector of actual labels.
#' @param prob Numeric predicted probabilities of the positive class.
#' @param positive_class The label of the positive class.
#' @return AUC in [0, 1]; `NA` if only one class is present.
#' @export
calc_auc <- function(actual, prob, positive_class = NULL) {
  if (is.null(positive_class)) positive_class <- levels(factor(actual))[1]
  pos <- actual == positive_class
  n_pos <- sum(pos); n_neg <- sum(!pos)
  if (n_pos == 0 || n_neg == 0 || anyNA(prob)) return(NA_real_)
  r <- rank(prob)                                   # average ranks (ties handled)
  (sum(r[pos]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

#' Percentile confidence interval from a vector of repeated measurements
#'
#' @param x Numeric vector (e.g. a metric across MCCV iterations).
#' @param probs Lower/upper probabilities (default 2.5\% / 97.5\%).
#' @return Named numeric `c(low, high)`; `NA`s if `x` has no finite values.
#' @export
pct_ci <- function(x, probs = c(0.025, 0.975)) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(c(low = NA_real_, high = NA_real_))
  q <- stats::quantile(x, probs = probs, names = FALSE, na.rm = TRUE)
  c(low = q[1], high = q[2])
}
