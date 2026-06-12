#' @title Permutation-based empirical significance of gene groups
#' @description Tests whether a group's score is better than expected from a
#'   random gene set of the same size. For each group, the observed score
#'   (training-only mean F1, the G-S-M Scoring criterion) is compared against a
#'   null distribution of scores from `n_perm` random gene sets drawn from a
#'   background universe. Reports an empirical p-value and a Benjamini-Hochberg
#'   adjusted p-value (FDR) across the tested groups.
#' @name permutation
NULL

#' Empirical p-values for gene groups via gene-set permutation
#'
#' Quantifies how surprising each group's score is relative to random gene sets
#' of the same size. This turns a bare ranking ("group X scored highest") into a
#' defensible claim ("group X scores higher than 199/200 random sets of the same
#' size, FDR < 0.01"). Scoring is leakage-free: the criterion is the mean F1 from
#' stratified k-fold CV on the **training** partition only (see
#' [score_group_train()]), so observed and null scores are directly comparable.
#'
#' @param dataset A `gsm_dataset` object.
#' @param train_idx Integer vector of training sample indices. Defaults to all
#'   samples (use a held-out split when you also want an unbiased test estimate).
#' @param group_names Optional character vector naming which groups to test.
#'   Defaults to every group in the dataset.
#' @param method Classifier used for scoring (`"rf"`, `"svm"`, `"knn"`,
#'   `"xgb"`). Default `"rf"`.
#' @param n_perm Number of random gene sets drawn per group. Default `200`.
#' @param universe Optional character vector of background genes to sample from.
#'   Defaults to all genes in the expression matrix (`rownames(dataset$expr)`).
#' @param cv Number of stratified CV folds for the scoring criterion. Default `3`.
#' @param score_args Named list of extra arguments forwarded to the classifier.
#' @param seed Optional integer base seed for reproducible permutations.
#' @param progress Optional callback `progress(i, n, group_name)` invoked once
#'   per group, for live UI updates.
#' @return A data.frame, one row per tested group, sorted by `p_value`:
#'   \describe{
#'     \item{Group}{Group name.}
#'     \item{size}{Number of group genes present in the universe.}
#'     \item{observed}{Observed group score (mean F1).}
#'     \item{null_mean}{Mean score of the random gene sets.}
#'     \item{p_value}{Empirical p-value `(1 + #{null >= observed}) / (n_perm + 1)`.}
#'     \item{p_adj}{Benjamini-Hochberg adjusted p-value (FDR) across groups.}
#'     \item{n_perm}{Permutations actually used for this group.}
#'   }
#' @export
permutation_test_groups <- function(dataset, train_idx = NULL,
                                     group_names = NULL, method = "rf",
                                     n_perm = 200, universe = NULL, cv = 3,
                                     score_args = list(), seed = NULL,
                                     progress = NULL) {
  stopifnot(inherits(dataset, "gsm_dataset"))
  stopifnot(is.numeric(n_perm) && n_perm >= 1)
  stopifnot(is.list(score_args))

  expr <- dataset$expr
  labels <- dataset$labels
  if (is.null(train_idx)) train_idx <- seq_len(ncol(expr))
  if (is.null(group_names)) group_names <- names(dataset$groups)
  group_names <- intersect(group_names, names(dataset$groups))
  if (length(group_names) == 0) stop("No matching groups to test.")

  if (is.null(universe)) universe <- rownames(expr)
  universe <- intersect(universe, rownames(expr))
  if (length(universe) < 2) stop("Background universe has fewer than 2 genes.")

  score_one <- function(genes) {
    do.call(score_group_train, c(
      list(expr = expr, labels = labels, genes = genes,
           train_idx = train_idx, method = method, cv = cv),
      score_args))
  }

  rows <- vector("list", length(group_names))
  for (j in seq_along(group_names)) {
    gname <- group_names[j]
    genes <- intersect(dataset$groups[[gname]], universe)
    size <- length(genes)

    if (size == 0) {
      rows[[j]] <- data.frame(Group = gname, size = 0L, observed = NA_real_,
                              null_mean = NA_real_, p_value = NA_real_,
                              p_adj = NA_real_, n_perm = 0L,
                              stringsAsFactors = FALSE)
      if (!is.null(progress)) progress(j, length(group_names), gname)
      next
    }

    observed <- score_one(genes)

    # Null: random same-size gene sets from the background universe. A fixed,
    # group-offset seed keeps the whole test reproducible without Math.random.
    if (!is.null(seed)) set.seed(seed + j)
    null_scores <- vapply(seq_len(n_perm), function(p) {
      rnd <- sample(universe, size)
      score_one(rnd)
    }, numeric(1))

    null_ok <- null_scores[is.finite(null_scores)]
    n_used <- length(null_ok)
    p_val <- if (!is.finite(observed) || n_used == 0) NA_real_
             else (1 + sum(null_ok >= observed)) / (n_used + 1)

    rows[[j]] <- data.frame(
      Group = gname, size = size, observed = observed,
      null_mean = if (n_used) mean(null_ok) else NA_real_,
      p_value = p_val, p_adj = NA_real_, n_perm = n_used,
      stringsAsFactors = FALSE)

    if (!is.null(progress)) progress(j, length(group_names), gname)
  }

  out <- do.call(rbind, rows)
  out$p_adj <- stats::p.adjust(out$p_value, method = "BH")
  out <- out[order(out$p_value, -out$observed, na.last = TRUE), , drop = FALSE]
  rownames(out) <- NULL
  out
}
