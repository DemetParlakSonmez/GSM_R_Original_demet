#' @title Pre-scoring: statistical group filtering
#' @description Abstract interface for pre-scoring gene groups using
#'   statistical methods. The pre-scoring step reduces the number of groups
#'   before the more expensive ML-based scoring step.
#' @name pre_score
NULL

#' Pre-score groups to select top candidates
#'
#' An abstract function that dispatches to a specific pre-scoring method.
#' The goal is to rank groups by statistical significance and retain only
#' the top-ranked fraction.
#'
#' @param dataset A `gsm_dataset` object.
#' @param train_idx Integer vector of training sample indices.
#' @param method Pre-scoring method. One of `"anova"` (default), `"limma"`, or
#'   `"none"`.
#' @param top_percent Percentage of top groups to retain. Default `20`.
#' @param progress Optional callback invoked once per group as
#'   `progress(i, n, group_name, score)`. Useful for live UI progress.
#' @param ... Additional arguments passed to the specific method.
#' @return A named list of character vectors — the subset of groups that
#'   passed pre-scoring, ordered by significance (most significant first).
#'   The full per-group scores (all groups, sorted descending) are attached
#'   as `attr(result, "scores")`, and the cutoff count as `attr(result, "n_keep")`.
#' @export
pre_score_groups <- function(dataset, train_idx, method = "anova",
                             top_percent = 20, progress = NULL, ...) {
  stopifnot(inherits(dataset, "gsm_dataset"))
  stopifnot(is.numeric(train_idx) && length(train_idx) > 0)
  stopifnot(top_percent > 0 && top_percent <= 100)

  fn <- switch(method,
    anova = pre_score_anova,
    limma = pre_score_limma,
    none  = pre_score_none,
    stop("Unknown pre-scoring method: '", method, "'. ",
         "Available methods: 'anova', 'limma', 'none'.")
  )
  fn(dataset, train_idx, top_percent = top_percent, progress = progress, ...)
}


#' Pre-score groups using one-way ANOVA
#'
#' For each group, computes a summary statistic from one-way ANOVA tests
#' on each gene in the group (training data only). Groups are ranked by
#' their aggregate p-value (Fisher's method).
#'
#' @param dataset A `gsm_dataset` object.
#' @param train_idx Integer vector of training sample indices.
#' @param top_percent Percentage of top groups to retain.
#' @param progress Optional callback `progress(i, n, group_name, score)`.
#' @param ... Ignored.
#' @return Named list of character vectors (top groups). The full sorted
#'   per-group scores are attached as `attr(., "scores")`.
#' @keywords internal
pre_score_anova <- function(dataset, train_idx, top_percent = 20,
                            progress = NULL, ...) {
  expr <- dataset$expr[, train_idx, drop = FALSE]
  labels <- dataset$labels[train_idx]
  groups <- dataset$groups

  gnames <- names(groups)
  n <- length(gnames)
  group_scores <- numeric(n)
  names(group_scores) <- gnames

  # Score each group
  for (j in seq_len(n)) {
    genes <- groups[[gnames[j]]]
    # Get expression for these genes (training only)
    gene_expr <- expr[genes, , drop = FALSE]

    # Run one-way ANOVA for each gene, collect p-values
    pvals <- vapply(seq_len(nrow(gene_expr)), \(i) {
      tryCatch({
        fit <- stats::oneway.test(gene_expr[i, ] ~ labels)
        fit$p.value
      }, error = \(e) 1.0)
    }, numeric(1))

    # Aggregate using Fisher's method: -2 * sum(log(p))
    # Higher value = more significant
    pvals[pvals == 0] <- .Machine$double.xmin  # avoid log(0)
    group_scores[j] <- -2 * sum(log(pvals))

    if (!is.null(progress)) progress(j, n, gnames[j], group_scores[j])
  }

  # Rank by score (descending — higher = more significant)
  sorted_scores <- sort(group_scores, decreasing = TRUE)

  # Keep top percent
  n_keep <- max(1, ceiling(n * top_percent / 100))
  top_names <- names(sorted_scores)[seq_len(n_keep)]

  result <- groups[top_names]
  attr(result, "scores") <- sorted_scores  # all groups, most significant first
  attr(result, "n_keep") <- n_keep
  result
}


#' Pre-score groups using limma moderated statistics
#'
#' Computes per-gene moderated test statistics with \pkg{limma} (empirical-Bayes
#' shrinkage of gene-wise variances) on the training data, then aggregates each
#' group's gene p-values with Fisher's method — the same group statistic as
#' [pre_score_anova()], but with limma's variance moderation, which is more
#' powerful for the small sample sizes typical of transcriptomic studies. For two
#' classes this is the moderated t-test; for more than two it is the moderated
#' F-test across all class contrasts.
#'
#' @param dataset A `gsm_dataset` object.
#' @param train_idx Integer vector of training sample indices.
#' @param top_percent Percentage of top groups to retain.
#' @param progress Optional callback `progress(i, n, group_name, score)`.
#' @param ... Ignored.
#' @return Named list of character vectors (top groups), with the full sorted
#'   per-group scores attached as `attr(., "scores")`.
#' @keywords internal
pre_score_limma <- function(dataset, train_idx, top_percent = 20,
                            progress = NULL, ...) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Package 'limma' is required for method 'limma'. ",
         "Install it with: BiocManager::install('limma')")
  }
  expr <- dataset$expr[, train_idx, drop = FALSE]
  labels <- droplevels(dataset$labels[train_idx])
  groups <- dataset$groups

  # One moderated fit over all genes; design = intercept + class effect(s).
  design <- stats::model.matrix(~ labels)
  fit <- limma::eBayes(limma::lmFit(expr, design))
  # F.p.value tests all non-intercept coefficients jointly (works for >2 classes).
  pvals <- fit$F.p.value
  names(pvals) <- rownames(expr)
  pvals[!is.finite(pvals)] <- 1
  pvals[pvals == 0] <- .Machine$double.xmin

  gnames <- names(groups)
  n <- length(gnames)
  group_scores <- numeric(n)
  names(group_scores) <- gnames
  for (j in seq_len(n)) {
    gp <- pvals[intersect(groups[[gnames[j]]], names(pvals))]
    # Fisher's method: -2 * sum(log(p)); higher = more significant.
    group_scores[j] <- if (length(gp)) -2 * sum(log(gp)) else 0
    if (!is.null(progress)) progress(j, n, gnames[j], group_scores[j])
  }

  sorted_scores <- sort(group_scores, decreasing = TRUE)
  n_keep <- max(1, ceiling(n * top_percent / 100))
  top_names <- names(sorted_scores)[seq_len(n_keep)]

  result <- groups[top_names]
  attr(result, "scores") <- sorted_scores
  attr(result, "n_keep") <- n_keep
  result
}


#' Pass-through pre-scoring (no filtering)
#'
#' Returns all groups unchanged. Useful when no pre-scoring is desired.
#'
#' @param dataset A `gsm_dataset` object.
#' @param train_idx Ignored.
#' @param ... Ignored.
#' @return The full groups list from the dataset.
#' @keywords internal
pre_score_none <- function(dataset, train_idx, progress = NULL, ...) {
  result <- dataset$groups
  attr(result, "scores") <- NULL
  attr(result, "n_keep") <- length(result)
  result
}
