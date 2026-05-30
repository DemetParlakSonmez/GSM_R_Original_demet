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
#' @param method Pre-scoring method. One of `"anova"` (default) or `"none"`.
#' @param top_percent Percentage of top groups to retain. Default `20`.
#' @param ... Additional arguments passed to the specific method.
#' @return A named list of character vectors — the subset of groups that
#'   passed pre-scoring, ordered by significance (most significant first).
#' @export
pre_score_groups <- function(dataset, train_idx, method = "anova",
                             top_percent = 20, ...) {
  stopifnot(inherits(dataset, "gsm_dataset"))
  stopifnot(is.numeric(train_idx) && length(train_idx) > 0)
  stopifnot(top_percent > 0 && top_percent <= 100)

  fn <- switch(method,
    anova = pre_score_anova,
    none  = pre_score_none,
    stop("Unknown pre-scoring method: '", method, "'. ",
         "Available methods: 'anova', 'none'.")
  )
  fn(dataset, train_idx, top_percent = top_percent, ...)
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
#' @param ... Ignored.
#' @return Named list of character vectors (top groups).
#' @keywords internal
pre_score_anova <- function(dataset, train_idx, top_percent = 20, ...) {
  expr <- dataset$expr[, train_idx, drop = FALSE]
  labels <- dataset$labels[train_idx]
  groups <- dataset$groups

  # Score each group
  group_scores <- vapply(names(groups), \(gname) {
    genes <- groups[[gname]]
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
    fisher_stat <- -2 * sum(log(pvals))
    fisher_stat
  }, numeric(1))

  # Rank by score (descending — higher = more significant)
  sorted_names <- names(sort(group_scores, decreasing = TRUE))

  # Keep top percent

  n_keep <- max(1, ceiling(length(sorted_names) * top_percent / 100))
  top_names <- sorted_names[seq_len(n_keep)]

  groups[top_names]
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
pre_score_none <- function(dataset, train_idx, ...) {
  dataset$groups
}
