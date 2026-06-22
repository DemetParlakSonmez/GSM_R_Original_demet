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
#' @param method Pre-scoring method. One of `"anova"` (default), `"limma"`,
#'   `"fisher_chisq"`, `"ebm"`, `"camera"`, or `"none"`. See *Details*.
#' @param top_percent Percentage of top groups to retain. Default `20`.
#' @param progress Optional callback invoked once per group as
#'   `progress(i, n, group_name, score)`. Useful for live UI progress.
#' @param ... Additional arguments passed to the specific method.
#' @return A named list of character vectors — the subset of groups that
#'   passed pre-scoring, ordered by significance (most significant first).
#'   The full per-group scores (all groups, sorted descending; higher = more
#'   significant) are attached as `attr(result, "scores")`, and the cutoff count
#'   as `attr(result, "n_keep")`. The newer p-value-based methods additionally
#'   attach the raw per-group p-values as `attr(result, "pvalues")` and their
#'   Benjamini-Hochberg adjusted values as `attr(result, "padj")`.
#' @details
#' The classic methods (`"anova"`, `"limma"`) aggregate per-gene p-values with
#' Fisher's combined statistic `-2 * sum(log(p))`, which grows with group size
#' and assumes independent genes. The newer methods address those biases:
#' \itemize{
#'   \item `"fisher_chisq"` — same per-gene ANOVA p-values, but the Fisher
#'     statistic is calibrated to a p-value via its chi-square(2k) null, removing
#'     the monotonic group-size inflation.
#'   \item `"ebm"` — Empirical Brown's Method: a correlation-corrected
#'     adaptation of Fisher's method that estimates inter-gene covariance from
#'     the data (Poole et al., Bioinformatics 2016; extends Brown 1975).
#'   \item `"camera"` — a competitive gene-set test that adjusts for inter-gene
#'     correlation (Wu & Smyth, NAR 2012; `limma::camera`), asking whether a
#'     group is *more* associated than the background. Two-class only; for >2
#'     classes it falls back to `"fisher_chisq"` with a warning.
#'   \item `"gsva"` — single-sample enrichment: each group is summarised into a
#'     per-sample activity score (ssGSEA/GSVA via the \pkg{GSVA} package when
#'     installed, otherwise a built-in combined z-score), then ranked by how
#'     well that activity separates the classes. This scores the same group-level
#'     signal the downstream ML classifier consumes (Hänzelmann et al. 2013).
#' }
#' @export
pre_score_groups <- function(dataset, train_idx, method = "anova",
                             top_percent = 20, progress = NULL, ...) {
  stopifnot(inherits(dataset, "gsm_dataset"))
  stopifnot(is.numeric(train_idx) && length(train_idx) > 0)
  stopifnot(top_percent > 0 && top_percent <= 100)

  fn <- switch(method,
    anova        = pre_score_anova,
    limma        = pre_score_limma,
    fisher_chisq = pre_score_fisher_chisq,
    ebm          = pre_score_ebm,
    camera       = pre_score_camera,
    gsva         = pre_score_gsva,
    none         = pre_score_none,
    stop("Unknown pre-scoring method: '", method, "'. ",
         "Available methods: 'anova', 'limma', 'fisher_chisq', 'ebm', ",
         "'camera', 'gsva', 'none'.")
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


# ──────────────────────────────────────────────────────────────────────────────
# Shared helpers for the p-value-based pre-scoring methods
# ──────────────────────────────────────────────────────────────────────────────

#' Per-gene one-way ANOVA p-values
#'
#' Welch one-way ANOVA p-value for every row (gene) of `expr`. Mirrors the
#' per-gene test used by [pre_score_anova()] so the calibrated/correlation-aware
#' methods build on the same gene-level statistics.
#'
#' @param expr Numeric matrix, genes (rows) x samples (columns), training only.
#' @param labels Factor of class labels, length `ncol(expr)`.
#' @return Named numeric vector of p-values, one per row of `expr`.
#' @keywords internal
.pre_score_pergene_anova <- function(expr, labels) {
  p <- vapply(seq_len(nrow(expr)), function(i) {
    tryCatch(stats::oneway.test(expr[i, ] ~ labels)$p.value, error = function(e) 1.0)
  }, numeric(1))
  names(p) <- rownames(expr)
  p
}

#' Finalise a p-value-based pre-scoring result
#'
#' Given one p-value per group, derives the significance score
#' (`-log10(p)`, so higher = more significant — matching the `"scores"` contract
#' of the Fisher-based methods), applies the top-percent cutoff and attaches the
#' raw p-values plus their Benjamini-Hochberg adjustment.
#'
#' @param groups The full named list of groups (for subsetting).
#' @param pvals Named numeric vector of per-group p-values (one per group).
#' @param top_percent Percentage of top groups to retain.
#' @return Named list of top groups with `scores`, `pvalues`, `padj` and
#'   `n_keep` attributes.
#' @keywords internal
.pre_score_finalize <- function(groups, pvals, top_percent) {
  pvals <- pvals[names(groups)]               # align order to groups
  pvals[!is.finite(pvals)] <- 1
  pvals[pvals <= 0] <- .Machine$double.xmin
  padj   <- stats::p.adjust(pvals, method = "BH")
  scores <- -log10(pvals)                     # higher = more significant

  ord    <- order(scores, decreasing = TRUE)
  scores <- scores[ord]; pvals <- pvals[ord]; padj <- padj[ord]

  n      <- length(scores)
  n_keep <- max(1, ceiling(n * top_percent / 100))
  top_names <- names(scores)[seq_len(n_keep)]

  result <- groups[top_names]
  attr(result, "scores")  <- scores   # all groups, most significant first
  attr(result, "pvalues") <- pvals
  attr(result, "padj")    <- padj
  attr(result, "n_keep")  <- n_keep
  result
}

#' Brown's empirical transform of a single gene's expression vector
#'
#' Maps each value to `-2 * log(F(z))` where `F` is the empirical CDF of the
#' standardised values — the transform underlying Empirical Brown's Method.
#'
#' @param x Numeric vector (one gene's training-sample values).
#' @return Numeric vector the same length as `x` (all zeros if `x` is constant).
#' @keywords internal
.ebm_transform <- function(x) {
  m <- mean(x); s <- stats::sd(x)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  z <- (x - m) / s
  -2 * log(stats::ecdf(z)(z))
}

#' Combine dependent p-values with Empirical Brown's Method
#'
#' Correlation-corrected adaptation of Fisher's method (Poole et al.,
#' Bioinformatics 2016). Inter-gene covariance is estimated from `data_matrix`;
#' when genes are uncorrelated this reduces exactly to Fisher's chi-square(2k).
#'
#' @param pvals Numeric vector of per-gene p-values.
#' @param data_matrix Numeric matrix, genes (rows) x samples (columns), the same
#'   data the p-values were computed from.
#' @return Single combined p-value.
#' @keywords internal
.ebm_combine <- function(pvals, data_matrix) {
  N <- length(pvals)
  if (N == 0) return(1)
  if (N == 1) return(unname(pvals[1]))

  wt    <- apply(data_matrix, 1, .ebm_transform)  # samples x genes
  covar <- stats::cov(wt)                          # genes x genes
  cov_sum  <- 2 * sum(covar[lower.tri(covar)])
  expected <- 2 * N
  varr     <- 4 * N + cov_sum
  cc       <- varr / (2 * expected)
  df_brown <- 2 * expected^2 / varr
  df_fisher <- 2 * N
  if (!is.finite(df_brown) || df_brown > df_fisher) { df_brown <- df_fisher; cc <- 1 }

  x <- -2 * sum(log(pvals))
  stats::pchisq(x / cc, df = df_brown, lower.tail = FALSE)
}


#' Pre-score groups using size-calibrated Fisher (chi-square)
#'
#' Computes per-gene Welch one-way ANOVA p-values (as [pre_score_anova()]) but
#' calibrates each group's Fisher statistic `-2 * sum(log(p))` to a p-value via
#' its chi-square(2k) null distribution, where k is the group size. This removes
#' the group-size inflation of the raw Fisher score: large and small groups are
#' placed on a common p-value scale. Groups are then ranked by `-log10(p)` and
#' Benjamini-Hochberg adjusted p-values are attached.
#'
#' @param dataset A `gsm_dataset` object.
#' @param train_idx Integer vector of training sample indices.
#' @param top_percent Percentage of top groups to retain.
#' @param progress Optional callback `progress(i, n, group_name, score)`.
#' @param ... Ignored.
#' @return Named list of character vectors (top groups), with `scores`,
#'   `pvalues`, `padj` and `n_keep` attributes.
#' @keywords internal
pre_score_fisher_chisq <- function(dataset, train_idx, top_percent = 20,
                                   progress = NULL, ...) {
  expr   <- dataset$expr[, train_idx, drop = FALSE]
  labels <- droplevels(dataset$labels[train_idx])
  groups <- dataset$groups

  universe <- intersect(unique(unlist(groups)), rownames(expr))
  gene_p   <- .pre_score_pergene_anova(expr[universe, , drop = FALSE], labels)

  gnames <- names(groups)
  n <- length(gnames)
  pg <- numeric(n); names(pg) <- gnames
  for (j in seq_len(n)) {
    gp <- gene_p[intersect(groups[[gnames[j]]], names(gene_p))]
    k  <- length(gp)
    if (k == 0) { pg[j] <- 1 } else {
      gp[gp <= 0] <- .Machine$double.xmin
      stat  <- -2 * sum(log(gp))
      pg[j] <- stats::pchisq(stat, df = 2 * k, lower.tail = FALSE)
    }
    if (!is.null(progress))
      progress(j, n, gnames[j], -log10(max(pg[j], .Machine$double.xmin)))
  }

  .pre_score_finalize(groups, pg, top_percent)
}


#' Pre-score groups using Empirical Brown's Method
#'
#' Like [pre_score_fisher_chisq()] but combines the per-gene p-values with
#' Empirical Brown's Method, which corrects the Fisher combination for the
#' inter-gene correlation present in co-regulated pathways (Poole et al.,
#' Bioinformatics 2016). Covariance is estimated per group from the training
#' expression of that group's genes. Returns one calibrated p-value per group;
#' ranking and BH adjustment match the other p-value-based methods.
#'
#' @param dataset A `gsm_dataset` object.
#' @param train_idx Integer vector of training sample indices.
#' @param top_percent Percentage of top groups to retain.
#' @param progress Optional callback `progress(i, n, group_name, score)`.
#' @param ... Ignored.
#' @return Named list of character vectors (top groups), with `scores`,
#'   `pvalues`, `padj` and `n_keep` attributes.
#' @keywords internal
pre_score_ebm <- function(dataset, train_idx, top_percent = 20,
                          progress = NULL, ...) {
  expr   <- dataset$expr[, train_idx, drop = FALSE]
  labels <- droplevels(dataset$labels[train_idx])
  groups <- dataset$groups

  universe <- intersect(unique(unlist(groups)), rownames(expr))
  gene_p   <- .pre_score_pergene_anova(expr[universe, , drop = FALSE], labels)

  gnames <- names(groups)
  n <- length(gnames)
  pg <- numeric(n); names(pg) <- gnames
  for (j in seq_len(n)) {
    g  <- intersect(groups[[gnames[j]]], names(gene_p))
    if (length(g) == 0) { pg[j] <- 1 } else {
      gp <- gene_p[g]
      gp[gp <= 0] <- .Machine$double.xmin
      pg[j] <- .ebm_combine(gp, expr[g, , drop = FALSE])
    }
    if (!is.null(progress))
      progress(j, n, gnames[j], -log10(max(pg[j], .Machine$double.xmin)))
  }

  .pre_score_finalize(groups, pg, top_percent)
}


#' Pre-score groups using CAMERA (competitive, correlation-aware)
#'
#' Wraps `limma::camera` (Wu & Smyth, NAR 2012): a *competitive* gene-set test
#' that asks whether a group's genes are more differentially associated with the
#' class than the genes outside it, while adjusting for inter-gene correlation.
#' Unlike the Fisher-based methods this is not self-contained, so a large group
#' of weakly-associated genes no longer ranks highly by size alone. The per-group
#' p-values are ranked by `-log10(p)` and BH-adjusted (note `camera` also reports
#' its own FDR, but BH over the retained groups is attached for consistency).
#'
#' Two-class comparisons only; for more than two classes it falls back to
#' [pre_score_fisher_chisq()] with a warning.
#'
#' @param dataset A `gsm_dataset` object.
#' @param train_idx Integer vector of training sample indices.
#' @param top_percent Percentage of top groups to retain.
#' @param progress Optional callback `progress(i, n, group_name, score)`.
#' @param ... Ignored.
#' @return Named list of character vectors (top groups), with `scores`,
#'   `pvalues`, `padj` and `n_keep` attributes.
#' @keywords internal
pre_score_camera <- function(dataset, train_idx, top_percent = 20,
                             progress = NULL, ...) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Package 'limma' is required for method 'camera'. ",
         "Install it with: BiocManager::install('limma')")
  }
  expr   <- dataset$expr[, train_idx, drop = FALSE]
  labels <- droplevels(dataset$labels[train_idx])
  groups <- dataset$groups

  if (nlevels(labels) > 2) {
    warning("method 'camera' supports two-class comparisons; ",
            "falling back to 'fisher_chisq' for >2 classes.")
    return(pre_score_fisher_chisq(dataset, train_idx, top_percent = top_percent,
                                  progress = progress))
  }

  design <- stats::model.matrix(~ labels)
  index  <- lapply(groups, function(g) which(rownames(expr) %in% g))
  keep   <- vapply(index, length, integer(1)) > 0

  res <- limma::camera(expr, index[keep], design, contrast = ncol(design))

  # camera reorders by significance; map p-values back onto every group
  # (groups with no genes present default to p = 1).
  pg <- stats::setNames(rep(1, length(groups)), names(groups))
  pg[rownames(res)] <- res$PValue

  if (!is.null(progress)) {
    gnames <- names(groups); n <- length(gnames)
    for (j in seq_len(n))
      progress(j, n, gnames[j], -log10(max(pg[gnames[j]], .Machine$double.xmin)))
  }

  .pre_score_finalize(groups, pg, top_percent)
}


# ──────────────────────────────────────────────────────────────────────────────
# Tier 2: single-sample enrichment pre-scoring
# ──────────────────────────────────────────────────────────────────────────────

#' Single-sample gene-set enrichment matrix
#'
#' Summarises each gene set into one activity score per sample. Uses the
#' \pkg{GSVA} package (ssGSEA or GSVA) when available — handling both the modern
#' parameter-object API and the legacy `gsva(expr, sets, method=)` API — and
#' otherwise falls back to a dependency-free combined z-score (per-gene
#' standardisation across samples, averaged over the set; Lee et al. 2008).
#'
#' @param expr Numeric matrix, genes (rows) x samples (columns), training only.
#' @param gene_sets Named list of character vectors (genes per set); entries are
#'   already restricted to genes present in `expr`.
#' @param method One of `"ssgsea"`, `"gsva"` (require \pkg{GSVA}) or `"zscore"`.
#' @return Numeric matrix, sets (rows) x samples (columns).
#' @keywords internal
.enrichment_matrix <- function(expr, gene_sets, method = "ssgsea") {
  gene_sets <- gene_sets[vapply(gene_sets, length, integer(1)) > 0]
  if (length(gene_sets) == 0) return(matrix(numeric(0), nrow = 0, ncol = ncol(expr)))

  if (method %in% c("ssgsea", "gsva")) {
    es <- tryCatch(
      suppressMessages(suppressWarnings(.run_gsva(expr, gene_sets, method))),
      error = function(e) NULL
    )
    if (!is.null(es)) return(as.matrix(es))
    warning("GSVA enrichment failed or unavailable; using built-in z-score.")
  }

  # Combined z-score fallback: standardise each gene, average over the set.
  mu  <- rowMeans(expr)
  sds <- sqrt(rowSums((expr - mu)^2) / max(1, ncol(expr) - 1))
  sds[!is.finite(sds) | sds == 0] <- 1
  z <- (expr - mu) / sds
  out <- t(vapply(gene_sets, function(g) colMeans(z[g, , drop = FALSE]),
                  numeric(ncol(expr))))
  rownames(out) <- names(gene_sets)
  out
}

#' Run GSVA across its modern and legacy APIs
#' @keywords internal
.run_gsva <- function(expr, gene_sets, method) {
  ns <- asNamespace("GSVA")
  if (exists("ssgseaParam", envir = ns)) {
    par <- if (method == "ssgsea") GSVA::ssgseaParam(expr, gene_sets)
           else                    GSVA::gsvaParam(expr, gene_sets)
    GSVA::gsva(par, verbose = FALSE)
  } else {
    GSVA::gsva(expr, gene_sets, method = method, verbose = FALSE)
  }
}

#' Two-sided p-value for a group-activity vector separating the classes
#'
#' Welch t-test for two classes, one-way ANOVA for more than two. Returns 1 when
#' the test is undefined (constant vector, a class with <2 samples, etc.).
#' @keywords internal
.activity_pvalue <- function(v, labels) {
  labels <- droplevels(labels)
  if (nlevels(labels) < 2 || any(table(labels) < 2)) return(1)
  p <- tryCatch({
    if (nlevels(labels) == 2) stats::t.test(v ~ labels)$p.value
    else stats::oneway.test(v ~ labels)$p.value
  }, error = function(e) 1)
  if (!is.finite(p)) 1 else p
}

#' Pre-score groups using single-sample enrichment (GSVA / ssGSEA)
#'
#' Transforms the training expression into a per-sample, per-group activity
#' matrix (see [.enrichment_matrix()]) and scores each group by how strongly its
#' activity separates the classes (Welch t / one-way ANOVA p-value). Unlike the
#' differential-expression methods, this evaluates the group as a single derived
#' feature — exactly the signal the downstream ML classifier consumes — so a
#' group only ranks highly when its aggregate behaviour is discriminative.
#'
#' @param dataset A `gsm_dataset` object.
#' @param train_idx Integer vector of training sample indices.
#' @param top_percent Percentage of top groups to retain.
#' @param progress Optional callback `progress(i, n, group_name, score)`.
#' @param enrichment Enrichment backend: `"ssgsea"`, `"gsva"` or `"zscore"`.
#'   Defaults to `"ssgsea"` when \pkg{GSVA} is installed, else `"zscore"`.
#' @param ... Ignored.
#' @return Named list of character vectors (top groups), with `scores`,
#'   `pvalues`, `padj` and `n_keep` attributes.
#' @keywords internal
pre_score_gsva <- function(dataset, train_idx, top_percent = 20,
                           progress = NULL, enrichment = NULL, ...) {
  expr   <- dataset$expr[, train_idx, drop = FALSE]
  labels <- droplevels(dataset$labels[train_idx])
  groups <- dataset$groups

  if (is.null(enrichment))
    enrichment <- if (requireNamespace("GSVA", quietly = TRUE)) "ssgsea" else "zscore"

  gl <- lapply(groups, function(g) intersect(g, rownames(expr)))
  es <- .enrichment_matrix(expr, gl, enrichment)

  gnames <- names(groups)
  n <- length(gnames)
  pg <- stats::setNames(rep(1, n), gnames)
  for (j in seq_len(n)) {
    nm <- gnames[j]
    if (nm %in% rownames(es)) pg[j] <- .activity_pvalue(es[nm, ], labels)
    if (!is.null(progress))
      progress(j, n, nm, -log10(max(pg[j], .Machine$double.xmin)))
  }

  .pre_score_finalize(groups, pg, top_percent)
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
