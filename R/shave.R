#' @title Training-only gene shaving (Welch t-test + BH FDR)
#' @description Phase-I statistical gene filter of the G-S-M framework: on the
#'   training partition only, run a per-gene Welch two-sample t-test, correct the
#'   p-values with Benjamini-Hochberg (FDR), and keep genes below the threshold.
#'   Surviving genes are then used for grouping, scoring, and modeling.
#' @name shave
NULL

#' Shave genes by a training-only Welch t-test with BH-FDR control
#'
#' Vectorised two-class Welch t-test computed on the training samples, with
#' Benjamini-Hochberg multiple-testing correction. Genes with adjusted p-value
#' below `alpha` are retained. Used once per train/test split before grouping.
#'
#' @param expr Numeric matrix (genes x samples).
#' @param labels Factor of class labels (exactly two levels used).
#' @param train_idx Integer vector of training sample indices.
#' @param alpha BH-FDR threshold (default `0.05`).
#' @return Character vector of surviving gene names. If fewer than two genes
#'   pass (or a class has <2 training samples), falls back to returning all
#'   genes so the pipeline can continue.
#' @export
shave_genes <- function(expr, labels, train_idx, alpha = 0.05) {
  genes <- rownames(expr)
  y <- droplevels(labels[train_idx])
  if (nlevels(y) != 2) return(genes)             # t-test needs two classes

  m <- expr[, train_idx, drop = FALSE]
  lev <- levels(y)
  a <- m[, y == lev[1], drop = FALSE]
  b <- m[, y == lev[2], drop = FALSE]
  na <- ncol(a); nb <- ncol(b)
  if (na < 2 || nb < 2) return(genes)

  ma <- rowMeans(a); mb <- rowMeans(b)
  va <- matrixStats_rowVars(a); vb <- matrixStats_rowVars(b)
  se2 <- va / na + vb / nb
  t <- (ma - mb) / sqrt(se2)
  # Welch-Satterthwaite degrees of freedom (per gene).
  df <- se2^2 / ((va / na)^2 / (na - 1) + (vb / nb)^2 / (nb - 1))
  p <- 2 * stats::pt(-abs(t), df = df)
  p[!is.finite(t) | !is.finite(df)] <- 1        # constant genes -> not significant

  p_adj <- stats::p.adjust(p, method = "BH")
  keep <- genes[is.finite(p_adj) & p_adj < alpha]
  if (length(keep) < 2) genes else keep          # safe fallback
}

# Per-row variance without a hard matrixStats dependency.
matrixStats_rowVars <- function(x) {
  n <- ncol(x)
  if (n < 2) return(rep(NA_real_, nrow(x)))
  mu <- rowMeans(x)
  rowSums((x - mu)^2) / (n - 1)
}

#' Intersect each group's genes with a surviving-gene set, dropping empties
#'
#' @param groups Named list of gene character vectors.
#' @param keep_genes Character vector of surviving gene names.
#' @return Named list of (non-empty) shaved groups.
#' @export
shave_groups <- function(groups, keep_genes) {
  out <- lapply(groups, function(g) intersect(g, keep_genes))
  out[vapply(out, length, integer(1)) > 0]
}
