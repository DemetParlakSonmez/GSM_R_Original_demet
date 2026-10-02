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
#' @param alpha BH-FDR threshold (default `0.05`). Used when `method = "fdr"`.
#' @param method Selection rule for surviving genes:
#'   * `"fdr"` (default) — keep genes with BH-adjusted p-value below `alpha`.
#'     The number of survivors varies from fold to fold with the data.
#'   * `"topk"` — keep the `top_k` genes with the strongest evidence (smallest
#'     raw p-value, ties broken by larger `|t|`), giving a **fixed** count per
#'     fold regardless of significance.
#' @param top_k Number of genes to keep when `method = "topk"`. Default `100`.
#' @return Character vector of surviving gene names. If fewer than two genes
#'   pass FDR (or a class has <2 training samples / data is not two-class), the
#'   function falls back to returning all genes so the pipeline can continue.
#' @export
shave_genes <- function(expr, labels, train_idx, alpha = 0.05,
                        method = c("fdr", "topk"), top_k = 100) {
  method <- match.arg(method)
  genes <- rownames(expr)
  wt <- welch_t_vector(expr, labels, train_idx)
  if (is.null(wt)) {
    return(genes)
  } # not two-class / too few samples

  p <- 2 * stats::pt(-abs(wt$t), df = wt$df)
  p[!is.finite(wt$t) | !is.finite(wt$df)] <- 1 # constant genes -> not significant

  if (method == "topk") {
    k <- max(1L, min(as.integer(top_k), length(genes)))
    ord <- order(p, -abs(wt$t)) # smallest p first, |t| tie-break
    return(genes[ord[seq_len(k)]])
  }

  p_adj <- stats::p.adjust(p, method = "BH")

  sig_idx <- which(is.finite(p_adj) & p_adj < alpha)

  keep <- genes[sig_idx]
  if (length(keep) < 2) genes else keep # safe fallback
}

# Vectorised two-class Welch t-test on the training columns.
# Returns list(t, df) named by gene, or NULL when not two-class / <2 per class.
welch_t_vector <- function(expr, labels, train_idx, genes = NULL) {
  if (!is.null(genes)) {
    genes <- intersect(genes, rownames(expr))
    if (length(genes) == 0) {
      return(NULL)
    }
    expr <- expr[genes, , drop = FALSE]
  }
  y <- droplevels(labels[train_idx])
  if (nlevels(y) != 2) {
    return(NULL)
  }
  m <- expr[, train_idx, drop = FALSE]
  lev <- levels(y)
  a <- m[, y == lev[1], drop = FALSE]
  b <- m[, y == lev[2], drop = FALSE]
  na <- ncol(a)
  nb <- ncol(b)
  if (na < 2 || nb < 2) {
    return(NULL)
  }
  ma <- rowMeans(a)
  mb <- rowMeans(b)
  va <- matrixStats_rowVars(a)
  vb <- matrixStats_rowVars(b)
  se2 <- va / na + vb / nb
  t <- (ma - mb) / sqrt(se2) # signed (class1 - class2)
  df <- se2^2 / ((va / na)^2 / (na - 1) + (vb / nb)^2 / (nb - 1))
  list(t = stats::setNames(t, rownames(m)), df = df)
}

#' Signed per-gene training statistic for direction-aware discovery
#'
#' Named signed Welch t-statistic (class1 − class2) over a bounded gene
#' universe, computed on the training partition only. Drives the UP/DOWN/ABS
#' modes of [gsm_detect()].
#'
#' @param expr Numeric matrix (genes x samples).
#' @param labels Factor of class labels (two levels used).
#' @param train_idx Integer training sample indices.
#' @param universe Optional character vector restricting the genes.
#' @return Named numeric signed t-statistic, or `NULL` if not two-class.
#' @export
gene_signed_stats <- function(expr, labels, train_idx, universe = NULL) {
  wt <- welch_t_vector(expr, labels, train_idx, genes = universe)
  if (is.null(wt)) {
    return(NULL)
  }
  wt$t[is.finite(wt$t)]
}

# Per-row variance without a hard matrixStats dependency.
matrixStats_rowVars <- function(x) {
  n <- ncol(x)
  if (n < 2) {
    return(rep(NA_real_, nrow(x)))
  }
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
