#' @title Label-blind expression filtering (denoising)
#' @description Reduce an expression matrix by dropping uninformative genes
#'   *without ever looking at the class labels*, so it can be applied once, on
#'   all samples, before any grouping / cross-validation without leaking class
#'   information. Three complementary, order-fixed steps: remove missing/constant
#'   genes, then low-expression genes, then keep the most variable genes.
#' @name filter
NULL

#' Filter genes from an expression matrix (label-blind)
#'
#' Applies, in order and only where the corresponding argument is active:
#' \enumerate{
#'   \item **missing / constant** — drop genes with more than `na_frac` missing
#'     values, or (when `drop_constant`) zero variance across samples;
#'   \item **low-expression** — drop genes whose mean is below `min_mean`;
#'   \item **low-variance** — keep the `max_genes` most variable genes, and/or
#'     those with variance at least `min_var`.
#' }
#' None of these inspect class labels, so the filter is leakage-free.
#'
#' @param expr Numeric matrix with genes as rows and samples as columns; row
#'   names (gene IDs) required.
#' @param drop_constant Logical; drop zero-variance genes. Default `TRUE`.
#' @param na_frac Drop genes with a greater fraction of missing values than this.
#'   Default `0.5`. Set `NULL` to skip.
#' @param min_mean Optional numeric; drop genes with row mean below this.
#'   `NULL` (default) skips the low-expression step.
#' @param min_var Optional numeric; drop genes with variance below this.
#'   `NULL` (default) skips this threshold.
#' @param max_genes Optional integer; keep at most this many genes, the highest
#'   by variance. `NULL` (default) = no cap.
#' @return The filtered matrix (a subset of the rows of `expr`), with an
#'   attribute `"filter_report"`: a named integer vector of how many genes
#'   remained after each step (`input`, `after_missing`, `after_lowexpr`,
#'   `after_variance`).
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(2000), nrow = 200,
#'             dimnames = list(paste0("G", 1:200), paste0("S", 1:10)))
#' X[1:20, ] <- 0            # constant genes
#' Xf <- filter_genes(X, max_genes = 100)
#' attr(Xf, "filter_report")
#' @export
filter_genes <- function(expr, drop_constant = TRUE, na_frac = 0.5,
                         min_mean = NULL, min_var = NULL, max_genes = NULL) {
  if (!is.matrix(expr) || !is.numeric(expr))
    stop("`expr` must be a numeric matrix (genes x samples).", call. = FALSE)
  if (is.null(rownames(expr)))
    stop("`expr` must have row names (gene IDs).", call. = FALSE)

  report <- c(input = nrow(expr))
  keep <- rep(TRUE, nrow(expr))

  rvar <- function(m) {
    n <- ncol(m); if (n < 2) return(rep(NA_real_, nrow(m)))
    mu <- rowMeans(m, na.rm = TRUE)
    rowSums((m - mu)^2, na.rm = TRUE) / (n - 1)
  }

  # 1) missing / constant
  if (!is.null(na_frac)) {
    na_share <- rowMeans(is.na(expr))
    keep <- keep & (na_share <= na_frac)
  }
  if (isTRUE(drop_constant)) {
    v <- rvar(expr)
    keep <- keep & is.finite(v) & v > 0
  }
  report["after_missing"] <- sum(keep)

  # 2) low-expression
  if (!is.null(min_mean)) {
    mu <- rowMeans(expr, na.rm = TRUE)
    keep <- keep & is.finite(mu) & (mu >= min_mean)
  }
  report["after_lowexpr"] <- sum(keep)

  # 3) low-variance threshold and/or top-N by variance
  if (!is.null(min_var) || !is.null(max_genes)) {
    v <- rvar(expr)
    v[!keep] <- -Inf                       # already-dropped genes can't be kept
    if (!is.null(min_var)) keep <- keep & is.finite(v) & (v >= min_var)
    if (!is.null(max_genes) && sum(keep) > max_genes) {
      cutoff_ranks <- order(v, decreasing = TRUE)[seq_len(max_genes)]
      top <- rep(FALSE, nrow(expr)); top[cutoff_ranks] <- TRUE
      keep <- keep & top
    }
  }
  report["after_variance"] <- sum(keep)

  if (sum(keep) < 2)
    stop("Filtering removed almost all genes (", sum(keep), " left). ",
         "Relax the thresholds.", call. = FALSE)

  out <- expr[keep, , drop = FALSE]
  attr(out, "filter_report") <- report
  out
}
