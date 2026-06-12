#' @title Consensus biomarker (gene) ranking from MCCV
#' @description Turns the per-fold, per-gene model importances produced by
#'   [run_mccv()] into a single stability-weighted biomarker list. A gene that is
#'   both frequently selected and consistently important across folds rises to the
#'   top, which is the gene-level claim a biologist ultimately wants from the
#'   group-level G-S-M result.
#' @name biomarkers
NULL

#' Aggregate per-gene importance across MCCV folds into a biomarker table
#'
#' Combines two independent signals across folds: how often a gene appears in the
#' best model (selection frequency) and how important it is when it does (mean
#' importance), plus a Robust Rank Aggregation consensus of the per-fold
#' importance rankings (see [rra_aggregate()]). Requires a scoring/modeling method
#' that emits per-gene importance — Random Forest (`MeanDecreaseGini`), linear-
#' kernel SVM (absolute weights), or XGBoost (gain). Returns `NULL` with a message
#' when no fold provided importances (e.g. KNN or radial SVM models).
#'
#' @param mccv_results The list returned by [run_mccv()] (must contain `$folds`).
#' @param top_n Optional integer; if set, keep only the top `top_n` genes.
#' @return A data.frame sorted by consensus strength, or `NULL`:
#'   \describe{
#'     \item{Gene}{Gene identifier.}
#'     \item{Frequency}{Fraction of successful folds whose best model used the gene.}
#'     \item{MeanImportance}{Mean importance over the folds where the gene appeared.}
#'     \item{SDImportance}{SD of that importance (NA when seen in one fold).}
#'     \item{RRA_Score}{Robust Rank Aggregation score (lower = stronger consensus).}
#'     \item{Rank}{Final rank (by RRA when available, else by frequency then importance).}
#'   }
#' @export
gsm_biomarkers <- function(mccv_results, top_n = NULL) {
  folds <- mccv_results$folds %||% mccv_results
  ok <- Filter(function(r) !is.null(r) && identical(r$status, "success"), folds)
  imp_lists <- Filter(function(x) !is.null(x) && length(x) > 0,
                      lapply(ok, `[[`, "importance"))
  if (length(imp_lists) == 0) {
    message("No per-gene importances available (need RF, linear SVM, or XGBoost ",
            "scoring/modeling). Returning NULL.")
    return(NULL)
  }
  n_folds <- length(ok)

  # Per-gene appearance count and importance values across folds.
  genes <- unique(unlist(lapply(imp_lists, names)))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  vals <- lapply(genes, function(g) {
    # Single-bracket by name yields NA (not an error) when the gene is absent
    # from a fold's importance vector.
    v <- vapply(imp_lists, function(imp) as.numeric(imp[g]), numeric(1))
    v[is.finite(v)]
  })
  names(vals) <- genes

  freq <- vapply(vals, length, integer(1)) / n_folds
  mean_imp <- vapply(vals, function(v) if (length(v)) mean(v) else NA_real_, numeric(1))
  sd_imp <- vapply(vals, function(v) if (length(v) > 1) stats::sd(v) else NA_real_, numeric(1))

  out <- data.frame(
    Gene = genes,
    Frequency = round(as.numeric(freq), 4),
    MeanImportance = signif(as.numeric(mean_imp), 4),
    SDImportance = signif(as.numeric(sd_imp), 4),
    stringsAsFactors = FALSE)

  # Consensus rank: prefer RRA over the per-fold importance rankings; fall back
  # to selection frequency then mean importance when RRA is unavailable.
  rra <- gsm_tbl_feature_rra(ok)
  if (!is.null(rra)) {
    out$RRA_Score <- rra$RRA_Score[match(out$Gene, rra$Gene)]
    out <- out[order(out$RRA_Score, -out$Frequency, -out$MeanImportance,
                     na.last = TRUE), , drop = FALSE]
  } else {
    out$RRA_Score <- NA_real_
    out <- out[order(-out$Frequency, -out$MeanImportance), , drop = FALSE]
  }
  out$Rank <- seq_len(nrow(out))
  rownames(out) <- NULL

  if (!is.null(top_n)) out <- utils::head(out, top_n)
  out
}
