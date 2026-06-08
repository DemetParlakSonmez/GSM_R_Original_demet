#' @title Robust Rank Aggregation (RRA)
#' @description Aggregate several per-iteration ranked lists (groups by F1, or
#'   genes by importance) into one consensus ranking. Uses the
#'   \pkg{RobustRankAggreg} package (Kolde et al. 2012) when available, and
#'   otherwise falls back to mean normalized rank so the pipeline still runs.
#' @name rra
NULL

#' Aggregate ranked lists into a consensus ranking
#'
#' @param rank_lists A list of character vectors, each ordered best-first
#'   (one per MCCV iteration).
#' @return A data.frame with `Name`, `RRA_Score` (lower = stronger consensus;
#'   a p-value-like score from RRA, or mean normalized rank in the fallback),
#'   `Rank`, and `Method` ("RRA" or "mean_rank").
#' @export
rra_aggregate <- function(rank_lists) {
  rank_lists <- Filter(function(x) length(x) > 0, rank_lists)
  if (length(rank_lists) == 0) return(NULL)

  if (requireNamespace("RobustRankAggreg", quietly = TRUE)) {
    N <- length(unique(unlist(rank_lists)))
    res <- tryCatch(
      RobustRankAggreg::aggregateRanks(glist = rank_lists, N = N,
                                       method = "RRA"),
      error = function(e) NULL)
    if (!is.null(res)) {
      out <- data.frame(Name = as.character(res$Name),
                        RRA_Score = signif(as.numeric(res$Score), 3),
                        stringsAsFactors = FALSE)
      out <- out[order(out$RRA_Score), , drop = FALSE]
      out$Rank <- seq_len(nrow(out))
      out$Method <- "RRA"
      return(out)
    }
  }

  # Fallback: mean of normalized ranks (rank/length) across lists.
  items <- unique(unlist(rank_lists))
  mean_nr <- vapply(items, function(it) {
    rs <- vapply(rank_lists, function(L) {
      p <- match(it, L)
      if (is.na(p)) 1 else p / length(L)          # absent -> worst (1.0)
    }, numeric(1))
    mean(rs)
  }, numeric(1))
  out <- data.frame(Name = items, RRA_Score = round(as.numeric(mean_nr), 6),
                    stringsAsFactors = FALSE)
  out <- out[order(out$RRA_Score), , drop = FALSE]
  out$Rank <- seq_len(nrow(out))
  out$Method <- "mean_rank"
  out
}

# Successful folds only (mirrors helper in export.R; kept local for use here).
.ok_folds_rra <- function(folds) {
  Filter(function(r) !is.null(r) && identical(r$status, "success"), folds)
}

#' RRA consensus of per-iteration GROUP rankings
#' @param folds List of fold results (each with `ranked_groups`).
#' @return data.frame (Name=Group, RRA_Score, Rank, Method) or NULL.
#' @keywords internal
gsm_tbl_group_rra <- function(folds) {
  ok <- .ok_folds_rra(folds)
  lists <- Filter(Negate(is.null), lapply(ok, `[[`, "ranked_groups"))
  res <- rra_aggregate(lists)
  if (is.null(res)) return(NULL)
  names(res)[1] <- "Group"
  res
}

#' RRA consensus of per-iteration FEATURE (gene) rankings by model importance
#' @param folds List of fold results (each with `importance`).
#' @return data.frame (Name=Gene, RRA_Score, Rank, Method) or NULL (when no
#'   fold provided importances, e.g. SVM/KNN scoring).
#' @keywords internal
gsm_tbl_feature_rra <- function(folds) {
  ok <- .ok_folds_rra(folds)
  lists <- lapply(ok, function(f) {
    imp <- f$importance
    if (is.null(imp) || length(imp) == 0) return(NULL)
    names(sort(imp, decreasing = TRUE))           # genes best-first by importance
  })
  lists <- Filter(Negate(is.null), lists)
  res <- rra_aggregate(lists)
  if (is.null(res)) return(NULL)
  names(res)[1] <- "Gene"
  res
}
