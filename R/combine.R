#' @title Combined cross-dataset report
#' @description Aggregate the per-dataset results from [run_batch()] into a
#'   single report: a per-dataset performance table plus consensus group and
#'   gene rankings across datasets (via Robust Rank Aggregation, the same engine
#'   used across folds within a dataset — it is dataset-agnostic).
#' @name combine
NULL

#' Per-dataset performance summary
#'
#' One row per dataset (mean/SD of each metric across its MCCV folds) plus a
#' final `ALL (mean)` meta-row averaging the per-dataset means.
#'
#' @param results A named list of per-dataset result objects, each with `name`,
#'   `metrics_summary` (a data.frame with metric rownames and `Mean`/`SD`
#'   columns) and `n_ok` — as produced by [run_batch()].
#' @return A data.frame, or `NULL` when `results` is empty.
#' @export
gsm_combined_metrics <- function(results) {
  if (length(results) == 0) return(NULL)
  metric_names <- unique(unlist(lapply(results, function(r) rownames(r$metrics_summary))))
  if (is.null(metric_names)) return(NULL)

  rows <- lapply(results, function(r) {
    mm <- r$metrics_summary
    row <- data.frame(
      Dataset      = r$name,
      Samples      = r$n_samples %||% NA_integer_,
      Genes        = r$n_genes   %||% NA_integer_,
      Groups       = r$n_groups  %||% NA_integer_,
      N_folds_ok   = r$n_ok      %||% NA_integer_,
      Mean_Best_k  = round(r$best_k_mean %||% NA_real_, 2),
      Elapsed_sec  = round(r$elapsed     %||% NA_real_, 1),
      stringsAsFactors = FALSE)
    for (m in metric_names) {
      row[[paste0("Mean_", m)]] <- round(tryCatch(mm[m, "Mean"], error = function(e) NA_real_), 4)
      row[[paste0("SD_", m)]]   <- round(tryCatch(mm[m, "SD"],   error = function(e) NA_real_), 4)
    }
    row
  })
  df <- do.call(rbind, rows)
  rownames(df) <- NULL

  # Meta row: sum the counts/time, average the rates.
  meta <- df[1, , drop = FALSE]
  meta$Dataset <- "ALL"
  for (cn in c("N_folds_ok", "Elapsed_sec"))
    meta[[cn]] <- sum(df[[cn]], na.rm = TRUE)
  for (cn in c("Samples", "Genes", "Groups", "Mean_Best_k",
               grep("^(Mean|SD)_", names(df), value = TRUE)))
    meta[[cn]] <- round(mean(df[[cn]], na.rm = TRUE), 4)
  rbind(df, meta)
}

#' Build a consensus ranking across datasets from per-dataset ranked lists
#' @keywords internal
.gsm_combine_consensus <- function(results, extract, item_col) {
  lists <- lapply(results, function(r) tryCatch(extract(r), error = function(e) NULL))
  lists <- Filter(function(x) !is.null(x) && length(x) > 0, lists)
  if (length(lists) == 0) return(NULL)
  agg <- rra_aggregate(lists)
  if (is.null(agg)) return(NULL)
  names(agg)[1] <- item_col

  freq <- table(unlist(lapply(lists, unique)))
  agg$Dataset_Frequency <- as.integer(freq[agg[[item_col]]])
  agg$Dataset_Frequency[is.na(agg$Dataset_Frequency)] <- 0L
  agg$N_Datasets <- length(lists)
  agg[order(agg$Rank), , drop = FALSE]
}

#' Consensus group ranking across datasets
#'
#' Aggregates each dataset's per-fold group consensus ([gsm_tbl_group_rra()])
#' into one cross-dataset ranking with `Dataset_Frequency` (how many datasets
#' ranked the group).
#'
#' @param results Per-dataset results from [run_batch()] (each with `mccv$folds`).
#' @return A data.frame (`Group`, `RRA_Score`, `Rank`, `Method`,
#'   `Dataset_Frequency`, `N_Datasets`) or `NULL`.
#' @export
gsm_combined_group_consensus <- function(results) {
  .gsm_combine_consensus(results, function(r) {
    tb <- gsm_tbl_group_rra(r$mccv$folds)
    if (is.null(tb)) return(NULL)
    tb$Group[order(tb$Rank)]
  }, "Group")
}

#' Consensus biomarkers (genes) across datasets
#'
#' Aggregates each dataset's [gsm_biomarkers()] gene ranking into one
#' cross-dataset ranking. Datasets whose scoring yields no per-gene importance
#' (e.g. radial SVM / KNN) contribute nothing and are skipped.
#'
#' @param results Per-dataset results from [run_batch()] (each with `mccv`).
#' @return A data.frame (`Gene`, `RRA_Score`, `Rank`, `Method`,
#'   `Dataset_Frequency`, `N_Datasets`) or `NULL`.
#' @export
gsm_combined_biomarkers <- function(results) {
  .gsm_combine_consensus(results, function(r) {
    bm <- gsm_biomarkers(r$mccv)
    if (is.null(bm) || nrow(bm) == 0) return(NULL)
    bm$Gene[order(bm$Rank)]
  }, "Gene")
}

#' Build and write the combined cross-dataset report
#'
#' Bundles the per-dataset performance table, consensus group ranking, consensus
#' biomarkers, and (optionally) the run manifest into a multi-sheet `.xlsx`
#' workbook or a zip of CSVs, reusing the single-dataset writers
#' ([gsm_write_xlsx()] / [gsm_write_csv_zip()]).
#'
#' @param results Per-dataset results from [run_batch()].
#' @param file Output path (`.xlsx` or `.zip`).
#' @param format `"xlsx"` or `"csv"`.
#' @param manifest Optional manifest data.frame to include as a sheet.
#' @param run_info Optional key/value data.frame of run parameters + timing
#'   (from `gsm_batch_run_info()`), included as the first sheet.
#' @return `file`, invisibly.
#' @export
gsm_export_combined_report <- function(results, file, format = c("xlsx", "csv"),
                                       manifest = NULL, run_info = NULL) {
  format <- match.arg(format)
  tables <- list(
    RunInfo            = run_info,
    PerDatasetMetrics  = gsm_combined_metrics(results),
    GroupConsensus     = gsm_combined_group_consensus(results),
    BiomarkerConsensus = gsm_combined_biomarkers(results)
  )
  if (!is.null(manifest) && nrow(manifest) > 0) tables$Manifest <- manifest
  tables <- Filter(function(t) !is.null(t) && nrow(t) > 0, tables)
  if (length(tables) == 0) stop("No combined tables to write (no successful datasets).")

  if (format == "xlsx") gsm_write_xlsx(tables, file)
  else gsm_write_csv_zip(tables, file, prefix = "Combined")
  invisible(file)
}
