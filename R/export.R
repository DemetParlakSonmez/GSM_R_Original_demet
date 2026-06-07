#' @title Export G-S-M scoring results
#' @description Build tidy result tables from a scoring run and write them to
#'   Excel (multi-sheet) or a CSV bundle (zip). Mirrors the output files of the
#'   original G-S-M framework: per-iteration results, ranked groups by
#'   frequency, group-feature associations over all iterations, ranked
#'   features, and a run-parameters summary.
#' @name export
NULL

# Internal null-coalescing helper (base R only ships %||% from 4.4.0).
`%|null|%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# Friendly classifier names (kept local to avoid cross-file coupling).
.clf_label <- function(m) {
  switch(as.character(m), rf = "RandomForest", svm = "SVM", knn = "KNN",
         as.character(m))
}

# Successful folds only.
.ok_folds <- function(folds) {
  Filter(function(r) !is.null(r) && identical(r$status, "success"), folds)
}

# Number of genes in a group (0 if unknown).
.n_genes <- function(dataset, g) {
  length(dataset$groups[[g]] %|null|% character(0))
}

# ── MCCV table builders ──────────────────────────────────────────────────────

#' Per-iteration, per-group scores (mirrors the original `res2` group sheet).
#' @keywords internal
gsm_tbl_iteration_group_scores <- function(folds, dataset) {
  ok <- .ok_folds(folds)
  rows <- lapply(ok, function(f) {
    accs <- f$group_accuracies
    if (length(accs) == 0) return(NULL)
    g <- names(accs)
    data.frame(
      Iteration    = f$fold,
      Group        = g,
      Accuracy     = round(as.numeric(accs), 6),
      Rank_in_fold = rank(-as.numeric(accs), ties.method = "min"),
      N_genes      = vapply(g, function(x) .n_genes(dataset, x), integer(1)),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Per-iteration final-model metrics (mirrors the original `res2` model sheet).
#' Error folds are kept (Status = "error", NA metrics) so nothing is dropped.
#' @keywords internal
gsm_tbl_iteration_metrics <- function(folds) {
  getm <- function(m, k) if (is.null(m[[k]])) NA_real_ else round(as.numeric(m[[k]]), 6)
  rows <- lapply(folds, function(f) {
    if (is.null(f)) return(NULL)
    ok <- identical(f$status, "success")
    m <- if (ok) f$metrics else list()
    data.frame(
      Iteration    = f$fold,
      Status       = f$status %|null|% "error",
      Accuracy     = getm(m, "Accuracy"),
      Sensitivity  = getm(m, "Sensitivity"),
      Specificity  = getm(m, "Specificity"),
      Precision    = getm(m, "Precision"),
      F_measure    = getm(m, "F_measure"),
      Cohens_kappa = getm(m, "Cohens_kappa"),
      Features     = if (ok) f$features_count %|null|% NA_integer_ else NA_integer_,
      Groups       = if (ok) f$n_groups %|null|% NA_integer_ else NA_integer_,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out[order(out$Iteration), , drop = FALSE]
}

#' Group ranking / stability across all iterations.
#' Mirrors `RankdGroups_by_Freq` + `ranked_Groups` + the top-N of `res1`.
#' @keywords internal
gsm_tbl_group_ranking <- function(folds, dataset, top_n = 10) {
  ok <- .ok_folds(folds)
  n_ok <- length(ok)
  long <- gsm_tbl_iteration_group_scores(folds, dataset)
  if (is.null(long) || nrow(long) == 0) return(NULL)

  groups <- sort(unique(long$Group))
  agg <- lapply(groups, function(g) {
    sub <- long[long$Group == g, , drop = FALSE]
    data.frame(
      Group           = g,
      Times_Selected  = nrow(sub),
      Selection_Rate  = round(100 * nrow(sub) / max(1, n_ok), 1),
      TopN_Appearances = sum(sub$Rank_in_fold <= top_n),
      TopN_Rate       = round(100 * sum(sub$Rank_in_fold <= top_n) / max(1, n_ok), 1),
      Mean_Accuracy   = round(mean(sub$Accuracy), 4),
      SD_Accuracy     = round(stats::sd(sub$Accuracy), 4),
      Mean_Rank       = round(mean(sub$Rank_in_fold), 2),
      N_genes         = sub$N_genes[1],
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, agg)
  out[order(-out$TopN_Appearances, -out$Mean_Accuracy, out$Mean_Rank), , drop = FALSE]
}

#' Aggregate metric summary (mirrors `res1`). Uses `result$metrics_summary`
#' when present, otherwise computes from folds.
#' @keywords internal
gsm_tbl_metric_summary <- function(result) {
  ms <- result$metrics_summary
  if (!is.null(ms)) {
    return(data.frame(Metric = rownames(ms),
                      Mean = round(ms$Mean, 6), SD = round(ms$SD, 6),
                      stringsAsFactors = FALSE))
  }
  im <- gsm_tbl_iteration_metrics(result$folds)
  cols <- c("Accuracy", "Sensitivity", "Specificity", "Precision",
            "F_measure", "Cohens_kappa")
  data.frame(
    Metric = cols,
    Mean = round(vapply(cols, function(c) mean(im[[c]], na.rm = TRUE), numeric(1)), 6),
    SD   = round(vapply(cols, function(c) stats::sd(im[[c]], na.rm = TRUE), numeric(1)), 6),
    stringsAsFactors = FALSE
  )
}

#' Group <-> feature associations over all iterations
#' (mirrors `Groups_Features_Associations_over_all_iterations.csv`).
#' @keywords internal
gsm_tbl_group_feature_assoc <- function(folds, dataset) {
  ok <- .ok_folds(folds)
  sel <- table(unlist(lapply(ok, function(f) names(f$group_accuracies))))
  if (length(sel) == 0) return(NULL)
  rows <- lapply(names(sel), function(g) {
    genes <- dataset$groups[[g]] %|null|% character(0)
    if (length(genes) == 0) return(NULL)
    data.frame(Group = g, Gene = genes,
               Times_Group_Selected = as.integer(sel[[g]]),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

#' Per-iteration pre-scoring significance (mirrors the KNIME significance file).
#' Returns NULL when no pre-scores were captured (e.g. method = "none").
#' @keywords internal
gsm_tbl_pre_score_significance <- function(folds) {
  ok <- .ok_folds(folds)
  rows <- lapply(ok, function(f) {
    ps <- f$pre_scores
    if (is.null(ps) || length(ps) == 0) return(NULL)
    data.frame(Iteration = f$fold, Group = names(ps),
               Fisher_Score = round(as.numeric(ps), 4),
               Pre_Rank = seq_along(ps), stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  if (is.null(out) || nrow(out) == 0) NULL else out
}

#' Ranked features over all iterations (mirrors `rankedFeatures`).
#' Uses true RandomForest importance (MeanDecreaseGini) when available;
#' always includes selection-frequency context.
#' @keywords internal
gsm_tbl_ranked_features <- function(folds, dataset) {
  ok <- .ok_folds(folds)
  if (length(ok) == 0) return(NULL)

  # Long (gene, group, accuracy) over every selected group in every fold.
  rows <- lapply(ok, function(f) {
    accs <- f$group_accuracies
    if (length(accs) == 0) return(NULL)
    do.call(rbind, lapply(names(accs), function(g) {
      genes <- dataset$groups[[g]] %|null|% character(0)
      if (length(genes) == 0) return(NULL)
      data.frame(Gene = genes, Group = g,
                 Accuracy = as.numeric(accs[[g]]), stringsAsFactors = FALSE)
    }))
  })
  long <- do.call(rbind, rows)
  if (is.null(long) || nrow(long) == 0) return(NULL)

  genes <- sort(unique(long$Gene))
  freq <- as.integer(table(long$Gene)[genes])
  ngrp <- vapply(genes, function(g)
    length(unique(long$Group[long$Gene == g])), integer(1))
  macc <- vapply(genes, function(g)
    mean(long$Accuracy[long$Gene == g]), numeric(1))

  # Mean RF importance across folds that provided it.
  imp_list <- Filter(Negate(is.null), lapply(ok, `[[`, "importance"))
  mean_imp <- rep(NA_real_, length(genes)); names(mean_imp) <- genes
  if (length(imp_list) > 0) {
    all_imp <- do.call(c, imp_list)
    by_gene <- tapply(as.numeric(all_imp), names(all_imp), mean)
    common <- intersect(genes, names(by_gene))
    mean_imp[common] <- by_gene[common]
  }

  out <- data.frame(
    Gene = genes,
    Mean_Importance = round(as.numeric(mean_imp), 5),
    Times_In_Selected_Group = freq,
    N_Groups_Containing = ngrp,
    Mean_Group_Accuracy = round(macc, 4),
    stringsAsFactors = FALSE
  )
  # Sort by importance (NA last), then by selection frequency.
  out[order(-replace(out$Mean_Importance, is.na(out$Mean_Importance), -Inf),
            -out$Times_In_Selected_Group), , drop = FALSE]
}

#' Run summary / parameters as a key-value table.
#' @keywords internal
gsm_tbl_run_summary <- function(summary, dataset = NULL) {
  s <- summary %|null|% list()
  if (length(s) == 0 && !is.null(dataset)) {
    s <- list(dataset = sprintf("%d samples x %d genes",
                                dataset$params$n_samples, dataset$params$n_genes))
  }
  flat <- vapply(s, function(v) {
    if (is.null(v) || length(v) == 0) "" else paste(as.character(v), collapse = ", ")
  }, character(1))
  out <- data.frame(Parameter = names(flat), Value = unname(flat),
                    stringsAsFactors = FALSE)
  ver <- tryCatch(as.character(utils::packageVersion("gsmR")), error = function(e) NA)
  rbind(out,
        data.frame(Parameter = c("gsmR_version", "exported_at"),
                   Value = c(ver, format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
                   stringsAsFactors = FALSE))
}

# ── Single-fold table builders ───────────────────────────────────────────────

.single_tables <- function(result, dataset, summary) {
  sr <- result$score_results
  fm <- result$final_model
  accs <- vapply(sr, function(x) as.numeric(x$accuracy), numeric(1))
  g <- names(sr)
  group_scores <- data.frame(
    Rank = rank(-accs, ties.method = "min"),
    Group = g, Accuracy = round(accs, 6),
    N_genes = vapply(g, function(x) .n_genes(dataset, x), integer(1)),
    stringsAsFactors = FALSE
  )
  group_scores <- group_scores[order(group_scores$Rank), , drop = FALSE]

  m <- fm$metrics %|null|% list()
  final_metrics <- data.frame(
    Metric = names(m), Value = round(as.numeric(unlist(m)), 6),
    stringsAsFactors = FALSE
  )

  feats <- fm$features %|null|% character(0)
  imp <- fm$importance
  feat_imp <- if (!is.null(imp)) as.numeric(imp[feats]) else rep(NA_real_, length(feats))
  ngrp <- vapply(feats, function(x)
    sum(vapply(dataset$groups, function(gg) x %in% gg, logical(1))), integer(1))
  feature_list <- data.frame(
    Gene = feats, Mean_Importance = round(feat_imp, 5),
    N_Groups_Containing = ngrp, stringsAsFactors = FALSE
  )

  list(
    RunSummary   = gsm_tbl_run_summary(summary, dataset),
    FinalMetrics = final_metrics,
    GroupScores  = group_scores,
    FeatureList  = feature_list
  )
}

# ── Public API ───────────────────────────────────────────────────────────────

#' Build all result tables for a scoring run
#'
#' @param result An MCCV result (list with `folds`, from [run_mccv()]) or a
#'   single-fold result `list(score_results=, final_model=)`.
#' @param dataset The `gsm_dataset` that was scored.
#' @param summary Optional named list of run metadata (classifier, params, ...).
#' @param top_n Top-N cutoff used for the group ranking table. Default 10.
#' @return A named list of data.frames (each becomes a sheet / CSV).
#' @export
gsm_export_tables <- function(result, dataset, summary = NULL, top_n = 10) {
  stopifnot(inherits(dataset, "gsm_dataset"))
  if (is.null(result$folds)) {
    return(.single_tables(result, dataset, summary))
  }
  folds <- result$folds
  tbls <- list(
    RunSummary           = gsm_tbl_run_summary(summary, dataset),
    MetricSummary        = gsm_tbl_metric_summary(result),
    IterationMetrics     = gsm_tbl_iteration_metrics(folds),
    IterationGroupScores = gsm_tbl_iteration_group_scores(folds, dataset),
    GroupRanking         = gsm_tbl_group_ranking(folds, dataset, top_n = top_n),
    GroupFeatureAssoc    = gsm_tbl_group_feature_assoc(folds, dataset),
    PreScoreSignificance = gsm_tbl_pre_score_significance(folds),
    RankedFeatures       = gsm_tbl_ranked_features(folds, dataset)
  )
  # Drop tables that are empty/NULL (e.g. pre-scores absent for method "none").
  Filter(function(x) !is.null(x) && nrow(x) > 0, tbls)
}

#' Write result tables to a multi-sheet Excel workbook
#' @param tables Named list of data.frames (from [gsm_export_tables()]).
#' @param file Output `.xlsx` path.
#' @return `file`, invisibly.
#' @export
gsm_write_xlsx <- function(tables, file) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package 'openxlsx' is required for Excel export. ",
         "Install it with: install.packages('openxlsx')")
  }
  wb <- openxlsx::createWorkbook()
  hdr <- openxlsx::createStyle(textDecoration = "bold",
                               fgFill = "#4a90d9", fontColour = "white",
                               border = "bottom")
  for (nm in names(tables)) {
    sheet <- substr(nm, 1, 31)
    openxlsx::addWorksheet(wb, sheet)
    openxlsx::writeData(wb, sheet, tables[[nm]], headerStyle = hdr)
    openxlsx::freezePane(wb, sheet, firstActiveRow = 2)
    openxlsx::setColWidths(wb, sheet, cols = seq_along(tables[[nm]]),
                           widths = "auto")
  }
  openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
  invisible(file)
}

# Map table keys to original-G-S-M-style CSV file names.
.csv_name_map <- c(
  RunSummary           = "RunSummary",
  MetricSummary        = "res1_MetricSummary",
  IterationMetrics     = "res2_IterationMetrics",
  IterationGroupScores = "res2_IterationGroupScores",
  GroupRanking         = "RankdGroups_by_Freq",
  GroupFeatureAssoc    = "Groups_Features_Associations_over_all_iterations",
  PreScoreSignificance = "PreScoreSignificance",
  RankedFeatures       = "rankedFeatures",
  FinalMetrics         = "FinalMetrics",
  GroupScores          = "GroupScores",
  FeatureList          = "FeatureList"
)

#' Write result tables to a zip bundle of CSV files
#' @param tables Named list of data.frames (from [gsm_export_tables()]).
#' @param file Output `.zip` path.
#' @param prefix Optional file-name prefix (e.g. dataset + classifier).
#' @return `file`, invisibly.
#' @export
gsm_write_csv_zip <- function(tables, file, prefix = "GSM") {
  dir <- file.path(tempdir(), paste0("gsm_csv_", as.integer(Sys.time())))
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  paths <- character(0)
  for (nm in names(tables)) {
    label <- .csv_name_map[[nm]] %|null|% nm
    p <- file.path(dir, sprintf("%s_%s.csv", prefix, label))
    utils::write.csv(tables[[nm]], p, row.names = FALSE)
    paths <- c(paths, p)
  }
  # Zip with relative paths so the archive has flat CSVs.
  old <- setwd(dir); on.exit(setwd(old), add = TRUE)
  utils::zip(zipfile = file, files = basename(paths), flags = "-q")
  invisible(file)
}

#' Suggested output file name, e.g. "GDS3257_GSMwithRandomForest_MCCV_2026-06-07.xlsx"
#' @param summary Run metadata list (uses `dataset_name`, `classifier`, `mode`).
#' @param dataset The `gsm_dataset` (fallback for the dataset token).
#' @param ext File extension without dot ("xlsx" or "zip").
#' @return A character file name.
#' @export
gsm_export_filename <- function(summary, dataset = NULL, ext = "xlsx") {
  s <- summary %|null|% list()
  clean <- function(x) gsub("[^A-Za-z0-9]+", "", as.character(x))
  dn <- clean(s$dataset_name %|null|% "dataset")
  # Accept either an internal code ("rf") or a friendly label ("Random Forest").
  clf <- clean(switch(as.character(s$classifier %|null|% "model"),
                      "Random Forest" = "RandomForest",
                      "Support Vector Machine" = "SVM",
                      "K-Nearest Neighbors" = "KNN",
                      .clf_label(s$classifier %|null|% "model")))
  mode <- if (isTRUE(grepl("MCCV|mccv", s$mode %|null|% ""))) "MCCV" else "Single"
  sprintf("%s_GSMwith%s_%s_%s.%s", dn, clf, mode, format(Sys.Date()), ext)
}

#' Build and write result tables in one call (console convenience)
#' @param result MCCV or single-fold result.
#' @param dataset The `gsm_dataset`.
#' @param file Output path (extension decides default format if `format` missing).
#' @param summary Optional run metadata.
#' @param format "xlsx" or "csv".
#' @param top_n Top-N cutoff for the group ranking table.
#' @return `file`, invisibly.
#' @export
gsm_export <- function(result, dataset, file, summary = NULL,
                       format = c("xlsx", "csv"), top_n = 10) {
  format <- match.arg(format)
  tbls <- gsm_export_tables(result, dataset, summary, top_n = top_n)
  if (format == "xlsx") gsm_write_xlsx(tbls, file)
  else gsm_write_csv_zip(tbls, file)
}
