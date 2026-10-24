#' @title Multi-dataset batch runner
#' @description Run the full G-S-M pipeline on many expression datasets that
#'   share a single grouping file, scheduling all `(dataset, fold)` work units
#'   over one shared pool of cores. Each dataset produces its own GUI-openable
#'   `.rds` session plus tables and a log; a combined cross-dataset report
#'   aggregates performance and consensus group/gene rankings.
#' @name batch
NULL

# ──────────────────────────────────────────────────────────────────────────────
# Scheduler design
#
# We flatten every (dataset, fold) pair into ONE queue scheduled over a single
# pool of `cores` workers. A fold (mccv_run_fold) is the atomic, single-threaded
# unit of work, and folds from different datasets are independent. This:
#   * keeps exactly `cores` folds running at once (no oversubscription),
#   * auto load-balances (a small dataset's folds free cores for a big one),
#   * bounds memory: datasets are loaded just-in-time and released once their
#     folds finish, with at most `max_resident` datasets in memory at a time,
#   * preserves per-fold seeding (seed + fold) so results match single runs.
#
# The state lives in an environment so the same stepper drives both the CLI
# (a while-loop) and the Shiny app (a non-blocking observe/invalidateLater loop).
# ──────────────────────────────────────────────────────────────────────────────

#' Create batch scheduler state
#'
#' @param expression_files Character vector of expression file paths (optionally
#'   named; names become dataset labels).
#' @param groups_file Path to the single shared grouping file.
#' @param params Named list of pipeline parameters (see [gsm_resolve_params()]).
#' @param outdir Output directory root.
#' @param cores Total worker cores for the shared fold pool. Default
#'   `detectCores() - 1`.
#' @param max_resident Maximum datasets held in memory at once.
#' @param discover Run module discovery per dataset.
#' @param groups_format `"auto"`, `"gmt"`, or `"csv"`.
#' @param progress Optional callback `progress(snapshot)` invoked each step.
#' @return An environment holding the scheduler state.
#' @keywords internal
gsm_batch_new_state <- function(expression_files, groups_file, params = list(),
                                outdir = "gsmR_batch", cores = NULL,
                                max_resident = 3, discover = TRUE,
                                groups_format = "auto", progress = NULL) {
  stopifnot(length(expression_files) >= 1)
  if (!file.exists(groups_file)) stop("Groups file not found: ", groups_file)
  missing <- expression_files[!file.exists(expression_files)]
  if (length(missing)) stop("Expression file(s) not found: ",
                            paste(missing, collapse = ", "))

  if (is.null(cores)) cores <- max(1L, parallel::detectCores() - 1L)

  nms <- names(expression_files)
  if (is.null(nms) || any(!nzchar(nms)))
    nms <- tools::file_path_sans_ext(basename(expression_files))
  nms <- make.unique(nms, sep = "_")
  names(expression_files) <- nms

  p <- gsm_resolve_params(params)
  fmt <- groups_format
  if (fmt == "auto")
    fmt <- if (tolower(tools::file_ext(groups_file)) == "gmt") "gmt" else "csv"
  groups <- if (fmt == "gmt") read_gmt(groups_file) else read_groups(groups_file)

  dir.create(file.path(outdir, "per_dataset"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(outdir, "combined"), showWarnings = FALSE, recursive = TRUE)

  st <- new.env(parent = emptyenv())
  st$expression_files <- expression_files
  st$names      <- nms
  st$groups     <- groups
  st$groups_file <- groups_file
  st$p          <- p
  st$cores      <- as.integer(cores)
  st$max_resident <- as.integer(max(1L, max_resident))
  st$discover   <- isTRUE(discover)
  st$outdir     <- outdir
  st$n_iter     <- p$mccv_iterations
  st$actual_model_method <- if (identical(p$model_method, "same")) p$score_method else p$model_method
  st$score_args <- .build_method_args(p$score_method, p)
  st$model_args <- .build_method_args(st$actual_model_method, p)
  st$progress   <- progress
  st$can_fork   <- .Platform$OS.type != "windows"

  st$status     <- stats::setNames(rep("pending", length(nms)), nms)
  st$loaded     <- list()   # name -> gsm_dataset
  st$exprdata   <- list()   # name -> list(expr, labels)
  st$folds      <- list()   # name -> list(length n_iter)
  st$done_count <- stats::setNames(integer(length(nms)), nms)
  st$ds_t0      <- stats::setNames(numeric(length(nms)), nms)
  st$window     <- character(0)
  st$next_ds    <- 1L
  st$queue      <- list()   # list of list(name=, fold=)
  st$running    <- list()   # key -> mcparallel job
  st$running_meta <- list() # key -> list(name=, fold=)
  st$results    <- list()   # name -> finalized per-dataset result
  st$manifest   <- list()
  st$errors     <- character(0)  # name -> failure message (for UI)
  st$log        <- character(0)
  st$logfile    <- file.path(outdir, "batch.log")
  st$active     <- TRUE
  st$all_done   <- FALSE
  st$outputs_written <- FALSE
  st$t0         <- proc.time()[["elapsed"]]
  st$started_at <- Sys.time()
  st$ended_at   <- NULL
  st
}

#' Build one manifest row with a fixed column set (so rows always rbind)
#' @keywords internal
.gsm_manifest_row <- function(dataset, status, n_ok = 0L, samples = NA_integer_,
                              genes = NA_integer_, groups = NA_integer_,
                              best_k = NA_real_, acc = NA_real_, fm = NA_real_,
                              auc = NA_real_, elapsed = NA_real_,
                              rds = NA_character_, error = NA_character_) {
  data.frame(
    Dataset = dataset, Status = status, N_folds_ok = as.integer(n_ok),
    Samples = samples, Genes = genes, Groups = groups,
    Mean_Best_k = round(best_k, 2), Mean_Accuracy = acc, Mean_F_measure = fm,
    Mean_AUC = auc, Elapsed_sec = round(elapsed, 1), RDS = rds, Error = error,
    stringsAsFactors = FALSE)
}

#' Append a line to the batch log (ring buffer + on-disk batch.log)
#' @keywords internal
gsm_batch_log <- function(st, msg) {
  line <- sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), msg)
  st$log <- utils::tail(c(st$log, line), 400L)
  tryCatch(cat(line, "\n", sep = "", file = st$logfile, append = TRUE),
           error = function(e) NULL)
  if (!is.null(st$progress)) tryCatch(st$progress(gsm_batch_snapshot(st)),
                                      error = function(e) NULL)
  invisible(NULL)
}

#' Small progress snapshot for UI / callbacks
#' @keywords internal
gsm_batch_snapshot <- function(st) {
  list(
    status     = st$status,
    done_count = st$done_count,
    n_iter     = st$n_iter,
    running    = length(st$running),
    queued     = length(st$queue),
    window     = st$window,
    active     = st$active,
    all_done   = st$all_done,
    errors     = st$errors,
    log_tail   = utils::tail(st$log, 25L)
  )
}

#' Admit the next dataset into the resident window
#' @keywords internal
.gsm_batch_admit <- function(st) {
  while (length(st$window) < st$max_resident && st$next_ds <= length(st$names)) {
    nm   <- st$names[st$next_ds]
    path <- st$expression_files[[nm]]
    st$next_ds <- st$next_ds + 1L
    ok <- tryCatch({
      st$status[nm] <- "loading"
      ed <- read_expression(path)
      ds <- gsm_dataset(ed$expr, ed$labels, st$groups,
                        min_group_size = st$p$min_group_size,
                        max_group_size = st$p$max_group_size)

      # Minimum group sayısı kontrolü: pre-score top% sonrası en az 3 grup kalmalı
      min_groups_needed <- max(3L, ceiling(ds$params$n_groups * st$p$top_percent / 100))
      if (ds$params$n_groups < min_groups_needed || ds$params$n_groups < 3L) {
        warn_msg <- sprintf(
          "Skipped '%s': only %d groups after size filtering (need >= %d for top_percent=%g%%)",
          nm, ds$params$n_groups, min_groups_needed, st$p$top_percent
        )
        st$status[nm] <- "skipped"
        st$errors[nm] <- warn_msg
        st$manifest[[length(st$manifest) + 1L]] <-
          .gsm_manifest_row(nm, "skipped", error = warn_msg)
        gsm_batch_log(st, sprintf("SKIPPED '%s': %s", nm, warn_msg))
        next
      }

      st$exprdata[[nm]] <- ed
      st$loaded[[nm]]   <- ds
      st$folds[[nm]]    <- vector("list", st$n_iter)
      st$ds_t0[nm]      <- proc.time()[["elapsed"]]
      st$window <- c(st$window, nm)
      st$status[nm] <- "running"
      # enqueue this dataset's folds
      for (f in seq_len(st$n_iter))
        st$queue[[length(st$queue) + 1L]] <- list(name = nm, fold = f)
      gsm_batch_log(st, sprintf("Loaded '%s' (%d samples, %d groups) -> queued %d folds",
                                nm, ds$params$n_samples, ds$params$n_groups, st$n_iter))
      TRUE
    }, error = function(e) {
      st$status[nm] <- "error"
      st$errors[nm] <- conditionMessage(e)
      st$manifest[[length(st$manifest) + 1L]] <-
        .gsm_manifest_row(nm, "error", error = conditionMessage(e))
      gsm_batch_log(st, sprintf("ERROR loading '%s': %s", nm, conditionMessage(e)))
      FALSE
    })
    invisible(ok)
  }
}


#' Launch queued fold jobs until the core budget is full
#' @keywords internal
.gsm_batch_launch <- function(st) {
  p <- st$p
  while (length(st$running) < st$cores && length(st$queue) > 0) {
    task <- st$queue[[1]]; st$queue[[1]] <- NULL
    nm <- task$name; f <- task$fold
    ds <- st$loaded[[nm]]
    key <- paste0(nm, "#", f)
    run_one <- function() {
      tryCatch(
        mccv_run_fold(ds, f, train_ratio = p$mccv_train_ratio,
                      pre_score_method = p$prescore_method,
                      pre_score_top_percent = p$top_percent,
                      pre_score_q_threshold = p$prescore_q_threshold %||% 0.05,
                      score_method = p$score_method,
                      model_method = st$actual_model_method,
                      score_args = st$score_args, model_args = st$model_args,
                      top_n = p$n_top_groups, shave = p$shave_genes,
                      fdr_alpha = p$fdr_alpha, normalize = p$normalize_z,
                      seed = p$score_seed),
        error = function(e) list(fold = f, status = "error",
                                 error_message = conditionMessage(e)))
    }
    if (st$can_fork) {
      st$running[[key]] <- parallel::mcparallel(run_one(), name = key)
      st$running_meta[[key]] <- task
    } else {
      # Windows / no-fork: run synchronously (one fold this tick)
      fr <- run_one()
      .gsm_batch_store_fold(st, nm, f, fr)
      break
    }
  }
  invisible(NULL)
}

#' Store a finished fold result into its dataset accumulator
#' @keywords internal
.gsm_batch_store_fold <- function(st, nm, f, fr) {
  st$folds[[nm]][[f]] <- fr
  st$done_count[nm] <- st$done_count[nm] + 1L
  stat <- if (is.list(fr) && identical(fr$status, "success")) {
    acc <- tryCatch(round(as.numeric(fr$metrics$Accuracy), 3), error = function(e) NA)
    sprintf("ok (acc=%s)", acc)
  } else "FAILED"
  gsm_batch_log(st, sprintf("  %s fold %d/%d %s", nm, f, st$n_iter, stat))
}

#' Reap finished fold jobs (non-blocking)
#' @keywords internal
.gsm_batch_reap <- function(st) {
  if (!st$can_fork || length(st$running) == 0) return(invisible(NULL))
  for (key in names(st$running)) {
    got <- tryCatch(parallel::mccollect(st$running[[key]], wait = FALSE, timeout = 0),
                    error = function(e) NULL)
    if (!is.null(got)) {
      meta <- st$running_meta[[key]]
      fr <- got[[1]]
      st$running[[key]] <- NULL
      st$running_meta[[key]] <- NULL
      .gsm_batch_store_fold(st, meta$name, meta$fold, fr)
    }
  }
  invisible(NULL)
}

#' Finalize any dataset whose folds are all complete
#' @keywords internal
.gsm_batch_finalize_ready <- function(st) {
  for (nm in st$window) {
    if (st$done_count[nm] >= st$n_iter) {
      .gsm_batch_finalize_one(st, nm)
    }
  }
  invisible(NULL)
}

#' Finalize one dataset: aggregate, save .rds/.yaml/tables/log, then release memory
#' @keywords internal
.gsm_batch_finalize_one <- function(st, nm) {
  p <- st$p
  ds <- st$loaded[[nm]]
  ed <- st$exprdata[[nm]]
  folds <- st$folds[[nm]]
  elapsed <- proc.time()[["elapsed"]] - st$ds_t0[nm]
  dsdir <- file.path(st$outdir, "per_dataset", nm)
  dir.create(dsdir, showWarnings = FALSE, recursive = TRUE)
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")

  agg <- tryCatch(aggregate_mccv_folds(folds), error = function(e) NULL)

  if (is.null(agg)) {
    # Surface a representative fold error so the UI can explain the failure.
    ferr <- NULL
    for (fr in folds) if (is.list(fr) && !is.null(fr$error_message)) {
      ferr <- fr$error_message; break
    }
    msg <- if (!is.null(ferr)) ferr else "all folds failed"
    st$status[nm] <- "error"
    st$errors[nm] <- msg
    st$manifest[[length(st$manifest) + 1L]] <- .gsm_manifest_row(
      nm, "error", samples = ds$params$n_samples, genes = ds$params$n_genes,
      groups = ds$params$n_groups, elapsed = elapsed, error = msg)
    gsm_batch_log(st, sprintf("Dataset '%s' FAILED: %s", nm, msg))
  } else {
    n_ok <- sum(vapply(folds, function(r) is.list(r) && identical(r$status, "success"),
                       logical(1)))
    score_summary <- list(
      mode = "mccv", n_folds = st$n_iter, n_successful = n_ok,
      score_method = p$score_method, model_method = st$actual_model_method,
      dataset = nm, seed = p$score_seed, prescore_method = p$prescore_method,
      prescore_top_percent = p$top_percent)

    detect <- NULL
    if (st$discover) {
      detect <- tryCatch(
        gsm_detect(agg, ds, mode = p$detect_mode, top_genes = p$detect_top_genes,
                   min_freq = p$detect_min_freq, resolution = p$detect_resolution,
                   n_boot = p$detect_nboot, seed = 1),
        error = function(e) { gsm_batch_log(st, sprintf("  %s discovery failed: %s",
                                                        nm, conditionMessage(e))); NULL })
    }

    session <- gsm_build_session(
      expr_data = ed, groups = st$groups, dataset = ds, dataset_name = nm, p = p,
      mccv_results = agg, score_summary = score_summary, time_score = elapsed,
      detect = detect, source = "batch")

    rds_path <- file.path(dsdir, sprintf("gsmR_results_%s_%s.rds", nm, ts))
    saveRDS(session, rds_path)

    # params yaml
    tryCatch({
      pe <- p; pe$expression_file <- st$expression_files[[nm]]; pe$groups_file <- st$groups_file
      writeLines(paste0("# gsmR batch params for ", nm, "\n", yaml::as.yaml(pe)),
                 file.path(dsdir, sprintf("gsmR_params_%s_%s.yaml", nm, ts)))
    }, error = function(e) NULL)

    # tables (xlsx + csv.zip); attach detect for export
    res_for_export <- agg; res_for_export$detect <- detect
    tryCatch(gsm_export(res_for_export, ds, file.path(dsdir, paste0(nm, "_tables.xlsx")),
                        summary = score_summary, format = "xlsx"),
             error = function(e) gsm_batch_log(st, sprintf("  %s xlsx export skipped: %s",
                                                          nm, conditionMessage(e))))
    tryCatch(gsm_export(res_for_export, ds, file.path(dsdir, paste0(nm, "_tables.zip")),
                        summary = score_summary, format = "csv"),
             error = function(e) gsm_batch_log(st, sprintf("  %s csv export skipped: %s",
                                                          nm, conditionMessage(e))))

    # per-dataset log
    mm <- agg$metrics_summary
    tryCatch(writeLines(c(
      sprintf("gsmR batch result: %s", nm),
      sprintf("expression: %s", st$expression_files[[nm]]),
      sprintf("groups: %s", st$groups_file),
      sprintf("classifier: %s | mccv_iterations: %d | folds ok: %d/%d",
              p$score_method, st$n_iter, n_ok, st$n_iter),
      sprintf("elapsed: %s", .fmt_time(elapsed)), "", "Mean metrics across folds:",
      utils::capture.output(print(mm))),
      file.path(dsdir, paste0(nm, ".log"))), error = function(e) NULL)

    getm <- function(k) tryCatch(mm[k, "Mean"], error = function(e) NA_real_)
    best_k_mean <- mean(vapply(folds, function(fr)
      if (is.list(fr) && identical(fr$status, "success"))
        as.numeric(fr$best_k %||% NA_real_) else NA_real_, numeric(1)), na.rm = TRUE)
    st$results[[nm]] <- list(name = nm, mccv = agg, summary = score_summary,
                             rds = rds_path, metrics_summary = mm, n_ok = n_ok,
                             elapsed = elapsed, detect = detect,
                             n_samples = ds$params$n_samples,
                             n_genes = ds$params$n_genes,
                             n_groups = ds$params$n_groups,
                             best_k_mean = best_k_mean)
    st$manifest[[length(st$manifest) + 1L]] <- .gsm_manifest_row(
      nm, "success", n_ok = n_ok, samples = ds$params$n_samples,
      genes = ds$params$n_genes, groups = ds$params$n_groups,
      best_k = best_k_mean, acc = getm("Accuracy"), fm = getm("F_measure"),
      auc = getm("AUC"), elapsed = elapsed, rds = rds_path)
    st$status[nm] <- "done"
    gsm_batch_log(st, sprintf("Dataset '%s' DONE: %d/%d folds, mean acc=%s, in %s",
                              nm, n_ok, st$n_iter, round(getm("Accuracy"), 3),
                              .fmt_time(elapsed)))
  }

  # release memory: drop the big expression matrix + dataset; keep folds-summary in results
  st$loaded[[nm]]   <- NULL
  st$exprdata[[nm]] <- NULL
  st$folds[[nm]]    <- NULL
  st$window <- setdiff(st$window, nm)
  gc(verbose = FALSE)
  invisible(NULL)
}

#' Advance the batch scheduler by one tick
#'
#' Reaps finished folds, finalizes completed datasets (freeing their memory),
#' admits new datasets into the resident window, and launches queued folds up to
#' the core budget. Mutates `st` in place. Sets `st$active <- FALSE` when all
#' work is complete.
#'
#' @param st A scheduler-state environment from [gsm_batch_new_state()].
#' @return `st`, invisibly.
#' @keywords internal
gsm_batch_step <- function(st) {
  if (!isTRUE(st$active)) return(invisible(st))
  .gsm_batch_reap(st)
  .gsm_batch_finalize_ready(st)
  .gsm_batch_admit(st)
  .gsm_batch_launch(st)

  if (length(st$queue) == 0 && length(st$running) == 0 &&
      st$next_ds > length(st$names) && length(st$window) == 0) {
    st$active <- FALSE
    st$all_done <- TRUE
  }
  invisible(st)
}

#' Full run-info: every chosen parameter + environment + timing (key/value)
#' @keywords internal
gsm_batch_run_info <- function(st) {
  p <- st$p
  n_ok <- sum(st$status == "done"); n_err <- sum(st$status == "error")
  total <- proc.time()[["elapsed"]] - st$t0
  ended <- if (!is.null(st$ended_at)) st$ended_at else Sys.time()
  kv <- function(k, v) data.frame(Parameter = k, Value = as.character(v),
                                  stringsAsFactors = FALSE)
  do.call(rbind, list(
    kv("Run started",            format(st$started_at, "%Y-%m-%d %H:%M:%S")),
    kv("Run finished",           format(ended, "%Y-%m-%d %H:%M:%S")),
    kv("Total wall-clock",       .fmt_time(total)),
    kv("Total wall-clock (sec)", round(total, 1)),
    kv("Datasets (total)",       length(st$names)),
    kv("Datasets (succeeded)",   n_ok),
    kv("Datasets (failed)",      n_err),
    kv("Datasets",               paste(st$names, collapse = ", ")),
    kv("Grouping file",          basename(st$groups_file)),
    kv("Classifier (scoring)",   p$score_method),
    kv("Final modeling",         st$actual_model_method),
    kv("Pre-scoring method",     p$prescore_method),
    kv("Pre-score top %",        p$top_percent),
    kv("MCCV iterations",        p$mccv_iterations),
    kv("Train ratio",            p$mccv_train_ratio),
    kv("Top-N groups",           p$n_top_groups),
    kv("Gene shaving",           p$shave_genes),
    kv("FDR alpha (shaving)",    p$fdr_alpha),
    kv("Z-normalize",            p$normalize_z),
    kv("Module discovery",       st$discover),
    kv("Seed",                   p$score_seed),
    kv("Group size range",       paste0(p$min_group_size, " - ", p$max_group_size)),
    kv("CPU cores",              st$cores),
    kv("Max resident datasets",  st$max_resident),
    kv("Output directory",       st$outdir),
    kv("Platform",               .Platform$OS.type),
    kv("R version",              R.version.string),
    kv("gsmR version",           tryCatch(as.character(utils::packageVersion("gsmR")),
                                          error = function(e) "unknown"))
  ))
}

#' Write the combined cross-dataset report + manifest + run-info (idempotent)
#' @keywords internal
gsm_write_batch_outputs <- function(st, format = c("xlsx", "csv")) {
  if (isTRUE(st$outputs_written)) return(invisible(NULL))
  format <- match.arg(format)
  st$ended_at <- Sys.time()
  manifest <- if (length(st$manifest))
    do.call(rbind, st$manifest) else data.frame()
  utils::write.csv(manifest, file.path(st$outdir, "batch_manifest.csv"),
                   row.names = FALSE)

  run_info <- gsm_batch_run_info(st)
  utils::write.csv(run_info, file.path(st$outdir, "batch_run_info.csv"),
                   row.names = FALSE)

  ok_results <- Filter(function(r) !is.null(r$mccv), st$results)
  if (length(ok_results) > 0) {
    file <- file.path(st$outdir, "combined",
                      if (format == "xlsx") "combined_report.xlsx" else "combined_tables.zip")
    tryCatch(gsm_export_combined_report(ok_results, file, format = format,
                                        manifest = manifest, run_info = run_info),
             error = function(e) gsm_batch_log(st,
               sprintf("Combined report failed: %s", conditionMessage(e))))
  }
  gsm_batch_log(st, sprintf("Batch complete: %d datasets ok, %d failed, total %s",
                            length(ok_results),
                            length(st$names) - length(ok_results),
                            .fmt_time(proc.time()[["elapsed"]] - st$t0)))
  st$outputs_written <- TRUE
  invisible(NULL)
}

#' Run the G-S-M pipeline on many datasets in parallel
#'
#' Schedules every `(dataset, fold)` work unit over one shared pool of `cores`
#' workers (see *Details* in [gsm_batch_new_state()] for the rationale). Each
#' dataset produces its own GUI-openable `.rds` session plus xlsx/csv tables and
#' a log under `outdir/per_dataset/<name>/`; a combined cross-dataset report
#' (per-dataset performance, consensus groups, consensus biomarkers) is written
#' under `outdir/combined/`.
#'
#' @inheritParams gsm_batch_new_state
#' @param report_format `"xlsx"` or `"csv"` for the combined report.
#' @return A manifest data.frame (one row per dataset), invisibly.
#' @export
#' @examples
#' \dontrun{
#' run_batch(
#'   expression_files = c("GDS1.csv", "GDS2.csv"),
#'   groups_file = "kegg.csv",
#'   params = list(score_method = "rf", mccv_iterations = 10),
#'   cores = 8, max_resident = 2
#' )
#' }
run_batch <- function(expression_files, groups_file, params = list(),
                      outdir = "gsmR_batch", cores = NULL, max_resident = 3,
                      discover = TRUE, groups_format = "auto", progress = NULL,
                      report_format = c("xlsx", "csv")) {
  report_format <- match.arg(report_format)
  st <- gsm_batch_new_state(expression_files, groups_file, params, outdir, cores,
                            max_resident, discover, groups_format, progress)
  gsm_batch_log(st, sprintf("Starting batch: %d datasets, %d cores, max_resident=%d, %d folds each",
                            length(st$names), st$cores, st$max_resident, st$n_iter))
  while (isTRUE(st$active)) {
    gsm_batch_step(st)
    # Yield while forked folds compute so we don't busy-spin a core.
    if (st$can_fork && length(st$running) > 0) Sys.sleep(0.05)
  }
  gsm_write_batch_outputs(st, format = report_format)
  invisible(if (length(st$manifest)) do.call(rbind, st$manifest) else data.frame())
}
