#' @title CLI Pipeline for gsmR
#' @description Run the full G-S-M analysis from the terminal using a YAML
#'   parameter file. Results are saved as `.rds` session files that the Shiny
#'   app can open directly.
#' @name pipeline
NULL

#' Run the full GSM-R pipeline
#'
#' Orchestrates the complete G-S-M analysis: load data, build dataset,
#' (optionally) collapse replicates, pre-score groups, run MCCV or single-fold
#' scoring, and (optionally) run module discovery. Results are saved as an
#' `.rds` session file compatible with the Shiny app's "Load Session" feature.
#'
#' @param params Path to a YAML parameter file (exported from the Shiny app's
#'   File > Export Settings), or a named list of parameters.
#' @param outdir Directory to write results. Created if it does not exist.
#'   Default: current directory.
#' @param discover Logical. Whether to run the discovery step. Default `TRUE`.
#' @param ... Override any parameter from the YAML file. These take precedence.
#'   Use the same names as the YAML keys (e.g.,
#'   `expression_file = "new.csv"`, `score_method = "svm"`).
#' @return Path to the saved `.rds` session file (invisibly).
#' @export
#'
#' @examples
#' \dontrun{
#'   # Run from an exported YAML file
#'   run_pipeline("params.yaml")
#'
#'   # Override the expression file
#'   run_pipeline("params.yaml", expression_file = "other_data.csv")
#'
#'   # Run without discovery, custom output directory
#'   run_pipeline("params.yaml", outdir = "results/", discover = FALSE)
#' }
run_pipeline <- function(params = NULL, outdir = ".", discover = TRUE, ...) {

  # ── 1. Load & merge parameters ──────────────────────────────────────────
  if (is.character(params) && length(params) == 1) {
    if (!file.exists(params)) stop("Parameter file not found: ", params)
    if (!requireNamespace("yaml", quietly = TRUE))
      stop("Package 'yaml' is required. Install with: install.packages('yaml')")
    cfg <- yaml::read_yaml(params)
    message("\u2139\ufe0f  Loaded parameters from: ", params)
  } else if (is.list(params)) {
    cfg <- params
  } else if (is.null(params)) {
    cfg <- list()
  } else {
    stop("`params` must be a file path (character), a named list, or NULL.")
  }

  # CLI overrides take precedence
  overrides <- list(...)
  for (nm in names(overrides)) {
    cfg[[nm]] <- overrides[[nm]]
  }

  # ── 2. Resolve parameters with defaults ─────────────────────────────────
  p <- gsm_resolve_params(cfg)

  # Validate required files
  if (is.null(p$expression_file) || !nzchar(p$expression_file))
    stop("'expression_file' is required. Specify it in the YAML or via run_pipeline(..., expression_file = '...').")
  if (is.null(p$groups_file) || !nzchar(p$groups_file))
    stop("'groups_file' is required. Specify it in the YAML or via run_pipeline(..., groups_file = '...').")
  if (!file.exists(p$expression_file))
    stop("Expression file not found: ", p$expression_file)
  if (!file.exists(p$groups_file))
    stop("Groups file not found: ", p$groups_file)

  # Resolve model_method
  actual_model_method <- if (identical(p$model_method, "same")) p$score_method else p$model_method

  # Build score_args based on method
  score_args <- .build_method_args(p$score_method, p)
  model_args <- .build_method_args(actual_model_method, p)

  # ── 3. Read data ────────────────────────────────────────────────────────
  message("\n\U0001F4C2 Reading expression file: ", p$expression_file)
  expr_data <- read_expression(p$expression_file)
  message("   ", ncol(expr_data$expr), " samples, ",
          nrow(expr_data$expr), " genes, ",
          nlevels(expr_data$labels), " classes (",
          paste(levels(expr_data$labels), collapse = ", "), ")")

  message("\U0001F4C2 Reading groups file: ", p$groups_file)
  fmt <- p$groups_format
  if (fmt == "auto") {
    ext <- tolower(tools::file_ext(p$groups_file))
    fmt <- if (ext == "gmt") "gmt" else "csv"
  }
  groups <- if (fmt == "gmt") read_gmt(p$groups_file) else read_groups(p$groups_file)
  message("   ", length(groups), " groups loaded")

  # ── 4. Build dataset ────────────────────────────────────────────────────
  message("\n\U0001F527 Building dataset...")
  dataset <- gsm_dataset(
    expr_data$expr, expr_data$labels, groups,
    min_group_size = p$min_group_size,
    max_group_size = p$max_group_size
  )
  message("   ", dataset$params$n_groups, " groups retained (",
          dataset$params$n_samples, " samples)")

  dataset_name <- tools::file_path_sans_ext(basename(p$expression_file))

  # ── 5. Run MCCV or single-fold scoring ──────────────────────────────────
  t0 <- proc.time()
  mccv_results <- NULL
  score_results <- NULL
  final_model <- NULL
  accumulative <- NULL
  score_summary <- NULL
  prescore_groups <- NULL
  prescore_scores <- NULL
  prescore_nkeep <- NULL
  train_idx <- NULL

  if (isTRUE(p$run_mccv_mode)) {
    message("\n\U0001F680 Running MCCV (",
            p$mccv_iterations, " iterations, ",
            p$mccv_cores, " cores, method: ", p$score_method, ")...")

    mccv_results <- run_mccv(
      dataset,
      n_iterations       = p$mccv_iterations,
      train_ratio        = p$mccv_train_ratio,
      pre_score_method   = p$prescore_method,
      pre_score_top_percent = p$top_percent,
      score_method       = p$score_method,
      model_method       = actual_model_method,
      score_args         = score_args,
      model_args         = model_args,
      top_n              = p$n_top_groups,
      shave              = p$shave_genes,
      fdr_alpha          = p$fdr_alpha,
      normalize          = p$normalize_z,
      cores              = p$mccv_cores,
      seed               = p$score_seed
    )

    # Build summary from MCCV results
    ok <- Filter(function(r) r$status == "success", mccv_results$folds)
    n_ok <- length(ok)
    message("   ", n_ok, "/", length(mccv_results$folds), " folds succeeded")
    if (n_ok > 0) {
      accs <- vapply(ok, function(r) as.numeric(r$metrics$Accuracy), numeric(1))
      message("   Mean accuracy: ", round(mean(accs), 4),
              " \u00b1 ", round(stats::sd(accs), 4))
    }
    score_summary <- list(
      mode = "mccv",
      n_folds = p$mccv_iterations,
      n_successful = n_ok,
      score_method = p$score_method,
      model_method = actual_model_method,
      dataset = dataset_name,
      seed = p$score_seed,
      prescore_method = p$prescore_method,
      prescore_top_percent = p$top_percent
    )
  } else {
    # Single-fold scoring
    message("\n\U0001F680 Running single-fold scoring (method: ", p$score_method, ")...")

    set.seed(p$prescore_seed)
    n <- ncol(dataset$expr)
    train_idx <- sort(sample(n, floor(n * p$train_ratio)))

    # Pre-score
    message("   Pre-scoring groups (", p$prescore_method, ", top ", p$top_percent, "%)...")
    ps_result <- pre_score_groups(dataset, train_idx,
                                   method = p$prescore_method,
                                   top_percent = p$top_percent)
    prescore_groups <- ps_result
    prescore_scores <- attr(ps_result, "scores")
    prescore_nkeep <- attr(ps_result, "n_keep")

    # Score each group
    groups_to_score <- if (identical(p$score_scope, "all")) {
      dataset$groups
    } else {
      head(ps_result, p$n_top_groups)
    }

    message("   Scoring ", length(groups_to_score), " groups...")
    score_results <- list()
    for (i in seq_along(groups_to_score)) {
      gname <- names(groups_to_score)[i]
      genes <- groups_to_score[[i]]
      res <- tryCatch(
        score_group(dataset$expr, dataset$labels, genes, train_idx,
                    method = p$score_method),
        error = function(e) list(accuracy = NA_real_, predictions = NULL,
                                  actual = NULL, error = conditionMessage(e))
      )
      res$group_name <- gname
      res$genes <- genes
      score_results[[gname]] <- res
      if (i %% 10 == 0 || i == length(groups_to_score))
        message("   ... ", i, "/", length(groups_to_score))
    }

    accs <- vapply(score_results, function(x) x$accuracy %||% NA_real_, numeric(1))
    message("   Mean accuracy: ", round(mean(accs, na.rm = TRUE), 4))

    score_summary <- list(
      mode = "single",
      score_method = p$score_method,
      n_groups_scored = length(score_results),
      dataset = dataset_name,
      seed = p$score_seed,
      prescore_method = p$prescore_method,
      prescore_top_percent = p$top_percent
    )
  }
  time_score <- (proc.time() - t0)[["elapsed"]]
  message("   Scoring completed in ", .fmt_time(time_score))

  # ── 6. Discovery ────────────────────────────────────────────────────────
  detect <- NULL
  if (isTRUE(discover) && !is.null(mccv_results)) {
    message("\n\U0001F30B Running module discovery (top ", p$detect_top_genes,
            " genes, mode: ", p$detect_mode, ")...")
    t1 <- proc.time()
    detect <- tryCatch(
      gsm_detect(
        mccv_results, dataset,
        mode       = p$detect_mode,
        top_genes  = p$detect_top_genes,
        min_freq   = p$detect_min_freq,
        resolution = p$detect_resolution,
        n_boot     = p$detect_nboot,
        seed       = 1
      ),
      error = function(e) {
        warning("Discovery failed: ", conditionMessage(e))
        NULL
      }
    )
    t1_elapsed <- (proc.time() - t1)[["elapsed"]]
    if (!is.null(detect)) {
      message("   Discovered ", nrow(detect$modules), " modules from ",
              detect$params$n_robust, " robust genes in ", .fmt_time(t1_elapsed))
    }
  } else if (isTRUE(discover) && is.null(mccv_results)) {
    message("\n\u26a0\ufe0f  Skipping discovery: requires MCCV results (run_mccv_mode: true).")
  }

  # ── 7. Save results ────────────────────────────────────────────────────
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")

  # Save session .rds (same format as Shiny)
  session_data <- gsm_build_session(
    expr_data = expr_data, groups = groups, dataset = dataset,
    dataset_name = dataset_name, p = p,
    mccv_results = mccv_results, score_results = score_results,
    final_model = final_model, accumulative = accumulative,
    prescore_groups = prescore_groups, prescore_scores = prescore_scores,
    prescore_nkeep = prescore_nkeep, train_idx = train_idx,
    score_summary = score_summary, time_score = time_score,
    detect = detect, source = "cli_pipeline"
  )

  rds_path <- file.path(outdir, sprintf("gsmR_results_%s_%s.rds", dataset_name, ts))
  saveRDS(session_data, rds_path)
  message("\n\U0001F4BE Session saved: ", rds_path)

  # Save timestamped parameter YAML
  yaml_path <- file.path(outdir, sprintf("gsmR_params_%s_%s.yaml", dataset_name, ts))
  # Write the actual parameters used (including overrides)
  p_export <- p
  # Convert logicals to yaml-friendly strings
  tryCatch({
    yaml_content <- yaml::as.yaml(p_export)
    # Add header comment
    header <- paste0(
      "# gsmR Pipeline Parameters\n",
      "# Generated: ", Sys.time(), "\n",
      "# Expression: ", p$expression_file, "\n",
      "# Groups: ", p$groups_file, "\n",
      "# Score method: ", p$score_method, "\n",
      "# MCCV mode: ", p$run_mccv_mode, "\n",
      "#\n"
    )
    writeLines(paste0(header, yaml_content), yaml_path)
    message("\U0001F4DD Parameters saved: ", yaml_path)
  }, error = function(e) {
    warning("Could not save parameters YAML: ", conditionMessage(e))
  })

  message("\n\u2705 Pipeline complete! Open the .rds file in the Shiny app to explore results.")
  invisible(rds_path)
}


# ── Internal helpers ──────────────────────────────────────────────────────────

#' Build classifier-specific argument list from pipeline parameters
#' @keywords internal
.build_method_args <- function(method, p) {
  switch(method,
    rf  = list(ntree = p$rf_ntree),
    svm = {
      gamma <- if (identical(p$svm_gamma, "auto") || is.null(p$svm_gamma)) NULL
               else as.numeric(p$svm_gamma)
      c(list(kernel = p$svm_kernel, cost = p$svm_cost),
        if (!is.null(gamma)) list(gamma = gamma) else list())
    },
    knn = list(k = p$knn_k),
    xgb = list(),
    mlp = list(hidden = .parse_hidden(p$mlp_hidden), dropout = p$nn_dropout,
               epochs = p$nn_epochs, lr = p$nn_lr),
    # `graph` is not set here: MCCV injects BRIDGE's per-fold differential
    # network at run time (see mccv_run_fold()). Absent that, the GNN builds a
    # train-only correlation graph itself.
    gnn = list(hidden = as.integer(p$gnn_hidden), dropout = p$nn_dropout,
               epochs = p$nn_epochs, lr = p$nn_lr,
               cor_threshold = p$gnn_cor_threshold),
    list()
  )
}

# "64,32" (YAML/GUI) or c(64, 32) -> integer vector of hidden-layer widths.
#' @keywords internal
.parse_hidden <- function(x) {
  if (is.null(x)) return(c(64L, 32L))
  if (is.numeric(x)) return(as.integer(x))
  v <- suppressWarnings(as.integer(trimws(strsplit(as.character(x), ",")[[1]])))
  v <- v[is.finite(v) & v > 0]
  if (length(v) == 0) c(64L, 32L) else v
}

#' Format elapsed seconds as human-readable string
#' @keywords internal
.fmt_time <- function(secs) {
  if (secs < 60) return(sprintf("%.1fs", secs))
  if (secs < 3600) return(sprintf("%.0fm %.0fs", secs %/% 60, secs %% 60))
  sprintf("%.0fh %.0fm", secs %/% 3600, (secs %% 3600) %/% 60)
}

#' Resolve gsmR parameters against their defaults
#'
#' Merges a (possibly partial) configuration list with the canonical defaults
#' used by [run_pipeline()]. Shared with the multi-dataset batch runner so both
#' resolve parameters identically. Does **not** validate input files — callers
#' that require them (e.g. [run_pipeline()]) check separately, since the batch
#' runner sets `expression_file` per dataset.
#'
#' @param cfg A named list of parameters (any subset; missing keys take defaults).
#' @return A fully-populated named list of resolved parameters.
#' @keywords internal
gsm_resolve_params <- function(cfg = list()) {
  list(
    expression_file    = cfg$expression_file    %||% NULL,
    groups_file        = cfg$groups_file        %||% NULL,
    groups_format      = cfg$groups_format      %||% "auto",
    min_group_size     = as.integer(cfg$min_group_size %||% 5),
    max_group_size     = as.integer(cfg$max_group_size %||% 500),
    collapse_method    = cfg$collapse_method    %||% "none",
    meta_conflict      = cfg$meta_conflict      %||% "first",
    keep_original      = as.logical(cfg$keep_original %||% FALSE),
    prescore_method    = cfg$prescore_method    %||% "anova",
    top_percent        = as.numeric(cfg$top_percent %||% 20),
    train_ratio        = as.numeric(cfg$train_ratio %||% 0.7),
    prescore_seed      = as.integer(cfg$prescore_seed %||% 42),
    score_method       = cfg$score_method       %||% "rf",
    rf_ntree           = as.integer(cfg$rf_ntree %||% 500),
    svm_kernel         = cfg$svm_kernel         %||% "radial",
    svm_cost           = as.numeric(cfg$svm_cost %||% 1),
    svm_gamma          = cfg$svm_gamma          %||% "auto",
    knn_k              = as.integer(cfg$knn_k %||% 5),
    mlp_hidden         = cfg$mlp_hidden         %||% "64,32",
    gnn_hidden         = as.integer(cfg$gnn_hidden %||% 32),
    gnn_cor_threshold  = as.numeric(cfg$gnn_cor_threshold %||% 0.3),
    nn_dropout         = as.numeric(cfg$nn_dropout %||% 0.3),
    nn_epochs          = as.integer(cfg$nn_epochs %||% 100),
    nn_lr              = as.numeric(cfg$nn_lr %||% 0.01),
    shave_genes        = as.logical(cfg$shave_genes %||% TRUE),
    fdr_alpha          = as.numeric(cfg$fdr_alpha %||% 0.05),
    normalize_z        = as.logical(cfg$normalize_z %||% FALSE),
    score_scope        = cfg$score_scope        %||% "selected",
    n_top_groups       = as.integer(cfg$n_top_groups %||% 10),
    score_seed         = as.integer(cfg$score_seed %||% 44),
    run_mccv_mode      = as.logical(cfg$run_mccv_mode %||% TRUE),
    mccv_iterations    = as.integer(cfg$mccv_iterations %||% 100),
    model_method       = cfg$model_method       %||% "same",
    mccv_train_ratio   = as.numeric(cfg$mccv_train_ratio %||% 0.9),
    mccv_cores         = as.integer(cfg$mccv_cores %||% 1),
    detect_mode        = cfg$detect_mode        %||% "abs",
    detect_top_genes   = as.integer(cfg$detect_top_genes %||% 200),
    detect_min_freq    = as.numeric(cfg$detect_min_freq %||% 0.5),
    detect_resolution  = as.numeric(cfg$detect_resolution %||% 1),
    detect_nboot       = as.integer(cfg$detect_nboot %||% 50)
  )
}

#' Assemble a gsmR session object (the `.rds` payload)
#'
#' Builds the session list saved by [run_pipeline()] and loaded by the Shiny app.
#' Factored out so the multi-dataset batch runner writes per-dataset `.rds` files
#' with an identical shape — making them openable in the existing results tabs.
#'
#' @param expr_data,groups,dataset,dataset_name Loaded inputs and built dataset.
#' @param p Resolved parameter list (see [gsm_resolve_params()]).
#' @param mccv_results,score_results,final_model,accumulative Scoring outputs.
#' @param prescore_groups,prescore_scores,prescore_nkeep,train_idx Pre-score outputs.
#' @param score_summary,time_score,detect Run summary, timing, discovery output.
#' @param source A short tag identifying the producer (e.g. "cli_pipeline",
#'   "batch").
#' @return A named list (the session payload) ready for [saveRDS()].
#' @keywords internal
gsm_build_session <- function(expr_data, groups, dataset, dataset_name, p,
                              mccv_results = NULL, score_results = NULL,
                              final_model = NULL, accumulative = NULL,
                              prescore_groups = NULL, prescore_scores = NULL,
                              prescore_nkeep = NULL, train_idx = NULL,
                              score_summary = NULL, time_score = NULL,
                              detect = NULL, source = "cli_pipeline") {
  list(
    # Data
    expr_data        = expr_data,
    groups           = groups,
    dataset          = dataset,
    dataset_name     = dataset_name,
    loaded           = TRUE,
    # Collapse (not run in CLI/batch — placeholder)
    collapse_result  = NULL,
    collapsed        = FALSE,
    time_collapse    = NULL,
    # Pre-score
    prescore_groups  = prescore_groups,
    prescore_scores  = prescore_scores,
    prescore_nkeep   = prescore_nkeep,
    prescored        = !is.null(prescore_groups),
    train_idx        = train_idx,
    time_prescore    = NULL,
    # Score
    score_results    = score_results,
    final_model      = final_model,
    accumulative     = accumulative,
    mccv_results     = mccv_results,
    scored           = TRUE,
    score_summary    = score_summary,
    time_score       = time_score,
    # Discovery
    detect           = detect,
    permtest         = NULL,
    mccv_group_summary = NULL,
    # Settings (for Shiny UI state restoration)
    settings = list(
      score_method      = p$score_method,
      run_mccv_mode     = p$run_mccv_mode,
      mccv_iterations   = p$mccv_iterations,
      model_method      = p$model_method,
      mccv_train_ratio  = p$mccv_train_ratio,
      score_scope       = p$score_scope,
      n_top_groups      = p$n_top_groups,
      score_seed        = p$score_seed,
      top_percent       = p$top_percent,
      normalize         = p$normalize_z,
      detect_mode       = p$detect_mode,
      detect_top_genes  = p$detect_top_genes,
      detect_min_freq   = p$detect_min_freq,
      detect_resolution = p$detect_resolution,
      detect_nboot      = p$detect_nboot
    ),
    # Metadata
    saved_at     = Sys.time(),
    r_version    = R.version.string,
    gsmR_version = tryCatch(as.character(utils::packageVersion("gsmR")),
                            error = function(e) "unknown"),
    source       = source
  )
}
