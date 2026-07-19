#' @title Parallel Monte Carlo Cross-Validation (MCCV)
#' @description Implementations for parallel execution of MCCV across available
#'   CPU cores in gsmR.
#' @name mccv
NULL

#' Run a single MCCV fold
#'
#' Executes one Monte Carlo CV iteration: stratified split, pre-scoring,
#' per-group scoring, and final joint model. Exported so callers (e.g. the
#' Shiny app) can run folds incrementally and stream live per-fold results.
#'
#' @param dataset A `gsm_dataset` object.
#' @param fold Integer fold index (also offsets the seed: `seed + fold`).
#' @param train_ratio,pre_score_method,pre_score_top_percent,score_method,model_method,score_args,model_args,seed
#'   See [run_mccv()].
#' @return A list with `fold`, `status` ("success"/"error"), and on success
#'   `metrics`, `features_count`, `n_groups`, `group_accuracies`, `pre_scores`
#'   (sorted pre-scoring scores for all groups), and `importance` (per-gene RF
#'   importance, `NULL` for SVM/KNN).
#' @export
mccv_run_fold <- function(dataset, fold, train_ratio = 0.9,
                          pre_score_method = "anova", pre_score_top_percent = 1,
                          score_method = "rf", model_method = "rf",
                          score_args = list(), model_args = list(),
                          top_n = 10, shave = TRUE, fdr_alpha = 0.05,
                          normalize = FALSE, seed = 42) {
  fold_seed <- seed + fold
  tryCatch(
    {
      split <- split_train_test(dataset$labels, ratio = train_ratio, seed = fold_seed)

      # Optional leakage-free z-score: fit mean/sd on TRAIN, apply to all columns.
      expr_use <- if (isTRUE(normalize)) {
        zscore_fit_apply(dataset$expr, split$train)
      } else {
        dataset$expr
      }

      # Phase I — training-only gene shaving (Welch t-test + BH FDR), then group.
      keep_genes <- if (isTRUE(shave)) {
        shave_genes(expr_use, dataset$labels, split$train, alpha = fdr_alpha)
      } else {
        rownames(expr_use)
      }

      # Candidate groups: optional ANOVA group pre-filter, then shave their genes.
      candidates <- pre_score_groups(dataset, split$train,
        method = pre_score_method,
        top_percent = pre_score_top_percent
      )
      pre_sc <- attr(candidates, "scores")
      candidates <- shave_groups(candidates, keep_genes)
      if (length(candidates) == 0) {
        stop("No groups survived gene shaving (FDR ", fdr_alpha, ").")
      }

      # Phase II — SCORE each candidate (mean Accuracy, k-fold CV, TRAIN only) -> rank.
      scores <- vapply(names(candidates), function(name) {
        do.call(score_group_train, c(
          list(
            expr = expr_use, labels = dataset$labels,
            genes = candidates[[name]], train_idx = split$train,
            method = score_method
          ),
          score_args
        ))
      }, numeric(1))
      ranked <- names(sort(scores, decreasing = TRUE))

      # Phase III — ACCUMULATIVE MODELING top-1..K; skip k that add no new genes.
      K <- min(top_n, length(ranked))
      acc_rows <- list()
      best_model <- NULL
      prev_features <- 0L
      for (k in seq_len(K)) {
        sub <- candidates[ranked[seq_len(k)]]
        n_feat <- length(unique(unlist(sub)))
        if (n_feat == prev_features) next # expand: no new features at this k
        prev_features <- n_feat
        m <- do.call(build_model, c(
          list(
            expr = expr_use, labels = dataset$labels, top_groups = sub,
            train_idx = split$train, method = model_method
          ),
          model_args
        ))
        acc_rows[[length(acc_rows) + 1L]] <- data.frame(
          k = k, group_added = ranked[k], n_features = length(m$features),
          Accuracy = m$metrics$Accuracy, Sensitivity = m$metrics$Sensitivity,
          Specificity = m$metrics$Specificity, Precision = m$metrics$Precision,
          F_measure = m$metrics$F_measure, Cohens_kappa = m$metrics$Cohens_kappa,
          AUC = m$metrics$AUC %||% NA_real_,
          stringsAsFactors = FALSE
        )
        m$.k <- k
        # Best model = highest Accuracy (the G-S-M selection metric).
        if (is.null(best_model) || m$metrics$Accuracy > best_model$metrics$Accuracy) {
          best_model <- m
        }
      }
      accumulative <- do.call(rbind, acc_rows)
      best_k <- best_model$.k

      list(
        fold = fold,
        status = "success",
        metrics = best_model$metrics, # performance at the best k (by Accuracy)
        best_k = best_k,
        features_count = length(best_model$features),
        n_groups = best_k,
        n_shaved_genes = length(keep_genes), # genes surviving Welch+BH
        accumulative = accumulative, # top-1..K performance curve
        ranked_groups = ranked,
        group_accuracies = scores, # train-only mean-Accuracy ranking scores
        pre_scores = pre_sc, # sorted ANOVA scores (all groups)
        importance = best_model$importance, # per-gene RF importance @ best k
        # Signed training Welch-t over the group-gene universe, for G-S-M detect
        # (de novo module discovery); classifier-independent, no leakage.
        gene_stats = gene_signed_stats(
          expr_use, dataset$labels, split$train,
          universe = intersect(unique(unlist(dataset$groups)), keep_genes)
        )
      )
    },
    error = function(e) {
      list(fold = fold, status = "error", error_message = conditionMessage(e))
    }
  )
}

#' Run Monte Carlo Cross-Validation (MCCV) in parallel
#'
#' Runs multiple iterations of stratified train/test splitting, statistical
#' pre-scoring, group scoring, and final model building. Folds are executed
#' in parallel across multiple CPU cores.
#'
#' @param dataset A `gsm_dataset` object.
#' @param n_iterations Integer. Number of MCCV iterations. Default `100`.
#' @param train_ratio Numeric. Proportion of samples for training. Default `0.9`.
#' @param pre_score_method Pre-scoring method. Default `"anova"`.
#' @param pre_score_top_percent Percentage of top groups to retain. Default `1`.
#' @param score_method Classifier method for group scoring. Default `"rf"`.
#' @param model_method Classifier method for final model. Default `"rf"`.
#' @param score_args Named list of additional arguments passed to the scoring classifier. Default `list()`.
#' @param model_args Named list of additional arguments passed to the final model classifier. Default `list()`.
#' @param top_n Integer. Number of accumulative top-k levels to evaluate per fold
#'   (the model is built on top-1, top-1..2, ..., top-1..`top_n` ranked groups).
#'   Default `10`.
#' @param cores Integer. Number of CPU cores to use. Default is all minus one.
#' @param seed Integer. Base random seed for reproducibility. Default `42`.
#' @param ... Additional arguments passed to model/scoring functions.
#' @return A list of MCCV results, including aggregated metrics across all folds.
#' @export
run_mccv <- function(dataset, n_iterations = 100, train_ratio = 0.9,
                     pre_score_method = "anova", pre_score_top_percent = 1,
                     score_method = "rf", model_method = "rf",
                     score_args = list(), model_args = list(),
                     top_n = 10, shave = TRUE, fdr_alpha = 0.05,
                     normalize = FALSE, cores = NULL, seed = 44, ...) {
  stopifnot(inherits(dataset, "gsm_dataset"))
  stopifnot(is.numeric(n_iterations) && n_iterations > 0)
  stopifnot(is.numeric(train_ratio) && train_ratio > 0 && train_ratio < 1)
  stopifnot(is.list(score_args))
  stopifnot(is.list(model_args))

  if (is.null(cores)) {
    cores <- max(1, parallel::detectCores() - 1)
  }

  message(sprintf("Starting MCCV with %d iterations on %d cores...", n_iterations, cores))

  # Single fold runner delegates to the exported helper (shared with the
  # Shiny app's live/streaming MCCV view).
  run_fold <- function(i) {
    mccv_run_fold(dataset, i,
      train_ratio = train_ratio,
      pre_score_method = pre_score_method,
      pre_score_top_percent = pre_score_top_percent,
      score_method = score_method, model_method = model_method,
      score_args = score_args, model_args = model_args,
      top_n = top_n, shave = shave, fdr_alpha = fdr_alpha,
      normalize = normalize, seed = seed
    )
  }

  # Run in parallel or sequentially based on cores and OS
  results <- list()
  if (cores > 1 && .Platform$OS.type != "windows") {
    results <- parallel::mclapply(seq_len(n_iterations), run_fold, mc.cores = cores)
  } else if (cores > 1 && .Platform$OS.type == "windows") {
    # Windows fallback socket cluster
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl))

    # Export package and variables
    parallel::clusterEvalQ(cl, {
      library(gsmR)
      library(randomForest)
      library(e1071)
    })
    parallel::clusterExport(cl, c(
      "dataset", "train_ratio", "pre_score_method",
      "pre_score_top_percent", "score_method",
      "model_method", "seed", "score_args",
      "model_args"
    ), envir = environment())

    results <- parallel::parLapply(cl, seq_len(n_iterations), run_fold)
  } else {
    results <- lapply(seq_len(n_iterations), run_fold)
  }

  aggregate_mccv_folds(results)
}


#' Aggregate MCCV fold results into a summary object
#'
#' Shared by [run_mccv()] and the multi-dataset batch runner ([run_batch()]) so
#' both produce an identically-shaped result. Filters to successful folds, then
#' computes per-metric mean/SD across them.
#'
#' @param results List of per-fold results (each with `$status` and `$metrics`),
#'   as returned by [mccv_run_fold()].
#' @return A list with `folds` (all fold results, including failures),
#'   `metrics_summary` (data.frame of `Mean`/`SD` per metric), and `raw_metrics`
#'   (one row per successful fold).
#' @keywords internal
aggregate_mccv_folds <- function(results) {
  successful_results <- Filter(function(res) res$status == "success", results)

  if (length(successful_results) == 0) {
    stop("All MCCV folds failed. Check errors.")
  }

  if (length(successful_results) < length(results)) {
    warning(sprintf(
      "%d folds failed during MCCV execution.",
      length(results) - length(successful_results)
    ))
  }

  # Aggregate metrics across all successful folds
  metrics_df <- do.call(rbind, lapply(successful_results, function(res) {
    as.data.frame(res$metrics)
  }))

  mean_metrics <- colMeans(metrics_df, na.rm = TRUE)
  sd_metrics <- apply(metrics_df, 2, stats::sd, na.rm = TRUE)

  list(
    folds = results,
    metrics_summary = data.frame(
      Mean = mean_metrics,
      SD = sd_metrics
    ),
    raw_metrics = metrics_df
  )
}
