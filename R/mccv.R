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
#'   `metrics`, `features_count`, and `group_accuracies`.
#' @export
mccv_run_fold <- function(dataset, fold, train_ratio = 0.7,
                          pre_score_method = "anova", pre_score_top_percent = 1,
                          score_method = "rf", model_method = "rf",
                          score_args = list(), model_args = list(), seed = 42) {
  fold_seed <- seed + fold
  tryCatch({
    split <- split_train_test(dataset$labels, ratio = train_ratio, seed = fold_seed)

    top_groups <- pre_score_groups(dataset, split$train,
                                   method = pre_score_method,
                                   top_percent = pre_score_top_percent)

    group_scores <- list()
    for (name in names(top_groups)) {
      score_args_combined <- c(
        list(expr = dataset$expr, labels = dataset$labels,
             genes = top_groups[[name]], train_idx = split$train,
             method = score_method),
        score_args
      )
      score <- do.call(score_group, score_args_combined)
      group_scores[[name]] <- score$accuracy
    }

    model_args_combined <- c(
      list(expr = dataset$expr, labels = dataset$labels, top_groups = top_groups,
           train_idx = split$train, method = model_method),
      model_args
    )
    final_model <- do.call(build_model, model_args_combined)

    list(
      fold = fold,
      status = "success",
      metrics = final_model$metrics,
      features_count = length(final_model$features),
      n_groups = length(top_groups),
      group_accuracies = group_scores
    )
  }, error = function(e) {
    list(fold = fold, status = "error", error_message = conditionMessage(e))
  })
}

#' Run Monte Carlo Cross-Validation (MCCV) in parallel
#'
#' Runs multiple iterations of stratified train/test splitting, statistical
#' pre-scoring, group scoring, and final model building. Folds are executed
#' in parallel across multiple CPU cores.
#'
#' @param dataset A `gsm_dataset` object.
#' @param n_iterations Integer. Number of MCCV iterations. Default `10`.
#' @param train_ratio Numeric. Proportion of samples for training. Default `0.7`.
#' @param pre_score_method Pre-scoring method. Default `"anova"`.
#' @param pre_score_top_percent Percentage of top groups to retain. Default `1`.
#' @param score_method Classifier method for group scoring. Default `"rf"`.
#' @param model_method Classifier method for final model. Default `"rf"`.
#' @param score_args Named list of additional arguments passed to the scoring classifier. Default `list()`.
#' @param model_args Named list of additional arguments passed to the final model classifier. Default `list()`.
#' @param cores Integer. Number of CPU cores to use. Default is all minus one.
#' @param seed Integer. Base random seed for reproducibility. Default `42`.
#' @param ... Additional arguments passed to model/scoring functions.
#' @return A list of MCCV results, including aggregated metrics across all folds.
#' @export
run_mccv <- function(dataset, n_iterations = 10, train_ratio = 0.7,
                     pre_score_method = "anova", pre_score_top_percent = 1,
                     score_method = "rf", model_method = "rf",
                     score_args = list(), model_args = list(),
                     cores = NULL, seed = 42, ...) {
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
    mccv_run_fold(dataset, i, train_ratio = train_ratio,
                  pre_score_method = pre_score_method,
                  pre_score_top_percent = pre_score_top_percent,
                  score_method = score_method, model_method = model_method,
                  score_args = score_args, model_args = model_args, seed = seed)
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
    parallel::clusterExport(cl, c("dataset", "train_ratio", "pre_score_method", 
                                  "pre_score_top_percent", "score_method", 
                                  "model_method", "seed", "score_args", 
                                  "model_args"), envir = environment())
    
    results <- parallel::parLapply(cl, seq_len(n_iterations), run_fold)
  } else {
    results <- lapply(seq_len(n_iterations), run_fold)
  }

  # Filter out error folds if any occurred
  successful_results <- Filter(function(res) res$status == "success", results)
  
  if (length(successful_results) == 0) {
    stop("All MCCV folds failed. Check errors.")
  }
  
  if (length(successful_results) < length(results)) {
    warning(sprintf("%d folds failed during MCCV execution.", length(results) - length(successful_results)))
  }

  # Aggregate metrics across all successful folds
  metrics_df <- do.call(rbind, lapply(successful_results, function(res) {
    as.data.frame(res$metrics)
  }))
  
  mean_metrics <- colMeans(metrics_df, na.rm = TRUE)
  sd_metrics <- apply(metrics_df, 2, stats::sd, na.rm = TRUE)
  
  agg_results <- list(
    folds = results,
    metrics_summary = data.frame(
      Mean = mean_metrics,
      SD = sd_metrics
    ),
    raw_metrics = metrics_df
  )
  
  agg_results
}
