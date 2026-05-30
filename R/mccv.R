#' @title Parallel Monte Carlo Cross-Validation (MCCV)
#' @description Implementations for parallel execution of MCCV across available
#'   CPU cores in gsmR.
#' @name mccv
NULL

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
#' @param cores Integer. Number of CPU cores to use. Default is all minus one.
#' @param seed Integer. Base random seed for reproducibility. Default `42`.
#' @param ... Additional arguments passed to model/scoring functions.
#' @return A list of MCCV results, including aggregated metrics across all folds.
#' @export
run_mccv <- function(dataset, n_iterations = 10, train_ratio = 0.7,
                     pre_score_method = "anova", pre_score_top_percent = 1,
                     score_method = "rf", model_method = "rf",
                     cores = NULL, seed = 42, ...) {
  stopifnot(inherits(dataset, "gsm_dataset"))
  stopifnot(is.numeric(n_iterations) && n_iterations > 0)
  stopifnot(is.numeric(train_ratio) && train_ratio > 0 && train_ratio < 1)

  if (is.null(cores)) {
    cores <- max(1, parallel::detectCores() - 1)
  }
  
  message(sprintf("Starting MCCV with %d iterations on %d cores...", n_iterations, cores))

  # Define the single fold runner
  run_fold <- function(i) {
    fold_seed <- seed + i
    
    tryCatch({
      # Run split
      split <- split_train_test(dataset$labels, ratio = train_ratio, seed = fold_seed)
      
      # Pre-scoring
      top_groups <- pre_score_groups(dataset, split$train, 
                                     method = pre_score_method, 
                                     top_percent = pre_score_top_percent)
      
      # Score individual groups
      group_scores <- list()
      for (name in names(top_groups)) {
        genes <- top_groups[[name]]
        score <- score_group(dataset$expr, dataset$labels, genes, split$train, method = score_method)
        group_scores[[name]] <- score$accuracy
      }
      
      # Build final model
      final_model <- build_model(dataset$expr, dataset$labels, top_groups, split$train, method = model_method)
      
      list(
        fold = i,
        status = "success",
        metrics = final_model$metrics,
        features_count = length(final_model$features),
        group_accuracies = group_scores
      )
    }, error = function(e) {
      list(
        fold = i,
        status = "error",
        error_message = e$message
      )
    })
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
    })
    parallel::clusterExport(cl, c("dataset", "train_ratio", "pre_score_method", 
                                  "pre_score_top_percent", "score_method", 
                                  "model_method", "seed"), envir = environment())
    
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
