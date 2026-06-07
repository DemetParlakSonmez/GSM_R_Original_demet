# Tests for parallel Monte Carlo Cross-Validation (MCCV)

test_that("run_mccv executes successfully on synthetic data", {
  f_expr <- test_data_path("test_expression.csv")
  f_groups <- test_data_path("test_groups.csv")

  expr_data <- read_expression(f_expr)
  groups <- read_groups(f_groups)

  # Build dataset
  ds <- gsm_dataset(expr_data$expr, expr_data$labels, groups, min_group_size = 5)

  # Run MCCV with 2 iterations, 1 core (sequential)
  mccv_seq <- run_mccv(ds, n_iterations = 2, cores = 1, seed = 42, pre_score_top_percent = 20)

  expect_type(mccv_seq, "list")
  expect_true("folds" %in% names(mccv_seq))
  expect_true("metrics_summary" %in% names(mccv_seq))
  expect_true("raw_metrics" %in% names(mccv_seq))

  expect_equal(length(mccv_seq$folds), 2)
  expect_equal(nrow(mccv_seq$raw_metrics), 2)

  # Check that each fold succeeded
  expect_equal(mccv_seq$folds[[1]]$status, "success")
  expect_equal(mccv_seq$folds[[2]]$status, "success")

  # Aggregated metrics checks
  expect_true("Accuracy" %in% rownames(mccv_seq$metrics_summary))
  expect_true(is.numeric(mccv_seq$metrics_summary$Mean))
  expect_true(all(mccv_seq$metrics_summary$Mean >= 0 & mccv_seq$metrics_summary$Mean <= 1))
})


test_that("run_mccv is reproducible with same seed", {
  f_expr <- test_data_path("test_expression.csv")
  f_groups <- test_data_path("test_groups.csv")

  expr_data <- read_expression(f_expr)
  groups <- read_groups(f_groups)
  ds <- gsm_dataset(expr_data$expr, expr_data$labels, groups, min_group_size = 5)

  # Run twice with same seed
  mccv1 <- run_mccv(ds, n_iterations = 2, cores = 1, seed = 123, pre_score_top_percent = 20)
  mccv2 <- run_mccv(ds, n_iterations = 2, cores = 1, seed = 123, pre_score_top_percent = 20)

  expect_equal(mccv1$raw_metrics, mccv2$raw_metrics)
  expect_equal(mccv1$metrics_summary, mccv2$metrics_summary)
})


test_that("run_mccv works in parallel", {
  f_expr <- test_data_path("test_expression.csv")
  f_groups <- test_data_path("test_groups.csv")

  expr_data <- read_expression(f_expr)
  groups <- read_groups(f_groups)
  ds <- gsm_dataset(expr_data$expr, expr_data$labels, groups, min_group_size = 5)

  # Run MCCV with 2 iterations on 2 cores
  mccv_par <- run_mccv(ds, n_iterations = 2, cores = 2, seed = 42, pre_score_top_percent = 20)

  expect_equal(length(mccv_par$folds), 2)
  expect_equal(mccv_par$folds[[1]]$status, "success")
  expect_equal(mccv_par$folds[[2]]$status, "success")
})


test_that("mccv_run_fold builds an accumulative top-k curve ranked by train score", {
  skip_if_not_installed("randomForest")
  expr_data <- read_expression(test_data_path("test_expression.csv"))
  groups <- read_groups(test_data_path("test_groups.csv"))
  ds <- gsm_dataset(expr_data$expr, expr_data$labels, groups,
                    min_group_size = 2, max_group_size = 500)

  f <- mccv_run_fold(ds, fold = 1, train_ratio = 0.7, pre_score_method = "none",
                     score_method = "rf", model_method = "rf", top_n = 4, seed = 1)
  expect_equal(f$status, "success")
  expect_s3_class(f$accumulative, "data.frame")
  expect_true(nrow(f$accumulative) >= 1 && nrow(f$accumulative) <= 4)
  expect_equal(f$accumulative$k, seq_len(nrow(f$accumulative)))
  expect_true(all(c("Accuracy", "Sensitivity", "Cohens_kappa") %in%
                    names(f$accumulative)))
  expect_false(is.unsorted(f$accumulative$n_features))     # genes accumulate
  expect_equal(f$best_k, which.max(f$accumulative$Accuracy))
  # ranked_groups follows the train-only group scores (descending).
  expect_equal(f$ranked_groups, names(sort(f$group_accuracies, decreasing = TRUE)))
  expect_true("GroupA_signal" %in% f$ranked_groups[1:2])   # signal group ranks top
})


test_that("run_mccv aggregates accumulative folds for export (res1 per k)", {
  skip_if_not_installed("randomForest")
  expr_data <- read_expression(test_data_path("test_expression.csv"))
  groups <- read_groups(test_data_path("test_groups.csv"))
  ds <- gsm_dataset(expr_data$expr, expr_data$labels, groups,
                    min_group_size = 2, max_group_size = 500)

  res <- run_mccv(ds, n_iterations = 3, top_n = 3, cores = 1, seed = 42,
                  pre_score_top_percent = 100)
  acc <- gsm_export_tables(res, ds)$AccumulativeSummary
  expect_equal(nrow(acc), 3L)            # one row per k = 1..3
  expect_true(all(c("Mean_Accuracy", "SD_Accuracy", "Mean_Features") %in% names(acc)))
})
