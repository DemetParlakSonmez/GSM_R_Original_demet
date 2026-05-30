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
