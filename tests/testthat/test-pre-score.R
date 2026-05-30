# Tests for pre-scoring functions

test_that("pre_score_anova ranks signal group highest", {
  # Build a dataset where GroupA has signal genes
  f_expr <- test_data_path("test_expression.csv")
  f_groups <- test_data_path("test_groups.gmt")

  expr_data <- read_expression(f_expr)
  groups <- read_gmt(f_groups)

  ds <- gsm_dataset(expr_data$expr, expr_data$labels, groups, min_group_size = 5)

  # Use first 14 samples as training
  train_idx <- 1:14

  # Pre-score with ANOVA, keep top 40% (should keep at least 2 groups)
  result <- pre_score_groups(ds, train_idx, method = "anova", top_percent = 40)

  expect_type(result, "list")
  expect_true(length(result) >= 1)
  # The signal group should be among the top
  expect_true("GroupA_signal" %in% names(result))
})


test_that("pre_score_none returns all groups", {
  expr <- matrix(rnorm(60), nrow = 6, ncol = 10)
  rownames(expr) <- paste0("G", 1:6)
  colnames(expr) <- paste0("S", 1:10)
  labels <- factor(rep(c("A", "B"), 5))
  names(labels) <- colnames(expr)
  groups <- list(
    g1 = paste0("G", 1:6),
    g2 = paste0("G", 1:6)
  )

  ds <- gsm_dataset(expr, labels, groups, min_group_size = 1)
  result <- pre_score_groups(ds, train_idx = 1:7, method = "none")

  expect_equal(length(result), length(ds$groups))
})


test_that("pre_score_groups errors on unknown method", {
  expr <- matrix(rnorm(60), nrow = 6, ncol = 10)
  rownames(expr) <- paste0("G", 1:6)
  colnames(expr) <- paste0("S", 1:10)
  labels <- factor(rep(c("A", "B"), 5))
  names(labels) <- colnames(expr)
  groups <- list(g1 = paste0("G", 1:6))

  ds <- gsm_dataset(expr, labels, groups, min_group_size = 1)
  expect_error(pre_score_groups(ds, 1:7, method = "bogus"), "Unknown")
})


test_that("pre_score_groups validates inputs", {
  expect_error(pre_score_groups("not_a_dataset", 1:5))
})
