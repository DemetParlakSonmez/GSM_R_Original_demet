# Tests for gsm_dataset construction and validation

test_that("gsm_dataset constructs from synthetic data", {
  f_expr <- test_data_path("test_expression.csv")
  f_groups <- test_data_path("test_groups.csv")

  expr_data <- read_expression(f_expr)
  groups <- read_groups(f_groups)

  # min_group_size = 5 to keep all groups (smallest has 10)
  ds <- gsm_dataset(expr_data$expr, expr_data$labels, groups, min_group_size = 5)

  expect_s3_class(ds, "gsm_dataset")
  expect_true(is.matrix(ds$expr))
  expect_s3_class(ds$labels, "factor")
  expect_type(ds$groups, "list")
  expect_true(length(ds$groups) >= 4)
  expect_true(ds$params$n_samples == 20)
  expect_true(ds$params$n_genes == 50)
})


test_that("gsm_dataset drops small groups", {
  expr <- matrix(rnorm(100), nrow = 10, ncol = 10)
  rownames(expr) <- paste0("G", 1:10)
  colnames(expr) <- paste0("S", 1:10)
  labels <- factor(rep(c("A", "B"), 5))
  names(labels) <- colnames(expr)

  groups <- list(
    big = paste0("G", 1:8),    # 8 genes — keeps
    small = paste0("G", 1:2)   # 2 genes — dropped at min_group_size=5
  )

  expect_message(
    ds <- gsm_dataset(expr, labels, groups, min_group_size = 5, max_group_size = 500),
    "Dropped"
  )
  expect_equal(length(ds$groups), 1)
  expect_true("big" %in% names(ds$groups))
})


test_that("gsm_dataset errors on sample mismatch", {
  expr <- matrix(1:12, nrow = 3, ncol = 4)
  rownames(expr) <- paste0("G", 1:3)
  colnames(expr) <- paste0("S", 1:4)
  labels <- factor(c("A", "B", "A"))  # Only 3 samples, expr has 4
  names(labels) <- paste0("S", 1:3)

  groups <- list(g1 = paste0("G", 1:3))

  expect_error(gsm_dataset(expr, labels, groups), "does not match")
})


test_that("gsm_dataset warns about missing genes in groups", {
  expr <- matrix(rnorm(60), nrow = 6, ncol = 10)
  rownames(expr) <- paste0("G", 1:6)
  colnames(expr) <- paste0("S", 1:10)
  labels <- factor(rep(c("A", "B"), 5))
  names(labels) <- colnames(expr)

  # Group references genes not in expression matrix
  groups <- list(
    g1 = c(paste0("G", 1:6), "MISSING1", "MISSING2", "MISSING3", "MISSING4")
  )

  expect_message(
    ds <- gsm_dataset(expr, labels, groups, min_group_size = 1),
    "not in the expression matrix"
  )
  # Group should have only the 6 available genes
  expect_equal(length(ds$groups[["g1"]]), 6)
})


test_that("validate_gsm_dataset catches invalid objects", {
  expect_error(validate_gsm_dataset(list()), "not a gsm_dataset")

  # Create a valid dataset then break it
  expr <- matrix(rnorm(60), nrow = 6, ncol = 10)
  rownames(expr) <- paste0("G", 1:6)
  colnames(expr) <- paste0("S", 1:10)
  labels <- factor(rep(c("A", "B"), 5))
  names(labels) <- colnames(expr)
  groups <- list(g1 = paste0("G", 1:6))

  ds <- gsm_dataset(expr, labels, groups, min_group_size = 1)
  expect_true(validate_gsm_dataset(ds))

  # Break it
  ds$labels <- ds$labels[1:5]
  expect_error(validate_gsm_dataset(ds), "does not match")
})


test_that("print.gsm_dataset produces output", {
  expr <- matrix(rnorm(60), nrow = 6, ncol = 10)
  rownames(expr) <- paste0("G", 1:6)
  colnames(expr) <- paste0("S", 1:10)
  labels <- factor(rep(c("A", "B"), 5))
  names(labels) <- colnames(expr)
  groups <- list(g1 = paste0("G", 1:6))

  ds <- gsm_dataset(expr, labels, groups, min_group_size = 1)
  output <- capture.output(print(ds))
  expect_true(any(grepl("GSM Dataset", output)))
  expect_true(any(grepl("Samples:", output)))
})
