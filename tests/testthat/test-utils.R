# Tests for utility functions

test_that("filter_groups_by_size filters correctly", {
  groups <- list(
    small = c("A", "B"),                    # 2 genes
    medium = paste0("G", 1:10),             # 10 genes
    big = paste0("G", 1:100),               # 100 genes
    huge = paste0("G", 1:600)               # 600 genes
  )

  # Default: min=5, max=500
  expect_message(result <- filter_groups_by_size(groups))
  expect_equal(length(result), 2)
  expect_true("medium" %in% names(result))
  expect_true("big" %in% names(result))
  expect_false("small" %in% names(result))
  expect_false("huge" %in% names(result))
})


test_that("filter_groups_by_size handles edge cases", {
  groups <- list(exact5 = paste0("G", 1:5))

  # Exactly at min boundary — should keep
  result <- filter_groups_by_size(groups, min_size = 5, max_size = 500)
  expect_equal(length(result), 1)

  # Below min — should drop
  expect_message(result <- filter_groups_by_size(groups, min_size = 6, max_size = 500))
  expect_equal(length(result), 0)
})


test_that("intersect_genes returns correct intersection", {
  available <- paste0("G", 1:10)
  groups <- list(
    g1 = c("G1", "G2", "G3", "MISSING1"),
    g2 = c("G5", "MISSING2"),
    g3 = c("MISSING3", "MISSING4")   # All missing — should be dropped
  )

  result <- intersect_genes(available, groups)
  expect_equal(length(result), 2)  # g3 dropped
  expect_equal(length(result[["g1"]]), 3)
  expect_equal(length(result[["g2"]]), 1)
})


test_that("normalize_expression zscore works", {
  set.seed(1)
  expr <- matrix(rnorm(50, mean = 10, sd = 3), nrow = 5, ncol = 10)
  rownames(expr) <- paste0("G", 1:5)

  normed <- normalize_expression(expr, method = "zscore")

  # Each row should have mean ≈ 0 and sd ≈ 1
  row_means <- rowMeans(normed)
  row_sds <- apply(normed, 1, sd)

  expect_equal(unname(row_means), rep(0, 5), tolerance = 1e-10)
  expect_equal(unname(row_sds), rep(1, 5), tolerance = 1e-10)
})


test_that("normalize_expression handles constant genes", {
  expr <- matrix(c(5, 5, 5, 5, 5, 1, 2, 3, 4, 5), nrow = 2, ncol = 5, byrow = TRUE)
  rownames(expr) <- c("constant", "varying")

  # Should not error on constant row
  normed <- normalize_expression(expr)
  expect_false(any(is.nan(normed)))
  expect_equal(normed["constant", ], rep(0, 5))  # (5-5)/1 = 0
})


test_that("normalize_expression none method passes through", {
  expr <- matrix(1:12, nrow = 3, ncol = 4)
  result <- normalize_expression(expr, method = "none")
  expect_identical(result, expr)
})


test_that("split_train_test produces correct ratio", {
  labels <- factor(c(rep("A", 30), rep("B", 70)))
  split <- split_train_test(labels, ratio = 0.7, seed = 42)

  # Total should equal n

  expect_equal(length(split$train) + length(split$test), 100)

  # Approximately 70% training
  expect_true(length(split$train) >= 60 && length(split$train) <= 80)

  # No overlap
  expect_equal(length(intersect(split$train, split$test)), 0)
})


test_that("split_train_test is reproducible with seed", {
  labels <- factor(c(rep("A", 50), rep("B", 50)))

  split1 <- split_train_test(labels, ratio = 0.7, seed = 123)
  split2 <- split_train_test(labels, ratio = 0.7, seed = 123)

  expect_equal(split1$train, split2$train)
  expect_equal(split1$test, split2$test)
})


test_that("split_train_test is stratified", {
  labels <- factor(c(rep("A", 20), rep("B", 80)))
  split <- split_train_test(labels, ratio = 0.7, seed = 42)

  # Both classes should be in training set
  train_labels <- labels[split$train]
  expect_true("A" %in% train_labels)
  expect_true("B" %in% train_labels)

  # Both classes should be in test set
  test_labels <- labels[split$test]
  expect_true("A" %in% test_labels)
  expect_true("B" %in% test_labels)
})


test_that("check_sample_match detects mismatch", {
  expr <- matrix(1:12, nrow = 3, ncol = 4)
  colnames(expr) <- paste0("S", 1:4)
  labels <- factor(c("A", "B", "A", "B"))
  names(labels) <- paste0("S", 1:4)

  # Should pass
  expect_true(gsmR:::check_sample_match(expr, labels))

  # Mismatch
  names(labels) <- paste0("X", 1:4)
  expect_error(gsmR:::check_sample_match(expr, labels), "do not match")
})
