# Tests for pluggable classifiers (RF, SVM, KNN)

test_that("score_group and build_model work with Random Forest", {
  f_expr <- test_data_path("test_expression.csv")
  f_groups <- test_data_path("test_groups.csv")
  expr_data <- read_expression(f_expr)
  groups <- read_groups(f_groups)
  ds <- gsm_dataset(expr_data$expr, expr_data$labels, groups, min_group_size = 5)

  train_idx <- 1:14
  genes <- ds$groups[[1]]

  # score_group RF
  score <- score_group(ds$expr, ds$labels, genes, train_idx, method = "rf", ntree = 50)
  expect_type(score, "list")
  expect_true("accuracy" %in% names(score))
  expect_true(score$accuracy >= 0 && score$accuracy <= 1)

  # build_model RF
  model <- build_model(ds$expr, ds$labels, ds$groups[1:2], train_idx, method = "rf", ntree = 50)
  expect_type(model, "list")
  expect_true("metrics" %in% names(model))
  expect_true(model$metrics$Accuracy >= 0 && model$metrics$Accuracy <= 1)
})


test_that("score_group and build_model work with SVM", {
  f_expr <- test_data_path("test_expression.csv")
  f_groups <- test_data_path("test_groups.csv")
  expr_data <- read_expression(f_expr)
  groups <- read_groups(f_groups)
  ds <- gsm_dataset(expr_data$expr, expr_data$labels, groups, min_group_size = 5)

  train_idx <- 1:14
  genes <- ds$groups[[1]]

  # score_group SVM (linear kernel, cost = 1)
  score <- score_group(ds$expr, ds$labels, genes, train_idx, method = "svm", kernel = "linear", cost = 1)
  expect_type(score, "list")
  expect_true("accuracy" %in% names(score))
  expect_true(score$accuracy >= 0 && score$accuracy <= 1)

  # build_model SVM
  model <- build_model(ds$expr, ds$labels, ds$groups[1:2], train_idx, method = "svm", kernel = "radial")
  expect_type(model, "list")
  expect_true("metrics" %in% names(model))
  expect_true(model$metrics$Accuracy >= 0 && model$metrics$Accuracy <= 1)
})


test_that("score_group and build_model work with KNN", {
  f_expr <- test_data_path("test_expression.csv")
  f_groups <- test_data_path("test_groups.csv")
  expr_data <- read_expression(f_expr)
  groups <- read_groups(f_groups)
  ds <- gsm_dataset(expr_data$expr, expr_data$labels, groups, min_group_size = 5)

  train_idx <- 1:14
  genes <- ds$groups[[1]]

  # score_group KNN (k = 3)
  score <- score_group(ds$expr, ds$labels, genes, train_idx, method = "knn", k = 3)
  expect_type(score, "list")
  expect_true("accuracy" %in% names(score))
  expect_true(score$accuracy >= 0 && score$accuracy <= 1)

  # build_model KNN
  model <- build_model(ds$expr, ds$labels, ds$groups[1:2], train_idx, method = "knn", k = 1)
  expect_type(model, "list")
  expect_true("metrics" %in% names(model))
  expect_true(model$metrics$Accuracy >= 0 && model$metrics$Accuracy <= 1)
})
