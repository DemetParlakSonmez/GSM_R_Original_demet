# Tests for added backends: AUC/importance wiring, XGBoost, limma pre-scoring

backend_ds <- function() {
  expr_data <- read_expression(test_data_path("test_expression.csv"))
  groups <- read_gmt(test_data_path("test_groups.gmt"))
  gsm_dataset(expr_data$expr, expr_data$labels, groups, min_group_size = 5)
}

test_that("single-group scorers now report AUC and probabilities", {
  ds <- backend_ds()
  genes <- ds$groups[["GroupA_signal"]]
  for (m in c("rf", "svm", "knn")) {
    args <- switch(m, rf = list(ntree = 50), svm = list(kernel = "linear"),
                   knn = list(k = 3))
    s <- do.call(score_group, c(list(ds$expr, ds$labels, genes, 1:14, method = m), args))
    expect_true(all(c("accuracy", "auc", "prob") %in% names(s)))
    if (!is.na(s$auc)) expect_true(s$auc >= 0 && s$auc <= 1)
  }
})

test_that("linear-SVM build_model returns per-gene importance", {
  ds <- backend_ds()
  m <- build_model(ds$expr, ds$labels, ds$groups[1:2], 1:14,
                   method = "svm", kernel = "linear")
  expect_false(is.null(m$importance))
  expect_true(all(m$importance >= 0, na.rm = TRUE))
  # Radial kernel has no linear weights -> NULL importance
  mr <- build_model(ds$expr, ds$labels, ds$groups[1:2], 1:14,
                    method = "svm", kernel = "radial")
  expect_null(mr$importance)
})

test_that("XGBoost scoring and modeling work", {
  skip_if_not_installed("xgboost")
  ds <- backend_ds()
  s <- score_group(ds$expr, ds$labels, ds$groups[["GroupA_signal"]], 1:14,
                   method = "xgb", nrounds = 20)
  expect_true(s$accuracy >= 0 && s$accuracy <= 1)
  expect_true("prob" %in% names(s))

  m <- build_model(ds$expr, ds$labels, ds$groups[1:2], 1:14,
                   method = "xgb", nrounds = 20)
  expect_true(m$metrics$Accuracy >= 0 && m$metrics$Accuracy <= 1)
  expect_false(is.null(m$importance))   # XGBoost emits per-gene Gain
})

test_that("limma pre-scoring ranks the signal group highly", {
  skip_if_not_installed("limma")
  ds <- backend_ds()
  res <- pre_score_groups(ds, train_idx = 1:14, method = "limma", top_percent = 40)
  expect_type(res, "list")
  expect_true("GroupA_signal" %in% names(res))
  expect_true(!is.null(attr(res, "scores")))
})
