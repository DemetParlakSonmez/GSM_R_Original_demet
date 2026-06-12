# Tests for consensus biomarker (gene) ranking from MCCV

bm_ds <- function() {
  expr_data <- read_expression(test_data_path("test_expression.csv"))
  groups <- read_gmt(test_data_path("test_groups.gmt"))
  gsm_dataset(expr_data$expr, expr_data$labels, groups, min_group_size = 5)
}

test_that("gsm_biomarkers aggregates RF importance into a ranked table", {
  ds <- bm_ds()
  res <- run_mccv(ds, n_iterations = 3, train_ratio = 0.7,
                  pre_score_method = "none", pre_score_top_percent = 100,
                  score_method = "rf", model_method = "rf",
                  model_args = list(ntree = 50), top_n = 3,
                  shave = FALSE, cores = 1, seed = 11)

  bm <- gsm_biomarkers(res)
  expect_s3_class(bm, "data.frame")
  expect_true(all(c("Gene", "Frequency", "MeanImportance", "Rank") %in% names(bm)))
  expect_true(nrow(bm) >= 1)
  expect_true(all(bm$Frequency > 0 & bm$Frequency <= 1))
  expect_equal(bm$Rank, seq_len(nrow(bm)))
  # top_n truncates
  expect_lte(nrow(gsm_biomarkers(res, top_n = 2)), 2)
})

test_that("gsm_biomarkers returns NULL when no importances are available (KNN)", {
  ds <- bm_ds()
  res <- run_mccv(ds, n_iterations = 2, train_ratio = 0.7,
                  pre_score_method = "none", pre_score_top_percent = 100,
                  score_method = "knn", model_method = "knn",
                  score_args = list(k = 3), model_args = list(k = 3),
                  top_n = 3, shave = FALSE, cores = 1, seed = 5)
  expect_message(out <- gsm_biomarkers(res), "No per-gene importances")
  expect_null(out)
})
