# Tests for the combined cross-dataset report

# Build two small per-dataset result objects (as run_batch stores them) by
# running run_mccv on two synthetic datasets.
make_combined_results <- function() {
  mk_ds <- function(seed, shift) {
    set.seed(seed)
    n <- 24; ng <- 30
    labels <- factor(rep(c("A", "B"), each = n / 2))
    m <- matrix(stats::rnorm(n * ng), nrow = ng,
                dimnames = list(paste0("G", seq_len(ng)), paste0("S", seq_len(n))))
    m[1:6, labels == "B"] <- m[1:6, labels == "B"] + shift
    names(labels) <- colnames(m)
    groups <- list(sig = paste0("G", 1:6), mix = paste0("G", 4:10),
                   noise = paste0("G", 20:27))
    gsm_dataset(m, labels, groups, min_group_size = 5)
  }
  run1 <- run_mccv(mk_ds(1, 2.0), n_iterations = 3, pre_score_top_percent = 100,
                   score_method = "rf", model_method = "rf",
                   score_args = list(ntree = 50), top_n = 3, cores = 1, seed = 44)
  run2 <- run_mccv(mk_ds(2, 2.5), n_iterations = 3, pre_score_top_percent = 100,
                   score_method = "rf", model_method = "rf",
                   score_args = list(ntree = 50), top_n = 3, cores = 1, seed = 44)
  list(
    dsA = list(name = "dsA", mccv = run1, metrics_summary = run1$metrics_summary,
               n_ok = sum(vapply(run1$folds, function(r) identical(r$status, "success"),
                                 logical(1)))),
    dsB = list(name = "dsB", mccv = run2, metrics_summary = run2$metrics_summary,
               n_ok = sum(vapply(run2$folds, function(r) identical(r$status, "success"),
                                 logical(1))))
  )
}

test_that("gsm_combined_metrics returns one row per dataset plus a meta row", {
  res <- make_combined_results()
  tb <- gsm_combined_metrics(res)
  expect_s3_class(tb, "data.frame")
  expect_equal(nrow(tb), 3)                       # dsA, dsB, ALL
  expect_equal(tb$Dataset[3], "ALL")
  expect_true(all(c("Mean_Accuracy", "SD_Accuracy", "Mean_F_measure") %in% names(tb)))
  # meta mean equals the average of the two dataset means (values rounded to 4 dp)
  expect_equal(tb$Mean_Accuracy[3], mean(tb$Mean_Accuracy[1:2]), tolerance = 1e-3)
})

test_that("consensus group ranking aggregates across datasets via RRA", {
  res <- make_combined_results()
  gc_tbl <- gsm_combined_group_consensus(res)
  expect_s3_class(gc_tbl, "data.frame")
  expect_true(all(c("Group", "RRA_Score", "Rank", "Dataset_Frequency",
                    "N_Datasets") %in% names(gc_tbl)))
  expect_equal(max(gc_tbl$N_Datasets), 2)
  expect_false(is.unsorted(gc_tbl$Rank))
})

test_that("consensus biomarkers aggregate gene rankings across datasets", {
  res <- make_combined_results()
  bm <- gsm_combined_biomarkers(res)            # RF → importances available
  expect_s3_class(bm, "data.frame")
  expect_true(all(c("Gene", "RRA_Score", "Rank", "Dataset_Frequency") %in% names(bm)))
  expect_true(nrow(bm) >= 1)
})

test_that("gsm_export_combined_report writes a csv.zip bundle", {
  res <- make_combined_results()
  f <- file.path(tempdir(), "combined_test.zip")
  gsm_export_combined_report(res, f, format = "csv")
  expect_true(file.exists(f))
  inside <- utils::unzip(f, list = TRUE)$Name
  expect_true(any(grepl("PerDatasetMetrics", inside)))
  expect_true(any(grepl("GroupConsensus", inside)))
})
