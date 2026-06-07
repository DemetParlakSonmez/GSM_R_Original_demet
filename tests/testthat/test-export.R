test_that("group ranking aggregates per-fold scores correctly", {
  # Minimal fake dataset: only $groups and $params are used by builders.
  ds <- structure(
    list(
      groups = list(GroupA = c("g1", "g2"), GroupB = c("g3"),
                    GroupC = c("g4", "g5", "g6")),
      params = list(n_samples = 20L, n_genes = 6L, n_groups = 3L)
    ),
    class = "gsm_dataset"
  )
  mk_metrics <- function(a) list(Accuracy = a, Sensitivity = a, Specificity = a,
                                 Precision = a, F_measure = a, Cohens_kappa = a)
  folds <- list(
    list(fold = 1, status = "success", metrics = mk_metrics(0.9),
         features_count = 3L, n_groups = 2L,
         group_accuracies = c(GroupA = 0.9, GroupB = 0.7),
         pre_scores = c(GroupA = 10, GroupB = 5), importance = c(g1 = 2, g2 = 1)),
    list(fold = 2, status = "success", metrics = mk_metrics(0.8),
         features_count = 3L, n_groups = 2L,
         group_accuracies = c(GroupA = 0.6, GroupC = 0.8),
         pre_scores = c(GroupC = 9, GroupA = 4), importance = c(g4 = 3))
  )

  gr <- gsmR:::gsm_tbl_group_ranking(folds, ds, top_n = 1)
  a <- gr[gr$Group == "GroupA", ]
  expect_equal(a$Times_Selected, 2L)             # selected in both folds
  expect_equal(a$Mean_Accuracy, 0.75)            # (0.9 + 0.6) / 2
  # GroupA ranks 1 in fold1, rank 2 in fold2 -> mean rank 1.5
  expect_equal(a$Mean_Rank, 1.5)
  # TopN_Appearances with top_n=1: GroupA is #1 only in fold1
  expect_equal(a$TopN_Appearances, 1L)

  # Group-feature associations: GroupC (3 genes) selected once
  assoc <- gsmR:::gsm_tbl_group_feature_assoc(folds, ds)
  cc <- assoc[assoc$Group == "GroupC", ]
  expect_equal(nrow(cc), 3L)
  expect_true(all(cc$Times_Group_Selected == 1L))
})

test_that("gsm_export_tables returns the expected MCCV table set", {
  ds <- structure(
    list(groups = list(GroupA = c("g1", "g2"), GroupB = c("g3")),
         params = list(n_samples = 10L, n_genes = 3L, n_groups = 2L)),
    class = "gsm_dataset"
  )
  mk_metrics <- function(a) list(Accuracy = a, Sensitivity = a, Specificity = a,
                                 Precision = a, F_measure = a, Cohens_kappa = a)
  folds <- list(
    list(fold = 1, status = "success", metrics = mk_metrics(0.9),
         features_count = 3L, n_groups = 2L,
         group_accuracies = c(GroupA = 0.9, GroupB = 0.7),
         pre_scores = c(GroupA = 10, GroupB = 5), importance = c(g1 = 2)),
    list(fold = 2, status = "error", error_message = "boom")
  )
  result <- list(folds = folds,
                 metrics_summary = data.frame(Mean = rep(0.9, 6), SD = rep(0, 6),
                                              row.names = names(mk_metrics(0))))
  tbls <- gsm_export_tables(result, ds)
  expect_true(all(c("RunSummary", "MetricSummary", "IterationMetrics",
                    "IterationGroupScores", "GroupRanking", "GroupFeatureAssoc",
                    "PreScoreSignificance", "RankedFeatures") %in% names(tbls)))
  # Error fold is retained in IterationMetrics
  expect_equal(nrow(tbls$IterationMetrics), 2L)
  expect_true("error" %in% tbls$IterationMetrics$Status)

  # CSV bundle writes a non-empty zip
  zf <- tempfile(fileext = ".zip")
  gsm_write_csv_zip(tbls, zf, prefix = "TEST")
  expect_true(file.exists(zf) && file.size(zf) > 0)

  skip_if_not_installed("openxlsx")
  xf <- tempfile(fileext = ".xlsx")
  gsm_write_xlsx(tbls, xf)
  expect_true(file.exists(xf))
  expect_setequal(openxlsx::getSheetNames(xf), substr(names(tbls), 1, 31))
})

test_that("single-fold export produces the reduced table set", {
  ds <- structure(
    list(groups = list(GroupA = c("g1", "g2"), GroupB = c("g3")),
         params = list(n_samples = 10L, n_genes = 3L, n_groups = 2L)),
    class = "gsm_dataset"
  )
  result <- list(
    score_results = list(GroupA = list(accuracy = 0.9), GroupB = list(accuracy = 0.6)),
    final_model = list(
      metrics = list(Accuracy = 0.85, Sensitivity = 0.8, Specificity = 0.9,
                     Precision = 0.83, F_measure = 0.81, Cohens_kappa = 0.7),
      features = c("g1", "g2", "g3"), importance = c(g1 = 2, g2 = 1, g3 = 0.5)
    )
  )
  tbls <- gsm_export_tables(result, ds)
  expect_setequal(names(tbls),
                  c("RunSummary", "FinalMetrics", "GroupScores", "FeatureList"))
  expect_equal(nrow(tbls$GroupScores), 2L)
  expect_equal(tbls$GroupScores$Group[1], "GroupA")  # sorted by accuracy
})
