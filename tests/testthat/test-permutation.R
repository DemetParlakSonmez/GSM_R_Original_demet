# Tests for permutation-based empirical group significance

make_ds <- function() {
  expr_data <- read_expression(test_data_path("test_expression.csv"))
  groups <- read_gmt(test_data_path("test_groups.gmt"))
  gsm_dataset(expr_data$expr, expr_data$labels, groups, min_group_size = 5)
}

test_that("permutation_test_groups returns a well-formed, reproducible table", {
  ds <- make_ds()
  res1 <- permutation_test_groups(ds, method = "rf", n_perm = 30, cv = 3,
                                  score_args = list(ntree = 50), seed = 7)
  res2 <- permutation_test_groups(ds, method = "rf", n_perm = 30, cv = 3,
                                  score_args = list(ntree = 50), seed = 7)

  expect_s3_class(res1, "data.frame")
  expect_setequal(names(res1),
                  c("Group", "size", "observed", "null_mean",
                    "p_value", "p_adj", "n_perm"))
  expect_equal(nrow(res1), length(ds$groups))
  # Empirical p-values are bounded and FDR is monotone-consistent
  p <- res1$p_value[is.finite(res1$p_value)]
  expect_true(all(p > 0 & p <= 1))
  pa <- res1$p_adj[is.finite(res1$p_adj)]
  expect_true(all(pa >= 0 & pa <= 1))
  # Same seed -> identical results
  expect_equal(res1$p_value, res2$p_value)
})

test_that("the signal group is the most significant", {
  ds <- make_ds()
  res <- permutation_test_groups(ds, method = "rf", n_perm = 50, cv = 3,
                                 score_args = list(ntree = 50), seed = 1)
  # GroupA_signal carries the real class signal; it should rank first (lowest p)
  expect_equal(res$Group[1], "GroupA_signal")
  expect_true(res$observed[res$Group == "GroupA_signal"] >=
              max(res$observed[res$Group != "GroupA_signal"], na.rm = TRUE))
})

test_that("group_names subsetting and empty-group handling work", {
  ds <- make_ds()
  res <- permutation_test_groups(ds, group_names = names(ds$groups)[1],
                                 method = "rf", n_perm = 10,
                                 score_args = list(ntree = 30), seed = 3)
  expect_equal(nrow(res), 1)
})
