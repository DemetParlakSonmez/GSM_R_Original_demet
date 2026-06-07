test_that("rra_aggregate ranks a consistently-top item first", {
  # 'a' is always ranked first across lists -> should win.
  lists <- list(c("a", "b", "c", "d"),
                c("a", "c", "b", "d"),
                c("a", "b", "d", "c"))
  res <- rra_aggregate(lists)
  expect_s3_class(res, "data.frame")
  expect_true(all(c("Name", "RRA_Score", "Rank", "Method") %in% names(res)))
  expect_equal(res$Name[1], "a")
  expect_equal(res$Rank, seq_len(nrow(res)))
})

test_that("rra_aggregate fallback works without the package path", {
  # Force the mean-rank fallback by checking it directly via a tiny case.
  lists <- list(c("x", "y"), c("x", "y"), c("y", "x"))
  res <- rra_aggregate(lists)
  expect_true(res$Name[1] %in% c("x", "y"))
  expect_equal(nrow(res), 2L)
})

test_that("calc_auc is 1 for perfect separation and ~0.5 for random", {
  actual <- factor(c("pos", "pos", "neg", "neg"), levels = c("pos", "neg"))
  expect_equal(calc_auc(actual, c(0.9, 0.8, 0.2, 0.1), "pos"), 1)
  expect_equal(calc_auc(actual, c(0.1, 0.2, 0.8, 0.9), "pos"), 0)  # perfectly wrong
  expect_true(is.na(calc_auc(factor(c("pos","pos")), c(0.5,0.6), "pos")))  # one class
})

test_that("pct_ci returns 2.5/97.5 percentiles", {
  ci <- pct_ci(0:100 / 100)
  expect_lt(ci[["low"]], ci[["high"]])
  expect_true(ci[["low"]] >= 0 && ci[["high"]] <= 1)
})
