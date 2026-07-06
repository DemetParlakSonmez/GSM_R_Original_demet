# Tests for label-blind expression filtering (R/filter.R).

.make_expr <- function(p = 200, n = 12, seed = 1) {
  set.seed(seed)
  X <- matrix(stats::rnorm(p * n), nrow = p,
              dimnames = list(paste0("G", seq_len(p)), paste0("S", seq_len(n))))
  X[1:10, ] <- 3                      # 10 constant genes
  X[11:20, ] <- X[11:20, ] * 0.001    # 10 near-zero-variance / low genes
  X
}

test_that("filter_genes drops constant genes and keeps names", {
  X <- .make_expr()
  Xf <- filter_genes(X)
  expect_false(any(paste0("G", 1:10) %in% rownames(Xf)))   # constants gone
  expect_true(all(rownames(Xf) %in% rownames(X)))          # names preserved
  expect_equal(ncol(Xf), ncol(X))                          # samples unchanged
  rep <- attr(Xf, "filter_report")
  expect_equal(rep[["input"]], nrow(X))
  expect_lt(rep[["after_missing"]], nrow(X))
})

test_that("max_genes keeps the most variable genes only", {
  X <- .make_expr()
  Xf <- filter_genes(X, max_genes = 50)
  expect_lte(nrow(Xf), 50)
  # The highest-variance genes should survive; the tiny-variance block should not.
  expect_false(any(paste0("G", 11:20) %in% rownames(Xf)))
})

test_that("min_mean drops low-expression genes", {
  X <- .make_expr()
  X[21:40, ] <- X[21:40, ] - 100      # very low mean
  Xf <- filter_genes(X, min_mean = -5)
  expect_false(any(paste0("G", 21:40) %in% rownames(Xf)))
})

test_that("filter is label-blind (independent of any y)", {
  X <- .make_expr()
  f1 <- rownames(filter_genes(X, max_genes = 80))
  # There is no y argument at all; a permutation of a hypothetical label vector
  # cannot change the result. Re-running must be identical.
  f2 <- rownames(filter_genes(X, max_genes = 80))
  expect_identical(f1, f2)
})

test_that("filter_genes errors clearly when nothing survives", {
  X <- matrix(5, nrow = 5, ncol = 4,
              dimnames = list(paste0("G", 1:5), paste0("S", 1:4)))  # all constant
  expect_error(filter_genes(X), "removed almost all genes")
})
