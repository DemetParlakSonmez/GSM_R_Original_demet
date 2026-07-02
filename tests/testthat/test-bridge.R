# Tests for the BRIDGE grouping engine (R/bridge.R) and dispatch (R/grouping.R).
#
# The essential contract: BRIDGE must return groups in the SAME shape the
# file readers return (a named list of character vectors, genes present in the
# expression matrix), so the Scoring/Modeling steps consume it unchanged.

# Simulate a tiny expression matrix with an embedded differential module.
.make_bridge_data <- function(p = 30, n = 60, module = paste0("G", 1:8),
                              seed = 42) {
  set.seed(seed)
  X <- matrix(stats::rnorm(p * n), nrow = p,
              dimnames = list(paste0("G", seq_len(p)),
                              paste0("S", seq_len(n))))
  y <- factor(rep(c("A", "B"), each = n / 2))
  # Shift + correlate the module in class B so it forms a differential block.
  shared <- stats::rnorm(sum(y == "B"))
  for (g in module) X[g, y == "B"] <- X[g, y == "B"] + 1.5 + shared
  names(y) <- colnames(X)
  list(X = X, y = y, module = module)
}

.has_backend <- function() {
  requireNamespace("glasso", quietly = TRUE) ||
    requireNamespace("JGL", quietly = TRUE)
}


test_that("bridge_prior_from_groups builds a valid [0,1] symmetric matrix", {
  skip_if_not_installed("Matrix")
  feats <- paste0("G", 1:10)
  groups <- list(m1 = paste0("G", 1:4), m2 = paste0("G", c(3, 4, 5, 6)))
  B <- bridge_prior_from_groups(groups, feats)

  expect_equal(dim(B), c(10, 10))
  expect_equal(rownames(B), feats)
  expect_true(isSymmetric(as.matrix(B)))
  expect_true(all(as.matrix(B) >= 0 & as.matrix(B) <= 1))
  expect_true(all(Matrix::diag(B) == 0))
  # Genes co-occurring in a group must have a positive prior edge.
  expect_gt(as.matrix(B)["G3", "G4"], 0)
  # Genes never co-occurring must not.
  expect_equal(as.matrix(B)["G1", "G10"], 0)
})


test_that("bridge_prior_from_edgelist aligns an edge list to features", {
  skip_if_not_installed("Matrix")
  feats <- paste0("G", 1:6)
  el <- data.frame(a = c("G1", "G2"), b = c("G2", "G3"), w = c(0.9, 0.4))
  B <- bridge_prior_from_edgelist(el, feats)

  expect_equal(dim(B), c(6, 6))
  expect_equal(as.matrix(B)["G1", "G2"], 0.9)
  expect_equal(as.matrix(B)["G2", "G3"], 0.4)
  expect_true(isSymmetric(as.matrix(B)))
})


test_that("bridge_grouping returns groups in the correct format", {
  skip_if_not_installed("Matrix")
  skip_if_not_installed("igraph")
  skip_if(!.has_backend(), "no graphical-lasso backend (glasso/JGL) installed")

  d <- .make_bridge_data()
  g <- bridge_grouping(d$X, d$y, B = list(known = d$module),
                       resolution = c(1, 1.5), tau = 0.4,
                       min_size = 3, verbose = FALSE)

  # ── The format contract ──────────────────────────────────────────────────
  expect_type(g, "list")
  if (length(g) > 0) {
    expect_false(is.null(names(g)))
    expect_true(all(nzchar(names(g))))
    expect_false(any(duplicated(names(g))))
    for (m in g) {
      expect_type(m, "character")
      expect_gt(length(m), 0)
      # Every gene must exist in the expression matrix.
      expect_true(all(m %in% rownames(d$X)))
    }
    # A grouping this shape must survive gsm_dataset() unchanged.
    ds <- gsm_dataset(d$X, d$y, g, min_group_size = 3, max_group_size = 500)
    expect_s3_class(ds, "gsm_dataset")
    expect_true(is.list(ds$groups))
  }
})


test_that("bridge_grouping works with B = NULL (data-driven) and warns", {
  skip_if_not_installed("igraph")
  skip_if(!.has_backend(), "no graphical-lasso backend installed")

  d <- .make_bridge_data()
  expect_warning(
    g <- bridge_grouping(d$X, d$y, B = NULL, resolution = 1, tau = 0.4,
                         min_size = 3, verbose = FALSE),
    "data-driven"
  )
  expect_type(g, "list")
})


test_that("bridge_grouping fails loudly on B / X feature mismatch", {
  skip_if_not_installed("Matrix")
  d <- .make_bridge_data(p = 10, n = 20)
  # Prior matrix over completely different feature names.
  Bad <- matrix(0, 5, 5,
                dimnames = list(paste0("Z", 1:5), paste0("Z", 1:5)))
  expect_error(
    bridge_grouping(d$X, d$y, B = Bad, verbose = FALSE),
    "no feature names|share no feature|do not overlap"
  )
})


test_that("bridge_grouping rejects non-binary y", {
  d <- .make_bridge_data(p = 10, n = 30)
  y3 <- factor(rep(c("A", "B", "C"), each = 10))
  names(y3) <- colnames(d$X)
  expect_error(bridge_grouping(d$X, y3, verbose = FALSE),
               "exactly two classes")
})


test_that("data_type = 'compositional' errors with SpiecEasi guidance", {
  d <- .make_bridge_data(p = 10, n = 30)
  expect_error(
    bridge_grouping(d$X, d$y, B = list(m = paste0("G", 1:4)),
                    data_type = "compositional", verbose = FALSE),
    "SpiecEasi"
  )
})


test_that("compute_grouping dispatches file method unchanged", {
  f <- test_data_path("test_groups.csv")
  via_dispatch <- compute_grouping("file", file = f, format = "csv")
  direct <- read_groups(f)
  expect_identical(via_dispatch, direct)
})


test_that("compute_grouping dispatches bridge method", {
  skip_if_not_installed("igraph")
  skip_if(!.has_backend(), "no graphical-lasso backend installed")

  d <- .make_bridge_data()
  g <- compute_grouping("bridge", X = d$X, y = d$y,
                        prior_source = "groups",
                        prior_groups = list(known = d$module),
                        bridge_args = list(resolution = 1, tau = 0.4,
                                           min_size = 3, verbose = FALSE))
  expect_type(g, "list")
})
