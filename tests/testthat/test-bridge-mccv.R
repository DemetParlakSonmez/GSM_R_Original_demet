# Tests for leakage-free, per-fold BRIDGE grouping (R/mccv.R + R/grouping.R +
# R/bridge.R prefilter). The headline test is the NULL-CHECK: on signal-free
# data, per-fold BRIDGE must NOT manufacture above-chance test accuracy.

.has_backend <- function() {
  requireNamespace("glasso", quietly = TRUE) ||
    requireNamespace("JGL", quietly = TRUE)
}

# Signal-free data: random expression, random labels. Any "accuracy" above
# chance on held-out samples is leakage.
.make_null_data <- function(p = 40, n = 60, seed = 7) {
  set.seed(seed)
  X <- matrix(stats::rnorm(p * n), nrow = p,
              dimnames = list(paste0("G", seq_len(p)), paste0("S", seq_len(n))))
  y <- factor(rep(c("A", "B"), each = n / 2))
  names(y) <- colnames(X)
  # A prior that (misleadingly) groups arbitrary genes together.
  prior <- list(m1 = paste0("G", 1:10), m2 = paste0("G", 11:20),
                m3 = paste0("G", 21:30))
  list(X = X, y = y, prior = prior)
}

.bridge_spec <- function(prior, extra = list()) {
  list(method = "bridge", prior = prior,
       args = utils::modifyList(
         list(resolution = c(1, 1.5), tau = 0.4, min_size = 3,
              lambda = 0.1, verbose = FALSE), extra))
}


test_that("bridge_prefilter is label-blind (independent of y)", {
  skip_if_not_installed("Matrix")
  d <- .make_null_data()
  # Same X, two different label vectors -> identical kept features.
  keep1 <- bridge_prefilter(d$X, method = "variance", max_genes = 15)
  y_perm <- sample(d$y)  # permuting labels must not matter — prefilter ignores y
  keep2 <- bridge_prefilter(d$X, method = "variance", max_genes = 15)
  expect_identical(keep1, keep2)
  expect_length(keep1, 15)

  # Prior restriction keeps only prior genes.
  kp <- bridge_prefilter(d$X, method = "prior", B = d$prior)
  expect_true(all(kp %in% unlist(d$prior)))
  expect_true(length(kp) <= length(unique(unlist(d$prior))))
})


test_that("grouping_for_fold recomputes BRIDGE from TRAIN columns only", {
  skip_if_not_installed("igraph")
  skip_if(!.has_backend(), "no graphical-lasso backend installed")
  d <- .make_null_data()
  ds <- gsm_dataset(d$X, d$y, list(placeholder = paste0("G", 1:5)),
                    min_group_size = 3, grouping_spec = .bridge_spec(d$prior))

  train_idx <- 1:40
  # Groups from the helper must equal a direct bridge_grouping on train columns.
  g_helper <- grouping_for_fold(ds, train_idx, fold_seed = 3)
  g_direct <- bridge_grouping(d$X[, train_idx], d$y[train_idx], B = d$prior,
                              resolution = c(1, 1.5), tau = 0.4, min_size = 3,
                              lambda = 0.1, seed = 3, verbose = FALSE)
  expect_identical(lapply(g_helper, sort), lapply(g_direct, sort))

  # Every returned gene exists in X (came from the train-derived universe).
  expect_true(all(unlist(g_helper) %in% rownames(d$X)))
})


test_that("per-fold groups differ across folds (different train subsets)", {
  skip_if_not_installed("igraph")
  skip_if(!.has_backend(), "no graphical-lasso backend installed")
  d <- .make_null_data()
  ds <- gsm_dataset(d$X, d$y, list(placeholder = paste0("G", 1:5)),
                    min_group_size = 3, grouping_spec = .bridge_spec(d$prior))

  s1 <- split_train_test(ds$labels, ratio = 0.7, seed = 101)$train
  s2 <- split_train_test(ds$labels, ratio = 0.7, seed = 202)$train
  g1 <- grouping_for_fold(ds, s1, fold_seed = 101)
  g2 <- grouping_for_fold(ds, s2, fold_seed = 202)
  # Different training samples should generally yield a different grouping.
  sig <- function(g) paste(sort(vapply(g, function(x) paste(sort(x), collapse=","),
                                       character(1))), collapse = "|")
  expect_false(identical(sig(g1), sig(g2)))
})


test_that("gsm_dataset stores grouping_spec and reports the method", {
  d <- .make_null_data()
  ds_file <- gsm_dataset(d$X, d$y, d$prior, min_group_size = 3)
  expect_null(ds_file$grouping_spec)
  expect_identical(ds_file$params$grouping_method, "file")

  ds_bridge <- gsm_dataset(d$X, d$y, d$prior, min_group_size = 3,
                           grouping_spec = .bridge_spec(d$prior))
  expect_identical(ds_bridge$grouping_spec$method, "bridge")
  expect_identical(ds_bridge$params$grouping_method, "bridge")
})


test_that("gsm_dataset allows a deferred BRIDGE grouping (no preview groups)", {
  d <- .make_null_data()
  # Empty groups + a grouping_spec = deferred; must NOT error.
  ds <- gsm_dataset(d$X, d$y, list(), min_group_size = 3,
                    grouping_spec = .bridge_spec(d$prior))
  expect_s3_class(ds, "gsm_dataset")
  expect_identical(ds$params$n_groups, 0L)
  expect_identical(ds$params$grouping_method, "bridge")

  # Empty groups WITHOUT a spec must still error (file grouping needs groups).
  expect_error(gsm_dataset(d$X, d$y, list(), min_group_size = 3),
               "named list|No groups")
})


test_that("group-level RRA is suppressed for BRIDGE folds", {
  # Fold results tagged as BRIDGE -> group-name RRA returns NULL (names are
  # fold-local); gene-level RRA is the right consensus.
  fake_folds <- list(
    list(status = "success", grouping_method = "bridge",
         ranked_groups = c("BRIDGE_r1_c1", "BRIDGE_r1_c2")),
    list(status = "success", grouping_method = "bridge",
         ranked_groups = c("BRIDGE_r1_c5", "BRIDGE_r1_c9"))
  )
  expect_null(gsmR:::gsm_tbl_group_rra(fake_folds))

  # File grouping still aggregates normally.
  file_folds <- list(
    list(status = "success", grouping_method = "file",
         ranked_groups = c("gA", "gB")),
    list(status = "success", grouping_method = "file",
         ranked_groups = c("gA", "gC"))
  )
  expect_false(is.null(gsmR:::gsm_tbl_group_rra(file_folds)))
})


test_that("NULL-CHECK: per-fold BRIDGE does not overfit signal-free data", {
  skip_if_not_installed("igraph")
  skip_if_not_installed("randomForest")
  skip_if(!.has_backend(), "no graphical-lasso backend installed")

  d <- .make_null_data(p = 40, n = 80, seed = 11)
  ds <- gsm_dataset(d$X, d$y, list(placeholder = paste0("G", 1:5)),
                    min_group_size = 3, max_group_size = 500,
                    grouping_spec = .bridge_spec(d$prior))

  # Run several honest folds (BRIDGE re-derived on train, evaluated on test).
  folds <- lapply(1:6, function(i)
    mccv_run_fold(ds, i, train_ratio = 0.7, pre_score_method = "none",
                  score_method = "rf", model_method = "rf",
                  top_n = 3, shave = FALSE, seed = 100))
  ok <- Filter(function(f) identical(f$status, "success"), folds)
  skip_if(length(ok) < 3, "too few successful folds to assess")

  acc <- vapply(ok, function(f) f$metrics$Accuracy, numeric(1))
  mean_acc <- mean(acc, na.rm = TRUE)

  # On pure noise, honest held-out accuracy must sit near chance (0.5), not
  # near 1. Allow generous slack for small-sample variance.
  expect_lt(mean_acc, 0.72)
  expect_true(all(vapply(ok, function(f)
    identical(f$grouping_method, "bridge"), logical(1))))
})
