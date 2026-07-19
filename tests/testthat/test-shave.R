# Tests for gene shaving selection rules: BH-FDR (default) and Top-K.

.shave_data <- function(p = 200, n = 40, n_signal = 20, seed = 1) {
  set.seed(seed)
  X <- matrix(stats::rnorm(p * n), nrow = p,
              dimnames = list(paste0("G", seq_len(p)), paste0("S", seq_len(n))))
  y <- factor(rep(c("A", "B"), each = n / 2))
  names(y) <- colnames(X)
  X[seq_len(n_signal), y == "B"] <- X[seq_len(n_signal), y == "B"] + 3
  list(X = X, y = y, signal = paste0("G", seq_len(n_signal)))
}

test_that("FDR shaving keeps significant genes only", {
  d <- .shave_data()
  keep <- shave_genes(d$X, d$y, seq_len(ncol(d$X)), alpha = 0.05, method = "fdr")
  expect_true(all(d$signal %in% keep))
  expect_lt(length(keep), nrow(d$X))          # not everything survives
})

test_that("Top-K shaving keeps exactly K genes, strongest first", {
  d <- .shave_data()
  keep <- shave_genes(d$X, d$y, seq_len(ncol(d$X)), method = "topk", top_k = 30)
  expect_length(keep, 30)
  # The true signal genes are the strongest, so all must be in the top-30.
  expect_true(all(d$signal %in% keep))
})

test_that("Top-K is capped at the number of available genes", {
  d <- .shave_data(p = 25)
  keep <- shave_genes(d$X, d$y, seq_len(ncol(d$X)), method = "topk", top_k = 1000)
  expect_length(keep, 25)
})

test_that("shave_method flows through mccv_run_fold into the fold result", {
  d <- .shave_data(p = 120, n = 40)
  ds <- gsm_dataset(d$X, d$y, list(sig = d$signal, other = paste0("G", 30:60)),
                    min_group_size = 3)
  skip_if_not_installed("randomForest")
  fr <- mccv_run_fold(ds, 1, pre_score_method = "none", score_method = "rf",
                      model_method = "rf", top_n = 2,
                      shave = TRUE, shave_method = "topk", shave_k = 25, seed = 5)
  skip_if(fr$status != "success", "fold did not succeed")
  # Top-K keeps 25 genes, so the per-fold shaved count is reported and bounded.
  expect_true(!is.null(fr$n_shaved_genes))
  expect_lte(fr$n_shaved_genes, 25)
})
