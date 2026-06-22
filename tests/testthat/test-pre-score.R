# Tests for pre-scoring functions

test_that("pre_score_anova ranks signal group highest", {
  # Build a dataset where GroupA has signal genes
  f_expr <- test_data_path("test_expression.csv")
  f_groups <- test_data_path("test_groups.gmt")

  expr_data <- read_expression(f_expr)
  groups <- read_gmt(f_groups)

  ds <- gsm_dataset(expr_data$expr, expr_data$labels, groups, min_group_size = 5)

  # Use first 14 samples as training
  train_idx <- 1:14

  # Pre-score with ANOVA, keep top 40% (should keep at least 2 groups)
  result <- pre_score_groups(ds, train_idx, method = "anova", top_percent = 40)

  expect_type(result, "list")
  expect_true(length(result) >= 1)
  # The signal group should be among the top
  expect_true("GroupA_signal" %in% names(result))
})


test_that("pre_score_none returns all groups", {
  expr <- matrix(rnorm(60), nrow = 6, ncol = 10)
  rownames(expr) <- paste0("G", 1:6)
  colnames(expr) <- paste0("S", 1:10)
  labels <- factor(rep(c("A", "B"), 5))
  names(labels) <- colnames(expr)
  groups <- list(
    g1 = paste0("G", 1:6),
    g2 = paste0("G", 1:6)
  )

  ds <- gsm_dataset(expr, labels, groups, min_group_size = 1)
  result <- pre_score_groups(ds, train_idx = 1:7, method = "none")

  expect_equal(length(result), length(ds$groups))
})


test_that("pre_score_groups errors on unknown method", {
  expr <- matrix(rnorm(60), nrow = 6, ncol = 10)
  rownames(expr) <- paste0("G", 1:6)
  colnames(expr) <- paste0("S", 1:10)
  labels <- factor(rep(c("A", "B"), 5))
  names(labels) <- colnames(expr)
  groups <- list(g1 = paste0("G", 1:6))

  ds <- gsm_dataset(expr, labels, groups, min_group_size = 1)
  expect_error(pre_score_groups(ds, 1:7, method = "bogus"), "Unknown")
})


test_that("pre_score_groups validates inputs", {
  expect_error(pre_score_groups("not_a_dataset", 1:5))
})


# ── New p-value-based methods ───────────────────────────────────────────────

# Build a two-class dataset with a small (10-gene) signal group and a large
# (200-gene) pure-noise group. Raw Fisher (-2*sum(log p)) is a SUM, so the big
# noise group accumulates a larger statistic than the small signal group and is
# ranked higher — the size bias the calibrated/competitive methods must fix.
make_size_bias_ds <- function(seed = 1) {
  set.seed(seed)
  n_samp <- 40
  labels <- factor(rep(c("A", "B"), each = n_samp / 2))
  n_signal <- 10; n_noise <- 200
  signal <- t(vapply(seq_len(n_signal), function(i) {
    # Moderate per-gene shift: each signal gene is detectable, but 10 of them
    # accumulate a SMALLER raw Fisher sum than 200 noise genes (~2 each) do.
    c(stats::rnorm(n_samp / 2, 0), stats::rnorm(n_samp / 2, 0.9))
  }, numeric(n_samp)))
  noise <- matrix(stats::rnorm(n_noise * n_samp), nrow = n_noise)
  expr <- rbind(signal, noise)
  rownames(expr) <- c(paste0("sig", seq_len(n_signal)),
                      paste0("noise", seq_len(n_noise)))
  colnames(expr) <- paste0("S", seq_len(n_samp))
  names(labels) <- colnames(expr)
  groups <- list(signal   = paste0("sig", seq_len(n_signal)),
                 noise_big = paste0("noise", seq_len(n_noise)))
  gsm_dataset(expr, labels, groups, min_group_size = 5, max_group_size = 500)
}

test_that("raw Fisher (anova) is biased toward the large group", {
  ds <- make_size_bias_ds()
  sc <- attr(pre_score_groups(ds, seq_len(40), method = "anova",
                              top_percent = 100), "scores")
  # The pure-noise 200-gene group outscores the 10-gene signal group.
  expect_gt(sc[["noise_big"]], sc[["signal"]])
})

test_that("fisher_chisq, ebm and camera rank the signal group above noise", {
  ds <- make_size_bias_ds()
  for (m in c("fisher_chisq", "ebm")) {
    res <- pre_score_groups(ds, seq_len(40), method = m, top_percent = 100)
    sc  <- attr(res, "scores")
    expect_gt(sc[["signal"]], sc[["noise_big"]])
    # contract: scores sorted descending, p-values + BH attached
    expect_false(is.unsorted(rev(sc)))
    expect_named(attr(res, "pvalues"))
    expect_named(attr(res, "padj"))
    expect_true(all(attr(res, "padj") >= attr(res, "pvalues") - 1e-9))
  }

  skip_if_not_installed("limma")
  res <- pre_score_groups(ds, seq_len(40), method = "camera", top_percent = 100)
  sc  <- attr(res, "scores")
  expect_gt(sc[["signal"]], sc[["noise_big"]])
})

test_that("camera falls back to fisher_chisq for >2 classes", {
  set.seed(2)
  expr <- matrix(stats::rnorm(30 * 18), nrow = 30)
  rownames(expr) <- paste0("G", seq_len(30))
  colnames(expr) <- paste0("S", seq_len(18))
  labels <- factor(rep(c("A", "B", "C"), each = 6))
  names(labels) <- colnames(expr)
  groups <- list(g1 = paste0("G", 1:10), g2 = paste0("G", 11:30))
  ds <- gsm_dataset(expr, labels, groups, min_group_size = 5)

  skip_if_not_installed("limma")
  expect_warning(
    res <- pre_score_groups(ds, seq_len(18), method = "camera", top_percent = 100),
    "fisher_chisq"
  )
  expect_type(res, "list")
  expect_true(length(res) >= 1)
})

test_that("gsva single-sample enrichment ranks the signal group above noise", {
  ds <- make_size_bias_ds()
  # Force the dependency-free z-score backend so the test never needs GSVA.
  res <- pre_score_groups(ds, seq_len(40), method = "gsva",
                          top_percent = 100, enrichment = "zscore")
  sc  <- attr(res, "scores")
  expect_gt(sc[["signal"]], sc[["noise_big"]])
  expect_named(attr(res, "padj"))
  expect_false(is.unsorted(rev(sc)))
})

test_that("gsva method runs via the dispatcher (GSVA backend when installed)", {
  ds <- make_size_bias_ds()
  res <- pre_score_groups(ds, seq_len(40), method = "gsva", top_percent = 100)
  expect_type(res, "list")
  expect_true(length(res) >= 1)
  # Whichever backend ran, the discriminative signal group should win.
  expect_gt(attr(res, "scores")[["signal"]], attr(res, "scores")[["noise_big"]])
})

test_that("EBM reduces to Fisher chi-square for a single-gene group", {
  set.seed(3)
  expr <- matrix(stats::rnorm(20 * 20), nrow = 20)
  rownames(expr) <- paste0("G", seq_len(20))
  colnames(expr) <- paste0("S", seq_len(20))
  labels <- factor(rep(c("A", "B"), 10))
  names(labels) <- colnames(expr)
  groups <- list(solo = "G1", g2 = paste0("G", 2:6))
  ds <- gsm_dataset(expr, labels, groups, min_group_size = 1)

  p_ebm <- attr(pre_score_groups(ds, seq_len(20), method = "ebm",
                                 top_percent = 100), "pvalues")[["solo"]]
  p_chi <- attr(pre_score_groups(ds, seq_len(20), method = "fisher_chisq",
                                 top_percent = 100), "pvalues")[["solo"]]
  # Single gene: both equal the gene's own ANOVA p-value.
  expect_equal(p_ebm, p_chi, tolerance = 1e-8)
})
