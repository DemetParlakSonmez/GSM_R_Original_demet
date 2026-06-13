# tests/testthat/test-pipeline.R
# ==============================================================================
# Unit and integration tests for run_pipeline()
# ==============================================================================

# ── Helpers ───────────────────────────────────────────────────────────────────

# Create minimal expression CSV and groups CSV files for testing.
make_test_files <- function(dir) {
  # Expression: 20 samples, 30 genes, 2 classes
  set.seed(123)
  n <- 20; g <- 30
  expr <- matrix(rnorm(n * g), nrow = n, ncol = g)
  colnames(expr) <- paste0("gene", seq_len(g))
  cls <- factor(rep(c("pos", "neg"), each = n / 2))
  df <- data.frame(class = cls, expr, check.names = FALSE)
  expr_file <- file.path(dir, "test_expr.csv")
  write.csv(df, expr_file, row.names = FALSE)

  # Groups: 6 groups of 5 genes each
  grp <- data.frame(
    gene  = paste0("gene", c(1:5, 6:10, 11:15, 16:20, 21:25, 26:30)),
    group = rep(paste0("G", 1:6), each = 5),
    stringsAsFactors = FALSE
  )
  groups_file <- file.path(dir, "test_groups.csv")
  write.csv(grp, groups_file, row.names = FALSE)

  list(expr_file = expr_file, groups_file = groups_file)
}

# ── Unit Tests: Parameter Parsing ─────────────────────────────────────────────

test_that("run_pipeline rejects missing expression_file", {
  expect_error(run_pipeline(params = list(groups_file = "x.csv")),
               "expression_file.*required")
})

test_that("run_pipeline rejects missing groups_file", {
  expect_error(run_pipeline(params = list(expression_file = "x.csv")),
               "groups_file.*required")
})

test_that("run_pipeline rejects non-existent expression file", {
  expect_error(
    run_pipeline(params = list(expression_file = "no_such_file.csv",
                               groups_file = "no_such_groups.csv")),
    "not found")
})

test_that("run_pipeline accepts a YAML file path", {
  dir <- tempdir()
  files <- make_test_files(dir)
  yaml_file <- file.path(dir, "params.yaml")
  writeLines(sprintf(
    "expression_file: '%s'\ngroups_file: '%s'\nrun_mccv_mode: false\nscore_method: rf\n",
    files$expr_file, files$groups_file), yaml_file)
  expect_error(
    run_pipeline(yaml_file, outdir = file.path(dir, "out"), discover = FALSE),
    NA)  # no error
  unlink(file.path(dir, "out"), recursive = TRUE)
  unlink(yaml_file)
})

test_that("run_pipeline CLI overrides take precedence over YAML", {
  dir <- tempdir()
  files <- make_test_files(dir)
  yaml_file <- file.path(dir, "params2.yaml")
  writeLines(sprintf(
    "expression_file: '%s'\ngroups_file: '%s'\nscore_method: svm\nrun_mccv_mode: false\n",
    files$expr_file, files$groups_file), yaml_file)
  # Override score_method from svm to rf
  outdir <- file.path(dir, "out2")
  result <- run_pipeline(yaml_file, outdir = outdir, discover = FALSE,
                         score_method = "rf")
  sess <- readRDS(result)
  expect_equal(sess$settings$score_method, "rf")  # override won
  unlink(outdir, recursive = TRUE)
  unlink(yaml_file)
})

# ── Unit Tests: Output Structure ──────────────────────────────────────────────

test_that("run_pipeline output .rds has correct structure", {
  dir <- tempdir()
  files <- make_test_files(dir)
  outdir <- file.path(dir, "out3")
  result <- run_pipeline(
    params = list(expression_file = files$expr_file,
                  groups_file = files$groups_file,
                  run_mccv_mode = FALSE,
                  score_method = "rf"),
    outdir = outdir, discover = FALSE)

  expect_true(file.exists(result))
  expect_match(result, "\\.rds$")

  sess <- readRDS(result)

  # Required fields (same as Shiny session save)
  required_fields <- c("expr_data", "groups", "dataset", "dataset_name",
                        "loaded", "scored", "settings", "saved_at",
                        "r_version", "gsmR_version", "source")
  for (f in required_fields) {
    expect_true(f %in% names(sess), info = paste("Missing field:", f))
  }

  expect_true(sess$loaded)
  expect_true(sess$scored)
  expect_equal(sess$source, "cli_pipeline")
  expect_s3_class(sess$saved_at, "POSIXct")

  unlink(outdir, recursive = TRUE)
})

test_that("run_pipeline saves a params YAML alongside results", {
  dir <- tempdir()
  files <- make_test_files(dir)
  outdir <- file.path(dir, "out4")
  result <- run_pipeline(
    params = list(expression_file = files$expr_file,
                  groups_file = files$groups_file,
                  run_mccv_mode = FALSE),
    outdir = outdir, discover = FALSE)

  yaml_files <- list.files(outdir, pattern = "\\.yaml$", full.names = TRUE)
  expect_length(yaml_files, 1)
  # Read and verify it's valid YAML
  cfg <- yaml::read_yaml(yaml_files[1])
  expect_true("expression_file" %in% names(cfg))
  expect_true("score_method" %in% names(cfg))

  unlink(outdir, recursive = TRUE)
})

# ── Integration Tests: Single-Fold Mode ───────────────────────────────────────

test_that("run_pipeline single-fold produces score_results", {
  dir <- tempdir()
  files <- make_test_files(dir)
  outdir <- file.path(dir, "out5")
  result <- run_pipeline(
    params = list(expression_file = files$expr_file,
                  groups_file = files$groups_file,
                  run_mccv_mode = FALSE,
                  score_method = "rf",
                  n_top_groups = 3,
                  top_percent = 50),
    outdir = outdir, discover = FALSE)

  sess <- readRDS(result)
  expect_false(is.null(sess$score_results))
  expect_true(length(sess$score_results) > 0)
  # Each result should have accuracy
  accs <- vapply(sess$score_results, function(x) x$accuracy, numeric(1))
  expect_true(all(is.finite(accs)))
  expect_true(all(accs >= 0 & accs <= 1))

  unlink(outdir, recursive = TRUE)
})

# ── Integration Tests: MCCV Mode ─────────────────────────────────────────────

test_that("run_pipeline MCCV mode produces mccv_results", {
  dir <- tempdir()
  files <- make_test_files(dir)
  outdir <- file.path(dir, "out6")
  result <- run_pipeline(
    params = list(expression_file = files$expr_file,
                  groups_file = files$groups_file,
                  run_mccv_mode = TRUE,
                  mccv_iterations = 2,
                  mccv_cores = 1,
                  score_method = "rf",
                  n_top_groups = 3,
                  top_percent = 50),
    outdir = outdir, discover = FALSE)

  sess <- readRDS(result)
  expect_false(is.null(sess$mccv_results))
  expect_true("folds" %in% names(sess$mccv_results))
  expect_equal(length(sess$mccv_results$folds), 2)

  unlink(outdir, recursive = TRUE)
})

# ── Integration Tests: Discovery ──────────────────────────────────────────────

test_that("run_pipeline with discover=TRUE runs module discovery", {
  dir <- tempdir()
  files <- make_test_files(dir)
  outdir <- file.path(dir, "out7")
  result <- run_pipeline(
    params = list(expression_file = files$expr_file,
                  groups_file = files$groups_file,
                  run_mccv_mode = TRUE,
                  mccv_iterations = 2,
                  mccv_cores = 1,
                  score_method = "rf",
                  n_top_groups = 3,
                  top_percent = 50,
                  detect_top_genes = 20,
                  detect_nboot = 5),
    outdir = outdir, discover = TRUE)

  sess <- readRDS(result)
  # detect may be NULL if discovery found no modules, but it should have run
  # without error — the field should exist
  expect_true("detect" %in% names(sess))

  unlink(outdir, recursive = TRUE)
})

test_that("run_pipeline skips discovery gracefully in single-fold mode", {
  dir <- tempdir()
  files <- make_test_files(dir)
  outdir <- file.path(dir, "out8")
  # Single-fold + discover=TRUE → should warn and skip (discovery needs MCCV)
  expect_message(
    run_pipeline(
      params = list(expression_file = files$expr_file,
                    groups_file = files$groups_file,
                    run_mccv_mode = FALSE),
      outdir = outdir, discover = TRUE),
    "Skipping discovery")

  unlink(outdir, recursive = TRUE)
})

# ── Integration Tests: Round-Trip Compatibility ───────────────────────────────

test_that("CLI .rds can be loaded by Shiny session loader", {
  dir <- tempdir()
  files <- make_test_files(dir)
  outdir <- file.path(dir, "out9")
  result <- run_pipeline(
    params = list(expression_file = files$expr_file,
                  groups_file = files$groups_file,
                  run_mccv_mode = FALSE),
    outdir = outdir, discover = FALSE)

  sess <- readRDS(result)

  # Simulate the Shiny load_session validation
  expect_true(is.list(sess))
  expect_false(is.null(sess$saved_at))
  expect_true(inherits(sess$dataset, "gsm_dataset"))
  expect_true(is.logical(sess$loaded) && sess$loaded)
  expect_true(is.logical(sess$scored) && sess$scored)

  # Settings should restore Shiny controls
  s <- sess$settings
  expect_true(is.list(s))
  expect_true("score_method" %in% names(s))
  expect_true("run_mccv_mode" %in% names(s))
  expect_true("n_top_groups" %in% names(s))

  unlink(outdir, recursive = TRUE)
})
