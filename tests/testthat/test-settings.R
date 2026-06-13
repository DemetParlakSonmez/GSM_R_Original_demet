# tests/testthat/test-settings.R
# ==============================================================================
# Tests for settings import/export, including input file paths
# ==============================================================================

test_that("gsmr_settings_template includes expression_file placeholder", {
  tmpl <- gsmr_settings_template()
  expect_match(tmpl, "expression_file.*\\{expression_file\\}")
})

test_that("gsmr_settings_template includes groups_file placeholder", {
  tmpl <- gsmr_settings_template()
  expect_match(tmpl, "groups_file.*\\{groups_file\\}")
})

test_that("export_settings writes valid YAML with file paths", {
  f <- tempfile(fileext = ".yaml")
  on.exit(unlink(f))

  current_values <- list(
    expression_file = "/data/expr.csv",
    groups_file     = "/data/groups.gmt",
    groups_format   = "auto",
    min_group_size  = 5,
    max_group_size  = 500,
    collapse_method = "mean",
    meta_conflict   = "first",
    keep_original   = FALSE,
    prescore_method = "anova",
    top_percent     = 20,
    train_ratio     = 0.7,
    prescore_seed   = 42,
    perm_method     = "rf",
    perm_nperm      = 100,
    score_method    = "rf",
    rf_ntree        = 500,
    svm_kernel      = "radial",
    svm_cost        = 1,
    svm_gamma       = "auto",
    knn_k           = 5,
    shave_genes     = TRUE,
    fdr_alpha       = 0.05,
    normalize_z     = FALSE,
    score_scope     = "selected",
    n_top_groups    = 10,
    score_seed      = 44,
    run_mccv_mode   = TRUE,
    mccv_iterations = 10,
    model_method    = "same",
    mccv_train_ratio = 0.7,
    mccv_cores      = 4,
    detect_mode     = "abs",
    detect_top_genes = 200,
    detect_min_freq  = 0.5,
    detect_resolution = 1,
    detect_nboot    = 50
  )

  export_settings(current_values, f)

  expect_true(file.exists(f))
  content <- paste(readLines(f), collapse = "\n")

  # File paths present
  expect_match(content, "expression_file: '/data/expr.csv'")
  expect_match(content, "groups_file: '/data/groups.gmt'")

  # Can be read back as valid YAML
  cfg <- read_settings(f)
  expect_equal(cfg$expression_file, "/data/expr.csv")
  expect_equal(cfg$groups_file, "/data/groups.gmt")
  expect_equal(cfg$score_method, "rf")
  expect_equal(cfg$mccv_iterations, 10)
})

test_that("read_settings returns a named list", {
  f <- tempfile(fileext = ".yaml")
  on.exit(unlink(f))
  writeLines("score_method: svm\nmccv_iterations: 5\n", f)
  cfg <- read_settings(f)
  expect_true(is.list(cfg))
  expect_equal(cfg$score_method, "svm")
  expect_equal(cfg$mccv_iterations, 5)
})

test_that("export_settings round-trips all parameters", {
  f <- tempfile(fileext = ".yaml")
  on.exit(unlink(f))

  vals <- list(
    expression_file = "data.csv", groups_file = "groups.gmt",
    groups_format = "csv", min_group_size = 10, max_group_size = 200,
    collapse_method = "median", meta_conflict = "NA", keep_original = TRUE,
    prescore_method = "limma", top_percent = 30, train_ratio = 0.8,
    prescore_seed = 99, perm_method = "svm", perm_nperm = 500,
    score_method = "svm", rf_ntree = 1000, svm_kernel = "linear",
    svm_cost = 10, svm_gamma = 0.01, knn_k = 7,
    shave_genes = FALSE, fdr_alpha = 0.1, normalize_z = TRUE,
    score_scope = "all", n_top_groups = 20, score_seed = 55,
    run_mccv_mode = FALSE, mccv_iterations = 25, model_method = "rf",
    mccv_train_ratio = 0.6, mccv_cores = 8,
    detect_mode = "up", detect_top_genes = 100,
    detect_min_freq = 0.3, detect_resolution = 2, detect_nboot = 100
  )

  export_settings(vals, f)
  cfg <- read_settings(f)

  expect_equal(cfg$expression_file, "data.csv")
  expect_equal(cfg$groups_file, "groups.gmt")
  expect_equal(cfg$score_method, "svm")
  expect_equal(cfg$mccv_iterations, 25)
  expect_equal(cfg$detect_mode, "up")
  expect_equal(cfg$svm_cost, 10)
})
