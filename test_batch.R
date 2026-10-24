pkgload::load_all(".")
expr_files <- list.files("example_input_files", pattern="\\.(csv|tsv|txt)$", full.names=TRUE)
expr_files <- setNames(expr_files, tools::file_path_sans_ext(basename(expr_files)))
cat("Found", length(expr_files), "expression files\n")
groups_path <- file.path("example_grouping_files", "cancer-Disgenet_gedinet.csv")
cat("Groups file exists:", file.exists(groups_path), "\n")
params <- list(
  score_method     = "rf",
  prescore_method  = "anova",
  mccv_iterations  = 2L,
  top_percent      = 20,
  mccv_train_ratio = 0.9,
  n_top_groups     = 10L,
  model_method     = "same",
  shave_genes      = TRUE,
  normalize_z      = TRUE,
  fdr_alpha        = 0.05,
  score_seed       = 44L
)
tryCatch({
  st <- gsm_batch_new_state(expr_files, groups_path, params=params, outdir=tempdir(), cores=1, max_resident=2)
  cat("Batch state created OK, datasets:", length(st$names), "\n")
}, error = function(e) {
  cat("ERROR:", conditionMessage(e), "\n")
  cat("TRACEBACK:\n")
  traceback()
})
