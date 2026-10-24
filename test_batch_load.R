pkgload::load_all(".", quiet = TRUE)

expr_files <- list.files("example_input_files", pattern = "\\.(csv|tsv|txt)$", full.names = TRUE)
expr_files <- setNames(expr_files, tools::file_path_sans_ext(basename(expr_files)))
cat("Toplam dataset:", length(expr_files), "\n\n")

groups_path <- file.path("example_grouping_files", "cancer-Disgenet_gedinet.csv")
cat("Grouping dosyası var mı:", file.exists(groups_path), "\n\n")

# Her dataset için hangi sorunlar var, test et
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

p <- gsm_resolve_params(params)
groups <- read_groups(groups_path)
cat("Toplam grup sayısı (grouping dosyasında):", length(groups), "\n\n")

for (nm in names(expr_files)) {
  path <- expr_files[[nm]]
  tryCatch({
    ed <- read_expression(path)
    ds <- gsm_dataset(ed$expr, ed$labels, groups,
                      min_group_size = p$min_group_size,
                      max_group_size = p$max_group_size)
    cat(sprintf("OK  | %-25s | samples=%d, groups=%d\n",
                nm, ds$params$n_samples, ds$params$n_groups))
  }, error = function(e) {
    cat(sprintf("ERR | %-25s | %s\n", nm, conditionMessage(e)))
  })
}
