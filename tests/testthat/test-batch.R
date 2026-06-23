# Tests for the multi-dataset batch runner

# Write a small two-class expression CSV (samples as rows, first col = class).
write_demo_csv <- function(path, n_per_class = 12, n_genes = 30, seed = 1,
                           signal_genes = 1:6, shift = 1.5) {
  set.seed(seed)
  n <- 2 * n_per_class
  labels <- rep(c("A", "B"), each = n_per_class)
  m <- matrix(stats::rnorm(n * n_genes), nrow = n,
              dimnames = list(paste0("S", seq_len(n)), paste0("G", seq_len(n_genes))))
  m[labels == "B", signal_genes] <- m[labels == "B", signal_genes] + shift
  df <- data.frame(class = labels, m, check.names = FALSE)
  utils::write.csv(df, path, row.names = FALSE)
  path
}

# A tiny grouping file (2-column CSV; read_groups expects gene, group order).
write_demo_groups <- function(path) {
  rows <- rbind(
    data.frame(gene = paste0("G", 1:6),   group = "sig"),
    data.frame(gene = paste0("G", 4:10),  group = "mix"),
    data.frame(gene = paste0("G", 20:27), group = "noise")
  )
  utils::write.csv(rows, path, row.names = FALSE)
  path
}

demo_batch_inputs <- function() {
  dir <- file.path(tempdir(), paste0("gsm_batch_", as.integer(stats::runif(1, 1, 1e6))))
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  list(
    dir = dir,
    files = c(dsA = write_demo_csv(file.path(dir, "dsA.csv"), seed = 11),
              dsB = write_demo_csv(file.path(dir, "dsB.csv"), seed = 22, shift = 2)),
    groups = write_demo_groups(file.path(dir, "groups.csv"))
  )
}

batch_params <- list(score_method = "rf", mccv_iterations = 3, mccv_train_ratio = 0.7,
                     prescore_method = "anova", top_percent = 100, n_top_groups = 3,
                     rf_ntree = 50)


test_that("run_batch produces per-dataset .rds with the pipeline session shape", {
  inp <- demo_batch_inputs()
  out <- file.path(inp$dir, "out")
  mf <- run_batch(inp$files, inp$groups, params = batch_params, outdir = out,
                  cores = 2, max_resident = 1, discover = FALSE, report_format = "csv")

  expect_s3_class(mf, "data.frame")
  expect_equal(nrow(mf), 2)
  expect_true(all(mf$Status == "success"))

  rds <- list.files(file.path(out, "per_dataset", "dsA"), pattern = "[.]rds$",
                    full.names = TRUE)
  expect_length(rds, 1)
  sess <- readRDS(rds)
  # Same key set a run_pipeline session has → opens in the Shiny results tabs.
  for (k in c("expr_data", "dataset", "dataset_name", "mccv_results", "scored",
              "score_summary", "settings", "source"))
    expect_true(k %in% names(sess))
  expect_identical(sess$source, "batch")
  expect_true(sess$scored)
  expect_true(!is.null(sess$mccv_results$metrics_summary))
})

test_that("run_batch writes a combined report + manifest", {
  inp <- demo_batch_inputs()
  out <- file.path(inp$dir, "out")
  run_batch(inp$files, inp$groups, params = batch_params, outdir = out,
            cores = 2, max_resident = 2, discover = FALSE, report_format = "csv")

  expect_true(file.exists(file.path(out, "batch_manifest.csv")))
  expect_true(file.exists(file.path(out, "combined", "combined_tables.zip")))
  for (nm in c("dsA", "dsB")) {
    dd <- file.path(out, "per_dataset", nm)
    expect_true(file.exists(file.path(dd, paste0(nm, ".log"))))
    expect_true(file.exists(file.path(dd, paste0(nm, "_tables.zip"))))
  }
})

test_that("batch fold results match a direct run_mccv (scheduling is deterministic)", {
  inp <- demo_batch_inputs()
  out <- file.path(inp$dir, "out_repro")
  run_batch(inp$files["dsA"], inp$groups, params = batch_params, outdir = out,
            cores = 2, max_resident = 1, discover = FALSE, report_format = "csv")
  sess <- readRDS(list.files(file.path(out, "per_dataset", "dsA"), "[.]rds$",
                             full.names = TRUE))
  batch_summary <- sess$mccv_results$metrics_summary

  # Same dataset, same seed, sequential single core.
  e <- read_expression(inp$files[["dsA"]])
  g <- read_groups(inp$groups)
  ds <- gsm_dataset(e$expr, e$labels, g, min_group_size = 5)
  direct <- run_mccv(ds, n_iterations = batch_params$mccv_iterations,
                     train_ratio = batch_params$mccv_train_ratio,
                     pre_score_method = "anova", pre_score_top_percent = 100,
                     score_method = "rf", model_method = "rf",
                     score_args = list(ntree = 50), model_args = list(ntree = 50),
                     top_n = 3, cores = 1, seed = 44)

  expect_equal(batch_summary$Mean, direct$metrics_summary$Mean, tolerance = 1e-8)
})

test_that("scheduler never exceeds max_resident datasets in memory", {
  inp <- demo_batch_inputs()
  # 2 datasets, window of 1 → at most one resident at any tick.
  st <- gsmR:::gsm_batch_new_state(
    inp$files, inp$groups, params = batch_params,
    outdir = file.path(inp$dir, "out_win"), cores = 2, max_resident = 1,
    discover = FALSE)
  max_seen <- 0L
  t0 <- Sys.time()
  while (isTRUE(st$active) &&
         as.numeric(difftime(Sys.time(), t0, units = "secs")) < 120) {
    gsmR:::gsm_batch_step(st)
    max_seen <- max(max_seen, length(st$window), length(st$loaded))
    Sys.sleep(0.02)   # yield so forked folds can finish
  }
  expect_lte(max_seen, 1L)
  expect_true(st$all_done)
})
