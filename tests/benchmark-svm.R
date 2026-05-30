#!/usr/bin/env Rscript

# Standalone Benchmark Test for G-S-M utilizing SVM Classifier in Parallel
#
# Runs:
# 1. Load GDS3257.csv (107 samples, 22,225 features).
# 2. Build GSM Dataset.
# 3. Execute 10-fold Parallel MCCV using Support Vector Machine (SVM) on 10 cores.
# 4. Report exact execution times and aggregate metrics.
#
# Run using:
#   Rscript tests/benchmark-svm.R

cat("======================================================\n")
cat("       gsmR: Standalone Parallel SVM Benchmark        \n")
cat("======================================================\n")

library(gsmR)
library(data.table)
library(e1071)

GDS_PATH <- "/Users/jallmer/Git/GSM/The-G-S-M-Grouping-Scoring-Modeling-Approach/example_input_files/GDS3257.csv"
GROUP_PATH <- "/Users/jallmer/Git/GSM/The-G-S-M-Grouping-Scoring-Modeling-Approach/example_grouping_files/cancer-DisGeNET_gedinet.csv"

if (!file.exists(GDS_PATH) || !file.exists(GROUP_PATH)) {
  stop("ERROR: Example files not found.")
}

# 1. Load Data
expr_data <- read_expression(GDS_PATH, class_col = "class")
groups <- read_groups(GROUP_PATH, gene_col = 1, group_col = 2)
dataset <- gsm_dataset(expr_data$expr, expr_data$labels, groups, min_group_size = 5, max_group_size = 500)

available_cores <- parallel::detectCores()
target_cores <- available_cores # Max power 10 cores!

cat(sprintf("\n[1/2] Preparing Parallel SVM MCCV (10 iterations on %d cores)...\n", target_cores))

start_time <- Sys.time()

# 2. Run MCCV with SVM and a linear kernel!
# We pass kernel = "linear" and cost = 1 dynamically!
mccv_results <- run_mccv(
  dataset = dataset,
  n_iterations = 10,
  train_ratio = 0.7,
  pre_score_method = "anova",
  pre_score_top_percent = 1,
  score_method = "svm",
  model_method = "svm",
  score_args = list(kernel = "linear", cost = 1),
  model_args = list(kernel = "linear", cost = 1),
  cores = target_cores,
  seed = 42
)

duration <- difftime(Sys.time(), start_time, units = "secs")

# 3. Display Aggregated Benchmark Results
cat("\n======================================================\n")
cat("                PARALLEL SVM BENCHMARK RESULTS         \n")
cat("======================================================\n")
cat(sprintf("  Iterations (Folds): %d\n", length(mccv_results$folds)))
cat(sprintf("  Cores Employed:     %d / %d available\n", target_cores, available_cores))
cat(sprintf("  Classifier:         SVM (Linear Kernel, Cost = 1)\n"))
cat(sprintf("  Total Execution:    %.2f seconds\n", as.numeric(duration)))
cat(sprintf("  Avg Time per Fold:  %.2f seconds (estimated raw)\n", as.numeric(duration) * target_cores / 10))
cat("------------------------------------------------------\n")
cat("  AGGREGATED METRICS (MEAN ± SD):\n")

sum_df <- mccv_results$metrics_summary
for (metric_name in rownames(sum_df)) {
  mean_val <- sum_df[metric_name, "Mean"]
  sd_val <- sum_df[metric_name, "SD"]
  cat(sprintf("  - %-16s : %6.2f%%  (± %5.2f%%)\n", 
              metric_name, mean_val * 100, sd_val * 100))
}
cat("======================================================\n")

# Assertions to verify correctness
mean_acc <- sum_df["Accuracy", "Mean"]
if (is.na(mean_acc) || mean_acc < 0.5) {
  stop("Verification Failed: Mean accuracy is unreasonably low or NaN.")
}

cat(sprintf("\n✓ SUCCESS: 10-fold Parallel SVM MCCV executed perfectly in %.2f seconds!\n", as.numeric(duration)))
