#!/usr/bin/env Rscript

# Standalone Benchmark and Integration Test for Parallel MCCV in gsmR
#
# Runs:
# 1. Load real published expression matrix (GDS3257.csv, 107 samples, 22,225 features).
# 2. Load groupings (cancer-DisGeNET_gedinet.csv, ~3,900 groups).
# 3. Construct and validate gsm_dataset.
# 4. Execute 10-fold Parallel MCCV on all available cores (minus one).
# 5. Report exact execution times and aggregated metrics (Mean ± SD).
#
# Run specifically using:
#   Rscript tests/benchmark-mccv.R

cat("======================================================\n")
cat("    gsmR: Standalone Parallel MCCV Benchmark Test     \n")
cat("======================================================\n")

# Load gsmR and parallel
library(gsmR)
library(data.table)
library(randomForest)

GDS_PATH <- "/Users/jallmer/Git/GSM/The-G-S-M-Grouping-Scoring-Modeling-Approach/example_input_files/GDS3257.csv"
GROUP_PATH <- "/Users/jallmer/Git/GSM/The-G-S-M-Grouping-Scoring-Modeling-Approach/example_grouping_files/cancer-DisGeNET_gedinet.csv"

if (!file.exists(GDS_PATH) || !file.exists(GROUP_PATH)) {
  stop("ERROR: Example GDS3257.csv or cancer-DisGeNET_gedinet.csv files not found. Ensure repository is cloned.")
}

# 1. Load Expression CSV
cat("\n[1/5] Parsing GDS3257.csv expression matrix...\n")
expr_data <- read_expression(GDS_PATH, class_col = "class")
cat(sprintf("      Loaded %d genes across %d samples (preserves raw headers).\n", nrow(expr_data$expr), ncol(expr_data$expr)))

# 2. Load Grouping CSV
cat("\n[2/5] Parsing cancer-DisGeNET_gedinet.csv groupings...\n")
groups <- read_groups(GROUP_PATH, gene_col = 1, group_col = 2)
cat(sprintf("      Loaded %d disease groups.\n", length(groups)))

# 3. Build GSM Dataset
cat("\n[3/5] Building and validating gsm_dataset...\n")
dataset <- gsm_dataset(expr_data$expr, expr_data$labels, groups, min_group_size = 5, max_group_size = 500)
cat(sprintf("      Intersecting groups remaining: %d\n", dataset$params$n_groups))

# Determine optimal cores (using detectCores() - 1 as standard, but customizable)
available_cores <- parallel::detectCores()
target_cores <- max(1, available_cores - 1)

cat(sprintf("\n[4/5] Preparing Parallel MCCV (10 iterations on %d cores)...\n", target_cores))

start_time <- Sys.time()

# 4. Run MCCV
mccv_results <- run_mccv(
  dataset = dataset,
  n_iterations = 10,
  train_ratio = 0.7,
  pre_score_method = "anova",
  pre_score_top_percent = 1, # keep top 1% for extremely fast runs
  score_method = "rf",
  model_method = "rf",
  cores = target_cores,
  seed = 42
)

duration <- difftime(Sys.time(), start_time, units = "secs")

# 5. Display Aggregated Benchmark Results
cat("\n======================================================\n")
cat("             PARALLEL MCCV BENCHMARK RESULTS          \n")
cat("======================================================\n")
cat(sprintf("  Iterations (Folds): %d\n", length(mccv_results$folds)))
cat(sprintf("  Cores Employed:     %d / %d available\n", target_cores, available_cores))
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

cat(sprintf("\n✓ SUCCESS: 10-fold Parallel MCCV executed perfectly in %.2f seconds!\n", as.numeric(duration)))
