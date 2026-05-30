#!/usr/bin/env Rscript

# Standalone Convergence and Stability Analysis for G-S-M
#
# Runs:
# 1. Load real published expression matrix and groupings.
# 2. Execute 10-fold Parallel MCCV.
# 3. For each fold, identify the highest-accuracy disease group.
# 4. Iteratively evaluate consensus from 1 fold up to 10 folds.
# 5. Display convergence table showing which group is on top and its majority vote.
#
# Run using:
#   Rscript tests/convergence-analysis.R

cat("======================================================\n")
cat("      gsmR: MCCV Convergence & Stability Analysis      \n")
cat("======================================================\n")

library(gsmR)
library(data.table)
library(randomForest)

GDS_PATH <- "/Users/jallmer/Git/GSM/The-G-S-M-Grouping-Scoring-Modeling-Approach/example_input_files/GDS3257.csv"
GROUP_PATH <- "/Users/jallmer/Git/GSM/The-G-S-M-Grouping-Scoring-Modeling-Approach/example_grouping_files/cancer-DisGeNET_gedinet.csv"

if (!file.exists(GDS_PATH) || !file.exists(GROUP_PATH)) {
  stop("ERROR: Example files not found.")
}

# 1. Load Data
expr_data <- read_expression(GDS_PATH, class_col = "class")
groups <- read_groups(GROUP_PATH, gene_col = 1, group_col = 2)
dataset <- gsm_dataset(expr_data$expr, expr_data$labels, groups, min_group_size = 5, max_group_size = 500)

# 2. Run 10-fold Parallel MCCV on all 10 cores
cat("\n[1/3] Executing 10-fold Parallel MCCV on 10 cores...\n")
start_time <- Sys.time()
mccv_results <- run_mccv(
  dataset = dataset,
  n_iterations = 10,
  train_ratio = 0.7,
  pre_score_method = "anova",
  pre_score_top_percent = 1,
  score_method = "rf",
  model_method = "rf",
  cores = 10,
  seed = 42
)
duration <- difftime(Sys.time(), start_time, units = "secs")
cat(sprintf("      Parallel execution finished in %.2f seconds.\n", as.numeric(duration)))

# 3. Identify top group for each fold
cat("\n[2/3] Extracting the top-scoring group for each fold...\n")
top_groups_per_fold <- character(10)

for (i in 1:10) {
  fold_res <- mccv_results$folds[[i]]
  if (fold_res$status != "success") next
  
  # group_accuracies is a named list of accuracies
  accs <- unlist(fold_res$group_accuracies)
  
  # Find the group with maximum accuracy. In case of ties, take the first one.
  max_idx <- which(accs == max(accs))[1]
  top_group <- names(accs)[max_idx]
  top_groups_per_fold[i] <- top_group
  
  cat(sprintf("      Fold %2d: Top Group = '%s' (Accuracy: %.1f%%)\n", 
              i, top_group, accs[max_idx] * 100))
}

# 4. Iteratively evaluate consensus stability
cat("\n[3/3] Analyzing convergence progression fold-by-fold...\n")
cat("=========================================================================================\n")
cat(sprintf("  %-9s | %-40s | %-12s | %-12s\n", "Folds 1:N", "Consensus Top Group (Majority)", "Vote Count", "Percentage"))
cat("=========================================================================================\n")

previous_consensus <- ""
stable_count <- 0

for (n in 1:10) {
  # Get top groups for first n folds
  sub_groups <- top_groups_per_fold[1:n]
  
  # Compute frequencies
  freqs <- table(sub_groups)
  sorted_freqs <- sort(freqs, decreasing = TRUE)
  
  consensus_group <- names(sorted_freqs)[1]
  vote_count <- as.integer(sorted_freqs[1])
  percentage <- (vote_count / n) * 100
  
  # Track stability (does it stay the same as previous?)
  if (n > 1) {
    if (consensus_group == previous_consensus) {
      stable_count <- stable_count + 1
    } else {
      stable_count <- 0
    }
  }
  
  previous_consensus <- consensus_group
  
  cat(sprintf("  Folds 1:%-2d | %-40s | %2d / %-8d | %6.1f%%\n", 
              n, 
              substring(consensus_group, 1, 40), 
              vote_count, 
              n, 
              percentage))
}
cat("=========================================================================================\n")

cat("\n✓ SUCCESS: Convergence analysis completed perfectly!\n")
