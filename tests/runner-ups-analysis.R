#!/usr/bin/env Rscript

# Standalone Runner-Ups and Stability Analysis for G-S-M
#
# Runs:
# 1. Load real published GDS3257 dataset.
# 2. Execute 10-fold Parallel MCCV.
# 3. For each fold, extract all scored groups and rank them by accuracy descending.
# 4. Across all 10 folds, calculate:
#    - How many times each group appeared in the Top 10.
#    - The average accuracy of each group.
#    - The average rank of each group (1 = best).
# 5. Display a beautiful, sorted summary table of the Top 10 stable groups.
#
# Run using:
#   Rscript tests/runner-ups-analysis.R

cat("======================================================\n")
cat("      gsmR: G-S-M Runner-Up & Stability Analysis      \n")
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

# 2. Run 10-fold Parallel MCCV on 10 cores
cat("\n[1/2] Executing 10-fold Parallel MCCV on 10 cores...\n")
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

# 3. Aggregate all group performances across folds
cat("\n[2/2] Aggregating stability, ranks, and accuracies across all 10 folds...\n")

# We will collect accuracy and rank for each group in each fold
group_performance <- list()

for (i in 1:10) {
  fold_res <- mccv_results$folds[[i]]
  if (fold_res$status != "success") next
  
  # group_accuracies is a named list of accuracies
  accs <- unlist(fold_res$group_accuracies)
  
  # Rank the groups (1 = highest accuracy, handles ties by average rank)
  ranks <- rank(-accs, ties.method = "min")
  
  for (gname in names(accs)) {
    if (is.null(group_performance[[gname]])) {
      group_performance[[gname]] <- list(
        accuracies = numeric(0),
        ranks = numeric(0),
        in_top_10 = 0
      )
    }
    
    group_performance[[gname]]$accuracies <- c(group_performance[[gname]]$accuracies, accs[gname])
    group_performance[[gname]]$ranks <- c(group_performance[[gname]]$ranks, ranks[gname])
    
    if (ranks[gname] <= 10) {
      group_performance[[gname]]$in_top_10 <- group_performance[[gname]]$in_top_10 + 1
    }
  }
}

# 4. Compute summaries
summary_list <- list()
for (gname in names(group_performance)) {
  perf <- group_performance[[gname]]
  
  # Folds in which this group was pre-scored and evaluated (since ANOVA selections can vary slightly)
  n_evals <- length(perf$accuracies)
  
  summary_list[[gname]] <- data.frame(
    Group = gname,
    Folds_Evaluated = n_evals,
    Top_10_Appearances = perf$in_top_10,
    Appearance_Rate = (perf$in_top_10 / 10) * 100,
    Mean_Accuracy = mean(perf$accuracies) * 100,
    Mean_Rank = mean(perf$ranks),
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, summary_list)

# Sort by:
# 1. Top 10 Appearances (descending)
# 2. Mean Accuracy (descending)
# 3. Mean Rank (ascending)
sorted_summary <- summary_df[order(-summary_df$Top_10_Appearances, -summary_df$Mean_Accuracy, summary_df$Mean_Rank), ]

# 5. Display beautiful ranked summary
cat("\n=========================================================================================\n")
cat("                       TOP 10 STABLE G-S-M GROUPS (RUNNER-UPS)                           \n")
cat("=========================================================================================\n")
cat(sprintf("  %-38s | %-13s | %-13s | %-10s\n", "Disease Group", "Top 10 Folds", "Mean Acc (%)", "Mean Rank"))
cat("-----------------------------------------------------------------------------------------\n")

# Display top 10
top_n_to_show <- min(10, nrow(sorted_summary))
for (k in 1:top_n_to_show) {
  row <- sorted_summary[k, ]
  cat(sprintf("  %2d. %-34s |    %2d / 10     |    %6.2f%%   |   %4.1f\n",
              k,
              substring(row$Group, 1, 34),
              row$Top_10_Appearances,
              row$Mean_Accuracy,
              row$Mean_Rank))
}
cat("=========================================================================================\n")

# Also list other groups that were scored but didn't make the top 10 stable list
if (nrow(sorted_summary) > 10) {
  cat("\n  Additional Groups Scored (Runner-Ups below Rank 10):\n")
  for (k in 11:min(20, nrow(sorted_summary))) {
    row <- sorted_summary[k, ]
    cat(sprintf("  - %-38s : Top 10 in %d/10 folds (Mean Acc: %.2f%%, Mean Rank: %.1f)\n",
                row$Group, row$Top_10_Appearances, row$Mean_Accuracy, row$Mean_Rank))
  }
}

cat("\n✓ SUCCESS: Runner-ups stability analysis completed perfectly!\n")
