#!/usr/bin/env Rscript

# Standalone Integration Test for gsmR using real published G-S-M example files.
#
# Runs:
# 1. Load raw expression matrix from GDS3257.csv (107 samples, 22,225 features).
# 2. Load gene disease groupings from cancer-DisGeNET_gedinet.csv (~3,900 groups).
# 3. Construct and validate gsm_dataset.
# 4. Split into training and testing partitions (stratified split).
# 5. Run ANOVA pre-scoring to select top groups.
# 6. Score the top groups using Random Forest classifier.
# 7. Build and evaluate the final integrated Random Forest classifier.
#
# Run specifically using:
#   Rscript tests/integration-test.js (or .R)

cat("======================================================\n")
cat("     gsmR: Standalone Real Data Integration Test      \n")
cat("======================================================\n")

# Load gsmR dependencies and source files directly in the workspace
cat("\n[0/7] Sourcing package in workspace...\n")
library(data.table)
library(randomForest)
source("R/input.R")
source("R/dataset.R")
source("R/utils.R")
source("R/pre_score.R")
source("R/score.R")
source("R/model.R")
source("R/metrics.R")

GDS_PATH <- "/Users/jallmer/Git/GSM/The-G-S-M-Grouping-Scoring-Modeling-Approach/example_input_files/GDS3257.csv"
GROUP_PATH <- "/Users/jallmer/Git/GSM/The-G-S-M-Grouping-Scoring-Modeling-Approach/example_grouping_files/cancer-DisGeNET_gedinet.csv"

if (!file.exists(GDS_PATH) || !file.exists(GROUP_PATH)) {
  stop("ERROR: Example GDS3257.csv or cancer-DisGeNET_gedinet.csv files not found. Ensure repository is cloned.")
}

start_time <- Sys.time()

# 1. Load Expression CSV
cat("\n[1/7] Parsing GDS3257.csv expression matrix...\n")
expr_data <- read_expression(GDS_PATH, class_col = "class")
cat(sprintf("      Loaded %d genes across %d samples.\n", nrow(expr_data$expr), ncol(expr_data$expr)))

# 2. Load Grouping CSV
cat("\n[2/7] Parsing cancer-DisGeNET_gedinet.csv groupings...\n")
groups <- read_groups(GROUP_PATH, gene_col = 1, group_col = 2)
cat(sprintf("      Loaded %d disease groups.\n", length(groups)))

# 3. Build GSM Dataset
cat("\n[3/7] Building and validating gsm_dataset object...\n")
dataset <- gsm_dataset(expr_data$expr, expr_data$labels, groups, min_group_size = 5, max_group_size = 500)
cat("      Successfully validated.\n")
cat(sprintf("      Intersecting groups remaining: %d\n", dataset$params$n_groups))

# 4. Split Train/Test
cat("\n[4/7] Performing reproducible stratified split (70% Train)...\n")
train_ratio <- 0.7
seed <- 42

split <- split_train_test(dataset$labels, ratio = train_ratio, seed = seed)
train_idx <- split$train
test_idx <- split$test

cat(sprintf("      Training samples: %d, Testing samples: %d\n", length(train_idx), length(test_idx)))

# 5. Pre-scoring
cat("\n[5/7] Running statistical pre-scoring (ANOVA)...\n")
# Select top 1% of groups to keep the test extremely fast while still validating accuracy
top_groups <- pre_score_groups(dataset, train_idx, method = "anova", top_percent = 1)
n_top <- length(top_groups)
cat(sprintf("      Selected top %d groups.\n", n_top))

# 6. Group Scoring
cat("\n[6/7] Scoring selected top groups with Random Forest...\n")
scored_results <- list()
for (name in names(top_groups)) {
  genes <- top_groups[[name]]
  score <- score_group(dataset$expr, dataset$labels, genes, train_idx, method = "rf")
  scored_results[[name]] <- score
  cat(sprintf("      Group scored: %s | Accuracy: %.1f%% | Genes: %d\n",
              name, score$accuracy * 100, length(genes)))
}

# 7. Final Model Building
cat("\n[7/7] Building integrated final Random Forest model on top groups...\n")
final_model <- build_model(dataset$expr, dataset$labels, top_groups, train_idx, method = "rf")

cat("\n======================================================\n")
cat("              INTEGRATION TEST RESULTS                \n")
cat("======================================================\n")
cat(sprintf("  Features Combined:  %d genes\n", length(final_model$features)))
cat(sprintf("  Accuracy:           %.2f%%\n", final_model$metrics$Accuracy * 100))
cat(sprintf("  Sensitivity:        %.2f%%\n", final_model$metrics$Sensitivity * 100))
cat(sprintf("  Specificity:        %.2f%%\n", final_model$metrics$Specificity * 100))
cat(sprintf("  F-Measure (F1):     %.2f%%\n", final_model$metrics$F_measure * 100))
cat(sprintf("  Cohen's Kappa:      %.4f\n", final_model$metrics$Cohens_kappa))
cat("======================================================\n")

# Verification assertions
if (length(final_model$features) == 0) {
  stop("Verification Failed: Final model has 0 features.")
}
if (is.na(final_model$metrics$Accuracy) || final_model$metrics$Accuracy < 0.5) {
  stop("Verification Failed: Model accuracy is unreasonably low or NaN.")
}

duration <- difftime(Sys.time(), start_time, units = "secs")
cat(sprintf("\n✓ SUCCESS: Integration test completed perfectly in %.2f seconds!\n", as.numeric(duration)))
