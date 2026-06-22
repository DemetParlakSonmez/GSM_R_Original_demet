#' @title Settings Management for gsmR
#' @description Import and Export configuration settings for the G-S-M framework.
#' @name settings
NULL

#' Generate YAML Settings Template
#'
#' @return A character string containing the YAML template with comments.
#' @export
gsmr_settings_template <- function() {
  "
# ==============================================================================
# gsmR Settings Configuration
# ==============================================================================
# This file stores parameters for the G-S-M framework so they can be easily
# shared, version-controlled, or reused across different datasets.
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Input Files
# ------------------------------------------------------------------------------

# Path to the expression CSV file (samples as rows, first column = class label).
expression_file: '{expression_file}'

# Path to the grouping file (GMT, CSV, or other supported format).
groups_file: '{groups_file}'

# ------------------------------------------------------------------------------
# 1. Data Input Parameters
# ------------------------------------------------------------------------------

# Format of the grouping file (e.g. 'auto', 'csv', 'gmt').
groups_format: '{groups_format}'

# Minimum number of genes required for a biological group to be considered.
min_group_size: {min_group_size}

# Maximum number of genes allowed for a biological group.
max_group_size: {max_group_size}

# ------------------------------------------------------------------------------
# 2. Replicate Collapsing Settings
# ------------------------------------------------------------------------------

# Method used to collapse multiple replicates into a single sample (e.g. 'mean', 'median', 'sum', 'min', 'max', 'first').
collapse_method: '{collapse_method}'

# How to handle metadata conflicts during collapsing (e.g. 'error', 'NA', 'concat', 'first').
meta_conflict: '{meta_conflict}'

# Whether to keep original uncollapsed samples in the dataset (true or false).
keep_original: {keep_original}

# ------------------------------------------------------------------------------
# 3. Pre-Scoring Parameters
# ------------------------------------------------------------------------------

# Method used to statistically pre-score and filter groups before ML ('anova', 'limma', 'none').
prescore_method: '{prescore_method}'

# Top percentage of groups to retain after pre-scoring (e.g., 20 for top 20%).
top_percent: {top_percent}

# Ratio of samples assigned to training split during pre-scoring (e.g., 0.7).
train_ratio: {train_ratio}

# Random seed for pre-scoring.
prescore_seed: {prescore_seed}

# Classifier used for empirical group significance permutation test.
perm_method: '{perm_method}'

# Number of permutations per group for significance testing.
perm_nperm: {perm_nperm}

# ------------------------------------------------------------------------------
# 4. Scoring & Modeling Parameters
# ------------------------------------------------------------------------------

# Classifier used to score individual groups ('rf', 'svm', 'knn', 'xgb').
score_method: '{score_method}'

# Number of trees for Random Forest (if score_method or model_method is 'rf').
rf_ntree: {rf_ntree}

# SVM Kernel type ('linear', 'radial', 'polynomial', 'sigmoid').
svm_kernel: '{svm_kernel}'

# SVM Cost parameter (C).
svm_cost: {svm_cost}

# SVM Gamma parameter (for non-linear kernels).
svm_gamma: {svm_gamma}

# K-Nearest Neighbors 'k' parameter.
knn_k: {knn_k}

# Enable Phase I gene shaving using Welch t-test and BH FDR filter (true or false).
shave_genes: {shave_genes}

# FDR threshold (alpha) for gene shaving.
fdr_alpha: {fdr_alpha}

# Z-score normalize features (true or false).
normalize_z: {normalize_z}

# Which groups to score ('selected' or 'all').
score_scope: '{score_scope}'

# Accumulative depth: Number of top-K groups to cumulatively model.
n_top_groups: {n_top_groups}

# Random seed for scoring split.
score_seed: {score_seed}

# ------------------------------------------------------------------------------
# 5. Parallel MCCV Parameters
# ------------------------------------------------------------------------------

# Whether to run in Multi-Fold Parallel MCCV mode (true or false).
run_mccv_mode: {run_mccv_mode}

# Number of cross-validation iterations for MCCV.
mccv_iterations: {mccv_iterations}

# Classifier used for building the final joint model on the top K groups ('same', 'rf', 'svm', 'knn', 'xgb').
model_method: '{model_method}'

# Ratio of samples assigned to the training split per iteration (e.g., 0.7).
mccv_train_ratio: {mccv_train_ratio}

# Number of CPU cores to use for MCCV (defaults to 1 if not set).
mccv_cores: {mccv_cores}

# ------------------------------------------------------------------------------
# 6. Discovery Parameters
# ------------------------------------------------------------------------------

# Direction for module discovery ('abs', 'up', 'down').
detect_mode: '{detect_mode}'

# Number of top genes to consider per iteration for discovery.
detect_top_genes: {detect_top_genes}

# Minimum appearance frequency across iterations.
detect_min_freq: {detect_min_freq}

# Resolution for Louvain clustering.
detect_resolution: {detect_resolution}

# Number of bootstraps for robustness testing.
detect_nboot: {detect_nboot}

# ------------------------------------------------------------------------------
# 7. Multi-Dataset Batch (optional; used by run_batch())
# ------------------------------------------------------------------------------
# When set, run_batch() processes EVERY listed expression file against the single
# groups_file above, sharing all parameters in this file. All other settings
# (pre-scoring, scoring, MCCV, discovery) apply to each dataset.

# List of expression files to process as a batch (uncomment and edit):
# expression_files:
#   - 'data/GDS1.csv'
#   - 'data/GDS2.csv'

# Output directory root for per-dataset results + the combined report.
batch_outdir: 'gsmR_batch'

# Total worker cores for the shared (dataset x fold) pool. Blank = detectCores()-1.
batch_cores: {mccv_cores}

# Maximum datasets held in memory at once (bounds peak RAM; folds still saturate cores).
max_resident: 3
"
}

#' Read settings from a YAML file
#'
#' @param file Path to the YAML file.
#' @return A named list of settings.
#' @export
read_settings <- function(file) {
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Package 'yaml' is required to read settings.")
  }
  yaml::read_yaml(file)
}

#' Export current settings to a YAML file
#'
#' @param current_values A named list of current settings values.
#' @param file Path to write the output YAML file.
#' @export
export_settings <- function(current_values, file) {
  template <- gsmr_settings_template()
  
  # Replace placeholders with actual values
  for (name in names(current_values)) {
    val <- current_values[[name]]
    if (is.null(val)) {
      val_str <- ""
    } else if (is.logical(val)) {
      val_str <- ifelse(val, "true", "false")
    } else if (is.numeric(val)) {
      val_str <- as.character(val)
    } else {
      val_str <- as.character(val)
    }
    
    # Simple search-and-replace for {name}
    pattern <- paste0("\\{", name, "\\}")
    template <- gsub(pattern, val_str, template)
  }
  
  # Write to file
  writeLines(trimws(template), file)
}
