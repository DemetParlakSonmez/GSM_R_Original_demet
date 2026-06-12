# gsmR — Grouping-Scoring-Modeling Framework for R

An optimized R implementation of the G-S-M framework for transcriptomic data analysis with biological prior knowledge (e.g. pathways, gene sets, or disease groups). This package extends the original KNIME workflow with high-performance concurrent processing, pluggable machine learning backends, and technical replicate aggregation.

---

## Key Features

- **📂 Replicate Collapsing**: Collapse technical replicates using dynamic aggregation methods (mean, median, sum, first, etc.) with automated sample mapping and provenance tracking.
- **📊 Statistical Pre-Scoring**: Filter out irrelevant biological groups using fast statistical tests (One-Way ANOVA) before executing expensive machine learning runs.
- **⚡ Parallel Monte Carlo Cross-Validation (MCCV)**: Run multi-core stratified CV folds concurrently utilizing multi-process socket clusters on Windows and process-forking (`mclapply`) on Unix/macOS.
- **🎯 Pluggable ML Backends**: Seamlessly train and score groups using **Random Forest (RF)**, **Support Vector Machines (SVM)** (via `e1071`), or **K-Nearest Neighbors (KNN)** (via `class`) with dynamic hyperparameter forwarding via R's `...` ellipsis.
- **🖥️ Upgraded Shiny GUI Dashboard**: An interactive, dark-themed dashboard to load data, collapse replicates, pre-score pathways, score individual groups, or launch parallel MCCV runs with visual fold boxplots and stability rankings.

---

## Installation

To install and compile the R package locally along with all required machine learning and visual dependencies, execute the following commands in your R console:

```r
# 1. Install devtools or remotes if not already present
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# 2. Install package dependencies
install.packages(c("shiny", "bslib", "DT", "data.table", "randomForest", "e1071"))

# 3. Install gsmR from local source
remotes::install_local("path/to/GSM4R", dependencies = TRUE, force = TRUE)
```

---

## Quick Start (R Console)

### 1. Load Data & Build Dataset
```r
library(gsmR)

# Read expression CSV (first column = class/labels, subsequent columns = genes/samples)
expr_data <- read_expression("example_input_files/GDS3257.csv", class_col = "class")

# Read biological groupings CSV (columns: gene, group)
groups <- read_groups("example_grouping_files/cancer-DisGeNET_gedinet.csv")

# Create validated G-S-M dataset object
ds <- gsm_dataset(
  expr   = expr_data$expr,
  labels = expr_data$labels,
  groups = groups,
  min_group_size = 5,
  max_group_size = 500
)
```

### 2. Pre-Score and Rank Groups
Filter out less significant groups using statistical pre-scoring (retaining the top 20%):
```r
train_idx <- 1:80
top_groups <- pre_score_groups(ds, train_idx = train_idx, method = "anova", top_percent = 20)
```

### 3. Score a Single Group (with Pluggable Classifiers)
```r
# A. Score using Random Forest (Default)
score_rf <- score_group(ds$expr, ds$labels, genes = ds$groups[["Glioblastoma"]], 
                        train_idx = train_idx, method = "rf", ntree = 100)

# B. Score using SVM (Linear Kernel, Cost = 5)
score_svm <- score_group(ds$expr, ds$labels, genes = ds$groups[["Glioblastoma"]], 
                         train_idx = train_idx, method = "svm", kernel = "linear", cost = 5)

# C. Score using K-Nearest Neighbors (k = 5)
score_knn <- score_group(ds$expr, ds$labels, genes = ds$groups[["Glioblastoma"]], 
                         train_idx = train_idx, method = "knn", k = 5)
```

### 4. Run High-Performance Parallel MCCV
Execute 10-fold Monte Carlo cross-validation concurrently on 10 cores using **SVM** for scoring and **Random Forest** for final joint model building:
```r
mccv_results <- run_mccv(
  dataset               = ds,
  n_iterations          = 10,
  train_ratio           = 0.7,
  pre_score_method      = "anova",
  pre_score_top_percent = 20,
  score_method          = "svm",
  model_method          = "rf",
  score_args            = list(kernel = "linear", cost = 1),
  model_args            = list(ntree = 500),
  cores                 = 10,
  seed                  = 42
)

# Access aggregate metrics across all folds
print(mccv_results$metrics_summary)
```

### 5. Empirical Significance of Groups (Permutation Test)
Ask *"is this group better than a random gene set of the same size?"* — the question
reviewers actually want answered. Reports an empirical p-value and BH-adjusted FDR:
```r
sig <- permutation_test_groups(
  ds,
  method     = "rf",
  n_perm     = 200,        # random same-size gene sets per group
  score_args = list(ntree = 100),
  seed       = 42
)
head(sig)   # Group, size, observed, null_mean, p_value, p_adj, n_perm
```

### 6. Consensus Biomarker (Gene) Ranking
Collapse the per-fold, per-gene model importances from an MCCV run into one
stability-weighted biomarker list (selection frequency × mean importance, plus an
RRA consensus). Works with RF, linear-kernel SVM, and XGBoost models:
```r
biomarkers <- gsm_biomarkers(mccv_results, top_n = 50)
head(biomarkers)   # Gene, Frequency, MeanImportance, RRA_Score, Rank
```

### 7. Additional Backends
- **XGBoost** gradient-boosted trees for scoring/modeling (`score_method`/`model_method = "xgb"`),
  returning per-gene `Gain` importance and probabilities for AUC.
- **limma** moderated-t / moderated-F pre-scoring (`pre_score_method = "limma"`),
  more powerful than ANOVA at the small sample sizes typical of transcriptomics.
- Probabilities are now wired through every classifier, so **AUC-ROC** is populated
  for `score_group()` and `build_model()` (not just `NA`).

---

## Interactive Shiny GUI Dashboard

gsmR contains a premium web-based interface built with modern bootstrap glassmorphism styling. Launch the application locally using:

```r
library(gsmR)
run_app()
```

To host the dashboard **always-on** (Docker, Hugging Face Spaces, Cloud Run, or shinyapps.io), see [DEPLOY.md](DEPLOY.md).

### Dashboard Tabs Overview

1. **📂 Data Tab**: Upload your expression matrix and prior biological knowledge files. Instantaneously explores class distributions and group sizes with clean charts.
2. **🔁 Replicates Tab**: Optionally upload metadata to resolve technical replicates. Automatically collapses duplicate entries via your chosen aggregation statistic with complete provenance tables.
3. **📊 Pre-Score Tab**: Run training-stratified ANOVA tests to prune the pathway list, increasing signal-to-noise ratio.
4. **🎯 Score Tab**: The central command center:
   - Toggle **Single Fold Mode** to rank top pathways and build a consolidated joint model using **RF**, **SVM**, or **KNN**.
   - Toggle **Parallel MCCV Mode** to execute multi-core CV runs in the background. Generates aggregate summary sheets, fold-performance boxplots, and ranked group stability list (runner-ups analysis).

---

## Developer Verification & Standalone Scripts

The `tests/` directory contains targeted R scripts to perform deep analytical investigations:
- **`benchmark-svm.R`**: Times and evaluates multi-core linear SVM performance on real GDS3257 data.
- **`convergence-analysis.R`**: Iteratively aggregates results from multiple MCCV folds to determine after how many iterations the top-performing groups reach consensus stability.
- **`runner-ups-analysis.R`**: Calculates group-by-group metrics (appearance rate in Top 10, average rank, average accuracy) across all cross-validation folds.

---

## Citation & References

- Khokhar, M., Bakir-Gungor, B., & Yousef, M. (2025). *Enhancing the Efficiency of the Grouping-Scoring-Modeling Framework with Statistical Pre-Scoring Component for Transcriptomic Data Analysis.* BIOINFORMATICS 2025.
  DOI: [10.5220/0013192600003911](https://doi.org/10.5220/0013192600003911)
- [Original G-S-M Repository](https://github.com/malikyousef/The-G-S-M-Grouping-Scoring-Modeling-Approach)
