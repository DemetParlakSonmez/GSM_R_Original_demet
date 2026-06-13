# gsmR — Grouping-Scoring-Modeling Framework for R

An optimized R implementation of the G-S-M framework for transcriptomic data analysis with biological prior knowledge (e.g. pathways, gene sets, or disease groups). This package extends the original KNIME workflow with high-performance concurrent processing, pluggable machine learning backends, and a bidirectional GUI ↔ Terminal pipeline.

---

## Key Features

- **📂 Replicate Collapsing**: Collapse technical replicates using dynamic aggregation methods (mean, median, sum, first, etc.) with automated sample mapping and provenance tracking.
- **📊 Statistical Pre-Scoring**: Filter out irrelevant biological groups using fast statistical tests (One-Way ANOVA, limma moderated statistics) before executing expensive machine learning runs.
- **⚡ Parallel Monte Carlo Cross-Validation (MCCV)**: Run multi-core stratified CV folds concurrently utilizing multi-process socket clusters on Windows and process-forking (`mclapply`) on Unix/macOS.
- **🎯 Pluggable ML Backends**: Seamlessly train and score groups using **Random Forest (RF)**, **Support Vector Machines (SVM)**, **K-Nearest Neighbors (KNN)**, or **XGBoost** with dynamic hyperparameter forwarding via R's `...` ellipsis.
- **🌋 De Novo Module Discovery**: Discover robust gene modules via RRA-ranked gene stability networks with Louvain clustering and bootstrap validation.
- **🖥️ Interactive Shiny Dashboard**: A dark-themed, modern GUI with glassmorphism styling, maximizable panels, contextual tooltips, and inline progress indicators.
- **🔁 GUI ↔ Terminal Pipeline**: Build your analysis parameters graphically, export them as YAML, run the analysis headlessly from the terminal, and load the results back into the GUI for interactive exploration.

---

## The GUI ↔ Terminal Workflow

A core design principle of gsmR is that **parameter configuration and result inspection are visual tasks, while computation belongs on the command line**. This allows you to:

1. **Design your experiment in the GUI** — load data, configure pre-scoring, choose classifiers, set MCCV iterations, tweak discovery parameters — all with immediate visual feedback.
2. **Export the configuration** — with one click, save a self-contained YAML parameter file that captures every setting including input file paths.
3. **Run the analysis in the terminal** — on your laptop, an HPC cluster, or a remote server. The YAML file is the single source of truth; command-line flags can override any parameter.
4. **Inspect results in the GUI** — the terminal produces `.rds` session files in the exact same format as the GUI's "Save Session". Load them back into the Shiny app and every tab (Overview, Performance, Groups, Genes, Discover, Iterations, Export) is fully populated.

```
┌─────────────────────────────────────────┐
│          Shiny GUI (Browser)            │
│                                         │
│  1. Load data, configure parameters     │
│  2. File → Export Settings → params.yaml│
│            │                            │
│            ▼                            │
│  ┌──────────────────────────────┐       │
│  │  params.yaml                 │       │
│  │  expression_file: data.csv   │       │
│  │  groups_file: groups.gmt     │       │
│  │  score_method: rf            │       │
│  │  mccv_iterations: 50         │       │
│  │  ...                         │       │
│  └──────────┬───────────────────┘       │
│             │                           │
│  4. File → Load Session ← results.rds  │
│     All tabs populated instantly        │
└─────────────┼───────────────────────────┘
              │
              ▼
┌─────────────────────────────────────────┐
│          Terminal / HPC / Server         │
│                                         │
│  3. Rscript -e '                        │
│       gsmR::run_pipeline("params.yaml", │
│         outdir = "results/")'           │
│                                         │
│  Writes:                                │
│   - gsmR_results_GDS3257_20260613.rds   │
│   - gsmR_params_GDS3257_20260613.yaml   │
└─────────────────────────────────────────┘
```

This separation is especially useful when:
- You have **large datasets** that take hours to process — configure once in the GUI, then run headlessly.
- You need **reproducibility** — the timestamped YAML captures every parameter, and the `.rds` contains all intermediate results.
- You want to **batch-process multiple datasets** — write a shell loop that calls `run_pipeline()` with different `--expr` files.
- You're running on a **remote server / HPC cluster** without a browser — the CLI pipeline requires no display.

---

## Installation

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

## CLI Pipeline (`run_pipeline`)

Run the full G-S-M analysis headlessly from the terminal. The pipeline reads a YAML parameter file (exported from the Shiny GUI or written by hand), executes all stages, and saves results as an `.rds` session file that the Shiny app can open directly.

### From R
```r
library(gsmR)

# Basic: run with a YAML config
run_pipeline("params.yaml")

# Override any parameter on the fly
run_pipeline("params.yaml",
  expression_file = "new_dataset.csv",
  mccv_iterations = 50,
  mccv_cores = 8,
  outdir = "results/"
)

# Skip module discovery (faster)
run_pipeline("params.yaml", discover = FALSE)
```

### From the Terminal (CLI Script)
```bash
# Basic usage
inst/bin/gsmr params.yaml

# Override expression file and output directory
inst/bin/gsmr params.yaml --expr new_data.csv --outdir results/

# Change classifier and increase iterations
inst/bin/gsmr params.yaml --method svm --iterations 50 --cores 8

# Skip discovery step
inst/bin/gsmr params.yaml --no-discover

# Print help
inst/bin/gsmr --help
```

### CLI Output
The pipeline produces two timestamped files:

| File | Contents |
|------|----------|
| `gsmR_results_{dataset}_{timestamp}.rds` | Full session file (loadable by Shiny's *Load Session*) |
| `gsmR_params_{dataset}_{timestamp}.yaml` | Exact parameters used (including any CLI overrides) |

### Parameter Override Precedence
```
CLI flags (--method svm)  >  YAML file  >  Built-in defaults
```

---

## Interactive Shiny GUI Dashboard

gsmR contains a premium web-based interface built with modern bootstrap glassmorphism styling. Launch the application locally using:

```r
library(gsmR)
run_app()
```

To host the dashboard **always-on** (Docker, Hugging Face Spaces, Cloud Run, or shinyapps.io), see [DEPLOY.md](DEPLOY.md).

### Dashboard Features

- **Dark theme** with carefully tuned contrast for extended use
- **Maximizable panels** — click the expand icon in each panel header to view any chart or table in full screen
- **Contextual tooltips** — hover over the ⓘ icons next to every parameter for explanations
- **Inline progress indicators** — analysis status banners appear directly in the content area (no distracting modal overlays)
- **Immediate button feedback** — buttons disable on click to prevent duplicate submissions during long computations
- **Session save/load** — save your entire analysis state (data, results, settings) to an `.rds` file and restore it later
- **Settings import/export** — export all parameters as YAML for reproducibility or CLI use

### Dashboard Tabs

| Tab | Purpose |
|-----|---------|
| **📂 Data** | Upload expression matrix and grouping files. Visualizes class distributions and group sizes. |
| **🔁 Replicates** | Optionally upload metadata to collapse technical replicates with full provenance tracking. |
| **📊 Pre-Score** | Run ANOVA/limma pre-scoring to prune the pathway list and increase signal-to-noise ratio. |
| **🎯 Score** | Central command: single-fold or parallel MCCV scoring with fold boxplots and stability rankings. |
| **🌋 Discover** | De novo module discovery: robust gene stability network visualization with interactive node tooltips. |
| **🔄 Iterations** | Per-fold breakdown with convergence analysis across MCCV iterations. |
| **📤 Export** | Download results as XLSX or CSV, export settings as YAML, or save/load full session files. |

### File Menu

The **File** menu (accessible from the sidebar) provides:
- **💾 Save Session** — save all data, results, and settings to an `.rds` file at a location you choose (native file dialog on macOS)
- **📂 Load Session** — restore a previously saved session (including CLI-generated `.rds` files)
- **📤 Export Settings** — download all current parameters as a YAML file for CLI use
- **📥 Import Settings** — load parameters from a YAML file to pre-configure the GUI

---

## Developer Verification & Testing

### Test Suite
The `tests/testthat/` directory contains 350+ automated tests covering:
- Data loading and dataset construction
- All classifier backends (RF, SVM, KNN, XGBoost)
- Pre-scoring (ANOVA, limma, pass-through)
- MCCV execution and fold aggregation
- Module discovery
- Export and settings round-trip
- **CLI pipeline** — parameter parsing, YAML loading, CLI overrides, output structure, single-fold mode, MCCV mode, discovery, and Shiny round-trip compatibility

Run the full suite:
```r
# From the package root
testthat::test_dir("tests/testthat")
```

### Standalone Analysis Scripts
The `tests/` directory also contains targeted R scripts for deep analytical investigations:
- **`benchmark-svm.R`**: Times and evaluates multi-core linear SVM performance on real GDS3257 data.
- **`convergence-analysis.R`**: Iteratively aggregates results from multiple MCCV folds to determine after how many iterations the top-performing groups reach consensus stability.
- **`runner-ups-analysis.R`**: Calculates group-by-group metrics (appearance rate in Top 10, average rank, average accuracy) across all cross-validation folds.

---

## Citation & References

- Khokhar, M., Bakir-Gungor, B., & Yousef, M. (2025). *Enhancing the Efficiency of the Grouping-Scoring-Modeling Framework with Statistical Pre-Scoring Component for Transcriptomic Data Analysis.* BIOINFORMATICS 2025.
  DOI: [10.5220/0013192600003911](https://doi.org/10.5220/0013192600003911)
- [Original G-S-M Repository](https://github.com/malikyousef/The-G-S-M-Grouping-Scoring-Modeling-Approach)

