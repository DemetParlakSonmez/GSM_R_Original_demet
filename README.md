# gsmR — Grouping-Scoring-Modeling Framework for Transcriptomic Data

An R implementation of the G-S-M framework for transcriptomic data analysis,
migrated from the original KNIME workflow.

## Installation

```r
# Install from local source
devtools::install("path/to/GSM4R")

# Or from GitHub
# remotes::install_github("jallmer/GSM4R")
```

## Quick Start

```r
library(gsmR)

# 1. Read data
expr_data <- read_expression("GDS3257.csv")
groups    <- read_groups("cancer-DisGeNET_gedinet.csv")

# 2. Construct validated dataset
ds <- gsm_dataset(
  expr   = expr_data$expr,
  labels = expr_data$labels,
  groups = groups
)

# 3. Pre-score groups (optional — select top 20% by ANOVA)
top_groups <- pre_score_groups(ds, train_idx = 1:80, method = "anova", top_percent = 20)

# 4. Score individual groups with Random Forest
score <- score_group(ds$expr, ds$labels, genes = ds$groups[["Glioblastoma"]],
                     train_idx = 1:80, method = "rf")
```

## Architecture

G-S-M uses a strategy pattern for pluggable backends:

| Step | Function | Default | Alternatives |
|------|----------|---------|--------------|
| Pre-Scoring | `pre_score_groups()` | ANOVA | limma (future), none |
| Scoring | `score_group()` | Random Forest | SVM (future), etc. |
| Modeling | `build_model()` | Random Forest | glmnet (future), etc. |

## Shiny App

gsmR includes an interactive web interface. Launch it with:

```r
library(gsmR)
run_app()
```

**Required packages:** `shiny`, `bslib`, `DT` (install with `install.packages(c("shiny", "bslib", "DT"))`).

The app has four tabs:

| Tab | Description |
|-----|-------------|
| **📂 Data** | Upload expression & grouping files, explore class distribution and group sizes |
| **🔁 Replicates** | Upload metadata, collapse technical replicates with provenance tracking |
| **📊 Pre-Score** | Run ANOVA pre-scoring to filter groups before ML scoring |
| **🎯 Score** | Score top groups with Random Forest and compare accuracies |

## References

- Khokhar, M., Bakir-Gungor, B., & Yousef, M. (2025). *Enhancing the Efficiency of the Grouping-Scoring-Modeling Framework with Statistical Pre-Scoring Component for Transcriptomic Data Analysis.* BIOINFORMATICS 2025.
  DOI: [10.5220/0013192600003911](https://doi.org/10.5220/0013192600003911)
- [Original G-S-M Repository](https://github.com/malikyousef/The-G-S-M-Grouping-Scoring-Modeling-Approach)
