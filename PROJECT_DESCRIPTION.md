# Project Description — gsmR & BRIDGE

## Summary

**gsmR** is an R package that implements the **Grouping–Scoring–Modeling (G-S-M)**
framework for analyzing high-dimensional biological data (e.g. gene-expression /
transcriptomics) with **prior biological knowledge**. Instead of treating each gene
independently, G-S-M works at the level of *biologically meaningful groups of genes*
(pathways, gene-ontology terms, disease-associated sets), which improves
interpretability, statistical power, and robustness.

The package modernizes the original KNIME G-S-M workflow into a high-performance R
implementation with parallel cross-validation, pluggable machine-learning backends,
a command-line pipeline, and an interactive Shiny dashboard.

**BRIDGE** (Biologically-Regularized Interaction-Differential Grouping Engine) is the
package's flagship **data-driven Grouping engine**. Where the classic framework reads a
*predefined* grouping from a file, BRIDGE **derives the groups directly from the data**,
regularized by prior biological knowledge — discovering condition-specific gene modules
that a fixed pathway database would miss.

---

## The G-S-M Framework

The framework runs in three coupled stages:

1. **Grouping (G)** — Partition genes into groups. Either from a knowledge source
   (KEGG, GO, DisGeNET, or any GMT / two-column file) **or** data-drivenly with BRIDGE.
2. **Scoring (S)** — Score how well each group discriminates the classes, using fast
   statistical pre-scoring (ANOVA, limma, Fisher, GSVA/ssGSEA, …) followed by a
   machine-learning classifier.
3. **Modeling (M)** — Build predictive models from the top-ranked groups and evaluate
   them under **Monte-Carlo cross-validation (MCCV)**.

Six classifiers are interchangeable in **both** the Scoring and Modeling steps:

| Method | Key | Backend | Per-gene importance |
|--------|-----|---------|---------------------|
| Random Forest | `rf` | `randomForest` | MeanDecreaseGini |
| Support Vector Machine | `svm` | `e1071` | \|w\| (linear kernel) |
| K-Nearest Neighbors | `knn` | `class` | — |
| XGBoost | `xgb` | `xgboost` | Gain |
| Neural Network (MLP) | `mlp` | `torch` | input-gradient saliency |
| **Graph Neural Network (GCN)** | `gnn` | `torch` | node saliency |

The **GNN** is the one that is not a drop-in commodity classifier. Each *sample* is a
graph: nodes are the group's genes, node features are that sample's expression values,
and the adjacency is **BRIDGE's differential network** — so the model convolves over
exactly the interaction structure that defined the groups, closing the loop between the
Grouping and Modeling steps. With a file-based grouping (KEGG/GO/DisGeNET) there is no
such network, so it falls back to a **training-only** correlation graph, keeping the
backend usable everywhere and leakage-free.

> **Caveat, stated plainly.** Transcriptomic studies typically have tens to a few hundred
> samples. The `mlp` and `gnn` backends have far more capacity than that supports; they
> are heavily regularized (dropout + weight decay + class-balanced loss) but should be
> treated as *comparators*, with Random Forest and XGBoost as the baselines to beat —
> not as automatic upgrades.

A core design principle is a **GUI ↔ Terminal split**: configure the experiment visually
in the Shiny app, export a self-contained YAML, run the heavy computation headlessly
(laptop, server, or HPC), then load the resulting `.rds` session back into the GUI for
inspection.

---

## What BRIDGE Does

Given an expression matrix, class labels, and an optional prior graph, BRIDGE:

1. **Pre-filters genes** (optional, *label-blind* so it is leakage-safe) to keep the
   estimation tractable on transcriptome-scale data.
2. **Estimates two class-conditional sparse networks** with a **knowledge-weighted
   graphical lasso** — edges known to be biologically related are penalized less.
3. **Builds a differential network** (`Δ = Θ1 − Θ0`) and runs **multi-resolution
   community detection** to propose candidate gene modules at multiple scales.
4. **Selects modules** by a Grouping–Scoring coupled objective
   `J(M) = D(M) + λ1·C(M) − λ2·|M|`, keeping only modules that pass a
   **bootstrap-stability** threshold and a size filter.

The output is a named list of gene groups, **identical in shape** to the file-based
readers, so Scoring and Modeling consume it unchanged. Inside cross-validation, BRIDGE
groups are **re-derived per training fold** so no test-set information leaks into grouping.

---

## Scientific Contribution

- **Data-driven, knowledge-regularized grouping.** BRIDGE bridges the gap between purely
  data-driven module discovery (which ignores biology) and fixed pathway databases
  (which ignore the data at hand).
- **Condition-specific interaction differences.** By contrasting class-conditional
  networks, it surfaces gene modules whose *interactions* change between conditions —
  signals invisible to single-gene differential-expression analysis.
- **Leakage-free evaluation.** Grouping is recomputed inside each CV fold on training
  data only, giving honest performance estimates.

---

## Inputs Required

| Input | Description |
|-------|-------------|
| **Expression matrix `X`** | Numeric, **genes × samples**, with gene IDs as row names. |
| **Labels `y`** | Factor/vector of length `ncol(X)`; **exactly two classes**. |
| **Prior graph `B`** *(optional)* | Named grouping list, 3-column edge list (`gene_i, gene_j, weight`), or symmetric weight matrix. Omit for pure data-driven grouping. |

**Environment:** R (≥ 4.1) with `Matrix` and a graphical-lasso backend (`glasso`;
optionally `JGL`, `huge`, `igraph`). In this dev environment R lives in the conda env
`gsmr` (`export PATH=/opt/conda/envs/gsmr/bin:$PATH`).

---

## How to Run & Collect Results

**1. Direct R call**
```r
library(gsmR)
groups <- bridge_grouping(X, y, B = prior_groups, details = TRUE)
attr(groups, "scores")       # J(M) per selected module
attr(groups, "diagnostics")  # per-module J / D / C / stability / size
```

**2. Full G-S-M pipeline (headless / HPC)**
```r
gsmR::run_pipeline("params.yaml", outdir = "results/")
# → timestamped .rds session + .yaml params; BRIDGE re-derived per fold.
```

**3. Shiny app** — the **Grouping (G) → BRIDGE** tab: choose a prior, tune
hyperparameters, run with a live progress bar, inspect modules/genes/diagnostics, then
apply the grouping to the scoring-and-modeling pipeline. Results load into every
dashboard tab (Overview, Performance, Groups, Genes, Iterations, Export).

**Outputs:** selected gene groups with per-module quality scores, plus (through the full
pipeline) cross-validated classification performance, ranked biomarker groups/genes, and
a reproducible session file.

---

## Key Features

- 📂 Replicate collapsing with provenance tracking
- 📊 Statistical pre-scoring to prune irrelevant groups before ML
- 🌉 **BRIDGE** data-driven, knowledge-regularized grouping
- ⚡ Parallel Monte-Carlo cross-validation (fork / socket clusters)
- 🎯 Pluggable ML backends (RF, SVM, KNN, XGBoost, MLP, Graph Neural Network)
- 🧮 Multi-dataset batch processing with combined reports
- 🖥️ Interactive Shiny dashboard + 🔁 GUI ↔ Terminal reproducible pipeline
