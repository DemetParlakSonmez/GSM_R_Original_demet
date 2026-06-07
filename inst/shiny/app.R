# ──────────────────────────────────────────────────────────────────────────────
# gsmR Shiny App — Interactive G-S-M Workflow
# ──────────────────────────────────────────────────────────────────────────────

library(shiny)
library(bslib)
library(DT)

# Load gsmR — either installed or via pkgload
if (requireNamespace("gsmR", quietly = TRUE)) {
  library(gsmR)
} else {
  # Dev mode: load from package root
  app_dir <- normalizePath(file.path(dirname(sys.frame(1)$ofile %||% "."), "..", ".."))
  if (file.exists(file.path(app_dir, "DESCRIPTION"))) {
    pkgload::load_all(app_dir)
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# Discover example datasets (repo-level example_*_files/ + bundled toy data)
# ══════════════════════════════════════════════════════════════════════════════
.launch_dir <- getOption("gsmR.launch_dir", getwd())
.example_input_dir  <- getOption("gsmR.example_input_dir",
                                 file.path(.launch_dir, "example_input_files"))
.example_groups_dir <- getOption("gsmR.example_groups_dir",
                                 file.path(.launch_dir, "example_grouping_files"))

.list_examples <- function(dir, pattern) {
  if (!dir.exists(dir)) return(character(0))
  sort(list.files(dir, pattern = pattern, ignore.case = TRUE))
}
.example_inputs <- .list_examples(.example_input_dir,  "\\.(csv|tsv|txt)$")
.example_groups <- .list_examples(.example_groups_dir, "\\.(csv|tsv|gmt|txt)$")

# Choice lists always include the bundled toy data as the first option.
.expr_choices   <- c("Bundled toy (test_expression.csv)" = "__toy__",
                     stats::setNames(.example_inputs, .example_inputs))
.groups_choices <- c("Bundled toy (test_groups.csv)" = "__toy__",
                     stats::setNames(.example_groups, .example_groups))

# ══════════════════════════════════════════════════════════════════════════════
# Helper: value box HTML
# ══════════════════════════════════════════════════════════════════════════════
value_box_ui <- function(value, label, icon = NULL) {
  tags$div(
    class = "value-box",
    if (!is.null(icon)) tags$span(style = "float:right; font-size:1.5rem; opacity:0.4;", icon),
    tags$div(class = "value", value),
    tags$div(class = "label", label)
  )
}

# Null-coalescing helper (base R only ships %||% from 4.4.0).
`%||%` <- function(a, b) if (is.null(a)) b else a

# Friendly classifier name + a "key = value" string of its hyperparameters.
classifier_label <- function(method) {
  switch(method, rf = "Random Forest", svm = "Support Vector Machine",
         knn = "K-Nearest Neighbors", method)
}
format_params <- function(args) {
  if (length(args) == 0) return("defaults")
  paste(sprintf("%s = %s", names(args), unlist(args)), collapse = ", ")
}

# Human-readable elapsed time, reported in minutes or hours (seconds for short runs).
format_duration <- function(secs) {
  if (is.null(secs) || is.na(secs)) return("—")
  if (secs < 60)     return(sprintf("%.1f sec", secs))
  mins <- secs / 60
  if (mins < 60)     return(sprintf("%.2f min (%.0f sec)", mins, secs))
  hrs <- mins / 60
  sprintf("%.2f hours (%.1f min)", hrs, mins)
}

# Full-width banner reporting how long a step took.
time_banner_ui <- function(secs, label = "Completed in") {
  tags$div(
    class = "value-box",
    style = "background: rgba(74,144,217,0.12); border-left: 4px solid #4a90d9;",
    tags$div(style = "font-size:1.4rem; font-weight:600;",
             paste0("⏱️ ", format_duration(secs))),
    tags$div(class = "label", label)
  )
}

# Numbered workflow tab title: "(n) <icon> <label>".
nav_step_title <- function(n, icon, label) {
  tags$span(
    class = "nav-step-title",
    tags$span(class = "nav-step", n),
    tags$span(class = "nav-ico", icon),
    tags$span(class = "nav-lbl", label)
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# UI
# ══════════════════════════════════════════════════════════════════════════════
ui <- page_navbar(
  title = tags$span(
    tags$strong("gsmR"),
    tags$span(style = "font-weight:300; opacity:0.6; margin-left:8px;",
              "Grouping \u2022 Scoring \u2022 Modeling")
  ),
  theme = bs_theme(
    version = 5,
    bootswatch = "darkly",
    primary = "#4a90d9",
    "font-size-base" = "0.92rem",
    "border-radius" = "0.75rem"
  ),
  header = tags$head(
    tags$link(rel = "stylesheet", href = "styles.css"),
    tags$style(HTML("
      /* ── Professional fixed top bar ───────────────────────────────────── */
      .navbar {
        position: sticky; top: 0; z-index: 1040;
        background: linear-gradient(135deg, #161a2e 0%, #11152a 100%) !important;
        border-bottom: 1px solid rgba(255,255,255,0.08);
        box-shadow: 0 2px 18px rgba(0,0,0,0.45);
        padding: 0.35rem 0.5rem;
        min-height: 60px;
      }
      .navbar-brand { color: #ffffff !important; font-size: 1.15rem; }
      .navbar-nav { gap: 4px; }

      /* Tabs: clear, tactile, always legible — pill with an active accent. */
      .navbar-nav .nav-link {
        color: rgba(255,255,255,0.72) !important;
        font-weight: 500;
        padding: 0.5rem 0.9rem !important;
        border-radius: 8px;
        display: flex; align-items: center;
        transition: background .15s ease, color .15s ease;
      }
      .navbar-nav .nav-link:hover {
        color: #ffffff !important;
        background: rgba(255,255,255,0.07);
      }
      .navbar-nav .nav-link.active {
        color: #ffffff !important;
        background: rgba(74,144,217,0.20);
        box-shadow: inset 0 -3px 0 #4a90d9;
      }

      /* Step badges (workflow order: 1 Data -> 4 Score). */
      .nav-step {
        display: inline-flex; align-items: center; justify-content: center;
        width: 22px; height: 22px; border-radius: 50%;
        background: rgba(255,255,255,0.14); color: #fff;
        font-size: 0.72rem; font-weight: 700; margin-right: 8px;
      }
      .nav-link.active .nav-step { background: #4a90d9; }
      .nav-ico { margin-right: 6px; font-size: 0.95rem; }
      /* Step connector arrows between tabs. */
      .nav-step-sep { color: rgba(255,255,255,0.25); align-self: center; padding: 0 2px; }
    "))
  ),

  # ── Tab 1: Data Upload ───────────────────────────────────────────────────
  nav_panel(
    value = "data",
    title = nav_step_title("1", "\U0001F4C2", "Data"),
    layout_sidebar(
      sidebar = sidebar(
        title = "Upload Files",
        width = 340,
        tags$div(
          class = "upload-area",
          fileInput("expr_file", "\U0001F9EC Expression Data (CSV)",
                    accept = c(".csv", ".tsv", ".txt")),
          tags$small(class = "text-muted",
                     "Samples as rows, first column = class labels")
        ),
        tags$hr(),
        tags$div(
          class = "upload-area",
          fileInput("groups_file", "\U0001F5C2\uFE0F Grouping File",
                    accept = c(".csv", ".tsv", ".gmt", ".txt")),
          tags$small(class = "text-muted",
                     "2-column CSV (gene, group) or GMT format")
        ),
        tags$hr(),
        selectInput("groups_format", "Grouping format",
                    choices = c("Auto-detect" = "auto",
                                "2-column CSV" = "csv",
                                "GMT" = "gmt")),
        numericInput("min_group_size", "Min group size", value = 5, min = 1),
        numericInput("max_group_size", "Max group size", value = 500, min = 10),
        tags$hr(),
        actionButton("btn_load", "\U0001F680 Load & Validate",
                     class = "btn-primary btn-lg w-100"),
        tags$div(class = "text-center text-muted my-2",
                 style = "font-size:0.8rem;", "— or pick an example —"),
        selectInput("example_expr", "Example expression dataset",
                    choices = .expr_choices),
        selectInput("example_groups", "Example grouping file",
                    choices = .groups_choices),
        actionButton("btn_example", "\U0001F4E6 Load Selected Example",
                     class = "btn-outline-light btn-sm w-100"),
        tags$small(class = "text-muted d-block mt-1 text-center",
                   "Large GDS files may take a few seconds to load")
      ),

      # Main panel
      conditionalPanel(
        condition = "output.data_loaded",
        layout_columns(
          col_widths = c(3, 3, 3, 3),
          uiOutput("vbox_samples"),
          uiOutput("vbox_genes"),
          uiOutput("vbox_groups"),
          uiOutput("vbox_classes")
        ),
        layout_columns(
          col_widths = c(6, 6),
          card(
            card_header("Class Distribution"),
            plotOutput("plot_class_dist", height = "260px")
          ),
          card(
            card_header("Group Size Distribution"),
            plotOutput("plot_group_sizes", height = "260px")
          )
        ),
        card(
          card_header("Expression Matrix Preview"),
          DTOutput("expr_preview")
        )
      ),
      conditionalPanel(
        condition = "!output.data_loaded",
        tags$div(
          style = "text-align:center; padding:80px 20px; opacity:0.5;",
          tags$h3("\U0001F4C1 Upload expression data and grouping file"),
          tags$p("Then click 'Load & Validate' to begin.")
        )
      )
    )
  ),

  # ── Tab 2: Replicate Collapsing ──────────────────────────────────────────
  nav_panel(
    value = "replicates",
    title = nav_step_title("2", "\U0001F501", "Replicates"),
    layout_sidebar(
      sidebar = sidebar(
        title = "Collapsing Settings",
        width = 340,
        fileInput("meta_file", "\U0001F4CB Metadata (CSV)",
                  accept = c(".csv", ".tsv", ".txt")),
        tags$small(class = "text-muted",
                   "Must include sample_id and replicate_group_id columns"),
        tags$hr(),
        textInput("sample_id_col", "Sample ID column", value = "sample_id"),
        textInput("replicate_group_col", "Replicate group column",
                  value = "replicate_group_id"),
        selectInput("collapse_method", "Aggregation method",
                    choices = c("mean", "median", "sum", "min", "max", "first")),
        selectInput("meta_conflict", "Metadata conflict handling",
                    choices = c("error", "NA", "concat", "first")),
        checkboxInput("keep_original", "Keep original samples", value = FALSE),
        tags$hr(),
        actionButton("btn_collapse", "\U0001F501 Collapse Replicates",
                     class = "btn-primary btn-lg w-100")
      ),
      conditionalPanel(
        condition = "output.collapse_done",
        uiOutput("time_collapse"),
        layout_columns(
          col_widths = c(4, 4, 4),
          uiOutput("vbox_collapsed_before"),
          uiOutput("vbox_collapsed_after"),
          uiOutput("vbox_collapsed_groups")
        ),
        card(
          card_header("Provenance"),
          DTOutput("provenance_table")
        ),
        card(
          card_header("Warnings"),
          verbatimTextOutput("collapse_warnings")
        )
      ),
      conditionalPanel(
        condition = "!output.collapse_done",
        tags$div(
          style = "text-align:center; padding:80px 20px; opacity:0.5;",
          tags$h3("\U0001F501 Optional: collapse replicate samples"),
          tags$p("Upload metadata with replicate grouping to combine replicates before analysis.")
        )
      )
    )
  ),

  # ── Tab 3: Pre-Scoring ──────────────────────────────────────────────────
  nav_panel(
    value = "prescore",
    title = nav_step_title("3", "\U0001F4CA", "Pre-Score"),
    layout_sidebar(
      sidebar = sidebar(
        title = "Pre-Scoring Settings",
        width = 340,
        selectInput("prescore_method", "Method",
                    choices = c("ANOVA" = "anova", "None" = "none")),
        sliderInput("top_percent", "Top % groups to keep",
                    min = 5, max = 100, value = 20, step = 5),
        sliderInput("train_ratio", "Training ratio",
                    min = 0.5, max = 0.9, value = 0.7, step = 0.05),
        numericInput("prescore_seed", "Random seed", value = 44),
        tags$hr(),
        checkboxInput("prescore_detail",
                      "Show detailed info (per-group scores + charts)",
                      value = TRUE),
        conditionalPanel(
          condition = "input.prescore_detail",
          numericInput("prescore_topn_plot",
                       "Groups to show in bar chart", value = 20, min = 5, max = 100)
        ),
        tags$hr(),
        actionButton("btn_prescore", "\U0001F4CA Run Pre-Scoring",
                     class = "btn-primary btn-lg w-100")
      ),
      # ── Live view while pre-scoring is running ─────────────────────────
      conditionalPanel(
        condition = "output.prescore_running",
        uiOutput("prescore_live_status"),
        layout_columns(
          col_widths = c(7, 5),
          card(
            card_header("Live: top groups so far"),
            plotOutput("plot_prescore_live", height = "340px")
          ),
          card(
            card_header("Live leaderboard (top 15 so far)"),
            tableOutput("prescore_live_table")
          )
        )
      ),
      conditionalPanel(
        condition = "output.prescore_done && !output.prescore_running",
        uiOutput("time_prescore"),
        layout_columns(
          col_widths = c(4, 4, 4),
          uiOutput("vbox_prescore_input"),
          uiOutput("vbox_prescore_output"),
          uiOutput("vbox_prescore_pct")
        ),
        # Detailed charts (only when the user opts in and ANOVA produced scores)
        conditionalPanel(
          condition = "input.prescore_detail && output.prescore_has_scores",
          layout_columns(
            col_widths = c(7, 5),
            card(
              card_header("Top Groups by Pre-Score (Fisher statistic)"),
              plotOutput("plot_prescore_bars", height = "320px")
            ),
            card(
              card_header("Score Distribution (all groups)"),
              plotOutput("plot_prescore_hist", height = "320px")
            )
          )
        ),
        card(
          card_header(textOutput("prescore_table_header", inline = TRUE)),
          DTOutput("prescore_table")
        )
      ),
      conditionalPanel(
        condition = "!output.prescore_done && !output.prescore_running",
        tags$div(
          style = "text-align:center; padding:80px 20px; opacity:0.5;",
          tags$h3("\U0001F4CA Pre-score groups"),
          tags$p("Load data first, then run pre-scoring to filter groups before ML scoring.")
        )
      )
    )
  ),

  # ── Tab 4: Scoring ─────────────────────────────────────────────────────
  nav_panel(
    value = "score",
    title = nav_step_title("4", "\U0001F3AF", "Score"),
    layout_sidebar(
      sidebar = sidebar(
        title = "Scoring Settings",
        width = 340,
        selectInput("score_method", "Classifier",
                    choices = c("Random Forest" = "rf",
                                "Support Vector Machine" = "svm",
                                "K-Nearest Neighbors" = "knn")),
        
        # Hyperparameters for Random Forest
        conditionalPanel(
          condition = "input.score_method == 'rf'",
          numericInput("rf_ntree", "Number of Trees (ntree)", value = 500, min = 1, step = 50)
        ),
        
        # Hyperparameters for SVM
        conditionalPanel(
          condition = "input.score_method == 'svm'",
          selectInput("svm_kernel", "SVM Kernel",
                      choices = c("Linear" = "linear",
                                  "Radial" = "radial",
                                  "Polynomial" = "polynomial",
                                  "Sigmoid" = "sigmoid"),
                      selected = "linear"),
          numericInput("svm_cost", "Cost (C)", value = 1, min = 0.001, step = 0.5),
          numericInput("svm_gamma", "Gamma (non-linear only)", value = 0.1, min = 0.0001, step = 0.01)
        ),
        
        # Hyperparameters for KNN
        conditionalPanel(
          condition = "input.score_method == 'knn'",
          numericInput("knn_k", "Neighbors (k)", value = 3, min = 1, max = 20, step = 1)
        ),
        
        # Phase I gene shaving: training-only Welch t-test + BH-FDR filter.
        checkboxInput("shave_genes",
                      "Gene shaving (Welch t-test + BH FDR, training only)",
                      value = TRUE),
        conditionalPanel(
          condition = "input.shave_genes",
          numericInput("fdr_alpha", "FDR threshold (alpha)",
                       value = 0.05, min = 0.001, max = 0.5, step = 0.01)
        ),
        checkboxInput("normalize_z",
                      "Z-score normalize (per-fold, train-fit)", value = FALSE),

        radioButtons("score_scope", "Which groups to score",
                     choices = c("Top pre-scored groups (recommended)" = "selected",
                                 "All groups (skip pre-score filter)" = "all"),
                     selected = "selected"),
        # Accumulative depth K: build models on top-1, top-1..2, ..., top-1..K
        # ranked groups and report performance at each level (both modes).
        numericInput("n_top_groups",
                     "Accumulative depth (top-K groups to model)",
                     value = 10, min = 1),
        numericInput("score_seed", "Random seed", value = 44),

        tags$hr(),
        checkboxInput("run_mccv_mode", "Run Multi-Fold Parallel MCCV?", value = FALSE),
        
        # Parallel MCCV options
        conditionalPanel(
          condition = "input.run_mccv_mode",
          numericInput("mccv_iterations", "MCCV Iterations", value = 10, min = 1, max = 500),
          sliderInput("mccv_cores", "CPU Cores", min = 1, max = parallel::detectCores(), 
                      value = max(1, parallel::detectCores() - 1), step = 1),
          selectInput("model_method", "Final Modeling Classifier",
                      choices = c("Same as Scoring" = "same",
                                  "Random Forest" = "rf",
                                  "Support Vector Machine" = "svm",
                                  "K-Nearest Neighbors" = "knn"),
                      selected = "same"),
          sliderInput("mccv_train_ratio", "Training Ratio", min = 0.5, max = 0.9, value = 0.7, step = 0.05)
        ),
        
        tags$hr(),
        actionButton("btn_score", "\U0001F680 Run Analysis",
                     class = "btn-primary btn-lg w-100")
      ),

      # ── Always-visible status alert (readiness / what's happening) ──────
      uiOutput("score_status_alert"),

      # ── Live view while Parallel MCCV is running ───────────────────────
      conditionalPanel(
        condition = "output.mccv_running",
        uiOutput("mccv_live_status"),
        layout_columns(
          col_widths = c(6, 6),
          card(
            card_header("Live: fold accuracy as folds complete"),
            plotOutput("plot_mccv_live", height = "340px")
          ),
          card(
            card_header("Live: running metric averages"),
            uiOutput("mccv_live_means")
          )
        ),
        card(
          card_header("Live: completed folds"),
          tableOutput("mccv_live_table")
        )
      ),

      # Main Panel Results Conditional Displays
      conditionalPanel(
        condition = "output.score_done && !output.mccv_running",
        uiOutput("time_score"),
        uiOutput("score_summary"),

        # ── Export results (Excel workbook + CSV bundle) ──────────────────
        card(
          card_header("\U0001F4BE Export results"),
          uiOutput("export_note"),
          tags$div(
            class = "d-flex gap-2 flex-wrap",
            downloadButton("dl_xlsx", "Download Excel workbook (.xlsx)",
                           class = "btn-primary"),
            downloadButton("dl_csv", "Download CSV bundle (.zip)",
                           class = "btn-outline-light"),
            actionButton("save_disk", "\U0001F4C1 Save to workspace folder",
                         class = "btn-outline-light",
                         icon = NULL)
          ),
          tags$small(class = "text-muted d-block mt-2",
                     "Tip: file downloads don't work in VS Code's Simple Browser — use ",
                     tags$b("Save to workspace folder"), " (writes files into the project), ",
                     "or open the app in a real browser tab."),
          uiOutput("save_disk_path")
        ),

        # Accumulative performance curve (top-1..K) — the hallmark G-S-M plot.
        card(
          card_header("Accumulative performance: metric vs. number of top groups (k)"),
          layout_sidebar(
            sidebar = sidebar(
              width = 220, position = "right",
              selectInput("acc_metric", "Metric",
                          choices = c("Accuracy", "AUC", "Sensitivity", "Specificity",
                                      "Precision", "F_measure", "Cohens_kappa"),
                          selected = "Accuracy")
            ),
            plotOutput("plot_accumulative", height = "320px")
          )
        ),

        # Display Single Fold Mode Results
        conditionalPanel(
          condition = "!input.run_mccv_mode",
          layout_columns(
            col_widths = c(6, 6),
            card(
              card_header("Group Scores (by Accuracy)"),
              DTOutput("score_table")
            ),
            card(
              card_header("Score Distribution"),
              plotOutput("plot_scores", height = "300px")
            )
          ),
          card(
            card_header("Final Joint Model Performance (Combined Top Features)"),
            layout_columns(
              col_widths = c(4, 4, 4),
              uiOutput("model_vbox_acc"),
              uiOutput("model_vbox_kappa"),
              uiOutput("model_vbox_sens")
            ),
            layout_columns(
              col_widths = c(4, 4, 4),
              uiOutput("model_vbox_spec"),
              uiOutput("model_vbox_prec"),
              uiOutput("model_vbox_f1")
            )
          )
        ),
        
        # Display Parallel MCCV Mode Results
        conditionalPanel(
          condition = "input.run_mccv_mode",
          card(
            card_header("Parallel MCCV Performance Summary (Aggregate)"),
            DTOutput("mccv_summary_table")
          ),
          layout_columns(
            col_widths = c(6, 6),
            card(
              card_header("MCCV Folds Performance Visualizer"),
              plotOutput("plot_mccv_folds", height = "350px")
            ),
            card(
              card_header("Average Group Accuracy (Top 10 Stable)"),
              plotOutput("plot_mccv_groups", height = "350px")
            )
          ),
          card(
            card_header("Group Stability & Accuracies (All Scored Groups)"),
            DTOutput("mccv_stability_table")
          )
        )
      ),
      
      conditionalPanel(
        condition = "!output.score_done && !output.mccv_running",
        tags$div(
          style = "text-align:center; padding:80px 20px; opacity:0.5;",
          tags$h3("\U0001F3AF Run scoring or parallel MCCV analysis"),
          tags$p("Select your classifier and parameters, then click 'Run Analysis' to see results.")
        )
      )
    )
  ),

  # ── Footer ─────────────────────────────────────────────────────────────
  nav_spacer(),
  nav_item(
    tags$a(href = "https://github.com/jallmer/GSM4R",
           target = "_blank",
           style = "color: #888; font-size: 0.8rem;",
           "GitHub \U0001F517")
  )
)


# ══════════════════════════════════════════════════════════════════════════════
# SERVER
# ══════════════════════════════════════════════════════════════════════════════
server <- function(input, output, session) {

  # ── Reactive values ────────────────────────────────────────────────────
  rv <- reactiveValues(
    expr_data  = NULL,   # list(expr, labels)
    groups     = NULL,   # named list
    dataset    = NULL,   # gsm_dataset object
    loaded     = FALSE,

    # Collapse
    collapse_result = NULL,
    collapsed       = FALSE,

    # Pre-score
    prescore_groups = NULL,
    prescored       = FALSE,
    train_idx       = NULL,

    # Score
    score_results = NULL,
    scored        = FALSE
  )

  # Chunked pre-scoring worker state (drives the live main-panel view).
  ps <- reactiveValues(
    active = FALSE, i = 0, n = 0, tick = 0,
    gnames = NULL, scores = NULL, groups = NULL,
    expr = NULL, labels = NULL,
    current = "", top_percent = 20, t0 = NULL
  )

  # Batched MCCV worker state (runs `cores` folds per cycle, streams results).
  mc <- reactiveValues(
    active = FALSE, done = 0, n = 0, tick = 0,
    cores = 1, folds = NULL, params = NULL, t0 = NULL
  )

  # ── Output flags for conditionalPanel ──────────────────────────────────
  output$data_loaded     <- reactive(rv$loaded)
  output$collapse_done   <- reactive(rv$collapsed)
  output$prescore_done   <- reactive(rv$prescored)
  output$prescore_running <- reactive(isTRUE(ps$active))
  output$score_done      <- reactive(rv$scored)
  output$mccv_running    <- reactive(isTRUE(mc$active))
  outputOptions(output, "data_loaded",      suspendWhenHidden = FALSE)
  outputOptions(output, "collapse_done",    suspendWhenHidden = FALSE)
  outputOptions(output, "prescore_done",    suspendWhenHidden = FALSE)
  outputOptions(output, "prescore_running", suspendWhenHidden = FALSE)
  outputOptions(output, "score_done",       suspendWhenHidden = FALSE)
  outputOptions(output, "mccv_running",     suspendWhenHidden = FALSE)

  # ══════════════════════════════════════════════════════════════════════
  # Tab 1: Load Data
  # ══════════════════════════════════════════════════════════════════════
  # Shared loader used by both the file upload and the example dataset buttons.
  load_into_rv <- function(expr_path, groups_path, groups_name, fmt) {
    tryCatch({
      withProgress(message = "Loading data...", value = 0.3, {
        # Read expression
        expr_data <- read_expression(expr_path)
        setProgress(0.5, detail = "Reading groups...")

        # Read groups (resolve auto-detect from the file extension)
        if (fmt == "auto") {
          ext <- tolower(tools::file_ext(groups_name))
          fmt <- if (ext == "gmt") "gmt" else "csv"
        }
        groups <- if (fmt == "gmt") {
          read_gmt(groups_path)
        } else {
          read_groups(groups_path)
        }

        setProgress(0.8, detail = "Building dataset...")

        # Build dataset
        ds <- gsm_dataset(
          expr_data$expr, expr_data$labels, groups,
          min_group_size = input$min_group_size,
          max_group_size = input$max_group_size
        )

        rv$expr_data <- expr_data
        rv$groups    <- groups
        rv$dataset   <- ds
        rv$loaded    <- TRUE
        rv$dataset_name <- tools::file_path_sans_ext(basename(expr_path))

        # Reset downstream
        rv$collapsed <- FALSE
        rv$prescored <- FALSE
        rv$scored    <- FALSE

        setProgress(1, detail = "Done!")
      })

      showNotification("\u2705 Data loaded successfully!", type = "message")
    }, error = function(e) {
      showNotification(paste("\u274C Error:", e$message), type = "error",
                       duration = 10)
    })
  }

  observeEvent(input$btn_load, {
    req(input$expr_file, input$groups_file)
    load_into_rv(input$expr_file$datapath, input$groups_file$datapath,
                 input$groups_file$name, input$groups_format)
  })

  # Load a selected example: either the bundled toy data (inst/extdata) or one
  # of the repo-level example_input_files/ + example_grouping_files/ datasets.
  observeEvent(input$btn_example, {
    expr_sel <- input$example_expr
    grp_sel  <- input$example_groups

    expr_path <- if (expr_sel == "__toy__") {
      system.file("extdata", "test_expression.csv", package = "gsmR")
    } else file.path(.example_input_dir, expr_sel)

    grp_name <- if (grp_sel == "__toy__") "test_groups.csv" else grp_sel
    grp_path <- if (grp_sel == "__toy__") {
      system.file("extdata", "test_groups.csv", package = "gsmR")
    } else file.path(.example_groups_dir, grp_sel)

    if (!nzchar(expr_path) || !file.exists(expr_path)) {
      showNotification(paste("\u274C Expression file not found:", expr_sel),
                       type = "error", duration = 10); return()
    }
    if (!nzchar(grp_path) || !file.exists(grp_path)) {
      showNotification(paste("\u274C Grouping file not found:", grp_sel),
                       type = "error", duration = 10); return()
    }
    # Respect the format dropdown; "auto" resolves from the file extension.
    load_into_rv(expr_path, grp_path, grp_name, input$groups_format)
  })

  # Value boxes
  output$vbox_samples <- renderUI({
    req(rv$dataset)
    value_box_ui(rv$dataset$params$n_samples, "Samples", "\U0001F9EA")
  })
  output$vbox_genes <- renderUI({
    req(rv$dataset)
    value_box_ui(format(rv$dataset$params$n_genes, big.mark = ","), "Genes", "\U0001F9EC")
  })
  output$vbox_groups <- renderUI({
    req(rv$dataset)
    value_box_ui(format(rv$dataset$params$n_groups, big.mark = ","), "Groups", "\U0001F5C2\uFE0F")
  })
  output$vbox_classes <- renderUI({
    req(rv$dataset)
    value_box_ui(rv$dataset$params$n_classes, "Classes", "\U0001F3F7\uFE0F")
  })

  # Plots
  output$plot_class_dist <- renderPlot({
    req(rv$dataset)
    tbl <- table(rv$dataset$labels)
    par(bg = "transparent", col.axis = "#aaaaaa", col.lab = "#aaaaaa",
        col.main = "#cccccc", fg = "#555555", mar = c(4, 4, 2, 1))
    bp <- barplot(tbl, col = c("#4a90d9", "#e06666", "#6aa84f", "#f1c232")[
      seq_along(tbl)
    ], border = NA, main = "", ylab = "Count", las = 1)
    text(bp, tbl / 2, labels = tbl, col = "white", font = 2, cex = 1.2)
  }, bg = "transparent")

  output$plot_group_sizes <- renderPlot({
    req(rv$dataset)
    sizes <- vapply(rv$dataset$groups, length, integer(1))
    par(bg = "transparent", col.axis = "#aaaaaa", col.lab = "#aaaaaa",
        col.main = "#cccccc", fg = "#555555", mar = c(4, 4, 2, 1))
    hist(sizes, breaks = 30, col = "#4a90d9", border = NA,
         main = "", xlab = "Genes per group", ylab = "Frequency", las = 1)
    abline(v = median(sizes), col = "#f1c232", lwd = 2, lty = 2)
    legend("topright", legend = paste("Median:", median(sizes)),
           col = "#f1c232", lwd = 2, lty = 2, bty = "n", text.col = "#aaaaaa")
  }, bg = "transparent")

  # Expression preview
  output$expr_preview <- renderDT({
    req(rv$dataset)
    # Show first 100 genes x all samples
    mat <- rv$dataset$expr[1:min(100, nrow(rv$dataset$expr)), , drop = FALSE]
    df <- data.frame(Gene = rownames(mat), round(mat, 3), check.names = FALSE)
    datatable(df, options = list(
      pageLength = 10, scrollX = TRUE,
      dom = "ftp",
      initComplete = JS("function(settings, json) {",
                         "$(this.api().table().header()).css({'color': '#ccc'});",
                         "}")
    ), rownames = FALSE, class = "compact stripe")
  })

  # ══════════════════════════════════════════════════════════════════════
  # Tab 2: Collapse Replicates
  # ══════════════════════════════════════════════════════════════════════
  observeEvent(input$btn_collapse, {
    req(rv$loaded, input$meta_file)

    tryCatch({
      t0 <- Sys.time()
      withProgress(message = "Collapsing replicates...", value = 0.3, {
        # Read metadata
        meta <- data.table::fread(input$meta_file$datapath, header = TRUE,
                                  check.names = FALSE)
        meta <- as.data.frame(meta)

        setProgress(0.6, detail = "Running collapse...")

        result <- collapse_replicates(
          expr              = rv$dataset$expr,
          metadata          = meta,
          sample_id_col     = input$sample_id_col,
          replicate_group_col = input$replicate_group_col,
          method            = input$collapse_method,
          metadata_conflict = input$meta_conflict,
          keep_original     = input$keep_original
        )

        rv$collapse_result <- result
        rv$collapsed <- TRUE
        setProgress(1, detail = "Done!")
      })
      rv$time_collapse <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      showNotification(paste0("\u2705 Replicates collapsed in ",
                              format_duration(rv$time_collapse), "!"),
                       type = "message")
    }, error = function(e) {
      showNotification(paste("\u274C Error:", e$message), type = "error",
                       duration = 10)
    })
  })

  output$time_collapse <- renderUI({
    req(rv$collapsed, rv$time_collapse)
    time_banner_ui(rv$time_collapse, "Replicate collapse runtime")
  })

  output$vbox_collapsed_before <- renderUI({
    req(rv$collapse_result)
    value_box_ui(ncol(rv$dataset$expr), "Samples Before", "\u23EA")
  })
  output$vbox_collapsed_after <- renderUI({
    req(rv$collapse_result)
    value_box_ui(ncol(rv$collapse_result$expr), "Samples After", "\u23E9")
  })
  output$vbox_collapsed_groups <- renderUI({
    req(rv$collapse_result, rv$collapse_result$provenance)
    value_box_ui(nrow(rv$collapse_result$provenance), "Groups Collapsed", "\U0001F501")
  })

  output$provenance_table <- renderDT({
    req(rv$collapse_result, rv$collapse_result$provenance)
    datatable(rv$collapse_result$provenance, options = list(
      pageLength = 10, scrollX = TRUE, dom = "ftp"
    ), rownames = FALSE, class = "compact stripe")
  })

  output$collapse_warnings <- renderPrint({
    req(rv$collapse_result)
    w <- rv$collapse_result$warnings
    if (length(w) == 0) cat("No warnings.") else cat(paste(w, collapse = "\n"))
  })

  # ══════════════════════════════════════════════════════════════════════
  # Tab 3: Pre-Scoring
  # ══════════════════════════════════════════════════════════════════════
  # How many groups to score per reactive cycle. Larger = faster overall but
  # less frequent repaints; smaller = smoother live updates.
  PRESCORE_CHUNK <- 25

  observeEvent(input$btn_prescore, {
    req(rv$loaded)
    if (isTRUE(ps$active)) return()  # ignore re-clicks while running

    tryCatch({
      ds <- rv$dataset
      split <- split_train_test(ds$labels,
                                ratio = input$train_ratio,
                                seed = input$prescore_seed)
      rv$train_idx <- split$train
      rv$prescored <- FALSE   # hide old results; the live panel takes over
      rv$scored <- FALSE

      # "None" method is instant \u2014 no need for the live chunked worker.
      if (input$prescore_method != "anova") {
        t0 <- Sys.time()
        top_groups <- pre_score_groups(ds, train_idx = split$train,
                                       method = input$prescore_method,
                                       top_percent = input$top_percent)
        rv$prescore_groups <- top_groups
        rv$prescore_scores <- attr(top_groups, "scores")
        rv$prescore_nkeep  <- attr(top_groups, "n_keep")
        rv$prescored <- TRUE
        rv$time_prescore <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
        showNotification(paste0("\u2705 Pre-scoring done in ",
                                format_duration(rv$time_prescore), "! ",
                                length(top_groups), " groups selected."),
                         type = "message")
        return()
      }

      # ANOVA: set up the chunked worker so the main panel can fill live.
      ps$expr        <- ds$expr[, split$train, drop = FALSE]
      ps$labels      <- ds$labels[split$train]
      ps$groups      <- ds$groups
      ps$gnames      <- names(ds$groups)
      ps$n           <- length(ps$gnames)
      ps$scores      <- stats::setNames(rep(NA_real_, ps$n), ps$gnames)
      ps$i           <- 0
      ps$current     <- ""
      ps$top_percent <- input$top_percent
      ps$t0          <- Sys.time()
      ps$active      <- TRUE
    }, error = function(e) {
      ps$active <- FALSE
      showNotification(paste("\u274C Error:", e$message), type = "error",
                       duration = 10)
    })
  })

  # Worker: scores PRESCORE_CHUNK groups per cycle, yielding between chunks so
  # Shiny can repaint the live leaderboard/chart. Finalizes when all are done.
  observe({
    req(ps$active)
    err <- tryCatch({ isolate({
      if (ps$i >= ps$n) {
        # \u2500\u2500 Finalize \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
        sorted <- sort(ps$scores, decreasing = TRUE)
        n_keep <- max(1, ceiling(ps$n * ps$top_percent / 100))
        top_names <- names(sorted)[seq_len(n_keep)]
        res <- ps$groups[top_names]
        attr(res, "scores") <- sorted
        attr(res, "n_keep") <- n_keep

        rv$prescore_groups <- res
        rv$prescore_scores <- sorted
        rv$prescore_nkeep  <- n_keep
        rv$prescored <- TRUE
        rv$scored <- FALSE
        rv$time_prescore <- as.numeric(difftime(Sys.time(), ps$t0, units = "secs"))
        ps$active <- FALSE
        showNotification(paste0("\u2705 Pre-scoring done in ",
                                format_duration(rv$time_prescore), "! ",
                                n_keep, " groups selected."),
                         type = "message")
        return()
      }

      # \u2500\u2500 Process one chunk \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
      chunk <- min(PRESCORE_CHUNK, ps$n - ps$i)
      sc <- ps$scores
      for (k in seq_len(chunk)) {
        j <- ps$i + k
        gene_expr <- ps$expr[ps$groups[[ps$gnames[j]]], , drop = FALSE]
        pvals <- vapply(seq_len(nrow(gene_expr)), function(r) {
          tryCatch(stats::oneway.test(gene_expr[r, ] ~ ps$labels)$p.value,
                   error = function(e) 1.0)
        }, numeric(1))
        pvals[pvals == 0] <- .Machine$double.xmin
        sc[j] <- -2 * sum(log(pvals))
      }
      ps$scores  <- sc
      ps$i       <- ps$i + chunk
      ps$current <- ps$gnames[ps$i]
      # Bump the throttled tick ~every 4 chunks so the heavy live chart/table
      # repaint periodically instead of on every single chunk.
      if ((ps$i %/% chunk) %% 4 == 0 || ps$i >= ps$n) ps$tick <- ps$tick + 1
      NULL
    }) }, error = function(e) conditionMessage(e))

    if (!is.null(err)) {
      ps$active <- FALSE
      showNotification(paste("⛔ Pre-scoring stopped:", err), type = "error", duration = 12)
      return()
    }
    if (isTRUE(ps$active)) invalidateLater(20)  # breather so the client repaints
  })

  # ── Live pre-scoring renderers (active only while ps$active) ────────────
  output$prescore_live_status <- renderUI({
    req(ps$active)
    pct <- if (ps$n > 0) round(100 * ps$i / ps$n) else 0
    done_scores <- ps$scores[!is.na(ps$scores)]
    best <- if (length(done_scores)) {
      b <- which.max(done_scores)
      sprintf("%s (score = %.1f)", names(done_scores)[b], done_scores[b])
    } else "—"
    elapsed <- as.numeric(difftime(Sys.time(), ps$t0, units = "secs"))
    eta <- if (ps$i > 0) elapsed / ps$i * (ps$n - ps$i) else NA
    tags$div(
      class = "value-box",
      style = "background: rgba(74,144,217,0.12); border-left: 4px solid #4a90d9;",
      tags$div(style = "display:flex; justify-content:space-between; align-items:baseline;",
        tags$div(style = "font-size:1.5rem; font-weight:700;",
                 sprintf("⏳ Scoring groups… %d / %d  (%d%%)", ps$i, ps$n, pct)),
        tags$div(class = "text-muted",
                 sprintf("elapsed %s · ETA %s",
                         format_duration(elapsed), format_duration(eta)))
      ),
      # progress bar
      tags$div(style = "height:10px; background:rgba(255,255,255,0.08); border-radius:6px; margin:10px 0;",
        tags$div(style = sprintf("height:10px; width:%d%%; background:#4a90d9; border-radius:6px; transition:width .2s;", pct))
      ),
      tags$div(class = "label",
               HTML(sprintf("Now: <b>%s</b> &nbsp;|&nbsp; Best so far: <b>%s</b>",
                            htmltools::htmlEscape(ps$current), htmltools::htmlEscape(best))))
    )
  })

  output$plot_prescore_live <- renderPlot({
    ps$tick  # throttled dependency
    req(isolate(ps$active))
    done <- isolate(ps$scores)[!is.na(isolate(ps$scores))]
    req(length(done) > 0)
    topn <- min(input$prescore_topn_plot %||% 20, length(done))
    s <- rev(sort(done, decreasing = TRUE)[seq_len(topn)])
    par(bg = "transparent", col.axis = "#aaaaaa", col.lab = "#aaaaaa",
        fg = "#555555", mar = c(4, 9, 1, 1))
    barplot(s, horiz = TRUE, names.arg = names(s), las = 1,
            col = "#4a90d9", border = NA, cex.names = 0.7,
            xlab = "Fisher score (computed so far)")
  }, bg = "transparent")

  output$prescore_live_table <- renderTable({
    ps$tick  # throttled dependency
    req(isolate(ps$active))
    done <- sort(isolate(ps$scores)[!is.na(isolate(ps$scores))], decreasing = TRUE)
    req(length(done) > 0)
    done <- utils::head(done, 15)
    data.frame(
      Rank  = as.character(seq_along(done)),
      Group = names(done),
      Score = round(as.numeric(done), 2),
      stringsAsFactors = FALSE
    )
  }, striped = TRUE, hover = TRUE, width = "100%", digits = 2)

  output$time_prescore <- renderUI({
    req(rv$prescored, rv$time_prescore)
    time_banner_ui(rv$time_prescore, "Pre-scoring runtime")
  })

  output$vbox_prescore_input <- renderUI({
    req(rv$dataset, rv$prescored)
    value_box_ui(rv$dataset$params$n_groups, "Input Groups", "\U0001F4E5")
  })
  output$vbox_prescore_output <- renderUI({
    req(rv$prescore_groups)
    value_box_ui(length(rv$prescore_groups), "Selected Groups", "\u2705")
  })
  output$vbox_prescore_pct <- renderUI({
    req(rv$prescore_groups, rv$dataset)
    pct <- round(100 * length(rv$prescore_groups) / rv$dataset$params$n_groups, 1)
    value_box_ui(paste0(pct, "%"), "Retained", "\U0001F4CA")
  })

  # Flag for the conditionalPanel: do we have numeric scores to chart?
  output$prescore_has_scores <- reactive({
    isTRUE(rv$prescored) && !is.null(rv$prescore_scores) &&
      length(rv$prescore_scores) > 0
  })
  outputOptions(output, "prescore_has_scores", suspendWhenHidden = FALSE)

  output$prescore_table_header <- renderText({
    if (isTRUE(input$prescore_detail) && !is.null(rv$prescore_scores)) {
      "All Groups Ranked by Pre-Score (selected = kept for scoring)"
    } else {
      "Top-Ranked Groups After Pre-Scoring"
    }
  })

  output$prescore_table <- renderDT({
    req(rv$prescore_groups)
    detailed <- isTRUE(input$prescore_detail) && !is.null(rv$prescore_scores)

    if (detailed) {
      scores <- rv$prescore_scores                       # all groups, sorted desc
      gnames <- names(scores)
      n_keep <- rv$prescore_nkeep %||% length(rv$prescore_groups)
      gene_counts <- vapply(gnames, function(g) length(rv$dataset$groups[[g]]),
                            integer(1))
      df <- data.frame(
        Rank     = seq_along(scores),
        Group    = gnames,
        Score    = round(as.numeric(scores), 2),
        Genes    = gene_counts,
        Selected = ifelse(seq_along(scores) <= n_keep, "✅ yes", "—"),
        stringsAsFactors = FALSE
      )
      dt <- datatable(df, options = list(
        pageLength = 20, scrollX = TRUE, dom = "ftp",
        order = list(list(0, "asc"))
      ), rownames = FALSE, class = "compact stripe") |>
        formatStyle("Score",
                    background = styleColorBar(range(df$Score), "#4a90d9"),
                    backgroundSize = "90% 70%",
                    backgroundRepeat = "no-repeat",
                    backgroundPosition = "center") |>
        formatStyle("Selected",
                    target = "row",
                    backgroundColor = styleEqual("—", "rgba(255,255,255,0.02)"))
      return(dt)
    }

    # Simple view (detail off, or method = none)
    df <- data.frame(
      Rank  = seq_along(rv$prescore_groups),
      Group = names(rv$prescore_groups),
      Genes = vapply(rv$prescore_groups, length, integer(1)),
      stringsAsFactors = FALSE
    )
    datatable(df, options = list(
      pageLength = 20, scrollX = TRUE, dom = "ftp"
    ), rownames = FALSE, class = "compact stripe")
  })

  # Bar chart: top-N groups by Fisher score, selected vs not.
  output$plot_prescore_bars <- renderPlot({
    req(rv$prescore_scores)
    scores <- rv$prescore_scores
    n_keep <- rv$prescore_nkeep %||% length(rv$prescore_groups)
    topn <- min(input$prescore_topn_plot %||% 20, length(scores))
    s <- rev(scores[seq_len(topn)])                       # rev: highest at top
    is_sel <- rev(seq_len(topn) <= n_keep)
    cols <- ifelse(is_sel, "#4a90d9", "#555555")

    par(bg = "transparent", col.axis = "#aaaaaa", col.lab = "#aaaaaa",
        fg = "#555555", mar = c(4, 9, 1, 1))
    bp <- barplot(s, horiz = TRUE, names.arg = names(s), las = 1,
                  col = cols, border = NA, cex.names = 0.7,
                  xlab = "Fisher score  (-2 * sum log p)")
    legend("bottomright", legend = c("selected", "filtered out"),
           fill = c("#4a90d9", "#555555"), border = NA, bty = "n",
           text.col = "#aaaaaa")
  }, bg = "transparent")

  # Histogram: distribution of all group scores with the selection cutoff.
  output$plot_prescore_hist <- renderPlot({
    req(rv$prescore_scores)
    scores <- as.numeric(rv$prescore_scores)
    n_keep <- rv$prescore_nkeep %||% length(rv$prescore_groups)
    cutoff <- scores[min(n_keep, length(scores))]         # score of last kept group

    par(bg = "transparent", col.axis = "#aaaaaa", col.lab = "#aaaaaa",
        fg = "#555555", mar = c(4, 4, 1, 1))
    hist(scores, breaks = min(40, max(10, length(scores) %/% 5)),
         col = "#4a90d9", border = NA, main = "",
         xlab = "Fisher score", ylab = "Number of groups", las = 1)
    abline(v = cutoff, col = "#f1c232", lwd = 2, lty = 2)
    legend("topright", legend = sprintf("cutoff = %.1f (top %d)", cutoff, n_keep),
           col = "#f1c232", lwd = 2, lty = 2, bty = "n", text.col = "#aaaaaa")
  }, bg = "transparent")

  # ══════════════════════════════════════════════════════════════════════
  # Tab 4: Scoring
  # ══════════════════════════════════════════════════════════════════════
  
  # Dynamic package requirements notice in sidebar
  # Prominent, color-coded readiness alert shown at the top of the Score panel.
  alert_box <- function(type, title, body = NULL) {
    pal <- switch(type,
      ok    = list(bg = "rgba(106,168,79,0.15)",  bd = "#6aa84f", ic = "✅"),
      warn  = list(bg = "rgba(241,194,50,0.15)",  bd = "#f1c232", ic = "⚠️"),
      error = list(bg = "rgba(224,102,102,0.18)", bd = "#e06666", ic = "⛔"),
      info  = list(bg = "rgba(74,144,217,0.15)",  bd = "#4a90d9", ic = "ℹ️"))
    tags$div(
      class = "value-box",
      style = sprintf("background:%s; border-left:5px solid %s;", pal$bg, pal$bd),
      tags$div(style = "font-size:1.15rem; font-weight:600;",
               paste0(pal$ic, " ", title)),
      if (!is.null(body)) tags$div(class = "label", style = "margin-top:4px;", HTML(body))
    )
  }

  output$score_status_alert <- renderUI({
    method_lbl <- switch(input$score_method,
                         rf = "Random Forest", svm = "Support Vector Machine",
                         knn = "K-Nearest Neighbors", input$score_method)
    req_pkg <- switch(input$score_method, rf = "randomForest", svm = "e1071", knn = "class")
    mode_lbl <- if (isTRUE(input$run_mccv_mode))
      sprintf("Parallel MCCV (%d folds)", input$mccv_iterations %||% 0) else "Single-fold scoring"

    if (isTRUE(mc$active)) {
      return(alert_box("info", "Analysis running…",
                       "Parallel MCCV is in progress — see the live panel below."))
    }
    if (!isTRUE(rv$loaded)) {
      return(alert_box("warn", "Step 1 — Load a dataset",
                       "Go to the <b>📂 Data</b> tab and load (or pick) an expression + grouping file."))
    }
    if (!requireNamespace(req_pkg, quietly = TRUE)) {
      return(alert_box("error", paste0("Missing package: ", req_pkg),
        sprintf("%s needs the R package <b>%s</b>, which isn't installed.<br>Install it in R: <code>install.packages(\"%s\")</code>, then click Run Analysis.",
                method_lbl, req_pkg, req_pkg)))
    }
    # Pre-Score is OPTIONAL — Score creates its own split and pre-scores on the
    # fly when needed. We just describe what will happen.
    n_pre <- length(rv$prescore_groups %||% list())
    scope_txt <- if (identical(input$score_scope, "all")) {
      sprintf("Will score <b>ALL %d groups</b> (pre-score filter skipped)%s.",
              rv$dataset$params$n_groups,
              if (isTRUE(input$run_mccv_mode)) " per fold" else "")
    } else if (isTRUE(input$run_mccv_mode)) {
      est <- max(1, ceiling(rv$dataset$params$n_groups * input$top_percent / 100))
      heavy <- if (est > 100)
        sprintf(" <span style='color:#f1c232'>⚠️ ~%d groups/fold × %d folds — this can take many minutes and the app is busy during each parallel batch.</span>",
                est, input$mccv_iterations %||% 0) else ""
      sprintf("Will pre-score &amp; score the <b>top %s%%</b> (~%d) groups per fold.%s",
              input$top_percent, est, heavy)
    } else if (isTRUE(rv$prescored) && n_pre > 0) {
      sprintf("Will score the <b>top %d</b> of %d pre-scored groups.",
              min(input$n_top_groups %||% 10, n_pre), n_pre)
    } else {
      sprintf("No pre-scoring run yet — Score will pre-score automatically (top %s%%) and score the <b>top %d</b>. Run the <b>📊 Pre-Score</b> tab first if you want its detailed view.",
              input$top_percent, input$n_top_groups %||% 10)
    }
    alert_box("ok", paste0("Ready — ", method_lbl, " · ", mode_lbl),
              paste0(scope_txt, " Click <b>🚀 Run Analysis</b>."))
  })

  observeEvent(input$btn_score, {
    if (isTRUE(mc$active)) {
      showNotification("⏳ Parallel MCCV is already running — please wait for it to finish.",
                       type = "warning", duration = 6)
      return()
    }
    # Only a loaded dataset is required. Pre-Score is optional — Score creates
    # its own train/test split and pre-scores on the fly when needed.
    if (!isTRUE(rv$loaded)) {
      showNotification("⚠️ Load a dataset on the Data tab first.",
                       type = "warning", duration = 8)
      return()
    }

    # 1. Package dependency validation
    methods_to_check <- unique(c(input$score_method, 
                                 if (input$run_mccv_mode && input$model_method != "same") input$model_method else input$score_method))
    for (m in methods_to_check) {
      if (m == "rf" && !requireNamespace("randomForest", quietly = TRUE)) {
        showNotification("\u274C Package 'randomForest' required. Install with: install.packages('randomForest')",
                         type = "error", duration = 10)
        return()
      }
      if (m == "svm" && !requireNamespace("e1071", quietly = TRUE)) {
        showNotification("\u274C Package 'e1071' required. Install with: install.packages('e1071')",
                         type = "error", duration = 10)
        return()
      }
      if (m == "knn" && !requireNamespace("class", quietly = TRUE)) {
        showNotification("\u274C Package 'class' required.",
                         type = "error", duration = 10)
        return()
      }
    }

    # 2. Extract classifier arguments
    score_args <- list()
    if (input$score_method == "rf") {
      score_args <- list(ntree = input$rf_ntree)
    } else if (input$score_method == "svm") {
      score_args <- list(kernel = input$svm_kernel, cost = input$svm_cost, gamma = input$svm_gamma)
    } else if (input$score_method == "knn") {
      score_args <- list(k = input$knn_k)
    }

    tryCatch({
      t0 <- Sys.time()
      # Reset previous state
      rv$score_results <- NULL
      rv$final_model <- NULL
      rv$mccv_results <- NULL
      rv$scored <- FALSE

      if (!input$run_mccv_mode) {
        # ── SINGLE FOLD SCORING ─────────────────────────────────────────────
        # Pre-Score is optional: make a train/test split now if none exists.
        if (is.null(rv$train_idx)) {
          split <- split_train_test(rv$dataset$labels,
                                    ratio = input$train_ratio,
                                    seed = input$prescore_seed)
          rv$train_idx <- split$train
        }

        # Candidate groups: all groups, or the pre-scored top set (auto if needed).
        scope_all <- identical(input$score_scope, "all")
        auto_prescored <- FALSE
        if (scope_all) {
          candidates <- rv$dataset$groups
        } else {
          candidates <- rv$prescore_groups
          if (is.null(candidates)) {
            candidates <- pre_score_groups(rv$dataset, rv$train_idx,
                                           method = input$prescore_method,
                                           top_percent = input$top_percent)
            auto_prescored <- TRUE
          }
        }

        results <- list()
        final_model <- NULL; best_k <- NA_integer_; n_shaved <- NA_integer_
        withProgress(message = "Gene shaving + scoring + accumulative modeling...", value = 0, {
          # Optional leakage-free z-score (fit on train, applied to all columns).
          expr_use <- if (isTRUE(input$normalize_z))
            zscore_fit_apply(rv$dataset$expr, rv$train_idx) else rv$dataset$expr

          # Phase I — training-only gene shaving (Welch t-test + BH FDR).
          if (isTRUE(input$shave_genes)) {
            keep_genes <- shave_genes(expr_use, rv$dataset$labels,
                                      rv$train_idx, alpha = input$fdr_alpha)
            candidates <- shave_groups(candidates, keep_genes)
            n_shaved <- length(keep_genes)
            validate(need(length(candidates) > 0,
                          "No groups survived gene shaving — raise the FDR threshold."))
          }
          K <- min(input$n_top_groups, length(candidates))

          # Phase II — score each candidate (mean F1, k-fold CV, TRAIN only) -> rank.
          scores <- numeric(length(candidates)); names(scores) <- names(candidates)
          for (i in seq_along(candidates)) {
            gname <- names(candidates)[i]
            setProgress(0.5 * i / length(candidates),
                        detail = sprintf("Scoring %s (%d/%d)", gname, i, length(candidates)))
            scores[i] <- do.call(score_group_train, c(
              list(expr = expr_use, labels = rv$dataset$labels,
                   genes = candidates[[gname]], train_idx = rv$train_idx,
                   method = input$score_method), score_args))
            results[[gname]] <- list(accuracy = scores[i])
          }
          ranked <- names(sort(scores, decreasing = TRUE))

          # Phase III — accumulative top-1..K; skip k that add no new features.
          acc_rows <- list(); best_model <- NULL; prev_features <- 0L
          for (k in seq_len(K)) {
            setProgress(0.5 + 0.5 * k / K,
                        detail = sprintf("Modeling top-%d / %d groups", k, K))
            sub <- candidates[ranked[seq_len(k)]]
            if (length(unique(unlist(sub))) == prev_features) next
            prev_features <- length(unique(unlist(sub)))
            m <- do.call(build_model, c(
              list(expr = expr_use, labels = rv$dataset$labels,
                   top_groups = sub, train_idx = rv$train_idx,
                   method = input$score_method), score_args))
            acc_rows[[length(acc_rows) + 1L]] <- data.frame(
              k = k, group_added = ranked[k], n_features = length(m$features),
              Accuracy = m$metrics$Accuracy, Sensitivity = m$metrics$Sensitivity,
              Specificity = m$metrics$Specificity, Precision = m$metrics$Precision,
              F_measure = m$metrics$F_measure, Cohens_kappa = m$metrics$Cohens_kappa,
              AUC = m$metrics$AUC %||% NA_real_,
              stringsAsFactors = FALSE)
            m$.k <- k
            if (is.null(best_model) || m$metrics$F_measure > best_model$metrics$F_measure)
              best_model <- m
          }
          final_model      <- best_model
          best_k           <- best_model$.k
          rv$accumulative  <- do.call(rbind, acc_rows)
          rv$final_model   <- final_model
          rv$score_results <- results
          rv$scored <- TRUE
          setProgress(1, detail = "Done!")
        })

        n <- length(candidates)
        rv$score_summary <- list(
          mode = "single",
          dataset_name = rv$dataset_name %||% "dataset",
          modeling = sprintf("accumulative top-1..%d (best k = %d)", K, best_k),
          best_k = best_k,
          gene_shaving = if (isTRUE(input$shave_genes))
            sprintf("Welch t-test + BH FDR %.3g (%d genes kept)", input$fdr_alpha, n_shaved)
            else "off",
          scope = if (scope_all) "All groups (pre-score filter skipped)"
                  else if (auto_prescored) "Top pre-scored groups (auto, on the fly)"
                  else "Top pre-scored groups",
          classifier = classifier_label(input$score_method),
          params = format_params(score_args),
          ranking = "group classifier score (train-only)",
          n_candidates = n,
          n_groups_scored = K,
          n_features = length(final_model$features),
          prescore_method = if (scope_all) "skipped" else input$prescore_method,
          prescore_top_percent = input$top_percent,
          seed = input$score_seed,
          dataset = sprintf("%d samples x %d genes",
                            rv$dataset$params$n_samples, rv$dataset$params$n_genes)
        )
        rv$time_score <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
        showNotification(sprintf("✅ Accumulative top-1..%d done (best k = %d) in %s!",
                                 K, best_k, format_duration(rv$time_score)),
                         type = "message")
      } else {
        # ── PARALLEL MULTI-FOLD MCCV ────────────────────────────────────────
        m_method <- if (input$model_method == "same") input$score_method else input$model_method
        # Determine modeling arguments
        model_args <- list()
        if (m_method == "rf") {
          model_args <- list(ntree = input$rf_ntree)
        } else if (m_method == "svm") {
          model_args <- list(kernel = input$svm_kernel, cost = input$svm_cost, gamma = input$svm_gamma)
        } else if (m_method == "knn") {
          model_args <- list(k = input$knn_k)
        }

        # Scope: "all" groups => keep 100% in each fold's pre-scoring (no filter).
        scope_all <- identical(input$score_scope, "all")
        fold_top_pct <- if (scope_all) 100 else input$top_percent
        mc$scope_all <- scope_all
        # Hand off to the batched worker so the panel fills live, fold by fold.
        mc$params <- list(
          dataset = rv$dataset, train_ratio = input$mccv_train_ratio,
          pre_score_method = input$prescore_method,
          pre_score_top_percent = fold_top_pct,
          score_method = input$score_method, model_method = m_method,
          score_args = score_args, model_args = model_args,
          top_n = input$n_top_groups, shave = isTRUE(input$shave_genes),
          fdr_alpha = input$fdr_alpha, normalize = isTRUE(input$normalize_z),
          seed = input$score_seed
        )
        mc$n      <- input$mccv_iterations
        mc$cores  <- input$mccv_cores
        mc$done   <- 0
        mc$folds  <- vector("list", mc$n)
        mc$tick   <- 0
        mc$primed <- FALSE   # first worker tick only paints the live panel
        mc$t0     <- Sys.time()
        mc$active <- TRUE
        showNotification(sprintf("🚀 Starting Parallel MCCV: %d folds on %d core(s)…",
                                 mc$n, mc$cores), type = "message", duration = 5)
      }
    }, error = function(e) {
      showNotification(paste("\u274C Error during execution:", e$message), type = "error",
                       duration = 10)
    })
  })

  # \u2500\u2500 MCCV worker: run `cores` folds per cycle, stream results, then aggregate
  observe({
    req(mc$active)
    # First tick: do NO heavy work \u2014 just yield so Shiny paints the live
    # "running" panel before the (blocking) parallel batch starts.
    if (!isTRUE(isolate(mc$primed))) {
      isolate(mc$primed <- TRUE)
      invalidateLater(60)
      return()
    }
    err <- tryCatch({ isolate({
      if (mc$done >= mc$n) {
        # \u2500\u2500 Finalize: aggregate exactly like run_mccv() \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
        results <- mc$folds
        ok <- Filter(function(r) !is.null(r) && r$status == "success", results)
        if (length(ok) == 0) {
          mc$active <- FALSE
          showNotification("\u274C All MCCV folds failed.", type = "error", duration = 10)
          return()
        }
        metrics_df <- do.call(rbind, lapply(ok, function(r) as.data.frame(r$metrics)))
        rv$mccv_results <- list(
          folds = results,
          metrics_summary = data.frame(
            Mean = colMeans(metrics_df, na.rm = TRUE),
            SD   = apply(metrics_df, 2, stats::sd, na.rm = TRUE)
          ),
          raw_metrics = metrics_df
        )
        p <- mc$params
        groups_per_fold <- vapply(ok, function(r) r$n_groups %||% NA_integer_, integer(1))
        feats_per_fold  <- vapply(ok, function(r) r$features_count %||% NA_integer_, integer(1))
        best_ks <- vapply(ok, function(r) r$best_k %||% NA_integer_, integer(1))
        best_k_mode <- as.integer(names(sort(table(best_ks), decreasing = TRUE))[1])
        rv$score_summary <- list(
          mode = "mccv",
          modeling = sprintf("accumulative top-1..%d per fold", p$top_n %||% 10),
          ranking = "group classifier score (train-only)",
          gene_shaving = if (isTRUE(p$shave))
            sprintf("Welch t-test + BH FDR %.3g (mean %.0f genes/fold)",
                    p$fdr_alpha %||% 0.05,
                    mean(vapply(ok, function(r) r$n_shaved_genes %||% NA_real_, numeric(1)),
                         na.rm = TRUE)) else "off",
          best_k = best_k_mode,
          dataset_name = rv$dataset_name %||% "dataset",
          scope = if (isTRUE(mc$scope_all)) "All groups (pre-score keeps 100%/fold)"
                  else "Top pre-scored groups",
          classifier = classifier_label(p$score_method),
          model_classifier = classifier_label(p$model_method),
          params = format_params(p$score_args),
          model_params = format_params(p$model_args),
          n_iterations = mc$n,
          n_folds_ok = length(ok),
          cores = mc$cores,
          train_ratio = p$train_ratio,
          prescore_method = p$pre_score_method,
          prescore_top_percent = p$pre_score_top_percent,
          groups_per_fold = round(mean(groups_per_fold, na.rm = TRUE), 1),
          features_per_fold = round(mean(feats_per_fold, na.rm = TRUE), 1),
          seed = p$seed,
          dataset = sprintf("%d samples x %d genes",
                            rv$dataset$params$n_samples, rv$dataset$params$n_genes)
        )
        rv$scored <- TRUE
        rv$time_score <- as.numeric(difftime(Sys.time(), mc$t0, units = "secs"))
        mc$active <- FALSE
        showNotification(sprintf("\u2705 Parallel MCCV finished %d/%d folds in %s!",
                                 length(ok), mc$n, format_duration(rv$time_score)),
                         type = "message")
        return()
      }

      # \u2500\u2500 Run the next batch of `cores` folds in parallel \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
      batch <- min(mc$cores, mc$n - mc$done)
      idx <- mc$done + seq_len(batch)
      p <- mc$params
      runner <- function(i) mccv_run_fold(
        p$dataset, i, train_ratio = p$train_ratio,
        pre_score_method = p$pre_score_method,
        pre_score_top_percent = p$pre_score_top_percent,
        score_method = p$score_method, model_method = p$model_method,
        score_args = p$score_args, model_args = p$model_args,
        top_n = p$top_n %||% 10, shave = p$shave %||% TRUE,
        fdr_alpha = p$fdr_alpha %||% 0.05, normalize = p$normalize %||% FALSE,
        seed = p$seed)
      res <- if (mc$cores > 1 && .Platform$OS.type != "windows") {
        parallel::mclapply(idx, runner, mc.cores = mc$cores)
      } else {
        lapply(idx, runner)
      }
      for (k in seq_along(idx)) mc$folds[[idx[k]]] <- res[[k]]
      mc$done <- mc$done + batch
      mc$tick <- mc$tick + 1
      NULL
    }) }, error = function(e) conditionMessage(e))

    # On any unexpected error, stop cleanly so the button isn't dead-locked.
    if (!is.null(err)) {
      mc$active <- FALSE
      showNotification(paste("⛔ MCCV stopped:", err), type = "error", duration = 12)
      return()
    }
    if (isTRUE(mc$active)) invalidateLater(20)
  })

  # Helper: accuracy vector from completed folds (in fold order).
  mccv_done_folds <- function() {
    Filter(function(r) !is.null(r) && r$status == "success", mc$folds)
  }

  output$mccv_live_status <- renderUI({
    req(mc$active)
    pct <- if (mc$n > 0) round(100 * mc$done / mc$n) else 0
    elapsed <- as.numeric(difftime(Sys.time(), mc$t0, units = "secs"))
    eta <- if (mc$done > 0) elapsed / mc$done * (mc$n - mc$done) else NA
    ok <- mccv_done_folds()
    accs <- vapply(ok, function(r) as.numeric(r$metrics$Accuracy), numeric(1))
    best <- if (length(accs)) sprintf("%.1f%%", 100 * max(accs)) else "\u2014"
    meanacc <- if (length(accs)) sprintf("%.1f%%", 100 * mean(accs)) else "\u2014"
    tags$div(
      class = "value-box",
      style = "background: rgba(74,144,217,0.12); border-left: 4px solid #4a90d9;",
      tags$div(style = "display:flex; justify-content:space-between; align-items:baseline;",
        tags$div(style = "font-size:1.5rem; font-weight:700;",
                 sprintf("\u23F3 Running MCCV fold %d / %d  (%d%%)", mc$done, mc$n, pct)),
        tags$div(class = "text-muted",
                 sprintf("%d cores \u00B7 elapsed %s \u00B7 ETA %s",
                         mc$cores, format_duration(elapsed), format_duration(eta)))
      ),
      tags$div(style = "height:10px; background:rgba(255,255,255,0.08); border-radius:6px; margin:10px 0;",
        tags$div(style = sprintf("height:10px; width:%d%%; background:#4a90d9; border-radius:6px; transition:width .2s;", pct))
      ),
      tags$div(class = "label",
               HTML(sprintf("Folds done: <b>%d</b> &nbsp;|&nbsp; Mean accuracy: <b>%s</b> &nbsp;|&nbsp; Best: <b>%s</b>",
                            length(ok), meanacc, best)))
    )
  })

  output$plot_mccv_live <- renderPlot({
    mc$tick
    req(isolate(mc$active))
    ok <- isolate(mccv_done_folds())
    req(length(ok) > 0)
    accs <- vapply(ok, function(r) as.numeric(r$metrics$Accuracy), numeric(1))
    par(bg = "transparent", col.axis = "#aaaaaa", col.lab = "#aaaaaa",
        fg = "#555555", mar = c(4, 4, 2, 1))
    if (length(accs) >= 3) {
      boxplot(accs, horizontal = TRUE, col = "#4a90d9", border = "#cccccc",
              ylim = c(0, 1), main = "", xlab = "Fold accuracy")
      points(accs, rep(1, length(accs)), pch = 19,
             col = grDevices::adjustcolor("#f1c232", 0.7))
    } else {
      barplot(accs, names.arg = paste0("F", seq_along(accs)), col = "#4a90d9",
              border = NA, ylim = c(0, 1), ylab = "Accuracy")
    }
  }, bg = "transparent")

  output$mccv_live_means <- renderUI({
    mc$tick
    req(isolate(mc$active))
    ok <- isolate(mccv_done_folds())
    if (length(ok) == 0) return(tags$p(class = "text-muted", "Waiting for first fold\u2026"))
    md <- do.call(rbind, lapply(ok, function(r) as.data.frame(r$metrics)))
    means <- colMeans(md, na.rm = TRUE)
    tagList(lapply(names(means), function(nm) {
      tags$div(style = "display:flex; justify-content:space-between; padding:4px 0; border-bottom:1px solid rgba(255,255,255,0.05);",
        tags$span(class = "text-muted", gsub("_", " ", nm)),
        tags$b(sprintf("%.3f", means[[nm]])))
    }))
  })

  output$mccv_live_table <- renderTable({
    mc$tick
    req(isolate(mc$active))
    rows <- Filter(Negate(is.null), isolate(mc$folds))
    req(length(rows) > 0)
    df <- do.call(rbind, lapply(rows, function(r) {
      if (r$status == "success") {
        data.frame(Fold = r$fold, Status = "OK",
                   Best_k = as.character(r$best_k %||% NA),
                   Accuracy = round(as.numeric(r$metrics$Accuracy), 4),
                   Kappa = round(as.numeric(r$metrics$Cohens_kappa), 4),
                   Sensitivity = round(as.numeric(r$metrics$Sensitivity), 4),
                   Features = as.character(r$features_count %||% NA),
                   stringsAsFactors = FALSE)
      } else {
        data.frame(Fold = r$fold, Status = "FAIL", Best_k = "—",
                   Accuracy = NA, Kappa = NA, Sensitivity = NA, Features = "—",
                   stringsAsFactors = FALSE)
      }
    }))
    df <- df[order(df$Fold), ]
    df$Fold <- as.character(df$Fold)
    df
  }, striped = TRUE, hover = TRUE, width = "100%", digits = 4)

  # Run summary card: what algorithm / parameters / how many groups were used.
  output$score_summary <- renderUI({
    req(rv$scored, rv$score_summary)
    s <- rv$score_summary
    row <- function(k, v) tags$div(
      style = "display:flex; justify-content:space-between; gap:1rem; padding:5px 0; border-bottom:1px solid rgba(255,255,255,0.06);",
      tags$span(class = "text-muted", k), tags$b(style = "text-align:right;", v))

    items <- if (identical(s$mode, "single")) {
      list(
        row("Mode", "Single-fold scoring"),
        row("Modeling", s$modeling %||% "accumulative top-k"),
        row("Gene shaving", s$gene_shaving %||% "off"),
        row("Ranking", s$ranking %||% "group classifier score (train-only)"),
        row("Groups source", s$scope),
        row("Scoring algorithm", s$classifier),
        row("Hyperparameters", s$params),
        row("Candidate groups", format(s$n_candidates %||% NA, big.mark = ",")),
        row("Accumulative depth (K)", s$n_groups_scored),
        row("Best k (max accuracy)", s$best_k %||% NA),
        row("Features at best k", format(s$n_features, big.mark = ",")),
        row("Pre-scoring", if (identical(s$prescore_method, "skipped")) "skipped"
                           else sprintf("%s · top %s%%", toupper(s$prescore_method),
                                        s$prescore_top_percent)),
        row("Pre-scoring", if (identical(s$prescore_method, "skipped")) "skipped"
                           else sprintf("%s · top %s%%", toupper(s$prescore_method),
                                        s$prescore_top_percent)),
        row("Dataset", s$dataset),
        row("Random seed", s$seed),
        row("Runtime", format_duration(rv$time_score))
      )
    } else {
      list(
        row("Mode", sprintf("Parallel MCCV · %d cores", s$cores)),
        row("Modeling", s$modeling %||% "accumulative top-k per fold"),
        row("Gene shaving", s$gene_shaving %||% "off"),
        row("Ranking", s$ranking %||% "group classifier score (train-only)"),
        row("Most frequent best k", s$best_k %||% NA),
        row("Groups source", s$scope),
        row("Iterations (folds)", sprintf("%d (%d succeeded)",
                                          s$n_iterations, s$n_folds_ok)),
        row("Scoring algorithm", sprintf("%s (%s)", s$classifier, s$params)),
        row("Final model algorithm", sprintf("%s (%s)", s$model_classifier, s$model_params)),
        row("Avg groups scored / fold", s$groups_per_fold),
        row("Avg features / fold", s$features_per_fold),
        row("Pre-scoring", sprintf("%s · top %s%%", toupper(s$prescore_method),
                                   s$prescore_top_percent)),
        row("Train ratio", s$train_ratio),
        row("Dataset", s$dataset),
        row("Random seed", s$seed),
        row("Runtime", format_duration(rv$time_score))
      )
    }
    tags$div(
      class = "value-box",
      style = "background: rgba(74,144,217,0.08); border-left: 4px solid #4a90d9;",
      tags$div(style = "font-size:1.1rem; font-weight:600; margin-bottom:6px;",
               "\U0001F4CB Run summary"),
      items
    )
  })

  # Accumulative metric-vs-k curve (single-fold: one line; MCCV: mean +/- SD band).
  output$plot_accumulative <- renderPlot({
    req(rv$scored)
    metric <- input$acc_metric %||% "Accuracy"
    par(bg = "transparent", col.axis = "#aaaaaa", col.lab = "#aaaaaa",
        col.main = "#cccccc", fg = "#555555", mar = c(4, 4, 2, 1))
    if (isTRUE(input$run_mccv_mode)) {
      ok <- Filter(function(r) !is.null(r) && r$status == "success",
                   rv$mccv_results$folds %||% list())
      accs <- Filter(Negate(is.null), lapply(ok, `[[`, "accumulative"))
      req(length(accs) > 0)
      long <- do.call(rbind, accs)
      ks <- sort(unique(long$k))
      mean_v <- vapply(ks, function(kk) mean(long[long$k == kk, metric], na.rm = TRUE), numeric(1))
      sd_v   <- vapply(ks, function(kk) stats::sd(long[long$k == kk, metric], na.rm = TRUE), numeric(1))
      sd_v[is.na(sd_v)] <- 0
      plot(ks, mean_v, type = "n", ylim = c(0, 1), xlab = "Number of top groups (k)",
           ylab = metric, las = 1, main = "")
      polygon(c(ks, rev(ks)), c(mean_v - sd_v, rev(mean_v + sd_v)),
              col = grDevices::adjustcolor("#4a90d9", 0.25), border = NA)
      lines(ks, mean_v, col = "#4a90d9", lwd = 2)
      points(ks, mean_v, pch = 19, col = "#4a90d9")
      bk <- ks[which.max(mean_v)]
      abline(v = bk, col = "#f1c232", lwd = 2, lty = 2)
      legend("bottomright", bty = "n", text.col = "#aaaaaa",
             legend = sprintf("best k = %d (mean %s = %.3f)", bk, metric, max(mean_v)),
             col = "#f1c232", lwd = 2, lty = 2)
    } else {
      a <- rv$accumulative
      req(!is.null(a))
      plot(a$k, a[[metric]], type = "o", pch = 19, col = "#4a90d9", lwd = 2,
           ylim = c(0, 1), xlab = "Number of top groups (k)", ylab = metric,
           las = 1, main = "")
      bk <- a$k[which.max(a[[metric]])]
      abline(v = bk, col = "#f1c232", lwd = 2, lty = 2)
      legend("bottomright", bty = "n", text.col = "#aaaaaa",
             legend = sprintf("best k = %d (%s = %.3f)", bk, metric, max(a[[metric]], na.rm = TRUE)),
             col = "#f1c232", lwd = 2, lty = 2)
    }
  }, bg = "transparent")

  # ── Results export (Excel + CSV bundle) ────────────────────────────────
  # Assemble the in-memory result object for the export functions.
  assemble_result <- reactive({
    if (isTRUE(input$run_mccv_mode)) {
      req(rv$mccv_results); rv$mccv_results
    } else {
      req(rv$score_results, rv$final_model)
      list(score_results = rv$score_results, final_model = rv$final_model)
    }
  })

  output$export_note <- renderUI({
    req(rv$scored)
    if (isTRUE(input$run_mccv_mode)) {
      ok <- Filter(function(r) !is.null(r) && r$status == "success",
                   rv$mccv_results$folds %||% list())
      n_groups <- length(unique(unlist(lapply(ok, function(f) names(f$group_accuracies)))))
      tags$p(class = "text-muted",
             sprintf("Results kept: %d iterations × %d distinct groups scored. Includes per-iteration scores, ranked groups, group–feature associations, ranked features, and pre-score significance.",
                     length(ok), n_groups))
    } else {
      tags$p(class = "text-muted",
             sprintf("Results kept: %d scored groups + final joint model.",
                     length(rv$score_results %||% list())))
    }
  })

  .do_export <- function(file, fmt) {
    tryCatch({
      tbls <- gsm_export_tables(assemble_result(), rv$dataset, rv$score_summary)
      if (fmt == "xlsx") {
        gsm_write_xlsx(tbls, file)
      } else {
        prefix <- sprintf("%s_GSMwith%s",
                          rv$score_summary$dataset_name %||% "dataset",
                          switch(rv$score_summary$classifier %||% "model",
                                 "Random Forest" = "RandomForest",
                                 "Support Vector Machine" = "SVM",
                                 "K-Nearest Neighbors" = "KNN",
                                 rv$score_summary$classifier))
        gsm_write_csv_zip(tbls, file, prefix = prefix)
      }
    }, error = function(e) {
      showNotification(paste("❌ Export failed:", conditionMessage(e)),
                       type = "error", duration = 12)
      # Write a small note so the download isn't an empty/corrupt file.
      writeLines(paste("Export failed:", conditionMessage(e)), file)
    })
  }

  output$dl_xlsx <- downloadHandler(
    filename = function() gsm_export_filename(rv$score_summary, rv$dataset, "xlsx"),
    content  = function(file) .do_export(file, "xlsx")
  )
  output$dl_csv <- downloadHandler(
    filename = function() gsm_export_filename(rv$score_summary, rv$dataset, "zip"),
    content  = function(file) .do_export(file, "csv")
  )

  # Save-to-disk fallback (works in any browser, incl. VS Code Simple Browser):
  # writes the workbook + CSV bundle into a "gsm_outputs/" folder in the project.
  observeEvent(input$save_disk, {
    req(rv$scored)
    out_dir <- file.path(getOption("gsmR.launch_dir", getwd()), "gsm_outputs")
    tryCatch({
      dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
      tbls   <- gsm_export_tables(assemble_result(), rv$dataset, rv$score_summary)
      xlsx_f <- file.path(out_dir, gsm_export_filename(rv$score_summary, rv$dataset, "xlsx"))
      zip_f  <- file.path(out_dir, gsm_export_filename(rv$score_summary, rv$dataset, "zip"))
      gsm_write_xlsx(tbls, xlsx_f)
      gsm_write_csv_zip(tbls, zip_f,
                        prefix = sub("\\.zip$", "", basename(zip_f)))
      rv$saved_paths <- c(xlsx_f, zip_f)
      showNotification(paste0("✅ Saved results to: ", out_dir),
                       type = "message", duration = 12)
    }, error = function(e) {
      showNotification(paste("❌ Save failed:", conditionMessage(e)),
                       type = "error", duration = 12)
    })
  })

  output$save_disk_path <- renderUI({
    req(rv$saved_paths)
    tags$div(
      class = "value-box", style = "margin-top:10px; background: rgba(106,168,79,0.12); border-left:4px solid #6aa84f;",
      tags$div(style = "font-weight:600;", "✅ Files written to workspace:"),
      tags$ul(style = "margin:6px 0 0 0;",
              lapply(rv$saved_paths, function(p) tags$li(tags$code(p))))
    )
  })

  output$time_score <- renderUI({
    req(rv$scored, rv$time_score)
    lbl <- if (isTRUE(input$run_mccv_mode)) "Parallel MCCV runtime" else "Scoring + model runtime"
    time_banner_ui(rv$time_score, lbl)
  })

  # Single Fold Renderers: Group Scores Table
  output$score_table <- renderDT({
    req(rv$score_results)
    df <- data.frame(
      Rank     = seq_along(rv$score_results),
      Group    = names(rv$score_results),
      Accuracy = vapply(rv$score_results, \(x) round(x$accuracy, 4), numeric(1)),
      Genes    = vapply(names(rv$score_results),
                        \(g) length(rv$prescore_groups[[g]]), integer(1)),
      stringsAsFactors = FALSE
    )
    df <- df[order(-df$Accuracy), ]
    df$Rank <- seq_len(nrow(df))
    datatable(df, options = list(
      pageLength = 20, scrollX = TRUE, dom = "ftp"
    ), rownames = FALSE, class = "compact stripe") |>
      formatStyle("Accuracy",
                  background = styleColorBar(c(0, 1), "#4a90d9"),
                  backgroundSize = "90% 70%",
                  backgroundRepeat = "no-repeat",
                  backgroundPosition = "center")
  })

  # Single Fold Renderers: Barplot distribution
  output$plot_scores <- renderPlot({
    req(rv$score_results)
    accs <- vapply(rv$score_results, \(x) x$accuracy, numeric(1))
    accs <- sort(accs, decreasing = TRUE)
    par(bg = "transparent", col.axis = "#aaaaaa", col.lab = "#aaaaaa",
        col.main = "#cccccc", fg = "#555555", mar = c(8, 4, 2, 1))
    bp <- barplot(accs, col = "#4a90d9", border = NA,
                  main = "", ylab = "Accuracy", las = 2, ylim = c(0, 1),
                  names.arg = names(accs), cex.names = 0.7)
    abline(h = 0.5, col = "#888888", lty = 2)
  }, bg = "transparent")

  # Single Fold Renderers: Final Joint Model Value Boxes
  output$model_vbox_acc <- renderUI({
    req(rv$final_model)
    value_box_ui(paste0(round(rv$final_model$metrics$Accuracy * 100, 2), "%"), "Model Accuracy", "\U0001F3AF")
  })
  output$model_vbox_kappa <- renderUI({
    req(rv$final_model)
    value_box_ui(round(rv$final_model$metrics$Cohens_kappa, 4), "Cohen's Kappa", "\u2705")
  })
  output$model_vbox_sens <- renderUI({
    req(rv$final_model)
    value_box_ui(paste0(round(rv$final_model$metrics$Sensitivity * 100, 2), "%"), "Sensitivity (TPR)", "\u23E9")
  })
  output$model_vbox_spec <- renderUI({
    req(rv$final_model)
    value_box_ui(paste0(round(rv$final_model$metrics$Specificity * 100, 2), "%"), "Specificity (TNR)", "\u23EA")
  })
  output$model_vbox_prec <- renderUI({
    req(rv$final_model)
    val <- rv$final_model$metrics$Precision
    if (is.null(val) || is.na(val)) val <- 0
    value_box_ui(paste0(round(val * 100, 2), "%"), "Precision", "\U0001F3F7\uFE0F")
  })
  output$model_vbox_f1 <- renderUI({
    req(rv$final_model)
    val <- rv$final_model$metrics$F_measure
    if (is.null(val) || is.na(val)) val <- 0
    value_box_ui(paste0(round(val * 100, 2), "%"), "F-Measure", "\U0001F4CA")
  })

  # Parallel MCCV Renderers: Summary statistics table
  output$mccv_summary_table <- renderDT({
    req(rv$mccv_results)
    df <- rv$mccv_results$metrics_summary
    df$Metric <- rownames(df)
    df <- df[, c("Metric", "Mean", "SD")]
    df$Mean <- round(df$Mean, 4)
    df$SD <- round(df$SD, 4)
    datatable(df, options = list(
      dom = "t",
      initComplete = JS("function(settings, json) {",
                         "$(this.api().table().header()).css({'color': '#ccc'});",
                         "}")
    ), rownames = FALSE, class = "compact stripe")
  })

  # Parallel MCCV Renderers: Metrics spread boxplot
  output$plot_mccv_folds <- renderPlot({
    req(rv$mccv_results)
    df <- rv$mccv_results$raw_metrics
    par(bg = "transparent", col.axis = "#aaaaaa", col.lab = "#aaaaaa",
        col.main = "#cccccc", fg = "#555555", mar = c(6, 4, 2, 1))
    boxplot(df, col = c("#4a90d9", "#e06666", "#6aa84f", "#f1c232", "#8e7cc3", "#e06666"),
            border = "#cccccc", las = 2, main = "", ylab = "Value", ylim = c(0, 1))
  }, bg = "transparent")

  # Parallel MCCV Renderers: Stability ranking table
  output$mccv_stability_table <- renderDT({
    req(rv$mccv_results)
    mccv_res <- rv$mccv_results
    
    group_performance <- list()
    successful_folds <- Filter(function(res) res$status == "success", mccv_res$folds)
    n_success <- length(successful_folds)
    
    if (n_success == 0) return(NULL)
    
    for (i in seq_along(successful_folds)) {
      fold_res <- successful_folds[[i]]
      accs <- unlist(fold_res$group_accuracies)
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
    
    summary_list <- list()
    for (gname in names(group_performance)) {
      perf <- group_performance[[gname]]
      summary_list[[gname]] <- data.frame(
        Group = gname,
        Folds_Evaluated = length(perf$accuracies),
        Top_10_Appearances = perf$in_top_10,
        Appearance_Rate = round((perf$in_top_10 / n_success) * 100, 1),
        Mean_Accuracy = round(mean(perf$accuracies) * 100, 2),
        Mean_Rank = round(mean(perf$ranks), 1),
        stringsAsFactors = FALSE
      )
    }
    
    summary_df <- do.call(rbind, summary_list)
    sorted_summary <- summary_df[order(-summary_df$Top_10_Appearances, -summary_df$Mean_Accuracy, summary_df$Mean_Rank), ]
    rv$mccv_group_summary <- sorted_summary
    
    datatable(sorted_summary, options = list(
      pageLength = 10, scrollX = TRUE, dom = "ftp"
    ), rownames = FALSE, class = "compact stripe") |>
      formatStyle("Mean_Accuracy",
                  background = styleColorBar(c(0, 100), "#4a90d9"),
                  backgroundSize = "90% 70%",
                  backgroundRepeat = "no-repeat",
                  backgroundPosition = "center")
  })

  # Parallel MCCV Renderers: Stability top groups barplot
  output$plot_mccv_groups <- renderPlot({
    req(rv$mccv_group_summary)
    df <- rv$mccv_group_summary
    top_n <- min(10, nrow(df))
    df_sub <- df[seq_len(top_n), ]
    
    accs <- df_sub$Mean_Accuracy
    names(accs) <- df_sub$Group
    
    par(bg = "transparent", col.axis = "#aaaaaa", col.lab = "#aaaaaa",
        col.main = "#cccccc", fg = "#555555", mar = c(8, 4, 2, 1))
    bp <- barplot(accs, col = "#4a90d9", border = NA,
                  main = "", ylab = "Mean Accuracy (%)", las = 2, ylim = c(0, 100))
    abline(h = 50, col = "#888888", lty = 2)
  }, bg = "transparent")
}

# ══════════════════════════════════════════════════════════════════════════════
# Run
# ══════════════════════════════════════════════════════════════════════════════
shinyApp(ui = ui, server = server)
