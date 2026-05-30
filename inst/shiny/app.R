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
      .navbar { border-bottom: 1px solid rgba(255,255,255,0.06); }
      .nav-link.active { border-bottom: 2px solid #4a90d9 !important; }
    "))
  ),

  # ── Tab 1: Data Upload ───────────────────────────────────────────────────
  nav_panel(
    title = "\U0001F4C2 Data",
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
                     class = "btn-primary btn-lg w-100")
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
    title = "\U0001F501 Replicates",
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
    title = "\U0001F4CA Pre-Score",
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
        numericInput("prescore_seed", "Random seed", value = 42),
        tags$hr(),
        actionButton("btn_prescore", "\U0001F4CA Run Pre-Scoring",
                     class = "btn-primary btn-lg w-100")
      ),
      conditionalPanel(
        condition = "output.prescore_done",
        layout_columns(
          col_widths = c(4, 4, 4),
          uiOutput("vbox_prescore_input"),
          uiOutput("vbox_prescore_output"),
          uiOutput("vbox_prescore_pct")
        ),
        card(
          card_header("Top-Ranked Groups After Pre-Scoring"),
          DTOutput("prescore_table")
        )
      ),
      conditionalPanel(
        condition = "!output.prescore_done",
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
    title = "\U0001F3AF Score",
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
        
        numericInput("n_top_groups", "Score top N groups", value = 10, min = 1),
        numericInput("score_seed", "Random seed", value = 42),
        
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
                     class = "btn-primary btn-lg w-100"),
        tags$hr(),
        uiOutput("tab4_dependencies")
      ),
      
      # Main Panel Results Conditional Displays
      conditionalPanel(
        condition = "output.score_done",
        
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
        condition = "!output.score_done",
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

  # ── Output flags for conditionalPanel ──────────────────────────────────
  output$data_loaded   <- reactive(rv$loaded)
  output$collapse_done <- reactive(rv$collapsed)
  output$prescore_done <- reactive(rv$prescored)
  output$score_done    <- reactive(rv$scored)
  outputOptions(output, "data_loaded",   suspendWhenHidden = FALSE)
  outputOptions(output, "collapse_done", suspendWhenHidden = FALSE)
  outputOptions(output, "prescore_done", suspendWhenHidden = FALSE)
  outputOptions(output, "score_done",    suspendWhenHidden = FALSE)

  # ══════════════════════════════════════════════════════════════════════
  # Tab 1: Load Data
  # ══════════════════════════════════════════════════════════════════════
  observeEvent(input$btn_load, {
    req(input$expr_file, input$groups_file)

    tryCatch({
      withProgress(message = "Loading data...", value = 0.3, {
        # Read expression
        expr_data <- read_expression(input$expr_file$datapath)
        setProgress(0.5, detail = "Reading groups...")

        # Read groups
        fmt <- input$groups_format
        if (fmt == "auto") {
          ext <- tolower(tools::file_ext(input$groups_file$name))
          fmt <- if (ext == "gmt") "gmt" else "csv"
        }
        groups <- if (fmt == "gmt") {
          read_gmt(input$groups_file$datapath)
        } else {
          read_groups(input$groups_file$datapath)
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
    par(bg = "transparent", col.axis = "#aaa", col.lab = "#aaa",
        col.main = "#ccc", fg = "#555", mar = c(4, 4, 2, 1))
    bp <- barplot(tbl, col = c("#4a90d9", "#e06666", "#6aa84f", "#f1c232")[
      seq_along(tbl)
    ], border = NA, main = "", ylab = "Count", las = 1)
    text(bp, tbl / 2, labels = tbl, col = "white", font = 2, cex = 1.2)
  }, bg = "transparent")

  output$plot_group_sizes <- renderPlot({
    req(rv$dataset)
    sizes <- vapply(rv$dataset$groups, length, integer(1))
    par(bg = "transparent", col.axis = "#aaa", col.lab = "#aaa",
        col.main = "#ccc", fg = "#555", mar = c(4, 4, 2, 1))
    hist(sizes, breaks = 30, col = "#4a90d9", border = NA,
         main = "", xlab = "Genes per group", ylab = "Frequency", las = 1)
    abline(v = median(sizes), col = "#f1c232", lwd = 2, lty = 2)
    legend("topright", legend = paste("Median:", median(sizes)),
           col = "#f1c232", lwd = 2, lty = 2, bty = "n", text.col = "#aaa")
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
      showNotification("\u2705 Replicates collapsed!", type = "message")
    }, error = function(e) {
      showNotification(paste("\u274C Error:", e$message), type = "error",
                       duration = 10)
    })
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
  observeEvent(input$btn_prescore, {
    req(rv$loaded)

    tryCatch({
      withProgress(message = "Pre-scoring groups...", value = 0.2, {
        ds <- rv$dataset

        # Train/test split
        split <- split_train_test(ds$labels,
                                  ratio = input$train_ratio,
                                  seed = input$prescore_seed)
        rv$train_idx <- split$train

        setProgress(0.4, detail = paste0("Scoring ", ds$params$n_groups, " groups..."))

        top_groups <- pre_score_groups(
          ds,
          train_idx   = split$train,
          method      = input$prescore_method,
          top_percent = input$top_percent
        )

        rv$prescore_groups <- top_groups
        rv$prescored <- TRUE
        rv$scored <- FALSE  # Reset downstream

        setProgress(1, detail = "Done!")
      })
      showNotification(paste0("\u2705 Pre-scoring done! ",
                              length(rv$prescore_groups), " groups selected."),
                       type = "message")
    }, error = function(e) {
      showNotification(paste("\u274C Error:", e$message), type = "error",
                       duration = 10)
    })
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

  output$prescore_table <- renderDT({
    req(rv$prescore_groups)
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

  # ══════════════════════════════════════════════════════════════════════
  # Tab 4: Scoring
  # ══════════════════════════════════════════════════════════════════════
  
  # Dynamic package requirements notice in sidebar
  output$tab4_dependencies <- renderUI({
    req_pkg <- if (input$score_method == "rf") {
      "randomForest"
    } else if (input$score_method == "svm") {
      "e1071"
    } else {
      "class (core)"
    }
    tags$small(class = "text-muted", paste("Requires package:", req_pkg))
  })

  observeEvent(input$btn_score, {
    req(rv$prescored, rv$train_idx)

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
      # Reset previous state
      rv$score_results <- NULL
      rv$final_model <- NULL
      rv$mccv_results <- NULL
      rv$scored <- FALSE

      if (!input$run_mccv_mode) {
        # ── SINGLE FOLD SCORING ─────────────────────────────────────────────
        groups_to_score <- rv$prescore_groups
        n <- min(input$n_top_groups, length(groups_to_score))
        groups_to_score <- groups_to_score[seq_len(n)]

        results <- list()
        withProgress(message = "Scoring groups...", value = 0, {
          for (i in seq_along(groups_to_score)) {
            gname <- names(groups_to_score)[i]
            setProgress(i / n, detail = paste0(gname, " (", i, "/", n, ")"))

            set.seed(input$score_seed)
            score_args_combined <- c(
              list(
                expr      = rv$dataset$expr,
                labels    = rv$dataset$labels,
                genes     = groups_to_score[[i]],
                train_idx = rv$train_idx,
                method    = input$score_method
              ),
              score_args
            )
            res <- do.call(score_group, score_args_combined)
            results[[gname]] <- res
          }
          
          setProgress(0.9, detail = "Building final joint model...")
          
          # Build final joint model on combined top features
          m_method <- input$score_method
          m_args <- score_args
          
          model_args_combined <- c(
            list(
              expr       = rv$dataset$expr,
              labels     = rv$dataset$labels,
              top_groups = groups_to_score,
              train_idx  = rv$train_idx,
              method     = m_method
            ),
            m_args
          )
          final_model <- do.call(build_model, model_args_combined)
          
          rv$final_model <- final_model
          rv$score_results <- results
          rv$scored <- TRUE
          setProgress(1.0, detail = "Done!")
        })
        
        showNotification(paste0("\u2705 Scored ", n, " groups and built final model!"), type = "message")
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

        withProgress(message = "Executing Parallel MCCV...", value = 0.1, {
          setProgress(0.3, detail = sprintf("Running %d iterations across %d cores...", 
                                             input$mccv_iterations, input$mccv_cores))
          
          mccv_res <- run_mccv(
            dataset               = rv$dataset,
            n_iterations          = input$mccv_iterations,
            train_ratio           = input$mccv_train_ratio,
            pre_score_method      = input$prescore_method,
            pre_score_top_percent = input$top_percent,
            score_method          = input$score_method,
            model_method          = m_method,
            score_args            = score_args,
            model_args            = model_args,
            cores                 = input$mccv_cores,
            seed                  = input$score_seed
          )
          
          rv$mccv_results <- mccv_res
          rv$scored <- TRUE
          setProgress(1.0, detail = "Done!")
        })
        
        showNotification(sprintf("\u2705 Parallel MCCV finished %d iterations successfully!", 
                                 sum(vapply(mccv_res$folds, \(x) x$status == "success", logical(1)))), 
                         type = "message")
      }
    }, error = function(e) {
      showNotification(paste("\u274C Error during execution:", e$message), type = "error",
                       duration = 10)
    })
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
    par(bg = "transparent", col.axis = "#aaa", col.lab = "#aaa",
        col.main = "#ccc", fg = "#555", mar = c(8, 4, 2, 1))
    bp <- barplot(accs, col = "#4a90d9", border = NA,
                  main = "", ylab = "Accuracy", las = 2, ylim = c(0, 1),
                  names.arg = names(accs), cex.names = 0.7)
    abline(h = 0.5, col = "#888", lty = 2)
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
    par(bg = "transparent", col.axis = "#aaa", col.lab = "#aaa",
        col.main = "#ccc", fg = "#555", mar = c(6, 4, 2, 1))
    boxplot(df, col = c("#4a90d9", "#e06666", "#6aa84f", "#f1c232", "#8e7cc3", "#e06666"),
            border = "#ccc", las = 2, main = "", ylab = "Value", ylim = c(0, 1))
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
    
    par(bg = "transparent", col.axis = "#aaa", col.lab = "#aaa",
        col.main = "#ccc", fg = "#555", mar = c(8, 4, 2, 1))
    bp <- barplot(accs, col = "#4a90d9", border = NA,
                  main = "", ylab = "Mean Accuracy (%)", las = 2, ylim = c(0, 100))
    abline(h = 50, col = "#888", lty = 2)
  }, bg = "transparent")
}

# ══════════════════════════════════════════════════════════════════════════════
# Run
# ══════════════════════════════════════════════════════════════════════════════
shinyApp(ui = ui, server = server)
