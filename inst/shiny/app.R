# ──────────────────────────────────────────────────────────────────────────────
# gsmR Shiny App — Interactive G-S-M Workflow
# ──────────────────────────────────────────────────────────────────────────────

library(shiny)
library(bslib)
library(DT)

# ── Fork safety: force single-threaded BLAS ──────────────────────────────────
# The app runs heavy estimators (BRIDGE's `mb`/huge, glasso) inside forked
# workers (mcparallel) for the Stop button and parallel MCCV. Multi-threaded
# OpenBLAS + fork() deadlocks the child (BRIDGE "runs forever"). Setting the
# thread env vars here — in the parent, before any fork — makes forks safe.
Sys.setenv(OPENBLAS_NUM_THREADS = "1", OMP_NUM_THREADS = "1",
           GOTO_NUM_THREADS = "1")

# Increase maximum file upload size to 500MB
options(shiny.maxRequestSize = 500 * 1024^2)

# Dev mode: locate the package root robustly.
# Strategy 1: caller set GSMR_PKG_DIR (e.g. run-live.sh or our launcher).
# Strategy 2: walk up from this file's own path (works when source()'d).
# Strategy 3: walk up from getwd() — works when runApp() is called from the repo root.
.find_pkg_root <- function() {
  candidates <- c(
    Sys.getenv("GSMR_PKG_DIR", unset = NA_character_),
    tryCatch({
      frames <- sys.frames()
      ofiles <- Filter(Negate(is.null),
                       lapply(frames, function(f) f$ofile))
      if (length(ofiles))
        normalizePath(file.path(dirname(ofiles[[length(ofiles)]]), "..", ".."))
      else NA_character_
    }, error = function(e) NA_character_),
    normalizePath(file.path(getwd(), "inst", "shiny", ".."), mustWork = FALSE),
    normalizePath(file.path(getwd(), "..", ".."), mustWork = FALSE),
    normalizePath(getwd(), mustWork = FALSE)
  )
  for (d in candidates) {
    if (!is.na(d) && nzchar(d) && file.exists(file.path(d, "DESCRIPTION"))) return(d)
  }
  NULL
}
app_dir <- .find_pkg_root()

# Load gsmR — either installed or via pkgload
if (!is.null(app_dir)) {
  pkgload::load_all(app_dir, quiet = TRUE)
} else if (requireNamespace("gsmR", quietly = TRUE)) {
  library(gsmR)
} else {
  warning("gsmR package not found and could not locate package root for dev mode. ",
          "Run the app from the package root directory.")
}

# ══════════════════════════════════════════════════════════════════════════════
# Discover example datasets (repo-level example_*_files/ + bundled toy data)
# ══════════════════════════════════════════════════════════════════════════════
.launch_dir <- getOption("gsmR.launch_dir", if (!is.null(app_dir)) app_dir else getwd())
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

# Default grouping file for the dropdowns: prefer the DisGeNET cancer set
# (real gene symbols that match the GDS example datasets), else first example.
.default_groups <- {
  hit <- grep("disgenet", .example_groups, ignore.case = TRUE, value = TRUE)
  if (length(hit)) hit[1]
  else if (length(.example_groups)) .example_groups[1]
  else "__toy__"
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

# Null-coalescing helper (base R only ships %||% from 4.4.0).
`%||%` <- function(a, b) if (is.null(a)) b else a

# Friendly classifier name + a "key = value" string of its hyperparameters.
classifier_label <- function(method) {
  switch(method, rf = "Random Forest", svm = "Support Vector Machine",
         knn = "K-Nearest Neighbors", xgb = "XGBoost",
         mlp = "Neural Network (MLP)", gnn = "Graph Neural Network (GCN)",
         method)
}

# The classifier dropdowns (scoring, modeling, permutation, batch) all offer the
# same backends, so they are built from one list rather than five copies.
CLASSIFIER_CHOICES <- c("Random Forest" = "rf",
                        "Support Vector Machine" = "svm",
                        "K-Nearest Neighbors" = "knn",
                        "XGBoost" = "xgb",
                        "Neural Network (MLP)" = "mlp",
                        "Graph Neural Network (GCN)" = "gnn")

# Hyperparameter args for the torch-backed classifiers. `graph` is deliberately
# absent: MCCV injects BRIDGE's per-fold differential network at run time, and
# without one the GNN builds a train-only correlation graph itself.
build_nn_args <- function(input, method) {
  common <- list(dropout = input$nn_dropout %||% 0.3,
                 epochs  = input$nn_epochs  %||% 100,
                 lr      = input$nn_lr      %||% 0.01)
  if (method == "mlp") {
    hid <- suppressWarnings(as.integer(trimws(strsplit(
      input$mlp_hidden %||% "64,32", ",")[[1]])))
    hid <- hid[is.finite(hid) & hid > 0]
    c(list(hidden = if (length(hid)) hid else c(64L, 32L)), common)
  } else {
    c(list(hidden = as.integer(input$gnn_hidden %||% 32),
           cor_threshold = input$gnn_cor_threshold %||% 0.3), common)
  }
}

# Assemble BRIDGE hyperparameters from prefixed inputs, so the Data tab
# (prefix "bridge") and the Grouping tab (prefix "g_bridge") build an identical
# args list. min/max group size are passed explicitly (their input ids differ).
build_bridge_args <- function(input, prefix, min_size, max_size) {
  g <- function(suf, default) input[[paste0(prefix, "_", suf)]] %||% default
  res_grid <- as.numeric(strsplit(g("resolution", "0.5,1,1.5,2"),
                                  "[,[:space:]]+")[[1]])
  res_grid <- res_grid[is.finite(res_grid) & res_grid > 0]
  if (!length(res_grid)) res_grid <- c(0.5, 1, 1.5, 2)
  mg <- suppressWarnings(as.integer(g("max_genes", NA)))
  if (length(mg) == 0 || is.na(mg) || mg <= 0) mg <- NULL
  list(
    lambda = g("lambda", 0.1), beta = g("beta", 0.5),
    lambda2 = g("lambda2", 0.05), lambda1 = g("lambda1", 0.25),
    resolution = res_grid, tau = g("tau", 0.6),
    data_type = g("data_type", "gaussian"),
    precision_method = g("precision", "glasso"),
    prefilter = g("prefilter", "none"), max_genes = mg,
    min_size = min_size, max_size = max_size, verbose = FALSE
  )
}

# A prominent "how to use this" info callout (for tab guidance).
gsm_help_callout <- function(title, body) {
  tags$div(class = "card",
    style = "margin-bottom:12px; padding:10px 14px; background:rgba(74,144,217,0.08); border-left:4px solid #4a90d9;",
    tags$div(style = "font-weight:700; margin-bottom:4px;", paste0("ℹ️ ", title)),
    tags$div(class = "text-muted", style = "font-size:0.85rem; line-height:1.5;", body))
}

# Snapshot the run settings (with concrete numbers) so the GSM tab can show
# exactly what a run used — grouping, shaving, pre-score, scoring, modeling.
build_run_config <- function(input, rv) {
  ds <- rv$dataset
  gm <- ds$params$grouping_method %||% "file"

  # -- Grouping section --
  grouping <- list("Method" = if (gm == "bridge") "BRIDGE (data-driven)" else "From file")
  if (gm == "bridge") {
    spec <- ds$grouping_spec
    a <- spec$args %||% list()
    prior_n <- length(unique(unlist(spec$prior)))
    grouping[["Prior graph"]] <- if (prior_n > 0)
      sprintf("%s (%d genes)", input$bridge_prior_source %||% "none", prior_n)
      else (input$bridge_prior_source %||% "none")
    grouping[["Speed pre-filter"]] <- a$prefilter %||% "none"
    if (!is.null(a$max_genes)) grouping[["Max genes"]] <- a$max_genes
    grouping[["Estimator"]] <- a$precision_method %||% "glasso"
    grouping[["Data type"]] <- a$data_type %||% "gaussian"
    grouping[["lambda / beta"]] <- sprintf("%.3g / %.3g", a$lambda %||% 0.1, a$beta %||% 0.5)
    grouping[["tau (stability)"]] <- a$tau %||% 0.6
    grouping[["Resolutions"]] <- paste(a$resolution %||% c(0.5,1,1.5,2), collapse = ", ")
    grouping[["Groups"]] <- "re-derived per training fold"
  } else {
    grouping[["Groups loaded"]] <- ds$params$n_groups
    grouping[["Group size range"]] <- sprintf("%d – %d",
      ds$params$min_group_size %||% 5, ds$params$max_group_size %||% 500)
  }

  # -- Gene shaving section --
  if (isTRUE(input$shave_genes)) {
    shaving <- list("Enabled" = "Yes (Welch t-test, train-only)")
    if (identical(input$shave_method, "topk"))
      shaving[["Rule"]] <- sprintf("Top-K = %s genes", input$shave_k %||% 100)
    else
      shaving[["Rule"]] <- sprintf("BH-FDR alpha = %s", input$fdr_alpha %||% 0.05)
  } else shaving <- list("Enabled" = "No")

  # -- Scoring / modeling section --
  mccv <- isTRUE(input$run_mccv_mode)
  scoring <- list(
    "Mode" = if (mccv) "Parallel MCCV" else "Single split",
    "Pre-score" = sprintf("%s (top %s%%)", input$prescore_method %||% "none",
                          if (identical(input$score_scope, "all")) "100" else input$top_percent),
    "Group scorer" = classifier_label(input$score_method),
    "Final model" = classifier_label(if (mccv && input$model_method != "same")
                                       input$model_method else input$score_method),
    "Accumulative K" = input$n_top_groups,
    "Z-score norm" = if (isTRUE(input$normalize_z)) "on (train-fit)" else "off"
  )
  if (mccv) {
    scoring[["Iterations"]] <- input$mccv_iterations
    scoring[["Train ratio"]] <- input$mccv_train_ratio
    scoring[["Cores"]] <- input$mccv_cores
  } else {
    scoring[["Train ratio"]] <- input$train_ratio
  }

  list(grouping = grouping, shaving = shaving, scoring = scoring,
       dataset = sprintf("%d samples × %d genes",
                         ds$params$n_samples, ds$params$n_genes))
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
    style = "background: rgba(74,144,217,0.12); border-left: 4px solid #4a90d9; position:relative;",
    tags$button(type = "button", class = "btn-close",
                style = "position:absolute; top:8px; right:8px; filter:invert(1); opacity:0.6;",
                onclick = "this.parentElement.style.display='none';",
                `aria-label` = "Close"),
    tags$div(style = "font-size:1.4rem; font-weight:600;",
             paste0("⏱️ ", format_duration(secs))),
    tags$div(class = "label", label)
  )
}

# Inline help tooltip: gray circle with "?" that shows a popup on hover.
help_tip <- function(text) {
  tags$span(
    class = "help-tip",
    `data-tip` = text,
    "?"
  )
}

# A field label with an inline info (?) icon; hover/focus the icon for `tip`.
# Use above an input given `label = NULL` to attach contextual help.
field_label <- function(text, tip) {
  tags$div(style = "display:flex; align-items:center; gap:4px; margin-bottom:4px;",
    tags$label(text, style = "margin:0;"),
    help_tip(tip)
  )
}

# Turn a raw batch failure message into a plain-language "what to do" suggestion.
batch_fail_hint <- function(msg) {
  m <- tolower(paste(msg, collapse = " "))
  if (grepl("no groups remain|gene id matching|genes in groups|after filtering", m))
    return(paste0("The grouping file's genes don't match these datasets. ",
                  "In “Grouping file” pick a real file (e.g. KEGG, mate, or ",
                  "DisGeNET) instead of the bundled toy, then re-run."))
  if (grepl("all folds failed|all mccv folds failed", m))
    return(paste0("Every cross-validation fold failed — usually too few samples ",
                  "per class or a single-class split. Check the dataset has ≥2 ",
                  "classes with enough samples, or lower MCCV iterations / train ratio."))
  if (grepl("2 distinct classes|distinct classes|class column|nlevels", m))
    return(paste0("The expression file needs a ‘class’ column with at least two ",
                  "classes. Verify the first column holds the labels."))
  if (grepl("not found|cannot open|no such file|cannot read", m))
    return("A file couldn't be read — re-select or re-upload the dataset / grouping file.")
  paste0("See the Batch log for details. Expected formats: expression = samples as rows ",
         "with a first ‘class’ column; grouping = two columns (gene, group) or GMT.")
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

# Utility tabs (Session, Batch): no step number — they are tools/modes, not
# sequential pipeline steps, so they read differently from the numbered flow.
nav_util_title <- function(icon, label) {
  tags$span(
    class = "nav-step-title nav-util",
    tags$span(class = "nav-ico", icon),
    tags$span(class = "nav-lbl", label)
  )
}

# A small chevron shown between the numbered pipeline steps to signal order.
nav_flow_sep <- function() {
  bslib::nav_item(tags$span(class = "nav-step-sep", "›"))
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
  # Fixed tab bar: never collapse into a hamburger; tabs stay visible at any width.
  navbar_options = navbar_options(collapsible = FALSE),
  theme = bs_theme(
    version = 5,
    bootswatch = "darkly",
    primary = "#4a90d9",
    secondary = "#5b6480",
    success = "#0f9b58",
    info = "#4a90d9",
    warning = "#f4a623",
    "font-size-base" = "0.92rem",
    "border-radius" = "0.75rem",
    "font-family-sans-serif" =
      "system-ui, -apple-system, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif"
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

      /* ── Fixed tab bar: one row, never collapse; scroll if truly too narrow ── */
      .navbar > .container-fluid { flex-wrap: nowrap; }
      .navbar .navbar-collapse { display: flex !important; flex-basis: auto; }
      .navbar-nav { flex-wrap: nowrap; align-items: center; }
      .navbar-nav.me-auto, .navbar .navbar-nav:first-of-type {
        overflow-x: auto; overflow-y: hidden; scrollbar-width: thin;
      }
      .navbar-nav .nav-lbl { white-space: nowrap; }

      /* Utility tabs (right group): muted, un-numbered, separated from the flow. */
      .nav-util { opacity: 0.92; }
      .nav-util .nav-lbl { font-weight: 500; }
      .nav-vrule {
        display: inline-block; width: 1px; height: 24px;
        background: rgba(255,255,255,0.15); margin: 0 8px; align-self: center;
      }
      
      /* Sidebar overflow visible for tooltips */
      .sidebar { overflow: visible !important; }

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

      /* -- Help-tip circle -- */
      .help-tip {
        display: inline-flex; align-items: center; justify-content: center;
        width: 18px; height: 18px; border-radius: 50%;
        background: rgba(255,255,255,0.18); color: rgba(255,255,255,0.6);
        font-size: 0.7rem; font-weight: 700; cursor: help;
        flex-shrink: 0; margin-left: 4px; vertical-align: middle;
        transition: background 0.15s;
      }
      .help-tip:hover { background: rgba(255,255,255,0.30); color: #fff; }
      /* Tooltip popup styled by JS */
      #help-tip-popup {
        position: fixed; display: none;
        background: #1e2235; color: #ddd;
        font-size: 0.78rem; font-weight: 400;
        padding: 8px 12px; border-radius: 8px;
        max-width: 300px;
        box-shadow: 0 4px 16px rgba(0,0,0,0.5);
        border: 1px solid rgba(255,255,255,0.12);
        z-index: 99999; line-height: 1.4;
        pointer-events: none;
      }

      /* -- visNetwork navigation buttons (dark theme override) -- */
      .vis-button {
        background-color: rgba(40,44,66,0.9) !important;
        border: 1px solid rgba(255,255,255,0.25) !important;
        border-radius: 6px !important;
        box-shadow: 0 2px 8px rgba(0,0,0,0.4) !important;
        filter: grayscale(1) invert(1) brightness(1.5) !important;
      }
      .vis-button:hover {
        background-color: rgba(74,144,217,0.5) !important;
        border-color: #4a90d9 !important;
        filter: grayscale(0) invert(1) brightness(2) !important;
      }

      /* -- visNetwork tooltip (dark theme) -- */
      .vis-tooltip {
        background: #1e2235 !important;
        color: #e0e0e0 !important;
        border: 1px solid rgba(255,255,255,0.15) !important;
        border-radius: 8px !important;
        padding: 8px 12px !important;
        font-size: 0.85rem !important;
        box-shadow: 0 4px 16px rgba(0,0,0,0.5) !important;
      }

      /* -- Live progress notification (withProgress): prominent + themed -- */
      #shiny-notification-panel {
        top: 74px !important; bottom: auto !important;
        right: 50% !important; transform: translateX(50%);
        width: min(560px, 92vw); max-width: 560px;
      }
      .shiny-notification {
        background: #161a2e !important; color: #e8e8e8 !important;
        border: 1px solid rgba(74,144,217,0.55) !important;
        border-left: 4px solid #4a90d9 !important;
        border-radius: 10px !important;
        box-shadow: 0 8px 28px rgba(0,0,0,0.55) !important;
        font-size: 0.92rem; padding: 14px 18px !important;
      }
      .shiny-progress .progress {
        height: 14px !important; margin-top: 8px !important;
        background: rgba(255,255,255,0.10) !important; border-radius: 7px;
      }
      .shiny-progress .progress-bar {
        background-color: #4a90d9 !important;
        transition: width .25s ease;
      }
      .shiny-progress .progress-message { font-weight: 700; font-size: 0.95rem; }
      .shiny-progress .progress-detail  { color: #9fb4d6; margin-top: 2px; display:block; }

      /* -- Card fullscreen button: move to header bar -- */
      .bslib-full-screen-enter {
        position: absolute !important;
        top: 6px !important; right: 8px !important;
        bottom: auto !important;
        opacity: 0.5; z-index: 10;
        background: rgba(255,255,255,0.1) !important;
        border: 1px solid rgba(255,255,255,0.2) !important;
        border-radius: 4px !important;
        padding: 2px 6px !important;
        transition: opacity 0.15s;
      }
      .bslib-full-screen-enter:hover { opacity: 1; }
      .card { position: relative; }
    ")),
    tags$div(id = "help-tip-popup"),
    tags$script(HTML("
      document.addEventListener('mouseover', function(e) {
        var tip = e.target.closest('.help-tip');
        var popup = document.getElementById('help-tip-popup');
        if (!tip) { popup.style.display = 'none'; return; }
        var text = tip.getAttribute('data-tip');
        if (!text) return;
        popup.textContent = text;
        popup.style.display = 'block';
        var r = tip.getBoundingClientRect();
        var pw = popup.offsetWidth, ph = popup.offsetHeight;
        var left = r.right + 8;
        var top = r.top - ph / 2 + r.height / 2;
        if (left + pw > window.innerWidth - 10) left = r.left - pw - 8;
        if (top < 5) top = 5;
        if (top + ph > window.innerHeight - 5) top = window.innerHeight - ph - 5;
        popup.style.left = left + 'px';
        popup.style.top = top + 'px';
      });
      document.addEventListener('mouseout', function(e) {
        if (e.target.closest('.help-tip')) document.getElementById('help-tip-popup').style.display = 'none';
      });

      /* Button disable/enable helpers called from R via session$sendCustomMessage */
      Shiny.addCustomMessageHandler('disableBtn', function(id) {
        var btn = document.getElementById(id);
        if (btn) { btn.disabled = true; btn.classList.add('disabled'); btn.style.opacity = '0.5'; }
      });
      Shiny.addCustomMessageHandler('enableBtn', function(id) {
        var btn = document.getElementById(id);
        if (btn) { btn.disabled = false; btn.classList.remove('disabled'); btn.style.opacity = '1'; }
      });

      /* Show an inline progress banner inside a Shiny output container */
      Shiny.addCustomMessageHandler('showBanner', function(msg) {
        var el = document.getElementById(msg.id);
        if (el) {
          el.innerHTML =
            '<div class=value-box style=background:rgba(74,144,217,0.12);border-left:4px_solid_#4a90d9;padding:12px>' +
            '<div style=font-weight:600;font-size:1.1rem>' + '\u23f3 ' + msg.title + '</div>' +
            (msg.detail ? '<div style=margin-top:4px;opacity:0.7>' + msg.detail + '</div>' : '') +
            '</div>';
        }
      });

      /* Make the bslib sidebar collapse toggle discoverable: give it an
         accessible label/tooltip and a 'Settings' hint shown when collapsed. */
      function gsmLabelToggles() {
        document.querySelectorAll('.bslib-sidebar-layout > .collapse-toggle').forEach(function(btn) {
          if (btn.dataset.gsmLabeled) return;
          btn.dataset.gsmLabeled = '1';
          btn.setAttribute('title', 'Show / hide the settings panel');
          btn.setAttribute('aria-label', 'Show or hide the settings panel');
          var layout = btn.parentElement;
          if (layout && !layout.querySelector('.gsm-sidebar-hint')) {
            var hint = document.createElement('span');
            hint.className = 'gsm-sidebar-hint';
            hint.textContent = 'Settings';
            layout.appendChild(hint);
          }
        });
      }
      document.addEventListener('DOMContentLoaded', gsmLabelToggles);
      setTimeout(gsmLabelToggles, 400);
      setTimeout(gsmLabelToggles, 1200);
      new MutationObserver(gsmLabelToggles).observe(document.body, { childList: true, subtree: true });
    "))
  ),

  # ── Tab 1: Data Upload ───────────────────────────────────────────────────
  nav_panel(
    value = "data",
    title = nav_step_title("1", "\U0001F4C2", "Data"),
    layout_sidebar(
      sidebar = sidebar(
        title = "\U2699\UFE0F Data — Setup",
        width = 340,
        # -- Expression file --
        tags$div(
          class = "upload-area",
          tags$div(style = "display:flex; align-items:center; gap:4px; ",
            tags$label("\U0001F9EC Expression Data"),
            help_tip("CSV/TSV with samples as rows. First column = class labels (e.g. tumor/normal). All other columns = gene expression values.")
          ),
          fileInput("expr_file", NULL, accept = c(".csv", ".tsv", ".txt")),
          uiOutput("expr_file_status")
        ),
        tags$hr(),

        # ── Data filtering (label-blind denoising, applied at load) ──────
        tags$div(style = "display:flex; align-items:center; gap:4px;",
          tags$label("\U0001F9F9 Data filtering", style = "margin:0; font-weight:600;"),
          help_tip("Reduce noise by dropping uninformative genes. Every filter is label-blind (it never looks at the class), so it does NOT cause overfitting. Applied to all samples when you load.")
        ),
        checkboxInput("flt_missing", "Drop constant / mostly-missing genes", value = TRUE),
        checkboxInput("flt_lowexpr", "Drop low-expression genes", value = FALSE),
        conditionalPanel(
          condition = "input.flt_lowexpr",
          numericInput("flt_min_mean", "Min mean expression", value = 0, step = 0.5)
        ),
        checkboxInput("flt_variance", "Keep only the most variable genes", value = FALSE),
        conditionalPanel(
          condition = "input.flt_variance",
          numericInput("flt_max_genes", "Max genes to keep (top variance)",
                       value = 5000, min = 50, step = 500)
        ),
        uiOutput("filter_summary"),
        tags$hr(),
        actionButton("btn_load", "\U0001F680 Load & Validate",
                     class = "btn-primary btn-lg w-100"),
        tags$div(class = "text-center text-muted my-2",
                 style = "font-size:0.8rem;", "— or pick an example —"),
        selectInput("example_expr", "Example expression dataset",
                    choices = .expr_choices),
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
          card(full_screen = TRUE,
            card_header("Class Distribution"),
            plotOutput("plot_class_dist", height = "260px")
          ),
          card(full_screen = TRUE,
            card_header("Group Size Distribution"),
            plotOutput("plot_group_sizes", height = "260px")
          )
        ),
        card(full_screen = TRUE,
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
  nav_flow_sep(),
  nav_panel(
    value = "replicates",
    title = nav_step_title("2", "\U0001F501", "Replicates"),
    layout_sidebar(
      sidebar = sidebar(
        title = "\U2699\UFE0F Replicates — Settings",
        width = 340,
        tags$div(style = "display:flex; align-items:center; gap:4px;",
          tags$label("\U0001F4CB Metadata (CSV)"),
          help_tip("CSV with sample_id and replicate_group_id columns. Maps samples to their biological replicate groups for averaging.")
        ),
        fileInput("meta_file", NULL, accept = c(".csv", ".tsv", ".txt")),
        tags$hr(),
        tags$div(style = "display:flex; align-items:center; gap:4px;",
          tags$label("Sample ID column"),
          help_tip("Column name in the metadata that matches the row names in the expression data.")
        ),
        textInput("sample_id_col", NULL, value = "sample_id"),
        tags$div(style = "display:flex; align-items:center; gap:4px;",
          tags$label("Replicate group column"),
          help_tip("Column in metadata that defines which samples are replicates of each other.")
        ),
        textInput("replicate_group_col", NULL, value = "replicate_group_id"),
        tags$div(style = "display:flex; align-items:center; gap:4px;",
          tags$label("Aggregation method"),
          help_tip("How to combine replicate expression values. Mean is most common.")
        ),
        selectInput("collapse_method", NULL,
                    choices = c("mean", "median", "sum", "min", "max", "first")),
        tags$div(style = "display:flex; align-items:center; gap:4px;",
          tags$label("Conflict handling"),
          help_tip("What to do if replicates have different class labels. Error = stop, NA = set missing, concat = combine, first = keep first.")
        ),
        selectInput("meta_conflict", NULL,
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

  # ── Tab 3: Grouping (G) — run a grouping method and inspect it ───────────
  nav_flow_sep(),
  nav_panel(
    value = "grouping",
    title = nav_step_title("3", "\U0001F9E9", "Grouping"),
    layout_sidebar(
      sidebar = sidebar(
        width = 320, open = TRUE,
        tags$div(class = "sidebar-title", "\U0001F9E9 Grouping (G)"),
        tags$p(class = "text-muted", style = "font-size:0.83rem;",
               "Self-contained: choose the expression data here (or reuse what's loaded), run a grouping method, inspect it, then apply it to the pipeline."),

        # ── Step A: expression data (shared with the whole pipeline) ──────
        tags$div(class = "sidebar-title", style = "font-size:0.9rem;", "① Expression data"),
        uiOutput("g_data_status"),
        selectInput("g_example_expr", "Example dataset",
                    choices = c("— none —" = "", .expr_choices)),
        fileInput("g_expr_upload", "…or upload expression CSV",
                  accept = c(".csv", ".tsv", ".txt")),
        actionButton("g_load_data", "\U0001F4C2 Load data",
                     class = "btn-outline-light btn-sm w-100"),
        tags$hr(),

        tags$div(class = "text-muted", style = "font-size:0.83rem;",
                 "② Choose a method tab on the right, set its options, and click Run."),
        tags$hr(),

        # ── Step C: apply the grouping you ran to the pipeline ────────────
        tags$div(class = "sidebar-title", style = "font-size:0.9rem;", "③ Use it downstream"),
        actionButton("g_apply", "✅ Apply active tab to pipeline",
                     class = "btn-success w-100"),
        uiOutput("g_active_status")
      ),
      navset_card_tab(
        id = "g_method_tabs",

        # ── Method: Prior-knowledge grouping (from a file) ─────────────
        nav_panel(
          title = "\U0001F4DA Prior knowledge", value = "g_file",
          layout_columns(
            col_widths = c(4, 8),
            card(
              card_header("Inputs"),
              tags$p(class = "text-muted", style = "font-size:0.82rem;",
                     "Use a predefined grouping from biological knowledge — genes grouped by pathway (KEGG), function (GO), or disease association (DisGeNET). Pick a bundled example or upload your own."),
              tags$div(style = "display:flex; align-items:center; gap:4px;",
                tags$label("Grouping source", style = "margin:0;"),
                help_tip("Choose a bundled example knowledge base, or upload your own gene→group file (2-column CSV, or GMT).")
              ),
              selectInput("g_file_source", NULL,
                          choices = c("Example knowledge base" = "example",
                                      "Upload my own file" = "upload")),
              conditionalPanel(
                condition = "input.g_file_source == 'example'",
                selectInput("g_file_example", "Example grouping",
                            choices = .groups_choices, selected = .default_groups)
              ),
              conditionalPanel(
                condition = "input.g_file_source == 'upload'",
                fileInput("g_file_upload", "Grouping file",
                          accept = c(".csv", ".tsv", ".gmt", ".txt")),
                selectInput("g_file_format", "Format",
                            choices = c("Auto-detect" = "auto",
                                        "2-column CSV" = "csv", "GMT" = "gmt"))
              ),
              actionButton("g_file_run", "\U0001F680 Load grouping",
                           class = "btn-primary w-100")
            ),
            tags$div(
              gsm_help_callout(
                "What is prior-knowledge grouping?",
                tagList(
                  "A curated file that maps each gene to one or more biological ",
                  tags$b("groups"), " — e.g. a KEGG pathway, a GO term, or a DisGeNET disease. ",
                  "It's ", tags$b("prior knowledge"), ": the groups come from biology, not from your data. ",
                  "Later, Scoring tests each group and Modeling keeps the most predictive ones.",
                  tags$br(), tags$br(),
                  tags$b("Steps: "), tags$b("①"), " load expression data (sidebar). ",
                  tags$b("②"), " pick an example knowledge base (or upload). ",
                  tags$b("③"), " click Load — the table shows each group, its size, and ",
                  tags$b("coverage"), " (how many of its genes are actually in your data). ",
                  tags$b("④"), " Apply to pipeline (sidebar) to score these groups in Pre-Score/GSM.")),
              gsm_help_callout(
                "Reading the results",
                tagList(
                  tags$b("Groups / Unique genes / Median size"), " — how many groups and genes the file defines. ",
                  tags$b("Gene coverage"), " — the share of your expression genes that appear in the grouping; low coverage means many groups reference genes you don't measure. ",
                  "The ", tags$b("Groups table"), " lists each group with its gene count and members; the ",
                  tags$b("histogram"), " shows the spread of group sizes.")),
              uiOutput("g_file_valueboxes"),
              card(full_screen = TRUE, card_header("Groups"),
                   DTOutput("g_file_groups_dt")),
              card(full_screen = TRUE, card_header("Group-size distribution"),
                   plotOutput("g_file_size_hist", height = "240px"))
            )
          )
        ),

        # ── Method: BRIDGE ─────────────────────────────────────────────
        nav_panel(
          title = "\U0001F309 BRIDGE", value = "g_bridge",
          layout_columns(
            col_widths = c(4, 8),
            card(
              card_header("Inputs"),
              tags$p(class = "text-muted", style = "font-size:0.82rem;",
                     "Data-driven grouping from expression + labels, regularized by a prior graph. This runs on ALL samples (exploratory); applied to the pipeline it is re-derived per training fold."),
              selectInput("g_bridge_prior_source", "Prior graph source",
                          choices = c("Example grouping file" = "example",
                                      "Upload edge list / grouping" = "upload",
                                      "Grouping loaded on Data tab" = "loaded",
                                      "None (data-driven)" = "none")),
              conditionalPanel(
                condition = "input.g_bridge_prior_source == 'example'",
                selectInput("g_bridge_prior_example", "Example prior",
                            choices = .groups_choices, selected = .default_groups)
              ),
              conditionalPanel(
                condition = "input.g_bridge_prior_source == 'upload'",
                fileInput("g_bridge_prior_file", "Prior file",
                          accept = c(".csv", ".tsv", ".gmt", ".txt"))
              ),
              selectInput("g_bridge_data_type", "Data type",
                          choices = c("Gaussian" = "gaussian",
                                      "Nonparanormal" = "nonparanormal",
                                      "Compositional (coming soon)" = "compositional")),
              selectInput("g_bridge_prefilter", "Speed pre-filter (label-blind)",
                          choices = c("None" = "none",
                                      "Top-variance genes" = "variance",
                                      "Prior-graph genes only" = "prior")),
              layout_columns(col_widths = c(6, 6),
                selectInput("g_bridge_precision", "Estimator",
                            choices = c("glasso (prior-weighted, \u2264250 genes)" = "glasso",
                                        "jgl (fused joint estimation)" = "jgl",
                                        "mb / fast (large p)" = "mb")),
                numericInput("g_bridge_max_genes", "Max genes", value = 1000, min = 10, step = 50)
              ),
              tags$details(
                tags$summary(class = "text-muted", style = "cursor:pointer; font-size:0.85rem;",
                             "⚙️ Advanced hyperparameters"),
                tags$div(style = "margin-top:6px;",
                  layout_columns(col_widths = c(6, 6),
                    numericInput("g_bridge_lambda", "lambda", value = 0.1, min = 0, step = 0.01),
                    numericInput("g_bridge_beta", "beta", value = 0.5, min = 0, max = 0.99, step = 0.05)
                  ),
                  layout_columns(col_widths = c(6, 6),
                    numericInput("g_bridge_lambda2", "lambda2", value = 0.05, min = 0, step = 0.01),
                    numericInput("g_bridge_lambda1", "lambda1", value = 0.25, min = 0, step = 0.05)
                  ),
                  textInput("g_bridge_resolution", "Resolution grid", value = "0.5, 1, 1.5, 2"),
                  sliderInput("g_bridge_tau", "tau (stability)", min = 0, max = 1, value = 0.6, step = 0.05),
                  layout_columns(col_widths = c(6, 6),
                    numericInput("g_bridge_min_size", "Min size", value = 5, min = 1),
                    numericInput("g_bridge_max_size", "Max size", value = 500, min = 10)
                  )
                )
              ),
              conditionalPanel(
                condition = "!output.g_bridge_running",
                actionButton("g_bridge_run", "\U0001F680 Run BRIDGE grouping",
                             class = "btn-primary w-100")
              ),
              conditionalPanel(
                condition = "output.g_bridge_running",
                tags$div(class = "text-muted text-center",
                         style = "font-size:0.82rem;", "Running… use Stop on the right →")
              )
            ),
            tags$div(
              gsm_help_callout(
                "What is BRIDGE?",
                tagList(
                  tags$b("BRIDGE discovers the groups from YOUR data"), " instead of taking them from a file. ",
                  "It learns how genes co-vary in each class, finds where those relationships ", tags$b("differ"),
                  " between the two classes (e.g. tumor vs normal), and turns those differential sub-networks into ", tags$b("modules"), " (groups). ",
                  "A ", tags$b("prior graph"), " (an optional pathway/disease file) nudges biologically related genes together — hence ",
                  tags$em("Biologically-Regularized Interaction-Differential Grouping Engine"), ".",
                  tags$br(), tags$br(),
                  tags$b("Use it when"), " you want data-driven, condition-specific groups rather than fixed pathway definitions. ",
                  "It runs here on ", tags$b("all samples"), " for exploration; once applied, GSM re-derives it ", tags$b("per training fold"), " so the evaluation stays honest (no leakage).")),
              gsm_help_callout(
                "How to run it",
                tagList(
                  tags$b("①"), " load expression data (sidebar). ",
                  tags$b("②"), " pick a ", tags$b("prior source"), " (example pathway/disease file, upload, or 'None'). ",
                  tags$b("③"), " keep ", tags$b("Speed pre-filter"), " = Top-variance and ", tags$b("Estimator"), " = mb / fast so it finishes in seconds. ",
                  tags$b("④"), " click ", tags$b("Run"), " — watch the live progress; ", tags$b("Stop"), " anytime. ",
                  tags$b("⑤"), " inspect, then ", tags$b("Apply to pipeline"), " (sidebar ③).")),
              gsm_help_callout(
                "Reading the results",
                tagList(
                  tags$b("Modules table — J = D + λ₁·C − λ₂·|M|"), " ranks each module by an overall score, where:",
                  tags$ul(style = "margin:4px 0 0 0;",
                    tags$li(tags$b("D"), " = discrimination — how well the module separates the two classes (higher = more predictive)."),
                    tags$li(tags$b("C"), " = prior cohesion — how much the module agrees with the prior graph (biological plausibility)."),
                    tags$li(tags$b("stability"), " = how often the module survives when the samples are resampled (1.0 = very robust)."),
                    tags$li(tags$b("size"), " = number of genes; |M| is penalized so huge modules don't win by size alone.")),
                  tags$b("Diagnostics"), " show how many genes/candidates went in and which backend ran. The ",
                  tags$b("network"), " draws the differential connections, coloured by module.")),
              uiOutput("g_bridge_progress"),
              uiOutput("g_bridge_valueboxes"),
              uiOutput("g_bridge_diag_card"),
              card(full_screen = TRUE, card_header("Modules (J = D + λ₁C − λ₂|M|)"),
                   DTOutput("g_bridge_modules_dt")),
              card(full_screen = TRUE, card_header("Module-size distribution"),
                   plotOutput("g_bridge_size_hist", height = "240px")),
              card(full_screen = TRUE,
                   card_header("Differential-network modules"),
                   visNetwork::visNetworkOutput("g_bridge_network", height = "420px"))
            )
          )
        )
      )
    )
  ),

  # ── Tab 4: Pre-Scoring ──────────────────────────────────────────────────
  nav_flow_sep(),
  nav_panel(
    value = "prescore",
    title = nav_step_title("4", "\U0001F4CA", "Pre-Score"),
    layout_sidebar(
      sidebar = sidebar(
        title = "\U2699\UFE0F Pre-Score — Settings",
        width = 340,
        tags$div(style = "display:flex; align-items:center; gap:4px;",
          tags$label("Method"),
          help_tip("Statistical test to rank groups by differential expression. ANOVA/limma combine per-gene p-values with Fisher's method (biased toward large groups). Fisher (chi-square) calibrates for group size; EBM also corrects for inter-gene correlation; CAMERA is a competitive, correlation-aware gene-set test (two-class only). GSVA scores each group's single-sample enrichment activity (the signal the ML model sees). None skips pre-scoring.")
        ),
        selectInput("prescore_method", NULL,
                    choices = c("ANOVA" = "anova",
                                "limma (moderated)" = "limma",
                                "Fisher (chi-square, size-calibrated)" = "fisher_chisq",
                                "Empirical Brown (correlation-aware)" = "ebm",
                                "CAMERA (competitive, 2-class)" = "camera",
                                "GSVA / ssGSEA enrichment" = "gsva",
                                "None" = "none")),
        tags$div(style = "display:flex; align-items:center; gap:4px;",
          tags$label("Top % groups to keep"),
          help_tip("Percentage of groups passing the pre-score filter. Lower = more selective, faster scoring. 20% is a good default.")
        ),
        sliderInput("top_percent", NULL,
                    min = 5, max = 100, value = 20, step = 5),
        tags$div(style = "display:flex; align-items:center; gap:4px;",
          tags$label("Training ratio"),
          help_tip("Fraction of samples used for training. The rest is held out for validation. 0.7 means 70% train, 30% test.")
        ),
        sliderInput("train_ratio", NULL,
                    min = 0.5, max = 0.9, value = 0.7, step = 0.05),
        tags$div(style = "display:flex; align-items:center; gap:4px;",
          tags$label("Random seed"),
          help_tip("Seed for reproducibility. Same seed = same train/test split.")
        ),
        numericInput("prescore_seed", NULL, value = 44),
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
                     class = "btn-primary btn-lg w-100"),
        tags$hr(),
        tags$strong("\U0001F3B2 Group significance (permutation)"),
        help_tip("Tests each group against random same-size gene sets. Reports empirical p-value + BH FDR."),
        selectInput("perm_method", "Classifier",
                    choices = c("Random Forest" = "rf",
                                "Support Vector Machine" = "svm",
                                "K-Nearest Neighbors" = "knn",
                                "XGBoost" = "xgb")),
        numericInput("perm_nperm", "Permutations per group",
                     value = 200, min = 20, max = 2000, step = 20),
        actionButton("btn_permtest", "\U0001F3B2 Test group significance",
                     class = "btn-outline-light w-100")
      ),
      # ── Live view while pre-scoring is running ─────────────────────────
      conditionalPanel(
        condition = "output.prescore_running",
        uiOutput("prescore_live_status"),
        layout_columns(
          col_widths = c(7, 5),
          card(full_screen = TRUE,
            card_header("Live: top groups so far"),
            plotOutput("plot_prescore_live", height = "340px")
          ),
          card(full_screen = TRUE,
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
            card(full_screen = TRUE,
              card_header("Top Groups by Pre-Score (Fisher statistic)"),
              plotOutput("plot_prescore_bars", height = "320px")
            ),
            card(full_screen = TRUE,
              card_header("Score Distribution (all groups)"),
              plotOutput("plot_prescore_hist", height = "320px")
            )
          )
        ),
        card(full_screen = TRUE,
          card_header(textOutput("prescore_table_header", inline = TRUE)),
          DTOutput("prescore_table")
        )
      ),
      conditionalPanel(
        condition = "output.permtest_done",
        card(full_screen = TRUE,
          card_header("\U0001F3B2 Empirical group significance (permutation p-value + BH FDR)"),
          tags$small(class = "text-muted",
                     "observed = train-only mean F1; null_mean = mean over random same-size gene sets; p_adj = BH-FDR across tested groups."),
          DTOutput("permtest_table")
        )
      ),
      conditionalPanel(
        condition = "!output.prescore_done && !output.prescore_running && !output.data_loaded",
        tags$div(
          style = "text-align:center; padding:80px 20px; opacity:0.5;",
          tags$h3("\U0001F4CA Pre-score groups"),
          tags$p("Load data first, then run pre-scoring to filter groups before ML scoring.")
        )
      ),
      conditionalPanel(
        condition = "!output.prescore_done && !output.prescore_running && output.data_loaded",
        tags$div(
          style = "text-align:center; padding:80px 20px; opacity:0.5;",
          tags$h3("\U0001F4CA Ready to pre-score"),
          tags$p("Data loaded. Click \U0001F4CA Run Pre-Scoring to filter groups before ML scoring.")
        )
      )
    )
  ),

  # ── Tab 4: Scoring ─────────────────────────────────────────────────────
  nav_flow_sep(),
  nav_panel(
    value = "score",
    title = nav_step_title("5", "\U0001F3AF", "GSM"),
    layout_sidebar(
      sidebar = sidebar(
        title = "\U2699\UFE0F GSM — Settings",
        width = 340,

        # Collapsible setting cards — click a header to open/close it.
        accordion(
          id = "gsm_settings_acc", multiple = TRUE, open = "Grouping",

          # ── Grouping approach (Prior knowledge / BRIDGE) ──────────────
          accordion_panel(
            "Grouping", icon = tags$span("\U0001F9E9"),
            tags$p(class = "text-muted", style = "font-size:0.8rem;",
                   "How genes are grouped. Configure and run it here, or explore in the Grouping tab."),
            radioButtons("gsm_grp_method", NULL,
                         choices = c("\U0001F4DA Prior knowledge (pathways/GO/disease)" = "file",
                                     "\U0001F309 BRIDGE (data-driven + prior)" = "bridge")),
            conditionalPanel(
              condition = "input.gsm_grp_method == 'file'",
              selectInput("gsm_grp_file_source", "Source",
                          choices = c("Example knowledge base" = "example",
                                      "Upload file" = "upload")),
              conditionalPanel("input.gsm_grp_file_source == 'example'",
                selectInput("gsm_grp_file_example", "Example grouping",
                            choices = .groups_choices, selected = .default_groups)),
              conditionalPanel("input.gsm_grp_file_source == 'upload'",
                fileInput("gsm_grp_file_upload", "Grouping file (CSV / GMT)",
                          accept = c(".csv", ".tsv", ".gmt", ".txt")))
            ),
            conditionalPanel(
              condition = "input.gsm_grp_method == 'bridge'",
              selectInput("gsm_grp_bridge_prior", "Prior graph",
                          choices = c("Example grouping file" = "example",
                                      "Upload" = "upload", "None (data-driven)" = "none")),
              conditionalPanel("input.gsm_grp_bridge_prior == 'example'",
                selectInput("gsm_grp_bridge_prior_example", "Example prior",
                            choices = .groups_choices, selected = .default_groups)),
              conditionalPanel("input.gsm_grp_bridge_prior == 'upload'",
                fileInput("gsm_grp_bridge_prior_file", "Prior file",
                          accept = c(".csv", ".tsv", ".gmt", ".txt"))),
              selectInput("gsm_grp_bridge_precision", "Estimator",
                          choices = c("mb / fast" = "mb",
                                      "glasso (prior-weighted, \u2264250 genes)" = "glasso",
                                      "jgl (fused joint)" = "jgl")),
              selectInput("gsm_grp_bridge_prefilter", "Speed pre-filter",
                          choices = c("Top-variance" = "variance",
                                      "Prior genes" = "prior", "None" = "none")),
              numericInput("gsm_grp_bridge_max_genes", "Max genes", value = 1000, min = 10, step = 100),
              tags$small(class = "text-muted d-block",
                         "Re-derived per training fold (leakage-free). mb + a pre-filter keeps folds fast.")
            ),
            actionButton("gsm_grp_run", "✅ Use this grouping", class = "btn-success w-100"),
            uiOutput("gsm_grouping_note")
          ),

          accordion_panel(
            "Classifier", icon = tags$span("\U0001F916"),
            selectInput("score_method", "Group scorer",
                        choices = CLASSIFIER_CHOICES),
            conditionalPanel("input.score_method == 'rf'",
              numericInput("rf_ntree", "Number of Trees (ntree)", value = 500, min = 1, step = 50)),
            conditionalPanel("input.score_method == 'svm'",
              selectInput("svm_kernel", "SVM Kernel",
                          choices = c("Linear" = "linear", "Radial" = "radial",
                                      "Polynomial" = "polynomial", "Sigmoid" = "sigmoid"),
                          selected = "linear"),
              numericInput("svm_cost", "Cost (C)", value = 1, min = 0.001, step = 0.5),
              numericInput("svm_gamma", "Gamma (non-linear only)", value = 0.1, min = 0.0001, step = 0.01)),
            conditionalPanel("input.score_method == 'knn'",
              numericInput("knn_k", "Neighbors (k)", value = 3, min = 1, max = 20, step = 1)),

            # ── Deep learning (torch): MLP and GNN ──────────────────────────
            conditionalPanel(
              "input.score_method == 'mlp' || input.score_method == 'gnn' ||
               input.model_method == 'mlp' || input.model_method == 'gnn'",
              tags$div(class = "text-muted small", style = "margin-bottom:6px;",
                       "Needs the 'torch' package. These models have far more ",
                       "capacity than typical transcriptomic sample sizes — ",
                       "treat Random Forest / XGBoost as the baseline to beat."),
              conditionalPanel("input.score_method == 'mlp' || input.model_method == 'mlp'",
                tags$div(style = "display:flex; align-items:center; gap:4px;",
                  textInput("mlp_hidden", "Hidden layers (widths)", value = "64,32"),
                  help_tip("Comma-separated layer widths. '64,32' = two hidden layers of 64 then 32 units."))),
              conditionalPanel("input.score_method == 'gnn' || input.model_method == 'gnn'",
                numericInput("gnn_hidden", "GCN hidden width", value = 32, min = 4, max = 256, step = 4),
                tags$div(style = "display:flex; align-items:center; gap:4px;",
                  numericInput("gnn_cor_threshold", "Fallback graph |cor| cutoff",
                               value = 0.3, min = 0, max = 0.99, step = 0.05),
                  help_tip("The GNN convolves over BRIDGE's differential network when a BRIDGE grouping is active. With a file-based grouping (KEGG/GO/DisGeNET) there is no such network, so it builds a gene-gene graph from training-sample correlation and keeps edges above this cutoff."))),
              numericInput("nn_epochs", "Training epochs", value = 100, min = 10, max = 2000, step = 10),
              numericInput("nn_lr", "Learning rate", value = 0.01, min = 0.0001, max = 1, step = 0.005),
              numericInput("nn_dropout", "Dropout", value = 0.3, min = 0, max = 0.9, step = 0.05)
            )
          ),

          accordion_panel(
            "Gene shaving", icon = tags$span("\U2702\UFE0F"),
            tags$div(style = "display:flex; align-items:center; gap:4px;",
              checkboxInput("shave_genes", "Enable gene shaving", value = TRUE),
              help_tip("Rank genes by a per-gene Welch t-test on TRAINING samples only (no leakage), then keep the significant ones.")),
            conditionalPanel("input.shave_genes",
              selectInput("shave_method", "Selection rule",
                          choices = c("BH-FDR threshold" = "fdr", "Top-K genes" = "topk")),
              conditionalPanel("input.shave_method == 'fdr'",
                numericInput("fdr_alpha", "FDR threshold (alpha)", value = 0.05, min = 0.001, max = 0.5, step = 0.01)),
              conditionalPanel("input.shave_method == 'topk'",
                numericInput("shave_k", "Top-K genes to keep", value = 100, min = 2, step = 10))),
            checkboxInput("normalize_z", "Z-score normalize (per-fold, train-fit)", value = FALSE)
          ),

          accordion_panel(
            "Scoring scope & depth", icon = tags$span("\U0001F3AF"),
            radioButtons("score_scope", "Groups to score",
                         choices = c("Top pre-scored (recommended)" = "selected",
                                     "All groups" = "all"), selected = "selected"),
            numericInput("n_top_groups", "Accumulative depth (K)", value = 10, min = 1),
            numericInput("score_seed", "Random seed", value = 44)
          ),

          accordion_panel(
            "Cross-validation", icon = tags$span("\U0001F501"),
            checkboxInput("run_mccv_mode", "Multi-Fold MCCV", value = FALSE),
            conditionalPanel("input.run_mccv_mode",
              numericInput("mccv_iterations", "MCCV Iterations", value = 10, min = 1, max = 500),
              sliderInput("mccv_cores", "CPU Cores", min = 1, max = parallel::detectCores(),
                          value = max(1, parallel::detectCores() - 1), step = 1),
              selectInput("model_method", "Final Modeling Classifier",
                          choices = c("Same as Scoring" = "same", CLASSIFIER_CHOICES),
                          selected = "same"),
              sliderInput("mccv_train_ratio", "Training Ratio", min = 0.5, max = 0.9, value = 0.7, step = 0.05))
          )
        ),

        tags$hr(),
        actionButton("btn_score", "\U0001F680 Run Analysis",
                     class = "btn-primary btn-lg w-100")
      ),

      # ── Always-visible status alert (readiness / what's happening) ──────
      uiOutput("score_status_alert"),

      # ── Run configuration: what THIS run used (grouping, shaving, …) ────
      uiOutput("gsm_run_config"),

      # ── Live view while Parallel MCCV is running ───────────────────────
      conditionalPanel(
        condition = "output.mccv_running",
        tags$div(class = "d-flex justify-content-end mb-2",
          actionButton("btn_stop_score", "⏹️ Stop run",
                       class = "btn-danger btn-sm")
        ),
        uiOutput("mccv_live_status"),
        card(
          card_header("Live: per-fold progress (updates every second)"),
          uiOutput("mccv_fold_progress")
        ),
        layout_columns(
          col_widths = c(6, 6),
          card(
            card_header("Live: fold accuracy as folds complete"),
            plotOutput("plot_mccv_live", height = "260px")
          ),
          card(
            card_header("Live: running metric averages"),
            uiOutput("mccv_live_means")
          )
        ),
        card(full_screen = TRUE,
          card_header("Groups created — most recent completed fold"),
          DTOutput("gsm_live_groups")
        )
      ),

      # Main Panel Results Conditional Displays
      conditionalPanel(
        condition = "output.score_done && !output.mccv_running",

        # ── Professional tabbed results dashboard ─────────────────────────
        navset_card_tab(
          id = "score_tabs",

          # ── Overview: run summary + headline performance ───────────────
          nav_panel(
            title = "\U0001F4CA Overview",
            uiOutput("score_summary"),
            conditionalPanel(
              condition = "!input.run_mccv_mode",
              card(full_screen = TRUE,
                card_header("Best-k joint model performance"),
                layout_columns(col_widths = c(4, 4, 4),
                  uiOutput("model_vbox_acc"), uiOutput("model_vbox_kappa"),
                  uiOutput("model_vbox_sens")),
                layout_columns(col_widths = c(4, 4, 4),
                  uiOutput("model_vbox_spec"), uiOutput("model_vbox_prec"),
                  uiOutput("model_vbox_f1"))
              )
            ),
            conditionalPanel(
              condition = "input.run_mccv_mode",
              card(full_screen = TRUE,
                card_header("Aggregate performance across folds (mean ± SD, 95% CI)"),
                DTOutput("mccv_summary_table")
              )
            )
          ),

          # ── Performance: accumulative metric-vs-k curve ────────────────
          nav_panel(
            title = "\U0001F4C8 Performance",
            card(full_screen = TRUE,
              card_header("Accumulative performance: metric vs. number of top groups (k)"),
              layout_sidebar(
                sidebar = sidebar(
                  width = 220, position = "right",
                  selectInput("acc_metric", "Metric",
                              choices = c("Accuracy", "AUC", "Sensitivity",
                                          "Specificity", "Precision", "F_measure",
                                          "Cohens_kappa"),
                              selected = "Accuracy"),
                  tags$small(class = "text-muted",
                             "MCCV shows mean ± SD across folds; the dashed line marks best k.")
                ),
                plotOutput("plot_accumulative", height = "340px")
              )
            ),
            conditionalPanel(
              condition = "input.run_mccv_mode",
              card(full_screen = TRUE,
                card_header("Per-fold metric distribution"),
                plotOutput("plot_mccv_folds", height = "320px")
              )
            )
          ),

          # ── Groups: scores, frequency ranking, and RRA consensus ───────
          nav_panel(
            title = "\U0001F5C2️ Groups",
            conditionalPanel(
              condition = "!input.run_mccv_mode",
              layout_columns(col_widths = c(7, 5),
                card(full_screen = TRUE, card_header("Group scores (train-only F1, ranked)"),
                     DTOutput("score_table")),
                card(full_screen = TRUE, card_header("Score distribution"),
                     plotOutput("plot_scores", height = "320px")))
            ),
            conditionalPanel(
              condition = "input.run_mccv_mode",
              layout_columns(col_widths = c(4, 4, 4),
                card(full_screen = TRUE, card_header("Average group accuracy (top stable)"),
                     plotOutput("plot_mccv_groups", height = "320px")),
                card(full_screen = TRUE, card_header("\U0001F3C5 RRA — groups over iterations"),
                     DTOutput("rra_group_table")),
                card(full_screen = TRUE, card_header("Group stability & accuracies"),
                     DTOutput("mccv_stability_table")))
            )
          ),

          # ── Genes: importance and RRA consensus across iterations ──────
          nav_panel(
            title = "\U0001F9EC Genes",
            conditionalPanel(
              condition = "input.run_mccv_mode",
              layout_columns(col_widths = c(4, 4, 4),
                card(full_screen = TRUE, card_header("\U0001F9EA Consensus biomarkers"),
                     tags$small(class = "text-muted",
                                "Selection freq × mean importance + RRA"),
                     DTOutput("biomarker_table")),
                card(full_screen = TRUE, card_header("\U0001F3C5 RRA — genes (importance)"),
                     DTOutput("rra_feature_table")),
                card(full_screen = TRUE, card_header("Ranked features"),
                     DTOutput("ranked_features_table")))
            ),
            conditionalPanel(
              condition = "!input.run_mccv_mode",
              card(full_screen = TRUE, card_header("Genes in the best-k joint model"),
                   DTOutput("single_feature_table"))
            )
          ),

          # ── Discover: de novo robust gene-module discovery (G-S-M detect) ──
          nav_panel(
            title = "\U0001F30B Discover",
            conditionalPanel(
              condition = "input.run_mccv_mode",
              layout_sidebar(
                sidebar = sidebar(
                  title = "Discovery Settings",
                  width = 280, open = TRUE,
                  tags$div(style = "display:flex; align-items:center; gap:4px;",
                    tags$label("Direction"),
                    help_tip("Which direction of differential expression to consider. Absolute = both up and down regulated genes.")
                  ),
                  radioButtons("detect_mode", NULL,
                               c("Absolute" = "abs", "Up" = "up", "Down" = "down"),
                               selected = "abs", inline = TRUE),
                  tags$div(style = "display:flex; align-items:center; gap:4px;",
                    tags$label("Top genes"),
                    help_tip("Number of most differentially expressed genes to include in the stability network.")
                  ),
                  numericInput("detect_top_genes", NULL, 200, 20, 2000, 10),
                  tags$div(style = "display:flex; align-items:center; gap:4px;",
                    tags$label("Min frequency"),
                    help_tip("Minimum fraction of MCCV folds a gene must appear in to be considered robust. Higher = more stringent.")
                  ),
                  numericInput("detect_min_freq", NULL, 0.5, 0, 1, 0.05),
                  layout_columns(
                    col_widths = c(6, 6),
                    tags$div(
                      tags$div(style = "display:flex; align-items:center; gap:4px;",
                        tags$label("Resolution"),
                        help_tip("Louvain clustering resolution. Higher = more, smaller modules.")
                      ),
                      numericInput("detect_resolution", NULL, 1, 0.1, 5, 0.1)
                    ),
                    tags$div(
                      tags$div(style = "display:flex; align-items:center; gap:4px;",
                        tags$label("Bootstraps"),
                        help_tip("Number of bootstrap iterations for network stability estimation. 0 = skip bootstrap.")
                      ),
                      numericInput("detect_nboot", NULL, 50, 0, 200, 10)
                    )
                  ),
                  actionButton("run_detect", "\U0001F30B Discover",
                               class = "btn-primary btn-lg w-100")
                ),
                uiOutput("detect_status"),
                layout_columns(
                  col_widths = c(6, 6),
                  card(full_screen = TRUE,
                    card_header("Robust genes (RRA + stability)"),
                    DTOutput("detect_robust_table")),
                  card(full_screen = TRUE,
                    card_header("Gene–gene stability network"),
                    visNetwork::visNetworkOutput("detect_network_vis", height = "360px"))
                ),
                layout_columns(
                  col_widths = c(6, 6),
                  card(full_screen = TRUE,
                    card_header("Discovered modules vs. known groups"),
                    DTOutput("detect_modules_table")),
                  card(full_screen = TRUE,
                    card_header("Module detail"),
                    selectInput("detect_module_pick", "Module", choices = NULL),
                    DTOutput("detect_module_detail"))
                )
              )
            ),
            conditionalPanel(
              condition = "!input.run_mccv_mode",
              tags$div(style = "padding:40px; opacity:0.6;",
                       "De novo discovery aggregates many iterations — run ",
                       tags$b("Parallel MCCV"), " first, then come back here.")
            )
          ),

          # ── Iterations: per-iteration detail (MCCV) ────────────────────
          nav_panel(
            title = "\U0001F501 Iterations",
            conditionalPanel(
              condition = "input.run_mccv_mode",
              card(full_screen = TRUE, card_header("Best-k model per iteration"),
                   DTOutput("iter_metrics_table")),
              card(full_screen = TRUE, card_header("Accumulative performance — every iteration × k"),
                   DTOutput("iter_accum_table"))
            ),
            conditionalPanel(
              condition = "!input.run_mccv_mode",
              card(full_screen = TRUE, card_header("Accumulative performance (this split)"),
                   DTOutput("iter_accum_table_single"))
            )
          ),

          # ── Export ─────────────────────────────────────────────────────
          nav_panel(
            title = "\U0001F4BE Export",
            card(
              uiOutput("export_note"),
              tags$div(
                class = "d-flex gap-2 flex-wrap",
                downloadButton("dl_xlsx", "Download Excel workbook (.xlsx)",
                               class = "btn-primary"),
                downloadButton("dl_csv", "Download CSV bundle (.zip)",
                               class = "btn-outline-light"),
                actionButton("save_disk", "\U0001F4C1 Save to workspace folder",
                             class = "btn-outline-light")
              ),
              tags$small(class = "text-muted d-block mt-2",
                         "Tip: downloads don't work in VS Code's Simple Browser — use ",
                         tags$b("Save to workspace folder"), ", or open the app in a real browser tab."),
              uiOutput("save_disk_path")
            ),
            conditionalPanel(
              condition = "output.has_timings",
              card(full_screen = TRUE,
                card_header("\U23F1\UFE0F Recorded Step Runtimes (Exported to 'Timings' Excel sheet)"),
                DTOutput("export_timings_table")
              )
            )
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

  # ── Tab 5: File ────────────────────────────────────────────────────
  # ── Utilities (right group): tools/modes, not sequential pipeline steps ──
  nav_spacer(),
  nav_item(tags$span(class = "nav-vrule")),
  nav_panel(
    value = "file_menu",
    title = nav_util_title("💾", "Session"),
    card(
      card_header("📝 Import & Export Parameters"),
      tags$p("Save your current analysis parameters to a YAML file to share or reuse later, or load an existing configuration."),
      layout_columns(
        col_widths = c(6, 6),
        card(
          card_header("Import Parameters"),
          fileInput("import_settings", "Upload YAML Settings", accept = c(".yaml", ".yml")),
          uiOutput("import_settings_status")
        ),
        card(
          card_header("Export Parameters"),
          downloadButton("export_settings", "Download YAML Settings", class = "btn-primary w-100")
        )
      )
    ),
    card(
      card_header("\U0001F4BE Save & Load Session"),
      tags$p(class = "text-muted",
             "Save the entire session (data, results, models, discovery) to a single file. ",
             "Reload it later to restore all results without recalculating."),
      layout_columns(
        col_widths = c(6, 6),
        card(
          card_header("Load Session"),
          fileInput("load_session_file", "Upload .rds Session File",
                    accept = c(".rds", ".RDS")),
          uiOutput("load_session_status")
        ),
        card(
          card_header("Save Session"),
          tags$div(style = "display:flex; gap:6px; align-items:end;",
            tags$div(style = "flex:1;",
              textInput("save_session_path", "Save to:",
                        value = file.path(path.expand("~"), "Downloads",
                                          paste0("gsmR_session_", format(Sys.time(), "%Y%m%d"), ".rds")))
            ),
            actionButton("btn_browse_save", "\U0001F4C2",
                         class = "btn-outline-light", style = "margin-bottom:16px;")
          ),
          actionButton("btn_save_session", "\U0001F4BE Save Session",
                       class = "btn-primary w-100"),
          uiOutput("save_session_status"),
          tags$small(class = "text-muted d-block mt-2",
                     "Includes: expression data, groups, dataset, pre-score results, ",
                     "scoring/MCCV results, models, discovery, and all settings.")
        )
      )
    )
  ),

  # ── Tab 6: Batch ───────────────────────────────────────────────────
  nav_panel(
    value = "batch",
    title = nav_util_title("\U0001F9EE", "Batch"),
    layout_sidebar(
      sidebar = sidebar(
        title = "\U2699\UFE0F Batch — Setup", width = 360,
        tags$div(style = "display:flex; align-items:center; gap:4px;",
          tags$label("\U0001F9EC Datasets"),
          help_tip("Pick built-in examples and/or upload your own expression files. Uploads and selected examples are combined. Each dataset runs against the single grouping file below, sharing the parameters here; all (dataset x fold) work units run over one shared core pool.")
        ),
        selectInput("batch_datasets", "Built-in examples",
                    choices = .example_inputs, multiple = TRUE),
        tags$div(style = "display:flex; gap:6px; margin-top:-6px; margin-bottom:8px;",
          actionButton("batch_select_all", "Select all",
                       class = "btn-outline-light btn-sm", style = "flex:1;"),
          actionButton("batch_clear_all", "Clear all",
                       class = "btn-outline-light btn-sm", style = "flex:1;")),
        fileInput("batch_upload", "…or upload expression files",
                  multiple = TRUE, accept = c(".csv", ".tsv", ".txt")),
        tags$hr(),
        tags$div(style = "display:flex; align-items:center; gap:4px;",
          tags$label("\U0001F5C2️ Grouping file"),
          help_tip("One grouping file shared by every dataset. Pick an example, or upload your own (2-column gene,group CSV or GMT). An upload overrides the dropdown.")),
        selectInput("batch_groups", "Built-in example", choices = .groups_choices,
                    selected = .default_groups),
        fileInput("batch_groups_upload", "…or upload grouping file",
                  accept = c(".csv", ".tsv", ".gmt", ".txt")),
        tags$hr(),
        field_label("Classifier", "Machine-learning algorithm used to score each group and build models. Random Forest is robust and gives gene importances; SVM suits small samples; KNN is simple; XGBoost handles complex patterns. The Neural Network (MLP) and Graph Neural Network (GCN) need the 'torch' package — the GCN convolves over BRIDGE's differential network, or a correlation graph with file-based groupings. Applied to every dataset."),
        selectInput("bp_score_method", NULL, choices = CLASSIFIER_CHOICES),
        field_label("Pre-scoring", "Fast statistical test that ranks groups so only the strongest go to ML scoring. Computed on the TRAINING split of each fold (no leakage). ANOVA/limma combine per-gene p-values; Fisher (chi-square) calibrates for group size; Empirical Brown also corrects gene-gene correlation; CAMERA is a competitive 2-class test; GSVA scores single-sample enrichment; None keeps all groups."),
        selectInput("bp_prescore_method", NULL,
                    choices = c("ANOVA" = "anova", "limma (moderated)" = "limma",
                                "Fisher (chi-square)" = "fisher_chisq",
                                "Empirical Brown" = "ebm",
                                "CAMERA (2-class)" = "camera",
                                "GSVA / ssGSEA" = "gsva", "None" = "none")),
        layout_columns(col_widths = c(6, 6),
          tags$div(
            field_label("MCCV iters", "Number of Monte-Carlo cross-validation folds per dataset. Each fold makes a fresh random train/test split and runs the full score->rank->model pipeline; results are averaged. More folds = more stable estimates but longer runtime."),
            numericInput("bp_mccv_iterations", NULL, value = 10, min = 1, max = 500)),
          tags$div(
            field_label("Pre-score top %", "Applied on the TRAINING part only, independently in each MCCV round (fold): groups are ranked by the pre-scoring test on that fold's training split, and only this top fraction is kept before the expensive ML scoring (the held-out test samples are never used). E.g. 20% of 327 groups keeps ~66. Lower = faster, more aggressive filtering; 100% keeps all groups."),
            numericInput("bp_top_percent", NULL, value = 100, min = 1, max = 100))),
        layout_columns(col_widths = c(6, 6),
          tags$div(
            field_label("Train ratio", "Fraction of samples used for training in each fold's random split (e.g. 0.7 = 70% train, 30% held-out test). Pre-scoring, gene shaving, scoring and modeling all use the training part only."),
            numericInput("bp_train_ratio", NULL, value = 0.7, min = 0.5, max = 0.95, step = 0.05)),
          tags$div(
            field_label("Top-N groups", "Accumulative modeling builds classifiers from the top-1, top-2, ... up to top-N ranked groups, reporting performance at each level. The best k (by F1) is chosen per fold."),
            numericInput("bp_n_top_groups", NULL, value = 10, min = 1))),
        field_label("Final modeling", "Classifier used to build the final accumulative model from the top-ranked groups. 'Same as scoring' reuses the classifier above; otherwise pick a different algorithm for the modeling step."),
        selectInput("bp_model_method", NULL,
                    choices = c("Same as scoring" = "same", CLASSIFIER_CHOICES)),
        layout_columns(col_widths = c(6, 6),
          tags$div(style = "display:flex; align-items:center; gap:4px;",
            checkboxInput("bp_shave", "Gene shaving", value = TRUE),
            help_tip("Before scoring, drop genes that are not individually differential (Welch t-test + BH FDR on the TRAINING split). Focuses each group on its informative genes.")),
          tags$div(style = "display:flex; align-items:center; gap:4px;",
            checkboxInput("bp_normalize", "Z-normalize", value = FALSE),
            help_tip("Standardize each gene to mean 0 / SD 1 (fit on training, applied to test). Helps distance-based classifiers (KNN, SVM); usually unnecessary for tree models."))),
        tags$div(style = "display:flex; align-items:center; gap:4px;",
          checkboxInput("bp_discover", "Module discovery (per dataset)", value = FALSE),
          help_tip("After MCCV, run de-novo module discovery (RRA-ranked robust genes -> co-selection network -> Louvain clustering with bootstrap validation) for each dataset. Adds time but yields data-driven gene modules.")),
        field_label("Speed / memory", "CPU cores = total workers shared across ALL datasets' folds (every (dataset, fold) unit runs over one pool). Max resident = how many datasets are held in memory at once (bounds RAM; folds still saturate all cores)."),
        layout_columns(col_widths = c(7, 5),
          tags$div(
            field_label("CPU cores", "Total parallel workers shared across all datasets and folds. Higher = faster, up to the machine's core count."),
            sliderInput("bp_cores", NULL, min = 1, max = parallel::detectCores(),
                      value = max(1, parallel::detectCores() - 1), step = 1)),
          tags$div(
            field_label("Max resident", "How many datasets are loaded in memory at the same time. Lower bounds peak RAM (datasets load/release just-in-time); cores stay busy as long as a couple datasets supply enough folds."),
            numericInput("bp_max_resident", NULL, value = 2, min = 1, max = 20))),
        layout_columns(col_widths = c(6, 6),
          tags$div(
            field_label("Seed", "Random seed for reproducibility. Fold f uses seed+f, so the same seed reproduces identical splits and results."),
            numericInput("bp_seed", NULL, value = 44, min = 1)),
          tags$div(
            field_label("Report", "Format of the combined cross-dataset report: a multi-sheet Excel workbook (.xlsx) or a zip of CSV files."),
            selectInput("bp_report_format", NULL,
                      choices = c("Excel (.xlsx)" = "xlsx", "CSV (.zip)" = "csv")))),
        tags$hr(),
        actionButton("btn_batch", "\U0001F680 Run Batch",
                     class = "btn-primary btn-lg w-100")
      ),
      # Main panel
      uiOutput("batch_overall"),
      uiOutput("batch_alert"),
      card(
        card_header("\U0001F4CA Per-dataset progress"),
        uiOutput("batch_status")
      ),
      conditionalPanel(
        condition = "output.batch_done",
        card(
          card_header("\U0001F9EC Combined cross-dataset report"),
          layout_columns(col_widths = c(4, 4, 4),
            downloadButton("dl_batch_report", "\U0001F4E5 Combined report",
                           class = "btn-primary w-100"),
            downloadButton("dl_batch_all", "\U0001F4E6 All results (zip)",
                           class = "btn-outline-light w-100"),
            tags$div(class = "text-muted small", textOutput("batch_outdir_note"))),
          tags$hr(),
          tags$strong("Run details (parameters + timing)"),
          DTOutput("batch_run_info"),
          tags$hr(),
          tags$strong("Per-dataset performance"),
          DTOutput("batch_combined_metrics")
        )
      ),
      card(
        card_header("\U0001F4DC Batch log"),
        tags$pre(style = "max-height:240px; overflow:auto; font-size:0.8rem;",
                 textOutput("batch_log"))
      )
    )
  ),

  # ── GitHub link (stays in the right utility cluster) ────────────────────
  nav_item(
    tags$a(href = "https://github.com/jallmer/GSM4R",
           target = "_blank",
           style = "color: #888; font-size: 0.8rem; align-self: center;",
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
    filter_report = NULL,
    groups     = NULL,   # named list
    dataset    = NULL,   # gsm_dataset object
    loaded     = FALSE,

    # Collapse
    collapse_result = NULL,
    collapsed       = FALSE,

    # Grouping (G) tab — standalone run + inspect results per method
    g_file   = NULL,   # list(groups)
    g_bridge = NULL,   # list(groups, diagnostics, summary, graph, spec)

    # Pre-score
    prescore_groups = NULL,
    prescored       = FALSE,
    train_idx       = NULL,

    # Score
    score_results = NULL,
    scored        = FALSE,
    detect        = NULL,

    # Group significance (permutation test)
    permtest      = NULL,

    # Operation timings (collected across the session, written to Timings sheet)
    # Named list with fields:
    #   data_load      : seconds for Data -> Load Selected Example
    #   grouping_load  : seconds for Grouping -> Prior knowledge -> Load grouping
    #   bridge_runs    : named list (estimator -> seconds) for BRIDGE runs
    #   gsm_run        : seconds for GSM -> Run Analysis
    timings       = list()
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
    cores = 1, folds = NULL, params = NULL, t0 = NULL,
    # Non-blocking parallel scheduler state:
    running = list(),       # named (fold idx) list of mcparallel job handles in flight
    queue = integer(0),     # fold indices not yet launched
    started = 0,            # how many folds have been launched
    fold_started = list(),  # per-fold launch time (secs since t0)
    summary_open = FALSE,   # keep the completion banner on screen until dismissed
    last_summary = NULL     # snapshot shown in the persistent "done" banner
  )

  # G-tab BRIDGE run: cancellable background job (mirrors the MCCV scheduler).
  gbr <- reactiveValues(active = FALSE, job = NULL, t0 = NULL, tick = 0)

  # ── Output flags for conditionalPanel ──────────────────────────────────
  output$data_loaded     <- reactive(rv$loaded)
  output$g_bridge_running <- reactive(isTRUE(gbr$active))
  output$collapse_done   <- reactive(rv$collapsed)
  output$prescore_done   <- reactive(rv$prescored)
  output$prescore_running <- reactive(isTRUE(ps$active))
  output$score_done      <- reactive(rv$scored)
  output$mccv_running    <- reactive(isTRUE(mc$active))
  output$permtest_done   <- reactive(!is.null(rv$permtest))
  outputOptions(output, "data_loaded",      suspendWhenHidden = FALSE)
  outputOptions(output, "collapse_done",    suspendWhenHidden = FALSE)
  outputOptions(output, "prescore_done",    suspendWhenHidden = FALSE)
  outputOptions(output, "prescore_running", suspendWhenHidden = FALSE)
  outputOptions(output, "score_done",       suspendWhenHidden = FALSE)
  outputOptions(output, "mccv_running",     suspendWhenHidden = FALSE)
  outputOptions(output, "g_bridge_running", suspendWhenHidden = FALSE)
  outputOptions(output, "permtest_done",    suspendWhenHidden = FALSE)

  # ══════════════════════════════════════════════════════════════════════
  # Tab 1: Load Data
  # ══════════════════════════════════════════════════════════════════════
  # Load EXPRESSION ONLY (+ label-blind filtering). Grouping is chosen later in
  # the Grouping / GSM tabs; rv$dataset is built when a grouping is committed.
  load_into_rv <- function(expr_path) {
    output$expr_file_status <- renderUI(NULL)
    tryCatch({
      withProgress(message = "Loading data...", value = 0.1, {

        # ── Quick-validate expression file ──────────────────────────────
        expr_peek <- tryCatch(
          data.table::fread(expr_path, nrows = 3, header = TRUE, check.names = FALSE),
          error = function(e) NULL)
        if (is.null(expr_peek) || ncol(expr_peek) < 3) {
          output$expr_file_status <- renderUI(tags$div(class = "text-danger mt-1",
            style = "font-weight:600;",
            "❌ Wrong file — expected an expression CSV with a 'class' column and gene columns."))
          stop("Expression file does not look like a valid expression CSV.")
        }
        if (!"class" %in% names(expr_peek)) {
          if (ncol(expr_peek) == 2) {
            output$expr_file_status <- renderUI(tags$div(class = "text-danger mt-1",
              style = "font-weight:600;",
              "❌ This looks like a grouping file, not expression data. Grouping is loaded in the Grouping tab."))
            stop("The uploaded expression file looks like a grouping file (2 columns).")
          }
          output$expr_file_status <- renderUI(tags$div(class = "text-danger mt-1",
            style = "font-weight:600;", "❌ Missing 'class' column. Found: ",
            paste(head(names(expr_peek), 5), collapse = ", "),
            if (ncol(expr_peek) > 5) "..." else ""))
          stop("Class column 'class' not found in expression file.")
        }
        setProgress(0.3, detail = "Reading expression...")

        # ── Parse expression ────────────────────────────────────────────
        expr_data <- read_expression(expr_path)
        n_before  <- nrow(expr_data$expr)

        # ── Label-blind filtering (denoising) ───────────────────────────
        setProgress(0.6, detail = "Filtering genes...")
        rep <- NULL
        if (isTRUE(input$flt_missing) || isTRUE(input$flt_lowexpr) ||
            isTRUE(input$flt_variance)) {
          mm <- if (isTRUE(input$flt_lowexpr)) input$flt_min_mean else NULL
          mg <- if (isTRUE(input$flt_variance)) {
            v <- suppressWarnings(as.integer(input$flt_max_genes))
            if (length(v) && !is.na(v) && v > 0) v else NULL
          } else NULL
          expr_data$expr <- filter_genes(expr_data$expr,
                                         drop_constant = isTRUE(input$flt_missing),
                                         min_mean = mm, max_genes = mg)
          rep <- attr(expr_data$expr, "filter_report")
        }
        n_after <- nrow(expr_data$expr)

        rv$expr_data <- expr_data
        rv$filter_report <- rep
        rv$loaded    <- TRUE
        rv$dataset   <- NULL      # no grouping chosen yet
        rv$groups    <- NULL
        rv$g_file    <- NULL
        rv$g_bridge  <- NULL
        rv$dataset_name <- tools::file_path_sans_ext(basename(expr_path))
        rv$expr_file_path <- normalizePath(expr_path, mustWork = FALSE)
        rv$collapsed <- FALSE; rv$prescored <- FALSE; rv$scored <- FALSE
        rv$prescore_groups <- NULL; rv$score_results <- NULL; rv$mccv_results <- NULL

        output$expr_file_status <- renderUI(tags$div(class = "text-success mt-1",
          sprintf("✅ %d samples × %d genes loaded%s", ncol(expr_data$expr),
                  n_after, if (n_after < n_before)
                    sprintf(" (%d after filtering)", n_after) else "")))
        setProgress(1, detail = "Done!")
      })
      showNotification("✅ Expression loaded. Now choose a grouping in the Grouping tab.",
                       type = "message")
    }, error = function(e) {
      showNotification(paste("❌ Error:", e$message), type = "error", duration = 10)
    })
  }

  # ── Instant validation on file upload ──────────────────────────────────
  observeEvent(input$expr_file, {
    req(input$expr_file)
    peek <- tryCatch(
      data.table::fread(input$expr_file$datapath, nrows = 3, header = TRUE, check.names = FALSE),
      error = function(e) NULL
    )
    if (is.null(peek) || ncol(peek) < 3) {
      output$expr_file_status <- renderUI(
        tags$div(class = "text-danger mt-1", style = "font-weight:600;",
                 "\u274c Wrong file type — expected an expression CSV with many gene columns.")
      )
    } else if (!"class" %in% names(peek)) {
      if (ncol(peek) == 2) {
        output$expr_file_status <- renderUI(
          tags$div(class = "text-danger mt-1", style = "font-weight:600;",
                   "\u274c Wrong file — this looks like a grouping file (2 columns). ",
                   "Upload it into the Grouping File box instead.")
        )
      } else {
        output$expr_file_status <- renderUI(
          tags$div(class = "text-danger mt-1", style = "font-weight:600;",
                   "\u274c Missing 'class' column. Found: ",
                   paste(head(names(peek), 5), collapse = ", "),
                   if (ncol(peek) > 5) "..." else "")
        )
      }
    } else {
      output$expr_file_status <- renderUI(
        tags$div(class = "text-success mt-1", "\u2705 Expression file looks valid")
      )
    }
  })

  observeEvent(input$groups_file, {
    req(input$groups_file)
    ext <- tolower(tools::file_ext(input$groups_file$name))
    if (ext == "gmt") {
      # GMT files are fine; can't easily peek
      output$groups_file_status <- renderUI(
        tags$div(class = "text-success mt-1", "\u2705 GMT file selected")
      )
      return()
    }
    peek <- tryCatch(
      data.table::fread(input$groups_file$datapath, nrows = 3, header = TRUE, check.names = FALSE),
      error = function(e) NULL
    )
    if (is.null(peek) || ncol(peek) < 2) {
      output$groups_file_status <- renderUI(
        tags$div(class = "text-danger mt-1", style = "font-weight:600;",
                 "\u274c Wrong file type — expected a 2-column CSV (gene, group) or GMT.")
      )
    } else if ("class" %in% names(peek) && ncol(peek) > 10) {
      output$groups_file_status <- renderUI(
        tags$div(class = "text-danger mt-1", style = "font-weight:600;",
                 "\u274c Wrong file — this looks like expression data (has 'class' column, ",
                 ncol(peek), " columns). Upload it into the Expression Data box instead.")
      )
    } else if (ncol(peek) == 2) {
      output$groups_file_status <- renderUI(
        tags$div(class = "text-success mt-1", "\u2705 Grouping file looks valid")
      )
    } else {
      output$groups_file_status <- renderUI(
        tags$div(class = "text-muted mt-1", "\u2139\ufe0f File has ", ncol(peek),
                 " columns — first two will be used as (gene, group).")
      )
    }
  })

  observeEvent(input$btn_load, {
    req(input$expr_file)
    t0_data <- Sys.time()
    load_into_rv(input$expr_file$datapath)
    elapsed_data <- as.numeric(difftime(Sys.time(), t0_data, units = "secs"))
    timings_cur <- rv$timings %||% list()
    timings_cur$data_load <- elapsed_data
    timings_cur$data_load_detail <- "Load & Validate (Uploaded expression)"
    rv$timings <- timings_cur
  })

  # ── Filter summary (before/after gene counts) ────────────────────────────
  output$filter_summary <- renderUI({
    rep <- rv$filter_report
    if (is.null(rep)) return(NULL)
    tags$div(class = "text-muted", style = "font-size:0.78rem; margin-top:2px;",
      sprintf("Genes: %d → %d kept (after missing/constant %d, low-expr %d, variance %d).",
              rep[["input"]], rep[["after_variance"]],
              rep[["after_missing"]], rep[["after_lowexpr"]], rep[["after_variance"]]))
  })

  # Load a selected example expression dataset (bundled toy or a repo GDS file).
  observeEvent(input$btn_example, {
    expr_sel <- input$example_expr
    expr_path <- if (expr_sel == "__toy__") {
      system.file("extdata", "test_expression.csv", package = "gsmR")
    } else file.path(.example_input_dir, expr_sel)
    if (!nzchar(expr_path) || !file.exists(expr_path)) {
      showNotification(paste("\u274C Expression file not found:", expr_sel),
                       type = "error", duration = 10); return()
    }
    t0_data <- Sys.time()
    load_into_rv(expr_path)
    elapsed_data <- as.numeric(difftime(Sys.time(), t0_data, units = "secs"))
    timings_cur <- rv$timings %||% list()
    timings_cur$data_load <- elapsed_data
    timings_cur$data_load_detail <- "Load Selected Example"
    rv$timings <- timings_cur
  })

  # Value boxes
  output$vbox_samples <- renderUI({
    req(rv$expr_data)
    value_box_ui(ncol(rv$expr_data$expr), "Samples", "\U0001F9EA")
  })
  output$vbox_genes <- renderUI({
    req(rv$expr_data)
    value_box_ui(format(nrow(rv$expr_data$expr), big.mark = ","), "Genes", "\U0001F9EC")
  })
  output$vbox_groups <- renderUI({
    req(rv$expr_data)
    if (is.null(rv$dataset))
      return(value_box_ui("\u2014", "Groups (choose in GSM)", "\U0001F5C2\uFE0F"))
    deferred <- identical(rv$dataset$params$grouping_method, "bridge") &&
                rv$dataset$params$n_groups == 0
    val <- if (deferred) "per-fold"
           else format(rv$dataset$params$n_groups, big.mark = ",")
    value_box_ui(val, if (deferred) "Groups (BRIDGE)" else "Groups", "\U0001F5C2\uFE0F")
  })
  output$vbox_classes <- renderUI({
    req(rv$expr_data)
    value_box_ui(nlevels(rv$expr_data$labels), "Classes", "\U0001F3F7\uFE0F")
  })

  # Plots
  output$plot_class_dist <- renderPlot({
    req(rv$expr_data)
    tbl <- table(rv$expr_data$labels)
    par(bg = "transparent", col.axis = "#aaaaaa", col.lab = "#aaaaaa",
        col.main = "#cccccc", fg = "#555555", mar = c(4, 4, 2, 1))
    bp <- barplot(tbl, col = c("#4a90d9", "#e06666", "#6aa84f", "#f1c232")[
      seq_along(tbl)
    ], border = NA, main = "", ylab = "Count", las = 1)
    text(bp, tbl / 2, labels = tbl, col = "white", font = 2, cex = 1.2)
  }, bg = "transparent")

  output$plot_group_sizes <- renderPlot({
    req(rv$expr_data)
    par(bg = "transparent", col.axis = "#aaaaaa", col.lab = "#aaaaaa",
        col.main = "#cccccc", fg = "#555555", mar = c(4, 4, 2, 1))
    sizes <- if (is.null(rv$dataset)) integer(0)
             else vapply(rv$dataset$groups, length, integer(1))
    if (length(sizes) == 0) {
      plot.new()
      text(0.5, 0.5, paste0(
        "No grouping selected yet.\n",
        "Prepare a grouping in the Grouping tab,\n",
        "then choose it in the GSM tab."),
        col = "#aaaaaa", cex = 1.05)
      return(invisible())
    }
    hist(sizes, breaks = 30, col = "#4a90d9", border = NA,
         main = "", xlab = "Genes per group", ylab = "Frequency", las = 1)
    abline(v = median(sizes), col = "#f1c232", lwd = 2, lty = 2)
    legend("topright", legend = paste("Median:", median(sizes)),
           col = "#f1c232", lwd = 2, lty = 2, bty = "n", text.col = "#aaaaaa")
  }, bg = "transparent")

  # Expression preview
  output$expr_preview <- renderDT({
    req(rv$expr_data)
    # Show first 100 genes x all samples
    mat <- rv$expr_data$expr[1:min(100, nrow(rv$expr_data$expr)), , drop = FALSE]
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
    if (is.null(rv$dataset)) {
      showNotification("⚠️ Choose a grouping first (Grouping tab → Apply, or the GSM grouping panel).",
                       type = "warning", duration = 10); return()
    }
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
                 sprintf("elapsed %s · ~%s remaining",
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
  alert_box <- function(type, title, body = NULL, closable = FALSE) {
    pal <- switch(type,
      ok    = list(bg = "rgba(106,168,79,0.15)",  bd = "#6aa84f", ic = "✅"),
      warn  = list(bg = "rgba(241,194,50,0.15)",  bd = "#f1c232", ic = "⚠️"),
      error = list(bg = "rgba(224,102,102,0.18)", bd = "#e06666", ic = "⛔"),
      info  = list(bg = "rgba(74,144,217,0.15)",  bd = "#4a90d9", ic = "ℹ️"))
    tags$div(
      class = "value-box",
      style = sprintf("background:%s; border-left:5px solid %s; position:relative;", pal$bg, pal$bd),
      if (closable) tags$button(type = "button", class = "btn-close",
                style = "position:absolute; top:8px; right:8px; filter:invert(1); opacity:0.6;",
                onclick = "this.parentElement.style.display='none';",
                `aria-label` = "Close"),
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
    # Show a dismissable "done" banner once results are ready
    if (isTRUE(rv$scored)) {
      return(alert_box("ok", "Analysis complete",
             paste0("Results are ready in the tabs below. ",
                    "Click \U0001F30B Discover to run de novo module discovery."),
             closable = TRUE))
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
    session$sendCustomMessage("disableBtn", "btn_score")
    session$sendCustomMessage("showBanner", list(
      id = "score_status_alert",
      title = "Analysis running…",
      detail = "Computing scores for each group. This may take a while."
    ))
    on.exit(session$sendCustomMessage("enableBtn", "btn_score"))
    # Only a loaded dataset is required. Pre-Score is optional — Score creates
    # its own train/test split and pre-scores on the fly when needed.
    if (is.null(rv$expr_data)) {
      showNotification("⚠️ Load expression data on the Data tab first.",
                       type = "warning", duration = 8)
      return()
    }
    if (is.null(rv$dataset)) {
      showNotification("⚠️ Choose a grouping first: open the 'Grouping' card in the settings, pick a method, and click 'Use this grouping'.",
                       type = "warning", duration = 10)
      return()
    }

    # ── Snapshot of what THIS run used (shown live in the config panel) ───────
    rv$run_config <- build_run_config(input, rv)

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
      if (m == "xgb" && !requireNamespace("xgboost", quietly = TRUE)) {
        showNotification("\u274C Package 'xgboost' required. Install with: install.packages('xgboost')",
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
    } else if (input$score_method == "xgb") {
      score_args <- list(nrounds = 100)
    } else if (input$score_method %in% c("mlp", "gnn")) {
      score_args <- build_nn_args(input, input$score_method)
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

        # Leakage-free grouping: for BRIDGE, RE-DERIVE groups on this run's
        # TRAINING samples only (rather than the load-time preview groups, which
        # saw all samples). File-based grouping is returned unchanged.
        ds_run <- rv$dataset
        if (!is.null(rv$dataset$grouping_spec) &&
            identical(rv$dataset$grouping_spec$method, "bridge")) {
          t0_bridge_single <- Sys.time()
          withProgress(message = "BRIDGE grouping (training samples only)...",
                       value = 0.2, {
            ds_run$groups <- grouping_for_fold(rv$dataset, rv$train_idx)
          })
          elapsed_bridge_single <- as.numeric(difftime(Sys.time(), t0_bridge_single, units = "secs"))
          timings_cur <- rv$timings %||% list()
          br <- timings_cur$bridge_runs %||% list()
          est_lbl <- rv$dataset$grouping_spec$args$precision %||% "bridge"
          if (is.null(br[[est_lbl]])) {
            br[[est_lbl]] <- elapsed_bridge_single
            timings_cur$bridge_runs <- br
            rv$timings <- timings_cur
          }
          validate(need(length(ds_run$groups) > 0,
                        "BRIDGE derived no groups on the training split. Lower tau or min group size."))
          ds_run$params$n_groups <- length(ds_run$groups)
          rv$prescore_groups <- NULL   # preview pre-scoring no longer applies
        }

        # Candidate groups: all groups, or the pre-scored top set (auto if needed).
        scope_all <- identical(input$score_scope, "all")
        auto_prescored <- FALSE
        if (scope_all) {
          candidates <- ds_run$groups
        } else {
          candidates <- rv$prescore_groups
          if (is.null(candidates)) {
            candidates <- pre_score_groups(ds_run, rv$train_idx,
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

          # Phase I — training-only gene shaving (Welch t-test; FDR or top-k).
          if (isTRUE(input$shave_genes)) {
            keep_genes <- shave_genes(expr_use, rv$dataset$labels,
                                      rv$train_idx, alpha = input$fdr_alpha,
                                      method = input$shave_method %||% "fdr",
                                      top_k = input$shave_k %||% 100)
            candidates <- shave_groups(candidates, keep_genes)
            n_shaved <- length(keep_genes)
            validate(need(length(candidates) > 0,
                          "No groups survived gene shaving — relax the FDR threshold or raise Top-K."))
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
            # Store gene list so export tables can report N_genes in BRIDGE mode
            # (dataset$groups is empty when grouping is BRIDGE per-fold)
            results[[gname]] <- list(accuracy = scores[i],
                                     genes = candidates[[gname]])
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
        session$sendCustomMessage("enableBtn", "btn_score")
          setProgress(1, detail = "Done!")
        })

        n <- length(candidates)
        rv$score_summary <- list(
          mode = "single",
          dataset_name = rv$dataset_name %||% "dataset",
          modeling = sprintf("accumulative top-1..%d (best k = %d)", K, best_k),
          best_k = best_k,
          gene_shaving = if (isTRUE(input$shave_genes))
            if (identical(input$shave_method, "topk"))
              sprintf("Welch t-test, top-%d genes (%d kept)", input$shave_k %||% 100L, n_shaved)
            else
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
        # Record GSM run timing.
        timings_cur <- rv$timings %||% list()
        timings_cur$gsm_run <- rv$time_score
        rv$timings <- timings_cur
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
        } else if (m_method == "xgb") {
          model_args <- list(nrounds = 100)
        } else if (m_method %in% c("mlp", "gnn")) {
          model_args <- build_nn_args(input, m_method)
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
          fdr_alpha = input$fdr_alpha,
          shave_method = input$shave_method %||% "fdr",
          shave_k = input$shave_k %||% 100,
          normalize = isTRUE(input$normalize_z),
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
        # Scheduler state: queue every fold, nothing running yet.
        mc$running      <- list()
        mc$queue        <- seq_len(mc$n)
        mc$started      <- 0
        mc$fold_started <- list()
        mc$summary_open <- FALSE
        mc$last_summary <- NULL
        # Per-fold live-progress files written by the forked workers.
        mc$progress_dir <- tempfile("gsm_mccv_")
        dir.create(mc$progress_dir, showWarnings = FALSE)
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
      # 1) Reap folds whose forked worker has finished — non-blocking, so the
      #    Shiny loop stays responsive and the live panel keeps repainting
      #    while the remaining folds keep running on their cores.
      if (length(mc$running) > 0) {
        running  <- mc$running
        folds    <- mc$folds
        done     <- mc$done
        finished <- character(0)
        for (k in names(running)) {
          got <- parallel::mccollect(running[[k]], wait = FALSE, timeout = 0)
          if (!is.null(got)) {
            idx <- as.integer(k)
            val <- got[[1]]
            folds[[idx]] <- if (is.null(val) || inherits(val, "try-error"))
              list(fold = idx, status = "error", error_message = "Fold worker crashed")
            else val
            finished <- c(finished, k)
            done <- done + 1L
          }
        }
        if (length(finished) > 0) {
          for (k in finished) running[[k]] <- NULL
          mc$running <- running; mc$folds <- folds; mc$done <- done
          mc$tick <- mc$tick + 1
        }
      }

      # 2) All folds launched AND reaped -> aggregate and finish.
      if (mc$done >= mc$n && length(mc$running) == 0) {
        # \u2500\u2500 Finalize: aggregate exactly like run_mccv() \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500
        results <- mc$folds
        ok <- Filter(function(r) !is.null(r) && r$status == "success", results)
        n_failed <- length(results) - length(ok)
        if (length(ok) == 0) {
          mc$active <- FALSE
          # Collect error messages from failed folds
          err_msgs <- unique(vapply(results, function(r) {
            if (is.null(r)) "Fold returned NULL" else r$error_message %||% "Unknown error"
          }, character(1)))
          showNotification(
            paste0("\u274c All ", length(results), " MCCV folds failed. Errors: ",
                   paste(err_msgs, collapse = "; ")),
            type = "error", duration = 15)
          return()
        }
        if (n_failed > 0) {
          showNotification(
            sprintf("\u26a0\ufe0f %d of %d folds failed. Results are based on %d successful folds.",
                    n_failed, length(results), length(ok)),
            type = "warning", duration = 10)
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
          gene_shaving = if (isTRUE(p$shave)) {
            mean_shaved <- mean(vapply(ok, function(r) r$n_shaved_genes %||% NA_real_,
                                       numeric(1)), na.rm = TRUE)
            if (identical(p$shave_method, "topk"))
              sprintf("Welch t-test, top-%d genes/fold (mean %.0f kept)",
                      p$shave_k %||% 100L, mean_shaved)
            else
              sprintf("Welch t-test + BH FDR %.3g (mean %.0f genes/fold)",
                      p$fdr_alpha %||% 0.05, mean_shaved)
          } else "off",
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
        # Record GSM run timing.
        timings_cur <- rv$timings %||% list()
        timings_cur$gsm_run <- rv$time_score
        rv$timings <- timings_cur
        # Snapshot for the persistent completion banner (stays until dismissed).
        mc_metrics <- colMeans(metrics_df, na.rm = TRUE)
        mc$last_summary <- list(
          n = mc$n, ok = length(ok), failed = n_failed, cores = mc$cores,
          elapsed = rv$time_score, mean_acc = mc_metrics[["Accuracy"]] %||% NA_real_,
          metrics = mc_metrics)
        mc$summary_open <- TRUE
        mc$active <- FALSE
        showNotification(sprintf("\u2705 Parallel MCCV finished %d/%d folds in %s!",
                                 length(ok), mc$n, format_duration(rv$time_score)),
                         type = "message")
        return()
      }

      # 3) Keep the core budget full: launch queued folds as background forks
      #    until `cores` are in flight. Each fold is independently seeded
      #    (seed + fold in mccv_run_fold), so parallel == sequential results.
      p     <- mc$params
      cores <- max(1L, as.integer(mc$cores %||% 1L))
      if (length(mc$queue) > 0 && .Platform$OS.type != "windows") {
        running <- mc$running
        queue   <- mc$queue
        started <- mc$fold_started
        pdir    <- mc$progress_dir     # read reactive HERE (parent), not in fork
        while (length(running) < cores && length(queue) > 0) {
          idx <- queue[1]; queue <- queue[-1]
          # Plain-string path built in the parent; the fork must not touch any
          # reactive object (would error: reactive context accessed cross-process).
          pf <- file.path(pdir, paste0("fold_", idx, ".txt"))
          job <- parallel::mcparallel(
            tryCatch(
              mccv_run_fold(
                p$dataset, idx, train_ratio = p$train_ratio,
                pre_score_method = p$pre_score_method,
                pre_score_top_percent = p$pre_score_top_percent,
                score_method = p$score_method, model_method = p$model_method,
                score_args = p$score_args, model_args = p$model_args,
                top_n = p$top_n %||% 10, shave = p$shave %||% TRUE,
                fdr_alpha = p$fdr_alpha %||% 0.05,
                shave_method = p$shave_method %||% "fdr", shave_k = p$shave_k %||% 100,
                normalize = p$normalize %||% FALSE,
                seed = p$seed,
                progress_file = pf),
              error = function(e) list(fold = idx, status = "error",
                                       error_message = conditionMessage(e))),
            name = as.character(idx))
          running[[as.character(idx)]] <- job
          started[[as.character(idx)]] <-
            as.numeric(difftime(Sys.time(), mc$t0, units = "secs"))
        }
        mc$running <- running; mc$queue <- queue; mc$fold_started <- started
        mc$started <- mc$n - length(queue)
      } else if (length(mc$queue) > 0) {
        # No fork available (Windows): run one fold synchronously per tick.
        idx <- mc$queue[1]; mc$queue <- mc$queue[-1]
        mc$started <- mc$n - length(mc$queue)
        fr <- tryCatch(
          mccv_run_fold(
            p$dataset, idx, train_ratio = p$train_ratio,
            pre_score_method = p$pre_score_method,
            pre_score_top_percent = p$pre_score_top_percent,
            score_method = p$score_method, model_method = p$model_method,
            score_args = p$score_args, model_args = p$model_args,
            top_n = p$top_n %||% 10, shave = p$shave %||% TRUE,
            fdr_alpha = p$fdr_alpha %||% 0.05,
            shave_method = p$shave_method %||% "fdr", shave_k = p$shave_k %||% 100,
            normalize = p$normalize %||% FALSE,
            seed = p$seed,
            progress_file = file.path(mc$progress_dir,
                                      paste0("fold_", idx, ".txt"))),
          error = function(e) list(fold = idx, status = "error",
                                   error_message = conditionMessage(e)))
        folds <- mc$folds; folds[[idx]] <- fr; mc$folds <- folds
        mc$done <- mc$done + 1L; mc$tick <- mc$tick + 1
      }
      NULL
    }) }, error = function(e) conditionMessage(e))

    # On any unexpected error, stop cleanly so the button isn't dead-locked.
    if (!is.null(err)) {
      mc$active <- FALSE
      showNotification(paste("⛔ MCCV stopped:", err), type = "error", duration = 12)
      return()
    }
    # Poll briskly while work is in flight so the live panel stays smooth.
    if (isTRUE(mc$active)) invalidateLater(200)
  })

  # ── Stop a running MCCV: kill in-flight forks, clear the queue, finish. ────
  observeEvent(input$btn_stop_score, {
    if (!isTRUE(mc$active)) return()
    # Terminate the forked fold workers still running on their cores.
    for (job in mc$running) {
      pid <- tryCatch(job$pid, error = function(e) NULL)
      if (!is.null(pid))
        tryCatch(tools::pskill(pid, tools::SIGKILL), error = function(e) NULL)
    }
    # Reap any now-dead workers so no zombies linger.
    if (length(mc$running) > 0)
      tryCatch(parallel::mccollect(mc$running, wait = FALSE), error = function(e) NULL)

    n_done <- mc$done
    mc$running     <- list()
    mc$queue       <- integer(0)
    mc$active      <- FALSE
    mc$summary_open <- FALSE     # no completion banner for a cancelled run
    rv$scored      <- FALSE
    session$sendCustomMessage("enableBtn", "btn_score")
    showNotification(
      sprintf("⏹️ MCCV stopped by user after %d/%d fold(s). No results saved.",
              n_done, mc$n %||% 0L),
      type = "warning", duration = 8)
  })

  # Helper: accuracy vector from completed folds (in fold order).
  mccv_done_folds <- function() {
    Filter(function(r) !is.null(r) && r$status == "success", mc$folds)
  }

  output$mccv_live_status <- renderUI({
    # -- Persistent completion banner: stays on screen until dismissed --
    if (!isTRUE(mc$active)) {
      if (!isTRUE(mc$summary_open) || is.null(mc$last_summary)) return(NULL)
      ls <- mc$last_summary; m <- ls$metrics
      metric_chip <- function(lab, val) tags$div(
        style = "flex:1 1 90px; text-align:center; padding:8px 6px; background:rgba(255,255,255,0.04); border-radius:8px;",
        tags$div(class = "text-muted", style = "font-size:.7rem; text-transform:uppercase; letter-spacing:.04em;", lab),
        tags$div(style = "font-size:1.15rem; font-weight:700;",
                 if (is.null(val) || is.na(val)) "—" else sprintf("%.3f", val)))
      return(tags$div(
        class = "value-box",
        style = "background: rgba(46,160,90,0.12); border-left: 4px solid #2ea05a;",
        tags$div(style = "display:flex; justify-content:space-between; align-items:center; gap:1rem;",
          tags$div(style = "font-size:1.35rem; font-weight:700;",
                   HTML(sprintf("✅ MCCV complete &mdash; %d/%d folds in %s",
                                ls$ok, ls$n, format_duration(ls$elapsed)))),
          actionButton("mccv_dismiss_summary", "Dismiss",
                       class = "btn btn-sm btn-outline-light")),
        tags$div(class = "text-muted", style = "margin:4px 0 10px;",
                 sprintf("Ran on %d CPU core(s)%s · mean accuracy %s",
                         ls$cores,
                         if (ls$failed > 0) sprintf(" · %d fold(s) failed", ls$failed) else "",
                         if (is.na(ls$mean_acc)) "—" else sprintf("%.1f%%", 100 * ls$mean_acc))),
        tags$div(style = "display:flex; gap:8px; flex-wrap:wrap;",
                 lapply(names(m), function(nm) metric_chip(gsub("_", " ", nm), m[[nm]])))
      ))
    }

    # -- Live running panel --
    invalidateLater(500)                 # keep the elapsed clock ticking
    n       <- mc$n; done <- mc$done
    running <- mc$running; n_run <- length(running); n_queue <- length(mc$queue)
    cores   <- max(1L, as.integer(mc$cores %||% 1L))
    pct     <- if (n > 0) round(100 * done / n) else 0
    elapsed <- as.numeric(difftime(Sys.time(), mc$t0, units = "secs"))
    eta     <- if (done > 0) elapsed / done * (n - done) else NA
    ok      <- mccv_done_folds()
    accs    <- vapply(ok, function(r) as.numeric(r$metrics$Accuracy), numeric(1))
    meanacc <- if (length(accs)) sprintf("%.1f%%", 100 * mean(accs)) else "—"

    run_idx  <- as.integer(names(running))
    fail_idx <- which(vapply(mc$folds, function(r)
                        !is.null(r) && identical(r$status, "error"), logical(1)))
    acc_by_fold <- stats::setNames(rep(NA_real_, n), as.character(seq_len(n)))
    for (r in ok) acc_by_fold[as.character(r$fold)] <- as.numeric(r$metrics$Accuracy)

    stat <- function(lab, val, col = "#eaeaea") tags$div(
      style = "flex:1 1 90px; text-align:center;",
      tags$div(style = sprintf("font-size:1.3rem; font-weight:700; color:%s;", col), val),
      tags$div(class = "text-muted",
               style = "font-size:.7rem; text-transform:uppercase; letter-spacing:.04em;", lab))

    chip <- function(i) {
      ik <- as.character(i)
      running_now <- (i %in% run_idx) && !(i %in% fail_idx) && is.na(acc_by_fold[ik])
      if (i %in% fail_idx) {
        bg <- "rgba(217,83,79,0.20)"; bd <- "#d9534f"; txt <- sprintf("F%d ✗", i)
      } else if (!is.na(acc_by_fold[ik])) {
        bg <- "rgba(46,160,90,0.20)"; bd <- "#2ea05a"
        txt <- sprintf("F%d %.0f%%", i, 100 * acc_by_fold[ik])
      } else if (running_now) {
        bg <- "rgba(74,144,217,0.30)"; bd <- "#4a90d9"; txt <- sprintf("F%d ▶", i)
      } else {
        bg <- "rgba(255,255,255,0.04)"; bd <- "rgba(255,255,255,0.12)"; txt <- sprintf("F%d", i)
      }
      tags$span(class = if (running_now) "gsm-chip gsm-chip-run" else "gsm-chip",
                style = sprintf("background:%s; border:1px solid %s;", bg, bd), txt)
    }
    grid <- if (n <= 120)
      tags$div(style = "display:flex; flex-wrap:wrap; gap:6px; margin-top:10px;",
               lapply(seq_len(n), chip))
    else
      tags$div(class = "text-muted", style = "margin-top:10px;",
               sprintf("%d folds — per-fold grid hidden for readability.", n))

    css <- paste0(
      ".gsm-chip{display:inline-block;min-width:54px;text-align:center;padding:5px 8px;",
      "border-radius:7px;font-size:.78rem;font-weight:600;font-variant-numeric:tabular-nums;}",
      "@keyframes gsmPulse{0%{opacity:.5}50%{opacity:1}100%{opacity:.5}}",
      ".gsm-chip-run{animation:gsmPulse 1s ease-in-out infinite;}")

    tags$div(
      class = "value-box",
      style = "background: rgba(74,144,217,0.12); border-left: 4px solid #4a90d9;",
      tags$style(HTML(css)),
      tags$div(style = "display:flex; justify-content:space-between; align-items:baseline; gap:1rem;",
        tags$div(style = "font-size:1.4rem; font-weight:700;",
                 HTML(sprintf("⚡ Parallel MCCV running &mdash; %d/%d folds (%d%%)", done, n, pct))),
        tags$div(class = "text-muted",
                 HTML(sprintf("elapsed <b>%s</b> &middot; ~%s remaining",
                              format_duration(elapsed), format_duration(eta))))),
      tags$div(style = "height:12px; background:rgba(255,255,255,0.08); border-radius:6px; margin:12px 0; overflow:hidden;",
        tags$div(style = sprintf("height:12px; width:%d%%; background:linear-gradient(90deg,#4a90d9,#5fd0c0); border-radius:6px; transition:width .3s;", pct))),
      tags$div(style = "display:flex; gap:6px; flex-wrap:wrap; margin:6px 0 2px;",
        stat("Cores busy", sprintf("%d / %d", n_run, cores), "#5fd0c0"),
        stat("Completed", sprintf("%d / %d", done, n)),
        stat("Running", n_run, "#4a90d9"),
        stat("Queued", n_queue, "#aaaaaa"),
        stat("Mean acc", meanacc, "#2ea05a"),
        stat("ETA", format_duration(eta))),
      grid
    )
  })

  # Dismiss the persistent MCCV completion banner.
  observeEvent(input$mccv_dismiss_summary, {
    mc$summary_open <- FALSE
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

  # ── Live per-fold progress: reads each forked worker's stage file ─────────
  output$mccv_fold_progress <- renderUI({
    req(mc$active)
    invalidateLater(1000)                       # tick the elapsed clocks
    n <- mc$n %||% 0L
    if (n == 0) return(tags$p(class = "text-muted", "Preparing folds…"))
    now <- as.numeric(difftime(Sys.time(), mc$t0, units = "secs"))
    running_names <- names(mc$running)
    started <- mc$fold_started
    dir <- mc$progress_dir %||% ""
    fmt <- function(s) { s <- max(0L, as.integer(round(s))); sprintf("%d:%02d", s %/% 60L, s %% 60L) }
    read_stage <- function(i) {
      f <- file.path(dir, paste0("fold_", i, ".txt"))
      if (nzchar(dir) && file.exists(f)) {
        v <- tryCatch(readLines(f, warn = FALSE), error = function(e) character(0))
        if (length(v)) return(utils::tail(v, 1))
      }
      ""
    }
    rows <- lapply(seq_len(n), function(i) {
      done <- !is.null(mc$folds[[i]])
      if (done) {
        ok <- identical(mc$folds[[i]]$status, "success")
        icon <- if (ok) "✅" else "❌"
        txt  <- if (ok) "done" else paste0("failed: ", mc$folds[[i]]$error_message %||% "error")
        el   <- ""
      } else if (as.character(i) %in% running_names) {
        icon <- "⏳"
        stg  <- read_stage(i); if (!nzchar(stg)) stg <- "starting"
        st   <- started[[as.character(i)]] %||% now
        el   <- fmt(now - st)
        txt  <- stg
      } else {
        icon <- "•"; txt <- "queued"; el <- ""
      }
      col <- if (grepl("^failed", txt)) "#e06666" else "#aaaaaa"
      tags$div(style = "display:flex; justify-content:space-between; gap:10px; padding:3px 0; border-bottom:1px solid rgba(255,255,255,0.05);",
        tags$span(sprintf("%s Fold %d", icon, i)),
        tags$span(style = sprintf("text-align:right; color:%s;", col),
                  paste0(txt, if (nzchar(el)) paste0("  ·  ", el) else "")))
    })
    tagList(rows)
  })

  # ── Run configuration panel: exactly what this run used, with numbers ─────
  output$gsm_run_config <- renderUI({
    cfg <- rv$run_config
    if (is.null(cfg)) return(NULL)
    kv_rows <- function(lst) tagList(lapply(names(lst), function(k)
      tags$div(style = "display:flex; justify-content:space-between; gap:12px; padding:3px 0; border-bottom:1px solid rgba(255,255,255,0.05);",
        tags$span(class = "text-muted", style = "white-space:nowrap;", k),
        tags$b(style = "text-align:right;", as.character(lst[[k]])))))
    section <- function(icon, title, lst) tags$div(
      style = "flex:1 1 220px; min-width:220px;",
      tags$div(style = "font-weight:700; margin-bottom:4px;", paste(icon, title)),
      kv_rows(lst))
    tags$div(class = "card", style = "margin-bottom:12px; padding:12px 14px; background:rgba(74,144,217,0.06); border-left:4px solid #4a90d9;",
      tags$div(style = "display:flex; justify-content:space-between; align-items:center; margin-bottom:8px;",
        tags$div(style = "font-weight:700; font-size:1.02rem;", "\U0001F9FE This run — selected configuration"),
        tags$span(class = "text-muted", style = "font-size:0.82rem;", cfg$dataset)),
      tags$div(style = "display:flex; flex-wrap:wrap; gap:24px;",
        section("\U0001F9E9", "Grouping (G)", cfg$grouping),
        section("\U2702\UFE0F", "Gene shaving", cfg$shaving),
        section("\U0001F3AF", "Scoring / Modeling", cfg$scoring)))
  })

  # ── Groups created in the most recent completed fold (live during MCCV) ───
  output$gsm_live_groups <- renderDT({
    mc$tick
    ok <- Filter(function(r) !is.null(r) && identical(r$status, "success"),
                 isolate(mc$folds))
    validate(need(length(ok) > 0, "Groups will appear here as folds complete…"))
    f  <- ok[[length(ok)]]                      # most recent successful fold
    gs <- f$group_sizes %||% integer(0)
    acc <- f$group_accuracies %||% numeric(0)
    nm <- f$ranked_groups %||% names(gs)
    nm <- nm[nm %in% names(gs)]
    validate(need(length(nm) > 0, "No groups recorded for this fold."))
    df <- data.frame(
      Rank       = seq_along(nm),
      Group      = nm,
      Genes      = as.integer(gs[nm]),
      TrainScore = round(as.numeric(acc[nm]), 4),
      stringsAsFactors = FALSE)
    .dt(df, page = 10)
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
      req(rv$mccv_results)
      res <- rv$mccv_results
      if (!is.null(rv$detect)) res$detect <- rv$detect   # include discovered modules
      res
    } else {
      req(rv$score_results, rv$final_model)
      list(score_results = rv$score_results, final_model = rv$final_model)
    }
  })

  # ── G-S-M detect: de novo module discovery (on demand) ──────────────────
  observeEvent(input$run_detect, {
    req(rv$mccv_results, rv$dataset)
    session$sendCustomMessage("disableBtn", "run_detect")
    session$sendCustomMessage("showBanner", list(
      id = "detect_status",
      title = "Discovering robust gene modules… please wait.",
      detail = "Building stability network, clustering modules, and running bootstraps."
    ))
    rv$detect_running <- TRUE
    on.exit({ rv$detect_running <- FALSE; session$sendCustomMessage("enableBtn", "run_detect") })
    rv$detect <- NULL
    tryCatch({
        rv$detect <- gsm_detect(
          rv$mccv_results, rv$dataset, mode = input$detect_mode,
          top_genes = input$detect_top_genes, min_freq = input$detect_min_freq,
          resolution = input$detect_resolution, n_boot = input$detect_nboot,
          seed = 1)
        updateSelectInput(session, "detect_module_pick",
                          choices = rv$detect$modules$Module)
      },
      error = function(e) {
        rv$detect <- NULL
        showNotification(paste("❌ Discovery failed:", conditionMessage(e)),
                         type = "error", duration = 12)
      })
  })

  output$detect_status <- renderUI({
    if (isTRUE(rv$detect_running)) {
      return(tags$div(
        class = "value-box",
        style = "background: rgba(74,144,217,0.12); border-left:4px solid #4a90d9;",
        tags$div(style = "font-weight:600; font-size:1.1rem;",
                 "\u23f3 Discovering robust gene modules\u2026 please wait."),
        tags$div(class = "label", "Building stability network, clustering modules, and running bootstraps.")
      ))
    }
    req(rv$detect)
    p <- rv$detect$params
    has_mods <- nrow(rv$detect$modules) > 0
    tags$div(id = "detect_status_bar",
      class = "value-box",
      style = sprintf("background: %s; border-left:4px solid %s; position:relative;",
                      if (has_mods) "rgba(106,168,79,0.12)" else "rgba(241,194,50,0.12)",
                      if (has_mods) "#6aa84f" else "#f1c232"),
      tags$button(type = "button", class = "btn-close",
                  style = "position:absolute; top:8px; right:8px; filter:invert(1); opacity:0.6;",
                  onclick = "this.parentElement.style.display='none';",
                  `aria-label` = "Close"),
      tags$div(style = "font-weight:600; padding-right:30px;",
               sprintf("%s %d modules \u00b7 %d robust genes \u00b7 mode = %s \u00b7 %d iterations",
                       if (has_mods) "\u2705" else "\u26a0\ufe0f",
                       nrow(rv$detect$modules), p$n_robust, p$mode, p$n_runs)),
      if (!is.null(p$note)) tags$div(class = "label", p$note),
      if (isTRUE(p$truncated))
        tags$div(class = "label",
                 sprintf("Graph capped to top %d genes.", p$top_genes)))
  })

  output$detect_robust_table <- renderDT({ req(rv$detect); .dt(rv$detect$robust_genes) })
  output$detect_modules_table <- renderDT({
    req(rv$detect)
    validate(need(nrow(rv$detect$modules) > 0,
                  "No multi-gene modules — see the note above."))
    .dt(rv$detect$modules) |>
      formatStyle("Label",
                  backgroundColor = styleEqual(
                    c("rediscovered", "fragment", "novel"),
                    c("rgba(106,168,79,0.25)", "rgba(241,194,50,0.25)",
                      "rgba(224,102,102,0.25)")))
  })
  output$detect_module_detail <- renderDT({
    req(rv$detect, input$detect_module_pick)
    m <- rv$detect$members
    .dt(m[m$Module == input$detect_module_pick, , drop = FALSE], dom = "tp")
  })
  output$detect_network_vis <- visNetwork::renderVisNetwork({
    req(rv$detect)
    validate(need(requireNamespace("igraph", quietly = TRUE),
                  "Install the 'igraph' package to view the network."))
    validate(need(requireNamespace("visNetwork", quietly = TRUE),
                  "Install the 'visNetwork' package for the interactive network."))
    adj <- rv$detect$graph$adjacency
    memb <- rv$detect$graph$membership
    validate(need(nrow(adj) >= 2 && sum(adj) > 0, "Network has no edges to draw."))
    g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected",
                                             weighted = TRUE, diag = FALSE)
    # Build visNetwork data from igraph
    vn_data <- visNetwork::toVisNetworkData(g)
    # Color nodes by module membership
    n_modules <- max(memb, na.rm = TRUE)
    palette <- c("#888888", grDevices::hcl.colors(n_modules, "Dark 3"))
    # memb is a named vector (gene -> module integer); match safely
    node_ids <- vn_data$nodes$id
    node_memb <- memb[node_ids]
    node_memb[is.na(node_memb)] <- 0L  # unassigned genes
    vn_data$nodes$group <- paste0("M", sprintf("%02d", node_memb))
    vn_data$nodes$color <- palette[node_memb + 1L]
    vn_data$nodes$label <- node_ids
    vn_data$nodes$title <- paste0("<b>", node_ids, "</b><br>Module: ",
                                   vn_data$nodes$group)
    vn_data$nodes$font.color <- "#dddddd"
    vn_data$nodes$font.size <- 12
    # Scale edge width by weight
    if (!is.null(vn_data$edges$weight) && nrow(vn_data$edges) > 0) {
      max_w <- max(vn_data$edges$weight, na.rm = TRUE)
      vn_data$edges$width <- 1 + 4 * (vn_data$edges$weight / max_w)
    }
    vn_data$edges$color <- "rgba(74,144,217,0.4)"
    visNetwork::visNetwork(vn_data$nodes, vn_data$edges,
                           background = "transparent") |>
      visNetwork::visPhysics(solver = "forceAtlas2Based",
                             forceAtlas2Based = list(gravitationalConstant = -30,
                                                     springLength = 80),
                             stabilization = list(iterations = 150)) |>
      visNetwork::visInteraction(navigationButtons = TRUE,
                                  zoomView = TRUE, dragView = TRUE,
                                  hover = TRUE,
                                  tooltipDelay = 100) |>
      visNetwork::visOptions(highlightNearest = list(enabled = TRUE, degree = 1,
                                                      hover = TRUE),
                              nodesIdSelection = FALSE) |>
      visNetwork::visLayout(randomSeed = 42)
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
      tbls <- gsm_export_tables(assemble_result(), rv$dataset, rv$score_summary,
                                timings = rv$timings %||% list())
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
    base_dir <- getOption("gsmR.launch_dir", NULL)
    if (is.null(base_dir)) {
      if (file.exists("DESCRIPTION")) base_dir <- getwd()
      else if (file.exists("../../DESCRIPTION")) base_dir <- normalizePath("../..")
      else base_dir <- getwd()
    }
    out_dir <- file.path(base_dir, "gsm_outputs")
    tryCatch({
      dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
      tbls   <- gsm_export_tables(assemble_result(), rv$dataset, rv$score_summary,
                                  timings = rv$timings %||% list())
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

  output$has_timings <- reactive({
    !is.null(rv$timings) && length(rv$timings) > 0
  })
  outputOptions(output, "has_timings", suspendWhenHidden = FALSE)

  output$export_timings_table <- renderDT({
    req(rv$timings)
    tbl <- gsm_tbl_timings(rv$timings)
    req(!is.null(tbl) && nrow(tbl) > 0)
    .dt(tbl, page = 10)
  })

  # Single Fold Renderers: Group Scores Table
  output$score_table <- renderDT({
    req(rv$score_results)
    .members <- function(g) {
      genes <- rv$score_results[[g]]$genes %||% rv$prescore_groups[[g]] %||%
               (if (!is.null(rv$dataset)) rv$dataset$groups[[g]] else NULL)
      genes
    }
    gnames <- names(rv$score_results)
    df <- data.frame(
      Rank     = seq_along(rv$score_results),
      Group    = gnames,
      Accuracy = vapply(rv$score_results, \(x) round(x$accuracy, 4), numeric(1)),
      Genes    = vapply(gnames, \(g) length(.members(g)), integer(1)),
      Members  = vapply(gnames, \(g) {
        gg <- .members(g); if (is.null(gg)) return("—")
        s <- paste(utils::head(gg, 20), collapse = ", ")
        if (length(gg) > 20) paste0(s, ", … (", length(gg), " total)") else s
      }, character(1)),
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
    # Best-k metric summary with mean, SD and bootstrap-style 95% CI (incl. AUC).
    df <- gsmR:::gsm_tbl_metric_summary(rv$mccv_results)
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
  # ── RRA + per-iteration tables for the results dashboard ────────────────
  .mccv_folds <- reactive({ req(rv$mccv_results); rv$mccv_results$folds })
  .dt <- function(df, page = 15, dom = "ftp")
    datatable(df, options = list(pageLength = page, dom = dom, scrollX = TRUE),
              rownames = FALSE, class = "compact stripe")

  # ══ Grouping (G) tab ════════════════════════════════════════════════════
  # Shared renderers for the "run a grouping method and inspect it" surface.
  .g_valueboxes <- function(groups, expr) {
    sizes <- vapply(groups, length, integer(1))
    genes <- unique(unlist(groups))
    cov <- if (nrow(expr) > 0)
      round(100 * length(intersect(genes, rownames(expr))) / nrow(expr), 1) else 0
    layout_columns(
      col_widths = c(3, 3, 3, 3),
      value_box_ui(length(groups), "Groups", "\U0001F9E9"),
      value_box_ui(format(length(genes), big.mark = ","), "Unique genes", "\U0001F9EC"),
      value_box_ui(if (length(sizes)) stats::median(sizes) else 0, "Median size", "\U0001F4CF"),
      value_box_ui(paste0(cov, "%"), "Gene coverage", "\U0001F4CA")
    )
  }
  .g_groups_df <- function(groups, scores = NULL) {
    if (length(groups) == 0) return(data.frame(Group = character(0)))
    df <- data.frame(
      Group = names(groups),
      Genes = vapply(groups, length, integer(1)),
      Members = vapply(groups, function(g) {
        s <- paste(utils::head(g, 12), collapse = ", ")
        if (length(g) > 12) paste0(s, ", …") else s
      }, character(1)),
      row.names = NULL, stringsAsFactors = FALSE)
    if (!is.null(scores) && length(scores) == nrow(df)) {
      df <- cbind(df[1:2], Score = round(as.numeric(scores), 4), Members = df$Members)
      df <- df[order(-df$Score), , drop = FALSE]
    } else {
      df <- df[order(-df$Genes), , drop = FALSE]
    }
    df
  }
  .g_size_hist <- function(groups, col = "#4a90d9") {
    sizes <- vapply(groups, length, integer(1))
    par(bg = "transparent", col.axis = "#aaaaaa", col.lab = "#aaaaaa",
        col.main = "#cccccc", fg = "#555555", mar = c(4, 4, 2, 1))
    if (length(sizes) == 0) { plot.new(); text(0.5, 0.5, "No groups yet.", col = "#aaaaaa"); return(invisible()) }
    hist(sizes, breaks = max(5, min(30, length(sizes))), col = col, border = NA,
         main = "", xlab = "Genes per group", ylab = "Frequency", las = 1)
    abline(v = stats::median(sizes), col = "#f1c232", lwd = 2, lty = 2)
  }
  .g_expr <- function() rv$expr_data$expr

  # Resolve an example grouping-file selection to a named list (reuses Data-tab
  # example dir + toy fallback).
  .g_read_example_groups <- function(sel) {
    if (is.null(sel) || !nzchar(sel))
      stop("Pick an example grouping file.")
    path <- if (identical(sel, "__toy__"))
      system.file("extdata", "test_groups.csv", package = "gsmR")
      else file.path(.example_groups_dir, sel)
    if (!file.exists(path)) stop("Example grouping file not found: ", sel)
    if (identical(tolower(tools::file_ext(path)), "gmt")) read_gmt(path)
    else read_groups(path)
  }

  # ── Step A: load expression data for grouping (shared with the pipeline) ──
  observeEvent(input$g_load_data, {
    src_upload <- !is.null(input$g_expr_upload)
    sel <- input$g_example_expr %||% ""
    if (!src_upload && !nzchar(sel))
      return(showNotification("⚠️ Pick an example dataset or upload an expression CSV.",
                              type = "warning"))
    path <- if (src_upload) input$g_expr_upload$datapath
            else if (identical(sel, "__toy__"))
              system.file("extdata", "test_expression.csv", package = "gsmR")
            else file.path(.example_input_dir, sel)
    if (is.null(path) || !nzchar(path) || !file.exists(path))
      return(showNotification("❌ Expression file not found.", type = "error", duration = 8))
    t0_data <- Sys.time()
    tryCatch({
      ed <- withProgress(message = "Loading expression…", value = 0.4,
                         read_expression(path))
      rv$expr_data <- ed
      elapsed_data <- as.numeric(difftime(Sys.time(), t0_data, units = "secs"))
      timings_cur <- rv$timings %||% list()
      timings_cur$data_load <- elapsed_data
      timings_cur$data_load_detail <- if (src_upload) "Uploaded expression (Grouping tab)" else "Load Selected Example (Grouping tab)"
      rv$timings <- timings_cur
      rv$g_file <- NULL; rv$g_bridge <- NULL           # clear stale results
      showNotification(sprintf("✅ Loaded %d samples × %d genes.",
                               ncol(ed$expr), nrow(ed$expr)), type = "message")
    }, error = function(e)
      showNotification(paste("❌", e$message), type = "error", duration = 10))
  })
  output$g_data_status <- renderUI({
    if (is.null(rv$expr_data))
      return(tags$div(class = "text-warning", style = "font-size:0.82rem;",
                      "No expression data yet — load one below or on the Data tab."))
    ed <- rv$expr_data
    tags$div(class = "text-success", style = "font-size:0.82rem;",
             sprintf("✅ %d samples × %d genes, %d classes", ncol(ed$expr),
                     nrow(ed$expr), nlevels(ed$labels)))
  })

  # ── Prior-knowledge (file) method: run + render ──────────────────────────
  observeEvent(input$g_file_run, {
    req(rv$expr_data)
    src <- input$g_file_source %||% "example"
    t0_gload <- Sys.time()
    tryCatch({
      groups <- if (src == "example") {
        .g_read_example_groups(input$g_file_example %||% "")
      } else {
        if (is.null(input$g_file_upload))
          stop("Upload a grouping file, or choose an example.")
        compute_grouping("file", file = input$g_file_upload$datapath,
                         format = input$g_file_format)
      }
      rv$g_file <- list(groups = groups)
      elapsed_gload <- as.numeric(difftime(Sys.time(), t0_gload, units = "secs"))
      timings_cur <- rv$timings %||% list()
      timings_cur$grouping_load <- elapsed_gload
      rv$timings <- timings_cur
      showNotification(sprintf("✅ Loaded %d group(s).", length(groups)), type = "message")
    }, error = function(e)
      showNotification(paste("❌", e$message), type = "error", duration = 10))
  })
  output$g_file_valueboxes <- renderUI({ req(rv$g_file); .g_valueboxes(rv$g_file$groups, .g_expr()) })
  output$g_file_groups_dt   <- renderDT({ req(rv$g_file); .dt(.g_groups_df(rv$g_file$groups)) })
  output$g_file_size_hist   <- renderPlot({ req(rv$g_file); .g_size_hist(rv$g_file$groups) }, bg = "transparent")

  # ── BRIDGE method: run (all samples, exploratory) + rich render ──────────
  # BRIDGE grouping runs in a BACKGROUND FORK so the UI stays responsive and the
  # Stop button works (a synchronous run would block Shiny and be uncancellable).
  observeEvent(input$g_bridge_run, {
    if (is.null(rv$expr_data)) {
      showNotification("⚠️ Load expression data first (① in the sidebar, or the Data tab).",
                       type = "warning", duration = 8); return()
    }
    if (isTRUE(gbr$active)) {
      showNotification("⏳ BRIDGE is already running — Stop it first.", type = "warning")
      return()
    }
    # Resolve the prior source (aborts with a clear message on any problem).
    src <- input$g_bridge_prior_source %||% "none"
    prior_obj <- NULL; ok <- TRUE
    if (src == "example") {
      sel <- input$g_bridge_prior_example %||% ""
      prior_obj <- tryCatch(.g_read_example_groups(sel),
        error = function(e) { showNotification(paste("❌", e$message),
                                type = "error", duration = 10); ok <<- FALSE; NULL })
    } else if (src == "loaded") {
      if (is.null(rv$groups) || length(rv$groups) == 0) {
        showNotification("❌ No grouping loaded on the Data tab. Pick 'Example', 'Upload', or 'None'.",
                         type = "error", duration = 10); ok <- FALSE
      } else prior_obj <- rv$groups
    } else if (src == "upload") {
      if (is.null(input$g_bridge_prior_file)) {
        showNotification("❌ Upload a prior file, or choose a different source.",
                         type = "error", duration = 10); ok <- FALSE
      } else prior_obj <- tryCatch(
        bridge_read_prior_file(input$g_bridge_prior_file$datapath),
        error = function(e) { showNotification(paste("❌", e$message),
                                type = "error", duration = 10); ok <<- FALSE; NULL })
    }
    if (!ok) return()

    args <- build_bridge_args(input, "g_bridge",
                              input$g_bridge_min_size, input$g_bridge_max_size)
    # Plain objects captured in the PARENT (never touch reactives inside a fork).
    X <- rv$expr_data$expr; y <- rv$expr_data$labels
    pfile <- tempfile("gbridge_prog_")   # plain string; fork writes stage markers
    call_args <- c(list(X = X, y = y, B = prior_obj, details = TRUE,
                        progress_file = pfile), args)
    rv$g_bridge <- NULL
    gbr$prior <- prior_obj; gbr$args <- args; gbr$progress_file <- pfile
    # Capture estimator label and start time for timing records.
    gbr$estimator <- input$g_bridge_precision %||% "unknown"
    gbr$t0 <- Sys.time(); gbr$tick <- 0; gbr$active <- TRUE

    if (.Platform$OS.type != "windows") {
      gbr$job <- parallel::mcparallel(
        tryCatch(do.call(bridge_grouping, call_args),
                 error = function(e) structure(list(), bridge_error = conditionMessage(e))),
        name = "g_bridge")
    } else {
      # No fork on Windows: run synchronously (blocking; Stop unavailable).
      res <- tryCatch(do.call(bridge_grouping, call_args),
                      error = function(e) structure(list(), bridge_error = conditionMessage(e)))
      gbr$active <- FALSE
      .g_bridge_finish(res, prior_obj, args)
    }
  })

  # Finalize a completed BRIDGE result (shared by fork poller + Windows path).
  .g_bridge_finish <- function(res, prior_obj, args) {
    err <- attr(res, "bridge_error")
    if (!is.null(err))
      return(showNotification(paste("❌", err), type = "error", duration = 12))
    # Record timing under the estimator label.
    elapsed_bridge <- as.numeric(difftime(Sys.time(), gbr$t0, units = "secs"))
    estimator_lbl  <- gbr$estimator %||% "unknown"
    timings_cur <- rv$timings %||% list()
    bridge_runs <- timings_cur$bridge_runs %||% list()
    bridge_runs[[estimator_lbl]] <- elapsed_bridge
    timings_cur$bridge_runs <- bridge_runs
    rv$timings <- timings_cur
    if (length(res) == 0)
      return(showNotification("BRIDGE produced no groups — lower tau/min size or adjust lambda.",
                              type = "warning", duration = 10))
    rv$g_bridge <- list(
      groups = res, diagnostics = attr(res, "diagnostics"),
      summary = attr(res, "summary"), graph = attr(res, "graph"),
      spec = list(method = "bridge", prior = prior_obj, args = args))
    showNotification(sprintf("✅ BRIDGE derived %d module(s).", length(res)), type = "message")
  }

  # Poller: reap the forked BRIDGE job without blocking the UI.
  observe({
    req(gbr$active)
    invalidateLater(600)
    isolate({
      if (is.null(gbr$job)) return()
      got <- tryCatch(parallel::mccollect(gbr$job, wait = FALSE, timeout = 0),
                      error = function(e) NULL)
      if (!is.null(got)) {
        res <- got[[1]]
        gbr$active <- FALSE; gbr$job <- NULL
        if (is.null(res)) {
          showNotification("❌ BRIDGE worker crashed.", type = "error", duration = 10)
        } else {
          .g_bridge_finish(res, gbr$prior, gbr$args)
        }
      }
    })
  })

  # Stop a running BRIDGE fork.
  observeEvent(input$g_bridge_stop, {
    req(gbr$active)
    if (!is.null(gbr$job)) {
      pid <- tryCatch(gbr$job$pid, error = function(e) NULL)
      if (!is.null(pid)) tryCatch(tools::pskill(pid, tools::SIGKILL), error = function(e) NULL)
      tryCatch(parallel::mccollect(gbr$job, wait = FALSE), error = function(e) NULL)
    }
    gbr$active <- FALSE; gbr$job <- NULL
    showNotification("⏹️ BRIDGE run stopped.", type = "warning", duration = 6)
  })

  # Live "running" banner: stage-based progress bar (steps 1..5) + elapsed clock.
  output$g_bridge_progress <- renderUI({
    req(gbr$active); invalidateLater(800)
    el <- as.integer(round(as.numeric(difftime(Sys.time(), gbr$t0, units = "secs"))))
    step <- 1L; label <- "starting…"
    pf <- gbr$progress_file %||% ""
    if (nzchar(pf) && file.exists(pf)) {
      v <- tryCatch(readLines(pf, warn = FALSE), error = function(e) character(0))
      if (length(v)) {
        parts <- strsplit(utils::tail(v, 1), "|", fixed = TRUE)[[1]]
        step <- suppressWarnings(as.integer(parts[1])); if (is.na(step)) step <- 1L
        label <- if (length(parts) > 1) parts[2] else "working…"
      }
    }
    pct <- max(8, min(100, round(step / 5 * 100)))
    tags$div(class = "value-box",
      style = "background: rgba(244,166,35,0.12); border-left:4px solid #f4a623;",
      tags$div(style = "display:flex; justify-content:space-between; align-items:center; margin-bottom:8px;",
        tags$div(tags$b("⏳ Running BRIDGE grouping"),
                 tags$div(class = "text-muted", style = "font-size:0.82rem;",
                          sprintf("Step %d/5 · %s", step, label))),
        tags$div(style = "display:flex; align-items:center; gap:10px;",
          tags$span(class = "text-muted", style = "font-size:0.82rem;",
                    sprintf("%d:%02d", el %/% 60L, el %% 60L)),
          actionButton("g_bridge_stop", "⏹️ Stop", class = "btn-danger btn-sm"))),
      tags$div(class = "progress", style = "height:14px; background:rgba(255,255,255,0.08);",
        tags$div(class = "progress-bar progress-bar-striped progress-bar-animated",
                 role = "progressbar",
                 style = sprintf("width:%d%%; background-color:#f4a623;", pct),
                 `aria-valuenow` = pct, `aria-valuemin` = 0, `aria-valuemax` = 100)),
      tags$div(class = "text-muted", style = "font-size:0.75rem; margin-top:4px;",
               "Runs on all samples (exploratory). The 'estimating networks' step is the slowest; a stronger pre-filter makes it faster."))
  })
  output$g_bridge_valueboxes <- renderUI({ req(rv$g_bridge); .g_valueboxes(rv$g_bridge$groups, .g_expr()) })
  output$g_bridge_size_hist  <- renderPlot({ req(rv$g_bridge); .g_size_hist(rv$g_bridge$groups, "#6aa84f") }, bg = "transparent")
  output$g_bridge_modules_dt <- renderDT({
    req(rv$g_bridge)
    d <- rv$g_bridge$diagnostics
    validate(need(!is.null(d) && nrow(d) > 0, "No modules."))
    gr <- rv$g_bridge$groups
    d$genes <- vapply(d$module, function(m) {
      g <- gr[[m]]; s <- paste(utils::head(g, 10), collapse = ", ")
      if (length(g) > 10) paste0(s, ", …") else s
    }, character(1))
    d[2:6] <- lapply(d[2:6], function(x) if (is.numeric(x)) round(x, 4) else x)
    .dt(d[order(-d$J), c("module", "size", "J", "D", "C", "stability", "genes")])
  })
  output$g_bridge_diag_card <- renderUI({
    req(rv$g_bridge); s <- rv$g_bridge$summary
    if (is.null(s)) return(NULL)
    chip <- function(lab, val) tags$div(
      style = "flex:1 1 130px; text-align:center; padding:8px; background:rgba(255,255,255,0.04); border-radius:8px;",
      tags$div(class = "text-muted", style = "font-size:.7rem; text-transform:uppercase;", lab),
      tags$div(style = "font-size:1.05rem; font-weight:700;", val))
    sweep <- paste(sprintf("r%s: %d", names(s$res_counts), s$res_counts), collapse = "  ·  ")
    tags$div(class = "card", style = "margin-bottom:12px; padding:12px;",
      tags$div(style = "font-weight:700; margin-bottom:8px;", "🔬 BRIDGE diagnostics"),
      tags$div(style = "display:flex; flex-wrap:wrap; gap:10px;",
        chip("Pre-filter kept", sprintf("%d / %d", s$prefilter_kept, s$n_features_total %||% s$prefilter_kept)),
        chip("Est. features", s$n_features),
        chip("Backend", s$backend),
        chip("Candidates → selected", sprintf("%d → %d", s$n_candidates, nrow(rv$g_bridge$diagnostics))),
        chip("Prior genes / edges", sprintf("%d / %d", s$prior_genes, s$prior_edges))),
      tags$div(class = "text-muted", style = "margin-top:8px; font-size:0.82rem;",
               paste0("Modules per resolution — ", sweep)))
  })
  output$g_bridge_network <- visNetwork::renderVisNetwork({
    req(rv$g_bridge)
    validate(need(requireNamespace("visNetwork", quietly = TRUE),
                  "Install 'visNetwork' for the interactive network."))
    edges <- rv$g_bridge$graph
    validate(need(!is.null(edges) && nrow(edges) > 0, "No differential-network edges to show."))
    # Map each gene to its (first) selected module for colouring.
    gr <- rv$g_bridge$groups
    g2m <- stats::setNames(rep(names(gr), lengths(gr)), unlist(gr))
    nodes_id <- union(edges$from, edges$to)
    nodes <- data.frame(id = nodes_id, label = nodes_id,
                        group = g2m[nodes_id] %||% "—", stringsAsFactors = FALSE)
     e <- data.frame(from = edges$from, to = edges$to,
                     value = edges$weight, stringsAsFactors = FALSE)
    visNetwork::visNetwork(nodes, e) |>
      visNetwork::visIgraphLayout(layout = "layout_with_fr") |>
      visNetwork::visOptions(highlightNearest = list(enabled = TRUE, degree = 1),
                             nodesIdSelection = FALSE) |>
      visNetwork::visInteraction(tooltipDelay = 100) |>
      visNetwork::visPhysics(stabilization = TRUE)
  })

  # ── Commit a grouping as the active pipeline dataset ─────────────────────
  # Generic: build rv$dataset from the loaded expression + (groups | spec).
  .commit_dataset <- function(groups, spec, label) {
    if (is.null(rv$expr_data)) {
      showNotification("⚠️ Load expression data on the Data tab first.",
                       type = "warning"); return(invisible(FALSE)) }
    ed <- rv$expr_data
    ds <- tryCatch(
      gsm_dataset(ed$expr, ed$labels, groups %||% list(),
                  min_group_size = 5, max_group_size = 500, grouping_spec = spec),
      error = function(e) { showNotification(paste("❌", e$message), type = "error", duration = 10); NULL })
    if (is.null(ds)) return(invisible(FALSE))
    rv$dataset <- ds; rv$groups <- ds$groups; rv$loaded <- TRUE
    rv$dataset_name <- rv$dataset_name %||% "dataset"
    rv$collapsed <- FALSE; rv$prescored <- FALSE; rv$scored <- FALSE
    rv$prescore_groups <- NULL; rv$score_results <- NULL; rv$mccv_results <- NULL
    showNotification(sprintf("✅ Using %s for scoring.", label), type = "message")
    invisible(TRUE)
  }

  # Commit a prepared G-tab result (used by the Grouping-tab "Apply" button).
  .commit_grouping <- function(which) {
    if (identical(which, "g_bridge")) {
      if (is.null(rv$g_bridge)) return(showNotification("⚠️ Run BRIDGE first.", type = "warning"))
      .commit_dataset(list(), rv$g_bridge$spec, "BRIDGE (per-fold)")
    } else {
      if (is.null(rv$g_file)) return(showNotification("⚠️ Run a Prior-knowledge grouping first.", type = "warning"))
      .commit_dataset(rv$g_file$groups, NULL,
                      sprintf("prior-knowledge (%d groups)", length(rv$g_file$groups)))
    }
  }
  observeEvent(input$g_apply, {
    .commit_grouping(input$g_method_tabs %||% "g_file")
  })
  output$g_active_status <- renderUI({
    val <- if (is.null(rv$dataset)) "— none applied yet —"
    else {
      gm <- rv$dataset$params$grouping_method %||% "file"
      n  <- rv$dataset$params$n_groups
      if (gm == "bridge" && n == 0) "BRIDGE (per-fold)"
      else sprintf("%s — %d groups", if (gm == "bridge") "BRIDGE" else "file", n)
    }
    tags$div(style = "margin-top:8px;",
      tags$div(class = "text-muted", style = "font-size:0.78rem; text-transform:uppercase;",
               "Active pipeline grouping"),
      tags$div(style = "font-weight:700; font-size:1.0rem; color:#4a90d9;", val))
  })

  # ── GSM tab: configure + run the grouping to score, right here ───────────
  observeEvent(input$gsm_grp_run, {
    if (is.null(rv$expr_data))
      return(showNotification("⚠️ Load expression data on the Data tab first.", type = "warning"))
    method <- input$gsm_grp_method %||% "file"

    if (method == "file") {
      t0_gload <- Sys.time()
      groups <- tryCatch({
        if ((input$gsm_grp_file_source %||% "example") == "example")
          .g_read_example_groups(input$gsm_grp_file_example %||% "")
        else {
          if (is.null(input$gsm_grp_file_upload)) stop("Upload a grouping file, or choose an example.")
          compute_grouping("file", file = input$gsm_grp_file_upload$datapath, format = "auto")
        }
      }, error = function(e) { showNotification(paste("❌", e$message), type = "error", duration = 10); NULL })
      if (is.null(groups)) return()
      elapsed_gload <- as.numeric(difftime(Sys.time(), t0_gload, units = "secs"))
      timings_cur <- rv$timings %||% list()
      timings_cur$grouping_load <- elapsed_gload
      rv$timings <- timings_cur
      .commit_dataset(groups, NULL, sprintf("prior-knowledge (%d groups)", length(groups)))

    } else {  # BRIDGE — commit the deferred spec (re-derived per fold)
      src <- input$gsm_grp_bridge_prior %||% "none"
      prior <- tryCatch(switch(src,
        none    = NULL,
        example = .g_read_example_groups(input$gsm_grp_bridge_prior_example %||% ""),
        upload  = {
          if (is.null(input$gsm_grp_bridge_prior_file)) stop("Upload a prior file, or choose a different source.")
          bridge_read_prior_file(input$gsm_grp_bridge_prior_file$datapath)
        }),
        error = function(e) { showNotification(paste("❌", e$message), type = "error", duration = 10); "ERR" })
      if (identical(prior, "ERR")) return()
      args <- build_bridge_args(input, "gsm_grp_bridge", 5, 500)
      spec <- list(method = "bridge", prior = prior, args = args)
      .commit_dataset(list(), spec, "BRIDGE (re-derived per training fold)")
    }
  })

  output$gsm_grouping_note <- renderUI({
    if (is.null(rv$dataset))
      return(tags$div(class = "text-warning", style = "font-size:0.82rem; margin-top:6px;",
                      "No grouping selected yet — configure one above and click 'Use this grouping'."))
    gm <- rv$dataset$params$grouping_method %||% "file"
    n  <- rv$dataset$params$n_groups
    txt <- if (gm == "bridge" && n == 0) "BRIDGE — re-derived per training fold (leakage-free)"
           else sprintf("%s — %d groups", if (gm == "bridge") "BRIDGE" else "Prior-knowledge", n)
    tags$div(class = "text-success", style = "font-size:0.85rem; margin-top:6px; font-weight:600;",
             paste0("✅ Active grouping: ", txt))
  })

  output$rra_group_table <- renderDT({
    tb <- gsmR:::gsm_tbl_group_rra(.mccv_folds())
    validate(need(!is.null(tb),
      if (identical(rv$dataset$params$grouping_method, "bridge"))
        "BRIDGE re-derives groups each fold, so group names are fold-local — see the gene-level consensus (RRA over genes) below."
      else
        "Run Parallel MCCV (≥ 2 iterations) to compute RRA over groups."))
    .dt(tb)
  })
  output$rra_feature_table <- renderDT({
    tb <- gsmR:::gsm_tbl_feature_rra(.mccv_folds())
    validate(need(!is.null(tb), "RRA over genes needs Random Forest importance from MCCV folds."))
    .dt(tb)
  })
  output$ranked_features_table <- renderDT({
    tb <- gsmR:::gsm_tbl_ranked_features(.mccv_folds(), rv$dataset)
    validate(need(!is.null(tb), "No ranked features yet."))
    .dt(tb)
  })

  # ── Consensus biomarkers (gsm_biomarkers) ──────────────────────────────
  output$biomarker_table <- renderDT({
    req(rv$mccv_results)
    bm <- gsm_biomarkers(rv$mccv_results)
    validate(need(!is.null(bm),
                  "Consensus biomarkers need per-gene importance — run MCCV with Random Forest, linear-kernel SVM, or XGBoost scoring/modeling."))
    .dt(bm)
  })

  # ── Group significance: permutation test (Pre-Score tab) ───────────────
  observeEvent(input$btn_permtest, {
    if (is.null(rv$dataset)) {
      showNotification("⚠️ Load a dataset first.", type = "warning", duration = 6)
      return()
    }
    m <- input$perm_method
    need_pkg <- switch(m, rf = "randomForest", svm = "e1071",
                       knn = "class", xgb = "xgboost", NULL)
    if (!is.null(need_pkg) && !requireNamespace(need_pkg, quietly = TRUE)) {
      showNotification(sprintf("❌ Package '%s' required for this classifier.", need_pkg),
                       type = "error", duration = 10)
      return()
    }
    sargs <- switch(m, rf = list(ntree = 100), xgb = list(nrounds = 100), list())
    withProgress(message = "Permutation testing groups...", value = 0.05, {
      res <- tryCatch(
        permutation_test_groups(
          rv$dataset, method = m, n_perm = input$perm_nperm,
          score_args = sargs, seed = input$prescore_seed,
          progress = function(i, n, g) setProgress(i / n, detail = g)),
        error = function(e) {
          showNotification(paste("Permutation test failed:", conditionMessage(e)),
                           type = "error", duration = 12)
          NULL
        })
      rv$permtest <- res
    })
    if (!is.null(rv$permtest))
      showNotification("✅ Group significance computed.", type = "message", duration = 5)
  })

  output$permtest_table <- renderDT({
    req(rv$permtest)
    df <- rv$permtest
    for (cl in intersect(c("observed", "null_mean", "p_value", "p_adj"), names(df)))
      df[[cl]] <- signif(df[[cl]], 3)
    .dt(df)
  })
  output$iter_metrics_table <- renderDT({
    tb <- gsmR:::gsm_tbl_iteration_metrics(.mccv_folds())
    validate(need(!is.null(tb), "No iterations yet."))
    .dt(tb)
  })
  output$iter_accum_table <- renderDT({
    tb <- gsmR:::gsm_tbl_iteration_accumulative(.mccv_folds())
    validate(need(!is.null(tb), "No iterations yet."))
    .dt(tb)
  })
  output$iter_accum_table_single <- renderDT({
    req(rv$accumulative)
    .dt(rv$accumulative, dom = "tp")
  })
  output$single_feature_table <- renderDT({
    req(rv$final_model)
    feats <- rv$final_model$features %||% character(0)
    validate(need(length(feats) > 0, "No features in the model."))
    imp <- rv$final_model$importance
    df <- data.frame(
      Gene = feats,
      Importance = if (!is.null(imp)) round(as.numeric(imp[feats]), 5) else NA_real_,
      stringsAsFactors = FALSE)
    if (!all(is.na(df$Importance))) df <- df[order(-df$Importance), ]
    .dt(df)
  })

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
      if (is.null(accs) || length(accs) == 0) next
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
    
    # Member genes per group. For fixed (prior-knowledge) groupings the group
    # name maps to a constant gene set; for BRIDGE the groups are re-derived per
    # fold (fold-local names), so a single member list is not defined.
    grp_members <- if (!is.null(rv$dataset) && length(rv$dataset$groups) > 0)
      rv$dataset$groups else NULL

    summary_list <- list()
    for (gname in names(group_performance)) {
      perf <- group_performance[[gname]]
      genes <- if (!is.null(grp_members)) grp_members[[gname]] else NULL
      gene_txt <- if (!is.null(genes)) {
        s <- paste(utils::head(genes, 20), collapse = ", ")
        if (length(genes) > 20) paste0(s, ", … (", length(genes), " total)") else s
      } else "(BRIDGE: re-derived per fold — see the gene consensus table)"
      summary_list[[gname]] <- data.frame(
        Group = gname,
        N_genes = if (!is.null(genes)) length(genes) else NA_integer_,
        Folds_Evaluated = length(perf$accuracies),
        Top_10_Appearances = perf$in_top_10,
        Appearance_Rate = round((perf$in_top_10 / n_success) * 100, 1),
        Mean_Accuracy = round(mean(perf$accuracies) * 100, 2),
        Mean_Rank = round(mean(perf$ranks), 1),
        Genes = gene_txt,
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

  # ══════════════════════════════════════════════════════════════════════
  # Settings Import / Export
  # ══════════════════════════════════════════════════════════════════════
  
  # Import Settings
  observeEvent(input$import_settings, {
    req(input$import_settings)
    tryCatch({
      settings <- gsmR::read_settings(input$import_settings$datapath)
      
      # Update UI inputs based on parsed YAML
      if (!is.null(settings$grouping_method)) updateSelectInput(session, "grouping_method", selected = settings$grouping_method)
      if (!is.null(settings$groups_format)) updateSelectInput(session, "groups_format", selected = settings$groups_format)
      # BRIDGE hyperparameters
      if (!is.null(settings$bridge_prior_source)) updateSelectInput(session, "bridge_prior_source", selected = settings$bridge_prior_source)
      if (!is.null(settings$bridge_data_type)) updateSelectInput(session, "bridge_data_type", selected = settings$bridge_data_type)
      if (!is.null(settings$bridge_lambda)) updateNumericInput(session, "bridge_lambda", value = as.numeric(settings$bridge_lambda))
      if (!is.null(settings$bridge_beta)) updateNumericInput(session, "bridge_beta", value = as.numeric(settings$bridge_beta))
      if (!is.null(settings$bridge_lambda2)) updateNumericInput(session, "bridge_lambda2", value = as.numeric(settings$bridge_lambda2))
      if (!is.null(settings$bridge_lambda1)) updateNumericInput(session, "bridge_lambda1", value = as.numeric(settings$bridge_lambda1))
      if (!is.null(settings$bridge_resolution)) updateTextInput(session, "bridge_resolution", value = settings$bridge_resolution)
      if (!is.null(settings$bridge_tau)) updateSliderInput(session, "bridge_tau", value = as.numeric(settings$bridge_tau))
      if (!is.null(settings$bridge_prefilter)) updateSelectInput(session, "bridge_prefilter", selected = settings$bridge_prefilter)
      if (!is.null(settings$bridge_precision)) updateSelectInput(session, "bridge_precision", selected = settings$bridge_precision)
      if (!is.null(settings$bridge_max_genes)) updateNumericInput(session, "bridge_max_genes", value = as.numeric(settings$bridge_max_genes))
      if (!is.null(settings$bridge_preview)) updateCheckboxInput(session, "bridge_preview", value = as.logical(settings$bridge_preview))
      if (!is.null(settings$min_group_size)) updateNumericInput(session, "min_group_size", value = as.numeric(settings$min_group_size))
      if (!is.null(settings$max_group_size)) updateNumericInput(session, "max_group_size", value = as.numeric(settings$max_group_size))
      
      if (!is.null(settings$collapse_method)) updateSelectInput(session, "collapse_method", selected = settings$collapse_method)
      if (!is.null(settings$meta_conflict)) updateSelectInput(session, "meta_conflict", selected = settings$meta_conflict)
      if (!is.null(settings$keep_original)) updateCheckboxInput(session, "keep_original", value = as.logical(settings$keep_original))
      
      if (!is.null(settings$prescore_method)) updateSelectInput(session, "prescore_method", selected = settings$prescore_method)
      if (!is.null(settings$top_percent)) updateSliderInput(session, "top_percent", value = as.numeric(settings$top_percent))
      if (!is.null(settings$train_ratio)) updateSliderInput(session, "train_ratio", value = as.numeric(settings$train_ratio))
      if (!is.null(settings$prescore_seed)) updateNumericInput(session, "prescore_seed", value = as.numeric(settings$prescore_seed))
      
      if (!is.null(settings$perm_method)) updateSelectInput(session, "perm_method", selected = settings$perm_method)
      if (!is.null(settings$perm_nperm)) updateNumericInput(session, "perm_nperm", value = as.numeric(settings$perm_nperm))
      
      if (!is.null(settings$score_method)) updateSelectInput(session, "score_method", selected = settings$score_method)
      if (!is.null(settings$rf_ntree)) updateNumericInput(session, "rf_ntree", value = as.numeric(settings$rf_ntree))
      if (!is.null(settings$svm_kernel)) updateSelectInput(session, "svm_kernel", selected = settings$svm_kernel)
      if (!is.null(settings$svm_cost)) updateNumericInput(session, "svm_cost", value = as.numeric(settings$svm_cost))
      if (!is.null(settings$svm_gamma)) updateNumericInput(session, "svm_gamma", value = as.numeric(settings$svm_gamma))
      if (!is.null(settings$knn_k)) updateNumericInput(session, "knn_k", value = as.numeric(settings$knn_k))
      
      if (!is.null(settings$shave_genes)) updateCheckboxInput(session, "shave_genes", value = as.logical(settings$shave_genes))
      if (!is.null(settings$shave_method)) updateSelectInput(session, "shave_method", selected = settings$shave_method)
      if (!is.null(settings$shave_k)) updateNumericInput(session, "shave_k", value = as.numeric(settings$shave_k))
      if (!is.null(settings$fdr_alpha)) updateNumericInput(session, "fdr_alpha", value = as.numeric(settings$fdr_alpha))
      if (!is.null(settings$normalize_z)) updateCheckboxInput(session, "normalize_z", value = as.logical(settings$normalize_z))
      
      if (!is.null(settings$score_scope)) updateRadioButtons(session, "score_scope", selected = settings$score_scope)
      if (!is.null(settings$n_top_groups)) updateNumericInput(session, "n_top_groups", value = as.numeric(settings$n_top_groups))
      if (!is.null(settings$score_seed)) updateNumericInput(session, "score_seed", value = as.numeric(settings$score_seed))
      
      if (!is.null(settings$run_mccv_mode)) updateCheckboxInput(session, "run_mccv_mode", value = as.logical(settings$run_mccv_mode))
      if (!is.null(settings$mccv_iterations)) updateNumericInput(session, "mccv_iterations", value = as.numeric(settings$mccv_iterations))
      if (!is.null(settings$model_method)) updateSelectInput(session, "model_method", selected = settings$model_method)
      if (!is.null(settings$mccv_train_ratio)) updateSliderInput(session, "mccv_train_ratio", value = as.numeric(settings$mccv_train_ratio))
      if (!is.null(settings$mccv_cores)) updateSliderInput(session, "mccv_cores", value = as.numeric(settings$mccv_cores))
      
      if (!is.null(settings$detect_mode)) updateRadioButtons(session, "detect_mode", selected = settings$detect_mode)
      if (!is.null(settings$detect_top_genes)) updateNumericInput(session, "detect_top_genes", value = as.numeric(settings$detect_top_genes))
      if (!is.null(settings$detect_min_freq)) updateNumericInput(session, "detect_min_freq", value = as.numeric(settings$detect_min_freq))
      if (!is.null(settings$detect_resolution)) updateNumericInput(session, "detect_resolution", value = as.numeric(settings$detect_resolution))
      if (!is.null(settings$detect_nboot)) updateNumericInput(session, "detect_nboot", value = as.numeric(settings$detect_nboot))

      output$import_settings_status <- renderUI({
        tags$div(class = "text-success mt-2", "✅ Settings imported successfully!")
      })
    }, error = function(e) {
      output$import_settings_status <- renderUI({
        tags$div(class = "text-danger mt-2", paste("❌ Error loading settings:", e$message))
      })
    })
  })

  # Export Settings
  output$export_settings <- downloadHandler(
    filename = function() {
      paste0("gsmr_settings_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".yaml")
    },
    content = function(file) {
      tryCatch({
        current_values <- list(
          expression_file = rv$expr_file_path %||% "",
          groups_file = rv$groups_file_path %||% "",
          grouping_method = input$grouping_method,
          groups_format = input$groups_format,
          bridge_prior_source = input$bridge_prior_source,
          bridge_data_type = input$bridge_data_type,
          bridge_lambda = input$bridge_lambda,
          bridge_beta = input$bridge_beta,
          bridge_lambda2 = input$bridge_lambda2,
          bridge_lambda1 = input$bridge_lambda1,
          bridge_resolution = input$bridge_resolution,
          bridge_tau = input$bridge_tau,
          bridge_prefilter = input$bridge_prefilter,
          bridge_precision = input$bridge_precision,
          bridge_max_genes = input$bridge_max_genes,
          bridge_preview = input$bridge_preview,
          min_group_size = input$min_group_size,
          max_group_size = input$max_group_size,
          collapse_method = input$collapse_method,
          meta_conflict = input$meta_conflict,
          keep_original = input$keep_original,
          prescore_method = input$prescore_method,
          top_percent = input$top_percent,
          train_ratio = input$train_ratio,
          prescore_seed = input$prescore_seed,
          perm_method = input$perm_method,
          perm_nperm = input$perm_nperm,
          score_method = input$score_method,
          rf_ntree = input$rf_ntree,
          svm_kernel = input$svm_kernel,
          svm_cost = input$svm_cost,
          svm_gamma = input$svm_gamma,
          knn_k = input$knn_k,
          shave_genes = input$shave_genes,
          shave_method = input$shave_method,
          shave_k = input$shave_k,
          fdr_alpha = input$fdr_alpha,
          normalize_z = input$normalize_z,
          score_scope = input$score_scope,
          n_top_groups = input$n_top_groups,
          score_seed = input$score_seed,
          run_mccv_mode = input$run_mccv_mode,
          mccv_iterations = input$mccv_iterations,
          model_method = input$model_method,
          mccv_train_ratio = input$mccv_train_ratio,
          mccv_cores = input$mccv_cores,
          detect_mode = input$detect_mode,
          detect_top_genes = input$detect_top_genes,
          detect_min_freq = input$detect_min_freq,
          detect_resolution = input$detect_resolution,
          detect_nboot = input$detect_nboot
        )
        gsmR::export_settings(current_values, file)
      }, error = function(e) {
        cat("EXPORT ERROR: ", e$message, "\n", file = stderr())
        writeLines(paste("Error generating YAML:", e$message), file)
      })
    }
  )

  # -- Browse button: open native OS file chooser (macOS osascript) --
  observeEvent(input$btn_browse_save, {
    init_dir <- dirname(input$save_session_path)
    init_file <- basename(input$save_session_path)
    chosen <- tryCatch({
      if (Sys.info()[["sysname"]] == "Darwin") {
        # macOS: use osascript for native Save dialog
        cmd <- sprintf(
          "osascript -e 'tell application \"System Events\" to set p to POSIX path of (choose file name with prompt \"Save session as:\" default name \"%s\" default location POSIX file \"%s\")'",
          init_file, init_dir)
        trimws(system(cmd, intern = TRUE))
      } else if (requireNamespace("tcltk", quietly = TRUE)) {
        as.character(tcltk::tkgetSaveFile(
          initialdir = init_dir, initialfile = init_file,
          filetypes = "{{RDS Files} {.rds}} {{All Files} {*}}"))
      } else ""
    }, error = function(e) "")
    if (nchar(chosen) > 0) {
      if (!grepl("\\.rds$", chosen, ignore.case = TRUE)) chosen <- paste0(chosen, ".rds")
      updateTextInput(session, "save_session_path", value = chosen)
    }
  })

  # -- Save button: write session to the chosen path --
  observeEvent(input$btn_save_session, {
    save_path <- input$save_session_path
    if (is.null(save_path) || nchar(trimws(save_path)) == 0) {
      output$save_session_status <- renderUI(
        tags$div(class = "text-danger mt-2", "\u274c Please specify a save path."))
      return()
    }
    tryCatch({
      dir.create(dirname(save_path), showWarnings = FALSE, recursive = TRUE)
      session_data <- list(
        expr_data = rv$expr_data, groups = rv$groups,
        dataset = rv$dataset, dataset_name = rv$dataset_name,
        loaded = rv$loaded,
        collapse_result = rv$collapse_result, collapsed = rv$collapsed,
        time_collapse = rv$time_collapse,
        prescore_groups = rv$prescore_groups, prescore_scores = rv$prescore_scores,
        prescore_nkeep = rv$prescore_nkeep, prescored = rv$prescored,
        train_idx = rv$train_idx, time_prescore = rv$time_prescore,
        score_results = rv$score_results, final_model = rv$final_model,
        accumulative = rv$accumulative, mccv_results = rv$mccv_results,
        scored = rv$scored, score_summary = rv$score_summary,
        time_score = rv$time_score,
        timings = rv$timings,
        detect = rv$detect, permtest = rv$permtest,
        mccv_group_summary = rv$mccv_group_summary,
        settings = list(
          score_method = input$score_method, run_mccv_mode = input$run_mccv_mode,
          mccv_iterations = input$mccv_iterations, model_method = input$model_method,
          mccv_train_ratio = input$mccv_train_ratio, score_scope = input$score_scope,
          n_top_groups = input$n_top_groups, score_seed = input$score_seed,
          top_percent = input$top_percent, normalize = input$normalize,
          detect_mode = input$detect_mode, detect_top_genes = input$detect_top_genes,
          detect_min_freq = input$detect_min_freq, detect_resolution = input$detect_resolution,
          detect_nboot = input$detect_nboot
        ),
        saved_at = Sys.time(), r_version = R.version.string,
        gsmR_version = tryCatch(as.character(utils::packageVersion("gsmR")),
                                error = function(e) "unknown")
      )
      saveRDS(session_data, save_path)
      output$save_session_status <- renderUI(
        tags$div(class = "text-success mt-2",
                 sprintf("\u2705 Session saved to %s", save_path)))
      showNotification(sprintf("\U0001F4BE Session saved to %s", save_path),
                       type = "message", duration = 5)
    }, error = function(e) {
      output$save_session_status <- renderUI(
        tags$div(class = "text-danger mt-2",
                 paste0("\u274c Save failed: ", conditionMessage(e))))
    })
  })

  observeEvent(input$load_session_file, {
    req(input$load_session_file)
    f <- input$load_session_file
    if (!grepl("\\.rds$", tolower(f$name))) {
      output$load_session_status <- renderUI(
        tags$div(class = "text-danger mt-2", "\u274c Wrong file type. Please upload a .rds session file.")
      )
      return()
    }
    tryCatch({
      sess <- readRDS(f$datapath)
      # Validate it's a gsmR session
      if (!is.list(sess) || is.null(sess$saved_at)) {
        output$load_session_status <- renderUI(
          tags$div(class = "text-danger mt-2",
                   "\u274c This file does not appear to be a valid gsmR session.")
        )
        return()
      }
      # Restore all reactive values
      rv$expr_data        <- sess$expr_data
      rv$groups           <- sess$groups
      rv$dataset          <- sess$dataset
      rv$dataset_name     <- sess$dataset_name
      rv$loaded           <- sess$loaded %||% FALSE
      rv$collapse_result  <- sess$collapse_result
      rv$collapsed        <- sess$collapsed %||% FALSE
      rv$time_collapse    <- sess$time_collapse
      rv$prescore_groups  <- sess$prescore_groups
      rv$prescore_scores  <- sess$prescore_scores
      rv$prescore_nkeep   <- sess$prescore_nkeep
      rv$prescored        <- sess$prescored %||% FALSE
      rv$train_idx        <- sess$train_idx
      rv$time_prescore    <- sess$time_prescore
      rv$score_results    <- sess$score_results
      rv$final_model      <- sess$final_model
      rv$accumulative     <- sess$accumulative
      rv$mccv_results     <- sess$mccv_results
      rv$scored           <- sess$scored %||% FALSE
      rv$score_summary    <- sess$score_summary
      rv$time_score       <- sess$time_score
      rv$timings          <- sess$timings %||% list()
      rv$detect           <- sess$detect
      rv$permtest         <- sess$permtest
      rv$mccv_group_summary <- sess$mccv_group_summary
      # Restore UI settings
      s <- sess$settings
      if (!is.null(s)) {
        if (!is.null(s$score_method))      updateRadioButtons(session, "score_method", selected = s$score_method)
        if (!is.null(s$run_mccv_mode))     updateCheckboxInput(session, "run_mccv_mode", value = s$run_mccv_mode)
        if (!is.null(s$mccv_iterations))   updateNumericInput(session, "mccv_iterations", value = s$mccv_iterations)
        if (!is.null(s$model_method))      updateSelectInput(session, "model_method", selected = s$model_method)
        if (!is.null(s$mccv_train_ratio))  updateSliderInput(session, "mccv_train_ratio", value = s$mccv_train_ratio)
        if (!is.null(s$score_scope))       updateRadioButtons(session, "score_scope", selected = s$score_scope)
        if (!is.null(s$n_top_groups))      updateNumericInput(session, "n_top_groups", value = s$n_top_groups)
        if (!is.null(s$score_seed))        updateNumericInput(session, "score_seed", value = s$score_seed)
        if (!is.null(s$normalize))         updateCheckboxInput(session, "normalize", value = s$normalize)
        if (!is.null(s$detect_mode))       updateRadioButtons(session, "detect_mode", selected = s$detect_mode)
        if (!is.null(s$detect_top_genes))  updateNumericInput(session, "detect_top_genes", value = s$detect_top_genes)
        if (!is.null(s$detect_min_freq))   updateNumericInput(session, "detect_min_freq", value = s$detect_min_freq)
        if (!is.null(s$detect_resolution)) updateNumericInput(session, "detect_resolution", value = s$detect_resolution)
        if (!is.null(s$detect_nboot))      updateNumericInput(session, "detect_nboot", value = s$detect_nboot)
      }
      # Update Discover module picker if detection results exist
      if (!is.null(sess$detect) && !is.null(sess$detect$modules)) {
        updateSelectInput(session, "detect_module_pick",
                          choices = sess$detect$modules$Module)
      }
      output$load_session_status <- renderUI(
        tags$div(class = "text-success mt-2",
                 sprintf("\u2705 Session restored! Saved %s with R %s.",
                         format(sess$saved_at, "%Y-%m-%d %H:%M"),
                         sess$r_version %||% "?"))
      )
      showNotification("\u2705 Session loaded \u2014 all results restored!", type = "message", duration = 6)
    }, error = function(e) {
      output$load_session_status <- renderUI(
        tags$div(class = "text-danger mt-2",
                 paste0("\u274c Failed to load session: ", conditionMessage(e)))
      )
    })
  })

  # \u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550
  # Multi-dataset batch processing
  #
  # Drives the shared scheduler (gsm_batch_new_state / gsm_batch_step) from a
  # non-blocking observe()/invalidateLater() loop \u2014 the same concurrency model as
  # the live MCCV view, but the work queue is every (dataset, fold) pair over one
  # shared core pool. Per-dataset .rds/tables/logs and a combined report are
  # written under a per-session temp directory.
  # \u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550
  bt <- reactiveValues(active = FALSE, st = NULL, snap = NULL, done = FALSE,
                       outdir = NULL, report = NULL, total = 0L,
                       report_format = "xlsx")

  output$batch_done <- reactive(isTRUE(bt$done))
  outputOptions(output, "batch_done", suspendWhenHidden = FALSE)

  observeEvent(input$batch_select_all, {
    updateSelectInput(session, "batch_datasets", selected = .example_inputs)
  })
  observeEvent(input$batch_clear_all, {
    updateSelectInput(session, "batch_datasets", selected = character(0))
  })

  observeEvent(input$btn_batch, {
    # Datasets = selected examples + uploaded files (mixed).
    sel <- input$batch_datasets
    example_files <- if (length(sel)) {
      stats::setNames(file.path(.example_input_dir, sel),
                      tools::file_path_sans_ext(basename(sel)))
    } else character(0)
    up <- input$batch_upload
    upload_files <- if (!is.null(up) && nrow(up)) {
      stats::setNames(up$datapath, tools::file_path_sans_ext(up$name))
    } else character(0)
    expr_files <- c(example_files, upload_files)
    if (length(expr_files) < 1) {
      showNotification("Select at least one example or upload one expression file.",
                       type = "warning"); return()
    }
    # Disambiguate duplicate dataset names (e.g. example + upload同名).
    names(expr_files) <- make.unique(names(expr_files), sep = "_")

    # Grouping: an upload overrides the example dropdown.
    gup <- input$batch_groups_upload
    if (!is.null(gup) && nrow(gup)) {
      groups_path <- gup$datapath
      groups_fmt  <- if (grepl("[.]gmt$", gup$name, ignore.case = TRUE)) "gmt" else "csv"
    } else {
      gsel <- input$batch_groups
      is_toy <- identical(gsel, "__toy__")
      groups_path <- if (is_toy)
        system.file("extdata", "test_groups.gmt", package = "gsmR")
      else file.path(.example_groups_dir, gsel)
      groups_fmt <- if (is_toy) "gmt" else "auto"
    }

    params <- list(
      score_method     = input$bp_score_method,
      prescore_method  = input$bp_prescore_method,
      mccv_iterations  = input$bp_mccv_iterations,
      top_percent      = input$bp_top_percent,
      mccv_train_ratio = input$bp_train_ratio,
      n_top_groups     = input$bp_n_top_groups,
      model_method     = input$bp_model_method,
      shave_genes      = input$bp_shave,
      normalize_z      = input$bp_normalize,
      score_seed       = input$bp_seed
    )
    outdir <- file.path(tempdir(),
                        paste0("gsmR_batch_", format(Sys.time(), "%Y%m%d_%H%M%S")))

    st <- tryCatch(
      gsmR:::gsm_batch_new_state(expr_files, groups_path, params = params, outdir = outdir,
                          cores = input$bp_cores, max_resident = input$bp_max_resident,
                          discover = isTRUE(input$bp_discover),
                          groups_format = groups_fmt),
      error = function(e) {
        showNotification(paste("Batch init failed:", conditionMessage(e)),
                         type = "error", duration = 10); NULL
      })
    if (is.null(st)) return()

    bt$st <- st; bt$active <- TRUE; bt$done <- FALSE; bt$outdir <- outdir
    bt$total <- length(st$names) * st$n_iter
    bt$report_format <- input$bp_report_format
    bt$snap <- gsmR:::gsm_batch_snapshot(st)
    showNotification(sprintf("\U0001F680 Batch started: %d datasets on %d cores",
                             length(st$names), st$cores), type = "message", duration = 5)
  })

  # Non-blocking scheduler loop
  observe({
    req(bt$active)
    st <- bt$st
    isolate({
      gsmR:::gsm_batch_step(st)
      bt$snap <- gsmR:::gsm_batch_snapshot(st)
      if (!isTRUE(st$active)) {
        fmt <- if (is.null(bt$report_format)) "xlsx" else bt$report_format
        tryCatch(gsmR:::gsm_write_batch_outputs(st, format = fmt),
                 error = function(e) NULL)
        cf <- file.path(st$outdir, "combined",
                        if (fmt == "xlsx") "combined_report.xlsx" else "combined_tables.zip")
        bt$report <- if (file.exists(cf)) cf else NULL
        bt$active <- FALSE; bt$done <- TRUE
        showNotification("\u2705 Batch complete!", type = "message")
      }
    })
    if (isTRUE(bt$active)) invalidateLater(400)
  })

  output$batch_overall <- renderUI({
    sn <- bt$snap
    if (is.null(sn)) return(NULL)
    done_folds <- sum(sn$done_count)
    pct <- if (bt$total > 0) round(100 * done_folds / bt$total) else 0
    n_done <- sum(sn$status == "done"); n_tot <- length(sn$status)
    n_err <- sum(sn$status == "error")
    tags$div(
      tags$h5(sprintf("%s \u2014 %d/%d datasets done%s | %d/%d folds (%d%%) | running %d, queued %d",
                      if (isTRUE(sn$all_done)) "\u2705 Complete" else "\u23f3 Running",
                      n_done, n_tot, if (n_err > 0) sprintf(", %d failed", n_err) else "",
                      done_folds, bt$total, pct, sn$running, sn$queued)),
      tags$div(class = "progress", style = "height:22px;",
        tags$div(class = "progress-bar bg-success", role = "progressbar",
                 style = sprintf("width:%d%%;", pct), sprintf("%d%%", pct)))
    )
  })

  output$batch_alert <- renderUI({
    sn <- bt$snap
    if (is.null(sn)) return(NULL)
    errs <- sn$errors
    n_err <- length(errs)
    if (n_err == 0) return(NULL)
    n_tot <- length(sn$status)
    all_failed <- n_err >= n_tot
    # Dominant failure reason → one shared suggestion.
    tab <- sort(table(unlist(unname(errs))), decreasing = TRUE)
    hint <- batch_fail_hint(names(tab)[1])
    cls <- if (all_failed) "alert-danger" else "alert-warning"
    icon <- if (all_failed) "\U274C" else "\U26A0\UFE0F"
    tags$div(class = paste("alert", cls), role = "alert",
      style = "border-radius:12px; margin-bottom:14px;",
      tags$div(style = "font-weight:700; font-size:1.02rem; margin-bottom:4px;",
        sprintf("%s %d of %d dataset%s failed", icon, n_err, n_tot,
                if (n_tot == 1) "" else "s")),
      tags$div(style = "margin-bottom:6px;", tags$strong("Likely fix: "), hint),
      tags$details(
        tags$summary(class = "small", style = "cursor:pointer;",
                     "Show per-dataset error messages"),
        tags$ul(class = "small mb-0", style = "margin-top:6px;",
          lapply(names(errs), function(nm)
            tags$li(tags$strong(nm), ": ", errs[[nm]])))))
  })

  output$batch_status <- renderUI({
    sn <- bt$snap
    if (is.null(sn)) return(tags$div(class = "gsm-empty",
      tags$span(class = "gsm-empty-icon", "\U0001F4CA"),
      tags$div(class = "gsm-empty-text",
               "No batch running yet. Pick datasets + a grouping file in the ",
               tags$strong("Settings"), " panel, then click ",
               tags$strong("\U0001F680 Run Batch"), ".")))
    errs <- sn$errors
    rows <- lapply(names(sn$status), function(nm) {
      stt <- sn$status[[nm]]; dc <- sn$done_count[[nm]]
      badge <- switch(stt, done = "success", running = "info",
                      error = "danger", loading = "warning", "secondary")
      frac <- if (sn$n_iter > 0) round(100 * dc / sn$n_iter) else 0
      err_msg <- if (!is.null(errs) && nm %in% names(errs)) errs[[nm]] else NULL
      tags$div(
        tags$div(class = "d-flex align-items-center gap-2 mb-1",
          tags$span(class = sprintf("badge bg-%s", badge),
                    style = "min-width:64px;", stt),
          tags$strong(style = "min-width:160px;", nm),
          tags$div(class = "progress flex-grow-1", style = "height:16px;",
            tags$div(class = "progress-bar", role = "progressbar",
                     style = sprintf("width:%d%%;", frac))),
          tags$span(class = "text-muted small", sprintf("%d/%d", dc, sn$n_iter))),
        if (!is.null(err_msg))
          tags$div(class = "small text-danger mb-2",
                   style = "margin-left:72px;",
                   tags$span(style = "opacity:0.85;", "\U26A0\UFE0F ", err_msg)))
    })
    do.call(tags$div, rows)
  })

  output$batch_log <- renderText({
    sn <- bt$snap
    if (is.null(sn)) return("")
    paste(sn$log_tail, collapse = "\n")
  })

  output$batch_outdir_note <- renderText({
    if (is.null(bt$outdir)) "" else
      paste0("Per-dataset .rds + tables + logs under: ",
             file.path(bt$outdir, "per_dataset"))
  })

  output$batch_run_info <- renderDT({
    req(bt$done, bt$st)
    info <- tryCatch(gsmR:::gsm_batch_run_info(bt$st), error = function(e) NULL)
    req(info)
    DT::datatable(info, rownames = FALSE,
                  options = list(scrollX = TRUE, dom = "t", pageLength = 100))
  })

  output$batch_combined_metrics <- renderDT({
    req(bt$done, bt$st)
    res <- bt$st$results
    if (length(res) == 0) return(NULL)
    tb <- tryCatch(gsm_combined_metrics(res), error = function(e) NULL)
    req(tb)
    num <- vapply(tb, is.numeric, logical(1))
    tb[num] <- lapply(tb[num], function(x) round(x, 3))
    DT::datatable(tb, rownames = FALSE,
                  options = list(scrollX = TRUE, dom = "t", pageLength = 50))
  })

  output$dl_batch_report <- downloadHandler(
    filename = function() {
      if (is.null(bt$report)) "combined_report.zip" else basename(bt$report)
    },
    content = function(file) {
      req(bt$report)
      file.copy(bt$report, file, overwrite = TRUE)
    }
  )

  # Zip the entire batch output tree (per-dataset .rds/tables/logs + combined + manifest).
  output$dl_batch_all <- downloadHandler(
    filename = function() paste0("gsmR_batch_results_",
                                 format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip"),
    content = function(file) {
      req(bt$outdir, dir.exists(bt$outdir))
      old <- setwd(dirname(bt$outdir)); on.exit(setwd(old), add = TRUE)
      utils::zip(zipfile = file, files = basename(bt$outdir), flags = "-r9Xq")
    }
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# Run
# ══════════════════════════════════════════════════════════════════════════════
shinyApp(ui = ui, server = server)
