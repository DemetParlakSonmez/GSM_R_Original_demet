#' @title BRIDGE grouping engine
#' @description BRIDGE (Biologically-Regularized Interaction-Differential Grouping
#'   Engine) is a data-driven grouping method for the G (Grouping) step. Unlike
#'   the file-based readers ([read_groups()], [read_gmt()]), which take a
#'   predefined grouping, BRIDGE *derives* groups from the expression matrix `X`,
#'   the class labels `y`, and an optional prior biological graph `B`.
#'
#'   It estimates two class-conditional sparse precision (partial-correlation)
#'   networks with a knowledge-weighted graphical lasso, forms a differential
#'   network `Delta = Theta1 - Theta0`, runs multi-resolution community
#'   detection on a graph that blends retained conditional edges with
#'   `|Delta|`, and selects the resulting modules by a Grouping-Scoring coupled
#'   objective `J(M)`. The output is the **same named list of character vectors**
#'   (`group name -> gene symbols`) that the file readers return, so the Scoring
#'   and Modeling steps consume it unchanged.
#' @name bridge
NULL

# ── Prior scaffold ───────────────────────────────────────────────────────────

#' Build a prior relatedness matrix `B` from a named grouping list
#'
#' Genes that co-occur in one or more groups get a prior edge. The raw
#' co-membership count is scaled to `[0, 1]` by dividing by its maximum so the
#' result is a valid `b_ij` weight matrix.
#'
#' @param groups Named list of character vectors (group name -> gene symbols),
#'   e.g. the output of [read_groups()].
#' @param features Character vector of feature (gene) IDs to build `B` over,
#'   typically `rownames(X)`. Genes not in `features` are ignored.
#' @return A symmetric `Matrix::dsCMatrix` (sparse) with `dimnames` set to
#'   `features`, values in `[0, 1]`, zero diagonal.
#' @export
bridge_prior_from_groups <- function(groups, features) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("BRIDGE needs the 'Matrix' package. Install it with ",
         "install.packages(\"Matrix\").", call. = FALSE)
  }
  stopifnot(is.character(features), length(features) > 0)
  p <- length(features)
  idx <- stats::setNames(seq_len(p), features)

  i_acc <- integer(0); j_acc <- integer(0)
  for (g in groups) {
    members <- unique(intersect(g, features))
    if (length(members) < 2) next
    mi <- idx[members]
    cmb <- utils::combn(sort(mi), 2)
    i_acc <- c(i_acc, cmb[1, ])
    j_acc <- c(j_acc, cmb[2, ])
  }
  if (length(i_acc) == 0) {
    B <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                              dims = c(p, p), symmetric = TRUE)
    dimnames(B) <- list(features, features)
    return(B)
  }
  # Count co-memberships per unordered pair.
  pair <- paste(i_acc, j_acc, sep = "_")
  tab <- table(pair)
  keys <- strsplit(names(tab), "_", fixed = TRUE)
  ii <- vapply(keys, function(k) as.integer(k[1]), integer(1))
  jj <- vapply(keys, function(k) as.integer(k[2]), integer(1))
  w <- as.numeric(tab)
  w <- w / max(w)  # scale to [0, 1]
  B <- Matrix::sparseMatrix(i = ii, j = jj, x = w, dims = c(p, p),
                            symmetric = TRUE)
  dimnames(B) <- list(features, features)
  B
}

#' Build a prior relatedness matrix `B` from an edge list or matrix
#'
#' Accepts either a 3-column edge list (`gene_i`, `gene_j`, `weight`) or an
#' already-formed symmetric weight matrix, and aligns it to `features`.
#'
#' @param x Either a `data.frame`/`matrix` edge list whose first two columns are
#'   gene IDs and whose third column is a numeric weight in `[0, 1]`, or a
#'   (possibly sparse) square symmetric matrix with gene IDs as `dimnames`.
#' @param features Character vector of feature IDs to align `B` to
#'   (typically `rownames(X)`).
#' @return A symmetric sparse matrix over `features`, values in `[0, 1]`.
#' @export
bridge_prior_from_edgelist <- function(x, features) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("BRIDGE needs the 'Matrix' package. Install it with ",
         "install.packages(\"Matrix\").", call. = FALSE)
  }
  stopifnot(is.character(features), length(features) > 0)
  p <- length(features)
  idx <- stats::setNames(seq_len(p), features)

  # Square matrix input -> subset/align to features.
  if ((is.matrix(x) || methods::is(x, "Matrix")) && nrow(x) == ncol(x) &&
      !is.null(rownames(x))) {
    common <- intersect(rownames(x), features)
    if (length(common) == 0) {
      stop("Prior matrix feature names do not overlap with X. ",
           "First few X features: ",
           paste(utils::head(features, 5), collapse = ", "),
           " | first few prior features: ",
           paste(utils::head(rownames(x), 5), collapse = ", "), call. = FALSE)
    }
    sub <- as.matrix(x[common, common, drop = FALSE])
    B <- Matrix::Matrix(0, p, p, sparse = TRUE)
    dimnames(B) <- list(features, features)
    B[common, common] <- sub
    B <- (B + Matrix::t(B)) / 2  # enforce symmetry
    Matrix::diag(B) <- 0
    B <- pmin(pmax(B, 0), 1)
    return(methods::as(B, "symmetricMatrix"))
  }

  # Edge-list input.
  df <- as.data.frame(x, stringsAsFactors = FALSE)
  if (ncol(df) < 3) {
    stop("Prior edge list must have at least 3 columns: gene_i, gene_j, weight.",
         call. = FALSE)
  }
  gi <- as.character(df[[1]]); gj <- as.character(df[[2]])
  w  <- as.numeric(df[[3]])
  keep <- gi %in% features & gj %in% features & is.finite(w)
  if (!any(keep)) {
    stop("No prior edges reference features present in X. Check that the edge ",
         "list gene IDs match rownames(X).", call. = FALSE)
  }
  gi <- gi[keep]; gj <- gj[keep]; w <- pmin(pmax(w[keep], 0), 1)
  B <- Matrix::sparseMatrix(i = idx[gi], j = idx[gj], x = w, dims = c(p, p),
                            symmetric = FALSE)
  B <- B + Matrix::t(B)                 # symmetrize
  Matrix::diag(B) <- 0
  B@x <- pmin(B@x, 1)
  dimnames(B) <- list(features, features)
  methods::as(B, "symmetricMatrix")
}

# ── Label-blind pre-filter (speed; no leakage) ───────────────────────────────

#' Reduce the feature set before BRIDGE estimation, without looking at labels
#'
#' The class-conditional graphical lasso scales ~cubically in the number of
#' features, so on transcriptome-scale data BRIDGE must run on a reduced gene
#' set. This filter is deliberately **label-blind** — it never inspects `y` — so
#' it can be applied once, globally, before cross-validation without leaking any
#' class information into the folds. Two complementary strategies:
#'
#' * `"variance"` keeps the `max_genes` most variable features (by variance).
#' * `"prior"` keeps only features that appear in the prior graph `B` (genes
#'   known to be biologically related), which is the most principled restriction
#'   for a knowledge-regularized method.
#'
#' `"prior"` and `"variance"` compose: with a prior supplied and `max_genes`
#' set, the prior restriction is applied first, then the variance cap.
#'
#' @param X Numeric expression matrix (genes x samples); row names required.
#' @param method One of `"none"`, `"variance"`, `"prior"`. Default `"none"`.
#' @param max_genes Optional integer cap on the number of features to keep
#'   (applied by the variance strategy). `NULL` = no cap.
#' @param B Optional prior over the features (named grouping list, edge list, or
#'   symmetric matrix) used by the `"prior"` strategy.
#' @return A character vector of the retained feature (row) names.
#' @export
bridge_prefilter <- function(X, method = c("none", "variance", "prior"),
                             max_genes = NULL, B = NULL) {
  method <- match.arg(method)
  stopifnot(is.matrix(X), !is.null(rownames(X)))
  feats <- rownames(X)

  if (method == "prior") {
    if (is.null(B)) {
      warning("bridge_prefilter(method = \"prior\") needs a prior B; ",
              "keeping all features.", call. = FALSE)
    } else {
      Bmat <- bridge_resolve_prior(B, feats, verbose = FALSE)
      if (!is.null(Bmat)) {
        deg <- Matrix::rowSums(Bmat != 0)
        keep <- feats[deg > 0]
        if (length(keep) >= 2) feats <- keep
      }
    }
  }

  if (method == "variance" || (!is.null(max_genes) && length(feats) > max_genes)) {
    v <- apply(X[feats, , drop = FALSE], 1, stats::var)
    ord <- feats[order(v, decreasing = TRUE)]
    cap <- if (is.null(max_genes)) length(ord) else min(max_genes, length(ord))
    feats <- ord[seq_len(cap)]
  }

  feats
}

# ── Class-conditional dependency estimation ──────────────────────────────────

#' Estimate two class-conditional sparse precision matrices
#'
#' Splits samples by class and estimates `Theta0`, `Theta1` with a
#' knowledge-weighted graphical lasso. The per-edge penalty is
#' `lambda_ij = lambda * (1 - beta * b_ij)` so biologically related edges are
#' penalized less. When JGL is available a fused penalty
#' `lambda2 * sum|theta0_ij - theta1_ij|` is applied jointly so the two networks
#' are assumed equal unless the data disagree.
#'
#' @param X0,X1 Sample-by-feature numeric matrices for class 0 and class 1
#'   (rows = samples, columns = features). Note this is the transpose of the
#'   `genes x samples` convention used elsewhere; [bridge_grouping()] handles the
#'   transpose for you.
#' @param B Symmetric prior matrix over the features (columns of `X0`/`X1`),
#'   values in `[0, 1]`. May be `NULL` for `B = 0`.
#' @param lambda Base sparsity penalty.
#' @param beta Knowledge-weighting strength in `[0, 1)`.
#' @param lambda2 Fused penalty coupling the two networks (JGL backend only).
#' @param data_type One of `"gaussian"`, `"nonparanormal"`, `"compositional"`.
#' @param precision_method Estimation backend preference: `"glasso"` (default;
#'   full graphical lasso, `O(p^3)`, uses the JGL fused backend when available)
#'   or `"mb"` (Meinshausen-Buhlmann neighborhood selection via `huge`,
#'   `~O(p^2)`, much faster for large `p`).
#' @return A list with `Theta0`, `Theta1` (dense numeric precision matrices) and
#'   the `backend` used.
#' @keywords internal
bridge_estimate_precision <- function(X0, X1, B, lambda = 0.1, beta = 0.5,
                                      lambda2 = 0.05,
                                      data_type = "gaussian",
                                      precision_method = "glasso") {
  data_type <- match.arg(data_type,
                         c("gaussian", "nonparanormal", "compositional"))
  precision_method <- match.arg(precision_method, c("glasso", "mb"))
  if (data_type == "compositional") {
    stop("data_type = \"compositional\" is not yet supported. It needs a CLR + ",
         "SPIEC-EASI estimator (the 'SpiecEasi' package, currently GitHub-only: ",
         "remotes::install_github(\"zdk123/SpiecEasi\")). Use \"gaussian\" or ",
         "\"nonparanormal\" for now.", call. = FALSE)
  }

  p <- ncol(X0)
  Bm <- if (is.null(B)) matrix(0, p, p) else as.matrix(B)

  if (data_type == "nonparanormal") {
    if (!requireNamespace("huge", quietly = TRUE)) {
      stop("data_type = \"nonparanormal\" needs the 'huge' package. Install it ",
           "with install.packages(\"huge\").", call. = FALSE)
    }
    X0 <- huge::huge.npn(X0, verbose = FALSE)
    X1 <- huge::huge.npn(X1, verbose = FALSE)
  }

  # Per-edge penalty weights: biologically related edges penalized less.
  Wpen <- lambda * (1 - beta * Bm)
  Wpen[Wpen < 0] <- 0

  # Fast backend: Meinshausen-Buhlmann neighborhood selection (O(p^2)), per
  # class, for large p. Knowledge weighting enters via the base penalty lambda;
  # the fused term is not available here (independent per-class estimation).
  if (precision_method == "mb") {
    if (!requireNamespace("huge", quietly = TRUE)) {
      stop("precision_method = \"mb\" needs the 'huge' package. Install it ",
           "with install.packages(\"huge\").", call. = FALSE)
    }
    est_mb <- function(Xc) {
      fit <- huge::huge(as.matrix(Xc), lambda = lambda, method = "mb",
                        verbose = FALSE)
      # huge returns a path; with a single lambda take the first (only) icov/path.
      P <- fit$path[[1]]
      Theta <- as.matrix(P) * 1.0
      Theta <- (Theta + t(Theta)) / 2          # symmetrize the adjacency
      diag(Theta) <- 1
      dimnames(Theta) <- list(colnames(Xc), colnames(Xc))
      Theta
    }
    return(list(Theta0 = est_mb(X0), Theta1 = est_mb(X1), backend = "huge-mb"))
  }

  # Preferred backend: JGL (joint, with fused penalty). Requires equal columns.
  if (requireNamespace("JGL", quietly = TRUE)) {
    fit <- tryCatch(
      JGL::JGL(list(as.matrix(X0), as.matrix(X1)),
               penalty = "fused", lambda1 = lambda, lambda2 = lambda2,
               weights = "equal", return.whole.theta = TRUE),
      error = function(e) NULL
    )
    if (!is.null(fit)) {
      return(list(Theta0 = fit$theta[[1]], Theta1 = fit$theta[[2]],
                  backend = "JGL"))
    }
    # fall through to glasso on failure
  }

  # Fallback: independent knowledge-weighted glasso per class (no fused term).
  if (!requireNamespace("glasso", quietly = TRUE)) {
    stop("BRIDGE needs a graphical-lasso backend. Install 'JGL' (preferred, ",
         "adds the fused penalty) or 'glasso' with ",
         "install.packages(c(\"JGL\", \"glasso\")).", call. = FALSE)
  }
  est_one <- function(Xc) {
    S <- stats::cov(Xc)
    fit <- glasso::glasso(S, rho = Wpen, penalize.diagonal = FALSE)
    Theta <- fit$wi
    dimnames(Theta) <- list(colnames(Xc), colnames(Xc))
    Theta
  }
  list(Theta0 = est_one(X0), Theta1 = est_one(X1), backend = "glasso")
}

# ── Differential network + module extraction ─────────────────────────────────

#' Build the blended differential graph and extract multi-resolution modules
#'
#' Forms `Delta = Theta1 - Theta0`, blends `|Delta|` with the retained
#' conditional edges (union of the two networks' support) and the prior `B`, then
#' runs Leiden community detection swept over a resolution grid. Falls back to
#' hierarchical clustering when igraph is unavailable (mirrors
#' [detect_communities()]).
#'
#' @param Theta0,Theta1 Precision matrices from [bridge_estimate_precision()].
#' @param B Symmetric prior matrix (or `NULL`).
#' @param resolution Numeric vector of Leiden resolutions to sweep.
#' @return A named list of character vectors: candidate modules across all
#'   resolutions (names encode the resolution and cluster id).
#' @keywords internal
bridge_extract_modules <- function(Theta0, Theta1, B, resolution) {
  feats <- rownames(Theta0)
  p <- length(feats)
  Delta <- abs(Theta1 - Theta0)

  retained <- (Theta0 != 0) | (Theta1 != 0)  # union support
  Bm <- if (is.null(B)) matrix(0, p, p) else as.matrix(B)

  # Blend: differential magnitude on retained edges, nudged by the prior.
  W <- (Delta * retained) + 0.25 * Bm
  W <- (W + t(W)) / 2  # enforce exact symmetry (glasso output can drift)
  diag(W) <- 0
  W[!is.finite(W)] <- 0
  # Normalize to [0, 1] for stable community detection.
  if (max(W) > 0) W <- W / max(W)
  dimnames(W) <- list(feats, feats)

  modules <- list()
  for (res in resolution) {
    memb <- detect_communities(W, method = "leiden", resolution = res)
    for (cl in unique(memb)) {
      genes <- names(memb)[memb == cl]
      if (length(genes) >= 2) {
        key <- sprintf("BRIDGE_r%s_c%s", format(res, trim = TRUE), cl)
        modules[[key]] <- sort(genes)
      }
    }
  }
  modules
}

# ── Grouping-Scoring coupling ────────────────────────────────────────────────

#' Module eigengene discrimination D(M)
#'
#' First principal component of the module's expression, evaluated by simple
#' k-fold cross-validated discrimination (mean absolute standardized difference
#' of the eigengene between classes across held-out folds).
#' @keywords internal
bridge_discrimination <- function(genes, X_gs, y, k = 3) {
  sub <- X_gs[genes, , drop = FALSE]
  if (nrow(sub) < 1 || ncol(sub) < 4) return(0)
  # Eigengene = first PC of the (samples x genes) module matrix.
  m <- t(sub)
  m <- scale(m, center = TRUE, scale = TRUE)
  m[!is.finite(m)] <- 0
  pc <- tryCatch(stats::prcomp(m, center = FALSE, scale. = FALSE)$x[, 1],
                 error = function(e) rowMeans(m))
  n <- length(pc)
  folds <- ((sample(n) - 1) %% k) + 1
  vals <- vapply(seq_len(k), function(f) {
    test <- folds == f
    if (length(unique(y[test])) < 2) return(NA_real_)
    lv <- levels(y)
    a <- pc[test & y == lv[1]]; b <- pc[test & y == lv[2]]
    if (length(a) < 1 || length(b) < 1) return(NA_real_)
    s <- stats::sd(pc[test]); if (!is.finite(s) || s == 0) s <- 1
    abs(mean(a) - mean(b)) / s
  }, numeric(1))
  mean(vals, na.rm = TRUE)
}

#' Mean prior cohesion C(M)
#' @keywords internal
bridge_cohesion <- function(genes, B) {
  if (is.null(B) || length(genes) < 2) return(0)
  sub <- as.matrix(B[genes, genes, drop = FALSE])
  m <- sub[upper.tri(sub)]
  if (length(m) == 0) return(0)
  mean(m)
}

#' Score, stability-filter and select candidate modules
#'
#' Computes `J(M) = D(M) + lambda1 * C(M) - lambda2 * |M|` for each module,
#' filters by a bootstrap-stability threshold `R(M) >= tau`, and keeps modules
#' with positive objective. Deduplicates identical gene sets, preferring the
#' higher-scoring one.
#'
#' @param modules Named list of candidate modules.
#' @param X_gs Genes-by-samples expression matrix.
#' @param y Class-label factor.
#' @param B Prior matrix (or `NULL`).
#' @param lambda1 Weight on prior cohesion in `J(M)`.
#' @param lambda2 Size penalty in `J(M)`.
#' @param tau Bootstrap-stability threshold in `[0, 1]`.
#' @param n_boot Number of bootstrap resamples for stability.
#' @return Named list of the selected modules (character vectors), with a
#'   `"scores"` attribute (the `J(M)` values).
#' @keywords internal
bridge_score_modules <- function(modules, X_gs, y, B, lambda1 = 0.25,
                                 lambda2 = 0.05, tau = 0.6, n_boot = 20) {
  if (length(modules) == 0) return(structure(list(), scores = numeric(0)))
  n <- ncol(X_gs)

  score_one <- function(genes) {
    D <- bridge_discrimination(genes, X_gs, y)
    C <- bridge_cohesion(genes, B)
    D + lambda1 * C - lambda2 * length(genes) / 100
  }

  J <- vapply(modules, score_one, numeric(1))

  # Bootstrap stability: fraction of resamples where D(M) stays positive.
  stability <- vapply(modules, function(genes) {
    ok <- vapply(seq_len(n_boot), function(b) {
      idx <- sample(n, n, replace = TRUE)
      if (length(unique(y[idx])) < 2) return(NA)
      bridge_discrimination(genes, X_gs[, idx, drop = FALSE], y[idx]) > 0
    }, logical(1))
    mean(ok, na.rm = TRUE)
  }, numeric(1))

  keep <- is.finite(J) & J > 0 & stability >= tau
  sel <- modules[keep]
  Jsel <- J[keep]
  if (length(sel) == 0) return(structure(list(), scores = numeric(0)))

  # Deduplicate identical gene sets, keeping the higher J.
  ord <- order(Jsel, decreasing = TRUE)
  sel <- sel[ord]; Jsel <- Jsel[ord]
  sigs <- vapply(sel, function(g) paste(sort(g), collapse = "\r"), character(1))
  dup <- duplicated(sigs)
  sel <- sel[!dup]; Jsel <- Jsel[!dup]

  structure(sel, scores = Jsel)
}

# ── Public entry point ───────────────────────────────────────────────────────

#' Derive feature groups with BRIDGE
#'
#' Runs the full BRIDGE pipeline: knowledge-weighted class-conditional precision
#' estimation, differential-network module extraction across resolutions, and
#' Grouping-Scoring coupled module selection. Returns groups in the exact format
#' the Scoring step consumes: a **named list of character vectors**
#' (`group name -> gene symbols`).
#'
#' @param X Numeric expression/abundance matrix with **genes as rows and samples
#'   as columns** (same orientation as `gsm_dataset()$expr`). Row names are
#'   required.
#' @param y Class-label vector or factor, length `ncol(X)`. BRIDGE currently
#'   requires exactly two classes.
#' @param B Prior biological graph over the features. Either `NULL` (pure
#'   data-driven, `B = 0`, with a warning), a symmetric weight matrix / sparse
#'   `Matrix`, a 3-column edge list (`gene_i`, `gene_j`, `weight`), or a named
#'   grouping list (as from [read_groups()]) which is converted with
#'   [bridge_prior_from_groups()]. Feature names must align with `rownames(X)`.
#' @param lambda Base graphical-lasso sparsity penalty. Default `0.1`.
#' @param beta Knowledge-weighting strength in `[0, 1)`; higher = prior edges
#'   penalized less. Default `0.5`.
#' @param lambda2 Fused penalty coupling the two class networks (JGL backend).
#'   Default `0.05`.
#' @param lambda1 Weight on prior cohesion `C(M)` in the selection objective
#'   `J(M)`. Default `0.25`.
#' @param resolution Numeric vector of Leiden resolutions to sweep for
#'   multi-scale modules. Default `c(0.5, 1, 1.5, 2)`.
#' @param tau Bootstrap-stability threshold `R(M) >= tau` in `[0, 1]`.
#'   Default `0.6`.
#' @param data_type One of `"gaussian"` (default), `"nonparanormal"`
#'   (applies [huge::huge.npn()] first), or `"compositional"` (not yet
#'   supported; errors with guidance to install SpiecEasi).
#' @param precision_method Estimation backend: `"glasso"` (default) or `"mb"`
#'   (faster neighborhood selection for large `p`). See
#'   [bridge_estimate_precision()].
#' @param prefilter Label-blind feature pre-filter applied first: `"none"`
#'   (default), `"variance"`, or `"prior"`. See [bridge_prefilter()]. This never
#'   inspects `y`, so it is safe to apply outside cross-validation.
#' @param max_genes Optional cap on the number of features kept by the variance
#'   pre-filter. `NULL` (default) = no cap.
#' @param harden If `TRUE` (default) return crisp groups (a plain named list of
#'   gene vectors) so the existing Scoring step works unchanged. Soft membership
#'   is not produced in this release; the argument is kept so that a future
#'   personalized-PageRank soft-membership mode can set `harden = FALSE` and
#'   return fuzzy assignments without changing this signature.
#' @param min_size,max_size Keep only modules whose size is in
#'   `[min_size, max_size]`. Defaults `5` and `500` (matching `gsm_dataset()`).
#' @param seed Random seed for reproducibility. Default `1`.
#' @param verbose Emit progress messages. Default `TRUE`.
#' @return A named list of character vectors (`group name -> gene symbols`),
#'   identical in shape to [read_groups()] output. Carries a `"scores"`
#'   attribute with each module's `J(M)`. If no module passes the filters, a
#'   warning is issued and an empty named list is returned.
#' @examples
#' \donttest{
#' set.seed(1)
#' # 40 genes x 60 samples, two classes; genes 1-8 are a differential module.
#' p <- 40; n <- 60
#' X <- matrix(stats::rnorm(p * n), nrow = p,
#'             dimnames = list(paste0("G", seq_len(p)), paste0("S", seq_len(n))))
#' y <- factor(rep(c("A", "B"), each = n / 2))
#' mod <- paste0("G", 1:8)
#' X[mod, y == "B"] <- X[mod, y == "B"] + 2  # shift the module in class B
#' # Prior: the true module is known to be related.
#' groups_prior <- list(known_module = mod)
#' if (requireNamespace("glasso", quietly = TRUE) ||
#'     requireNamespace("JGL", quietly = TRUE)) {
#'   g <- bridge_grouping(X, y, B = groups_prior, resolution = c(1, 1.5))
#'   str(g)
#' }
#' }
#' @export
bridge_grouping <- function(X, y, B = NULL,
                            lambda = 0.1, beta = 0.5, lambda2 = 0.05,
                            lambda1 = 0.25,
                            resolution = c(0.5, 1, 1.5, 2),
                            tau = 0.6,
                            data_type = c("gaussian", "nonparanormal",
                                          "compositional"),
                            precision_method = c("glasso", "mb"),
                            prefilter = c("none", "variance", "prior"),
                            max_genes = NULL,
                            harden = TRUE, min_size = 5, max_size = 500,
                            seed = 1, verbose = TRUE) {
  data_type <- match.arg(data_type)
  precision_method <- match.arg(precision_method)
  prefilter <- match.arg(prefilter)
  if (!is.null(seed)) set.seed(seed)

  # ── Validate X / y ─────────────────────────────────────────────────────────
  if (!is.matrix(X) || !is.numeric(X)) {
    stop("`X` must be a numeric matrix (genes x samples).", call. = FALSE)
  }
  if (is.null(rownames(X))) {
    stop("`X` must have row names (feature/gene IDs).", call. = FALSE)
  }
  y <- if (is.factor(y)) y else factor(y)
  if (length(y) != ncol(X)) {
    stop("length(y) (", length(y), ") must equal ncol(X) (", ncol(X), ").",
         call. = FALSE)
  }
  if (nlevels(y) != 2) {
    stop("BRIDGE currently supports exactly two classes; found ", nlevels(y),
         ".", call. = FALSE)
  }
  # ── Label-blind pre-filter (speed; never inspects y) ───────────────────────
  if (prefilter != "none" || !is.null(max_genes)) {
    kept <- bridge_prefilter(X, method = prefilter, max_genes = max_genes, B = B)
    if (verbose && length(kept) < nrow(X)) {
      message("BRIDGE: pre-filter kept ", length(kept), "/", nrow(X),
              " features (label-blind).")
    }
    X <- X[kept, , drop = FALSE]
  }
  feats <- rownames(X)

  # ── Resolve the prior B into a matrix aligned to feats ─────────────────────
  Bmat <- bridge_resolve_prior(B, feats, verbose = verbose)

  # ── p >> n: restrict to connected components of B when it is informative ───
  active <- bridge_active_features(X, Bmat, verbose = verbose)
  Xa <- X[active, , drop = FALSE]
  Ba <- if (is.null(Bmat)) NULL else Bmat[active, active, drop = FALSE]

  # ── Split by class; estimate precision matrices ────────────────────────────
  lv <- levels(y)
  X0 <- t(Xa[, y == lv[1], drop = FALSE])   # samples x features
  X1 <- t(Xa[, y == lv[2], drop = FALSE])
  if (verbose) message("BRIDGE: estimating class-conditional networks over ",
                       length(active), " features...")
  est <- bridge_estimate_precision(X0, X1, Ba, lambda = lambda, beta = beta,
                                   lambda2 = lambda2, data_type = data_type,
                                   precision_method = precision_method)
  # Ensure dimnames are the active features.
  dimnames(est$Theta0) <- list(active, active)
  dimnames(est$Theta1) <- list(active, active)
  if (verbose) message("BRIDGE: precision backend = ", est$backend)

  # ── Differential network + multi-resolution modules ────────────────────────
  modules <- bridge_extract_modules(est$Theta0, est$Theta1, Ba, resolution)
  if (verbose) message("BRIDGE: ", length(modules),
                       " candidate modules across resolutions.")

  # ── Score / stability-filter / select ──────────────────────────────────────
  selected <- bridge_score_modules(modules, Xa, y, Ba, lambda1 = lambda1,
                                   lambda2 = lambda2, tau = tau)

  # Size filter (matches gsm_dataset conventions).
  if (length(selected) > 0) {
    sizes <- vapply(selected, length, integer(1))
    keep <- sizes >= min_size & sizes <= max_size
    sc <- attr(selected, "scores")[keep]
    selected <- selected[keep]
    attr(selected, "scores") <- sc
  }

  if (length(selected) == 0) {
    warning("BRIDGE selected no modules passing the J(M)/stability/size ",
            "filters. Try lowering `tau`, `min_size`, or `lambda`.",
            call. = FALSE)
    out <- structure(list(), scores = numeric(0))
    return(out)
  }

  if (verbose) message("BRIDGE: selected ", length(selected),
                       " module(s) as final groups.")

  # `harden = FALSE` is reserved for a future soft-membership mode; today the
  # engine always yields crisp modules, so hardening is a no-op either way.
  selected
}

# ── Internal glue ────────────────────────────────────────────────────────────

#' Resolve a user-supplied prior into a feature-aligned matrix (or NULL)
#' @keywords internal
bridge_resolve_prior <- function(B, feats, verbose = TRUE) {
  if (is.null(B)) {
    warning("No prior graph B supplied; BRIDGE runs data-driven (B = 0).",
            call. = FALSE)
    return(NULL)
  }
  # Named grouping list -> co-membership prior.
  if (is.list(B) && !is.data.frame(B) && !methods::is(B, "Matrix")) {
    return(bridge_prior_from_groups(B, feats))
  }
  # Square symmetric matrix already over features.
  if ((is.matrix(B) || methods::is(B, "Matrix")) && nrow(B) == ncol(B) &&
      !is.null(rownames(B))) {
    if (length(intersect(rownames(B), feats)) == 0) {
      stop("Prior matrix and X share no feature names. ",
           "Check that B's dimnames match rownames(X).", call. = FALSE)
    }
    return(bridge_prior_from_edgelist(B, feats))
  }
  # Edge list.
  bridge_prior_from_edgelist(B, feats)
}

#' Restrict to features in non-trivial connected components of B (p >> n aid)
#' @keywords internal
bridge_active_features <- function(X, B, verbose = TRUE) {
  feats <- rownames(X)
  p <- length(feats); n <- ncol(X)
  # Only prune when p >> n and B is informative enough to define components.
  if (is.null(B) || p <= 2 * n) return(feats)
  if (!requireNamespace("igraph", quietly = TRUE)) return(feats)

  g <- igraph::graph_from_adjacency_matrix(as.matrix(B > 0), mode = "undirected",
                                           diag = FALSE)
  comp <- igraph::components(g)
  big <- which(comp$csize >= 2)
  active <- feats[comp$membership %in% big]
  if (length(active) < 2) return(feats)  # nothing useful to restrict to
  if (verbose && length(active) < p) {
    message("BRIDGE: p >> n; restricting estimation to ", length(active),
            "/", p, " features in connected components of B.")
  }
  active
}
