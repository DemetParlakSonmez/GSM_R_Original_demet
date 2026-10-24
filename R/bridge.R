#' @title BRIDGE grouping engine
#' @description BRIDGE (Biologically-Regularized Interaction-Differential Grouping
#'   Engine) is a data-driven grouping method for the G (Grouping) step. Unlike
#'   the file-based readers ([read_groups()], [read_gmt()]), which take a
#'   predefined grouping, BRIDGE *derives* groups from the expression matrix `X`,
#'   the class labels `y`, and an optional prior biological graph `B`.
#'
#'   It estimates two class-conditional sparse precision (partial-correlation)
#'   networks and forms a differential network `Delta = Theta1 - Theta0`.
#'
#'   **How the prior `B` is used — two distinct stages:**
#'
#'   *Stage A – Network estimation (glasso backend only):*
#'   When the `glasso` backend is active (JGL not installed or JGL fails), `B`
#'   directly modifies the per-edge penalty as
#'   `lambda_ij = lambda * (1 - beta * b_ij)`, so biologically related edges are
#'   penalized less during sparse precision estimation. This is the
#'   "knowledge-weighted graphical lasso" step.
#'
#'   *Stage B – Network integration / module selection (always):*
#'   Regardless of the estimation backend, `B` enters the blended differential
#'   graph `W = |Delta| * retained_edges + 0.25 * B` and the cohesion term
#'   `C(M)` in the module-selection objective `J(M) = D(M) + lambda1 * C(M)`.
#'   These two mechanisms keep the knowledge prior active even when JGL is the
#'   estimation backend.
#'
#'   *JGL backend caveat:* JGL's fused graphical lasso API does not accept
#'   per-edge penalty matrices; it applies a scalar `lambda1` uniformly to both
#'   networks and a scalar `lambda2` fused coupling. Therefore, when JGL is the
#'   active backend, the `B`-weighted penalty (Stage A) is **not** applied during
#'   estimation. The prior still influences the result through Stage B.
#'   If you require knowledge-weighted estimation, use `precision_method = "glasso"`.
#'
#'   After network estimation, multi-resolution community detection is run on the
#'   blended differential graph, and modules are selected by the
#'   Grouping-Scoring coupled objective `J(M)`. The output is the **same named
#'   list of character vectors** (`group name -> gene symbols`) that the file
#'   readers return, so the Scoring and Modeling steps consume it unchanged.
#' @name bridge
NULL

# ── Prior scaffold ───────────────────────────────────────────────────────────

#' Build a prior relatedness matrix `B` from a named grouping list
#'
#' Genes that co-occur in one or more groups get a prior edge. The raw
#' co-membership count is scaled to `[0, 1]` by dividing by its maximum so the
#' result is a valid `b_ij` weight matrix.
#'
#' **Important Conceptual Note:**
#' The resulting values in `B` represent a **normalized prior co-membership score**,
#' reflecting purely how often two genes appear together across the given lists.
#' They must **not** be interpreted as a "biological probability", "statistical
#' significance", or "confidence probability" of true interaction. If your source
#' knowledge base provides true quantitative biological confidence scores
#' (e.g., STRING PPI interaction probabilities), do not use this grouping fallback.
#' Instead, supply those exact scores as a 3-column edge list directly to
#' [bridge_prior_from_edgelist()].
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
      "install.packages(\"Matrix\").",
      call. = FALSE
    )
  }
  stopifnot(is.character(features), length(features) > 0)
  p <- length(features)
  idx <- stats::setNames(seq_len(p), features)

  # HIZ VE ÖLÇEKLENEBİLİRLİK OPTİMİZASYONU:
  # Eskiden kullanılan combn() ve c() ile vektör büyütme (pair explosion)
  # devasa pathway'lerde (örn. 5000 gen = 12.5M çift) saatlerce sürebilir
  # ve R'ın hafızasını tüketirdi. Bunun yerine bir O(n) "incident matrix" (A)
  # kurulup A %*% t(A) ile co-membership hesaplanır (saf C düzeyinde sparse çarpım).

  # Her grup için, `features` içinde var olan üyelerin indexlerini al.
  valid_members <- lapply(groups, function(g) {
    v <- idx[g]
    as.integer(v[!is.na(v)])
  })

  lens <- vapply(valid_members, length, integer(1))
  keep <- lens >= 2

  if (!any(keep)) {
    B <- Matrix::sparseMatrix(
      i = integer(0), j = integer(0), x = numeric(0),
      dims = c(p, p), symmetric = TRUE
    )
    dimnames(B) <- list(features, features)
    return(B)
  }

  valid_members <- valid_members[keep]
  lens <- lens[keep]
  K <- length(valid_members)

  # A: Features x Groups (p x K) incident matrisi
  i_idx <- unlist(valid_members, use.names = FALSE)
  j_idx <- rep(seq_len(K), times = lens)
  A <- Matrix::sparseMatrix(i = i_idx, j = j_idx, x = 1.0, dims = c(p, K))

  # Co-membership matrisi: b_ij = i ve j'nin birlikte bulunduğu grup sayısı
  B_raw <- A %*% Matrix::t(A)

  # Köşegeni sıfırla (kendi kendisiyle co-membership sayılmaz)
  Matrix::diag(B_raw) <- 0
  B_raw <- Matrix::drop0(B_raw)

  if (length(B_raw@x) == 0) {
    B <- Matrix::sparseMatrix(
      i = integer(0), j = integer(0), x = numeric(0),
      dims = c(p, p), symmetric = TRUE
    )
    dimnames(B) <- list(features, features)
    return(B)
  }

  # Ağırlıkları [0, 1] aralığına ölçekle
  B_raw@x <- B_raw@x / max(B_raw@x)

  dimnames(B_raw) <- list(features, features)
  methods::as(B_raw, "symmetricMatrix")
}

#' Build a prior relatedness matrix `B` from an edge list or matrix
#'
#' Accepts either a 3-column edge list (`gene_i`, `gene_j`, `weight`) or an
#' already-formed symmetric weight matrix, and aligns it to `features`.
#'
#' **Duplicate and Reciprocal Edges:** If the input edge list contains reciprocal
#' directed edges (e.g., `A -> B: 0.4` and `B -> A: 0.7`) or duplicate entries,
#' they are collapsed into a single undirected edge by taking the **maximum**
#' weight (e.g., `A - B: 0.7`). This prevents artificial confidence inflation
#' caused by silently summing weights.
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
      "install.packages(\"Matrix\").",
      call. = FALSE
    )
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
        paste(utils::head(rownames(x), 5), collapse = ", "),
        call. = FALSE
      )
    }
    sub <- as.matrix(x[common, common, drop = FALSE])
    B <- Matrix::Matrix(0, p, p, sparse = TRUE)
    dimnames(B) <- list(features, features)
    B[common, common] <- sub
    B <- (B + Matrix::t(B)) / 2 # enforce symmetry
    Matrix::diag(B) <- 0
    B <- pmin(pmax(B, 0), 1)
    return(methods::as(B, "symmetricMatrix"))
  }

  # Edge-list input.
  df <- as.data.frame(x, stringsAsFactors = FALSE)
  if (ncol(df) < 3) {
    stop("Prior edge list must have at least 3 columns: gene_i, gene_j, weight.",
      call. = FALSE
    )
  }
  gi <- as.character(df[[1]])
  gj <- as.character(df[[2]])
  w <- as.numeric(df[[3]])
  keep <- gi %in% features & gj %in% features & is.finite(w)
  if (!any(keep)) {
    stop("No prior edges reference features present in X. Check that the edge ",
      "list gene IDs match rownames(X).",
      call. = FALSE
    )
  }
  gi <- gi[keep]
  gj <- gj[keep]
  w <- pmin(pmax(w[keep], 0), 1)

  # Map to integer indices
  gi_idx <- idx[gi]
  gj_idx <- idx[gj]

  # Standardize edges to be undirected (i <= j) to handle reciprocals safely
  swap <- gi_idx > gj_idx
  temp <- gi_idx[swap]
  gi_idx[swap] <- gj_idx[swap]
  gj_idx[swap] <- temp

  # Resolve duplicates and reciprocal edges by taking the MAXIMUM weight.
  # (sparseMatrix naturally sums duplicate (i,j) pairs, which artificially
  # inflates confidence and requires an arbitrary cap. Max is biologically sound).
  ord <- order(w, decreasing = TRUE)
  gi_idx <- gi_idx[ord]
  gj_idx <- gj_idx[ord]
  w <- w[ord]

  edge_id <- as.numeric(gi_idx) + (as.numeric(gj_idx) - 1) * p
  first_occur <- !duplicated(edge_id)

  B <- Matrix::sparseMatrix(
    i = gi_idx[first_occur], j = gj_idx[first_occur],
    x = w[first_occur], dims = c(p, p), symmetric = TRUE
  )
  Matrix::diag(B) <- 0
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
        "keeping all features.",
        call. = FALSE
      )
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
#' Splits samples by class and estimates `Theta0`, `Theta1`.
#'
#' **Three distinct backends:**
#'
#' - **`"glasso"` — Knowledge-weighted graphical lasso** (`O(p^3)`):
#'   Uses a per-edge penalty `lambda_ij = lambda * (1 - beta * b_ij)` so that
#'   biologically related edges (high `b_ij` in the prior `B`) are penalized
#'   less. This is the *only* backend where `B` directly shapes the estimated
#'   networks during estimation.  Because of the cubic cost, a **safety guard**
#'   caps `p` at 250 features; if `p > 250` an automatic fallback to `"mb"` is
#'   triggered with a warning.
#'
#' - **`"jgl"` — Joint fused graphical lasso** (requires the `JGL` package):
#'   Applies scalar penalties `lambda1 = sparsity_penalty` and
#'   `lambda2 = fused_penalty` to jointly estimate both class networks with a
#'   fused coupling term.  The JGL API does *not* accept per-edge penalty
#'   matrices, so `B` is **not** used during estimation; it influences the
#'   result only through the downstream network-integration and module-scoring
#'   stages.  A scalability guard skips JGL when `p > 500` or `p > 2*n`
#'   (falls back to `"glasso"` or `"mb"`).
#'
#' - **`"mb"` — Meinshausen-Buhlmann neighborhood selection** (`~O(p^2)`,
#'   via `huge`): The fastest backend, suitable for large `p`.  Uses a scalar
#'   penalty; `B` is not applied during estimation.
#'
#' @param X0,X1 Sample-by-feature numeric matrices for class 0 and class 1
#'   (rows = samples, columns = features). Note this is the transpose of the
#'   `genes x samples` convention used elsewhere; [bridge_grouping()] handles the
#'   transpose for you.
#' @param B Symmetric prior matrix over the features (columns of `X0`/`X1`),
#'   values in `[0, 1]`. May be `NULL` for `B = 0`. Used in per-edge penalty
#'   weighting for the `glasso` backend only; passed through to
#'   `bridge_extract_modules` and `bridge_score_modules` for all backends.
#' @param sparsity_penalty Base sparsity penalty (previously `lambda`).
#' @param beta Knowledge-weighting strength in `[0, 1)`. Only effective for
#'   the `glasso` backend; ignored by JGL and mb.
#' @param fused_penalty Fused penalty coupling the two networks (JGL backend only, previously `lambda2`).
#' @param data_type One of `"gaussian"`, `"nonparanormal"`, `"compositional"`.
#' @param precision_method Estimation backend: `"glasso"` (knowledge-weighted
#'   graphical lasso, `O(p^3)`, max ~250 genes), `"jgl"` (joint fused graphical
#'   lasso, requires JGL package), or `"mb"` (Meinshausen-Buhlmann neighborhood
#'   selection via `huge`, `~O(p^2)`, fastest for large `p`).
#' @param glasso_max_p Safety cap for the `glasso` backend. If `p` exceeds this
#'   value, the backend automatically falls back to `"mb"` with a warning.
#'   Default `250`.
#' @return A list with `net0`, `net1` (dense numeric matrices), the `backend`
#'   used, and `type` (`"precision"` or `"adjacency"`).
#' @keywords internal
bridge_estimate_precision <- function(X0, X1, B, sparsity_penalty = 0.1, beta = 0.5,
                                      fused_penalty = 0.05,
                                      data_type = "gaussian",
                                      precision_method = "glasso",
                                      glasso_max_p = 250L,
                                      progress = NULL) {
  pg <- function(lbl) if (!is.null(progress)) progress(lbl)
  data_type <- match.arg(
    data_type,
    c("gaussian", "nonparanormal", "compositional")
  )
  precision_method <- match.arg(precision_method, c("glasso", "jgl", "mb"))
  if (data_type == "compositional") {
    stop("data_type = \"compositional\" is not yet supported. It needs a CLR + ",
      "SPIEC-EASI estimator (the 'SpiecEasi' package, currently GitHub-only: ",
      "remotes::install_github(\"zdk123/SpiecEasi\")). Use \"gaussian\" or ",
      "\"nonparanormal\" for now.",
      call. = FALSE
    )
  }

  p <- ncol(X0)
  n0 <- nrow(X0)
  n1 <- nrow(X1)
  n <- n0 + n1

  # ── Dimension / scalability report ─────────────────────────────────────────
  # Always emit this so the user can diagnose hangs without verbose=TRUE in the
  # parent call. p/n >> 1 is the primary cause of estimation timeouts.
  message(sprintf(
    "BRIDGE [estimation]: p=%d  n0=%d  n1=%d  n=%d  p/n=%.1f  backend_pref=%s",
    p, n0, n1, n, p / n, precision_method
  ))

  if (data_type == "nonparanormal") {
    if (!requireNamespace("huge", quietly = TRUE)) {
      stop("data_type = \"nonparanormal\" needs the 'huge' package. Install it ",
        "with install.packages(\"huge\").",
        call. = FALSE
      )
    }
    X0 <- huge::huge.npn(X0, verbose = FALSE)
    X1 <- huge::huge.npn(X1, verbose = FALSE)
  }

  # ══════════════════════════════════════════════════════════════════════════
  # Backend 1: mb — Meinshausen-Buhlmann neighborhood selection (O(p^2))
  # Fastest backend, suitable for large p. Scalar penalty; B is NOT used.
  # ══════════════════════════════════════════════════════════════════════════
  if (precision_method == "mb") {
    if (!requireNamespace("huge", quietly = TRUE)) {
      stop("precision_method = \"mb\" needs the 'huge' package. Install it ",
        "with install.packages(\"huge\").",
        call. = FALSE
      )
    }
    est_mb <- function(Xc) {
      fit <- huge::huge(as.matrix(Xc),
        lambda = sparsity_penalty, method = "mb",
        verbose = FALSE
      )
      # huge returns a path (adjacency); with a single lambda take the first (only) path.
      P <- fit$path[[1]]
      Adj <- as.matrix(P) * 1.0
      Adj <- (Adj + t(Adj)) / 2 # symmetrize the adjacency
      diag(Adj) <- 1
      dimnames(Adj) <- list(colnames(Xc), colnames(Xc))
      Adj
    }
    pg("class 1/2 network (mb)…")
    A0 <- est_mb(X0)
    pg("class 2/2 network (mb)…")
    A1 <- est_mb(X1)
    return(list(net0 = A0, net1 = A1, backend = "huge-mb", type = "adjacency"))
  }

  # ══════════════════════════════════════════════════════════════════════════
  # Backend 2: jgl — Joint Fused Graphical Lasso (requires JGL package)
  # Scalar penalties lambda1/lambda2; B is NOT used during estimation.
  # Scalability guard: skipped when p > 500 or p > 2*n (O(p^3) SDP hangs).
  # ══════════════════════════════════════════════════════════════════════════
  if (precision_method == "jgl") {
    if (!requireNamespace("JGL", quietly = TRUE)) {
      stop("precision_method = \"jgl\" requires the 'JGL' package. Install it ",
        "with install.packages(\"JGL\").",
        call. = FALSE
      )
    }
    jgl_too_large <- (p > 500L) || (p > 2L * n)
    if (jgl_too_large) {
      warning(
        sprintf(
          "BRIDGE: JGL requested but p=%d exceeds safe threshold (500 or 2*n=%d). ",
          p, 2L * n
        ),
        "Falling back to 'mb' to avoid hangs.",
        call. = FALSE
      )
      # Recursive call with mb
      return(bridge_estimate_precision(X0, X1, B,
        sparsity_penalty = sparsity_penalty, beta = beta,
        fused_penalty = fused_penalty, data_type = data_type,
        precision_method = "mb", glasso_max_p = glasso_max_p,
        progress = progress
      ))
    }
    pg("joint graphical lasso (JGL)…")
    fit <- tryCatch(
      JGL::JGL(list(as.matrix(X0), as.matrix(X1)),
        penalty = "fused", lambda1 = sparsity_penalty, lambda2 = fused_penalty,
        weights = "equal", return.whole.theta = TRUE
      ),
      error = function(e) {
        warning("BRIDGE: JGL failed with error: ", conditionMessage(e),
          ". Falling back to 'mb'.",
          call. = FALSE
        )
        NULL
      }
    )
    if (!is.null(fit)) {
      return(list(
        net0 = fit$theta[[1]], net1 = fit$theta[[2]],
        backend = "JGL", type = "precision"
      ))
    }
    # JGL failed; fall back to mb
    return(bridge_estimate_precision(X0, X1, B,
      sparsity_penalty = sparsity_penalty, beta = beta,
      fused_penalty = fused_penalty, data_type = data_type,
      precision_method = "mb", glasso_max_p = glasso_max_p,
      progress = progress
    ))
  }

  # ══════════════════════════════════════════════════════════════════════════
  # Backend 3: glasso — Knowledge-weighted graphical lasso (O(p^3))
  # Per-edge penalty: lambda_ij = lambda * (1 - beta * b_ij)
  # This is the ONLY backend where B directly shapes estimation.
  # SAFETY GUARD: p > glasso_max_p triggers automatic fallback to mb.
  # ══════════════════════════════════════════════════════════════════════════
  if (p > glasso_max_p) {
    warning(
      sprintf(
        "BRIDGE: glasso requested but p=%d exceeds safety cap (glasso_max_p=%d). ",
        p, glasso_max_p
      ),
      "The O(p^3) glasso solver would freeze or run for hours at this scale. ",
      "Falling back to 'mb' (fast). Reduce gene count with prefilter/max_genes ",
      "to stay within the glasso cap.",
      call. = FALSE
    )
    return(bridge_estimate_precision(X0, X1, B,
      sparsity_penalty = sparsity_penalty, beta = beta,
      fused_penalty = fused_penalty, data_type = data_type,
      precision_method = "mb", glasso_max_p = glasso_max_p,
      progress = progress
    ))
  }
  if (!requireNamespace("glasso", quietly = TRUE)) {
    stop("precision_method = \"glasso\" requires the 'glasso' package. ",
      "Install it with install.packages(\"glasso\").",
      call. = FALSE
    )
  }
  message(sprintf("BRIDGE [glasso]: building knowledge-weighted penalty matrix (p=%d, beta=%.2f)", p, beta))
  Bm <- if (is.null(B)) matrix(0, p, p) else as.matrix(B) # dense only here
  Wpen <- sparsity_penalty * (1 - beta * Bm)
  Wpen[Wpen < 0] <- 0
  est_one <- function(Xc) {
    S <- stats::cov(Xc)
    fit <- glasso::glasso(S, rho = Wpen, penalize.diagonal = FALSE)
    Theta <- fit$wi
    dimnames(Theta) <- list(colnames(Xc), colnames(Xc))
    Theta
  }
  pg("class 1/2 network (glasso)…")
  T0 <- est_one(X0)
  pg("class 2/2 network (glasso)…")
  T1 <- est_one(X1)
  list(net0 = T0, net1 = T1, backend = "glasso", type = "precision")
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
#' **Role of `B` here (network integration stage):**
#' `B` contributes to the blended weight matrix `W`. To prevent scale imbalances
#' between the data-driven differential signal `|Delta|` and the prior `B`, both
#' are max-normalized before being mixed via a convex combination controlled by
#' `prior_weight`. This ensures the prior neither dominates nor is ignored
#' due to arbitrary numeric scales.
#'
#' @param Theta0,Theta1 Precision matrices from [bridge_estimate_precision()].
#' @param B Symmetric prior matrix (or `NULL`).
#' @param resolution Numeric vector of Leiden resolutions to sweep.
#' @param prior_weight Numeric in `[0, 1]`. Weight of the normalized prior graph
#'   in the blended differential network. Default `0.25`.
#' @param edge_tol Numeric threshold to suppress dense solver artifacts. Absolute
#'   precision matrix entries below this are treated as zero. Default `1e-4`.
#' @param details If `TRUE`, attach the normalized blended graph `W` (attribute
#'   `"graph_W"`) and the per-resolution module counts (attribute `"res_counts"`)
#'   for inspection.
#' @return A named list of character vectors: candidate modules across all
#'   resolutions (names encode the resolution and cluster id).
#' @keywords internal
bridge_extract_modules <- function(net0, net1, B, resolution,
                                   prior_weight = 0.25, edge_tol = 1e-4,
                                   details = FALSE, progress = NULL) {
  feats <- rownames(net0)
  p <- length(feats)
  Delta <- abs(net1 - net0)

  # Suppress numerical solver artifacts (glasso/JGL typical tol is 1e-4)
  retained <- (abs(net0) > edge_tol) | (abs(net1) > edge_tol) # union support
  Bm <- if (is.null(B)) matrix(0, p, p) else as.matrix(B)

  # Blend: scale-normalized convex combination of differential signal and prior
  D_eff <- Delta * retained
  D_scale <- if (max(D_eff) > 0) max(D_eff) else 1
  B_scale <- if (max(Bm) > 0) max(Bm) else 1

  W <- (1 - prior_weight) * (D_eff / D_scale) + prior_weight * (Bm / B_scale)

  W <- (W + t(W)) / 2 # enforce exact symmetry (glasso output can drift)
  diag(W) <- 0
  W[!is.finite(W)] <- 0
  # Normalize to [0, 1] for stable community detection.
  if (max(W) > 0) W <- W / max(W)
  dimnames(W) <- list(feats, feats)

  modules <- list()
  res_counts <- integer(0)
  for (ri in seq_along(resolution)) {
    res <- resolution[ri]
    if (!is.null(progress)) {
      progress(sprintf(
        "community detection: resolution %s (%d/%d)…",
        format(res, trim = TRUE), ri, length(resolution)
      ))
    }
    memb <- detect_communities(W, method = "leiden", resolution = res)
    n_here <- 0L
    for (cl in unique(memb)) {
      genes <- names(memb)[memb == cl]
      if (length(genes) >= 2) {
        prefix <- if (is.null(B)) "DataDrivenBRIDGE" else "PriorInformedBRIDGE"
        key <- sprintf("%s_r%s_c%s", prefix, format(res, trim = TRUE), cl)
        modules[[key]] <- sort(genes)
        n_here <- n_here + 1L
      }
    }
    res_counts[format(res, trim = TRUE)] <- n_here
  }
  if (details) {
    attr(modules, "graph_W") <- W
    attr(modules, "res_counts") <- res_counts
  }
  modules
}

# ── Grouping-Scoring coupling ────────────────────────────────────────────────

#' @keywords internal
bridge_eigengene <- function(genes, X_gs) {
  sub <- X_gs[genes, , drop = FALSE]
  if (nrow(sub) < 1 || ncol(sub) < 4) {
    return(rep(0, ncol(X_gs)))
  }
  m <- t(sub)
  m <- scale(m, center = TRUE, scale = TRUE)
  m[!is.finite(m)] <- 0
  tryCatch(stats::prcomp(m, center = FALSE, scale. = FALSE)$x[, 1],
    error = function(e) rowMeans(m)
  )
}

#' @keywords internal
fast_discrimination <- function(pc, y, k = 3) {
  n <- length(pc)
  folds <- ((sample(n) - 1) %% k) + 1
  vals <- vapply(seq_len(k), function(f) {
    test <- folds == f
    if (length(unique(y[test])) < 2) {
      return(NA_real_)
    }
    lv <- levels(y)
    a <- pc[test & y == lv[1]]
    b <- pc[test & y == lv[2]]
    if (length(a) < 1 || length(b) < 1) {
      return(NA_real_)
    }
    s <- stats::sd(pc[test])
    if (!is.finite(s) || s == 0) s <- 1
    abs(mean(a) - mean(b)) / s
  }, numeric(1))
  mean(vals, na.rm = TRUE)
}

#' Module eigengene discrimination D(M)
#'
#' Evaluates module discrimination without data leakage (leakage-free CV):
#' For each fold, PCA (mean, variance, and PC1 loading vector) is fit strictly
#' on the training samples. Held-out test samples are then projected onto that
#' learned eigenspace, and class discrimination is computed on the out-of-fold test samples.
#'
#' @keywords internal
bridge_discrimination <- function(genes, X_gs, y, k = 3) {
  sub <- X_gs[genes, , drop = FALSE]
  p_genes <- nrow(sub)
  n <- ncol(sub)
  if (p_genes < 1 || n < 4) {
    return(0)
  }

  y <- as.factor(y)
  lv <- levels(y)
  if (length(lv) < 2) {
    return(0)
  }

  tab <- table(y)
  if (min(tab) < 2) {
    return(0)
  }

  k_eff <- min(k, min(tab))
  if (k_eff < 2) {
    return(0)
  }

  # Stratified fold partition: her katmanda her sınıftan dengeli örnek bulunsun
  folds <- integer(n)
  for (l in lv) {
    idx_l <- which(y == l)
    folds[idx_l] <- ((sample(length(idx_l)) - 1) %% k_eff) + 1
  }

  vals <- vapply(seq_len(k_eff), function(f) {
    test_idx <- folds == f
    train_idx <- !test_idx

    y_train <- y[train_idx]
    y_test <- y[test_idx]

    if (length(unique(y_train)) < 2 || length(unique(y_test)) < 2) {
      return(NA_real_)
    }

    m_train <- t(sub[, train_idx, drop = FALSE])
    m_test  <- t(sub[, test_idx, drop = FALSE])

    # 1. TRAIN: Parametreleri (mean, sd) SADECE train'den öğren (Data leakage önleme)
    train_means <- colMeans(m_train)
    train_sds <- apply(m_train, 2, stats::sd)
    train_sds[!is.finite(train_sds) | train_sds == 0] <- 1

    m_train_sc <- scale(m_train, center = train_means, scale = train_sds)
    m_train_sc[!is.finite(m_train_sc)] <- 0

    # 2. PCA FIT: SADECE train üzerinde PC1 özvektörünü (rotation) çıkar
    rot <- tryCatch({
      fit <- stats::prcomp(m_train_sc, center = FALSE, scale. = FALSE)
      fit$rotation[, 1]
    }, error = function(e) NULL)

    if (is.null(rot) || length(rot) != p_genes) {
      rot <- rep(1 / sqrt(max(1, p_genes)), p_genes)
    }

    # 3. TEST PROJECTION: Test verisini TRAIN parametreleriyle ölçekle ve projekte et
    m_test_sc <- scale(m_test, center = train_means, scale = train_sds)
    m_test_sc[!is.finite(m_test_sc)] <- 0

    pc_test <- as.vector(m_test_sc %*% rot)

    # 4. DISCRIMINATION: Out-of-fold test projeksiyonu üzerinde standardize sınıf farkı
    a <- pc_test[y_test == lv[1]]
    b <- pc_test[y_test == lv[2]]
    if (length(a) < 1 || length(b) < 1) {
      return(NA_real_)
    }

    s <- stats::sd(pc_test)
    if (!is.finite(s) || s == 0) s <- 1
    abs(mean(a) - mean(b)) / s
  }, numeric(1))

  res <- mean(vals, na.rm = TRUE)
  if (!is.finite(res)) 0 else res
}

#' Mean prior cohesion C(M)
#'
#' Measures how strongly a module's genes are connected in the prior graph `B`.
#' This is the module-selection (network integration) use of `B`: modules whose
#' genes share high prior co-membership receive a higher `C(M)`, boosting their
#' `J(M)` score proportionally to `lambda1`.
#' @keywords internal
bridge_cohesion <- function(genes, B) {
  if (is.null(B) || length(genes) < 2) {
    return(0)
  }
  sub <- as.matrix(B[genes, genes, drop = FALSE])
  m <- sub[upper.tri(sub)]
  if (length(m) == 0) {
    return(0)
  }
  mean(m)
}

#' Score, stability-filter and select candidate modules
#'
#' Computes `J(M) = D(M) + lambda1 * C(M) - lambda2 * |M|` for each module,
#' filters by a bootstrap-stability threshold `R(M) >= tau`, and keeps modules
#' with positive objective. Deduplicates identical gene sets, preferring the
#' higher-scoring one.
#'
#' **Role of `B` here (network integration stage):**
#' `C(M)` (prior cohesion) measures the mean prior edge weight among a module's
#' genes. High `C(M)` boosts `J(M)`, so biologically coherent modules (according
#' to `B`) are preferred regardless of which estimation backend was used.
#'
#' @param modules Named list of candidate modules.
#' @param X_gs Genes-by-samples expression matrix.
#' @param y Class-label factor.
#' @param B Prior matrix (or `NULL`). Used to compute `C(M)` (network
#'   integration stage, not estimation).
#' @param cohesion_weight Weight on prior cohesion in `J(M)` (previously `lambda1`).
#' @param size_penalty Size penalty in `J(M)` (previously `lambda2`).
#' @param tau Bootstrap-stability threshold in `[0, 1]`.
#' @param n_boot Number of bootstrap resamples for stability.
#' @param max_overlap Maximum allowed Jaccard similarity between modules before deduplication.
#' @param details If `TRUE`, attach a per-module `"diagnostics"` data.frame
#'   (`module, size, J, D, C, stability`) alongside `"scores"`.
#' @return Named list of the selected modules (character vectors), with a
#'   `"scores"` attribute (the `J(M)` values) and, when `details = TRUE`, a
#'   `"diagnostics"` data.frame.
#' @keywords internal
bridge_score_modules <- function(modules, X_gs, y, B, cohesion_weight = 0.25,
                                 size_penalty = 0.05, tau = 0.6, n_boot = 20,
                                 max_overlap = 0.85, details = FALSE) {
  empty <- structure(list(), scores = numeric(0))
  if (details) {
    attr(empty, "diagnostics") <-
      data.frame(
        module = character(0), size = integer(0), J = numeric(0),
        D = numeric(0), C = numeric(0), stability = numeric(0),
        stringsAsFactors = FALSE
      )
  }
  if (length(modules) == 0) {
    return(empty)
  }
  n <- ncol(X_gs)

  # D(M) and C(M) computed ONCE per module; J is derived from them (so any
  # diagnostics table is consistent with the J it explains).
  D <- vapply(
    modules, function(genes) bridge_discrimination(genes, X_gs, y),
    numeric(1)
  )
  C <- vapply(modules, function(genes) bridge_cohesion(genes, B), numeric(1))
  sizes <- vapply(modules, length, integer(1))
  J <- D + cohesion_weight * C - size_penalty * sizes / 100

  # HIZLI VE GEÇERLİ STABİLİTE ALTERNATİFİ:
  # JGL/glasso tabanlı tüm ağ tahminini (network estimation) bootstrap içinde
  # yüzlerce kez tekrarlayıp Jaccard similarity ölçmek (computational complexity)
  # açısından çok maliyetlidir (saatler/günler sürer).
  # Eski koddaki D > 0 kontrolü ise 'D' mutlak (absolute) bir fark olduğu için
  # her zaman TRUE dönüyor ve sahte bir stabilite (%100) üretiyordu.
  #
  # YENİ METODOLOJİ (Permutation-Anchored Bootstrap Stability):
  # Bootstrap sırasında ağ yeniden oluşturulmaz, ancak modülün ayırıcılık gücünün
  # (discrimination) o bootstrap örneğinde rastgele bir şans olup olmadığı test edilir.
  # Modülün o örneklemdeki gücü (D_boot), etiketleri rastgele karıştırılmış
  # halinden (D_null) büyükse modül 'stabil/güvenilir' kabul edilir.
  # RUNTIME OPTİMİZASYONU: PCA işlemini bootstrap döngüsünden çıkar.
  # PCA projeksiyonunu (eigengene) bir kez hesapla, bootstrap'lerde yeniden kullan.
  # Bu optimizasyon karmaşıklığı O(M * B * genes^2) seviyesinden O(M * B * N) seviyesine
  # indirir (dakikalar sürecek bir işlemi milisaniyelere çeker).
  PCs <- lapply(modules, function(genes) bridge_eigengene(genes, X_gs))

  stability <- vapply(seq_along(modules), function(i) {
    pc_full <- PCs[[i]]
    if (all(pc_full == 0)) {
      return(0)
    }

    ok <- vapply(seq_len(n_boot), function(b) {
      idx <- sample(n, n, replace = TRUE)
      y_boot <- y[idx]
      if (length(unique(y_boot)) < 2) {
        return(NA)
      }

      pc_boot <- pc_full[idx]

      # Orijinal etiketler ile diskriminasyon
      D_boot <- fast_discrimination(pc_boot, y_boot)

      # Rastgele etiketler ile null (şans eseri) diskriminasyon
      y_null <- sample(y_boot)
      D_null <- fast_discrimination(pc_boot, y_null)

      # Gerçek sinyal, şans eseri oluşan gürültüyü geçiyor mu?
      D_boot > D_null
    }, logical(1))
    mean(ok, na.rm = TRUE)
  }, numeric(1))

  keep <- is.finite(J) & J > 0 & stability >= tau
  if (!any(keep)) {
    return(empty)
  }
  sel <- modules[keep]
  Jk <- J[keep]
  Dk <- D[keep]
  Ck <- C[keep]
  Sk <- stability[keep]

  # Deduplicate highly similar/nested modules from different resolutions.
  # We iterate in decreasing order of J(M) score (highest quality first).
  # A candidate is dropped if it has > max_overlap Jaccard similarity with an already-kept module.
  ord <- order(Jk, decreasing = TRUE)
  sel <- sel[ord]
  Jk <- Jk[ord]
  Dk <- Dk[ord]
  Ck <- Ck[ord]
  Sk <- Sk[ord]

  if (length(sel) > 1) {
    keep_idx <- c(1L)
    for (i in seq_len(length(sel))[-1]) {
      g1 <- sel[[i]]
      is_dup <- FALSE
      for (j in keep_idx) {
        g2 <- sel[[j]]
        intersect_len <- length(intersect(g1, g2))
        union_len <- length(g1) + length(g2) - intersect_len
        if (intersect_len / union_len > max_overlap) {
          is_dup <- TRUE
          break
        }
      }
      if (!is_dup) keep_idx <- c(keep_idx, i)
    }
    sel <- sel[keep_idx]
    Jk <- Jk[keep_idx]
    Dk <- Dk[keep_idx]
    Ck <- Ck[keep_idx]
    Sk <- Sk[keep_idx]
  }

  out <- structure(sel, scores = Jk)
  if (details) {
    # Extract resolution from module names (e.g., "DataDrivenBRIDGE_r1.5_c3" -> 1.5)
    res_vals <- as.numeric(gsub("^.*BRIDGE_r([0-9.]+)_c.*$", "\\1", names(sel)))
    attr(out, "diagnostics") <- data.frame(
      module = names(sel), resolution = res_vals, size = vapply(sel, length, integer(1)),
      J = Jk, D = Dk, C = Ck, stability = Sk,
      row.names = NULL, stringsAsFactors = FALSE
    )
  }
  out
}

# ── Public entry point ───────────────────────────────────────────────────────

#' Derive feature groups with BRIDGE
#'
#' Runs the full BRIDGE pipeline: class-conditional precision estimation,
#' differential-network module extraction across resolutions, and
#' Grouping-Scoring coupled module selection. Returns groups in the exact format
#' the Scoring step consumes: a **named list of character vectors**
#' (`group name -> gene symbols`).
#'
#' **How prior `B` is used across the pipeline:**
#'
#' `B` participates in two distinct stages:
#'
#' 1. **Network estimation** (glasso backend only): Per-edge penalty
#'    `lambda_ij = lambda * (1 - beta * b_ij)` downweights penalization of
#'    biologically related edges. This is the knowledge-weighted graphical lasso.
#'    *Not available in the JGL or mb backends* (their APIs accept only scalar
#'    penalties). Use `precision_method = "glasso"` and ensure JGL is not
#'    installed if you need this behaviour.
#'
#' 2. **Network integration and module selection** (all backends): `B` enters
#'    the blended differential graph `W` (via a scale-normalized convex combination
#'    controlled by `prior_weight` in [bridge_extract_modules()]) and the cohesion
#'    term `C(M)` of the module objective `J(M)` (in [bridge_score_modules()]).
#'    These effects are backend-independent.
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
#'   `B` is used in the network-integration and module-selection stages for all
#'   backends; it additionally shapes the estimation penalty when
#'   `precision_method = "glasso"`.
#' @param sparsity_penalty Base graphical-lasso sparsity penalty (previously `lambda`). Default `0.1`.
#'   **Note:** We use a scalar value rather than a sequence/path (e.g., `huge.select`)
#'   because model selection in BRIDGE is driven by the downstream differential module
#'   objective `J(M)` and stability `R(M)`, not by unsupervised network-level criteria
#'   like EBIC or StARS. This is a tunable methodological hyper-parameter.
#' @param beta Knowledge-weighting strength in `[0, 1)`; higher = prior edges
#'   penalized less during glasso estimation. Only effective when the `"glasso"`
#'   backend is active. Default `0.5`.
#' @param fused_penalty Fused penalty coupling the two class networks (`"jgl"` backend
#'   only, previously `lambda2`). Default `0.05`.
#' @param cohesion_weight Weight on prior cohesion `C(M)` in the selection objective
#'   `J(M)` (previously `lambda1`). Default `0.25` is an implementation heuristic.
#' @param size_penalty Size penalty in the selection objective `J(M)` (previously overloaded with `lambda2`).
#'   Default `0.05` is an implementation heuristic.
#' @param prior_weight Weight of the normalized prior in the blended differential
#'   graph `W = (1 - w)*Delta_norm + w*B_norm` during network integration.
#'   Default `0.25`. This value is a tunable hyper-parameter, not a biological threshold.
#' @param edge_tol Numeric threshold for edge selection. Absolute precision matrix
#'   entries below this are treated as zero. Graphical lasso solvers (glasso, JGL)
#'   often leave tiny numerical artifacts (e.g., `1e-6`) due to `1e-4` convergence
#'   tolerances, which can cause spurious dense networks if checked via `!= 0`.
#'   Default `1e-4`.
#' @param resolution Numeric vector of Leiden resolutions to sweep for
#'   multi-scale modules. Default `c(0.5, 1, 1.5, 2)`. **Note:** This default grid is
#'   NOT an automatically 'optimal' biological set; it is merely a computational
#'   starting point. Users should change these values to sweep different granularities
#'   and compare the resulting module scores (`J(M)`) in the diagnostics output.
#' @param tau Bootstrap-stability threshold `R(M) >= tau` in `[0, 1]`.
#'   Default `0.6`. This is a tunable heuristic, not an absolute truth.
#' @param max_overlap Maximum allowed Jaccard similarity `[0, 1]` between candidate modules.
#'   Modules exceeding this overlap are deduplicated to prevent redundant nested modules,
#'   keeping only the one with the highest `J(M)` score. Default `0.85`.
#' @param n_boot Number of permutation-anchored bootstrap resamples to compute
#'   module stability. Default `20`. Use a small value (e.g., `20` or `50`) for fast
#'   interactive exploration or debugging. For final publication-quality
#'   analysis, increase to `100` or `500`. Note that higher `n_boot` linearly
#'   increases the runtime of the scoring phase.
#' @param data_type One of `"gaussian"` (default), `"nonparanormal"`
#'   (applies [huge::huge.npn()] first), or `"compositional"` (not yet
#'   supported; errors with guidance to install SpiecEasi).
#' @param precision_method Estimation backend: `"glasso"` (knowledge-weighted
#'   graphical lasso, uses prior `B`, max ~250 genes), `"jgl"` (joint fused
#'   graphical lasso, requires JGL package), or `"mb"` (fast neighborhood
#'   selection via `huge`, best for large `p`). See
#'   [bridge_estimate_precision()].
#' @param prefilter Label-blind feature pre-filter applied first: `"none"`
#'   (default), `"variance"`, or `"prior"`. See [bridge_prefilter()]. This never
#'   inspects `y`, so it is safe to apply outside cross-validation.
#' @param max_genes Optional cap on the number of features kept by the variance
#'   pre-filter. `NULL` (default) = no cap.

#' @param min_size,max_size Keep only modules whose size is in
#'   `[min_size, max_size]`. Defaults `5` and `500`. These are implementation constraints
#'   to prevent downstream crashes (e.g., matching `gsm_dataset()`) and filter out
#'   dense artifactual giant components; they are not strictly biological thresholds.
#' @param seed Random seed for reproducibility. Default `1`.
#' @param verbose Emit progress messages. Default `TRUE`.
#' @param progress_file Optional path; when set, a `"step|label"` stage marker
#'   (steps 1-5) is written at each phase so a watcher (the Shiny Grouping tab)
#'   can show a live progress bar. `NULL` (default) is a no-op.
#' @param details If `TRUE`, attach inspection attributes for the interactive
#'   Grouping tab: `"diagnostics"` (per-module `J/D/C/stability/size`),
#'   `"summary"` (candidate/prefilter/backend/prior counts, per-resolution
#'   sweep), and `"graph"` (a capped edge list of the differential network among
#'   selected genes). Default `FALSE` keeps the per-fold cross-validation path
#'   fast and its output unchanged.
#' @param graph_max_nodes,graph_max_edges Caps on the `"graph"` edge list emitted
#'   under `details = TRUE`. The defaults keep the interactive network plot
#'   responsive. The `"gnn"` classifier raises them (via
#'   [grouping_for_fold()]`(with_graph = TRUE)`) because it convolves over the
#'   graph and a truncated one would leave genes with no edges.
#' @return A named list of character vectors (`group name -> gene symbols`),
#'   identical in shape to [read_groups()] output. Carries a `"scores"`
#'   attribute with each module's `J(M)`. If no module passes the filters, a
#'   warning is issued and an empty named list is returned.
#' @examples
#' \donttest{
#' set.seed(1)
#' # 40 genes x 60 samples, two classes; genes 1-8 are a differential module.
#' p <- 40
#' n <- 60
#' X <- matrix(stats::rnorm(p * n),
#'   nrow = p,
#'   dimnames = list(paste0("G", seq_len(p)), paste0("S", seq_len(n)))
#' )
#' y <- factor(rep(c("A", "B"), each = n / 2))
#' mod <- paste0("G", 1:8)
#' X[mod, y == "B"] <- X[mod, y == "B"] + 2 # shift the module in class B
#' # Prior: the true module is known to be related.
#' groups_prior <- list(known_module = mod)
#' if (requireNamespace("glasso", quietly = TRUE) ||
#'   requireNamespace("JGL", quietly = TRUE)) {
#'   g <- bridge_grouping(X, y, B = groups_prior, resolution = c(1, 1.5))
#'   str(g)
#' }
#' }
#' @export
bridge_grouping <- function(X, y, B = NULL,
                            sparsity_penalty = 0.1, fused_penalty = 0.05,
                            cohesion_weight = 0.25, size_penalty = 0.05,
                            beta = 0.5, prior_weight = 0.25,
                            edge_tol = 1e-4,
                            resolution = c(0.5, 1, 1.5, 2),
                            tau = 0.6, n_boot = 20, max_overlap = 0.85,
                            data_type = c(
                              "gaussian", "nonparanormal",
                              "compositional"
                            ),
                            precision_method = c("glasso", "jgl", "mb"),
                            prefilter = c("none", "variance", "prior"),
                            max_genes = NULL,
                            min_size = 5, max_size = 500,
                            seed = 1, verbose = TRUE, details = FALSE,
                            graph_max_nodes = 150L, graph_max_edges = 2000L,
                            progress_file = NULL,
                            # Backward compatibility for old naming:
                            lambda = NULL, lambda1 = NULL, lambda2 = NULL) {
  if (!is.null(lambda)) {
    warning("BRIDGE: 'lambda' is deprecated; mapping to 'sparsity_penalty'.", call. = FALSE)
    sparsity_penalty <- lambda
  }
  if (!is.null(lambda1)) {
    warning("BRIDGE: 'lambda1' is deprecated; mapping to 'cohesion_weight'.", call. = FALSE)
    cohesion_weight <- lambda1
  }
  if (!is.null(lambda2)) {
    warning("BRIDGE: 'lambda2' is deprecated; mapping to 'fused_penalty' AND 'size_penalty'. ",
      "Please use the new separate arguments.",
      call. = FALSE
    )
    fused_penalty <- lambda2
    size_penalty <- lambda2
  }

  data_type <- match.arg(data_type)
  precision_method <- match.arg(precision_method)
  prefilter <- match.arg(prefilter)
  if (!is.null(seed)) {
    old_seed <- globalenv()$.Random.seed
    on.exit(
      {
        if (is.null(old_seed)) {
          rm(list = ".Random.seed", envir = globalenv(), inherits = FALSE)
        } else {
          assign(".Random.seed", old_seed, envir = globalenv())
        }
      },
      add = TRUE
    )
    set.seed(seed)
  }
  n_features_total <- nrow(X)

  # Stream a "step|label" stage marker so a watcher (the Shiny G tab) can show a
  # live progress bar even though this may run in a forked worker. Steps are
  # out of 5. No-op when progress_file is NULL (CLI / per-fold MCCV).
  mp <- function(step, label) {
    if (is.null(progress_file)) {
      return(invisible())
    }
    try(writeLines(paste0(step, "|", label), progress_file), silent = TRUE)
  }

  t_total_start <- Sys.time()
  t_phase <- Sys.time()
  profile_times <- list()
  profile_msg <- function(phase, details = "") {
    elapsed <- as.numeric(difftime(Sys.time(), t_phase, units = "secs"))
    profile_times[[phase]] <<- elapsed
    if (verbose) {
      msg <- sprintf("BRIDGE Profile: [%s] %.2f sec", phase, elapsed)
      if (nzchar(details)) msg <- paste0(msg, " | ", details)
      message(msg)
    }
    t_phase <<- Sys.time()
  }

  mp(1, "validating inputs")

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
      call. = FALSE
    )
  }
  if (nlevels(y) != 2) {
    stop("BRIDGE currently supports exactly two classes; found ", nlevels(y),
      ".",
      call. = FALSE
    )
  }
  # ── Label-blind pre-filter (speed; never inspects y) ───────────────────────
  mp(1, "pre-filtering genes")
  if (prefilter != "none" || !is.null(max_genes)) {
    kept <- bridge_prefilter(X, method = prefilter, max_genes = max_genes, B = B)
    if (verbose && length(kept) < nrow(X)) {
      message(
        "BRIDGE: pre-filter kept ", length(kept), "/", nrow(X),
        " features (label-blind)."
      )
    }
    X <- X[kept, , drop = FALSE]
  }
  feats <- rownames(X)
  prefilter_kept <- length(feats)
  profile_msg("Feature Filtering", sprintf("%d genes input, %d kept", n_features_total, prefilter_kept))

  # ── Resolve the prior B into a matrix aligned to feats ─────────────────────
  Bmat <- bridge_standardize_prior(B, feats, verbose = verbose)

  # ── p >> n: restrict to connected components of B when it is informative ───
  active <- bridge_active_features(X, Bmat, verbose = verbose)
  Xa <- X[active, , drop = FALSE]
  Ba <- if (is.null(Bmat)) NULL else Bmat[active, active, drop = FALSE]

  pipeline_mode <- if (is.null(Ba)) "Data-Driven Differential Grouping" else "Prior-Informed BRIDGE"
  if (verbose) message(sprintf("\n=== BRIDGE EXECUTION: %s ===", pipeline_mode))

  profile_msg("Prior Construction & Active Set", sprintf("%d active genes", length(active)))

  # ── Split by class; estimate precision matrices ────────────────────────────
  lv <- levels(y)
  X0 <- t(Xa[, y == lv[1], drop = FALSE]) # samples x features
  X1 <- t(Xa[, y == lv[2], drop = FALSE])
  n0_act <- nrow(X0)
  n1_act <- nrow(X1)
  p_act <- ncol(X0)

  # ── Sample size safety checks ──────────────────────────────────────────────
  if (n0_act < 3 || n1_act < 3) {
    stop(sprintf(
      "BRIDGE requires at least 3 samples per class for network estimation (found n0=%d, n1=%d).",
      n0_act, n1_act
    ), call. = FALSE)
  }
  if (n0_act < 10 || n1_act < 10) {
    warning(sprintf(
      "BRIDGE: Sample sizes are very small (n0=%d, n1=%d). Network estimation may be highly unstable.",
      n0_act, n1_act
    ), call. = FALSE)
  }

  if (verbose) {
    message(sprintf(
      "BRIDGE: estimating class-conditional networks — p=%d  n0=%d  n1=%d  p/n=%.1f",
      p_act, n0_act, n1_act, p_act / (n0_act + n1_act)
    ))
  }
  if (p_act > 500L || p_act > 2L * (n0_act + n1_act)) {
    warning(sprintf("BRIDGE: p is very large relative to n (p=%d, n=%d). ", p_act, n0_act + n1_act),
      "Consider prefilter='variance', max_genes=<200, or precision_method='mb' to avoid long runtimes. ",
      "glasso is safe up to ~250 genes; jgl up to ~500.",
      call. = FALSE
    )
  }

  mp(2, sprintf(
    "estimating networks: p=%d n=%d (slowest step)…",
    p_act, n0_act + n1_act
  ))
  est <- bridge_estimate_precision(X0, X1, Ba,
    sparsity_penalty = sparsity_penalty, beta = beta,
    fused_penalty = fused_penalty, data_type = data_type,
    precision_method = precision_method,
    progress = function(lbl) {
      mp(2, paste0("estimating ", lbl))
    }
  )
  profile_msg("Network Estimation", sprintf("Backend: %s", est$backend))
  # Ensure dimnames are the active features.
  dimnames(est$net0) <- list(active, active)
  dimnames(est$net1) <- list(active, active)
  if (verbose) message("BRIDGE: network backend = ", est$backend, " (", est$type, ")")

  # ── Differential network + multi-resolution modules ────────────────────────
  mp(3, "extracting modules (community detection)…")
  modules <- bridge_extract_modules(est$net0, est$net1, Ba, resolution,
    prior_weight = prior_weight,
    edge_tol = edge_tol,
    details = details,
    progress = function(lbl) mp(3, lbl)
  )
  n_candidates <- length(modules)
  profile_msg("Community Detection (Leiden)", sprintf("%d candidate modules", n_candidates))
  if (verbose) {
    message(
      "BRIDGE: ", n_candidates,
      " candidate modules across resolutions."
    )
  }

  # ── Score / stability-filter / select ──────────────────────────────────────
  mp(4, sprintf("scoring %d modules (bootstrap stability)…", n_candidates))
  selected <- bridge_score_modules(modules, Xa, y, Ba,
    cohesion_weight = cohesion_weight,
    size_penalty = size_penalty, tau = tau,
    n_boot = n_boot, max_overlap = max_overlap,
    details = details
  )
  profile_msg("Module Scoring & Bootstrap Stability", sprintf("%d modules passed initial stability", length(selected)))
  mp(5, "finalizing")
  diag_all <- attr(selected, "diagnostics") # NULL unless details

  # Size filter (matches gsm_dataset conventions).
  if (length(selected) > 0) {
    sizes <- vapply(selected, length, integer(1))
    keep <- sizes >= min_size & sizes <= max_size
    sc <- attr(selected, "scores")[keep]
    selected <- selected[keep]
    attr(selected, "scores") <- sc
    if (details && !is.null(diag_all)) {
      diag_all <- diag_all[keep, , drop = FALSE]
    }
  }

  if (length(selected) == 0) {
    warning("BRIDGE selected no modules passing the J(M)/stability/size ",
      "filters. Try lowering `tau`, `min_size`, or `lambda`.",
      call. = FALSE
    )
    out <- structure(list(), scores = numeric(0))
    if (details) {
      profile_times[["Total Runtime"]] <- as.numeric(difftime(Sys.time(), t_total_start, units = "secs"))
      attr(out, "summary") <- list(
        input_gene_count = n_features_total,
        retained_modules = 0L,
        candidate_modules = n_candidates,
        runtime_profiler = profile_times
      )
    }
    return(out)
  }

  if (verbose) {
    message(
      "BRIDGE: selected ", length(selected),
      " module(s) as final groups."
    )
  }

  # ── Optional rich diagnostics for the G-tab (default off; keeps MCCV fast) ──
  if (details) {
    raw_prior_genes <- if (is.null(Bmat)) 0L else as.integer(sum(Matrix::rowSums(Bmat != 0) > 0))
    prior_genes <- if (is.null(Ba)) {
      0L
    } else {
      as.integer(sum(Matrix::rowSums(methods::as(Ba, "CsparseMatrix") != 0) > 0))
    }
    prior_edges <- if (is.null(Ba)) {
      0L
    } else {
      as.integer(length(methods::as(Ba, "CsparseMatrix")@x) / 2)
    }

    edges_c0 <- sum(abs(est$net0) > edge_tol) / 2
    edges_c1 <- sum(abs(est$net1) > edge_tol) / 2
    edges_diff <- sum(abs(est$net1 - est$net0) > edge_tol) / 2

    attr(selected, "diagnostics") <- diag_all

    profile_times[["Total Runtime"]] <- as.numeric(difftime(Sys.time(), t_total_start, units = "secs"))

    attr(selected, "summary") <- list(
      pipeline_mode = pipeline_mode,
      input_gene_count = n_features_total,
      n_expression = length(feats),
      n_prior = if (!is.null(Bmat)) attr(Bmat, "n_prior") else 0L,
      n_intersection = if (!is.null(Bmat)) attr(Bmat, "n_intersection") else 0L,
      n_unmatched = if (!is.null(Bmat)) attr(Bmat, "n_unmatched") else 0L,
      raw_prior_gene_count = raw_prior_genes,
      prefilter_kept_count = prefilter_kept,
      active_gene_count = length(active),
      class_0_samples = n0_act,
      class_1_samples = n1_act,
      network_backend = est$backend,
      network_type = est$type,
      parameters = list(
        sparsity = sparsity_penalty, fused = fused_penalty,
        cohesion = cohesion_weight, size = size_penalty,
        beta = beta, prior_weight = prior_weight, edge_tol = edge_tol
      ),
      edges_class0 = edges_c0,
      edges_class1 = edges_c1,
      edges_differential = edges_diff,
      leiden_communities_found = sum(attr(modules, "res_counts")),
      candidate_modules = n_candidates,
      retained_modules = length(selected),
      resolutions = resolution,
      res_counts = attr(modules, "res_counts"),
      active_prior_genes = prior_genes,
      active_prior_edges = prior_edges,
      runtime_profiler = profile_times
    )
    attr(selected, "graph") <- bridge_module_edges(
      attr(modules, "graph_W"), unique(unlist(selected)),
      max_nodes = graph_max_nodes, max_edges = graph_max_edges
    )
  }

  if (verbose) {
    total_elapsed <- as.numeric(difftime(Sys.time(), t_total_start, units = "secs"))
    message(sprintf("BRIDGE Profile: [Total Runtime] %.2f sec", total_elapsed))
  }

  selected
}

# Edge list among selected-module genes from the blended graph W, capped for
# interactive display (top edges by weight, bounded node set).
#' @keywords internal
bridge_module_edges <- function(W, genes, max_nodes = 150L, max_edges = 2000L) {
  empty <- structure(
    data.frame(
      from = character(0), to = character(0), weight = numeric(0),
      stringsAsFactors = FALSE
    ),
    capped = FALSE
  )
  if (is.null(W) || length(genes) < 2) {
    return(empty)
  }
  genes <- intersect(genes, rownames(W))
  if (length(genes) < 2) {
    return(empty)
  }
  capped <- FALSE
  if (length(genes) > max_nodes) { # keep the highest-degree nodes
    deg <- rowSums(abs(W[genes, genes, drop = FALSE]))
    genes <- names(sort(deg, decreasing = TRUE))[seq_len(max_nodes)]
    capped <- TRUE
  }
  sub <- W[genes, genes, drop = FALSE]
  ut <- which(upper.tri(sub) & sub > 0, arr.ind = TRUE)
  if (nrow(ut) == 0) {
    return(empty)
  }
  df <- data.frame(
    from = genes[ut[, 1]], to = genes[ut[, 2]],
    weight = sub[ut], stringsAsFactors = FALSE
  )
  df <- df[order(df$weight, decreasing = TRUE), , drop = FALSE]
  if (nrow(df) > max_edges) {
    df <- df[seq_len(max_edges), , drop = FALSE]
    capped <- TRUE
  }
  structure(df, capped = capped)
}

# ── Internal glue ────────────────────────────────────────────────────────────

#' Standardize and resolve a user-supplied prior into a feature-aligned matrix
#'
#' This function acts as an explicit preprocessing and diagnostic step before
#' running BRIDGE. It aligns the prior graph `B` with the expression features,
#' reports the intersection/overlap percentage, and deduplicates the prior edges.
#'
#' @param B A biological prior (edge list dataframe, matrix, or named list of groups).
#' @param feats Character vector of expression feature names (e.g. `rownames(X)`).
#' @param verbose Logical; if `TRUE`, prints mapping statistics and overlap percentages.
#' @return A symmetric sparse matrix perfectly aligned to `feats` (with missing genes dropped),
#'   and attaching mapping statistics as attributes.
#' @export
bridge_standardize_prior <- function(B, feats, verbose = TRUE) {
  if (is.null(B)) {
    warning("No prior graph B supplied; BRIDGE runs data-driven (B = 0).",
      call. = FALSE
    )
    return(NULL)
  }

  # ── Validate ID mapping (Prior vs Expression Data) ─────────────────────────
  b_genes <- character(0)
  if (is.list(B) && !is.data.frame(B) && !methods::is(B, "Matrix")) {
    b_genes <- unique(as.character(unlist(B)))
  } else if ((is.matrix(B) || methods::is(B, "Matrix")) && nrow(B) == ncol(B) && !is.null(rownames(B))) {
    b_genes <- rownames(B)
  } else if (is.data.frame(B)) {
    b_genes <- unique(c(as.character(B[[1]]), as.character(B[[2]])))
  } else if (is.matrix(B) && ncol(B) >= 2) {
    b_genes <- unique(c(as.character(B[, 1]), as.character(B[, 2])))
  }

  if (length(b_genes) > 0) {
    overlap <- sum(b_genes %in% feats)
    match_pct <- (overlap / length(b_genes)) * 100
    if (verbose) {
      message(sprintf(
        "BRIDGE [prior]: Identifier mapping — %d prior genes supplied. %d matched with expression data (%.1f%%).",
        length(b_genes), overlap, match_pct
      ))
    }
    if (overlap == 0) {
      stop("BRIDGE [prior]: FATAL ERROR - 0% identifier overlap between biological prior and expression matrix. Check ID formats (e.g., Ensembl vs Gene Symbol mismatch).", call. = FALSE)
    }
    if (match_pct < 5.0) {
      warning(sprintf("BRIDGE [prior]: DANGER - Only %.1f%% of biological prior genes matched the expression data. Are ID types correct? This will cause the algorithm to largely ignore the biological knowledge and act as a purely data-driven method.", match_pct), call. = FALSE)
    }
  }

  # Named grouping list -> co-membership prior.
  ret <- NULL
  if (is.list(B) && !is.data.frame(B) && !methods::is(B, "Matrix")) {
    ret <- bridge_prior_from_groups(B, feats)
  } else if ((is.matrix(B) || methods::is(B, "Matrix")) && nrow(B) == ncol(B) && !is.null(rownames(B))) {
    if (length(intersect(rownames(B), feats)) == 0) {
      stop("Prior matrix and X share no feature names. ",
        "Check that B's dimnames match rownames(X).",
        call. = FALSE
      )
    }
    ret <- bridge_prior_from_edgelist(B, feats)
  } else {
    ret <- bridge_prior_from_edgelist(B, feats)
  }

  if (!is.null(ret) && length(b_genes) > 0) {
    attr(ret, "n_prior") <- length(b_genes)
    attr(ret, "n_intersection") <- overlap
    attr(ret, "n_unmatched") <- length(b_genes) - overlap
  }

  ret
}

#' Restrict to features in non-trivial connected components of B (p >> n aid)
#' @keywords internal
bridge_active_features <- function(X, B, verbose = TRUE) {
  feats <- rownames(X)
  p <- length(feats)
  n <- ncol(X)
  # Only prune when p >> n and B is informative enough to define components.
  if (is.null(B) || p <= 2 * n) {
    return(feats)
  }
  if (!requireNamespace("igraph", quietly = TRUE)) {
    return(feats)
  }

  g <- igraph::graph_from_adjacency_matrix(B,
    mode = "undirected",
    weighted = TRUE, diag = FALSE
  )
  comp <- igraph::components(g)
  big <- which(comp$csize >= 2)
  active <- feats[comp$membership %in% big]
  if (length(active) < 2) {
    return(feats)
  } # nothing useful to restrict to

  if (verbose) {
    n_comp <- length(big)
    if (n_comp > 0) {
      max_comp <- max(comp$csize[big])
      message(sprintf(
        "BRIDGE [prior]: Found %d informative components (size >= 2). Largest component has %d genes.",
        n_comp, max_comp
      ))

      if (max_comp > 500L) {
        message(
          "BRIDGE [prior]: WARNING — The prior graph contains a giant component (", max_comp, " genes). ",
          "Connected-component filtering alone will not reduce 'p' enough to make estimation fast. ",
          "Consider using prefilter='variance' with max_genes to strictly cap estimation size."
        )
      }
    }

    if (length(active) < p) {
      message(
        "BRIDGE [prior]: p >> n; restricted estimation from ", p, " to ", length(active),
        " features in connected components of B."
      )
    } else {
      message("BRIDGE [prior]: Connected-component filtering did not reduce the feature set (p=", p, ").")
    }
  }

  active
}
