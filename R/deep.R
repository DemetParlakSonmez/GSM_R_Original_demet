#' @title Deep-learning classifiers (MLP and Graph Neural Network)
#' @description Two `torch`-backed classifiers that plug into the Scoring (S)
#'   and Modeling (M) steps alongside the existing `rf`/`svm`/`knn`/`xgb`
#'   backends:
#'
#'   * `"mlp"` — a multi-layer perceptron over the group's genes.
#'   * `"gnn"` — a graph convolutional network (GCN). Each *sample* is a graph
#'     whose nodes are the group's genes, whose node feature is that gene's
#'     expression in that sample, and whose adjacency is a gene-gene graph
#'     shared across samples. The adjacency is BRIDGE's differential network
#'     when a BRIDGE grouping is active (see [bridge_grouping()]), which couples
#'     the Grouping and Modeling steps; otherwise it falls back to a
#'     training-only correlation graph so the GNN also works with file-based
#'     groupings (KEGG, GO, DisGeNET).
#'
#'   Both are two-class only, and both report per-gene importance (input-gradient
#'   saliency) so they feed the consensus biomarker ranking in [gsm_biomarkers()]
#'   exactly like the Random Forest and XGBoost backends do.
#'
#'   Note these models have far more capacity than the sample sizes typical of
#'   transcriptomic studies (tens to a few hundred samples). Defaults are
#'   deliberately small and heavily regularized (dropout + weight decay), but
#'   expect them to be competitive with, not automatically better than, `rf` and
#'   `xgb`.
#' @name deep
NULL

# ── torch backend plumbing ───────────────────────────────────────────────────

# Fail with actionable guidance rather than a bare "there is no package" error.
# Also pins torch to a single thread: MCCV folds run inside forked workers
# (mclapply), and a multi-threaded LibTorch inherited across fork() can deadlock
# the child — the same failure mode this package already guards against for
# OpenBLAS and sets `nthread = 1` for in XGBoost.
.torch_require <- function(method) {
  if (!requireNamespace("torch", quietly = TRUE)) {
    stop("Package 'torch' is required for method '", method, "'. Install it with:\n",
         "  install.packages(\"torch\"); torch::install_torch()", call. = FALSE)
  }
  if (!torch::torch_is_installed()) {
    stop("The 'torch' R package is present but its LibTorch backend is not ",
         "installed. Run: torch::install_torch()", call. = FALSE)
  }
  try(torch::torch_set_num_threads(1L), silent = TRUE)
  invisible(TRUE)
}

# Resolve the binary problem: levels, positive class, 0/1 label vector, and the
# positive-class weight that offsets class imbalance in the BCE loss.
.binary_setup <- function(y_tr, positive_class = NULL) {
  y_tr <- droplevels(as.factor(y_tr))
  lev <- levels(y_tr)
  if (length(lev) != 2)
    stop("Deep-learning backends support two-class problems only; found ",
         length(lev), " classes.", call. = FALSE)
  if (is.null(positive_class) || !positive_class %in% lev)
    positive_class <- if ("pos" %in% lev) "pos" else lev[1]
  neg <- setdiff(lev, positive_class)
  label <- as.numeric(y_tr == positive_class)
  n_pos <- sum(label); n_neg <- length(label) - n_pos
  pos_weight <- if (n_pos > 0) n_neg / n_pos else 1
  list(levels = lev, pos = positive_class, neg = neg, label = label,
       pos_weight = pos_weight)
}

# Standardize columns using TRAIN-only mean/sd, then apply to both matrices.
# Zero-variance genes would otherwise produce NaN and poison the whole network.
.standardize_train <- function(x_tr, x_te) {
  mu <- colMeans(x_tr)
  sdv <- apply(x_tr, 2, stats::sd)
  sdv[!is.finite(sdv) | sdv == 0] <- 1
  list(tr = scale(x_tr, center = mu, scale = sdv),
       te = scale(x_te, center = mu, scale = sdv))
}

# Full-batch Adam training loop with BCE-on-logits. Sample counts here are small
# (tens to hundreds), so full batch is both faster and more stable than
# mini-batching, and it keeps the loop free of a dataloader dependency.
.torch_train_loop <- function(model, forward, y, epochs, lr, weight_decay,
                              pos_weight) {
  opt <- torch::optim_adam(model$parameters, lr = lr, weight_decay = weight_decay)
  crit <- torch::nn_bce_with_logits_loss(
    pos_weight = torch::torch_tensor(pos_weight, dtype = torch::torch_float()))
  model$train()
  for (e in seq_len(epochs)) {
    opt$zero_grad()
    loss <- crit(forward(), y)
    loss$backward()
    opt$step()
  }
  model$eval()
  invisible(model)
}

# Per-feature saliency: mean |d logit / d input| over the training samples.
# For the MLP the features are genes directly; for the GNN they are graph nodes,
# which are also genes. Either way this yields a named per-gene importance, so
# `mlp` and `gnn` feed gsm_biomarkers() like `rf` (Gini) and `xgb` (Gain) do.
.saliency <- function(model, x, forward_on, feature_names) {
  tryCatch({
    xg <- x$clone()$detach()$requires_grad_(TRUE)
    out <- forward_on(xg)$sum()
    out$backward()
    g <- xg$grad$abs()
    # Reduce to one value per feature. MLP input is [samples, genes]; GNN input is
    # [samples, genes, node_feats]. Average over samples (the leading dim) first,
    # then collapse any TRAILING node-feature dims — averaging dim 1 twice would
    # collapse the genes themselves and leave nothing to name.
    g <- g$mean(dim = 1)
    while (length(dim(g)) > 1) g <- g$mean(dim = length(dim(g)))
    stats::setNames(as.numeric(as.array(g)), feature_names)
  }, error = function(e) NULL)
}

# ── Gene-gene adjacency for the GNN ──────────────────────────────────────────

#' Build the gene-gene adjacency a GNN convolves over
#'
#' Resolves the graph the `"gnn"` backend needs, in priority order:
#'
#' 1. `graph` — an edge list (`from`, `to`, `weight`), as produced by BRIDGE and
#'    attached to its groups (the `"graph"` attribute of [bridge_grouping()]
#'    output). This is the differential network `|Theta1 - Theta0|` blended with
#'    the prior, so the GNN convolves over exactly the interaction structure that
#'    defined the groups.
#' 2. A **training-only** absolute-correlation graph over `genes`, thresholded at
#'    `cor_threshold`. Used when no BRIDGE graph is available (e.g. a KEGG/GO
#'    file grouping), so the GNN is not restricted to BRIDGE runs. Computed on
#'    `train_idx` only — never the test samples — so it introduces no leakage.
#'
#' The result is symmetrically normalized with self-loops,
#' `A_hat = D^-1/2 (A + I) D^-1/2`, the standard GCN propagation matrix.
#'
#' @param genes Character vector of genes (the graph's nodes), in the column
#'   order the feature matrix uses.
#' @param graph Optional edge-list `data.frame` (`from`, `to`, `weight`).
#' @param expr Numeric matrix (genes x samples), used for the correlation
#'   fallback.
#' @param train_idx Integer vector of training-sample columns; the correlation
#'   fallback sees only these.
#' @param cor_threshold Absolute-correlation cutoff for the fallback graph.
#' @return A dense symmetric numeric matrix (`length(genes)` square), normalized
#'   for GCN propagation, with `genes` as dimnames.
#' @export
gsm_gene_adjacency <- function(genes, graph = NULL, expr = NULL,
                               train_idx = NULL, cor_threshold = 0.3) {
  g <- length(genes)
  A <- matrix(0, g, g, dimnames = list(genes, genes))

  # 1. BRIDGE differential network, restricted to this group's genes.
  if (!is.null(graph) && is.data.frame(graph) && nrow(graph) > 0 &&
      all(c("from", "to") %in% names(graph))) {
    e <- graph[graph$from %in% genes & graph$to %in% genes, , drop = FALSE]
    if (nrow(e) > 0) {
      w <- if ("weight" %in% names(e)) as.numeric(e$weight) else rep(1, nrow(e))
      w[!is.finite(w)] <- 0
      A[cbind(match(e$from, genes), match(e$to, genes))] <- w
      A[cbind(match(e$to, genes), match(e$from, genes))] <- w
    }
  }

  # 2. Correlation fallback (train-only) when the BRIDGE graph gave us nothing.
  if (sum(A) == 0 && !is.null(expr) && g >= 2) {
    cols <- if (is.null(train_idx)) seq_len(ncol(expr)) else train_idx
    sub <- expr[intersect(genes, rownames(expr)), cols, drop = FALSE]
    if (nrow(sub) >= 2 && ncol(sub) >= 3) {
      C <- suppressWarnings(abs(stats::cor(t(sub))))
      C[!is.finite(C)] <- 0
      diag(C) <- 0
      C[C < cor_threshold] <- 0
      keep <- rownames(sub)
      A[keep, keep] <- C
    }
  }

  diag(A) <- 0
  # GCN normalization with self-loops: D^-1/2 (A + I) D^-1/2.
  A_sl <- A + diag(1, g)
  d <- rowSums(A_sl)
  d[d <= 0] <- 1
  dinv <- 1 / sqrt(d)
  A_hat <- A_sl * outer(dinv, dinv)
  dimnames(A_hat) <- list(genes, genes)
  A_hat
}

# ── MLP ──────────────────────────────────────────────────────────────────────

# Fit an MLP and predict positive-class probabilities on x_te. Shared by
# score_group_mlp(), .cv_predict() and build_model_mlp(), mirroring .xgb_fit().
# Returns list(prob, pos, neg, levels, importance).
.mlp_fit <- function(x_tr, y_tr, x_te, positive_class = NULL,
                     hidden = c(64L, 32L), dropout = 0.3, epochs = 100L,
                     lr = 0.01, weight_decay = 1e-4, seed = NULL, ...) {
  .torch_require("mlp")
  if (!is.null(seed)) torch::torch_manual_seed(as.integer(seed))
  bs <- .binary_setup(y_tr, positive_class)

  sc <- .standardize_train(as.matrix(x_tr), as.matrix(x_te))
  p <- ncol(sc$tr)
  xt <- torch::torch_tensor(sc$tr, dtype = torch::torch_float())
  xv <- torch::torch_tensor(sc$te, dtype = torch::torch_float())
  yt <- torch::torch_tensor(bs$label, dtype = torch::torch_float())

  # Sequential MLP: [linear -> relu -> dropout] x len(hidden), then a logit head.
  layers <- list(); prev <- p
  for (h in as.integer(hidden)) {
    layers <- c(layers, list(torch::nn_linear(prev, h), torch::nn_relu(),
                             torch::nn_dropout(dropout)))
    prev <- h
  }
  layers <- c(layers, list(torch::nn_linear(prev, 1L)))
  model <- do.call(torch::nn_sequential, layers)

  fwd <- function(x) model(x)$squeeze(2)
  .torch_train_loop(model, function() fwd(xt), yt, epochs, lr, weight_decay,
                    bs$pos_weight)

  prob <- as.numeric(as.array(torch::torch_sigmoid(fwd(xv))))
  imp <- .saliency(model, xt, fwd, colnames(sc$tr))
  list(prob = prob, pos = bs$pos, neg = bs$neg, levels = bs$levels,
       importance = imp)
}

# ── GNN (graph convolutional network) ────────────────────────────────────────

# Two-layer GCN over a fixed gene-gene graph. Each sample is a graph: node
# features are that sample's expression values (one scalar per gene), the
# adjacency A_hat is shared across samples. Propagation is H' = relu(A_hat H W);
# after two layers a global mean-pool over nodes gives a graph-level embedding,
# and a linear head produces the logit.
.gcn_module <- function() {
  torch::nn_module(
    "gsm_gcn",
    initialize = function(n_feat, hidden, dropout) {
      self$l1 <- torch::nn_linear(n_feat, hidden)
      self$l2 <- torch::nn_linear(hidden, hidden)
      self$drop <- torch::nn_dropout(dropout)
      self$head <- torch::nn_linear(hidden, 1L)
    },
    forward = function(x, A) {
      # x: [batch, nodes, n_feat]; A: [nodes, nodes] (broadcasts over batch).
      h <- torch::nnf_relu(torch::torch_matmul(A, self$l1(x)))
      h <- self$drop(h)
      h <- torch::nnf_relu(torch::torch_matmul(A, self$l2(h)))
      h <- h$mean(dim = 2)          # global mean pool over nodes -> [batch, hidden]
      self$head(h)                  # [batch, 1]
    }
  )
}

# Fit a GCN and predict positive-class probabilities on x_te.
# x_tr/x_te are samples x genes; `adj` is the normalized gene-gene adjacency
# from gsm_gene_adjacency(), aligned to the columns of x_tr.
.gnn_fit <- function(x_tr, y_tr, x_te, adj, positive_class = NULL,
                     hidden = 32L, dropout = 0.3, epochs = 150L,
                     lr = 0.01, weight_decay = 1e-4, seed = NULL, ...) {
  .torch_require("gnn")
  if (!is.null(seed)) torch::torch_manual_seed(as.integer(seed))
  bs <- .binary_setup(y_tr, positive_class)

  sc <- .standardize_train(as.matrix(x_tr), as.matrix(x_te))
  genes <- colnames(sc$tr)

  # Node features: one scalar (the standardized expression) per gene per sample.
  # [samples, genes] -> [samples, genes, 1]
  to_nodes <- function(m) torch::torch_tensor(
    array(as.numeric(m), dim = c(nrow(m), ncol(m), 1L)),
    dtype = torch::torch_float())
  xt <- to_nodes(sc$tr)
  xv <- to_nodes(sc$te)
  yt <- torch::torch_tensor(bs$label, dtype = torch::torch_float())
  A <- torch::torch_tensor(adj, dtype = torch::torch_float())

  model <- .gcn_module()(n_feat = 1L, hidden = as.integer(hidden),
                         dropout = dropout)

  fwd <- function(x) model(x, A)$squeeze(2)
  .torch_train_loop(model, function() fwd(xt), yt, epochs, lr, weight_decay,
                    bs$pos_weight)

  prob <- as.numeric(as.array(torch::torch_sigmoid(fwd(xv))))
  imp <- .saliency(model, xt, fwd, genes)
  list(prob = prob, pos = bs$pos, neg = bs$neg, levels = bs$levels,
       importance = imp)
}

# ── Scoring-step entry points ────────────────────────────────────────────────

#' Score a group using a multi-layer perceptron
#'
#' @param expr Numeric matrix (genes x samples).
#' @param labels Factor of class labels.
#' @param genes Character vector of gene names.
#' @param train_idx Integer vector of training sample indices.
#' @param ... Hyperparameters forwarded to the MLP: `hidden` (integer vector of
#'   layer widths), `dropout`, `epochs`, `lr`, `weight_decay`, `seed`.
#' @return A list with `accuracy`, `auc`, `predictions`, `prob`, and `actual`.
#' @keywords internal
score_group_mlp <- function(expr, labels, genes, train_idx, ...) {
  test_idx <- setdiff(seq_len(ncol(expr)), train_idx)
  if (length(test_idx) == 0) stop("No test samples.")

  x_tr <- t(expr[genes, train_idx, drop = FALSE])
  x_te <- t(expr[genes, test_idx, drop = FALSE])
  test_labels <- labels[test_idx]
  positive_class <- if ("pos" %in% levels(labels)) "pos" else levels(labels)[1]

  res <- .mlp_fit(x_tr, labels[train_idx], x_te, positive_class, ...)
  predictions <- factor(ifelse(res$prob > 0.5, res$pos, res$neg),
                        levels = levels(labels))
  list(
    accuracy    = calc_accuracy(test_labels, predictions),
    auc         = calc_auc(test_labels, res$prob, positive_class),
    predictions = predictions,
    prob        = res$prob,
    actual      = test_labels
  )
}

#' Score a group using a graph convolutional network
#'
#' Convolves over the gene-gene graph induced on this group's genes — BRIDGE's
#' differential network when available, otherwise a training-only correlation
#' graph. See [gsm_gene_adjacency()].
#'
#' @param expr Numeric matrix (genes x samples).
#' @param labels Factor of class labels.
#' @param genes Character vector of gene names.
#' @param train_idx Integer vector of training sample indices.
#' @param graph Optional edge list (`from`, `to`, `weight`) — BRIDGE's
#'   differential network for this fold. `NULL` falls back to a train-only
#'   correlation graph.
#' @param cor_threshold Absolute-correlation cutoff for the fallback graph.
#' @param ... Hyperparameters forwarded to the GCN: `hidden`, `dropout`,
#'   `epochs`, `lr`, `weight_decay`, `seed`.
#' @return A list with `accuracy`, `auc`, `predictions`, `prob`, and `actual`.
#' @keywords internal
score_group_gnn <- function(expr, labels, genes, train_idx, graph = NULL,
                            cor_threshold = 0.3, ...) {
  test_idx <- setdiff(seq_len(ncol(expr)), train_idx)
  if (length(test_idx) == 0) stop("No test samples.")

  x_tr <- t(expr[genes, train_idx, drop = FALSE])
  x_te <- t(expr[genes, test_idx, drop = FALSE])
  test_labels <- labels[test_idx]
  positive_class <- if ("pos" %in% levels(labels)) "pos" else levels(labels)[1]

  adj <- gsm_gene_adjacency(genes, graph = graph, expr = expr,
                            train_idx = train_idx,
                            cor_threshold = cor_threshold)
  res <- .gnn_fit(x_tr, labels[train_idx], x_te, adj, positive_class, ...)
  predictions <- factor(ifelse(res$prob > 0.5, res$pos, res$neg),
                        levels = levels(labels))
  list(
    accuracy    = calc_accuracy(test_labels, predictions),
    auc         = calc_auc(test_labels, res$prob, positive_class),
    predictions = predictions,
    prob        = res$prob,
    actual      = test_labels
  )
}

# ── Modeling-step entry points ───────────────────────────────────────────────

#' Build an MLP model from top groups
#'
#' @param expr Numeric matrix (genes x samples).
#' @param labels Factor of class labels.
#' @param top_groups Named list of character vectors.
#' @param train_idx Integer vector of training sample indices.
#' @param ... Hyperparameters forwarded to the MLP.
#' @return A list with `metrics`, `predictions`, `actual`, `features`,
#'   `n_groups`, and `importance` (per-gene input-gradient saliency).
#' @keywords internal
build_model_mlp <- function(expr, labels, top_groups, train_idx, ...) {
  prep <- .deep_model_prep(expr, labels, top_groups, train_idx)
  if (is.null(prep$genes)) return(prep$empty)

  res <- .mlp_fit(prep$x_tr, labels[train_idx], prep$x_te,
                  prep$positive_class, ...)
  .deep_model_out(prep, res, labels, top_groups)
}

#' Build a graph convolutional network model from top groups
#'
#' Nodes are the genes of the top-k groups; the adjacency is BRIDGE's
#' differential network when a BRIDGE grouping is active, otherwise a
#' training-only correlation graph. See [gsm_gene_adjacency()].
#'
#' @param expr Numeric matrix (genes x samples).
#' @param labels Factor of class labels.
#' @param top_groups Named list of character vectors.
#' @param train_idx Integer vector of training sample indices.
#' @param graph Optional edge list (`from`, `to`, `weight`) — BRIDGE's
#'   differential network for this fold.
#' @param cor_threshold Absolute-correlation cutoff for the fallback graph.
#' @param ... Hyperparameters forwarded to the GCN.
#' @return A list with `metrics`, `predictions`, `actual`, `features`,
#'   `n_groups`, and `importance` (per-gene node saliency).
#' @keywords internal
build_model_gnn <- function(expr, labels, top_groups, train_idx, graph = NULL,
                            cor_threshold = 0.3, ...) {
  prep <- .deep_model_prep(expr, labels, top_groups, train_idx)
  if (is.null(prep$genes)) return(prep$empty)

  adj <- gsm_gene_adjacency(prep$genes, graph = graph, expr = expr,
                            train_idx = train_idx,
                            cor_threshold = cor_threshold)
  res <- .gnn_fit(prep$x_tr, labels[train_idx], prep$x_te, adj,
                  prep$positive_class, ...)
  .deep_model_out(prep, res, labels, top_groups)
}

# Shared front half of the two build_model_* functions: resolve the feature set,
# split train/test, and pre-build the empty result for the no-usable-genes case.
.deep_model_prep <- function(expr, labels, top_groups, train_idx) {
  test_idx <- setdiff(seq_len(ncol(expr)), train_idx)
  if (length(test_idx) == 0) stop("No test samples.")
  all_genes <- intersect(unique(unlist(top_groups)), rownames(expr))

  if (length(all_genes) == 0) {
    return(list(genes = NULL, empty = list(
      metrics     = calc_all_metrics(labels[test_idx], labels[test_idx]),
      predictions = factor(rep(NA, length(test_idx)), levels = levels(labels)),
      actual      = labels[test_idx],
      features    = character(0),
      n_groups    = length(top_groups),
      importance  = NULL
    )))
  }
  list(
    genes = all_genes, test_idx = test_idx,
    x_tr = t(expr[all_genes, train_idx, drop = FALSE]),
    x_te = t(expr[all_genes, test_idx, drop = FALSE]),
    test_labels = labels[test_idx],
    positive_class = if ("pos" %in% levels(labels)) "pos" else levels(labels)[1]
  )
}

# Shared back half: threshold probabilities, score metrics, package the result.
.deep_model_out <- function(prep, res, labels, top_groups) {
  predictions <- factor(ifelse(res$prob > 0.5, res$pos, res$neg),
                        levels = levels(labels))
  list(
    metrics     = calc_all_metrics(prep$test_labels, predictions,
                                   prep$positive_class, prob = res$prob),
    predictions = predictions,
    actual      = prep$test_labels,
    features    = prep$genes,
    n_groups    = length(top_groups),
    importance  = res$importance
  )
}
