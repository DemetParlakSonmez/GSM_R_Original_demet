test_that("gsm_gene_adjacency builds a normalized graph from a BRIDGE edge list", {
  genes <- c("G1", "G2", "G3", "G4")
  el <- data.frame(from = c("G1", "G2"), to = c("G2", "G3"),
                   weight = c(0.9, 0.5), stringsAsFactors = FALSE)
  A <- gsm_gene_adjacency(genes, graph = el)

  expect_equal(dim(A), c(4L, 4L))
  expect_equal(dimnames(A), list(genes, genes))
  expect_equal(A, t(A))                    # symmetric
  expect_true(all(is.finite(A)))
  expect_true(all(diag(A) > 0))            # self-loops survive normalization
  # G4 is isolated in the edge list: self-loop only.
  expect_equal(unname(A["G4", setdiff(genes, "G4")]), rep(0, 3))
})

test_that("gsm_gene_adjacency falls back to a TRAIN-ONLY correlation graph", {
  set.seed(1)
  X <- matrix(rnorm(6 * 30), nrow = 6,
              dimnames = list(paste0("G", 1:6), paste0("S", 1:30)))
  genes <- rownames(X)
  train <- 1:20

  A_tr  <- gsm_gene_adjacency(genes, graph = NULL, expr = X, train_idx = train)
  A_all <- gsm_gene_adjacency(genes, graph = NULL, expr = X,
                              train_idx = seq_len(ncol(X)))

  expect_equal(dim(A_tr), c(6L, 6L))
  expect_equal(A_tr, t(A_tr))
  # Held-out samples must not influence the graph, so the two differ.
  expect_false(isTRUE(all.equal(A_tr, A_all)))
})

test_that("a supplied BRIDGE graph takes priority over the correlation fallback", {
  set.seed(2)
  X <- matrix(rnorm(4 * 20), nrow = 4,
              dimnames = list(paste0("G", 1:4), paste0("S", 1:20)))
  genes <- rownames(X)
  el <- data.frame(from = "G1", to = "G2", weight = 1, stringsAsFactors = FALSE)

  A_bridge <- gsm_gene_adjacency(genes, graph = el, expr = X, train_idx = 1:14)
  A_cor    <- gsm_gene_adjacency(genes, graph = NULL, expr = X, train_idx = 1:14)

  expect_false(isTRUE(all.equal(A_bridge, A_cor)))
  # Only the G1-G2 edge exists, so G3/G4 stay isolated.
  expect_equal(unname(A_bridge["G3", "G4"]), 0)
})

test_that("score_group and build_model dispatch to mlp and gnn", {
  skip_if_not_installed("torch")
  skip_if_not(torch::torch_is_installed(), "LibTorch backend not installed")

  set.seed(3)
  p <- 20; n <- 40
  X <- matrix(rnorm(p * n), nrow = p,
              dimnames = list(paste0("G", 1:p), paste0("S", 1:n)))
  y <- factor(rep(c("neg", "pos"), each = n / 2))
  names(y) <- colnames(X)
  genes <- paste0("G", 1:8)
  X[genes, y == "pos"] <- X[genes, y == "pos"] + 2   # a real signal to find
  train <- c(1:14, 21:34)

  for (m in c("mlp", "gnn")) {
    r <- score_group(X, y, genes, train, method = m, epochs = 30)
    expect_true(is.finite(r$accuracy))
    expect_true(r$accuracy >= 0 && r$accuracy <= 1)
    expect_length(r$predictions, n - length(train))
    expect_true(all(levels(r$predictions) == levels(y)))

    f1 <- score_group_train(X, y, genes, train, method = m, cv = 3, epochs = 30)
    expect_true(is.finite(f1))

    mo <- build_model(X, y, list(real = genes), train, method = m, epochs = 30)
    expect_true(is.finite(mo$metrics$F_measure))
    expect_equal(mo$features, genes)
    # Both backends must report per-gene importance so they feed gsm_biomarkers().
    expect_false(is.null(mo$importance))
    expect_named(mo$importance, genes)
    expect_true(all(is.finite(mo$importance)))
  }
})

test_that("the GNN accepts BRIDGE's differential network as its adjacency", {
  skip_if_not_installed("torch")
  skip_if_not(torch::torch_is_installed(), "LibTorch backend not installed")

  set.seed(4)
  p <- 20; n <- 40
  X <- matrix(rnorm(p * n), nrow = p,
              dimnames = list(paste0("G", 1:p), paste0("S", 1:n)))
  y <- factor(rep(c("neg", "pos"), each = n / 2))
  names(y) <- colnames(X)
  genes <- paste0("G", 1:8)
  X[genes, y == "pos"] <- X[genes, y == "pos"] + 2
  train <- c(1:14, 21:34)

  # An explicit gene-gene graph, as BRIDGE would hand over.
  el <- data.frame(from = genes[-length(genes)], to = genes[-1],
                   weight = 1, stringsAsFactors = FALSE)
  r <- build_model(X, y, list(real = genes), train, method = "gnn",
                   graph = el, epochs = 30)
  expect_true(is.finite(r$metrics$F_measure))
  expect_named(r$importance, genes)
})

test_that("grouping_for_fold attaches BRIDGE's graph only when asked", {
  skip_if_not_installed("glasso")
  skip_if_not_installed("igraph")

  set.seed(5)
  p <- 30; n <- 40
  X <- matrix(rnorm(p * n), nrow = p,
              dimnames = list(paste0("G", 1:p), paste0("S", 1:n)))
  y <- factor(rep(c("neg", "pos"), each = n / 2))
  names(y) <- colnames(X)
  mod <- paste0("G", 1:10)
  X[mod, y == "pos"] <- X[mod, y == "pos"] + 1.5
  groups <- list(real = mod, other = paste0("G", 11:20))

  ds <- gsm_dataset(X, y, groups)
  ds$grouping_spec <- list(
    method = "bridge", prior = groups,
    args = list(resolution = c(1, 1.5), tau = 0.3, min_size = 3, verbose = FALSE))
  train <- c(1:14, 21:34)

  g_off <- grouping_for_fold(ds, train, expr = X, fold_seed = 1)
  expect_null(attr(g_off, "graph"))         # default stays fast

  g_on <- grouping_for_fold(ds, train, expr = X, fold_seed = 1, with_graph = TRUE)
  gr <- attr(g_on, "graph")
  skip_if(is.null(gr) || nrow(gr) == 0, "BRIDGE selected no modules on this seed")
  expect_true(all(c("from", "to", "weight") %in% names(gr)))
  # The display cap (150 nodes) must not truncate the graph the GNN convolves over.
  expect_false(isTRUE(attr(gr, "capped")))
})

test_that("unknown classifier methods are rejected with the full backend list", {
  X <- matrix(rnorm(20), nrow = 4,
              dimnames = list(paste0("G", 1:4), paste0("S", 1:5)))
  y <- factor(c("a", "a", "b", "b", "b"))
  names(y) <- colnames(X)
  expect_error(score_group(X, y, rownames(X), 1:3, method = "nope"), "gnn")
  expect_error(build_model(X, y, list(g = rownames(X)), 1:3, method = "nope"), "gnn")
})
