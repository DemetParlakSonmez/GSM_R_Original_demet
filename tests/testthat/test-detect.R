# Tests for G-S-M detect (de novo robust gene-module discovery)

# Build a synthetic folds list: two co-occurring gene sets with stable, opposite
# directions, plus one planted noise gene that flips around.
make_folds <- function(n = 12, seed = 1) {
  set.seed(seed)
  setA <- c("A1", "A2", "A3", "A4")          # stable, UP (strong positive t)
  setB <- c("B1", "B2", "B3")                # stable, DOWN (strong negative t)
  noise <- paste0("N", 1:10)                 # many weak genes (fill the bottom)
  lapply(seq_len(n), function(i) {
    a  <- stats::setNames(runif(length(setA), 6, 8), setA)       # strong +
    b  <- stats::setNames(-runif(length(setB), 6, 8), setB)      # strong -
    nz <- stats::setNames(rnorm(length(noise), 0, 0.4), noise)   # weak, random sign
    list(status = "success", gene_stats = c(a, b, nz),
         importance = c(a, abs(b), abs(nz)))
  })
}

fake_dataset <- function() {
  structure(list(
    groups = list(Known_A = c("A1", "A2", "A3", "A4", "A5"),  # superset of setA
                  Known_B = c("B1", "B2", "B3")),
    params = list(n_samples = 20L, n_genes = 8L, n_groups = 2L)),
    class = "gsm_dataset")
}

test_that("rank normalization gives [0,1] ranks, presence and signs", {
  nz <- gsmR:::detect_normalize_ranks(make_folds(), mode = "abs")
  expect_true(all(nz$rank_mat >= 0 & nz$rank_mat <= 1))
  expect_equal(dim(nz$presence), dim(nz$rank_mat))
  expect_true(all(nz$signs %in% c(-1, 0, 1)))
  # UP mode keeps only positive-t genes in the ranked lists.
  nz_up <- gsmR:::detect_normalize_ranks(make_folds(), mode = "up")
  expect_true(all(c("A1", "A2") %in% nz_up$universe))
  expect_false(any(c("B1", "B2") %in% nz_up$universe))
})

test_that("robust genes recover the stable sets with high direction consistency", {
  nz <- gsmR:::detect_normalize_ranks(make_folds(), mode = "abs")
  rg <- gsmR:::detect_robust_genes(nz, top_frac = 0.5, min_freq = 0.5)
  expect_true(all(c("A1", "B1") %in% rg$Gene))
  expect_true(all(rg$Direction_Consistency[rg$Gene %in% c("A1","A2","B1","B2")] == 1))
  expect_false("N1" %in% rg$Gene)              # noise gene filtered out
})

test_that("community detection recovers >= 2 modules (hclust fallback)", {
  nz <- gsmR:::detect_normalize_ranks(make_folds(), mode = "abs")
  rg <- gsmR:::detect_robust_genes(nz, top_frac = 0.5, min_freq = 0.5)
  adj <- gsmR:::detect_stability_graph(rg$Gene, nz, top_k = 0.5)
  memb <- gsmR:::detect_communities(adj, method = "hclust")
  expect_gte(length(unique(memb)), 2)
})

test_that("annotation maps a module to its known group", {
  ann <- gsmR:::detect_annotate(list(c("B1", "B2", "B3")), fake_dataset()$groups)
  expect_equal(ann$Best_Known_Group, "Known_B")
  expect_equal(ann$Jaccard, 1)
  expect_equal(ann$Label, "rediscovered")
})

test_that("gsm_detect end-to-end returns the documented shape", {
  d <- gsm_detect(list(folds = make_folds()), fake_dataset(),
                  mode = "abs", top_genes = 50, top_k_frac = 0.5,
                  min_freq = 0.5, n_boot = 5)
  expect_true(all(c("robust_genes", "modules", "members", "graph", "params") %in% names(d)))
  expect_true(all(c("Module", "Size", "Best_Known_Group", "Jaccard", "Label",
                    "Recovery_Freq") %in% names(d$modules)))
  expect_gte(nrow(d$modules), 1)
})
