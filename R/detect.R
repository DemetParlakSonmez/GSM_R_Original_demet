#' @title G-S-M detect: de novo Robust Rank Group Discovery
#' @description Discover gene modules de novo from the per-iteration ranked gene
#'   lists produced by an MCCV run, then annotate them against the dataset's
#'   known groups. Where G-S-M uses predefined groupings, this spin-off infers
#'   stable gene modules directly from repeated stochastic rankings.
#' @name detect
NULL

# Successful folds with a usable per-gene ranking signal.
.detect_ok_folds <- function(folds) {
  Filter(function(r) !is.null(r) && identical(r$status, "success") &&
           (!is.null(r$gene_stats) || !is.null(r$importance)), folds)
}

# Best-first ranked gene vector for one fold under a direction mode.
.fold_ranked <- function(f, mode) {
  s <- f$gene_stats
  if (is.null(s) || length(s) == 0) {                 # fall back to RF importance (abs)
    imp <- f$importance
    if (is.null(imp) || length(imp) == 0) return(NULL)
    return(names(sort(imp, decreasing = TRUE)))
  }
  if (mode == "up")        s <- s[s > 0]
  else if (mode == "down") s <- -s[s < 0]             # flip so larger = stronger down
  else                     s <- abs(s)                # abs
  if (length(s) == 0) return(NULL)
  names(sort(s, decreasing = TRUE))
}

# ── Step 1 — normalize ranks per run ─────────────────────────────────────────
#' @keywords internal
detect_normalize_ranks <- function(folds, mode = "abs", gene_universe = NULL,
                                    missing_penalty = 1) {
  ok <- .detect_ok_folds(folds)
  ranked <- Filter(Negate(is.null), lapply(ok, .fold_ranked, mode = mode))
  signed <- lapply(ok, function(f) f$gene_stats)
  if (length(ranked) < 2)
    stop("Need >= 2 successful iterations with gene rankings for discovery.")

  universe <- sort(unique(unlist(ranked)))
  if (!is.null(gene_universe)) universe <- intersect(universe, gene_universe)
  R <- length(ranked)
  rank_mat <- matrix(missing_penalty, nrow = length(universe), ncol = R,
                     dimnames = list(universe, NULL))
  presence <- matrix(0L, nrow = length(universe), ncol = R, dimnames = list(universe, NULL))
  signs    <- matrix(0,  nrow = length(universe), ncol = R, dimnames = list(universe, NULL))
  for (r in seq_len(R)) {
    L <- ranked[[r]]; n <- length(L)
    inu <- L[L %in% universe]
    pr <- (match(inu, L)) / n                          # percentile rank in [0,1]
    rank_mat[inu, r] <- pr
    presence[inu, r] <- 1L
    st <- signed[[r]]
    if (!is.null(st)) {
      common <- intersect(inu, names(st))
      signs[common, r] <- sign(st[common])
    }
  }
  list(rank_mat = rank_mat, presence = presence, signs = signs,
       universe = universe, n_runs = R)
}

# ── Step 2 — robust genes (RRA + stability metrics) ──────────────────────────
#' @keywords internal
detect_robust_genes <- function(norm, top_frac = 0.1, min_freq = 0.5) {
  rm <- norm$rank_mat; pres <- norm$presence; sg <- norm$signs; R <- norm$n_runs
  # Reconstruct per-run best-first lists for RRA from the percentile ranks.
  lists <- lapply(seq_len(R), function(r) {
    present <- rownames(rm)[pres[, r] == 1]
    present[order(rm[present, r])]
  })
  rra <- rra_aggregate(lists)
  if (is.null(rra)) return(NULL)
  names(rra)[1] <- "Gene"

  # Top region = top fraction, but at least 3 genes (so small universes work).
  thr <- max(top_frac, 3 / nrow(rm))
  freq_topk <- rowMeans(rm <= thr)                     # fraction of runs in top region
  med_rank  <- apply(rm, 1, stats::median)
  rank_var  <- apply(rm, 1, stats::var)
  dir_cons  <- apply(sg, 1, function(s) { s <- s[s != 0]; if (!length(s)) 0 else abs(mean(s)) })

  g <- rra$Gene
  out <- data.frame(
    Gene = g, RRA_Score = rra$RRA_Score, Rank = rra$Rank,
    Freq_topK = round(freq_topk[g], 3), Median_Rank = round(med_rank[g], 3),
    Rank_Var = round(rank_var[g], 4), Direction_Consistency = round(dir_cons[g], 3),
    stringsAsFactors = FALSE)
  out <- out[out$Freq_topK >= min_freq, , drop = FALSE]
  out[order(out$Rank), , drop = FALSE]
}

# ── Step 3 — gene-gene stability graph ───────────────────────────────────────
#' @keywords internal
detect_stability_graph <- function(genes, norm, top_k = 0.1, tau = 0.25) {
  rm <- norm$rank_mat[genes, , drop = FALSE]
  pres <- norm$presence[genes, , drop = FALSE]
  sg <- norm$signs[genes, , drop = FALSE]
  top_thr <- max(top_k, 3 / nrow(norm$rank_mat))       # floor for small universes
  topset <- rm <= top_thr                              # in top-K region per run
  n <- length(genes)
  W <- matrix(0, n, n, dimnames = list(genes, genes))
  for (i in seq_len(n - 1)) for (j in (i + 1):n) {
    ti <- topset[i, ]; tj <- topset[j, ]
    inter <- sum(ti & tj); uni <- sum(ti | tj)
    jacc <- if (uni == 0) 0 else inter / uni           # 3.1 co-occurrence
    both <- pres[i, ] == 1 & pres[j, ] == 1
    if (sum(both) < 2) next
    prox <- exp(-stats::median(abs(rm[i, both] - rm[j, both])) / tau)   # 3.2
    # 3.3 rank-profile similarity. Correlation is undefined for near-constant
    # (very stable) profiles — treat those as similar (stability is the point);
    # rank proximity above already encodes how close their levels are.
    si <- stats::sd(rm[i, both]); sj <- stats::sd(rm[j, both])
    if (si < 1e-9 || sj < 1e-9) {
      cc <- 1
    } else {
      cc <- suppressWarnings(stats::cor(rm[i, both], rm[j, both]))
      cc <- if (is.na(cc)) 0 else max(0, cc)
    }
    dagr <- mean(sg[i, both] == sg[j, both])                            # 3.4
    w <- jacc * prox * cc * dagr
    if (w > 1e-6) { W[i, j] <- w; W[j, i] <- w }
  }
  W
}

# ── Step 4 — community detection ─────────────────────────────────────────────
#' @keywords internal
detect_communities <- function(adj, method = c("leiden", "louvain", "hclust"),
                               resolution = 1) {
  method <- match.arg(method)
  genes <- rownames(adj)
  if (length(genes) < 2) return(stats::setNames(rep(1L, length(genes)), genes))

  if (method %in% c("leiden", "louvain") && requireNamespace("igraph", quietly = TRUE)) {
    g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected",
                                             weighted = TRUE, diag = FALSE)
    cl <- tryCatch({
      if (method == "leiden")
        igraph::cluster_leiden(g, objective_function = "modularity",
                               resolution_parameter = resolution, weights = igraph::E(g)$weight)
      else
        igraph::cluster_louvain(g, weights = igraph::E(g)$weight, resolution = resolution)
    }, error = function(e) NULL)
    if (!is.null(cl))
      return(stats::setNames(as.integer(igraph::membership(cl)), genes))
  }
  # Dependency-free fallback: hierarchical clustering on 1 - similarity.
  d <- stats::as.dist(1 - adj)
  hc <- stats::hclust(d, method = "average")
  k <- max(2, min(length(genes) - 1, round(sqrt(length(genes)))))
  stats::setNames(stats::cutree(hc, k = k), genes)
}

# Group-stability score per run = median signed normalized score of members.
.group_stability <- function(module_genes, norm) {
  rm <- norm$rank_mat[module_genes, , drop = FALSE]
  sg <- norm$signs[module_genes, , drop = FALSE]
  score <- (1 - rm) * sg                               # high & directional -> large |score|
  per_run <- apply(score, 2, stats::median)
  stats::median(abs(per_run))
}

# ── Step 5 — leave-one-gene-out refinement ───────────────────────────────────
#' @keywords internal
detect_refine_loo <- function(module_genes, norm) {
  full <- .group_stability(module_genes, norm)
  if (length(module_genes) < 2 || full == 0)
    return(data.frame(Gene = module_genes, Contribution = 0,
                      Class = "core", stringsAsFactors = FALSE))
  contrib <- vapply(module_genes, function(g) {
    (full - .group_stability(setdiff(module_genes, g), norm)) / full
  }, numeric(1))
  pos <- contrib[contrib > 0]
  thr <- if (length(pos)) stats::median(pos) else 0
  cls <- ifelse(contrib < 0, "discarded",
                ifelse(contrib >= thr, "core", "supporting"))
  data.frame(Gene = module_genes, Contribution = round(contrib, 4),
             Class = cls, stringsAsFactors = FALSE)
}

# ── Step 7 — split / merge ───────────────────────────────────────────────────
#' @keywords internal
detect_split_merge <- function(membership, adj, norm,
                               merge_thr = 0.15, density_thr = 0.05) {
  mods <- split(names(membership), membership)
  mods <- Filter(function(m) length(m) >= 2, mods)
  if (length(mods) <= 1) return(membership)

  # Merge pairs whose between-module mean edge weight exceeds threshold.
  repeat {
    nm <- length(mods); merged <- FALSE
    if (nm < 2) break
    for (i in seq_len(nm - 1)) {
      for (j in (i + 1):nm) {
        bw <- mean(adj[mods[[i]], mods[[j]], drop = FALSE])
        if (bw > merge_thr) {
          mods[[i]] <- union(mods[[i]], mods[[j]]); mods[[j]] <- NULL
          merged <- TRUE; break
        }
      }
      if (merged) break
    }
    if (!merged) break
  }
  # Split modules with very low internal edge density via sub-clustering.
  final <- list()
  for (m in mods) {
    sub <- adj[m, m, drop = FALSE]
    dens <- if (length(m) < 2) 1 else mean(sub[upper.tri(sub)])
    if (length(m) >= 4 && dens < density_thr) {
      memb <- detect_communities(sub, method = "hclust")
      for (s in split(names(memb), memb)) final[[length(final) + 1L]] <- s
    } else final[[length(final) + 1L]] <- m
  }
  out <- integer(0)
  for (i in seq_along(final)) out[final[[i]]] <- i
  out
}

# ── Step 6 — validation (leave-one-run-out + bootstrap) ──────────────────────
#' @keywords internal
detect_validate <- function(modules, folds, mode, gene_universe, params,
                            n_boot = 100, max_boot = 100, seed = 1) {
  ok <- .detect_ok_folds(folds)
  R <- length(ok)
  best_jacc <- function(set, ref_sets) {
    if (!length(ref_sets)) return(0)
    max(vapply(ref_sets, function(s)
      length(intersect(set, s)) / max(1, length(union(set, s))), numeric(1)))
  }
  rerun <- function(sub_folds) {
    res <- tryCatch({
      nz <- detect_normalize_ranks(sub_folds, mode, gene_universe)
      rg <- detect_robust_genes(nz, params$top_k_frac, params$min_freq)
      if (is.null(rg) || nrow(rg) < 2) return(list())
      g <- utils::head(rg$Gene, params$top_genes)
      adj <- detect_stability_graph(g, nz, top_k = params$top_k_frac, tau = params$tau)
      memb <- detect_communities(adj, resolution = params$resolution)
      split(names(memb), memb)
    }, error = function(e) list())
    res
  }
  # Leave-one-run-out recovery.
  loo_sets <- lapply(seq_len(R), function(r) rerun(ok[-r]))
  recovery <- vapply(modules, function(mset) {
    mean(vapply(loo_sets, function(sets) best_jacc(mset, sets) >= 0.5, logical(1)))
  }, numeric(1))
  # Bootstrap mean Jaccard.
  boot_j <- rep(NA_real_, length(modules))
  nb <- min(n_boot, max_boot)
  if (nb > 0) {
    set.seed(seed)
    boots <- lapply(seq_len(nb), function(b) rerun(ok[sample.int(R, R, replace = TRUE)]))
    boot_j <- vapply(modules, function(mset)
      mean(vapply(boots, function(sets) best_jacc(mset, sets), numeric(1))), numeric(1))
  }
  data.frame(Recovery_Freq = round(recovery, 3),
             Boot_Mean_Jaccard = round(boot_j, 3))
}

# ── Step 8 — annotation vs known groups ──────────────────────────────────────
#' @keywords internal
detect_annotate <- function(module_sets, groups) {
  do.call(rbind, lapply(module_sets, function(mset) {
    if (!length(groups))
      return(data.frame(Best_Known_Group = NA, Jaccard = 0, Overlap = 0,
                        Label = "novel", stringsAsFactors = FALSE))
    jac <- vapply(groups, function(g)
      length(intersect(mset, g)) / max(1, length(union(mset, g))), numeric(1))
    cont <- vapply(groups, function(g)
      length(intersect(mset, g)) / max(1, length(mset)), numeric(1))
    bi <- which.max(jac)
    overlap <- length(intersect(mset, groups[[bi]]))
    label <- if (jac[bi] >= 0.5) "rediscovered"
             else if (max(cont) >= 0.6) "fragment" else "novel"
    data.frame(Best_Known_Group = names(groups)[bi], Jaccard = round(jac[bi], 3),
               Overlap = overlap, Label = label, stringsAsFactors = FALSE)
  }))
}

#' Discover robust gene modules de novo from an MCCV result
#'
#' Runs the full G-S-M detect pipeline (normalize ranks, robust-gene detection
#' via RRA, gene-gene stability graph, community detection, leave-one-gene-out
#' refinement, split/merge, bootstrap validation, and annotation against the
#' dataset's known groups).
#'
#' @param result An MCCV result with `$folds` (from [run_mccv()] or the app).
#' @param dataset The scored `gsm_dataset` (its `$groups` are the known reference).
#' @param mode Direction mode: "abs" (default), "up", or "down".
#' @param top_genes Cap on robust genes entering the graph (default 200).
#' @param top_k_frac Top fraction defining the "top region" (default 0.1).
#' @param min_freq Minimum top-region frequency to call a gene robust (0.5).
#' @param tau Rank-proximity scale (default 0.25).
#' @param resolution Community-detection resolution (default 1).
#' @param n_boot Bootstrap iterations for validation (default 100; 0 to skip).
#' @param seed Random seed for the bootstrap.
#' @return A list: `robust_genes`, `modules`, `members`, `graph`, `params`.
#' @export
gsm_detect <- function(result, dataset, mode = c("abs", "up", "down"),
                       top_genes = 200, top_k_frac = 0.1, min_freq = 0.5,
                       tau = 0.25, resolution = 1, n_boot = 100, seed = 1) {
  mode <- match.arg(mode)
  stopifnot(inherits(dataset, "gsm_dataset"))
  if (is.null(result$folds)) stop("gsm_detect() needs an MCCV result with $folds.")
  universe <- unique(unlist(dataset$groups))

  norm <- detect_normalize_ranks(result$folds, mode, gene_universe = universe)
  robust <- detect_robust_genes(norm, top_frac = top_k_frac, min_freq = min_freq)
  if (is.null(robust) || nrow(robust) < 2)
    stop("Too few robust genes — lower 'min frequency' or run more iterations.")

  truncated <- nrow(robust) > top_genes
  genes <- utils::head(robust$Gene, top_genes)
  if (truncated) message(sprintf("Capping graph to top %d of %d robust genes.",
                                  top_genes, nrow(robust)))

  adj <- detect_stability_graph(genes, norm, top_k = top_k_frac, tau = tau)
  membership <- detect_communities(adj, resolution = resolution)
  membership <- detect_split_merge(membership, adj, norm)

  params0 <- list(mode = mode, top_genes = top_genes, top_k_frac = top_k_frac,
                  min_freq = min_freq, tau = tau, resolution = resolution,
                  n_boot = n_boot, n_runs = norm$n_runs, n_robust = nrow(robust),
                  truncated = truncated)

  mod_ids <- sort(unique(membership))
  members_rows <- list(); mod_rows <- list(); mod_sets <- list()
  for (mi in mod_ids) {
    mg <- names(membership)[membership == mi]
    if (length(mg) < 2) next
    loo <- detect_refine_loo(mg, norm)
    keep <- loo$Gene[loo$Class != "discarded"]          # drop genes that hurt the module
    if (length(keep) < 2) keep <- mg
    nm <- sprintf("M%02d", length(mod_sets) + 1L)
    loo$Module <- nm
    members_rows[[nm]] <- loo[, c("Module", "Gene", "Contribution", "Class")]
    mod_sets[[nm]] <- keep
    mod_rows[[nm]] <- data.frame(
      Module = nm, Size = length(keep),
      Core_Genes = paste(loo$Gene[loo$Class == "core"], collapse = ", "),
      Supporting_Genes = paste(loo$Gene[loo$Class == "supporting"], collapse = ", "),
      Mean_Stability = round(.group_stability(keep, norm), 4),
      stringsAsFactors = FALSE)
  }
  modules <- do.call(rbind, mod_rows)
  members <- do.call(rbind, members_rows)
  if (is.null(modules)) {
    # No multi-gene module survived (e.g. too few/sparse robust genes). Return
    # the robust-gene evidence anyway with an explanatory note (no hard error).
    params0$note <- "No multi-gene modules detected — too few or sparsely co-occurring robust genes. Try more iterations, a higher 'top genes', a lower 'min frequency', or turn gene shaving off."
    return(list(robust_genes = robust,
                modules = data.frame(), members = data.frame(),
                graph = list(adjacency = adj, membership = membership),
                params = params0))
  }

  params <- params0
  val <- detect_validate(mod_sets, result$folds, mode, universe, params,
                         n_boot = n_boot, seed = seed)
  ann <- detect_annotate(mod_sets, dataset$groups)
  modules <- cbind(modules, val, ann)
  rownames(modules) <- NULL; rownames(members) <- NULL

  list(robust_genes = robust, modules = modules, members = members,
       graph = list(adjacency = adj, membership = membership), params = params)
}
