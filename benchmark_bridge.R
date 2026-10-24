source("R/bridge.R")
source("R/utils.R")
source("R/detect.R")

suppressPackageStartupMessages(library(igraph))
options(error = traceback)

run_scaling_test <- function(p_list) {
  cat("\n=== BRIDGE SCALING BENCHMARK ===\n")
  
  results <- data.frame(
    Genes_p = integer(),
    Runtime_sec = numeric(),
    Edges = integer(),
    Modules = integer(),
    stringsAsFactors = FALSE
  )
  
  for (p in p_list) {
    n <- 40 
    set.seed(42)
    X <- matrix(rnorm(p * n), nrow = p, dimnames = list(paste0("G", 1:p), paste0("S", 1:n)))
    y <- factor(rep(c("A", "B"), each = n/2))
    
    edge_df <- data.frame(
      gene_i = paste0("G", 1:(p-1)),
      gene_j = paste0("G", 2:p),
      weight = 0.5
    )
    
    cat(sprintf("\nRunning BRIDGE with p = %d...\n", p))
    t0 <- Sys.time()
    
    res <- tryCatch({
      bridge_grouping(
        X = X, 
        y = y, 
        B = edge_df, 
        prefilter = "none",  
        precision_method = "glasso", # use glasso since huge might not be installed
        resolution = 1,
        n_boot = 5,
        tau = 0.5,
        verbose = FALSE,
        details = TRUE
      )
    }, error = function(e) {
      cat("Error:", e$message, "\n")
      NULL
    })
    
    t1 <- Sys.time()
    runtime <- as.numeric(difftime(t1, t0, units = "secs"))
    
    if (!is.null(res)) {
      mod_count <- length(res)
      results <- rbind(results, data.frame(
        Genes_p = p,
        Runtime_sec = round(runtime, 2),
        Edges = NA,
        Modules = mod_count
      ))
      cat(sprintf("DONE. Runtime: %.2f sec | Modules: %d\n", runtime, mod_count))
    }
  }
  
  cat("\n=== FINAL SCALING RESULTS ===\n")
  print(results)
}

run_scaling_test(c(100, 300))
