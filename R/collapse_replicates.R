#' @title Replicate Collapsing / Preprocessing
#'
#' @description
#' Collapse selected replicate samples into single representative samples
#' using a user-chosen aggregation method. This is a clean preprocessing step
#' that runs **before** any G-S-M scoring or modelling.
#'
#' ## Technical vs Biological Replicates
#'
#' - **Technical replicates** arise from repeated measurement of the same
#'   biological sample (e.g. multiple sequencing runs of the same library).
#'   These should typically be collapsed *before* analysis. For RNA-seq count
#'   data, **summing** is recommended because counts from the same library are
#'   additive.
#'
#' - **Biological replicates** arise from distinct biological samples (e.g.
#'   different patients, different mice). These should generally **not** be
#'   collapsed. Instead, use grouped cross-validation so that all samples from
#'   the same subject stay in the same fold.
#'
#' ## Recommendations
#'
#' | Data type | Replicate type | Recommended method |
#' |-----------|---------------|-------------------|
#' | RNA-seq raw counts | Technical | `"sum"` |
#' | Normalised expression | Technical | `"mean"` or `"median"` |
#' | Microarray intensity | Technical | `"mean"` or `"median"` |
#' | Any | Biological | Do **not** collapse; use grouped CV |
#'
#' @name collapse_replicates
NULL


#' Collapse replicate samples in an expression matrix
#'
#' For each replicate group, the corresponding sample columns are aggregated
#' into a single new column. Samples not belonging to any replicate group (or
#' not selected for collapsing) are kept unchanged.
#'
#' @param expr Numeric matrix with features/genes as rows and samples as
#'   columns. Row names are feature IDs; column names are sample IDs.
#' @param metadata A `data.frame` with one row per sample. Must include a
#'   column matching `sample_id_col` whose values correspond to `colnames(expr)`.
#' @param sample_id_col Name of the column in `metadata` that holds sample
#'   identifiers. Default `"sample_id"`.
#' @param replicate_group_col Name of the column in `metadata` that defines
#'   which samples are replicates of each other. Samples sharing the same
#'   value are collapsed together. Default `"replicate_group_id"`.
#' @param method Aggregation method: one of `"mean"`, `"median"`, `"sum"`,
#'   `"min"`, `"max"`, `"first"`. Default `"mean"`.
#' @param selected_groups Optional character vector of replicate group IDs to
#'   collapse. If `NULL` (default), all groups are collapsed.
#' @param replicate_type_col Optional column name indicating replicate type
#'   (e.g. `"technical"`, `"biological"`).
#' @param replicate_type_filter Optional value(s) to filter on in
#'   `replicate_type_col`. Only samples matching this filter are collapsed.
#' @param metadata_conflict How to handle conflicting metadata values within a
#'   replicate group. One of:
#'   - `"error"` (default): stop with an error if any non-ID metadata values
#'     differ within a group.
#'   - `"NA"`: set conflicting values to `NA`.
#'   - `"concat"`: concatenate conflicting values with `"; "`.
#'   - `"first"`: keep the value from the first sample in the group.
#' @param keep_original If `TRUE`, retain the original uncollapsed samples
#'   alongside the new collapsed columns. Default `FALSE`.
#' @param na_rm Logical. Should `NA` values be removed when aggregating?
#'   Default `TRUE`.
#' @param output_provenance If `TRUE` (default), include a provenance table
#'   in the output describing which samples were collapsed and how.
#' @return A list with components:
#'   \describe{
#'     \item{expr}{Numeric matrix with collapsed samples.}
#'     \item{metadata}{Updated `data.frame` with one row per output sample.}
#'     \item{provenance}{A `data.frame` recording: collapsed sample ID,
#'       original sample IDs, method, replicate group, and count.
#'       `NULL` if `output_provenance = FALSE`.}
#'     \item{warnings}{Character vector of warnings generated during
#'       collapsing.}
#'   }
#' @export
#'
#' @examples
#' # --- Example 1: Collapse technical RNA-seq count replicates by sum ---
#' expr <- matrix(c(10, 20, 12, 22, 5, 15), nrow = 2,
#'                dimnames = list(c("GeneA", "GeneB"), c("S1", "S2", "S3")))
#' meta <- data.frame(sample_id = c("S1", "S2", "S3"),
#'                    replicate_group_id = c("Patient1", "Patient1", "Patient2"),
#'                    condition = c("treat", "treat", "ctrl"))
#' result <- collapse_replicates(expr, meta, method = "sum")
#'
#' # --- Example 2: Collapse normalised expression replicates by mean ---
#' result <- collapse_replicates(expr, meta, method = "mean")
#'
#' # --- Example 3: Collapse only selected replicate groups ---
#' result <- collapse_replicates(expr, meta, method = "mean",
#'                               selected_groups = "Patient1")
#'
#' # --- Example 4: Error on conflicting labels ---
#' meta_bad <- meta
#' meta_bad$condition <- c("treat", "ctrl", "ctrl")
#' try(collapse_replicates(expr, meta_bad, method = "mean"))
collapse_replicates <- function(
    expr,
    metadata,
    sample_id_col       = "sample_id",
    replicate_group_col = "replicate_group_id",
    method              = c("mean", "median", "sum", "min", "max", "first"),
    selected_groups     = NULL,
    replicate_type_col  = NULL,
    replicate_type_filter = NULL,
    metadata_conflict   = c("error", "NA", "concat", "first"),
    keep_original       = FALSE,
    na_rm               = TRUE,
    output_provenance   = TRUE
) {
  # --- Argument matching ---
  method <- match.arg(method)
  metadata_conflict <- match.arg(metadata_conflict)
  warn_msgs <- character(0)

  # --- Validate matrix ---
  if (!is.matrix(expr) || !is.numeric(expr)) {
    stop("`expr` must be a numeric matrix (features x samples).")
  }
  if (is.null(rownames(expr))) {
    stop("`expr` must have row names (feature IDs).")
  }
  if (is.null(colnames(expr))) {
    stop("`expr` must have column names (sample IDs).")
  }

  # --- Validate metadata ---
  if (!is.data.frame(metadata)) {
    stop("`metadata` must be a data.frame.")
  }
  if (!sample_id_col %in% names(metadata)) {
    stop("Column '", sample_id_col, "' not found in metadata.")
  }
  if (!replicate_group_col %in% names(metadata)) {
    stop("Column '", replicate_group_col, "' not found in metadata.")
  }

  sample_ids <- as.character(metadata[[sample_id_col]])

  # Check for duplicate sample IDs in metadata
  if (anyDuplicated(sample_ids)) {
    stop("Duplicate sample IDs in metadata column '", sample_id_col, "'.")
  }

  # --- Cross-validate matrix and metadata ---
  mat_samples <- colnames(expr)

  missing_from_matrix <- setdiff(sample_ids, mat_samples)
  if (length(missing_from_matrix) > 0) {
    stop("Metadata contains sample IDs not present in the expression matrix: ",
         paste(head(missing_from_matrix, 5), collapse = ", "))
  }
  missing_from_metadata <- setdiff(mat_samples, sample_ids)
  if (length(missing_from_metadata) > 0) {
    stop("Expression matrix contains samples not present in metadata: ",
         paste(head(missing_from_metadata, 5), collapse = ", "))
  }

  # --- Determine which groups to collapse ---
  group_vec <- as.character(metadata[[replicate_group_col]])

  # Apply replicate type filter if specified
  collapse_mask <- rep(TRUE, nrow(metadata))
  if (!is.null(replicate_type_col)) {
    if (!replicate_type_col %in% names(metadata)) {
      stop("Column '", replicate_type_col, "' not found in metadata.")
    }
    if (!is.null(replicate_type_filter)) {
      type_vec <- as.character(metadata[[replicate_type_col]])
      collapse_mask <- type_vec %in% replicate_type_filter
    }
  }

  # Apply selected_groups filter
  if (!is.null(selected_groups)) {
    collapse_mask <- collapse_mask & (group_vec %in% selected_groups)
  }

  # Identify groups to collapse (groups with >0 samples selected)
  groups_to_collapse <- unique(group_vec[collapse_mask])

  # Samples that are NOT being collapsed
  non_collapse_idx <- which(!collapse_mask | !(group_vec %in% groups_to_collapse))

  # --- Warn about singleton groups ---
  for (g in groups_to_collapse) {
    g_idx <- which(group_vec == g & collapse_mask)
    if (length(g_idx) == 1) {
      msg <- paste0("Replicate group '", g,
                     "' contains only 1 sample — no collapsing needed.")
      warn_msgs <- c(warn_msgs, msg)
    }
  }

  # --- Data-type warnings ---
  all_integer <- all(expr == floor(expr), na.rm = TRUE) &&
    !any(expr < 0, na.rm = TRUE)
  if (method == "sum" && !all_integer) {
    msg <- paste0("Using method='sum' on data that may not be integer counts. ",
                  "Sum is recommended for raw RNA-seq count data.")
    warn_msgs <- c(warn_msgs, msg)
  }
  if (method %in% c("mean", "median") && all_integer &&
      max(expr, na.rm = TRUE) > 100) {
    msg <- paste0("Using method='", method, "' on data that appears to be ",
                  "raw counts. For technical RNA-seq count replicates, ",
                  "consider method='sum' instead.")
    warn_msgs <- c(warn_msgs, msg)
  }

  # --- Get the aggregation function ---
  agg_fn <- switch(method,
    mean   = \(x, na.rm) mean(x, na.rm = na.rm),
    median = \(x, na.rm) stats::median(x, na.rm = na.rm),
    sum    = \(x, na.rm) sum(x, na.rm = na.rm),
    min    = \(x, na.rm) min(x, na.rm = na.rm),
    max    = \(x, na.rm) max(x, na.rm = na.rm),
    first  = \(x, na.rm) x[1]
  )

  # --- Build output ---
  # Metadata columns to resolve (excluding sample_id and replicate_group)
  meta_resolve_cols <- setdiff(names(metadata),
                               c(sample_id_col, replicate_group_col))

  # Provenance tracking
  prov_rows <- list()

  # Collected output columns (expression) and metadata rows
  out_expr_cols <- list()
  out_meta_rows <- list()
  out_sample_order <- character(0)

  # Process each unique group
  processed_samples <- character(0)

  for (g in groups_to_collapse) {
    g_mask <- group_vec == g & collapse_mask
    g_sample_ids <- sample_ids[g_mask]

    if (length(g_sample_ids) == 0) next

    # --- Collapse expression ---
    g_expr <- expr[, g_sample_ids, drop = FALSE]

    if (length(g_sample_ids) == 1) {
      collapsed_col <- g_expr[, 1, drop = TRUE]
    } else {
      collapsed_col <- apply(g_expr, 1, agg_fn, na.rm = na_rm)
    }

    # New sample ID
    new_id <- paste0(g, "_collapsed_", method)
    names(collapsed_col) <- rownames(expr)
    out_expr_cols[[new_id]] <- collapsed_col
    out_sample_order <- c(out_sample_order, new_id)

    # --- Resolve metadata ---
    g_meta <- metadata[g_mask, , drop = FALSE]
    new_meta_vals <- list()
    new_meta_vals[[sample_id_col]] <- new_id
    new_meta_vals[[replicate_group_col]] <- g

    for (col in meta_resolve_cols) {
      vals <- as.character(g_meta[[col]])
      unique_vals <- unique(vals[!is.na(vals)])

      if (length(unique_vals) <= 1) {
        # All agree (or all NA)
        new_meta_vals[[col]] <- if (length(unique_vals) == 1) {
          unique_vals
        } else {
          NA_character_
        }
      } else {
        # Conflict
        switch(metadata_conflict,
          error = stop(
            "Conflicting metadata values in replicate group '", g,
            "' for column '", col, "': ",
            paste(unique_vals, collapse = " vs "), ".\n",
            "Use metadata_conflict='NA', 'concat', or 'first' to resolve."
          ),
          "NA" = {
            new_meta_vals[[col]] <- NA_character_
          },
          concat = {
            new_meta_vals[[col]] <- paste(unique_vals, collapse = "; ")
          },
          first = {
            new_meta_vals[[col]] <- vals[1]
          }
        )
      }
    }
    new_meta_row <- as.data.frame(new_meta_vals, stringsAsFactors = FALSE)
    out_meta_rows[[new_id]] <- new_meta_row

    # --- Provenance ---
    if (output_provenance) {
      prov_rows[[new_id]] <- data.frame(
        collapsed_sample_id = new_id,
        original_sample_ids = paste(g_sample_ids, collapse = "; "),
        replicate_group     = g,
        method              = method,
        n_replicates        = length(g_sample_ids),
        stringsAsFactors    = FALSE
      )
    }

    processed_samples <- c(processed_samples, g_sample_ids)
  }

  # --- Non-collapsed samples ---
  uncollapsed_ids <- setdiff(sample_ids, processed_samples)
  for (sid in uncollapsed_ids) {
    out_expr_cols[[sid]] <- expr[, sid, drop = TRUE]
    out_sample_order <- c(out_sample_order, sid)

    row_idx <- which(sample_ids == sid)
    out_meta_rows[[sid]] <- metadata[row_idx, , drop = FALSE]
  }

  # --- Keep originals if requested ---
  if (keep_original && length(processed_samples) > 0) {
    for (sid in processed_samples) {
      orig_id <- paste0(sid, "_original")
      out_expr_cols[[orig_id]] <- expr[, sid, drop = TRUE]
      out_sample_order <- c(out_sample_order, orig_id)

      row_idx <- which(sample_ids == sid)
      orig_meta <- metadata[row_idx, , drop = FALSE]
      orig_meta[[sample_id_col]] <- orig_id
      out_meta_rows[[orig_id]] <- orig_meta
    }
  }

  # --- Assemble output matrix ---
  out_expr <- do.call(cbind, out_expr_cols[out_sample_order])
  colnames(out_expr) <- out_sample_order
  rownames(out_expr) <- rownames(expr)

  # --- Assemble output metadata ---
  out_metadata <- do.call(rbind, out_meta_rows[out_sample_order])
  rownames(out_metadata) <- NULL

  # Restore original column types where possible
  for (col in names(out_metadata)) {
    orig_col <- metadata[[col]]
    if (is.numeric(orig_col)) {
      out_metadata[[col]] <- as.numeric(out_metadata[[col]])
    } else if (is.logical(orig_col)) {
      out_metadata[[col]] <- as.logical(out_metadata[[col]])
    } else if (is.factor(orig_col)) {
      out_metadata[[col]] <- factor(out_metadata[[col]], levels = levels(orig_col))
    }
  }

  # --- Assemble provenance ---
  provenance <- NULL
  if (output_provenance && length(prov_rows) > 0) {
    provenance <- do.call(rbind, prov_rows)
    rownames(provenance) <- NULL

    # Attach original metadata as a nested attribute
    orig_meta_by_group <- list()
    for (g in groups_to_collapse) {
      g_mask <- group_vec == g & collapse_mask
      if (any(g_mask)) {
        orig_meta_by_group[[g]] <- metadata[g_mask, , drop = FALSE]
      }
    }
    attr(provenance, "original_metadata") <- orig_meta_by_group
  }

  # --- Issue collected warnings ---
  for (w in warn_msgs) {
    warning(w, call. = FALSE)
  }

  list(
    expr       = out_expr,
    metadata   = out_metadata,
    provenance = provenance,
    warnings   = warn_msgs
  )
}
