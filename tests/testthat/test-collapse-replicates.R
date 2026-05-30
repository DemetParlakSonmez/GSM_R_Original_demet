# Tests for collapse_replicates()

# ──────────────────────────────────────────────────────────────────────────────
# Test fixtures
# ──────────────────────────────────────────────────────────────────────────────

make_count_fixture <- function() {
  # Two technical replicates for Patient1, one sample for Patient2

  expr <- matrix(
    c(10, 20,   # S1: GeneA=10, GeneB=20
      12, 22,   # S2: GeneA=12, GeneB=22
       5, 15),  # S3: GeneA= 5, GeneB=15
    nrow = 2, byrow = FALSE,
    dimnames = list(c("GeneA", "GeneB"), c("S1", "S2", "S3"))
  )
  meta <- data.frame(
    sample_id          = c("S1", "S2", "S3"),
    replicate_group_id = c("Patient1", "Patient1", "Patient2"),
    condition          = c("treat", "treat", "ctrl"),
    stringsAsFactors   = FALSE
  )
  list(expr = expr, meta = meta)
}


make_three_rep_fixture <- function() {
  expr <- matrix(
    c(10, 20,    # S1
      12, 22,    # S2
      14, 18,    # S3
       5, 15),   # S4
    nrow = 2, byrow = FALSE,
    dimnames = list(c("GeneA", "GeneB"), c("S1", "S2", "S3", "S4"))
  )
  meta <- data.frame(
    sample_id          = c("S1", "S2", "S3", "S4"),
    replicate_group_id = c("G1", "G1", "G1", "G2"),
    condition          = c("treat", "treat", "treat", "ctrl"),
    stringsAsFactors   = FALSE
  )
  list(expr = expr, meta = meta)
}


# ──────────────────────────────────────────────────────────────────────────────
# Core aggregation tests
# ──────────────────────────────────────────────────────────────────────────────

test_that("sum correctly collapses two technical count replicates", {
  fix <- make_count_fixture()
  result <- suppressWarnings(
    collapse_replicates(fix$expr, fix$meta, method = "sum")
  )

  # Patient1: S1+S2 = (10+12, 20+22) = (22, 42)
  expect_equal(result$expr["GeneA", "Patient1_collapsed_sum"], 22)
  expect_equal(result$expr["GeneB", "Patient1_collapsed_sum"], 42)

  # Patient2: singleton, (5, 15) — collapsed but values unchanged
  expect_equal(result$expr["GeneA", "Patient2_collapsed_sum"], 5)
  expect_equal(result$expr["GeneB", "Patient2_collapsed_sum"], 15)

  # Output has 2 columns

  expect_equal(ncol(result$expr), 2)
})


test_that("mean correctly collapses two replicates", {
  fix <- make_count_fixture()
  result <- suppressWarnings(
    collapse_replicates(fix$expr, fix$meta, method = "mean")
  )

  # Patient1: mean(10,12)=11, mean(20,22)=21
  expect_equal(result$expr["GeneA", "Patient1_collapsed_mean"], 11)
  expect_equal(result$expr["GeneB", "Patient1_collapsed_mean"], 21)
})


test_that("median works with three replicates (odd count)", {
  fix <- make_three_rep_fixture()
  result <- suppressWarnings(
    collapse_replicates(fix$expr, fix$meta, method = "median")
  )

  # G1: median(10,12,14)=12, median(20,22,18)=20
  expect_equal(result$expr["GeneA", "G1_collapsed_median"], 12)
  expect_equal(result$expr["GeneB", "G1_collapsed_median"], 20)
})


test_that("median works with two replicates (even count)", {
  fix <- make_count_fixture()
  result <- suppressWarnings(
    collapse_replicates(fix$expr, fix$meta, method = "median")
  )

  # Patient1: median(10,12)=11, median(20,22)=21
  expect_equal(result$expr["GeneA", "Patient1_collapsed_median"], 11)
  expect_equal(result$expr["GeneB", "Patient1_collapsed_median"], 21)
})


test_that("min selects minimum across replicates", {
  fix <- make_count_fixture()
  result <- suppressWarnings(
    collapse_replicates(fix$expr, fix$meta, method = "min")
  )

  expect_equal(result$expr["GeneA", "Patient1_collapsed_min"], 10)
  expect_equal(result$expr["GeneB", "Patient1_collapsed_min"], 20)
})


test_that("max selects maximum across replicates", {
  fix <- make_count_fixture()
  result <- suppressWarnings(
    collapse_replicates(fix$expr, fix$meta, method = "max")
  )

  expect_equal(result$expr["GeneA", "Patient1_collapsed_max"], 12)
  expect_equal(result$expr["GeneB", "Patient1_collapsed_max"], 22)
})


test_that("first keeps the first replicate's values", {
  fix <- make_count_fixture()
  result <- suppressWarnings(
    collapse_replicates(fix$expr, fix$meta, method = "first")
  )

  # First sample in Patient1 is S1
  expect_equal(result$expr["GeneA", "Patient1_collapsed_first"], 10)
  expect_equal(result$expr["GeneB", "Patient1_collapsed_first"], 20)
})


# ──────────────────────────────────────────────────────────────────────────────
# Selective collapsing
# ──────────────────────────────────────────────────────────────────────────────

test_that("selected_groups collapses only specified groups", {
  fix <- make_count_fixture()
  result <- suppressWarnings(
    collapse_replicates(fix$expr, fix$meta, method = "mean",
                        selected_groups = "Patient1")
  )

  # Patient1 collapsed, Patient2 kept as-is (original sample S3)
  expect_true("Patient1_collapsed_mean" %in% colnames(result$expr))
  expect_true("S3" %in% colnames(result$expr))
  expect_false("Patient2_collapsed_mean" %in% colnames(result$expr))
  expect_equal(ncol(result$expr), 2)
})


test_that("replicate_type_filter collapses only matching types", {
  expr <- matrix(
    c(10, 20, 12, 22, 5, 15, 7, 17),
    nrow = 2,
    dimnames = list(c("GA", "GB"), c("S1", "S2", "S3", "S4"))
  )
  meta <- data.frame(
    sample_id          = c("S1", "S2", "S3", "S4"),
    replicate_group_id = c("G1", "G1", "G2", "G2"),
    replicate_type     = c("technical", "technical", "biological", "biological"),
    condition          = c("A", "A", "B", "B"),
    stringsAsFactors   = FALSE
  )

  result <- suppressWarnings(
    collapse_replicates(expr, meta, method = "mean",
                        replicate_type_col = "replicate_type",
                        replicate_type_filter = "technical")
  )

  # G1 collapsed (technical), G2 samples kept as-is (biological)
  expect_true("G1_collapsed_mean" %in% colnames(result$expr))
  expect_true("S3" %in% colnames(result$expr))
  expect_true("S4" %in% colnames(result$expr))
  expect_equal(ncol(result$expr), 3)
})


# ──────────────────────────────────────────────────────────────────────────────
# Samples without replicate group (singletons / non-collapsed)
# ──────────────────────────────────────────────────────────────────────────────

test_that("samples not in selected groups remain unchanged", {
  fix <- make_count_fixture()
  result <- suppressWarnings(
    collapse_replicates(fix$expr, fix$meta, method = "mean",
                        selected_groups = "Patient1")
  )

  # S3 should be unchanged
  expect_equal(result$expr["GeneA", "S3"], 5)
  expect_equal(result$expr["GeneB", "S3"], 15)

  # Metadata for S3 should be preserved
  s3_row <- result$metadata[result$metadata$sample_id == "S3", ]
  expect_equal(nrow(s3_row), 1)
  expect_equal(s3_row$condition, "ctrl")
})


# ──────────────────────────────────────────────────────────────────────────────
# Metadata conflict handling
# ──────────────────────────────────────────────────────────────────────────────

test_that("metadata conflicts trigger error by default", {
  fix <- make_count_fixture()
  fix$meta$condition <- c("treat", "ctrl", "ctrl")  # Conflict in Patient1

  expect_error(
    suppressWarnings(
      collapse_replicates(fix$expr, fix$meta, method = "mean")
    ),
    "Conflicting metadata"
  )
})


test_that("metadata_conflict='NA' sets conflicts to NA", {
  fix <- make_count_fixture()
  fix$meta$condition <- c("treat", "ctrl", "ctrl")

  result <- suppressWarnings(
    collapse_replicates(fix$expr, fix$meta, method = "mean",
                        metadata_conflict = "NA")
  )

  p1_row <- result$metadata[result$metadata$sample_id == "Patient1_collapsed_mean", ]
  expect_true(is.na(p1_row$condition))
})


test_that("metadata_conflict='concat' concatenates conflicting values", {
  fix <- make_count_fixture()
  fix$meta$condition <- c("treat", "ctrl", "ctrl")

  result <- suppressWarnings(
    collapse_replicates(fix$expr, fix$meta, method = "mean",
                        metadata_conflict = "concat")
  )

  p1_row <- result$metadata[result$metadata$sample_id == "Patient1_collapsed_mean", ]
  expect_true(grepl("treat", p1_row$condition))
  expect_true(grepl("ctrl", p1_row$condition))
})


test_that("metadata_conflict='first' keeps the first value", {
  fix <- make_count_fixture()
  fix$meta$condition <- c("treat", "ctrl", "ctrl")

  result <- suppressWarnings(
    collapse_replicates(fix$expr, fix$meta, method = "mean",
                        metadata_conflict = "first")
  )

  p1_row <- result$metadata[result$metadata$sample_id == "Patient1_collapsed_mean", ]
  expect_equal(p1_row$condition, "treat")
})


test_that("agreeing metadata values are preserved without conflict", {
  fix <- make_count_fixture()

  result <- suppressWarnings(
    collapse_replicates(fix$expr, fix$meta, method = "mean")
  )

  p1_row <- result$metadata[result$metadata$sample_id == "Patient1_collapsed_mean", ]
  expect_equal(p1_row$condition, "treat")
})


# ──────────────────────────────────────────────────────────────────────────────
# Provenance
# ──────────────────────────────────────────────────────────────────────────────

test_that("provenance records original sample IDs", {
  fix <- make_count_fixture()
  result <- suppressWarnings(
    collapse_replicates(fix$expr, fix$meta, method = "mean")
  )

  expect_false(is.null(result$provenance))
  expect_true(is.data.frame(result$provenance))

  p1_prov <- result$provenance[result$provenance$replicate_group == "Patient1", ]
  expect_equal(p1_prov$n_replicates, 2)
  expect_true(grepl("S1", p1_prov$original_sample_ids))
  expect_true(grepl("S2", p1_prov$original_sample_ids))
  expect_equal(p1_prov$method, "mean")
})


test_that("provenance stores original metadata as attribute", {
  fix <- make_count_fixture()
  result <- suppressWarnings(
    collapse_replicates(fix$expr, fix$meta, method = "sum")
  )

  orig_meta <- attr(result$provenance, "original_metadata")
  expect_type(orig_meta, "list")
  expect_true("Patient1" %in% names(orig_meta))
  expect_equal(nrow(orig_meta[["Patient1"]]), 2)
})


test_that("provenance can be disabled", {
  fix <- make_count_fixture()
  result <- suppressWarnings(
    collapse_replicates(fix$expr, fix$meta, method = "mean",
                        output_provenance = FALSE)
  )
  expect_null(result$provenance)
})


# ──────────────────────────────────────────────────────────────────────────────
# Validation errors
# ──────────────────────────────────────────────────────────────────────────────

test_that("error if sample IDs in metadata are missing from matrix", {
  fix <- make_count_fixture()
  fix$meta$sample_id[1] <- "MISSING"

  expect_error(
    collapse_replicates(fix$expr, fix$meta, method = "mean"),
    "not present in the expression matrix"
  )
})


test_that("error if matrix columns are missing from metadata", {
  fix <- make_count_fixture()
  # Add an extra column to the matrix that has no metadata row
  extra_col <- matrix(1:2, nrow = 2, dimnames = list(c("GeneA", "GeneB"), "EXTRA"))
  fix$expr <- cbind(fix$expr, extra_col)

  expect_error(
    collapse_replicates(fix$expr, fix$meta, method = "mean"),
    "not present in metadata"
  )
})


test_that("error if expr is not a numeric matrix", {
  fix <- make_count_fixture()
  expect_error(
    collapse_replicates(as.data.frame(fix$expr), fix$meta, method = "mean"),
    "numeric matrix"
  )
})


test_that("error if sample_id_col is missing", {
  fix <- make_count_fixture()
  expect_error(
    collapse_replicates(fix$expr, fix$meta, sample_id_col = "wrong",
                        method = "mean"),
    "not found"
  )
})


test_that("error if replicate_group_col is missing", {
  fix <- make_count_fixture()
  expect_error(
    collapse_replicates(fix$expr, fix$meta, replicate_group_col = "wrong",
                        method = "mean"),
    "not found"
  )
})


test_that("error on duplicate sample IDs in metadata", {
  fix <- make_count_fixture()
  fix$meta$sample_id[2] <- "S1"

  expect_error(
    collapse_replicates(fix$expr, fix$meta, method = "mean"),
    "Duplicate"
  )
})


# ──────────────────────────────────────────────────────────────────────────────
# Warnings
# ──────────────────────────────────────────────────────────────────────────────

test_that("warning for singleton replicate group", {
  fix <- make_count_fixture()
  expect_warning(
    collapse_replicates(fix$expr, fix$meta, method = "mean"),
    "only 1 sample"
  )
})


test_that("warning when sum used on non-integer data", {
  fix <- make_count_fixture()
  fix$expr["GeneA", "S1"] <- 10.5

  # Use selected_groups to avoid the singleton warning interfering
  expect_warning(
    collapse_replicates(fix$expr, fix$meta, method = "sum",
                        selected_groups = "Patient1"),
    "may not be integer"
  )
})


test_that("warning when mean used on count-like data", {
  fix <- make_count_fixture()
  fix$expr <- fix$expr * 100  # make counts large

  # Use selected_groups to avoid the singleton warning interfering
  expect_warning(
    collapse_replicates(fix$expr, fix$meta, method = "mean",
                        selected_groups = "Patient1"),
    "raw counts"
  )
})


# ──────────────────────────────────────────────────────────────────────────────
# NA handling
# ──────────────────────────────────────────────────────────────────────────────

test_that("NA values are handled with na_rm=TRUE", {
  fix <- make_count_fixture()
  fix$expr["GeneA", "S1"] <- NA

  result <- suppressWarnings(
    collapse_replicates(fix$expr, fix$meta, method = "mean", na_rm = TRUE)
  )

  # mean(NA, 12) with na.rm=TRUE = 12
  expect_equal(result$expr["GeneA", "Patient1_collapsed_mean"], 12)
})


test_that("NA values propagate with na_rm=FALSE", {
  fix <- make_count_fixture()
  fix$expr["GeneA", "S1"] <- NA

  result <- suppressWarnings(
    collapse_replicates(fix$expr, fix$meta, method = "mean", na_rm = FALSE)
  )

  # mean(NA, 12) with na.rm=FALSE = NA
  expect_true(is.na(result$expr["GeneA", "Patient1_collapsed_mean"]))
})


# ──────────────────────────────────────────────────────────────────────────────
# Output structure and synchronisation
# ──────────────────────────────────────────────────────────────────────────────

test_that("output matrix and metadata stay synchronized", {
  fix <- make_count_fixture()
  result <- suppressWarnings(
    collapse_replicates(fix$expr, fix$meta, method = "mean")
  )

  # Same number of columns and metadata rows
  expect_equal(ncol(result$expr), nrow(result$metadata))

  # Column names match sample_id column
  expect_equal(colnames(result$expr),
               as.character(result$metadata$sample_id))
})


test_that("feature row names are preserved exactly", {
  fix <- make_count_fixture()
  result <- suppressWarnings(
    collapse_replicates(fix$expr, fix$meta, method = "mean")
  )
  expect_equal(rownames(result$expr), c("GeneA", "GeneB"))
})


test_that("result is a list with expected components", {
  fix <- make_count_fixture()
  result <- suppressWarnings(
    collapse_replicates(fix$expr, fix$meta, method = "mean")
  )

  expect_type(result, "list")
  expect_true(all(c("expr", "metadata", "provenance", "warnings") %in%
                    names(result)))
  expect_true(is.matrix(result$expr))
  expect_true(is.data.frame(result$metadata))
  expect_true(is.character(result$warnings))
})


# ──────────────────────────────────────────────────────────────────────────────
# keep_original option
# ──────────────────────────────────────────────────────────────────────────────

test_that("keep_original=TRUE retains original samples", {
  fix <- make_count_fixture()
  result <- suppressWarnings(
    collapse_replicates(fix$expr, fix$meta, method = "mean",
                        keep_original = TRUE)
  )

  # Should have collapsed + uncollapsed + originals
  # Patient1: collapsed + S1_original + S2_original
  # Patient2: collapsed
  expect_true("S1_original" %in% colnames(result$expr))
  expect_true("S2_original" %in% colnames(result$expr))
  expect_equal(result$expr["GeneA", "S1_original"], 10)
})


# ──────────────────────────────────────────────────────────────────────────────
# Edge case: many genes, many samples
# ──────────────────────────────────────────────────────────────────────────────

test_that("works with larger data", {
  set.seed(42)
  n_genes <- 100
  n_samples <- 20

  expr <- matrix(rnorm(n_genes * n_samples), nrow = n_genes,
                 dimnames = list(paste0("G", 1:n_genes),
                                 paste0("S", 1:n_samples)))

  # Each pair of replicates shares the same condition
  meta <- data.frame(
    sample_id          = paste0("S", 1:n_samples),
    replicate_group_id = rep(paste0("Group", 1:10), each = 2),
    condition          = rep(c("A", "B", "A", "B", "A",
                               "B", "A", "B", "A", "B"), each = 2),
    stringsAsFactors   = FALSE
  )

  result <- collapse_replicates(expr, meta, method = "mean")

  expect_equal(ncol(result$expr), 10)
  expect_equal(nrow(result$expr), 100)
  expect_equal(nrow(result$metadata), 10)
})
