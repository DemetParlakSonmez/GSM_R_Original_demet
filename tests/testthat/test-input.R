# Tests for input reader functions

test_that("read_expression reads synthetic test data correctly", {
  f <- test_data_path("test_expression.csv")
  result <- read_expression(f)

  expect_type(result, "list")
  expect_true("expr" %in% names(result))
  expect_true("labels" %in% names(result))

  # Dimensions: 50 genes x 20 samples

  expect_equal(nrow(result$expr), 50)
  expect_equal(ncol(result$expr), 20)

  # Labels
  expect_s3_class(result$labels, "factor")
  expect_equal(length(result$labels), 20)
  expect_true(all(c("pos", "neg") %in% levels(result$labels)))
  expect_equal(sum(result$labels == "pos"), 10)
  expect_equal(sum(result$labels == "neg"), 10)

  # Numeric content
  expect_true(is.numeric(result$expr))
  expect_false(any(is.na(result$labels)))
})


test_that("read_expression errors on non-existent file", {
  expect_error(read_expression("nonexistent.csv"))
})


test_that("read_expression errors on missing class column", {
  # Create a temporary file without a 'class' column
  tmp <- tempfile(fileext = ".csv")
  write.csv(data.frame(x = 1:5, y = 6:10), tmp, row.names = FALSE)
  expect_error(read_expression(tmp, class_col = "class"), "Class column")
  unlink(tmp)
})


test_that("read_groups reads 2-column grouping file correctly", {
  f <- test_data_path("test_groups.csv")
  groups <- read_groups(f)

  expect_type(groups, "list")
  expect_true(length(groups) >= 4)  # At least 4 groups
  expect_true("GroupA_signal" %in% names(groups))
  expect_true("GroupB_noise" %in% names(groups))

  # GroupA has 10 unique genes (duplicates removed)
  expect_equal(length(groups[["GroupA_signal"]]), 10)

  # All values are character
  expect_true(all(vapply(groups, is.character, logical(1))))
})


test_that("read_groups errors on non-existent file", {
  expect_error(read_groups("nonexistent.csv"))
})


test_that("read_gmt reads GMT file correctly", {
  f <- test_data_path("test_groups.gmt")
  groups <- read_gmt(f)

  expect_type(groups, "list")
  expect_equal(length(groups), 5)
  expect_true("GroupA_signal" %in% names(groups))

  # GroupA has 10 genes
  expect_equal(length(groups[["GroupA_signal"]]), 10)
  # GroupE has 15 genes
  expect_equal(length(groups[["GroupE_noise"]]), 15)
})


test_that("read_gmt and read_groups produce consistent results", {
  f_gmt <- test_data_path("test_groups.gmt")
  f_csv <- test_data_path("test_groups.csv")

  groups_gmt <- read_gmt(f_gmt)
  groups_csv <- read_groups(f_csv)

  # Same group names
  expect_setequal(names(groups_gmt), names(groups_csv))

  # Same genes in each group (both are sorted)
  for (nm in names(groups_gmt)) {
    expect_equal(groups_gmt[[nm]], groups_csv[[nm]],
                 info = paste("Group:", nm))
  }
})


test_that("read_expression_matrix reads a transposed matrix", {
  # Create a temporary transposed matrix file
  tmp <- tempfile(fileext = ".tsv")
  mat <- matrix(rnorm(30), nrow = 6, ncol = 5)
  rownames(mat) <- paste0("GENE", 1:6)
  colnames(mat) <- paste0("S", 1:5)
  write.table(mat, tmp, sep = "\t", quote = FALSE, col.names = NA)

  result <- read_expression_matrix(tmp, sep = "\t")
  expect_true(is.matrix(result))
  expect_true(is.numeric(result))
  expect_equal(nrow(result), 6)
  expect_equal(ncol(result), 5)
  unlink(tmp)
})
