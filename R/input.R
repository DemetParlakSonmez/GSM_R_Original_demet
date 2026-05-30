#' @title Input readers for GSM expression and grouping data
#' @description Functions to read expression matrices, grouping files (2-column
#'   CSV or GMT format), and pre-transposed expression matrices.
#' @import data.table
#' @name input
NULL

#' Read an expression CSV in GSM format
#'
#' Reads a CSV where rows are samples, the first column contains class labels,
#' and remaining columns are gene expression values. Returns a transposed
#' expression matrix (genes x samples) and a named factor of class labels.
#'
#' @param file Path to the expression CSV file.
#' @param class_col Name or index of the class label column. Default `"class"`.
#' @param sep Column separator. Default `","`.
#' @return A list with components:
#'   \describe{
#'     \item{expr}{Numeric matrix with genes as rows and samples as columns.}
#'     \item{labels}{Named factor of class labels (names = sample IDs).}
#'   }
#' @export
read_expression <- function(file, class_col = "class", sep = ",", sample_id_col = NULL) {
  stopifnot(file.exists(file))

  # Use data.table::fread for speed on large files
  dt <- data.table::fread(file, sep = sep, header = TRUE, check.names = FALSE)

  first_col <- names(dt)[1]

  # Extract class labels
  if (is.character(class_col)) {
    if (!class_col %in% names(dt)) {
      stop("Class column '", class_col, "' not found in file. ",
           "Available columns: ", paste(head(names(dt), 5), collapse = ", "), "...")
    }
    label_vec <- dt[[class_col]]
    dt[, (class_col) := NULL]
  } else {
    label_vec <- dt[[class_col]]
    dt[, (class_col) := NULL, with = FALSE]
  }

  # Detect optional sample ID column
  if (is.null(sample_id_col) && first_col != class_col) {
    is_sample_col <- grepl("sample|id|name|patient|subject", first_col, ignore.case = TRUE) ||
      all(vapply(head(dt[[first_col]], 5), function(val) {
        (is.character(val) || is.factor(val)) && is.na(suppressWarnings(as.numeric(as.character(val))))
      }, logical(1)))
    
    if (is_sample_col) {
      sample_id_col <- first_col
    }
  }

  if (!is.null(sample_id_col)) {
    sample_ids <- as.character(dt[[sample_id_col]])
    dt[, (sample_id_col) := NULL]
  } else {
    sample_ids <- paste0("S", seq_len(nrow(dt)))
  }

  # Convert to matrix and transpose (genes x samples)
  mat <- as.matrix(dt)

  # Ensure numeric
  if (!is.numeric(mat)) {
    # Try coercing — some files quote numbers
    mat <- apply(mat, 2, as.numeric)
  }

  # Transpose: genes x samples
  expr <- t(mat)
  colnames(expr) <- sample_ids

  # Labels
  labels <- factor(label_vec)
  names(labels) <- sample_ids

  # Validate
  if (length(unique(labels)) < 2) {
    stop("Expression data must have at least 2 distinct class labels. Found: ",
         paste(unique(labels), collapse = ", "))
  }
  if (any(is.na(labels))) {
    stop("Class labels contain NA values.")
  }

  list(expr = expr, labels = labels)
}


#' Read a 2-column grouping file
#'
#' Reads a CSV/TSV with gene-group associations (e.g. gene-disease from
#' DisGeNET, gene-pathway from KEGG, gene-miRNA). Returns a named list where
#' each element is a character vector of gene symbols belonging to that group.
#'
#' @param file Path to the grouping file.
#' @param gene_col Name or index of the gene column. Default `1`.
#' @param group_col Name or index of the group column. Default `2`.
#' @param sep Column separator. Default `","`.
#' @return Named list of character vectors (group name -> gene symbols).
#' @export
read_groups <- function(file, gene_col = 1, group_col = 2, sep = ",") {
  stopifnot(file.exists(file))

  dt <- data.table::fread(file, sep = sep, header = TRUE, check.names = FALSE)

  if (ncol(dt) < 2) {
    stop("Grouping file must have at least 2 columns. Found: ", ncol(dt))
  }

  # Resolve column references
  gene_col_name <- if (is.numeric(gene_col)) names(dt)[gene_col] else gene_col
  group_col_name <- if (is.numeric(group_col)) names(dt)[group_col] else group_col

  if (!gene_col_name %in% names(dt)) {
    stop("Gene column '", gene_col_name, "' not found.")
  }
  if (!group_col_name %in% names(dt)) {
    stop("Group column '", group_col_name, "' not found.")
  }

  genes <- as.character(dt[[gene_col_name]])
  groups_vec <- as.character(dt[[group_col_name]])

  # Remove duplicates
  pairs <- unique(data.frame(gene = genes, group = groups_vec,
                             stringsAsFactors = FALSE))

  # Split into named list
  groups <- split(pairs$gene, pairs$group)

  # Sort genes within each group for reproducibility
  groups <- lapply(groups, sort)

  groups
}


#' Read a GMT (Gene Matrix Transposed) file
#'
#' Reads standard GMT format where each line is:
#' `group_name<TAB>description<TAB>gene1<TAB>gene2<TAB>...`
#'
#' @param file Path to the GMT file.
#' @return Named list of character vectors (group name -> gene symbols).
#' @export
read_gmt <- function(file) {
  stopifnot(file.exists(file))

  lines <- readLines(file, warn = FALSE)
  # Remove empty lines
  lines <- lines[nchar(trimws(lines)) > 0]

  groups <- list()
  for (line in lines) {
    fields <- strsplit(line, "\t")[[1]]
    if (length(fields) < 3) next
    name <- fields[1]
    # fields[2] is the description, skip it
    genes <- fields[3:length(fields)]
    genes <- genes[nchar(trimws(genes)) > 0]
    if (length(genes) > 0) {
      groups[[name]] <- sort(genes)
    }
  }

  groups
}


#' Read a pre-transposed expression matrix
#'
#' Reads a matrix file where genes are rows and samples are columns.
#' No class label column is expected.
#'
#' @param file Path to the file.
#' @param sep Column separator. Default `"\t"`.
#' @param row_names Column to use as row names. Default `1`.
#' @return Numeric matrix with genes as rows and samples as columns.
#' @export
read_expression_matrix <- function(file, sep = "\t", row_names = 1) {
  stopifnot(file.exists(file))

  dt <- data.table::fread(file, sep = sep, header = TRUE, check.names = FALSE)

  # Use first column as row names if specified
  rn <- NULL
  if (!is.null(row_names)) {
    rn <- as.character(dt[[row_names]])
    dt[, (row_names) := NULL, with = is.numeric(row_names)]
  }

  mat <- as.matrix(dt)
  if (!is.numeric(mat)) {
    mat <- apply(mat, 2, as.numeric)
  }

  if (!is.null(rn)) {
    rownames(mat) <- rn
  }

  mat
}
