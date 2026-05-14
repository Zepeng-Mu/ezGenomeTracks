# ezGenomeTracks - Helper functions (split from helpers.R)
#' Filter gene labels based on priority and maximum count
#'
#' This helper function filters gene labels to reduce overlap in crowded plots.
#' It prioritizes genes based on specified criteria (e.g., length, name) and
#' returns the top N genes for labeling.
#'
#' @param label_data Data frame containing gene information to be labeled
#' @param max_labels Maximum number of labels to show. If NULL, all labels are kept.
#' @param label_priority Priority criterion for filtering. Options:
#'   - "length": Prioritize longer genes (calculated as end - start)
#'   - "name": Sort alphabetically by gene name
#'   - A column name in label_data: Sort by that column (descending for numeric, alphabetical for character)
#' @param start_col Name of the start coordinate column. Default: "start"
#' @param end_col Name of the end coordinate column. Default: "end"
#'
#' @return A filtered data frame with at most max_labels rows
#'
#' @keywords internal
#' @noRd
filter_labels <- function(
  label_data,
  max_labels = NULL,
  label_priority = "length",
  start_col = "start",
  end_col = "end"
) {
  # If max_labels is NULL or greater than nrow, return all labels
  if (is.null(max_labels) || max_labels >= nrow(label_data)) {
    return(label_data)
  }

  # Ensure max_labels is positive
  if (max_labels < 1) {
    warning("max_labels must be >= 1. Using 1.")
    max_labels <- 1
  }

  # Calculate priority scores
  if (label_priority == "length") {
    # Prioritize longer genes
    if (!start_col %in% names(label_data) || !end_col %in% names(label_data)) {
      stop(sprintf("Columns '%s' and '%s' required for priority='length'", start_col, end_col))
    }
    label_data$priority_score <- label_data[[end_col]] - label_data[[start_col]]
    label_data <- label_data[order(label_data$priority_score, decreasing = TRUE), ]
  } else if (label_priority == "name") {
    # Sort alphabetically (ascending)
    # Assume a column named "gene_name" or similar exists
    name_col <- if ("gene_name" %in% names(label_data)) "gene_name" else names(label_data)[1]
    label_data <- label_data[order(label_data[[name_col]]), ]
  } else if (label_priority %in% names(label_data)) {
    # Use specified column
    priority_col <- label_data[[label_priority]]
    if (is.numeric(priority_col)) {
      # Sort descending for numeric columns
      label_data <- label_data[order(priority_col, decreasing = TRUE), ]
    } else {
      # Sort ascending for character/factor columns
      label_data <- label_data[order(priority_col), ]
    }
  } else {
    warning(sprintf("Priority '%s' not recognized or not a column in data. Using default order.", label_priority))
  }

  # Return top N labels
  label_data[seq_len(min(max_labels, nrow(label_data))), ]
}


#' Auto-detect OrgDb package for gene symbol mapping
#'
#' This internal helper function attempts to auto-detect and load a suitable
#' OrgDb package for mapping ENTREZ gene IDs to gene symbols. It checks
#' common OrgDb packages for human, mouse, rat, and other model organisms.
#'
#' @return An OrgDb object if one is found, NULL otherwise
#' @keywords internal
#' @noRd
.detect_org_db <- function() {
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    return(NULL)
  }

  # Common OrgDb package names to check
  orgdb_packages <- c(
    "org.Hs.eg.db", # Human
    "org.Mm.eg.db", # Mouse
    "org.Rn.eg.db", # Rat
    "org.Dm.eg.db", # Drosophila
    "org.Ce.eg.db", # C. elegans
    "org.Sc.sgd.db", # Yeast
    "org.Dr.eg.db", # Zebrafish
    "org.At.tair.db" # Arabidopsis
  )

  for (pkg_name in orgdb_packages) {
    # Check if package is available and try to load the OrgDb object
    if (requireNamespace(pkg_name, quietly = TRUE)) {
      # Use package::package pattern to access the OrgDb object
      candidate <- tryCatch(
        getExportedValue(pkg_name, pkg_name),
        error = function(e) NULL
      )
      if (!is.null(candidate) && methods::is(candidate, "OrgDb")) {
        message("Auto-detected OrgDb: ", pkg_name)
        return(candidate)
      }
    }
  }

  return(NULL)
}


#' Convert gene name to genomic region
#'
#' This function looks up a gene by name in a TxDb object and returns a region
#' string representing the gene body with optional padding. This allows users
#' to specify genes by name instead of coordinates when using ez_* functions.
#'
#' @param gene_name Character string. Gene symbol to look up (e.g., "PTPRC", "TP53").
#' @param txdb A TxDb object containing gene annotations
#'   (e.g., TxDb.Hsapiens.UCSC.hg38.knownGene).
#' @param org_db Optional OrgDb object for mapping gene symbols to ENTREZ IDs.
#'   If NULL (default), auto-detects available OrgDb packages.
#' @param extend Numeric. Amount to extend the region beyond the gene body.
#'   Interpretation depends on `extend_type`. Default: 0.1 (10% of gene length).
#' @param extend_type Character. How to interpret `extend`:
#'   - "proportion": `extend` is a proportion of gene length (default)
#'   - "bp": `extend` is an absolute number of base pairs
#'
#' @return A character string in region format "chr:start-end" suitable for
#'   use with ez_* functions.
#'
#' @details
#' The function maps gene symbols to ENTREZ IDs using the OrgDb, then queries
#' the TxDb for gene coordinates. If multiple genes match (e.g., same symbol
#' on different chromosomes), a warning is issued and the first match
#' (sorted by chromosome, then start position) is used.
#'
#' The padding extends the region on both sides of the gene body. For example,
#' with `extend = 0.1` (default) and a gene of length 10kb, the region will
#' extend 1kb upstream and 1kb downstream (total region: 12kb).
#'
#' @export
#' @importFrom GenomicRanges start end seqnames width
#' @importFrom methods is
#' @examples
#' \dontrun{
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' library(org.Hs.eg.db)
#'
#' # Get region for PTPRC gene with 10% padding (default)
#' region <- gene_to_region("PTPRC", TxDb.Hsapiens.UCSC.hg38.knownGene)
#'
#' # Use in ez_coverage
#' ez_coverage(signal_data, region)
#'
#' # With 5kb fixed padding on each side
#' region <- gene_to_region("TP53", TxDb.Hsapiens.UCSC.hg38.knownGene,
#'                          extend = 5000, extend_type = "bp")
#' }
gene_to_region <- function(
  gene_name,
  txdb,
  org_db = NULL,
  extend = 0.1,
  extend_type = c("proportion", "bp")
) {
  extend_type <- match.arg(extend_type)

  # Validate inputs
  if (!is.character(gene_name) || length(gene_name) != 1 || nchar(gene_name) == 0) {
    stop("gene_name must be a single non-empty character string")
  }

  if (!methods::is(txdb, "TxDb")) {
    stop("txdb must be a TxDb object (e.g., TxDb.Hsapiens.UCSC.hg38.knownGene)")
  }

  if (!is.numeric(extend) || extend < 0) {
    stop("extend must be a non-negative number")
  }

  # Check required packages
  if (!requireNamespace("GenomicFeatures", quietly = TRUE)) {
    stop(
      "Package 'GenomicFeatures' is required for TxDb support. ",
      "Install it with: BiocManager::install('GenomicFeatures')"
    )
  }
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    stop(
      "Package 'AnnotationDbi' is required for gene symbol lookup. ",
      "Install it with: BiocManager::install('AnnotationDbi')"
    )
  }

  # Auto-detect OrgDb if not provided
  if (is.null(org_db)) {
    org_db <- .detect_org_db()
  }

  if (is.null(org_db)) {
    stop(
      "No OrgDb package found for gene symbol mapping. ",
      "Install an appropriate OrgDb package (e.g., org.Hs.eg.db for human) ",
      "or provide one via the org_db argument."
    )
  }

  # Map gene symbol to ENTREZ ID
  gene_info <- tryCatch(
    {
      AnnotationDbi::select(
        org_db,
        keys = gene_name,
        columns = c("ENTREZID", "SYMBOL"),
        keytype = "SYMBOL"
      )
    },
    error = function(e) {
      stop(sprintf("Failed to look up gene '%s': %s", gene_name, e$message))
    }
  )

  if (nrow(gene_info) == 0 || all(is.na(gene_info$ENTREZID))) {
    stop(sprintf("Gene '%s' not found in the OrgDb database", gene_name))
  }

  # Get unique ENTREZ IDs (filter out NAs)
  entrez_ids <- unique(gene_info$ENTREZID[!is.na(gene_info$ENTREZID)])

  # Get all genes from TxDb
  all_genes <- GenomicFeatures::genes(txdb)

  # Find matching genes
  matching_genes <- all_genes[names(all_genes) %in% entrez_ids]

  if (length(matching_genes) == 0) {
    stop(sprintf(
      "Gene '%s' (ENTREZ ID: %s) not found in the TxDb database",
      gene_name,
      paste(entrez_ids, collapse = ", ")
    ))
  }

  # Handle multiple matches
  if (length(matching_genes) > 1) {
    # Sort by chromosome, then by start position
    matching_df <- data.frame(
      entrez_id = names(matching_genes),
      seqnames = as.character(GenomicRanges::seqnames(matching_genes)),
      start = GenomicRanges::start(matching_genes),
      end = GenomicRanges::end(matching_genes),
      stringsAsFactors = FALSE
    )
    matching_df <- matching_df[order(matching_df$seqnames, matching_df$start), ]

    warning(sprintf(
      "Multiple genomic locations found for gene '%s' (%d matches). Using first match at %s:%d-%d",
      gene_name,
      length(matching_genes),
      matching_df$seqnames[1],
      matching_df$start[1],
      matching_df$end[1]
    ))

    # Use first match
    matching_genes <- matching_genes[matching_df$entrez_id[1]]
  }

  # Extract coordinates
  chr <- as.character(GenomicRanges::seqnames(matching_genes))
  gene_start <- GenomicRanges::start(matching_genes)
  gene_end <- GenomicRanges::end(matching_genes)
  gene_length <- GenomicRanges::width(matching_genes)

  # Calculate padding
  if (extend_type == "proportion") {
    padding <- round(gene_length * extend)
  } else {
    padding <- extend
  }

  # Calculate extended region (ensure start is not negative)
  region_start <- max(1, gene_start - padding)
  region_end <- gene_end + padding

  # Return region string
  sprintf("%s:%d-%d", chr, region_start, region_end)
}

