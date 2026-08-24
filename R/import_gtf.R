#' Import and parse GTF/GFF files into standardized data frame
#'
#' This function reads GTF/GFF annotation files (e.g., from GENCODE, Ensembl, RefSeq)
#' and converts them to a standardized data frame format compatible with ez_* functions.
#'
#' The function automatically extracts gene identifiers, gene names, and transcript IDs
#' from the GTF attributes column, handling missing values gracefully.
#'
#' @param file Path to GTF/GFF file. Can be gzipped (.gtf.gz, .gff.gz) or uncompressed.
#' @param region Optional genomic region in format "chr:start-end" to filter features.
#'   If NULL (default), imports all features in the file. Example: "chr19:16000000-16100000"
#' @param features Character vector of feature types to retain (e.g., c("gene", "exon")).
#'   Default: c("gene", "exon", "transcript", "CDS", "start_codon", "stop_codon").
#'   Set to "all" to keep all features found in the file.
#' @param gene_id_col Column name/attribute to extract as gene_id. Default: "gene_id".
#'   This is typically "gene_id" in GENCODE/Ensembl GTF files.
#' @param gene_name_col Column name/attribute to extract as gene_name. Default: "gene_name".
#'   Falls back to gene_id if gene_name is not available in the annotation.
#' @param transcript_id_col Column name/attribute for transcript IDs. Default: "transcript_id".
#' @param verbose If TRUE, print import progress messages. Default: FALSE
#'
#' @return A data frame with standardized columns:
#'   \describe{
#'     \item{seqnames}{Chromosome name (e.g., "chr19", "1")}
#'     \item{start}{Start coordinate (1-based, as in GTF)}
#'     \item{end}{End coordinate (1-based, inclusive, as in GTF)}
#'     \item{width}{Feature width (end - start)}
#'     \item{strand}{Strand information: "+" (forward), "-" (reverse), or "*" (unknown)}
#'     \item{type}{Feature type: "gene", "exon", "transcript", "CDS", etc.}
#'     \item{gene_id}{Gene identifier (e.g., "ENSG00000000003")}
#'     \item{gene_name}{Gene symbol/name (e.g., "PTPRC"), fallback to gene_id if unavailable}
#'     \item{transcript_id}{Transcript identifier (e.g., "ENST00000000005"), may contain NA}
#'     \item{...}{Additional columns from GTF attributes (varies by source)}
#'   }
#'
#' @export
#' @importFrom rtracklayer import
#' @importFrom methods is
#' @importFrom GenomicRanges seqnames start end strand
#' @importFrom S4Vectors mcols
#'
#' @examples
#' \dontrun{
#' # Download GENCODE v47 (replace with actual URL/file):
#' # url <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/
#' #        gencode.v47.primary_assembly.annotation.gtf.gz"
#'
#' # Import entire GTF
#' gtf_df <- import_gtf("gencode.v47.annotation.gtf.gz")
#' head(gtf_df)
#' table(gtf_df$type)  # View feature type distribution
#'
#' # Import only genes and exons (reduces memory usage)
#' gtf_df <- import_gtf("gencode.v47.annotation.gtf.gz",
#'                      features = c("gene", "exon"))
#'
#' # Import region of interest (e.g., PTPRC locus)
#' gtf_df <- import_gtf("gencode.v47.annotation.gtf.gz",
#'                      region = "chr19:16000000-16100000")
#' subset(gtf_df, gene_name == "PTPRC")
#'
#' # Use imported GTF with ez_gene
#' library(ezGenomeTracks)
#' ez_gene(gtf_df, "chr19:16000000-16100000", label = "gene_name")
#' }
import_gtf <- function(
  file,
  region = NULL,
  features = c("gene", "exon", "transcript", "CDS", "start_codon", "stop_codon"),
  gene_id_col = "gene_id",
  gene_name_col = "gene_name",
  transcript_id_col = "transcript_id",
  verbose = FALSE
) {
  # Validate file exists
  if (!file.exists(file)) {
    stop("File not found: ", file)
  }

  # Parse region if provided
  if (!is.null(region)) {
    region_gr <- parse_region(region)
  } else {
    region_gr <- NULL
  }

  if (verbose) cat("Importing GTF/GFF: ", file, "\n", sep = "")

  # Import using rtracklayer
  # rtracklayer handles both .gtf and .gff, compressed and uncompressed
  tryCatch(
    {
      gene_gr <- rtracklayer::import(
        file,
        which = region_gr,
        format = "gtf"
      )
    },
    error = function(e) {
      # Try as GFF if GTF fails
      stop(
        "Failed to import GTF/GFF file. Ensure file is valid format. Error: ",
        conditionMessage(e)
      )
    }
  )

  if (length(gene_gr) == 0) {
    stop("No features found in GTF file. Check file format or region specification.")
  }

  if (verbose) cat("Loaded ", length(gene_gr), " features\n", sep = "")

  # Filter features by type
  if (!("all" %in% features)) {
    feature_mask <- as.character(gene_gr$type) %in% features
    gene_gr <- gene_gr[feature_mask]
    if (verbose) cat("Filtered to ", length(gene_gr), " features\n", sep = "")
  }

  # Convert GRanges to data frame
  gtf_df <- granges_to_df(gene_gr)

  # Extract and standardize GTF attributes
  # rtracklayer automatically extracts common attributes into mcols
  attributes_df <- as.data.frame(S4Vectors::mcols(gene_gr), stringsAsFactors = FALSE)

  # Extract gene_id
  gene_id_vals <- .extract_gtf_attribute(
    attributes_df,
    gene_id_col,
    default = names(gene_gr)  # Fallback to GRanges names
  )

  # Extract gene_name (fallback to gene_id if not available)
  gene_name_vals <- .extract_gtf_attribute(
    attributes_df,
    gene_name_col,
    default = gene_id_vals
  )

  # Extract transcript_id (allow NA)
  transcript_id_vals <- .extract_gtf_attribute(
    attributes_df,
    transcript_id_col,
    default = NA_character_
  )

  # Build standardized output data frame
  gtf_df$gene_id <- gene_id_vals
  gtf_df$gene_name <- gene_name_vals
  gtf_df$transcript_id <- transcript_id_vals

  # Reorder columns to put key columns first
  key_cols <- c("seqnames", "start", "end", "width", "strand", "type",
                "gene_id", "gene_name", "transcript_id")
  key_cols <- key_cols[key_cols %in% names(gtf_df)]
  other_cols <- setdiff(names(gtf_df), key_cols)
  gtf_df <- gtf_df[, c(key_cols, other_cols)]

  if (verbose) {
    cat("Success! Data frame with:\n")
    cat("  Rows: ", nrow(gtf_df), "\n", sep = "")
    cat("  Columns: ", ncol(gtf_df), "\n", sep = "")
    cat("  Feature types: ", paste(unique(gtf_df$type), collapse = ", "), "\n", sep = "")
  }

  return(gtf_df)
}
