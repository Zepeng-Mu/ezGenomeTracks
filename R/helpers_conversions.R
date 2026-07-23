# ezGenomeTracks - Helper functions (split from helpers.R)
#' Convert a GRanges object to a tidy data frame
#'
#' This function converts a GRanges object to a tidy data frame suitable for use with ggplot2.
#' It preserves all metadata columns and adds columns for chromosome, start, end, and width.
#'
#' @param gr A GRanges object
#' @param keep.mcols Logical indicating whether to keep metadata columns (default: TRUE)
#' @return A tidy data frame
#' @export
#' @importFrom GenomicRanges seqnames start end width
#' @importFrom S4Vectors mcols
#' @importFrom methods is
#' @examples
#' \dontrun{
#' library(GenomicRanges)
#' gr <- GRanges(seqnames = c("chr1", "chr1", "chr2"),
#'               ranges = IRanges(start = c(1, 100, 200), end = c(50, 150, 250)),
#'               score = c(0.1, 0.5, 0.9))
#' gr_df <- granges_to_df(gr)
#' }
granges_to_df <- function(gr, keep.mcols = TRUE) {
  if (!methods::is(gr, "GRanges")) {
    stop("Input must be a GRanges object")
  }

  # Extract coordinates

  start_vals <- GenomicRanges::start(gr)
  end_vals <- GenomicRanges::end(gr)

  # Adjust start when start equals end
  start_vals <- ifelse(start_vals == end_vals, start_vals - 1, start_vals)

  # Create base data frame with coordinates
  df <- data.frame(
    seqnames = as.character(GenomicRanges::seqnames(gr)),
    start = start_vals,
    end = end_vals,
    width = end_vals - start_vals,
    strand = as.character(GenomicRanges::strand(gr))
  )

  # Add metadata columns if requested
  if (keep.mcols && ncol(S4Vectors::mcols(gr)) > 0) {
    df <- cbind(df, as.data.frame(S4Vectors::mcols(gr)))
  }

  return(df)
}

#' Convert a data frame to a GRanges object
#'
#' This function converts a data frame to a GRanges object. The data frame must have
#' columns for chromosome (seqnames), start, and end positions.
#'
#' @param df A data frame with at least seqnames, start, and end columns
#' @param seqnames Column name for chromosome (default: "seqnames")
#' @param start Column name for start position (default: "start")
#' @param end Column name for end position (default: "end")
#' @param strand Column name for strand (default: "strand")
#' @return A GRanges object
#' @export
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @examples
#' \dontrun{
#' df <- data.frame(
#'   seqnames = c("chr1", "chr1", "chr2"),
#'   start = c(1, 100, 200),
#'   end = c(50, 150, 250),
#'   score = c(0.1, 0.5, 0.9)
#' )
#' gr <- df_to_granges(df)
#' }
df_to_granges <- function(
  df,
  seqnames = "seqnames",
  start = "start",
  end = "end",
  strand = "strand"
) {
  if (!all(c(seqnames, start, end) %in% colnames(df))) {
    stop("Data frame must contain columns for seqnames, start, and end")
  }

  # Extract metadata columns (all columns except coordinate columns)
  coord_cols <- c(seqnames, start, end)
  if (strand %in% colnames(df)) {
    coord_cols <- c(coord_cols, strand)
    strand_values <- df[[strand]]
  } else {
    strand_values <- "*"
  }

  mcols_df <- df[, !colnames(df) %in% coord_cols, drop = FALSE]

  # Create GRanges object
  gr <- GenomicRanges::GRanges(
    seqnames = df[[seqnames]],
    ranges = IRanges::IRanges(start = df[[start]], end = df[[end]]),
    strand = strand_values
  )

  # Add metadata columns if any exist
  if (ncol(mcols_df) > 0) {
    S4Vectors::mcols(gr) <- mcols_df
  }

  return(gr)
}

#' Import genomic data from a file into a tidy data frame
#'
#' This function imports genomic data from various file formats (BED, bigWig, GFF, etc.)
#' using rtracklayer and converts it to a tidy data frame suitable for ggplot2.
#'
#' @param file Path to the genomic data file
#' @param format File format (default: NULL, auto-detected from file extension)
#' @param which GRanges object specifying the genomic region to import (default: NULL, import all)
#' @param n_bins Number of bins used for megadepth BigWig queries.
#'   Only used for BigWig inputs when `which` is provided. Default: 2000.
#' @param verbose Logical. If `TRUE`, print megadepth progress output for BigWig
#'   imports. If `FALSE` (default), suppress routine messages.
#' @return A tidy data frame
#' @export
#' @importFrom rtracklayer import
#' @examples
#' \dontrun{
#' # Import a BED file
#' peaks_df <- import_genomic_data("peaks.bed")
#'
#' # Import a specific region from a bigWig file
#' library(GenomicRanges)
#' region <- GRanges("chr1", IRanges(1000000, 2000000))
#' signal_df <- import_genomic_data("signal.bw", which = region)
#' }
import_genomic_data <- function(file, which = NULL, n_bins = 2000L, verbose = FALSE) {
  # Handle non-character input (data frames, etc.)
  if (!is.character(file)) {
    # For non-file inputs, just use rtracklayer/granges_to_df conversion
    if (is.data.frame(file)) {
      return(file)
    }
    # For GRanges or other objects
    gr <- rtracklayer::import(file, which = which)
    df <- granges_to_df(gr)
    return(df)
  }

  # Determine file type and dispatch to appropriate importer
  file_ext <- tolower(tools::file_ext(file))
  is_bigwig <- (file_ext == "bw" || file_ext == "bigwig")
  is_remote <- grepl("^https?://", file)

  # Use megadepth for BigWig files (local or remote) for improved performance
  if ((is_bigwig || is_remote) && !is.null(which)) {
    tryCatch(
      {
        gr <- import_bigwig_megadepth(
          file = file,
          which = which,
          op = "mean",
          n_bins = n_bins,
          verbose = verbose
        )
        # Ensure output has score column for consistency
        if (!("score" %in% names(S4Vectors::mcols(gr)))) {
          S4Vectors::mcols(gr)$score <- 0
        }
        df <- granges_to_df(gr)
        return(df)
      },
      error = function(e) {
        # Fallback to rtracklayer if megadepth fails
        warning(
          "megadepth processing failed: ",
          conditionMessage(e),
          "\nFalling back to rtracklayer for: ", file
        )
        gr <- rtracklayer::import(file, which = which)
        df <- granges_to_df(gr)
        return(df)
      }
    )
  }

  # For other file types or when no region is specified, use rtracklayer
  gr <- rtracklayer::import(file, which = which)

  # Convert to data frame
  df <- granges_to_df(gr)

  return(df)
}

# Internal helper: Detect and return the correct column name from a data frame
# given a list of candidate column names. Used for flexible column name handling
# in functions that accept multiple column naming conventions (e.g., GWAS vs GRanges style).
#
# @param data A data frame to search
# @param candidates Character vector of candidate column names
# @param param_name Name of the parameter (used in error message)
# @param required Logical indicating if the column is required (default: TRUE)
# @return The name of the detected column, or NULL if not found and not required