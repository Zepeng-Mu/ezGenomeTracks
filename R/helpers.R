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

  # Create base data frame with coordinates
  df <- data.frame(
    seqnames = as.character(GenomicRanges::seqnames(gr)),
    start = GenomicRanges::start(gr),
    end = GenomicRanges::end(gr),
    width = GenomicRanges::width(gr),
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
df_to_granges <- function(df, seqnames = "seqnames", start = "start", end = "end", strand = "strand") {
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
import_genomic_data <- function(file, format = NULL, which = NULL) {
  # Import data using rtracklayer
  gr <- rtracklayer::import(file, format = format, which = which)

  # Convert to data frame
  df <- granges_to_df(gr)

  return(df)
}

#' Extract signal data for a single input element
#'
#' This function extracts genomic signal data for a specified region from either
#' a data frame or a file. It filters the data to only include features within
#' the specified region and optionally adds a track name.
#'
#' @param input Either a data frame with genomic coordinates or a character string
#'   specifying the path to a genomic data file (BED, bigWig, GFF, etc.)
#' @param region A genomic region string in the format "chr:start-end" or a GRanges object
#' @param name Optional name to assign to the track (default: NULL)
#' @return A data frame containing the filtered genomic data with an optional name column
#' @export
#' @importFrom dplyr filter mutate
#' @examples
#' \dontrun{
#' # Extract data from a data frame
#' df <- data.frame(
#'   seqnames = c("chr1", "chr1", "chr2"),
#'   start = c(1, 100, 200),
#'   end = c(50, 150, 250),
#'   score = c(0.1, 0.5, 0.9)
#' )
#' region_data <- get_single_signal(df, "chr1:50-150", name = "track1")
#'
#' # Extract data from a file
#' file_data <- get_single_signal("peaks.bed", "chr1:1000000-2000000", name = "peaks")
#' }
get_single_signal <- function(input, region, name = NULL) {
  region_gr <- parse_region(region = region)
  if (is(input, "data.frame")) {
    # Single track, data frame
    track_data <- input %>%
      dplyr::filter(seqnames == as.character(region_gr@seqnames),
                    start >= region_gr@start,
                    end <= region_gr@end) %>%
      dplyr::mutate(name = name)
  } else if (is(input, "character")) {
    # Single track, file name
    track_data <- import_genomic_data(input, region) %>%
      dplyr::mutate(name = name)
  }

  return(track_data)
}



#' Parse a genomic region string into a GRanges object
#'
#' This function parses a genomic region string in the format "chr:start-end" into a GRanges object.
#'
#' @param region A string specifying a genomic region (e.g., "chr1:1000000-2000000")
#' @return A GRanges object representing the specified region
#' @export
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @examples
#' \dontrun{
#' region_gr <- parse_region("chr1:1000000-2000000")
#' }
parse_region <- function(region) {
  if (!is.character(region) || length(region) != 1) {
    stop("Region must be a single character string")
  }

  # Parse the region string using regular expression
  # Match pattern: chromosome followed by any separator (: _ -) then start and end positions
  matches <- regexpr("^(.+?)[:_-](\\d+)[:_-](\\d+)$", region, perl = TRUE)
  if (matches == -1) {
    stop("Region must be in the format 'chr*start*end' where * can be :, _, or -")
  }

  # Extract the matched groups
  parsed <- regmatches(region, matches)
  groups <- stringr::str_match(parsed, "^(.+?)[:_-](\\d+)[:_-](\\d+)$")

  chr <- groups[, 2]
  start <- as.numeric(groups[, 3])
  end <- as.numeric(groups[, 4])

  if (is.na(start) || is.na(end)) {
    stop("Start and end positions must be numeric")
  }

  # Create GRanges object
  gr <- GenomicRanges::GRanges(
    seqnames = chr,
    ranges = IRanges::IRanges(start = start, end = end)
  )

  return(gr)
}