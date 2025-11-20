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
import_genomic_data <- function(file, which = NULL) {
  # Import data using rtracklayer
  gr <- rtracklayer::import(file, which = which)

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
#' @importFrom dplyr filter mutate bind_rows
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
    track_data <- import_genomic_data(input, which = region_gr) %>%
      dplyr::mutate(name = name)
  }

  return(track_data)
}

#' Process signal input into standardized data frame
#'
#' This function converts any input type (data.frame, character vector, or list)
#' into a standardized data frame with consistent columns for signal visualization.
#'
#' @param input Either a data frame, character vector of file paths, or named list
#' @param region A genomic region string in the format "chr:start-end"
#' @param track_labels Optional vector of track labels (used for character vector input)
#' @return A data frame with standardized columns: seqnames, start, end, score,
#'   and optionally track and group columns
#' @export
#' @importFrom dplyr bind_rows mutate filter
#' @examples
#' \dontrun{
#' # Data frame input
#' df <- data.frame(seqnames = "chr1", start = 1:100, end = 1:100, score = rnorm(100))
#' process_signal_input(df, "chr1:1-100")
#'
#' # Character vector input
#' files <- c("file1.bw", "file2.bw")
#' process_signal_input(files, "chr1:1-100", track_labels = c("Sample1", "Sample2"))
#'
#' # List input
#' data_list <- list("Track1" = df, "Track2" = files)
#' process_signal_input(data_list, "chr1:1-100")
#' }
process_signal_input <- function(input, region, track_labels = NULL) {
  if (is.data.frame(input)) {
    # Case 1: Data frame input
    # Validate required columns
    required_cols <- c("seqnames", "start", "end", "score")
    if (!all(required_cols %in% colnames(input))) {
      stop("Data frame must contain columns: ", paste(required_cols, collapse = ", "))
    }

    # Filter by region if needed
    region_gr <- parse_region(region)
    filtered_data <- input %>%
      dplyr::filter(seqnames == as.character(region_gr@seqnames),
                    start >= region_gr@start,
                    end <= region_gr@end)

    return(filtered_data)

  } else if (is.character(input)) {
    # Case 2: Character vector input (file paths)
    if (length(input) == 1) {
      # Single file
      track_name <- ifelse(is.null(track_labels), "Track 1", track_labels[1])
      return(get_single_signal(input, region, name = track_name))
    } else {
      # Multiple files
      track_data_list <- list()
      for (i in seq_along(input)) {
        track_name <- ifelse(is.null(track_labels), paste0("Track ", i), track_labels[i])
        track_data <- get_single_signal(input[i], region, name = track_name)
        track_data$group <- track_name
        track_data_list[[i]] <- track_data
      }
      return(dplyr::bind_rows(track_data_list))
    }

  } else if (is.list(input)) {
    # Case 3: List input
    if (is.null(names(input)) && is.null(track_labels)) {
      names(input) <- paste0("Track ", seq_along(input))
    } else if (is.null(names(input)) && !is.null(track_labels)) {
      names(input) <- track_labels
    }

    track_data_list <- list()
    for (i in seq_along(input)) {
      track_name <- names(input)[i]
      track_element <- input[[i]]

      if (is.data.frame(track_element)) {
        # Data frame element
        processed_data <- process_signal_input(track_element, region)
        processed_data$track <- track_name
        track_data_list[[i]] <- processed_data

      } else if (is.character(track_element)) {
        # Character vector element (multiple files for this track)
        if (length(track_element) == 1) {
          # Single file
          processed_data <- get_single_signal(track_element, region, name = track_name)
        } else {
          # Multiple files within this track
          file_data_list <- list()
          for (j in seq_along(track_element)) {
            file_data <- get_single_signal(track_element[j], region, name = paste0(track_name, "_", j))
            file_data$track <- track_name
            file_data$group <- paste0(track_name, "_", j)
            file_data_list[[j]] <- file_data
          }
          processed_data <- dplyr::bind_rows(file_data_list)
        }
        processed_data$track <- track_name
        track_data_list[[i]] <- processed_data

      } else {
        stop("List elements must be data frames or character vectors")
      }
    }

    names(track_data_list) <- names(input)
    return(dplyr::bind_rows(track_data_list))

  } else {
    stop("Input must be a data frame, character vector, or named list")
  }
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

#' Calculate y-axis limits for link tracks
#'
#' This function calculates appropriate y-axis limits for link/arc tracks based on
#' the maximum genomic distance span and height factor. This ensures curves are not
#' clipped and provides consistent spacing for multi-track plots.
#'
#' @param data A data frame with link/interaction data containing start1 and start2 columns
#' @param height_factor Height of curves as proportion of genomic distance span
#' @param direction Direction of curves: "down" (negative y) or "up" (positive y)
#' @return A numeric vector of length 2 with y-axis limits c(ymin, ymax)
#' @export
#' @examples
#' \dontrun{
#' df <- data.frame(
#'   start1 = c(1000, 2000, 3000),
#'   start2 = c(5000, 6000, 7000)
#' )
#' ylim <- calculate_link_ylim(df, height_factor = 0.15, direction = "down")
#' }
calculate_link_ylim <- function(data, height_factor = 0.15, direction = "down") {
  # Calculate maximum span across all links
  if (!all(c("start1", "start2") %in% colnames(data))) {
    stop("Data must contain 'start1' and 'start2' columns")
  }
  
  max_span <- max(abs(data$start2 - data$start1), na.rm = TRUE)
  max_height <- max_span * height_factor
  
  # Add 20% padding for visual spacing
  padding <- max_height * 0.2
  
  # Set limits based on direction
  # Curves start at y=0 and extend up or down
  if (direction == "down") {
    # Curves extend downward (negative y)
    # Y-axis needs to go from negative (bottom) to slightly positive (top)
    return(c(-(max_height + padding), padding))
  } else {
    # Curves extend upward (positive y)
    # Y-axis needs to go from slightly negative (bottom) to positive (top)
    return(c(-padding, max_height + padding))
  }
}
