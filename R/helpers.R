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
    track_data <- input |>
      dplyr::filter(
        seqnames == as.character(region_gr@seqnames),
        start >= region_gr@start,
        end <= region_gr@end
      ) |>
      dplyr::mutate(name = name)
  } else if (is(input, "character")) {
    # Single track, file name
    track_data <- import_genomic_data(input, which = region_gr) |>
      dplyr::mutate(name = name)
  }

  return(track_data)
}

#' Process signal input into standardized data frame
#'
#' This function converts any input type (GRanges, data.frame, character vector, or list)
#' into a standardized data frame with consistent columns for signal visualization.
#'
#' @param input A GRanges object, data frame, character vector of file paths, or named list.
#'   GRanges objects are automatically converted to data frames.
#' @param region A genomic region string in the format "chr:start-end"
#' @param track_labels Optional vector of track labels (used for character vector input)
#' @return A data frame with standardized columns: seqnames, start, end, score,
#'   and optionally track and group columns
#' @export
#' @importFrom dplyr bind_rows mutate filter
#' @importFrom methods is
#' @examples
#' \dontrun{
#' # GRanges input
#' library(GenomicRanges)
#' gr <- GRanges(seqnames = "chr1", ranges = IRanges(1:100, 1:100), score = rnorm(100))
#' process_signal_input(gr, "chr1:1-100")
#'
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
  if (methods::is(input, "GRanges")) {
    # Case 0: GRanges input - convert to data frame and process
    input <- granges_to_df(input)
  }

  if (is.data.frame(input)) {
    # Case 1: Data frame input
    # Validate required columns
    required_cols <- c("seqnames", "start", "end", "score")
    if (!all(required_cols %in% colnames(input))) {
      stop(
        "Data frame must contain columns: ",
        paste(required_cols, collapse = ", ")
      )
    }

    # Filter by region if needed
    region_gr <- parse_region(region)
    filtered_data <- input |>
      dplyr::filter(
        seqnames == as.character(region_gr@seqnames),
        start >= region_gr@start,
        end <= region_gr@end
      )

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
        track_name <- ifelse(
          is.null(track_labels),
          paste0("Track ", i),
          track_labels[i]
        )
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
          processed_data <- get_single_signal(
            track_element,
            region,
            name = track_name
          )
        } else {
          # Multiple files within this track
          file_data_list <- list()
          for (j in seq_along(track_element)) {
            file_data <- get_single_signal(
              track_element[j],
              region,
              name = paste0(track_name, "_", j)
            )
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

#' Process Manhattan plot input into standardized data frame
#'
#' This function converts any input type (data.frame or named list) into a
#' standardized data frame with consistent columns for Manhattan plot visualization.
#' Supports both GWAS-style (CHR, BP, P) and GRanges-style (seqnames, start, pvalue)
#' column naming conventions with auto-detection.
#'
#' @param input Either a data frame or named list of data frames
#' @param chr Column name for chromosome. If NULL, auto-detects from common names
#'   (CHR, chr, seqnames, chrom, chromosome). Default: NULL
#' @param bp Column name for base pair position. If NULL, auto-detects from common names
#'   (BP, bp, start, pos, position, POS). Default: NULL
#' @param p Column name for p-value. If NULL, auto-detects from common names
#'   (P, p, pvalue, p.value, pval). Default: NULL
#' @param snp Column name for SNP identifier. If NULL, auto-detects from common names
#'   (SNP, snp, rsid, id, variant_id, marker). Default: NULL (optional column)
#' @param track_labels Optional vector of track labels (used for unnamed list input)
#' @return A data frame with standardized columns and optionally track column
#' @export
#' @importFrom dplyr bind_rows mutate
#' @examples
#' # Data frame input with GWAS-style columns
#' df <- data.frame(CHR = 1, BP = 1:100, P = runif(100), SNP = paste0("rs", 1:100))
#' process_manhattan_input(df)
#'
#' # Data frame input with GRanges-style columns
#' df2 <- data.frame(seqnames = "chr1", start = 1:100, pvalue = runif(100))
#' process_manhattan_input(df2)
#'
#' # List input
#' data_list <- list("GWAS1" = df1, "GWAS2" = df2)
#' process_manhattan_input(data_list)
process_manhattan_input <- function(
  input,
  chr = NULL,
  bp = NULL,
  p = NULL,
  snp = NULL,
  track_labels = NULL
) {
  # Helper function to auto-detect column names
  detect_column <- function(data, candidates, param_name, required = TRUE) {
    for (col in candidates) {
      if (col %in% colnames(data)) return(col)
    }
    if (required) {
      stop(paste0(
        "Could not find ",
        param_name,
        " column. Expected one of: ",
        paste(candidates, collapse = ", ")
      ))
    }
    return(NULL)
  }

  # Define candidate column names
  chr_candidates <- c("CHR", "chr", "seqnames", "chrom", "chromosome")
  bp_candidates <- c("BP", "bp", "start", "pos", "position", "POS")
  p_candidates <- c("P", "p", "pvalue", "p.value", "pval", "P.value")
  snp_candidates <- c("SNP", "snp", "rsid", "id", "variant_id", "marker")

  if (is.data.frame(input)) {
    # Case 1: Data frame input
    # Auto-detect column names if not provided
    if (is.null(chr)) {
      chr <- detect_column(input, chr_candidates, "chromosome")
    }
    if (is.null(bp)) {
      bp <- detect_column(input, bp_candidates, "position")
    }
    if (is.null(p)) {
      p <- detect_column(input, p_candidates, "p-value")
    }
    if (is.null(snp)) {
      snp <- detect_column(input, snp_candidates, "SNP", required = FALSE)
    }

    return(input)
  } else if (is.list(input)) {
    # Case 2: List input (named or unnamed)
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
        # Auto-detect column names for each track element if not provided
        local_chr <- if (is.null(chr)) {
          detect_column(track_element, chr_candidates, "chromosome")
        } else {
          chr
        }
        local_bp <- if (is.null(bp)) {
          detect_column(track_element, bp_candidates, "position")
        } else {
          bp
        }
        local_p <- if (is.null(p)) {
          detect_column(track_element, p_candidates, "p-value")
        } else {
          p
        }

        # Add track column
        track_element$track <- track_name
        track_data_list[[i]] <- track_element
      } else {
        stop("List elements must be data frames")
      }
    }

    names(track_data_list) <- names(input)
    return(dplyr::bind_rows(track_data_list))
  } else {
    stop("Input must be a data frame or named list of data frames")
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

  # Remove commas (thousand separators) from the region string
  region <- gsub(",", "", region)

  # Parse the region string using regular expression
  # Match pattern: chromosome followed by any separator (: _ -) then start and end positions
  matches <- regexpr("^(.+?)[:_-](\\d+)[:_-](\\d+)$", region, perl = TRUE)
  if (matches == -1) {
    stop(
      "Region must be in the format 'chr*start*end' where * can be `:`, `_`, or `-`"
    )
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

#' Process sashimi plot data for visualization
#'
#' This function processes coverage and junction data for sashimi plot visualization.
#' It calculates y-positions for arc endpoints based on direction: "up" arcs start
#' at coverage height, "down" arcs start at zero, and "alternate" assigns alternating
#' directions to junctions sorted by start position.
#'
#' @param coverage_data Coverage input: a data frame with seqnames, start, end, score columns,
#'   or a file path to a bedGraph/bigWig file
#' @param junction_data Junction input: a data frame with seqnames, start, end, score columns,
#'   or a file path to a BED file. start/end represent the splice junction coordinates.
#' @param region A genomic region string in the format "chr:start-end"
#' @param junction_direction Direction of junction arcs: "alternate" (default), "up", or "down"
#' @return A list with two data frames: coverage_df and junction_df.
#'   junction_df has additional columns: y_start, y_end (arc endpoint y-values),
#'   arc_direction ("up" or "down" for each junction)
#' @export
#' @importFrom dplyr filter mutate arrange row_number
#' @examples
#' \dontrun{
#' coverage <- data.frame(
#'   seqnames = "chr1", start = 1:100, end = 2:101, score = runif(100, 0, 10)
#' )
#' junctions <- data.frame(
#'   seqnames = "chr1", start = c(20, 40, 60), end = c(35, 55, 80), score = c(5, 10, 3)
#' )
#' result <- process_sashimi_data(coverage, junctions, "chr1:1-100")
#' }
process_sashimi_data <- function(
  coverage_data,
  junction_data,
  region,
  junction_direction = c("alternate", "up", "down")
) {
  junction_direction <- match.arg(junction_direction)
  region_gr <- parse_region(region)

  # Process coverage data
  if (is.character(coverage_data) && length(coverage_data) == 1) {
    coverage_df <- import_genomic_data(coverage_data, which = region_gr)
  } else if (is.data.frame(coverage_data)) {
    required_cols <- c("seqnames", "start", "end", "score")
    if (!all(required_cols %in% colnames(coverage_data))) {
      stop(
        "Coverage data frame must contain columns: ",
        paste(required_cols, collapse = ", ")
      )
    }
    coverage_df <- coverage_data |>
      dplyr::filter(
        seqnames == as.character(region_gr@seqnames),
        end >= GenomicRanges::start(region_gr),
        start <= GenomicRanges::end(region_gr)
      )
  } else {
    stop("coverage_data must be a data frame or file path")
  }

  # Process junction data
  if (is.character(junction_data) && length(junction_data) == 1) {
    junction_df <- import_genomic_data(junction_data, which = region_gr)
  } else if (is.data.frame(junction_data)) {
    required_cols <- c("seqnames", "start", "end", "score")
    if (!all(required_cols %in% colnames(junction_data))) {
      stop(
        "Junction data frame must contain columns: ",
        paste(required_cols, collapse = ", ")
      )
    }
    junction_df <- junction_data |>
      dplyr::filter(
        seqnames == as.character(region_gr@seqnames),
        end >= GenomicRanges::start(region_gr),
        start <= GenomicRanges::end(region_gr)
      )
  } else {
    stop("junction_data must be a data frame or file path")
  }

  # Handle empty junction data

  if (nrow(junction_df) == 0) {
    junction_df$y_start <- numeric(0)
    junction_df$y_end <- numeric(0)
    junction_df$arc_direction <- character(0)
    junction_df$start1 <- numeric(0)
    junction_df$end1 <- numeric(0)
    junction_df$start2 <- numeric(0)
    junction_df$end2 <- numeric(0)
    return(list(coverage_df = coverage_df, junction_df = junction_df))
  }

  # Sort junctions by start position and assign direction
  junction_df <- junction_df |>
    dplyr::arrange(start) |>
    dplyr::mutate(
      arc_direction = if (junction_direction == "alternate") {
        ifelse(dplyr::row_number() %% 2 == 1, "up", "down")
      } else {
        junction_direction
      }
    )

  # Helper function to get coverage value at a position
  get_coverage_at_pos <- function(pos, cov_df) {
    # Find overlapping coverage bin
    idx <- which(cov_df$start <= pos & cov_df$end >= pos)
    if (length(idx) == 0) {
      return(0)
    }
    # Return the score of the first matching bin (or max if multiple)
    return(max(cov_df$score[idx], na.rm = TRUE))
  }

  # Calculate y-positions for arc endpoints
  junction_df <- junction_df |>
    dplyr::mutate(
      y_start = ifelse(
        arc_direction == "up",
        sapply(start, get_coverage_at_pos, cov_df = coverage_df),
        0
      ),
      y_end = ifelse(
        arc_direction == "up",
        sapply(end, get_coverage_at_pos, cov_df = coverage_df),
        0
      ),
      # Add start1/end1/start2/end2 columns for compatibility with geom_link
      start1 = start,
      end1 = start,
      start2 = end,
      end2 = end
    )

  return(list(coverage_df = coverage_df, junction_df = junction_df))
}

#' Process sashimi plot input into standardized data frames
#'
#' This function converts coverage and junction inputs of various types into
#' standardized data frames for sashimi plot visualization. It supports the same
#' input types as \code{\link{process_signal_input}} for coverage data (GRanges,
#' data.frame, character vector of file paths, or named list), and processes
#' junction data in parallel.
#'
#' @param coverage_data A GRanges object, data frame, character vector of file paths,
#'   or named list. Same input types supported as \code{\link{ez_coverage}}.
#' @param junction_data Junction input that must match the structure of coverage_data:
#'   - If coverage_data is a single source (GRanges, data.frame, or single file path),
#'     junction_data must also be a single source.
#'   - If coverage_data is a named list, junction_data must be a named list with
#'     the same names.
#' @param region A genomic region string in the format "chr:start-end"
#' @param track_labels Optional vector of track labels (used for character vector input)
#' @param junction_direction Direction of junction arcs: "alternate" (default), "up", or "down"
#' @return A list with two data frames: coverage_df and junction_df, both with
#'   a 'track' column when input is a list. junction_df has additional columns:
#'   y_start, y_end, arc_direction, start1, start2
#' @export
#' @importFrom dplyr bind_rows filter mutate arrange row_number
#' @importFrom methods is
#' @examples
#' \dontrun{
#' # Single data frame input
#' coverage <- data.frame(seqnames = "chr1", start = 1:100, end = 2:101, score = runif(100))
#' junctions <- data.frame(seqnames = "chr1", start = c(20, 60), end = c(40, 80), score = c(5, 10))
#' result <- process_sashimi_input(coverage, junctions, "chr1:1-100")
#'
#' # Named list input (multiple tracks)
#' cov_list <- list("Sample1" = cov_df1, "Sample2" = cov_df2)
#' junc_list <- list("Sample1" = junc_df1, "Sample2" = junc_df2)
#' result <- process_sashimi_input(cov_list, junc_list, "chr1:1-100")
#' }
process_sashimi_input <- function(
  coverage_data,
  junction_data,
  region,
  track_labels = NULL,
  junction_direction = c("alternate", "up", "down")
) {
  junction_direction <- match.arg(junction_direction)
  region_gr <- parse_region(region)

  # Helper function to get coverage value at a position
  get_coverage_at_pos <- function(pos, cov_df) {
    idx <- which(cov_df$start <= pos & cov_df$end >= pos)
    if (length(idx) == 0) return(0)
    return(max(cov_df$score[idx], na.rm = TRUE))
  }

  # Helper function to process a single junction source
  process_single_junction <- function(junc_input, cov_df, track_name = NULL) {
    # Import/filter junction data
    if (methods::is(junc_input, "GRanges")) {
      junc_input <- granges_to_df(junc_input)
    }

    if (is.character(junc_input) && length(junc_input) == 1) {
      junction_df <- import_genomic_data(junc_input, which = region_gr)
    } else if (is.data.frame(junc_input)) {
      required_cols <- c("seqnames", "start", "end", "score")
      if (!all(required_cols %in% colnames(junc_input))) {
        stop(
          "Junction data frame must contain columns: ",
          paste(required_cols, collapse = ", ")
        )
      }
      junction_df <- junc_input |>
        dplyr::filter(
          seqnames == as.character(region_gr@seqnames),
          end >= GenomicRanges::start(region_gr),
          start <= GenomicRanges::end(region_gr)
        )
    } else {
      stop("junction_data element must be a GRanges, data frame, or file path")
    }

    # Handle empty junction data
    if (nrow(junction_df) == 0) {
      junction_df$y_start <- numeric(0)
      junction_df$y_end <- numeric(0)
      junction_df$arc_direction <- character(0)
      junction_df$start1 <- numeric(0)
      junction_df$end1 <- numeric(0)
      junction_df$start2 <- numeric(0)
      junction_df$end2 <- numeric(0)
      if (!is.null(track_name)) junction_df$track <- character(0)
      return(junction_df)
    }

    # Sort junctions and assign arc direction
    junction_df <- junction_df |>
      dplyr::arrange(start) |>
      dplyr::mutate(
        arc_direction = if (junction_direction == "alternate") {
          ifelse(dplyr::row_number() %% 2 == 1, "up", "down")
        } else {
          junction_direction
        }
      )

    # Calculate y-positions based on coverage
    junction_df <- junction_df |>
      dplyr::mutate(
        y_start = ifelse(
          arc_direction == "up",
          sapply(start, get_coverage_at_pos, cov_df = cov_df),
          0
        ),
        y_end = ifelse(
          arc_direction == "up",
          sapply(end, get_coverage_at_pos, cov_df = cov_df),
          0
        ),
        start1 = start,
        end1 = start,
        start2 = end,
        end2 = end
      )

    if (!is.null(track_name)) {
      junction_df$track <- track_name
    }

    return(junction_df)
  }

  # Determine input type and process accordingly
  coverage_is_list <- is.list(coverage_data) && !is.data.frame(coverage_data) &&
    !methods::is(coverage_data, "GRanges")
  junction_is_list <- is.list(junction_data) && !is.data.frame(junction_data) &&
    !methods::is(junction_data, "GRanges")

  if (coverage_is_list) {
    # Multi-track mode: coverage_data is a list
    if (!junction_is_list) {
      stop(
        "When coverage_data is a list, junction_data must also be a list ",
        "with matching names."
      )
    }

    # Assign names if not present
    if (is.null(names(coverage_data)) && is.null(track_labels)) {
      names(coverage_data) <- paste0("Track ", seq_along(coverage_data))
    } else if (is.null(names(coverage_data)) && !is.null(track_labels)) {
      names(coverage_data) <- track_labels
    }

    if (is.null(names(junction_data))) {
      names(junction_data) <- names(coverage_data)
    }

    # Validate matching names
    if (!setequal(names(coverage_data), names(junction_data))) {
      stop(
        "coverage_data and junction_data lists must have the same names. ",
        "Coverage names: ", paste(names(coverage_data), collapse = ", "), ". ",
        "Junction names: ", paste(names(junction_data), collapse = ", "), "."
      )
    }

    # Process each track
    coverage_list <- list()
    junction_list <- list()

    for (track_name in names(coverage_data)) {
      cov_element <- coverage_data[[track_name]]
      junc_element <- junction_data[[track_name]]

      # Process coverage for this track using process_signal_input logic
      if (methods::is(cov_element, "GRanges")) {
        cov_element <- granges_to_df(cov_element)
      }

      if (is.data.frame(cov_element)) {
        required_cols <- c("seqnames", "start", "end", "score")
        if (!all(required_cols %in% colnames(cov_element))) {
          stop(
            "Coverage data frame must contain columns: ",
            paste(required_cols, collapse = ", ")
          )
        }
        cov_df <- cov_element |>
          dplyr::filter(
            seqnames == as.character(region_gr@seqnames),
            end >= GenomicRanges::start(region_gr),
            start <= GenomicRanges::end(region_gr)
          )
      } else if (is.character(cov_element) && length(cov_element) == 1) {
        cov_df <- import_genomic_data(cov_element, which = region_gr)
      } else {
        stop("List elements must be GRanges, data frames, or single file paths")
      }

      cov_df$track <- track_name
      coverage_list[[track_name]] <- cov_df

      # Process junction for this track
      junction_list[[track_name]] <- process_single_junction(
        junc_element, cov_df, track_name
      )
    }

    coverage_df <- dplyr::bind_rows(coverage_list)
    junction_df <- dplyr::bind_rows(junction_list)

    return(list(coverage_df = coverage_df, junction_df = junction_df))

  } else {
    # Single-track mode: use existing process_sashimi_data logic
    if (junction_is_list) {
      stop(
        "When coverage_data is a single source (not a list), ",
        "junction_data must also be a single source."
      )
    }

    # Process coverage data
    if (methods::is(coverage_data, "GRanges")) {
      coverage_data <- granges_to_df(coverage_data)
    }

    if (is.character(coverage_data) && length(coverage_data) == 1) {
      coverage_df <- import_genomic_data(coverage_data, which = region_gr)
    } else if (is.data.frame(coverage_data)) {
      required_cols <- c("seqnames", "start", "end", "score")
      if (!all(required_cols %in% colnames(coverage_data))) {
        stop(
          "Coverage data frame must contain columns: ",
          paste(required_cols, collapse = ", ")
        )
      }
      coverage_df <- coverage_data |>
        dplyr::filter(
          seqnames == as.character(region_gr@seqnames),
          end >= GenomicRanges::start(region_gr),
          start <= GenomicRanges::end(region_gr)
        )
    } else {
      stop("coverage_data must be a GRanges, data frame, file path, or named list")
    }

    # Process junction data
    junction_df <- process_single_junction(junction_data, coverage_df)

    return(list(coverage_df = coverage_df, junction_df = junction_df))
  }
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
calculate_link_ylim <- function(
  data,
  height_factor = 0.15,
  direction = "down"
) {
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

#' Process interaction input into standardized data frame
#'
#' This function converts any input type (GRanges, GInteractions, data.frame,
#' character vector or list) into a standardized data frame with consistent
#' columns for interaction/link visualization.
#'
#' @param input A GRanges, GInteractions, data frame, character vector of file paths,
#'   or named list of data sources.
#' @param region A genomic region string in the format "chr:start-end"
#' @param track_labels Optional vector of track labels (used for character vector input)
#' @return A data frame with standardized columns: start1, end1, start2, end2,
#'   and optionally score, track, and group columns
#' @export
#' @importFrom dplyr bind_rows mutate filter
#' @importFrom methods is
#' @importFrom rtracklayer import
#' @examples
#' \dontrun{
#' # Data frame input
#' df <- data.frame(
#'   start1 = c(1000, 2000), end1 = c(1500, 2500),
#'   start2 = c(5000, 6000), end2 = c(5500, 6500),
#'   score = c(10, 20)
#' )
#' process_interaction_input(df, "chr1:1-10000")
#'
#' # Character vector input
#' files <- c("file1.bedpe", "file2.bedpe")
#' process_interaction_input(files, "chr1:1-10000", track_labels = c("Sample1", "Sample2"))
#'
#' # List input
#' data_list <- list("Track1" = df, "Track2" = "file.bedpe")
#' process_interaction_input(data_list, "chr1:1-10000")
#' }
process_interaction_input <- function(input, region, track_labels = NULL) {
  region_gr <- parse_region(region)

  # Helper function to import a single file
  import_interaction_file <- function(file, region_gr) {
    gr <- tryCatch(
      {
        rtracklayer::import(file, which = region_gr)
      },
      error = function(e) {
        # Fallback: import all and filter
        gr_all <- rtracklayer::import(file)
        IRanges::subsetByOverlaps(gr_all, region_gr)
      }
    )
    # Process the imported data
    process_interaction_data(gr)
  }

  # Helper to standardize column names (chr1/seqnames1 -> consistent format)
  standardize_interaction_df <- function(df) {
    # Common column name mappings
    if ("seqnames1" %in% colnames(df) && !"chr1" %in% colnames(df)) {
      df$chr1 <- df$seqnames1
    }
    if ("seqnames2" %in% colnames(df) && !"chr2" %in% colnames(df)) {
      df$chr2 <- df$seqnames2
    }
    df
  }

  if (methods::is(input, "GRanges") || methods::is(input, "GInteractions")) {
    # Case 0: GRanges or GInteractions input - process directly
    return(process_interaction_data(input))
  }

  if (is.data.frame(input)) {
    # Case 1: Data frame input
    # Validate required columns
    input <- standardize_interaction_df(input)
    required_cols <- c("start1", "start2")
    if (!all(required_cols %in% colnames(input))) {
      stop(
        "Data frame must contain columns: ",
        paste(required_cols, collapse = ", ")
      )
    }
    return(input)
  } else if (is.character(input)) {
    # Case 2: Character vector input (file paths)
    if (length(input) == 1) {
      # Single file
      track_name <- ifelse(is.null(track_labels), "Track 1", track_labels[1])
      df <- import_interaction_file(input[1], region_gr)
      df$track <- track_name
      return(df)
    } else {
      # Multiple files
      track_data_list <- list()
      for (i in seq_along(input)) {
        track_name <- ifelse(
          is.null(track_labels),
          paste0("Track ", i),
          track_labels[i]
        )
        track_data <- import_interaction_file(input[i], region_gr)
        track_data$track <- track_name
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
        processed_data <- standardize_interaction_df(track_element)
        processed_data$track <- track_name
        track_data_list[[i]] <- processed_data
      } else if (methods::is(track_element, "GRanges") ||
                 methods::is(track_element, "GInteractions")) {
        processed_data <- process_interaction_data(track_element)
        processed_data$track <- track_name
        track_data_list[[i]] <- processed_data
      } else if (is.character(track_element)) {
        # Character vector element (file paths for this track)
        if (length(track_element) == 1) {
          # Single file
          processed_data <- import_interaction_file(track_element, region_gr)
        } else {
          # Multiple files within this track
          file_data_list <- list()
          for (j in seq_along(track_element)) {
            file_data <- import_interaction_file(track_element[j], region_gr)
            file_data$track <- track_name
            file_data$group <- paste0(track_name, "_", j)
            file_data_list[[j]] <- file_data
          }
          processed_data <- dplyr::bind_rows(file_data_list)
        }
        processed_data$track <- track_name
        track_data_list[[i]] <- processed_data
      } else {
        stop("List elements must be data frames, GRanges, GInteractions, or character vectors")
      }
    }

    names(track_data_list) <- names(input)
    return(dplyr::bind_rows(track_data_list))
  } else {
    stop("Input must be a data frame, GRanges, GInteractions, character vector, or named list")
  }
}


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

