# ezGenomeTracks - Helper functions (split from helpers.R)
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
    validate_required_columns(input, required_cols, "Data frame")
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

