# ezGenomeTracks - Helper functions (split from helpers.R)
#' Process signal input into standardized data frame
#'
#' This function converts any input type (GRanges, data.frame, character vector, or list)
#' into a standardized data frame with consistent columns for signal visualization.
#'
#' @param input A GRanges object, data frame, character vector of file paths, or named list.
#'   GRanges objects are automatically converted to data frames.
#' @param region A genomic region string in the format "chr:start-end"
#' @param track_labels Optional vector of track labels (used for character vector input)
#' @param n_bins Maximum number of bins used for BigWig queries
#'   through megadepth (default: 2000).
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
#' process_coverage_input(gr, "chr1:1-100")
#'
#' # Data frame input
#' df <- data.frame(seqnames = "chr1", start = 1:100, end = 1:100, score = rnorm(100))
#' process_coverage_input(df, "chr1:1-100")
#'
#' # Character vector input
#' files <- c("file1.bw", "file2.bw")
#' process_coverage_input(files, "chr1:1-100", track_labels = c("Sample1", "Sample2"))
#'
#' # List input
#' data_list <- list("Track1" = df, "Track2" = files)
#' process_coverage_input(data_list, "chr1:1-100")
#' }
process_coverage_input <- function(
  input,
  region,
  track_labels = NULL,
  n_bins = 2000L
) {
  if (methods::is(input, "GRanges")) {
    # Case 0: GRanges input - convert to data frame and process
    input <- granges_to_df(input)
  }

  if (is.data.frame(input)) {
    # Case 1: Data frame input
    # Validate required columns
    required_cols <- c("seqnames", "start", "end", "score")
    validate_required_columns(input, required_cols, "Data frame")

    # Filter by region if needed
    region_gr <- parse_region(region)
    filtered_data <- input |>
      dplyr::filter(
        seqnames == as.character(GenomicRanges::seqnames(region_gr)),
        start >= GenomicRanges::start(region_gr),
        end <= GenomicRanges::end(region_gr)
      )

    return(filtered_data)
  } else if (is.character(input)) {
    # Case 2: Character vector input (file paths)
    if (length(input) == 1) {
      # Single file
      track_name <- ifelse(is.null(track_labels), "Track 1", track_labels[1])
      return(get_single_coverage(
        input,
        region,
        name = track_name,
        n_bins = n_bins
      ))
    } else {
      # Multiple files
      track_data_list <- list()
      for (i in seq_along(input)) {
        track_name <- ifelse(
          is.null(track_labels),
          paste0("Track ", i),
          track_labels[i]
        )
        track_data <- get_single_coverage(
          input[i],
          region,
          name = track_name,
          n_bins = n_bins
        )
        track_data$group <- track_name
        track_data_list[[i]] <- track_data
      }
      return(dplyr::bind_rows(track_data_list))
    }
  } else if (is.list(input)) {
    # Case 3: List input
    input <- ensure_track_names(input, track_labels)

    track_data_list <- list()
    for (i in seq_along(input)) {
      track_name <- names(input)[i]
      track_element <- input[[i]]

      if (is.data.frame(track_element)) {
        # Data frame element
        processed_data <- process_coverage_input(
          track_element,
          region,
          n_bins = n_bins
        )
        processed_data$track <- track_name
        track_data_list[[i]] <- processed_data
      } else if (is.character(track_element)) {
        # Character vector element (multiple files for this track)
        if (length(track_element) == 1) {
          # Single file
          processed_data <- get_single_coverage(
            track_element,
            region,
            name = track_name,
            n_bins = n_bins
          )
        } else {
          # Multiple files within this track
          file_data_list <- list()
          for (j in seq_along(track_element)) {
            file_data <- get_single_coverage(
              track_element[j],
              region,
              name = paste0(track_name, "_", j),
              n_bins = n_bins
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
