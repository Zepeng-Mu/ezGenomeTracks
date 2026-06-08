# ezGenomeTracks - Helper functions (split from helpers.R)

# Internal helper: Get coverage value at a genomic position
# Finds the first overlapping coverage bin and returns max(score), or 0 if none.
#
# @param pos Genomic position (integer)
# @param cov_df Data frame with start, end, score columns
# @return Numeric coverage value
get_coverage_at_pos <- function(pos, cov_df) {
  idx <- which(cov_df$start <= pos & cov_df$end >= pos)
  if (length(idx) == 0) return(0)
  return(max(cov_df$score[idx], na.rm = TRUE))
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
  process_sashimi_input(
    coverage_data = coverage_data,
    junction_data = junction_data,
    region = region,
    junction_direction = junction_direction
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

  # Helper function to process a single junction source
  # direction_lookup: optional data frame with columns start, end, arc_direction
  #   providing a globally consistent direction assignment across tracks
  process_single_junction <- function(junc_input, cov_df, track_name = NULL, direction_lookup = NULL) {
    # Import/filter junction data
    if (methods::is(junc_input, "GRanges")) {
      junc_input <- granges_to_df(junc_input)
    }

    if (is.character(junc_input) && length(junc_input) == 1) {
      junction_df <- import_genomic_data(junc_input, which = region_gr)
    } else if (is.data.frame(junc_input)) {
      required_cols <- c("seqnames", "start", "end", "score")
      validate_required_columns(junc_input, required_cols, "Junction data frame")
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
        arc_direction = if (!is.null(direction_lookup)) {
          # Use global lookup for consistency across tracks
          lookup_key <- paste(start, end)
          global_key <- paste(direction_lookup$start, direction_lookup$end)
          dir <- direction_lookup$arc_direction[match(lookup_key, global_key)]
          ifelse(is.na(dir), "up", dir)
        } else if (junction_direction == "alternate") {
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
    coverage_data <- ensure_track_names(coverage_data, track_labels)

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

    # Build a global arc direction lookup so the same junction has the same
    # direction across all tracks (only needed for "alternate" mode)
    direction_lookup <- NULL
    if (junction_direction == "alternate") {
      all_coords <- lapply(names(coverage_data), function(nm) {
        junc_el <- junction_data[[nm]]
        if (methods::is(junc_el, "GRanges")) junc_el <- granges_to_df(junc_el)
        if (is.character(junc_el) && length(junc_el) == 1) {
          junc_el <- import_genomic_data(junc_el, which = region_gr)
        }
        if (is.data.frame(junc_el)) {
          junc_el <- junc_el |>
            dplyr::filter(
              seqnames == as.character(region_gr@seqnames),
              end >= GenomicRanges::start(region_gr),
              start <= GenomicRanges::end(region_gr)
            )
          return(junc_el[, c("start", "end"), drop = FALSE])
        }
        return(data.frame(start = integer(0), end = integer(0)))
      })
      all_coords_df <- unique(dplyr::bind_rows(all_coords))
      if (nrow(all_coords_df) > 0) {
        all_coords_df <- all_coords_df[order(all_coords_df$start), , drop = FALSE]
        all_coords_df$arc_direction <- ifelse(
          seq_len(nrow(all_coords_df)) %% 2 == 1, "up", "down"
        )
        direction_lookup <- all_coords_df[, c("start", "end", "arc_direction")]
      }
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
        validate_required_columns(cov_element, required_cols, "Coverage data frame")
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
        junc_element, cov_df, track_name, direction_lookup = direction_lookup
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
      validate_required_columns(coverage_data, required_cols, "Coverage data frame")
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
