# ezGenomeTracks - Helper functions (split from helpers.R)
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
        seqnames == as.character(GenomicRanges::seqnames(region_gr)),
        start >= GenomicRanges::start(region_gr),
        end <= GenomicRanges::end(region_gr)
      ) |>
      dplyr::mutate(name = name)
  } else if (is(input, "character")) {
    # Single track, file name
    # import_genomic_data() already handles BigWig/megadepth vs rtracklayer dispatch
    track_data <- import_genomic_data(input, which = region_gr) |>
      dplyr::mutate(name = name)
  } else {
    stop("Input must be a data.frame or file path (character), got: ",
         paste(class(input), collapse = "/"))
  }

  return(track_data)
}
