# ezGenomeTracks - Megadepth helpers for fast BigWig processing

#' Import BigWig coverage using megadepth
#'
#' Fast region-based BigWig queries using the megadepth command-line tool
#' wrapped through the megadepth R package. This provides significant performance
#' improvements over rtracklayer for BigWig file access, especially for multiple
#' region queries and remote files.
#'
#' @param file Path to BigWig file (local or remote URL)
#' @param which GRanges object specifying the genomic region to query
#' @param op Summary operation to compute across the region:
#'   - "mean": average coverage (default)
#'   - "sum": total coverage
#'   - "min": minimum coverage
#'   - "max": maximum coverage
#' @param n_bins Number of bins used to tile `which` before querying.
#'   Higher values preserve more detail; lower values are faster. Default: 2000.
#'
#' @return GRanges object with metadata column "score" containing the computed
#'   coverage values for the input region
#'
#' @details
#' This function uses megadepth::get_coverage() internally by:
#' 1. Converting the input GRanges region to a temporary BED file
#' 2. Calling megadepth to compute coverage statistics
#' 3. Returning the result as a GRanges object compatible with existing code
#'
#' @noRd
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors mcols
#' @keywords internal
import_bigwig_megadepth <- function(file, which, op = "mean", n_bins = 2000L) {
  # Validate op parameter
  op <- match.arg(op, c("sum", "mean", "min", "max"))
  stopifnot(
    "which must be a non-empty GRanges object" = methods::is(which, "GRanges") &&
      length(which) > 0,
    "n_bins must be a positive integer" = is.numeric(n_bins) &&
      length(n_bins) == 1 && n_bins >= 1
  )

  # Create temporary BED file from the region (required by megadepth)
  temp_bed <- tempfile(fileext = ".bed")
  on.exit(unlink(temp_bed), add = TRUE)

  # Build a tiled BED annotation across the requested range.
  # Using only one BED interval returns one summary score (flat signal).
  # Tiling preserves local shape while still leveraging megadepth speed.
  which_reduced <- GenomicRanges::reduce(which)
  total_width <- sum(GenomicRanges::width(which_reduced))
  bin_width <- max(1L, as.integer(ceiling(total_width / as.integer(n_bins))))

  bed_list <- lapply(seq_along(which_reduced), function(i) {
    chr_i <- as.character(GenomicRanges::seqnames(which_reduced)[i])
    start_i <- GenomicRanges::start(which_reduced)[i]
    end_i <- GenomicRanges::end(which_reduced)[i]
    starts <- seq.int(start_i, end_i, by = bin_width)
    ends <- pmin(starts + bin_width - 1L, end_i)
    data.frame(
      chr = rep(chr_i, length(starts)),
      start = starts - 1L, # BED uses 0-based starts
      end = ends
    )
  })
  bed_data <- do.call(rbind, bed_list)

  utils::write.table(
    bed_data,
    file = temp_bed,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

  # Call megadepth::get_coverage() with the annotation file
  # Use a temp prefix to avoid collisions when called repeatedly.
  out_prefix <- tempfile("megadepth_cov_")
  on.exit(unlink(Sys.glob(paste0(out_prefix, "*"))), add = TRUE)
  result_gr <- megadepth::get_coverage(
    bigwig_file = file,
    op = op,
    annotation = temp_bed,
    prefix = out_prefix
  )

  return(result_gr)
}
