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
import_bigwig_megadepth <- function(file, which, op = "mean") {
  # Validate op parameter
  op <- match.arg(op, c("sum", "mean", "min", "max"))

  # Create temporary BED file from the region (required by megadepth)
  temp_bed <- tempfile(fileext = ".bed")
  on.exit(unlink(temp_bed), add = TRUE)

  # Convert GRanges region to BED format (3 columns: chr, start, end)
  bed_data <- data.frame(
    chr = as.character(GenomicRanges::seqnames(which)),
    start = GenomicRanges::start(which) - 1,  # BED uses 0-based coordinates
    end = GenomicRanges::end(which)
  )

  utils::write.table(
    bed_data,
    file = temp_bed,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

  # Call megadepth::get_coverage() with the annotation file
  result_gr <- megadepth::get_coverage(
    bw_file = file,
    op = op,
    annotation = temp_bed,
    verbose = FALSE
  )

  return(result_gr)
}
