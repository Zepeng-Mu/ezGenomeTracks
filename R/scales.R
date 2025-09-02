#' Scale for genomic coordinates on the x-axis
#'
#' This function creates a continuous scale for genomic coordinates on the x-axis.
#' It formats the axis labels in a genomic coordinate style (e.g., 1Mb, 500kb).
#'
#' @param ... Additional arguments passed to scale_x_continuous
#' @param unit_suffix Suffix to use for the unit (default: "b" for base pairs)
#' @param breaks Breaks for the axis (default: waiver())
#' @param labels Labels for the axis (default: waiver())
#' @return A ggplot2 scale object
#' @export
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom scales label_number_si
#' @examples
#' \dontrun{
#' library(ggplot2)
#' p <- ggplot(data, aes(x = start)) + geom_point() + scale_x_genome()
#' }
scale_x_genome <- function(..., unit_suffix = "b", breaks = waiver(), labels = waiver()) {
  if (is.function(labels)) {
    # Use provided label function
    label_fun <- labels
  } else if (is.waive(labels)) {
    # Create a function to format genomic coordinates
    label_fun <- function(x) {
      scales::label_number_si(unit = unit_suffix)(x)
    }
  } else {
    # Use provided labels
    label_fun <- labels
  }
  
  ggplot2::scale_x_continuous(..., breaks = breaks, labels = label_fun)
}

#' Format genomic coordinates
#'
#' This function formats genomic coordinates in a human-readable way (e.g., 1Mb, 500kb).
#'
#' @param x Numeric vector of genomic coordinates
#' @param unit_suffix Suffix to use for the unit (default: "b" for base pairs)
#' @return Character vector of formatted coordinates
#' @export
#' @importFrom scales label_number_si
#' @examples
#' \dontrun{
#' format_genomic_coord(c(1000, 1000000, 1500000))
#' # [1] "1kb" "1Mb" "1.5Mb"
#' }
format_genomic_coord <- function(x, unit_suffix = "b") {
  scales::label_number_si(unit = unit_suffix)(x)
}

#' Create a genomic position scale for a specific chromosome region
#'
#' This function creates a continuous scale for a specific chromosome region.
#' It sets the limits to the start and end positions of the region and formats
#' the axis labels in a genomic coordinate style.
#'
#' @param region A string specifying a genomic region (e.g., "chr1:1000000-2000000") or a GRanges object
#' @param ... Additional arguments passed to scale_x_genome
#' @return A ggplot2 scale object
#' @export
#' @importFrom GenomicRanges start end
#' @importFrom methods is
#' @examples
#' \dontrun{
#' library(ggplot2)
#' p <- ggplot(data, aes(x = start)) + geom_point() + 
#'      scale_x_genome_region("chr1:1000000-2000000")
#' }
scale_x_genome_region <- function(region, ...) {
  # Parse region if it's a string
  if (is.character(region)) {
    region_gr <- parse_region(region)
  } else if (methods::is(region, "GRanges")) {
    region_gr <- region
  } else {
    stop("Region must be a character string or GRanges object")
  }
  
  # Extract start and end positions
  start_pos <- GenomicRanges::start(region_gr)
  end_pos <- GenomicRanges::end(region_gr)
  
  # Create scale with limits set to the region
  scale_x_genome(..., limits = c(start_pos, end_pos))
}