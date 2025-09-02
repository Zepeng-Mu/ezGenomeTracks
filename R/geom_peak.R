#' Geom for genomic peak tracks
#'
#' This function creates a geom for genomic peak tracks, such as ChIP-seq peaks,
#' ATAC-seq peaks, or any interval-based features. It displays the peaks as rectangles.
#'
#' @inheritParams ggplot2::geom_rect
#' @param height Height of the peaks (default: 0.8)
#' @param color Border color of the peaks (default: "black")
#' @param fill Fill color of the peaks (default: "gray70")
#' @param alpha Transparency (default: 0.7)
#' @return A ggplot2 layer
#' @export
#' @importFrom ggplot2 geom_rect aes
#' @examples
#' \dontrun{
#' library(ggplot2)
#' p <- ggplot(peak_data) + geom_peak(aes(xmin = start, xmax = end, fill = score))
#' }
geom_peak <- function(mapping = NULL, data = NULL, stat = "identity",
                       position = "identity", ..., height = 0.8,
                       color = "black", fill = "gray70", alpha = 0.7,
                       show.legend = NA, inherit.aes = TRUE) {

  # Default mapping for peaks
  if (is.null(mapping)) {
    mapping <- ggplot2::aes(xmin = start, xmax = end, ymin = 0, ymax = height)
  } else {
    # If mapping is provided but doesn't include y positions, add them
    if (!all(c("ymin", "ymax") %in% names(mapping))) {
      mapping$ymin <- 0
      mapping$ymax <- height
    }
  }

  # Create the rectangle geom
  ggplot2::geom_rect(
    mapping = mapping, data = data, stat = stat,
    position = position, color = color, fill = fill, alpha = alpha, ...,
    show.legend = show.legend, inherit.aes = inherit.aes
  )
}

#' Create a peak track from a BED file
#'
#' This function creates a peak track from a BED file. It imports the data
#' for a specific region and creates a ggplot2 layer for visualization.
#'
#' @param file Path to the BED file
#' @param region Genomic region to display (e.g., "chr1:1000000-2000000")
#' @param color Border color of the peaks (default: "black")
#' @param fill Fill color of the peaks (default: "gray70")
#' @param alpha Transparency (default: 0.7)
#' @param height Height of the peaks (default: 0.8)
#' @param use_score Use the score column for fill color (default: FALSE)
#' @param ... Additional arguments passed to geom_peak
#' @return A ggplot2 layer
#' @export
#' @importFrom ggplot2 ggplot aes scale_fill_gradient
#' @examples
#' \dontrun{
#' p <- peak_track("peaks.bed", "chr1:1000000-2000000", use_score = TRUE)
#' }
peak_track <- function(file, region, color = "black", fill = "gray70",
                       alpha = 0.7, height = 0.8, use_score = FALSE, ...) {
  # Parse the region
  region_gr <- parse_region(region)

  # Import the data
  peak_data <- import_genomic_data(file, which = region_gr)

  # Create the plot
  if (use_score && "score" %in% colnames(peak_data)) {
    p <- ggplot2::ggplot(peak_data) +
      geom_peak(ggplot2::aes(xmin = start, xmax = end, fill = score),
                color = color, alpha = alpha, height = height, ...) +
      ggplot2::scale_fill_gradient(low = "white", high = fill)
  } else {
    p <- ggplot2::ggplot(peak_data) +
      geom_peak(ggplot2::aes(xmin = start, xmax = end),
                color = color, fill = fill, alpha = alpha, height = height, ...)
  }

  # Apply the appropriate theme and scale
  p <- p + ez_peak_theme() +
    scale_x_genome_region(region) +
    ggplot2::ylim(0, 1)  # Fixed y-axis for peaks

  return(p)
}