#' Stack multiple genome tracks
#'
#' This function stacks multiple genome tracks vertically, ensuring they share
#' a common x-axis. It uses the aplot package for track stacking.
#'
#' @param ... ggplot2 objects representing genome tracks
#' @param region Genomic region to display (e.g., "chr1:1000000-2000000")
#' @param heights Relative heights of the tracks (default: NULL, equal heights)
#' @param ncol Number of columns for the legend (default: 1)
#' @param align_tracks Align tracks by x-axis (default: TRUE)
#' @return A composite plot with stacked tracks
#' @export
#' @importFrom aplot plot_list
#' @examples
#' \dontrun{
#' track1 <- ez_coverage("signal.bw", "chr1:1000000-2000000")
#' track2 <- ez_peak("peaks.bed", "chr1:1000000-2000000")
#' track3 <- ez_gene("genes.gtf", "chr1:1000000-2000000")
#' p <- genome_plot(track1, track2, track3, region = "chr1:1000000-2000000")
#' }
genome_plot <- function(..., region = NULL, heights = NULL, ncol = 1, align_tracks = TRUE) {
  # Collect the tracks
  tracks <- list(...)

  # Check if we have any tracks
  if (length(tracks) == 0) {
    stop("No tracks provided")
  }

  # If region is provided, apply it to all tracks
  if (!is.null(region)) {
    tracks <- lapply(tracks, function(track) {
      track + scale_x_genome_region(region)
    })
  }

  # Stack the tracks using aplot
  if (align_tracks) {
    p <- aplot::plot_list(
      plotlist = tracks,
      heights = heights,
      ncol = 1,
      guides = "collect",
      legend_ncol = ncol
    )
  } else {
    # If not aligning tracks, just use aplot::plot_list without alignment
    p <- aplot::plot_list(
      plotlist = tracks,
      heights = heights,
      ncol = 1,
      guides = "collect",
      legend_ncol = ncol,
      align = "none"
    )
  }

  return(p)
}

#' Add a vertical line to a genome track
#'
#' This function adds a vertical line to a genome track at a specific position.
#' It is useful for highlighting specific genomic positions.
#'
#' @param plot A ggplot2 object representing a genome track
#' @param position Genomic position for the vertical line
#' @param color Color of the line (default: "red")
#' @param size Size of the line (default: 0.5)
#' @param linetype Line type (default: "dashed")
#' @param alpha Transparency (default: 0.7)
#' @return A ggplot2 object with the vertical line added
#' @export
#' @importFrom ggplot2 geom_vline
#' @examples
#' \dontrun{
#' track <- ez_coverage("signal.bw", "chr1:1000000-2000000")
#' track <- add_vline(track, 1500000)
#' }
add_vline <- function(plot, position, color = "red", linewidth = 0.5,
                      linetype = "dashed", alpha = 0.7) {
  plot + ggplot2::geom_vline(
    xintercept = position,
    linewidth = linewidth,
    linetype = linetype,
    alpha = alpha
  )
}

#' Add a horizontal line to a genome track
#'
#' This function adds a horizontal line to a genome track at a specific y-value.
#' It is useful for highlighting specific signal thresholds.
#'
#' @param plot A ggplot2 object representing a genome track
#' @param y Y-value for the horizontal line
#' @param color Color of the line (default: "blue")
#' @param size Size of the line (default: 0.5)
#' @param linetype Line type (default: "dashed")
#' @param alpha Transparency (default: 0.7)
#' @return A ggplot2 object with the horizontal line added
#' @export
#' @importFrom ggplot2 geom_hline
#' @examples
#' \dontrun{
#' track <- ez_coverage("signal.bw", "chr1:1000000-2000000")
#' track <- add_hline(track, 10)
#' }
add_hline <- function(plot, y, color = "blue", size = 0.5,
                      linetype = "dashed", alpha = 0.7) {
  plot + ggplot2::geom_hline(
    yintercept = y,
    color = color,
    size = size,
    linetype = linetype,
    alpha = alpha
  )
}

#' Add a rectangle to a genome track
#'
#' This function adds a rectangle to a genome track between two genomic positions.
#' It is useful for highlighting specific genomic regions.
#'
#' @param plot A ggplot2 object representing a genome track
#' @param xmin Minimum x-value (start position)
#' @param xmax Maximum x-value (end position)
#' @param ymin Minimum y-value (default: -Inf)
#' @param ymax Maximum y-value (default: Inf)
#' @param fill Fill color of the rectangle (default: "yellow")
#' @param alpha Transparency (default: 0.2)
#' @param color Border color of the rectangle (default: NA)
#' @return A ggplot2 object with the rectangle added
#' @export
#' @importFrom ggplot2 geom_rect
#' @examples
#' \dontrun{
#' track <- ez_coverage("signal.bw", "chr1:1000000-2000000")
#' track <- add_rect(track, 1200000, 1400000)
#' }
add_rect <- function(plot, xmin, xmax, ymin = -Inf, ymax = Inf,
                      fill = "yellow", alpha = 0.2, color = NA) {
  plot + ggplot2::geom_rect(
    xmin = xmin,
    xmax = xmax,
    ymin = ymin,
    ymax = ymax,
    fill = fill,
    alpha = alpha,
    color = color
  )
}

#' Add text annotation to a genome track
#'
#' This function adds text annotation to a genome track at a specific position.
#' It is useful for labeling specific genomic features.
#'
#' @param plot A ggplot2 object representing a genome track
#' @param x X-value (genomic position)
#' @param y Y-value
#' @param label Text label
#' @param color Text color (default: "black")
#' @param size Text size (default: 3)
#' @param angle Text angle (default: 0)
#' @param hjust Horizontal justification (default: 0.5)
#' @param vjust Vertical justification (default: 0.5)
#' @return A ggplot2 object with the text annotation added
#' @export
#' @importFrom ggplot2 geom_text
#' @examples
#' \dontrun{
#' track <- ez_coverage("signal.bw", "chr1:1000000-2000000")
#' track <- add_text(track, 1500000, 20, "Peak")
#' }
add_text <- function(plot, x, y, label, color = "black", size = 3,
                      angle = 0, hjust = 0.5, vjust = 0.5) {
  plot + ggplot2::geom_text(
    x = x,
    y = y,
    label = label,
    color = color,
    size = size,
    angle = angle,
    hjust = hjust,
    vjust = vjust
  )
}