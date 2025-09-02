#' Geom for continuous signal tracks
#'
#' This function creates a geom for continuous signal tracks, such as bigWig files,
#' RNA-seq coverage, or ATAC-seq signal. It can display the signal as a line, area,
#' or heatmap.
#'
#' @inheritParams ggplot2::geom_line
#' @param type Type of signal visualization: "line", "area", or "heatmap" (default: "area")
#' @param baseline Baseline value for area plots (default: 0)
#' @param fill Fill color for area plots (default: "steelblue")
#' @param color Line color (default: "steelblue")
#' @param alpha Transparency (default: 0.5)
#' @return A ggplot2 layer
#' @export
#' @importFrom ggplot2 geom_line geom_area geom_tile aes
#' @examples
#' \dontrun{
#' library(ggplot2)
#' p <- ggplot(signal_data, aes(x = start, y = score)) + geom_signal()
#' }
geom_signal <- function(mapping = NULL, data = NULL, stat = "identity",
                         position = "identity", ..., type = "area",
                         baseline = 0, fill = "steelblue", color = "steelblue",
                         alpha = 0.5, show.legend = NA, inherit.aes = TRUE) {
  
  # Create the appropriate geom based on the type
  if (type == "line") {
    return(ggplot2::geom_line(
      mapping = mapping, data = data, stat = stat,
      position = position, color = color, ...,
      show.legend = show.legend, inherit.aes = inherit.aes
    ))
  } else if (type == "area") {
    return(ggplot2::geom_area(
      mapping = mapping, data = data, stat = stat,
      position = position, fill = fill, color = color, alpha = alpha, ...,
      show.legend = show.legend, inherit.aes = inherit.aes,
      baseline = baseline
    ))
  } else if (type == "heatmap") {
    # For heatmap, we need to ensure we have the right mapping
    if (is.null(mapping)) {
      mapping <- ggplot2::aes()
    }
    if (!"fill" %in% names(mapping)) {
      mapping$fill <- substitute(score)
    }
    
    return(ggplot2::geom_tile(
      mapping = mapping, data = data, stat = stat,
      position = position, alpha = alpha, ...,
      show.legend = show.legend, inherit.aes = inherit.aes
    ))
  } else {
    stop("Type must be one of 'line', 'area', or 'heatmap'")
  }
}

#' Stat for binning genomic signal data
#'
#' This function creates a stat for binning genomic signal data. It is useful for
#' reducing the size of large datasets for visualization.
#'
#' @inheritParams ggplot2::stat_bin
#' @param binwidth Width of bins in base pairs
#' @param bins Number of bins
#' @param summary_fun Function to summarize values within each bin (default: mean)
#' @return A ggplot2 layer
#' @export
#' @importFrom ggplot2 stat_summary_bin
#' @examples
#' \dontrun{
#' library(ggplot2)
#' p <- ggplot(signal_data, aes(x = start, y = score)) + 
#'      stat_bin_signal(binwidth = 1000)
#' }
stat_bin_signal <- function(mapping = NULL, data = NULL, geom = "line",
                            position = "identity", ..., binwidth = NULL,
                            bins = 30, summary_fun = mean,
                            show.legend = NA, inherit.aes = TRUE) {
  
  ggplot2::stat_summary_bin(
    mapping = mapping, data = data, geom = geom,
    position = position, fun = summary_fun, ...,
    binwidth = binwidth, bins = bins,
    show.legend = show.legend, inherit.aes = inherit.aes
  )
}

#' Create a signal track from a bigWig file
#'
#' This function creates a signal track from a bigWig file. It imports the data
#' for a specific region and creates a ggplot2 layer for visualization.
#'
#' @param file Path to the bigWig file
#' @param region Genomic region to display (e.g., "chr1:1000000-2000000")
#' @param type Type of signal visualization: "line", "area", or "heatmap" (default: "area")
#' @param color Line color (default: "steelblue")
#' @param fill Fill color for area plots (default: "steelblue")
#' @param alpha Transparency (default: 0.5)
#' @param binwidth Width of bins in base pairs (default: NULL)
#' @param ... Additional arguments passed to geom_signal
#' @return A ggplot2 layer
#' @export
#' @importFrom ggplot2 ggplot aes
#' @examples
#' \dontrun{
#' p <- signal_track("signal.bw", "chr1:1000000-2000000")
#' }
signal_track <- function(file, region, type = "area", color = "steelblue",
                        fill = "steelblue", alpha = 0.5, binwidth = NULL, ...) {
  # Parse the region
  region_gr <- parse_region(region)
  
  # Import the data
  signal_data <- import_genomic_data(file, which = region_gr)
  
  # Create the plot
  p <- ggplot2::ggplot(signal_data, ggplot2::aes(x = start, y = score)) +
    geom_signal(type = type, color = color, fill = fill, alpha = alpha, ...)
  
  # Apply binning if requested
  if (!is.null(binwidth)) {
    p <- p + stat_bin_signal(binwidth = binwidth)
  }
  
  # Apply the appropriate theme and scale
  p <- p + ez_signal_theme() +
    scale_x_genome_region(region)
  
  return(p)
}