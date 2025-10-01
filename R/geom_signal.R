#' Geom for continuous signal tracks
#'
#' This function creates a geom for continuous signal tracks, such as bigWig files,
#' RNA-seq coverage, or ATAC-seq signal. It can display the signal as a line, area,
#' or heatmap.
#'
#' @inheritParams ggplot2::geom_line
#' @param type Type of signal visualization: "line", "area", or "heatmap" (default: "area").
#'   "line" and "area" types expect `start`, `end`, and `score` columns in the data.
#'   "heatmap" expects `start`, `end`, and `score` columns in the data, with `score` mapped to `fill`.
#' @param fill Fill color for area plots (default: "purple2").
#' @param color Line color (default: "purple2").
#' @param alpha Transparency (default: 0.5).
#' @return A ggplot2 layer.
#' @export
#' @importFrom ggplot2 GeomSegment geom_tile aes geom_segment geom_ribbon
#' @examples
#' \dontrun{
#' library(ggplot2)
#' p <- ggplot(signal_data, aes(x = start, y = score)) + geom_signal()
#' }
#' Create a signal track for genomic data visualization
#'
#' This function creates a signal track visualization using different geom types.
#' It supports three visualization types: 'line', 'area', and 'heatmap', each with
#' customizable parameters.
#'
#' @param mapping Aesthetic mapping created with aes()
#' @param data The data to be displayed
#' @param stat The statistical transformation to use
#' @param position Position adjustment
#' @param ... Common parameters passed to all geom types
#' @param type Type of visualization: "line", "area", or "heatmap"
#' @param line.params List of parameters specific to line type visualization
#' @param area.params List of parameters specific to area type visualization
#' @param heatmap.params List of parameters specific to heatmap type visualization
#' @param show.legend Logical. Should this layer be included in the legends?
#' @param inherit.aes If FALSE, overrides the default aesthetics
#'
#' @return A ggplot2 layer or list of layers
geom_signal <- function(mapping = NULL, data = NULL, stat = "identity",
                        position = "identity", type = "area",
                        plot.params = list(), ...,
                        show.legend = NA, inherit.aes = TRUE) {

  # Validate that mapping is created by aes()
  if (!is.null(mapping) && !ggplot2::is.ggproto(mapping) && !inherits(mapping, "uneval")) {
    stop("`mapping` must be created by `aes()`.")
  }

  if (type == "line") {
    # Base aesthetics for line type
    base_aes <- ggplot2::aes(x = .data$start, xend = .data$end, y = 0, yend = .data$score)

    # Combine with user-provided mapping if it exists
    if (!is.null(mapping)) {
      mapping <- modifyList(base_aes, mapping)
    } else {
      mapping <- base_aes
    }

    ggplot2::layer(
      data = data,
      mapping = mapping,
      stat = stat,
      geom = GeomSegment,
      position = position,
      show.legend = show.legend,
      inherit.aes = inherit.aes,
      params = list(...)
    )
  } else if (type == "area") {
    # Base aesthetics for area type
    base_aes <- ggplot2::aes(x = .data$start, y = .data$score)

    # Combine with user-provided mapping if it exists
    if (!is.null(mapping)) {
      mapping <- modifyList(base_aes, mapping)
    } else {
      mapping <- base_aes
    }

    ggplot2::layer(
      data = data,
      mapping = mapping,
      stat = stat,
      geom = GeomArea,
      position = position,
      show.legend = show.legend,
      inherit.aes = inherit.aes,
      params = list(...)
    )

  } else {
    stop("Type must be one of 'line' or 'area'")
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
