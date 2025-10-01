#' Geom for continuous signal tracks
#'
#' This function creates a geom for continuous signal tracks, such as bigWig files,
#' RNA-seq coverage, or ATAC-seq signal. It can display the signal as a line, area,
#' or heatmap.
#'
#' @inheritParams ggplot2::layer
#' @param type Type of signal visualization: "line", "area", or "heatmap" (default: "area").
#'   "line" and "area" types expect `start`, `end`, and `score` columns in the data.
#'   "heatmap" expects `start`, `end`, and `score` columns in the data, with `score` mapped to `fill`.
#' @param na.rm If `TRUE`, silently removes `NA` values.
#' @param ... Other arguments passed on to \code{\link[ggplot2]{layer}}. These are
#'   often aesthetics, used to set an aesthetic to a fixed value, like
#'   `color = "red"` or `linewidth = 3`.
#' @return A ggplot2 layer.
#' @export
#' @importFrom ggplot2 GeomSegment GeomArea GeomTile layer aes ggproto Geom
#' @examples
#' \dontrun{
#' library(ggplot2)
#' signal_data <- data.frame(
#'   start = seq(1, 100, by = 10),
#'   end = seq(1, 100, by = 10) + 9,
#'   score = rnorm(10)
#' )
#' ggplot(signal_data) + geom_signal(aes(x = start, y = score, xend = end))
#' ggplot(signal_data) + geom_signal(aes(x = start, y = score), type = "area")
#' ggplot(signal_data) +
#'   geom_signal(aes(x = start, fill = score), type = "heatmap", y = 1, height = 1) +
#'   scale_fill_viridis_c()
#' }
geom_signal <- function(mapping = NULL, data = NULL, stat = "identity",
                        position = "identity", type = "area", ...,
                        na.rm = TRUE, show.legend = NA, inherit.aes = TRUE) {

  type <- match.arg(type, c("area", "line", "heatmap"))
  if (type == "area") {
     default_aes <- aes(x = .data$start, y = .data$score)
  } else if (type == "line") {
    default_aes <- aes(x = .data$start, y = 0, xend = .data$end, yend = .data$score)
  } else if (type == "heatmap") {
    default_aes <- aes(x = (.data$start + .data$end) / 2, fill = .data$score, y = 1, height = 1)
  }

  if (is.null(mapping)) {
    mapping <- default_aes
  } else {
    mapping <- utils::modifyList(default_aes, as.list(mapping))
    mapping <- do.call(aes, mapping)
  }

  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomSignal,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      type = type,
      na.rm = na.rm,
      ...
    )
  )
}

#' @rdname geom_signal
#' @format NULL
#' @usage NULL
GeomSignal <- ggproto("GeomSignal", Geom,
  required_aes = c("x", "y"),
  setup_params = function(data, params) {
    params$type <- match.arg(params$type, c("area", "line", "heatmap"))
    params
  },

  draw_panel = function(data, panel_params, coord, type = "area", na.rm = FALSE) {
    Geom <- switch(type,
      line = GeomSegment,
      area = GeomArea,
      heatmap = GeomTile
    )
    Geom$draw_panel(data, panel_params, coord)
  },

  default_aes = aes(
    colour = "purple2", fill = "purple2", linewidth = 0.5, linetype = 1,
    alpha = 0.7
  )
)

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
