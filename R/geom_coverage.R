#' Coverage visualization geom
#'
#' Visualize quantitative genomic signal as line, area, or heatmap tiles.
#' Input data must contain genomic coordinates (`start`, `end`) and a numeric
#' signal value (`score`). The geom automatically maps these columns to the
#' required aesthetics for the chosen `type`.
#'
#' @inheritParams ggplot2::layer
#' @param type Visualization style: `"area"` (default), `"line"`, or `"heatmap"`.
#'   - `"area"`/`"line"`: `score` is mapped to `y` (height).
#'   - `"heatmap"`: `score` is mapped to `fill`, producing colored tiles.
#' @param na.rm If `TRUE`, silently drop `NA` values.
#' @param ... Additional arguments passed to [ggplot2::layer()], e.g.
#'   `color = "black"`, `linewidth = 0.8`, or `alpha = 0.6`.
#'
#' @return A ggplot2 layer that can be added to a plot.
#' @export
#' @importFrom ggplot2 GeomRect GeomTile layer aes ggproto Geom
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' df <- data.frame(
#'   start = seq(1, 100, 10),
#'   end   = seq(10, 100, 10),
#'   score = rnorm(10)
#' )
#'
#' # Area plot (default)
#' ggplot(df) +
#'   geom_coverage()
#'
#' # Line plot
#' ggplot(df) +
#'   geom_coverage(type = "line")
#'
#' # Heatmap tiles
#' ggplot(df) +
#'   geom_coverage(type = "heatmap") +
#'   scale_fill_viridis_c()
#' }
geom_coverage <- function(mapping = NULL, data = NULL, stat = "identity",
                          position = "identity", type = "area", ...,
                          na.rm = TRUE, show.legend = NA, inherit.aes = TRUE) {
  type <- match.arg(type, c("area", "line", "heatmap"))
  if (type == "area") {
    default_aes <- aes(xmin = .data$start, xmax = .data$end,
                       ymin = 0, ymax = .data$score)
  } else if (type == "line") {
    default_aes <- aes(xmin = .data$start, xmax = .data$end,
                       ymin = .data$score, ymax = .data$score)
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

#' @rdname geom_coverage
#' @format NULL
#' @usage NULL
GeomSignal <- ggproto("GeomSignal", Geom,
  required_aes = c("xmin", "xmax", "ymin", "ymax"),
  setup_params = function(data, params) {
    params$type <- match.arg(params$type, c("area", "line", "heatmap"))
    params
  },
  draw_panel = function(data, panel_params, coord, type = "area", na.rm = FALSE) {
    if (type == "heatmap") {
      # For heatmap, transform xmin/xmax to x and keep y/height
      data$x <- (data$xmin + data$xmax) / 2
      data$y <- 1
      data$height <- 1
      GeomTile$draw_panel(data, panel_params, coord)
    } else {
      # Use GeomRect for area and line (which accepts xmin, xmax, ymin, ymax)
      GeomRect$draw_panel(data, panel_params, coord)
    }
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
#'   stat_bin_signal(binwidth = 1000)
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
