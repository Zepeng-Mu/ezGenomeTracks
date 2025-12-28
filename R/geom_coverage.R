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
#'   `linewidth = 0.8`, or `alpha = 0.6`.
#'
#' @return A ggplot2 layer that can be added to a plot.
#' @export
#' @importFrom ggplot2 GeomRect GeomTile GeomPath layer aes ggproto Geom
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
geom_coverage <- function(
  mapping = NULL,
  data = NULL,
  stat = "identity",
  position = "identity",
  type = "area",
  ...,
  na.rm = TRUE,
  show.legend = NA,
  inherit.aes = TRUE
) {
  type <- match.arg(type, c("area", "line", "heatmap"))
  if (type == "area") {
    default_aes <- aes(
      xmin = .data$start,
      xmax = .data$end,
      ymin = 0,
      ymax = .data$score
    )
  } else if (type == "line") {
    default_aes <- aes(
      xmin = .data$start,
      xmax = .data$end,
      ymin = .data$score,
      ymax = .data$score
    )
  } else if (type == "heatmap") {
    default_aes <- aes(
      x = (.data$start + .data$end) / 2,
      fill = .data$score,
      y = 1,
      height = 1
    )
  }

  if (is.null(mapping)) {
    mapping <- default_aes
  } else {
    mapping <- utils::modifyList(default_aes, as.list(mapping))
    mapping <- do.call(aes, mapping)
  }

  params_list <- utils::modifyList(
    list(
      type = type,
      na.rm = na.rm
    ),
    list(...)
  )

  geom_layer <- layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomCoverage,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = params_list
  )

  geom_layer
}

#' @rdname geom_coverage
#' @format NULL
#' @usage NULL
#' @importFrom ggplot2 GeomPath
GeomCoverage <- ggproto(
  "GeomCoverage",
  Geom,
  required_aes = c("xmin", "xmax", "ymin", "ymax"),
  setup_params = function(data, params) {
    params$type <- match.arg(params$type, c("area", "line", "heatmap"))
    params
  },
  draw_panel = function(
    data,
    panel_params,
    coord,
    type = "area",
    na.rm = FALSE
  ) {
    if (type == "heatmap") {
      # For heatmap, transform xmin/xmax to x and keep y/height
      data$x <- (data$xmin + data$xmax) / 2
      data$y <- 1
      data$height <- 1
      GeomTile$draw_panel(data, panel_params, coord)
    } else if (type == "line") {
      # For line type, create step-like path connecting score values
      # Use 'fill' aesthetic as the line colour (since fill is the user-facing "color" param)
      n <- nrow(data)
      if (n == 0) return(grid::nullGrob())

      # Determine the grouping variable
      if ("group" %in% names(data)) {
        groups <- unique(data$group)
      } else {
        groups <- 1
        data$group <- 1
      }

      # Build path data for each group
      path_data_list <- lapply(groups, function(g) {
        gdata <- data[data$group == g, , drop = FALSE]
        gdata <- gdata[order(gdata$xmin), , drop = FALSE]
        gn <- nrow(gdata)
        if (gn == 0) return(NULL)

        # Create step pattern: start -> end at ymax level for each bin
        x_coords <- numeric(2 * gn)
        y_coords <- numeric(2 * gn)
        for (i in seq_len(gn)) {
          x_coords[2 * i - 1] <- gdata$xmin[i]
          x_coords[2 * i] <- gdata$xmax[i]
          y_coords[2 * i - 1] <- gdata$ymax[i]
          y_coords[2 * i] <- gdata$ymax[i]
        }

        # Use fill as the line colour (fill is the user-facing "color" of the track)
        line_colour <- if ("fill" %in% names(gdata) &&
                           length(gdata$fill) > 0 &&
                           !all(is.na(gdata$fill))) {
          gdata$fill[1]
        } else {
          "steelblue"
        }

        data.frame(
          x = x_coords,
          y = y_coords,
          PANEL = gdata$PANEL[1],
          group = g,
          colour = line_colour,
          linewidth = gdata$linewidth[1],
          linetype = gdata$linetype[1],
          alpha = gdata$alpha[1]
        )
      })

      path_data <- do.call(rbind, path_data_list)

      GeomPath$draw_panel(path_data, panel_params, coord)
    } else {
      # Use GeomRect for area type - suppress borders by setting colour = NA
      data$colour <- NA
      GeomRect$draw_panel(data, panel_params, coord)
    }
  },
  default_aes = aes(
    colour = NA,
    fill = "purple2",
    linewidth = 0.5,
    linetype = 1,
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
stat_bin_signal <- function(
  mapping = NULL,
  data = NULL,
  geom = "line",
  position = "identity",
  ...,
  binwidth = NULL,
  bins = 30,
  summary_fun = mean,
  show.legend = NA,
  inherit.aes = TRUE
) {
  ggplot2::stat_summary_bin(
    mapping = mapping,
    data = data,
    geom = geom,
    position = position,
    fun = summary_fun,
    ...,
    binwidth = binwidth,
    bins = bins,
    show.legend = show.legend,
    inherit.aes = inherit.aes
  )
}
