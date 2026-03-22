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
#' @param area_border Logical; if `TRUE` (default), draws thin borders on area rectangles
#'   to eliminate white-line rendering artifacts. Set to `FALSE` for cleaner appearance
#'   when artifacts are not visible. Only affects `type = "area"`.
#' @param na.rm If `TRUE`, silently drop `NA` values.
#' @param ... Additional arguments passed to [ggplot2::layer()], e.g.
#'   `linewidth = 0.8`, or `alpha = 0.6`.
#'
#' @return A ggplot2 layer that can be added to a plot.
#' @export
#' @importFrom ggplot2 GeomRect GeomTile GeomPath layer aes ggproto Geom
#' @importFrom rlang `%||%`
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
  area_border = TRUE,
  ...,
  na.rm = TRUE,
  show.legend = NA,
  inherit.aes = TRUE
) {
  type <- match.arg(type, c("area", "line", "heatmap"))
  if (type == "area") {
    default_aes <- aes(
      start = .data$start,
      end = .data$end,
      ymin = 0,
      ymax = .data$score
    )
  } else if (type == "line") {
    default_aes <- aes(
      start = .data$start,
      end = .data$end,
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
      area_border = area_border,
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
  required_aes = c("start", "end", "ymin", "ymax"),
  setup_params = function(data, params) {
    params$type <- match.arg(params$type, c("area", "line", "heatmap"))
    params
  },
  draw_key = function(data, params, size) {
    # Determine the type from params (set in setup_params)
    type <- params$type %||% "area"

    if (type == "line") {
      # For line type: draw a line using fill as the colour
      # This matches the logic in draw_panel where fill is converted to line colour
      data$colour <- data$fill %||% "steelblue"
      ggplot2::draw_key_path(data, params, size)
    } else {
      # For area/heatmap type: draw a filled rectangle
      ggplot2::draw_key_polygon(data, params, size)
    }
  },
  draw_panel = function(
    data,
    panel_params,
    coord,
    type = "area",
    area_border = TRUE,
    na.rm = FALSE
  ) {
    # Transform start/end to xmin/xmax for drawing
    data$xmin <- data$start
    data$xmax <- data$end

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
      if (n == 0) {
        return(grid::nullGrob())
      }

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
        if (gn == 0) {
          return(NULL)
        }

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
        line_colour <- if (
          "fill" %in%
            names(gdata) &&
            length(gdata$fill) > 0 &&
            !all(is.na(gdata$fill))
        ) {
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
      # Use GeomRect for area type
      if (area_border) {
        # Bake alpha into fill & border colours directly, then draw at full
        # opacity so the border never double-composites over the fill.
        # This eliminates white-line artifacts between adjacent rectangles.
        a <- ifelse(is.na(data$alpha), 1, data$alpha)
        data$fill   <- scales::alpha(data$fill, a)
        data$colour <- scales::alpha(data$fill, 1)
        data$alpha   <- NA
        data$linewidth <- 0.1
      } else {
        # No border - may show rendering artifacts on some devices
        data$colour <- NA
      }
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
