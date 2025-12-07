#' A minimal theme for genome track visualization
#'
#' This function creates a minimal, clean theme for genome track visualization.
#' It removes most of the background elements and grid lines to focus on the data.
#'
#' @param base_size Base font size
#' @param base_family Base font family
#' @param base_line_size Base line size
#' @param base_rect_size Base rect size
#' @param show_grid Show grid lines (default: FALSE)
#' @param show_ticks Show axis ticks (default: TRUE)
#' @param show_x_axis Show x-axis line (default: TRUE)
#' @param show_y_axis Show y-axis line (default: FALSE)
#' @return A ggplot2 theme object
#' @export
#' @importFrom ggplot2 theme element_blank element_line element_text margin
#' @examples
#' \dontrun{
#' library(ggplot2)
#' p <- ggplot(data, aes(x = start, y = score)) +
#'   geom_line() +
#'   ez_theme()
#' }
ez_theme <- function(base_size = 10, base_family = "", base_line_size = 0.2,
                     base_rect_size = 0.2, show_grid = FALSE,
                     show_ticks = TRUE, show_x_axis = TRUE, show_y_axis = TRUE) {
  # Start with a minimal theme
  theme <- ggplot2::theme_minimal(
    base_size = base_size, base_family = base_family,
    base_line_size = base_line_size, base_rect_size = base_rect_size
  )

  # Modify the theme to be even more minimal
  theme <- theme + ggplot2::theme(
    # Remove panel grid
    panel.grid.major = if (show_grid) ggplot2::element_line(color = "gray90", linewidth = 0.2) else ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),

    # Remove panel border
    panel.border = ggplot2::element_blank(),

    # Customize axis
    axis.line.x = if (show_x_axis) ggplot2::element_line(color = "black", linewidth = 0.2) else ggplot2::element_blank(),
    axis.line.y = if (show_y_axis) ggplot2::element_line(color = "black", linewidth = 0.2) else ggplot2::element_blank(),
    axis.ticks = if (show_ticks) ggplot2::element_line(color = "black", linewidth = 0.2) else ggplot2::element_blank(),
    axis.ticks.length = ggplot2::unit(2, "pt"),

    # Reduce plot margins
    plot.margin = ggplot2::margin(0, 5, 0, 5),

    # Customize text
    # axis.title = ggplot2::element_text(size = ggplot2::rel(0.8)),
    # axis.text = ggplot2::element_text(size = ggplot2::rel(0.7)),

    # Remove legend background
    legend.background = ggplot2::element_blank(),
    legend.key = ggplot2::element_blank()
  )

  return(theme)
}

#' A theme for coverage tracks
#'
#' This function creates a theme specifically for coverage/signal tracks.
#' It removes the y-axis text and title, and makes the plot more compact.
#'
#' @param y_axis_style Style of the y-axis. Options are "none", "simple", and "full".
#' @param ... Additional arguments passed to ez_theme
#' @return A ggplot2 theme object
#' @export
#' @importFrom ggplot2 theme element_blank
#' @examples
#' \dontrun{
#' library(ggplot2)
#' p <- ggplot(data, aes(x = start, y = score)) +
#'   geom_line() +
#'   ez_coverage_theme()
#' }
ez_coverage_theme <- function(y_axis_style = c("none", "simple", "full"), ...) {
  y_axis_style <- match.arg(y_axis_style)

  if (y_axis_style == "none") {
    theme <- ez_theme(...) +
      ggplot2::theme(
        axis.title.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.line.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(5, 5, 5, 5)
      )
  } else if (y_axis_style == "simple") {
    theme <- ez_theme(...) +
      ggplot2::theme(
        axis.title.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.line.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(5, 5, 5, 5)
      )
  } else if (y_axis_style == "full") {
    theme <- ez_theme(...) +
      ggplot2::theme(
        axis.title.y = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(5, 5, 5, 5)
      )
  }
}

#' A theme for gene tracks
#'
#' This function creates a theme specifically for gene tracks.
#' It removes the y-axis text and title, and makes the plot more compact.
#'
#' @param ... Additional arguments passed to ez_theme
#' @return A ggplot2 theme object
#' @export
#' @importFrom ggplot2 theme element_blank
#' @examples
#' \dontrun{
#' library(ggplot2)
#' p <- ggplot(data, aes(x = start, y = gene_name)) +
#'   geom_gene() +
#'   ez_gene_theme()
#' }
ez_gene_theme <- function(...) {
  ez_theme(...) + ggplot2::theme(
    axis.text.y = ggplot2::element_text(size = ggplot2::rel(1)),
    axis.line.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(5, 5, 5, 5)
  )
}

#' A theme for feature tracks
#'
#' This function creates a theme specifically for feature tracks.
#' It removes the y-axis text and title, and makes the plot more compact.
#'
#' @param ... Additional arguments passed to ez_theme
#' @return A ggplot2 theme object
#' @export
#' @importFrom ggplot2 theme element_blank
#' @examples
#' \dontrun{
#' library(ggplot2)
#' p <- ggplot(data, aes(x = start, y = 1)) +
#'   geom_feature() +
#'   ez_feature_theme()
#' }
ez_feature_theme <- function(...) {
  ez_theme(...) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
    )
}

#' A theme for sashimi plots
#'
#' This function creates a theme specifically for sashimi plots that combine
#' coverage tracks with junction arcs. It provides extra vertical margin for
#' arc labels and handles y-axis display options.
#'
#' @param y_axis_style Style of the y-axis. Options are "none", "simple", and "full".
#'   Default: "none"
#' @param ... Additional arguments passed to ez_theme
#' @return A ggplot2 theme object
#' @export
#' @importFrom ggplot2 theme element_blank margin
#' @examples
#' \dontrun{
#' library(ggplot2)
#' p <- ggplot(data, aes(x = start, y = score)) +
#'   geom_area() +
#'   ez_sashimi_theme()
#' }
ez_sashimi_theme <- function(y_axis_style = c("none", "simple", "full"), ...) {
  y_axis_style <- match.arg(y_axis_style)

  if (y_axis_style == "none") {
    theme <- ez_theme(...) +
      ggplot2::theme(
        axis.title.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.line.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        # Extra margins for arc labels
        plot.margin = ggplot2::margin(10, 5, 10, 5)
      )
  } else if (y_axis_style == "simple") {
    theme <- ez_theme(...) +
      ggplot2::theme(
        axis.title.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.line.y = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(10, 5, 10, 5)
      )
  } else if (y_axis_style == "full") {
    theme <- ez_theme(...) +
      ggplot2::theme(
        axis.title.y = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(10, 5, 10, 5)
      )
  }

  return(theme)
}
