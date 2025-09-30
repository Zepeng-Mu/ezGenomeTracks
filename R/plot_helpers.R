#' Plot signal tracks from a data frame with grouping
#'
#' This function creates signal tracks from a data frame with grouping capabilities.
#' It can separate data into individual tracks using facets (`track_by`) and create
#' overlapping tracks within each facet using different colors (`group_by`).
#'
#' @param data A data frame containing signal data with columns: start, end, score,
#'   track_by (optional), and group_by (optional)
#' @param region Genomic region to display (e.g., "chr1:1000000-2000000")
#' @param track_by Column name in data for separating tracks into facets (optional)
#' @param group_by Column name in data for grouping overlapping tracks (optional)
#' @param color_by Column name in data for coloring tracks (defaults to group_by if not specified)
#' @param colors Vector of colors for different groups (default: NULL, uses default ggplot2 colors)
#' @param type Type of signal visualization: "line", "area", or "heatmap" (default: "area")
#' @param fill Fill color for area plots (default: "purple2", ignored if group_by is specified)
#' @param color Line color (default: "purple2", ignored if group_by is specified)
#' @param alpha Transparency (default: 0.8)
#' @param binwidth Width of bins in base pairs (default: NULL)
#' @param facet_scales Scale parameter for facet_wrap (default: "free_y")
#' @param ... Additional arguments passed to geom_signal
#' @return A ggplot2 object
#' @export
#' @importFrom ggplot2 ggplot aes facet_wrap scale_fill_manual scale_color_manual
#' @examples
#' \dontrun{
#' # Basic usage with track separation
#' p <- plot_signal_df(signal_data, "chr1:1000000-2000000", track_by = "sample")
#'
#' # With both track separation and grouping
#' p <- plot_signal_df(signal_data, "chr1:1000000-2000000",
#'   track_by = "sample", group_by = "condition"
#' )
#' }
plot_signal_df <- function(data, region, track_by = NULL, group_by = NULL,
                           color_by = NULL, colors = NULL, type = "area",
                           alpha = 0.8, binwidth = NULL,
                           facet_scales = "free_y", ...) {
  # Check required columns
  required_cols <- c("start", "end", "score")
  if (!all(required_cols %in% colnames(data))) {
    stop("Data must have columns: ", paste(required_cols, collapse = ", "))
  }

  # Create base plot
  p <- ggplot2::ggplot(data, ggplot2::aes(x = start, y = score))

  # Add group aesthetics if group_by is specified
  if (!is.null(group_by)) {
    # Create mapping with grouping
    if (!is.null(color_by)) {
      mapping <- ggplot2::aes(
        fill = .data[[color_by]],
        color = .data[[color_by]],
        group = .data[[group_by]]
      )
    } else {
      mapping <- ggplot2::aes(
        group = .data[[group_by]]
      )
    }

    # Add geom_signal with group mapping
    p <- p + geom_signal(mapping = mapping, type = type, alpha = alpha, ...)

    # Add color scales if colors are specified
    if (!is.null(colors)) {
      p <- p +
        ggplot2::scale_fill_manual(values = colors) +
        ggplot2::scale_color_manual(values = colors)
    }
  } else {
    # Add geom_signal without grouping
    p <- p + geom_signal(
      type = type,
      alpha = alpha, ...
    )
  }

  # Apply binning if requested
  if (!is.null(binwidth)) {
    p <- p + stat_bin_signal(binwidth = binwidth)
  }

  # Add faceting if track_by is specified
  if (!is.null(track_by)) {
    p <- p + ggplot2::facet_wrap(
      ggplot2::vars(.data[[track_by]]),
      ncol = 1,
      scales = facet_scales,
      strip.position = "left"
    )
  }

  # Apply the appropriate theme and scale
  p <- p + ez_signal_theme() +
    scale_x_genome_region(region)

  return(p)
}
