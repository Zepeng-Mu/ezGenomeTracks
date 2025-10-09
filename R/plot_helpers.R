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
  p <- ggplot2::ggplot(data)

  # Set default color_by to group_by if not specified
  if (is.null(color_by) && !is.null(group_by)) {
    color_by <- group_by
  }

  # Create the mapping
  if (!is.null(group_by)) {
    # Base mapping with grouping
    mapping_list <- list(
      x = rlang::expr(.data$start),
      y = rlang::expr(.data$score),
      group = rlang::expr(.data[[group_by]])
    )

    # Add color aesthetics if color_by is specified
    if (!is.null(color_by)) {
      mapping_list$fill <- rlang::expr(.data[[color_by]])
      mapping_list$colour <- rlang::expr(.data[[color_by]])
    }

    # Convert list to an aes object
    mapping <- do.call(ggplot2::aes, mapping_list)

    # Add geom_signal with the mapping
    p <- p + geom_coverage(
      mapping = mapping,
      type = type,
      alpha = alpha,
      ...
    )

    # Add color scales if colors are specified
    if (!is.null(colors)) {
      legend_name <- color_by
      p <- p +
        ggplot2::scale_fill_manual(values = colors, name = legend_name) +
        ggplot2::scale_colour_manual(values = colors, name = legend_name)
    }
  } else {
    # Add geom_signal without grouping or color mapping
    p <- p + geom_coverage(
      mapping = ggplot2::aes(x = .data$start, y = .data$score),
      type = type,
      alpha = alpha,
      ...
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
