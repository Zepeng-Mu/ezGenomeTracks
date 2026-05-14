# ezGenomeTracks - Link functions (split from ez_wrappers.R)
#' Easy interaction track visualization
#'
#' This function creates an interaction track visualization from various input types.
#' It provides a flexible interface with support for grouping and multiple tracks,
#' similar to \code{\link{ez_coverage}}.
#'
#' @param input A GRanges, GInteractions, data frame, character vector of file paths,
#'   or named list of data sources with interaction data.
#' @param region Genomic region to display (e.g., "chr1:1000000-2000000").
#'   Either `region` or `gene` (with `gene_db`) must be provided.
#' @param gene Gene name/symbol to look up (e.g., "PTPRC", "TP53").
#'   When provided, the region is automatically determined from the gene coordinates
#'   in `gene_db`. Either `region` or `gene` must be provided.
#' @param gene_db TxDb object for gene coordinate lookup when using `gene` parameter.
#' @param org_db Optional OrgDb object for gene symbol mapping. If NULL (default),
#'   auto-detects available OrgDb packages.
#' @param extend Numeric. Amount to extend the region beyond the gene body when
#'   using `gene` parameter. Default: 0.1 (10% of gene length on each side).
#' @param extend_type How to interpret `extend`: "proportion" (relative to gene
#'   length) or "bp" (absolute base pairs). Default: "proportion".
#' @param track_labels Optional vector of track labels (used for character vector input)
#' @param group_var Column name for grouping data within a single data frame (default: NULL)
#' @param color_by Whether colors distinguish "group" or "track" (default: "group")
#' @param colors Color(s) for the arcs. Can be a single color (e.g., "gray50") or
#'   a vector of colors for multiple tracks/groups (e.g., c("blue", "orange", "green")).
#'   If fewer colors than tracks/groups are provided, colors will be recycled.
#'   Default is "gray50".
#' @param curvature Numeric value controlling the arc curvature (0-1).
#'   Higher values create more pronounced curves. Default: 0.5
#' @param height_factor Height of curves as proportion of genomic distance span.
#'   Higher values create taller arcs. Default: 0.15
#' @param direction Direction of curve arcs: "down" (negative y, default) or "up" (positive y).
#'   Default: "down"
#' @param size Line width of the arcs. Default: 0.5
#' @param alpha Transparency level of the arcs (0 = transparent, 1 = opaque).
#'   Default: 0.7
#' @param use_score Logical indicating whether to use the 'score' column for
#'   arc coloring. If TRUE, a color gradient will be applied based on the
#'   interaction scores. Default: FALSE
#' @param facet_label_position Position of facet labels: "top" or "left" (default: "top")
#' @param border Logical. If `TRUE`, adds a black border around the plotting panel (default: FALSE)
#' @param show_legend Logical. If `TRUE`, displays the legend (default: FALSE)
#' @param ... Additional arguments passed to `geom_link()`
#'
#' @return A ggplot2 object representing the link track.
#'
#' @details
#' The function automatically handles different input types:
#' \itemize{
#'   \item For BEDPE files: Reads and processes the interaction data
#'   \item For data frames: Expects columns for interaction coordinates (start1, end1, start2, end2)
#'     and optionally a 'score' column if `use_score = TRUE`
#'   \item For named lists: Creates multiple faceted tracks
#'   \item With group_var: Groups interactions within a single data frame for overlaid visualization
#' }
#'
#' The visualization includes:
#' \itemize{
#'   \item Arcs connecting interaction anchors
#'   \item Optional score-based coloring
#'   \item Support for grouping and faceting
#'   \item Automatic scaling to fit the specified genomic region
#' }
#'
#' @export
#' @importFrom ggplot2 ggplot aes scale_color_gradient scale_color_manual facet_wrap
#' @importFrom methods is
#'
#' @examples
#' # From a data frame with uniform coloring
#' data(example_interactions)
#' ez_link(
#'   example_interactions,
#'   "chr1:1-15000",
#'   colors = "darkblue",
#'   size = 1,
#'   alpha = 0.8
#' )
#'
#' # Single data frame with grouping
#' df <- data.frame(
#'   start1 = c(1000, 2000, 5000, 6000),
#'   end1 = c(1500, 2500, 5500, 6500),
#'   start2 = c(3000, 4000, 8000, 9000),
#'   end2 = c(3500, 4500, 8500, 9500),
#'   sample = c("A", "A", "B", "B")
#' )
#' ez_link(df, "chr1:1-10000", group_var = "sample", colors = c("blue", "orange"))
ez_link <- function(
  input,
  region = NULL,
  gene = NULL,
  gene_db = NULL,
  org_db = NULL,
  extend = 0.1,
  extend_type = c("proportion", "bp"),
  track_labels = NULL,
  group_var = NULL,
  color_by = c("group", "track"),
  colors = "gray50",
  curvature = 0.5,
  height_factor = 0.15,
  direction = c("down", "up"),
  size = 0.5,
  alpha = 0.7,
  use_score = FALSE,
  facet_label_position = c("top", "left"),
  border = FALSE,
  show_legend = FALSE,
  ...
) {
  # Resolve region from either region string or gene name
  region <- .resolve_region(
    region = region,
    gene = gene,
    gene_db = gene_db,
    org_db = org_db,
    extend = extend,
    extend_type = extend_type
  )

  # Validate inputs
  direction <- match.arg(direction)
  color_by <- match.arg(color_by)
  facet_label_position <- match.arg(facet_label_position)

  stopifnot(
    "alpha must be between 0 and 1" = alpha >= 0 && alpha <= 1
  )

  # Process input using helper function (handles GRanges, GInteractions, data.frame, character, list)
  plotDt <- process_interaction_input(input, region, track_labels)

  # Determine plotting strategy
  has_track <- "track" %in% colnames(plotDt)

  # Auto-detect group column when color_by = "group" and group_var not specified
  if (
    is.null(group_var) && color_by == "group" && "group" %in% colnames(plotDt)
  ) {
    group_var <- "group"
  }

  has_group <- !is.null(group_var) && group_var %in% colnames(plotDt)

  # If score is requested but not present, warn and disable
  if (use_score && !"score" %in% colnames(plotDt)) {
    warning("Score column not found, disabling score-based coloring")
    use_score <- FALSE
  }

  # Calculate y-axis limits based on maximum curve height
  y_limits <- calculate_link_ylim(plotDt, height_factor, direction)

  # Create base plot with appropriate aesthetics
  if (use_score) {
    # Score-based coloring takes precedence
    p <- ggplot2::ggplot(plotDt) +
      geom_link(
        ggplot2::aes(
          start1 = start1,
          end1 = end1,
          start2 = start2,
          end2 = end2,
          color = score
        ),
        curvature = curvature,
        height_factor = height_factor,
        direction = direction,
        linewidth = size,
        alpha = alpha,
        ...
      ) +
      ggplot2::scale_color_gradient(
        low = "gray80",
        high = colors[1],
        guide = if (show_legend) "colorbar" else "none"
      )
  } else if (has_group) {
    # Data has grouping - use color mapping
    if (has_track) {
      # Multiple tracks with grouping
      if (color_by == "group") {
        aes_mapping <- ggplot2::aes(
          start1 = start1,
          end1 = end1,
          start2 = start2,
          end2 = end2,
          color = .data[[group_var]]
        )
        color_values <- unique(plotDt[[group_var]])
        legend_name <- group_var
      } else {
        aes_mapping <- ggplot2::aes(
          start1 = start1,
          end1 = end1,
          start2 = start2,
          end2 = end2,
          color = track
        )
        color_values <- unique(plotDt$track)
        legend_name <- "Track"
      }
    } else {
      # Single track with grouping
      aes_mapping <- ggplot2::aes(
        start1 = start1,
        end1 = end1,
        start2 = start2,
        end2 = end2,
        color = .data[[group_var]]
      )
      color_values <- unique(plotDt[[group_var]])
      legend_name <- group_var
    }

    p <- ggplot2::ggplot(plotDt, aes_mapping) +
      geom_link(
        curvature = curvature,
        height_factor = height_factor,
        direction = direction,
        linewidth = size,
        alpha = alpha,
        ...
      )

    # Apply color scales
    plot_colors <- resolve_plot_colors(colors, color_values)

    p <- p +
      ggplot2::scale_color_manual(
        values = plot_colors,
        name = legend_name,
        guide = if (show_legend) "legend" else "none"
      )
  } else {
    # No grouping - use single color or track-based colors
    if (has_track) {
      # Multiple tracks without grouping
      aes_mapping <- ggplot2::aes(
        start1 = start1,
        end1 = end1,
        start2 = start2,
        end2 = end2,
        color = track
      )
      color_values <- unique(plotDt$track)
      legend_name <- "Track"

      p <- ggplot2::ggplot(plotDt, aes_mapping) +
        geom_link(
          curvature = curvature,
          height_factor = height_factor,
          direction = direction,
          linewidth = size,
          alpha = alpha,
          ...
        )

      # Apply color scales
      plot_colors <- resolve_plot_colors(colors, color_values)

      p <- p +
        ggplot2::scale_color_manual(
          values = plot_colors,
          name = legend_name,
          guide = if (show_legend) "legend" else "none"
        )
    } else {
      # Single track without grouping
      p <- ggplot2::ggplot(
        plotDt,
        ggplot2::aes(start1 = start1, end1 = end1, start2 = start2, end2 = end2)
      ) +
        geom_link(
          curvature = curvature,
          height_factor = height_factor,
          direction = direction,
          color = colors[1],
          linewidth = size,
          alpha = alpha,
          ...
        )
    }
  }

  # Add faceting if multiple tracks
  if (has_track) {
    p <- apply_facet_position(p, facet_label_position, strip_placement = "inside")
  }

  # Apply the appropriate theme and scale with automatic clipping prevention
  p <- p +
    scale_x_genome_region(region) +
    ggplot2::scale_y_continuous(limits = y_limits, expand = c(0, 0)) +
    ggplot2::coord_cartesian(clip = "off") +
    ez_theme() +
    ggplot2::theme(
      axis.line.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank()
    )

  # Apply border after theme (so it doesn't get overwritten)
  if (border) {
    p <- apply_border_theme(p)
  }

  # Apply facet label theming last (after all other themes, so it doesn't get overwritten)
  if (has_track) {
    if (facet_label_position == "left") {
      p <- p + ggplot2::theme(
        strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 1),
        strip.placement = "inside",
        panel.spacing.y = ggplot2::unit(0, "pt")
      )
    } else {
      p <- p + ggplot2::theme(panel.spacing.y = ggplot2::unit(0, "pt"))
    }
  }

  return(p)
}
