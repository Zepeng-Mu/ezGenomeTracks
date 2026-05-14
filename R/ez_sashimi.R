# ezGenomeTracks - Sashimi functions (split from ez_wrappers.R)
#' Easy sashimi plot visualization
#'
#' This function creates a sashimi plot combining coverage tracks with splice junction arcs.
#' It is commonly used for RNA-seq data to visualize both read coverage and splicing events.
#' Supports the same flexible input types as \code{\link{ez_coverage}}.
#'
#' @param coverage_data Coverage input: a GRanges object, data frame, character vector
#'   of file paths, or named list of data sources. Same input types as \code{\link{ez_coverage}}.
#'   When a named list is provided, each element becomes a separate faceted track.
#' @param junction_data Junction input that must match the structure of coverage_data:
#'   - If coverage_data is a single source (GRanges, data.frame, or file path),
#'     junction_data must also be a single source.
#'   - If coverage_data is a named list, junction_data must be a named list with
#'     the same names.
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
#' @param track_labels Optional vector of track labels (used for unnamed list input)
#' @param colors Color(s) for coverage fill and junction arcs. Can be a single color
#'   or a vector of colors for multiple tracks. If fewer colors than tracks are
#'   provided, colors will be recycled. When NULL (default), uses vibeColors palette
#'   if available, otherwise a built-in default.
#' @param coverage_fill Deprecated. Use `colors` instead.
#' @param junction_direction Direction of junction arcs: "alternate" (default, odd junctions up,
#'   even down), "up" (all arcs above coverage), or "down" (all arcs below zero line)
#' @param junction_curvature Curvature of junction arcs (0-1). Higher values create more
#'   pronounced curves. Default: 0.05
#' @param height_factor Height of arcs as proportion of genomic distance span.
#'   Higher values create taller arcs. Default: 0.05
#' @param alpha Transparency for both coverage and junction arcs (0-1). Default: 1
#' @param score_transform Transformation for junction scores when mapping to line width:
#'   "identity" (no transformation), "log10", or "sqrt". Default: "identity"
#' @param linewidth_range Range of line widths for junction arcs, mapped from score.
#'   Default: c(0.1, 1.5)
#' @param show_labels Logical, whether to show score labels at arc centers. Default: TRUE
#' @param label_size Font size for score labels. Default: 3
#' @param label_color Color for score labels. Default: "black"
#' @param y_axis_style Y-axis style: "none", "simple", or "full". Default: "none"
#' @param facet_label_position Position of facet labels: "top" or "left" (default: "top").
#'   Only applies when using multiple tracks via named list input.
#' @param border Logical. If `TRUE`, adds a black border around the plotting panel (default: FALSE)
#' @param show_legend Logical. If `TRUE`, displays the legend (default: FALSE)
#' @param ... Additional arguments passed to geom_coverage
#'
#' @return A ggplot2 object representing the sashimi plot
#'
#' @details
#' Sashimi plots combine two visual elements:
#' \itemize{
#'   \item Coverage track: Shows read depth across the region as a filled area
#'   \item Junction arcs: Connect splice donor and acceptor sites, with line width
#'     proportional to the number of supporting reads (score)
#' }
#'
#' For "up" direction, arc endpoints are positioned at the coverage height at those
#' positions. For "down" direction, arcs extend below the zero line. The "alternate"
#' mode assigns alternating directions to junctions sorted by start position, which
#' helps reduce visual overlap.
#'
#' Junction scores are mapped to line width using the specified transformation.
#' Use "log10" or "sqrt" for data with wide score ranges to prevent very thick lines.
#'
#' For multi-track visualization, provide named lists for both coverage_data and
#' junction_data with matching names. Each track will be displayed as a separate
#' facet with its own coverage and junction arcs.
#'
#' @export
#' @importFrom ggplot2 ggplot aes geom_text scale_linewidth_continuous coord_cartesian labs facet_wrap theme element_text unit scale_fill_manual element_rect
#' @importFrom dplyr mutate group_by summarise
#'
#' @examples
#' # Basic sashimi plot
#' coverage <- data.frame(
#'   seqnames = "chr1", start = 1:100, end = 2:101,
#'   score = c(runif(30, 5, 10), rep(0, 20), runif(50, 5, 10))
#' )
#' junctions <- data.frame(
#'   seqnames = "chr1", start = c(30), end = c(50), score = c(15)
#' )
#' ez_sashimi(coverage, junctions, "chr1:1-100")
#'
#' # Multi-track sashimi plot with named lists
#' cov_list <- list("Sample1" = coverage1, "Sample2" = coverage2)
#' junc_list <- list("Sample1" = junctions1, "Sample2" = junctions2)
#' ez_sashimi(cov_list, junc_list, "chr1:1-100", colors = c("purple", "orange"))
ez_sashimi <- function(
  coverage_data,
  junction_data,
  region = NULL,
  gene = NULL,
  gene_db = NULL,
  org_db = NULL,
  extend = 0.1,
  extend_type = c("proportion", "bp"),
  track_labels = NULL,
  colors = NULL,
  coverage_fill = NULL,
  junction_direction = c("alternate", "up", "down"),
  junction_curvature = 0.05,
  height_factor = 0.05,
  alpha = 1,
  score_transform = c("identity", "log10", "sqrt"),
  linewidth_range = c(0.1, 1.5),
  show_labels = TRUE,
  label_size = 3,
  label_color = "black",
  y_axis_style = c("none", "simple", "full"),
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
  junction_direction <- match.arg(junction_direction)
  score_transform <- match.arg(score_transform)
  y_axis_style <- match.arg(y_axis_style)
  facet_label_position <- match.arg(facet_label_position)

  stopifnot(
    "alpha must be between 0 and 1" = alpha >= 0 && alpha <= 1,
    "linewidth_range must be a numeric vector of length 2" = is.numeric(
      linewidth_range
    ) &&
      length(linewidth_range) == 2
  )

  # Resolve default colors
  if (is.null(colors)) {
    colors <- ez_default_single_color("sashimi")
  }

  # Handle deprecated coverage_fill parameter
  if (!is.null(coverage_fill)) {
    warning(
      "coverage_fill is deprecated. Use 'colors' instead.",
      call. = FALSE
    )
    colors <- coverage_fill
  }

  chr <- stringr::str_remove(stringr::str_split(region, ":")[[1]][1], "chr")

  # Process data using new helper function
  processed <- process_sashimi_input(
    coverage_data = coverage_data,
    junction_data = junction_data,
    region = region,
    track_labels = track_labels,
    junction_direction = junction_direction
  )

  coverage_df <- processed$coverage_df
  junction_df <- processed$junction_df

  # Determine if we have multiple tracks
  has_track <- "track" %in% colnames(coverage_df)

  # Set up colors for tracks
  if (has_track) {
    track_names <- unique(coverage_df$track)
    plot_colors <- resolve_plot_colors(colors, track_names)
  } else {
    plot_colors <- colors[1]
  }

  # Apply score transformation for linewidth mapping
  if (nrow(junction_df) > 0) {
    junction_df <- junction_df |>
      dplyr::mutate(
        score_transformed = switch(
          score_transform,
          "identity" = score,
          "log10" = log10(score + 1),
          "sqrt" = sqrt(score)
        )
      )
  }

  # Build plot differently based on whether we have tracks
  if (has_track) {
    # Multi-track mode: use fill aesthetic mapped to track
    p <- ggplot2::ggplot() +
      geom_coverage(
        data = coverage_df,
        ggplot2::aes(
          start = start,
          end = end,
          ymin = 0,
          ymax = score,
          fill = track
        ),
        alpha = alpha,
        ...
      ) +
      ggplot2::scale_fill_manual(
        values = plot_colors,
        name = "Track",
        guide = if (show_legend) "legend" else "none"
      )

    # Add junction arcs for each track with matching colors
    if (nrow(junction_df) > 0) {
      for (track_name in track_names) {
        track_junctions <- junction_df[junction_df$track == track_name, ]
        track_color <- plot_colors[track_name]

        junctions_up <- track_junctions[track_junctions$arc_direction == "up", ]
        junctions_down <- track_junctions[
          track_junctions$arc_direction == "down",
        ]

        if (nrow(junctions_up) > 0) {
          p <- p +
            geom_link(
              data = junctions_up,
              ggplot2::aes(
                start1 = start1,
                end1 = end1,
                start2 = start2,
                end2 = end2,
                y = y_start,
                yend = y_end,
                linewidth = score_transformed
              ),
              curvature = junction_curvature,
              height_factor = height_factor,
              direction = "up",
              color = track_color,
              alpha = alpha
            )
        }

        if (nrow(junctions_down) > 0) {
          p <- p +
            geom_link(
              data = junctions_down,
              ggplot2::aes(
                start1 = start1,
                end1 = end1,
                start2 = start2,
                end2 = end2,
                y = y_start,
                yend = y_end,
                linewidth = score_transformed
              ),
              curvature = junction_curvature,
              height_factor = height_factor,
              direction = "down",
              color = track_color,
              alpha = alpha
            )
        }
      }
    }

    # Add faceting for multiple tracks
    p <- apply_facet_position(p, facet_label_position)
  } else {
    # Single track mode
    p <- ggplot2::ggplot() +
      geom_coverage(
        data = coverage_df,
        ggplot2::aes(start = start, end = end, ymin = 0, ymax = score),
        fill = plot_colors,
        alpha = alpha,
        ...
      )

    # Add junction arcs
    if (nrow(junction_df) > 0) {
      junctions_up <- junction_df[junction_df$arc_direction == "up", ]
      junctions_down <- junction_df[junction_df$arc_direction == "down", ]

      if (nrow(junctions_up) > 0) {
        p <- p +
          geom_link(
            data = junctions_up,
            ggplot2::aes(
              start1 = start1,
              end1 = end1,
              start2 = start2,
              end2 = end2,
              y = y_start,
              yend = y_end,
              linewidth = score_transformed
            ),
            curvature = junction_curvature,
            height_factor = height_factor,
            direction = "up",
            color = plot_colors,
            alpha = alpha
          )
      }

      if (nrow(junctions_down) > 0) {
        p <- p +
          geom_link(
            data = junctions_down,
            ggplot2::aes(
              start1 = start1,
              end1 = end1,
              start2 = start2,
              end2 = end2,
              y = y_start,
              yend = y_end,
              linewidth = score_transformed
            ),
            curvature = junction_curvature,
            height_factor = height_factor,
            direction = "down",
            color = plot_colors,
            alpha = alpha
          )
      }
    }
  }

  # Add linewidth scale if there are junctions
  if (nrow(junction_df) > 0) {
    p <- p +
      ggplot2::scale_linewidth_continuous(
        range = linewidth_range,
        name = if (score_transform == "identity") {
          "Junction\nReads"
        } else {
          paste0("Junction\nReads\n(", score_transform, ")")
        },
        guide = if (show_legend) "legend" else "none"
      )

    # Add score labels at arc centers if requested
    if (show_labels) {
      junction_df <- junction_df |>
        dplyr::mutate(
          label_x = (start1 + start2) / 2,
          arc_span = abs(start2 - start1),
          arc_peak = arc_span * height_factor,
          label_y = ifelse(
            arc_direction == "up",
            (y_start + y_end) / 2 + arc_peak,
            -arc_peak
          ),
          label_vjust = ifelse(arc_direction == "up", -0.3, 1.3)
        )

      p <- p +
        ggplot2::geom_text(
          data = junction_df,
          ggplot2::aes(x = label_x, y = label_y, label = score),
          color = label_color,
          size = label_size,
          vjust = junction_df$label_vjust
        )
    }
  }

  # Calculate y-axis limits to accommodate both coverage and arcs
  max_coverage <- max(coverage_df$score, na.rm = TRUE)
  if (nrow(junction_df) > 0) {
    max_arc_height <- max(
      abs(junction_df$start2 - junction_df$start1) * height_factor,
      na.rm = TRUE
    )
    has_up_arcs <- any(junction_df$arc_direction == "up")
    has_down_arcs <- any(junction_df$arc_direction == "down")

    y_upper <- max_coverage +
      if (has_up_arcs) max_arc_height * 1.5 else max_coverage * 0.1
    y_lower <- if (has_down_arcs) -(max_arc_height * 1.5) else 0
  } else {
    y_upper <- max_coverage * 1.1
    y_lower <- 0
  }

  # Apply theme and scales
  p <- p +
    ez_sashimi_theme(y_axis_style = y_axis_style) +
    scale_x_genome_region(region) +
    ggplot2::coord_cartesian(ylim = c(y_lower, y_upper), clip = "off") +
    ggplot2::labs(x = paste0("Chr", chr))

  # Apply border after theme
  if (border) {
    p <- apply_border_theme(p)
  }

  # Remove spacing between facet panels (applied last)
  if (has_track) {
    p <- p + ggplot2::theme(panel.spacing.y = ggplot2::unit(0, "pt"))
  }

  return(p)
}
