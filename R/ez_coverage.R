# ezGenomeTracks - Coverage functions (split from ez_wrappers.R)
# Internal helper to resolve region from either region string or gene name
# This avoids code duplication across ez_* functions
#' @noRd
.resolve_region <- function(
  region = NULL,
  gene = NULL,
  gene_db = NULL,
  org_db = NULL,
  extend = 0.1,
  extend_type = c("proportion", "bp")
) {
  extend_type <- match.arg(extend_type)

  # Case 1: region provided directly
  if (!is.null(region)) {
    if (!is.null(gene)) {
      warning(
        "Both 'region' and 'gene' provided. Using 'region' and ignoring 'gene'."
      )
    }
    return(region)
  }

  # Case 2: gene name provided - look up coordinates
  if (!is.null(gene)) {
    if (is.null(gene_db)) {
      stop(
        "When using 'gene' parameter, 'gene_db' (TxDb object) must be provided ",
        "for coordinate lookup."
      )
    }
    return(gene_to_region(
      gene_name = gene,
      txdb = gene_db,
      org_db = org_db,
      extend = extend,
      extend_type = extend_type
    ))
  }

  # Case 3: neither provided
  stop("Either 'region' or 'gene' (with 'gene_db') must be provided.")
}


#' Easy coverage track visualization
#'
#' This function creates a coverage track visualization from various input types.
#' It provides a flexible interface with support for grouping and multiple tracks.
#'
#' @param input A GRanges object, data frame, character vector of file paths, or named list of data sources.
#' @param region Genomic region to display (e.g., "chr1:1000000-2000000").
#'   Either `region` or `gene` (with `gene_db`) must be provided.
#' @param gene Gene name/symbol to look up (e.g., "PTPRC", "TP53").
#'   When provided, the region is automatically determined from the gene coordinates
#'   in `gene_db`. Either `region` or `gene` must be provided.
#' @param gene_db TxDb object for gene coordinate lookup when using `gene` parameter
#'   (e.g., TxDb.Hsapiens.UCSC.hg38.knownGene).
#' @param org_db Optional OrgDb object for gene symbol mapping. If NULL (default),
#'   auto-detects available OrgDb packages.
#' @param extend Numeric. Amount to extend the region beyond the gene body when
#'   using `gene` parameter. Default: 0.1 (10% of gene length on each side).
#' @param extend_type How to interpret `extend`: "proportion" (relative to gene
#'   length) or "bp" (absolute base pairs). Default: "proportion".
#' @param track_labels Optional vector of track labels (used for character vector input)
#' @param type Type of signal visualization: "line", "area", or "heatmap" (default: "area")
#' @param group_var Column name for grouping data within a single data frame (default: NULL)
#' @param color_by Whether colors distinguish "group" or "track" (default: "group")
#' @param colors Color(s) for the coverage track. Can be a single color or
#'   a vector of colors for multiple tracks/groups (e.g., c("blue", "orange", "green")).
#'   If fewer colors than tracks/groups are provided, colors will be recycled.
#'   When NULL (default), uses vibeColors palette if available, otherwise "steelblue".
#'   For multiple tracks with a single color, automatically uses a colorblind-safe palette.
#' @param y_axis_style Y-axis style: "none", "simple", "minmax", or "full" (default: "none")
#'   - "none": No y-axis displayed
#'   - "simple": Shows y-range as \[min - max\] label at top-left
#'   - "minmax": Shows only min and max values on y-axis with ticks
#'   - "full": Full y-axis with all ticks and labels
#' @param y_range Y-axis range limits. When NULL (default), uses the global maximum
#'   across all tracks so all tracks share the same y-scale.
#' @param alpha Transparency (default: 0.5)
#' @param bin_width Width of bins in base pairs (default: NULL)
#' @param area_border Logical; if `TRUE` (default), draws thin borders on area rectangles
#'   to eliminate white-line rendering artifacts. Only affects `type = "area"`.
#' @param facet_label_position Position of facet labels: "top" or "left" (default: "top")
#' @param average Logical. If `TRUE`, averages overlapping tracks into a single
#'   track before plotting. Applies when `input` is a character vector of
#'   multiple file paths (overlapping tracks), or when `input` is a named list
#'   whose elements are character vectors with multiple files (averages within
#'   each list element separately, keeping tracks independent).
#'   Uses [average_coverage()] internally. (default: FALSE)
#' @param summary_fun Summary function used when `average = TRUE`. One of
#'   `"mean"`, `"median"`, `"max"`, `"min"`, `"sum"`. (default: `"mean"`)
#' @param n_bins Maximum number of bins used when reading BigWig
#'   files via megadepth (default: 2000). Increase for finer detail, decrease
#'   for faster plotting.
#' @param border Logical. If `TRUE`, adds a black border around the plotting panel (default: FALSE)
#' @param show_legend Logical. If `TRUE`, displays the legend (default: FALSE)
#' @param label_chr Logical. If `TRUE` (default), labels the x-axis with the chromosome name
#'   (e.g., "Chr1"). Set to `FALSE` to suppress the x-axis label.
#' @param ... Additional arguments passed to geom_coverage
#' @return A ggplot2 object
#' @export
#' @importFrom ggplot2 ggplot aes scale_y_continuous coord_cartesian labs facet_wrap scale_color_manual scale_fill_manual theme element_text
#' @importFrom dplyr filter mutate bind_rows
#' @examples
#' # From a GRanges object
#' library(GenomicRanges)
#' gr <- GRanges(
#'   seqnames = "chr1",
#'   ranges = IRanges(start = 1:100, end = 1:100),
#'   score = rnorm(100)
#' )
#' ez_coverage(gr, "chr1:1-100")
#'
#' # Single data frame with grouping
#' df <- data.frame(
#'   seqnames = "chr1", start = 1:100, end = 1:100,
#'   score = rnorm(100), sample = rep(c("A", "B"), 50)
#' )
#' ez_coverage(df, "chr1:1-100", group_var = "sample", colors = c("blue", "orange"))
#'
#' \dontrun{
#' # Using gene name instead of coordinates
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' ez_coverage(signal_data, gene = "PTPRC", gene_db = TxDb.Hsapiens.UCSC.hg38.knownGene)
#' }
ez_coverage <- function(
  input,
  region = NULL,
  gene = NULL,
  gene_db = NULL,
  org_db = NULL,
  extend = 0.1,
  extend_type = c("proportion", "bp"),
  track_labels = NULL,
  type = c("area", "line", "heatmap"),
  group_var = NULL,
  color_by = c("group", "track"),
  colors = NULL,
  y_axis_style = c("none", "simple", "minmax", "full"),
  y_range = NULL,
  alpha = 0.5,
  bin_width = NULL,
  area_border = TRUE,
  facet_label_position = c("top", "left"),
  border = FALSE,
  show_legend = FALSE,
  label_chr = TRUE,
  average = FALSE,
  summary_fun = c("mean", "median", "max", "min", "sum"),
  n_bins = 2000L,
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
  type <- match.arg(type)
  y_axis_style <- match.arg(y_axis_style)
  color_by <- match.arg(color_by)
  facet_label_position <- match.arg(facet_label_position)

  summary_fun <- match.arg(summary_fun)

  stopifnot(
    "alpha must be between 0 and 1" = alpha >= 0 && alpha <= 1,
    "bin_width must be positive integer" = is.null(bin_width) ||
      (bin_width > 0 && is.numeric(bin_width)),
    "n_bins must be a positive integer" =
      is.numeric(n_bins) &&
      length(n_bins) == 1 &&
      n_bins >= 1
  )

  # Resolve default colors
  if (is.null(colors)) {
    colors <- ez_default_single_color("coverage")
  }

  chr <- stringr::str_remove(stringr::str_split(region, ":")[[1]][1], "chr")

  # --- Average signal across samples if requested ---
  if (average) {
    if (is.character(input) && length(input) > 1) {
      # Character vector of multiple files (overlapping tracks): average into one
      plotDt <- average_coverage(
        inputs = input,
        region = region,
        n_bins = n_bins,
        summary_fun = summary_fun
      )
    } else if (is.list(input) && !is.data.frame(input)) {
      # Named list: each element is an independent track.
      # Average only within elements that are character vectors with >1 file.
      track_data_list <- list()
      for (i in seq_along(input)) {
        # Priority: track_labels > names(input) > default "Track i"
        if (!is.null(track_labels) && length(track_labels) >= i) {
          track_name <- track_labels[i]
        } else if (!is.null(names(input)[i]) && names(input)[i] != "") {
          track_name <- names(input)[i]
        } else {
          track_name <- paste0("Track ", i)
        }
        track_element <- input[[i]]

        if (is.character(track_element) && length(track_element) > 1) {
          avg_df <- average_coverage(
            inputs = track_element,
            region = region,
            n_bins = n_bins,
            summary_fun = summary_fun
          )
          avg_df$track <- track_name
          track_data_list[[i]] <- avg_df
        } else {
          # Single file, data frame, or GRanges — process normally
          single_df <- process_coverage_input(
            track_element,
            region,
            n_bins = n_bins
          )
          single_df$track <- track_name
          track_data_list[[i]] <- single_df
        }
      }
      plotDt <- dplyr::bind_rows(track_data_list)
    } else {
      warning("average = TRUE has no effect for a single input. Proceeding normally.")
      plotDt <- process_coverage_input(
        input,
        region,
        track_labels,
        n_bins = n_bins
      )
    }
  } else {
    # Process input using helper function (handles GRanges, data.frame, character, list)
    plotDt <- process_coverage_input(
      input,
      region,
      track_labels,
      n_bins = n_bins
    )
  }

  # Determine plotting strategy
  has_track <- "track" %in% colnames(plotDt)

  # Auto-detect group column when color_by = "group" and group_var not specified

  # This handles the case where process_coverage_input adds a "group" column
  # for character vector inputs (multiple files → overlapping tracks)
  if (
    is.null(group_var) && color_by == "group" && "group" %in% colnames(plotDt)
  ) {
    group_var <- "group"
  }

  has_group <- !is.null(group_var) && group_var %in% colnames(plotDt)

  # Create base plot
  # Always use 'fill' aesthetic - GeomCoverage handles converting fill to line colour for type="line"
  if (has_group) {
    # Data has grouping - use fill mapping
    if (has_track) {
      # Multiple tracks with grouping
      if (color_by == "group") {
        aes_mapping <- ggplot2::aes(
          start = start,
          end = end,
          ymin = 0,
          ymax = score,
          fill = .data[[group_var]]
        )
        color_values <- unique(plotDt[[group_var]])
        legend_name <- group_var
      } else {
        aes_mapping <- ggplot2::aes(
          start = start,
          end = end,
          ymin = 0,
          ymax = score,
          fill = track
        )
        color_values <- unique(plotDt$track)
        legend_name <- "Track"
      }
    } else {
      # Single track with grouping
      aes_mapping <- ggplot2::aes(
        start = start,
        end = end,
        ymin = 0,
        ymax = score,
        fill = .data[[group_var]]
      )
      color_values <- unique(plotDt[[group_var]])
      legend_name <- group_var
    }

    p <- ggplot2::ggplot(plotDt, aes_mapping) +
      geom_coverage(type = type, area_border = area_border, alpha = alpha, ...)

    # Apply color scales
    plot_colors <- resolve_plot_colors(colors, color_values)

    p <- p +
      ggplot2::scale_fill_manual(
        values = plot_colors,
        name = legend_name,
        guide = if (show_legend) "legend" else "none"
      )
  } else {
    # No grouping - use single color or track-based colors
    if (has_track) {
      # Multiple tracks without grouping
      aes_mapping <- ggplot2::aes(
        start = start,
        end = end,
        ymin = 0,
        ymax = score,
        fill = track
      )
      color_values <- unique(plotDt$track)
      legend_name <- "Track"

      p <- ggplot2::ggplot(plotDt, aes_mapping) +
        geom_coverage(type = type, area_border = area_border, alpha = alpha, ...)

      # Apply color scales
      plot_colors <- resolve_plot_colors(colors, color_values)

      p <- p +
        ggplot2::scale_fill_manual(
          values = plot_colors,
          name = legend_name,
          guide = if (show_legend) "legend" else "none"
        )
    } else {
      # Single track without grouping
      p <- ggplot2::ggplot(
        plotDt,
        ggplot2::aes(start = start, end = end, ymin = 0, ymax = score)
      ) +
        geom_coverage(
          type = type,
          area_border = area_border,
          fill = colors[1],
          alpha = alpha,
          ...
        )
    }
  }

  # Add faceting if multiple tracks
  if (has_track) {
    p <- apply_facet_position(p, facet_label_position, strip_placement = "inside")
  }

  # Apply binning if requested
  if (!is.null(bin_width)) {
    p <- p + stat_bin_signal(binwidth = bin_width)
  }

  # Parse region for x positioning
  region_gr <- parse_region(region)
  x_min <- GenomicRanges::start(region_gr)
  x_max <- GenomicRanges::end(region_gr)

  # Keep the x-axis baseline independent from the track fill colour.
  p <- p + ggplot2::geom_hline(yintercept = 0, colour = "black", linewidth = 0.2)

  # Calculate y-axis limits for annotations (track-specific if multiple tracks)
  if (is.null(y_range)) {
    y_min <- 0
    y_max <- max(plotDt$score, na.rm = TRUE)
    if (has_track) {
      # Use the global maximum across all tracks as the shared y-range
      y_range_df <- data.frame(
        track = unique(plotDt$track),
        y_min = y_min,
        y_max = y_max,
        y_label = paste0("[", round(y_min, 1), " - ", round(y_max, 1), "]"),
        x = x_min
      )
    }
  } else {
    y_min <- y_range[1]
    y_max <- y_range[2]
    if (has_track) {
      # Use the same y_range for all tracks
      y_range_df <- data.frame(
        track = unique(plotDt$track),
        y_min = y_min,
        y_max = y_max,
        y_label = paste0("[", round(y_min, 1), " - ", round(y_max, 1), "]"),
        x = x_min
      )
    }
  }

  # Apply the appropriate theme and scale
  if (y_axis_style == "minmax") {
    if (has_track) {
      # For multiple tracks with free_y scales, add per-track min/max labels
      # Create label data for both min and max values
      x_offset <- (x_max - x_min) * 0.025
      y_labels_df <- y_range_df |>
        dplyr::mutate(
          min_label = round(y_min, 1),
          max_label = round(y_max, 1),
          x_label = x_min - x_offset
        )
      p <- p +
        ez_coverage_theme(y_axis_style = y_axis_style) +
        scale_x_genome_region(region, oob = scales::oob_keep) +
        ggplot2::scale_y_continuous(
          expand = c(0, 0),
          breaks = c(y_min, y_max)
        ) +
        ggplot2::coord_cartesian(ylim = c(y_min, y_max), clip = "off") +
        ggplot2::labs(x = if (label_chr) paste0("Chr", chr) else NULL) +
        # Max label at top
        ggplot2::geom_text(
          data = y_labels_df,
          ggplot2::aes(
            x = .data$x_label,
            y = .data$y_max,
            label = .data$max_label
          ),
          hjust = 0.5,
          vjust = 0.5,
          size = 2.5,
          inherit.aes = FALSE
        ) +
        # Min label at bottom
        ggplot2::geom_text(
          data = y_labels_df,
          ggplot2::aes(
            x = .data$x_label,
            y = .data$y_min,
            label = .data$min_label
          ),
          hjust = 0.5,
          vjust = 0.5,
          size = 2.5,
          inherit.aes = FALSE
        )
    } else {
      # Single track - use annotate for min/max labels
      x_offset <- (x_max - x_min) * 0.025
      p <- p +
        ez_coverage_theme(y_axis_style = y_axis_style) +
        scale_x_genome_region(region, oob = scales::oob_keep) +
        ggplot2::scale_y_continuous(
          expand = c(0, 0),
          breaks = c(y_min, y_max)
        ) +
        ggplot2::coord_cartesian(ylim = c(y_min, y_max), clip = "off") +
        ggplot2::labs(x = if (label_chr) paste0("Chr", chr) else NULL) +
        # Max label at top
        ggplot2::annotate(
          "text",
          x = x_min - x_offset,
          y = y_max,
          label = round(y_max, 1),
          hjust = 0.5,
          vjust = 0.5,
          size = 2.5
        ) +
        # Min label at bottom
        ggplot2::annotate(
          "text",
          x = x_min - x_offset,
          y = y_min,
          label = round(y_min, 1),
          hjust = 0.5,
          vjust = 0.5,
          size = 2.5
        )
    }
  } else if (y_axis_style == "simple") {
    if (has_track) {
      # Use geom_text for per-track labels
      p <- p +
        ez_coverage_theme(y_axis_style = y_axis_style) +
        scale_x_genome_region(region) +
        ggplot2::scale_y_continuous(expand = c(0, 0)) +
        ggplot2::coord_cartesian(ylim = c(y_min, y_max), clip = "off") +
        ggplot2::labs(x = if (label_chr) paste0("Chr", chr) else NULL) +
        ggplot2::geom_text(
          data = y_range_df,
          ggplot2::aes(x = .data$x, y = .data$y_max, label = .data$y_label),
          hjust = 0,
          vjust = 1,
          size = 3,
          nudge_x = (x_max - x_min) * 0.01,
          nudge_y = -(max(y_range_df$y_max) * 0.05),
          inherit.aes = FALSE
        )
    } else {
      # Single track - use annotate
      y_label <- paste0("[", round(y_min, 1), " - ", round(y_max, 1), "]")
      p <- p +
        ez_coverage_theme(y_axis_style = y_axis_style) +
        scale_x_genome_region(region) +
        ggplot2::scale_y_continuous(expand = c(0, 0)) +
        ggplot2::coord_cartesian(ylim = c(y_min, y_max), clip = "off") +
        ggplot2::labs(x = if (label_chr) paste0("Chr", chr) else NULL) +
        ggplot2::annotate(
          "text",
          x = x_min + (x_max - x_min) * 0.01,
          y = y_max * 0.95,
          label = y_label,
          hjust = 0,
          vjust = 0,
          size = 3
        )
    }
  } else {
    p <- p +
      ez_coverage_theme(y_axis_style = y_axis_style) +
      scale_x_genome_region(region) +
      ggplot2::scale_y_continuous(expand = c(0, 0)) +
      ggplot2::coord_cartesian(ylim = c(y_min, y_max)) +
      ggplot2::labs(x = if (label_chr) paste0("Chr", chr) else NULL)
  }

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
