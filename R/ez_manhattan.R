# ezGenomeTracks - Manhattan functions (genome-wide Manhattan plots)

#' Easy genome-wide Manhattan plot visualization
#'
#' This function creates a genome-wide Manhattan plot for GWAS results across
#' multiple chromosomes. It is a wrapper around geom_manhattan that provides
#' a flexible interface with support for grouping and multiple tracks.
#'
#' @description
#' This function creates a Manhattan plot from GWAS (Genome-Wide Association Study)
#' data, visualizing p-values across the entire genome with chromosomes displayed
#' along the x-axis using cumulative positions and alternating colors.
#'
#' For regional association plots (LocusZoom-style) focused on a single locus,
#' use `ez_locusZoom()` instead.
#'
#' @param input A data frame or named list of data frames containing GWAS results
#'   with columns for chromosome, position, p-values, and optionally SNP names.
#'   Must contain data from multiple chromosomes. Supports both GWAS-style
#'   (CHR, BP, P) and GRanges-style (seqnames, start, pvalue) column naming.
#' @param chr Character string specifying the column name for chromosome numbers.
#'   Default: auto-detect from "CHR", "seqnames", "chrom", etc.
#' @param bp Character string specifying the column name for base pair positions.
#'   Default: auto-detect from "BP", "start", "pos", etc.
#' @param p Character string specifying the column name for p-values.
#'   Default: auto-detect from "P", "pvalue", "p.value", etc.
#' @param snp Character string specifying the column name for SNP identifiers.
#'   Default: auto-detect from "SNP", "rsid", "variant_id", etc.
#' @param track_labels Optional vector of track labels (used for unnamed list input).
#'   Default: NULL.
#' @param group_var Column name for grouping data within a single data frame.
#'   Default: NULL.
#' @param logp Logical indicating whether to plot -log10(p-values).
#'   Default: TRUE.
#' @param size Numeric value for point size in the plot. Default: 0.5.
#' @param colors Vector of colors for alternating chromosomes.
#'   Default: `c("grey", "skyblue")`.
#' @param highlight_snps Character vector of SNP IDs to highlight.
#'   Default: NULL.
#' @param highlight_color Color for highlighting significant SNPs.
#'   Default: "purple".
#' @param threshold_p Numeric p-value threshold for drawing a significance line
#'   (e.g., 5e-8 for genome-wide significance). If NULL, no line is drawn.
#'   Default: NULL.
#' @param threshold_color Color for the significance threshold line.
#'   Default: "red".
#' @param threshold_linetype Linetype for the significance threshold line.
#'   Default: 2 (dashed).
#' @param facet_label_position Position of facet labels for multi-track plots:
#'   "top" or "left". Default: "top".
#' @param region Deprecated. Use `ez_locusZoom()` for regional plots.
#' @param gene Deprecated. Use `ez_locusZoom()` for regional plots.
#' @param ... Additional arguments passed to `geom_manhattan()`.
#'
#' @return A ggplot2 object containing the Manhattan plot.
#'
#' @details
#' The function creates a genome-wide Manhattan plot with:
#' - Chromosomes displayed along the x-axis with cumulative positions
#' - Alternating colors for adjacent chromosomes (customizable via `colors`)
#' - -log10(p-values) on the y-axis
#' - Optional significance threshold line
#' - Optional SNP highlighting
#'
#' For multiple tracks (via named list), plots are stacked vertically using facets.
#' For grouped data (via group_var), colors distinguish different groups.
#'
#' @section Regional plots:
#' For LocusZoom-style regional association plots with LD coloring, gene lookup,
#' and stackability with other tracks, use `ez_locusZoom()` instead.
#'
#' @seealso
#' \code{\link{ez_locusZoom}} for regional association plots,
#' \code{\link{geom_manhattan}} for the underlying geom
#'
#' @export
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual scale_x_continuous
#'   scale_y_continuous geom_hline labs theme_minimal theme element_blank facet_wrap element_text
#' @importFrom dplyr mutate arrange group_by ungroup filter
#' @importFrom rlang .data
#'
#' @examples
#' # Basic genome-wide Manhattan plot
#' df <- data.frame(
#'   CHR = rep(1:3, each = 20),
#'   BP = rep(1:20, 3) * 1000,
#'   P = runif(60, 0.0001, 1),
#'   SNP = paste0("rs", 1:60)
#' )
#' ez_manhattan(df)
#'
#' # With alternating chromosome colors
#' ez_manhattan(df, colors = c("navy", "orange"))
#'
#' # With significance threshold
#' ez_manhattan(df, threshold_p = 5e-8)
ez_manhattan <- function(
    input,
    chr = NULL,
    bp = NULL,
    p = NULL,
    snp = NULL,
    track_labels = NULL,
    group_var = NULL,
    logp = TRUE,
    size = 0.5,
    colors = c("grey", "skyblue"),
    highlight_snps = NULL,
    highlight_color = "purple",
    threshold_p = NULL,
    threshold_color = "red",
    threshold_linetype = 2,
    facet_label_position = c("top", "left"),
    # Deprecated parameters (kept for backward compatibility)
    region = NULL,
    gene = NULL,
    ...) {
  # Handle deprecated parameters
  if (!is.null(region) || !is.null(gene)) {
    warning(
      "The 'region' and 'gene' parameters are deprecated in ez_manhattan().\n",
      "For regional association plots (LocusZoom-style), please use ez_locusZoom() instead.\n",
      "Redirecting to ez_locusZoom()...",
      call. = FALSE
    )
    # Extract additional args that ez_locusZoom accepts
    dots <- list(...)

    return(ez_locusZoom(
      input = input,
      region = region,
      gene = gene,
      gene_db = dots$gene_db,
      org_db = dots$org_db,
      extend = dots$extend %||% 0.1,
      extend_type = dots$extend_type %||% "proportion",
      chr = chr,
      bp = bp,
      p = p,
      snp = snp,
      logp = logp,
      size = size,
      color = dots$color %||% "grey50",
      lead_snp = dots$lead_snp,
      r2 = dots$r2,
      colors = if (identical(colors, c("grey", "skyblue"))) NULL else colors,
      highlight_snps = highlight_snps,
      highlight_color = highlight_color,
      threshold_p = threshold_p,
      threshold_color = threshold_color,
      threshold_linetype = threshold_linetype,
      y_axis_style = dots$y_axis_style %||% "none",
      y_axis_label = dots$y_axis_label %||% expression(paste("-log"[10], "(P)"))
    ))
  }

  # Validate inputs
  facet_label_position <- match.arg(facet_label_position)

  # Default color palette function for tracks/groups
  get_default_colors <- function(n) {
    ez_default_palette(n)
  }

  # Process input using helper function
  if (is.data.frame(input)) {
    plotDt <- input
  } else {
    plotDt <- process_manhattan_input(
      input,
      chr = chr,
      bp = bp,
      p = p,
      snp = snp,
      track_labels = track_labels
    )
  }

  # Auto-detect column names if not provided
  if (is.null(chr)) {
    chr <- detect_column(
      plotDt,
      c("CHR", "chr", "seqnames", "chrom", "chromosome"),
      "chromosome"
    )
  }
  if (is.null(bp)) {
    bp <- detect_column(
      plotDt,
      c("BP", "bp", "start", "pos", "position", "POS"),
      "position"
    )
  }
  if (is.null(p)) {
    p <- detect_column(
      plotDt,
      c("P", "p", "pvalue", "p.value", "pval", "P.value"),
      "p-value"
    )
  }
  if (is.null(snp)) {
    snp <- detect_column(
      plotDt,
      c("SNP", "snp", "rsid", "id", "variant_id", "marker"),
      "SNP",
      required = FALSE
    )
  }

  # Check for single-chromosome data (should use ez_locusZoom instead)
  n_chr <- length(unique(plotDt[[chr]]))
  if (n_chr == 1) {
    stop(
      "ez_manhattan() requires data from multiple chromosomes.\n",
      "Your data contains only one chromosome. For single-region plots, ",
      "please use ez_locusZoom() instead:\n",
      "  ez_locusZoom(data, region = \"chr", unique(plotDt[[chr]])[1], ":start-end\")",
      call. = FALSE
    )
  }

  # Determine plotting strategy
  has_track <- "track" %in% colnames(plotDt)
  has_group <- !is.null(group_var) && group_var %in% colnames(plotDt)

  # Create base plot based on grouping/tracking
  if (has_group || has_track) {
    # Determine which variable to use for coloring
    # Priority: group_var if present, otherwise track
    if (has_group) {
      color_var <- group_var
      color_values <- unique(plotDt[[group_var]])
    } else {
      # Only track
      color_var <- "track"
      color_values <- unique(plotDt$track)
    }

    # Get colors for groups/tracks
    n_colors <- length(color_values)
    if (length(colors) >= n_colors) {
      plot_colors <- colors[1:n_colors]
    } else {
      plot_colors <- get_default_colors(n_colors)
    }
    names(plot_colors) <- color_values

    # For each group/track, create a Manhattan plot layer
    plot_obj <- ggplot2::ggplot(plotDt)

    # Add Manhattan layers for each color group
    for (i in seq_along(color_values)) {
      val <- color_values[i]
      subset_data <- plotDt[plotDt[[color_var]] == val, ]

      plot_obj <- plot_obj +
        geom_manhattan(
          data = subset_data,
          mode = "genome_wide",
          chr = chr,
          bp = bp,
          p = p,
          snp = snp,
          logp = logp,
          size = size,
          color = plot_colors[val],
          colors = rep(plot_colors[val], 2),
          highlight_snps = highlight_snps,
          highlight_color = highlight_color,
          threshold_p = threshold_p,
          threshold_color = threshold_color,
          threshold_linetype = threshold_linetype,
          color_by = "none",
          ...
        )
    }

    # Add faceting if multiple tracks
    if (has_track) {
      if (facet_label_position == "left") {
        plot_obj <- plot_obj +
          ggplot2::facet_wrap(
            ~track,
            ncol = 1,
            scales = "free_y",
            strip.position = "left"
          )
      } else {
        plot_obj <- plot_obj +
          ggplot2::facet_wrap(~track, ncol = 1, scales = "free_y")
      }
    }

    plot_obj <- plot_obj + ez_theme()

    # Apply facet label theming last (after ez_theme, so it doesn't get overwritten)
    if (has_track && facet_label_position == "left") {
      plot_obj <- plot_obj +
        ggplot2::theme(
          strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 1),
          strip.placement = "outside"
        )
    }
  } else {
    # Standard single-track genome-wide Manhattan plot
    plot_obj <- ggplot2::ggplot(plotDt) +
      geom_manhattan(
        data = plotDt,
        mode = "genome_wide",
        chr = chr,
        bp = bp,
        p = p,
        snp = snp,
        logp = logp,
        size = size,
        colors = colors,
        color_by = chr, # Color by chromosome for alternating colors
        highlight_snps = highlight_snps,
        highlight_color = highlight_color,
        threshold_p = threshold_p,
        threshold_color = threshold_color,
        threshold_linetype = threshold_linetype,
        ...
      )

    plot_obj <- plot_obj + ez_theme()
  }

  return(plot_obj)
}
