# ezGenomeTracks - Manhattan functions (split from ez_wrappers.R)
#' Easy Manhattan plot visualization
#'
#' This function creates a Manhattan plot for GWAS results.
#' It is a wrapper around geom_manhattan that provides a flexible interface
#' with support for grouping and multiple tracks.
#'
#' @description
#' This function creates a Manhattan plot from GWAS (Genome-Wide Association Study)
#' data, which is a standard way to visualize p-values across the genome.
#' Supports both genome-wide and regional (LocusZoom-style) modes with automatic
#' detection based on data content or explicit region specification.
#'
#' @param input A data frame or named list of data frames containing GWAS results
#'   with columns for chromosome, position, p-values, and optionally SNP names.
#'   Supports both GWAS-style (CHR, BP, P) and GRanges-style (seqnames, start, pvalue)
#'   column naming conventions.
#' @param region Optional genomic region string (e.g., "chr1:1000000-2000000") to
#'   force regional mode. When provided, data is filtered to this region and the
#'   plot uses coordinate-based x-axis consistent with `ez_coverage` and `ez_gene`.
#'   Either `region` or `gene` can be used to specify regional mode.
#' @param gene Optional gene name/symbol to look up (e.g., "PTPRC", "TP53").
#'   When provided, the region is automatically determined from the gene coordinates
#'   in `gene_db`. Forces regional mode.
#' @param gene_db TxDb object for gene coordinate lookup when using `gene` parameter.
#' @param org_db Optional OrgDb object for gene symbol mapping. If NULL (default),
#'   auto-detects available OrgDb packages.
#' @param extend Numeric. Amount to extend the region beyond the gene body when
#'   using `gene` parameter. Default: 0.1 (10% of gene length on each side).
#' @param extend_type How to interpret `extend`: "proportion" (relative to gene
#'   length) or "bp" (absolute base pairs). Default: "proportion".
#' @param chr Character string specifying the column name for chromosome numbers.
#'   Default: "CHR". Also accepts "seqnames", "chrom", etc.
#' @param bp Character string specifying the column name for base pair positions.
#'   Default: "BP". Also accepts "start", "pos", "position", etc.
#' @param p Character string specifying the column name for p-values.
#'   Default: "P". Also accepts "pvalue", "p.value", etc.
#' @param snp Character string specifying the column name for SNP identifiers.
#'   Default: "SNP". Also accepts "rsid", "variant_id", etc.
#' @param track_labels Optional vector of track labels (used for unnamed list input).
#'   Default: NULL.
#' @param group_var Column name for grouping data within a single data frame.
#'   Default: NULL.
#' @param logp Logical indicating whether to plot -log10(p-values).
#'   Default: TRUE.
#' @param size Numeric value for point size in the plot.
#'   Default: 0.5.
#' @param color Default point color for regional mode when color_by is not "r2".
#'   Default: "grey50".
#' @param lead_snp Character string or vector of SNP IDs to highlight as the lead variant(s).
#'   Default: NULL.
#' @param r2 Numeric vector of r² values for coloring points by linkage
#'   disequilibrium (LD) with lead variant. Must be same length as number of
#'   rows in data. Default: NULL.
#' @param colors Vector of colors for coloring points. Usage depends on `color_by`:
#'   - For discrete columns: colors are recycled/mapped to factor levels
#'   - For continuous columns: colors define a gradient (default: viridis-like palette)
#'   - For multi-track or grouped plots: colors for each track/group
#'   Default: NULL (appropriate defaults are chosen automatically).
#' @param highlight_snps Character vector of SNP IDs to highlight.
#'   Default: NULL.
#' @param highlight_color Color for highlighting significant or lead SNPs.
#'   Default: "purple".
#' @param threshold_p Numeric p-value threshold for drawing a significance line.
#'   If NULL, no line is drawn. Default: NULL.
#' @param threshold_color Color for the significance threshold line.
#'   Default: "red".
#' @param threshold_linetype Linetype for the significance threshold line.
#'   Default: 2 (dashed).
#' @param color_by How points should be colored. Can be:
#'   - A column name in the data (e.g., "CHR", "gene", "maf"): Colors by that column.
#'     Use `colors` to specify a custom palette. For chromosome coloring, use
#'     `color_by = "CHR"` (or your chr column name) with `colors = c("grey", "skyblue")`.
#'   - "r2": LD-based gradient coloring (requires `r2` parameter)
#'   - "none": Single color specified by `color` parameter
#'   - "auto" (default): Uses "r2" if `r2` is provided, otherwise "none"
#'   Note: In grouped/multi-track plots, color_by is handled differently.
#' @param y_axis_style Y-axis style: "none", "simple", or "full" (default: "none").
#'   Only applies in regional mode.
#' @param y_axis_label Label for the y-axis. Default: `expression(paste("-log"[10], "(P)"))`.
#' @param facet_label_position Position of facet labels: "top" or "left" (default: "top")
#' @param ... Additional arguments passed to `geom_manhattan()`.
#'
#' @return A ggplot2 object containing the Manhattan plot.
#'
#' @details
#' The function creates a Manhattan plot with chromosomes on the x-axis and
#' -log10(p-values) on the y-axis. The plot mode is automatically determined:
#'
#' - **Regional mode**: When `region` or `gene` is provided OR when data contains only one
#'   chromosome, the plot uses genomic coordinate formatting consistent with
#'   `ez_coverage` and `ez_gene`, making it suitable for stacking with other
#'   tracks via `vstack_plot()`. This is ideal for LocusZoom-style regional
#'   association plots.
#'
#' - **Genome-wide mode**: When data contains multiple chromosomes and no region
#'   is specified, chromosomes are displayed with alternating colors and
#'   cumulative positions.
#'
#' For LD-based coloring (LocusZoom style), provide `r2` values and set
#' `color_by = "r2"`.
#'
#' For multiple tracks (via named list), plots are stacked vertically using facets.
#' For grouped data (via group_var), colors distinguish different groups within tracks.
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
ez_manhattan <- function(
  input,
  region = NULL,
  gene = NULL,
  gene_db = NULL,
  org_db = NULL,
  extend = 0.1,
  extend_type = c("proportion", "bp"),
  chr = NULL,
  bp = NULL,
  p = NULL,
  snp = NULL,
  track_labels = NULL,
  group_var = NULL,
  logp = TRUE,
  size = 0.5,
  color = "grey50",
  lead_snp = NULL,
  r2 = NULL,
  colors = NULL,
  highlight_snps = NULL,
  highlight_color = "purple",
  threshold_p = NULL,
  threshold_color = "red",
  threshold_linetype = 2,
  color_by = "auto",
  y_axis_style = c("none", "simple", "full"),
  y_axis_label = expression(paste("-log"[10], "(P)")),
  facet_label_position = c("top", "left"),
  ...
) {
  # Resolve region from gene name if provided
  if (!is.null(gene)) {
    if (is.null(gene_db)) {
      stop("When using 'gene' parameter, 'gene_db' (TxDb object) must be provided.")
    }
    region <- gene_to_region(
      gene_name = gene,
      txdb = gene_db,
      org_db = org_db,
      extend = extend,
      extend_type = match.arg(extend_type)
    )
  }

  # Validate inputs
  y_axis_style <- match.arg(y_axis_style)
  facet_label_position <- match.arg(facet_label_position)

  # Default color palette function
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

  # Filter data to region if provided
  if (!is.null(region)) {
    region_gr <- parse_region(region)
    region_chr <- as.character(GenomicRanges::seqnames(region_gr))
    region_start <- GenomicRanges::start(region_gr)
    region_end <- GenomicRanges::end(region_gr)

    # Normalize chromosome names for comparison (handle chr prefix)
    region_chr_normalized <- gsub("^chr", "", region_chr, ignore.case = TRUE)

    plotDt <- plotDt |>
      dplyr::filter(
        gsub("^chr", "", as.character(.data[[chr]]), ignore.case = TRUE) ==
          region_chr_normalized,
        .data[[bp]] >= region_start,
        .data[[bp]] <= region_end
      )

    if (nrow(plotDt) == 0) {
      warning("No data points found in the specified region.")
    }
  }

  # Determine plotting strategy
  has_track <- "track" %in% colnames(plotDt)
  has_group <- !is.null(group_var) && group_var %in% colnames(plotDt)

  # Determine if regional mode (for theme application)
  n_chr <- length(unique(plotDt[[chr]]))
  is_regional <- !is.null(region) || n_chr == 1

  # Create base plot based on grouping/tracking
  if (has_group || has_track) {
    # When we have grouping or multiple tracks, we need to modify the color_by approach
    # Note: r2 coloring and group/track coloring are mutually exclusive
    if (!is.null(r2)) {
      warning(
        "r2 coloring is not supported with group_var or multiple tracks. Ignoring r2 parameter."
      )
      r2 <- NULL
    }

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
    if (!is.null(colors) && length(colors) >= n_colors) {
      # User provided custom colors
      plot_colors <- colors[1:n_colors]
    } else {
      # Use default palette
      plot_colors <- get_default_colors(n_colors)
    }
    names(plot_colors) <- color_values

    # For each group/track, create a Manhattan plot layer
    # We'll use geom_manhattan but need to handle coloring differently
    plot_obj <- ggplot2::ggplot(plotDt)

    # Add Manhattan layers for each color group
    for (i in seq_along(color_values)) {
      val <- color_values[i]
      subset_data <- plotDt[plotDt[[color_var]] == val, ]

      plot_obj <- plot_obj +
        geom_manhattan(
          data = subset_data,
          region = region,
          chr = chr,
          bp = bp,
          p = p,
          snp = snp,
          logp = logp,
          size = size,
          color = plot_colors[val],
          lead_snp = lead_snp,
          r2 = NULL, # r2 is disabled for grouped plots
          colors = rep(plot_colors[val], 2), # Use single color for this group
          highlight_snps = highlight_snps,
          highlight_color = highlight_color,
          threshold_p = threshold_p,
          threshold_color = threshold_color,
          threshold_linetype = threshold_linetype,
          color_by = if (color_by == "auto") "none" else color_by,
          y_axis_label = y_axis_label,
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
          ) +
          ggplot2::theme(
            strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 1),
            strip.placement = "outside"
          )
      } else {
        plot_obj <- plot_obj +
          ggplot2::facet_wrap(~track, ncol = 1, scales = "free_y")
      }
    }

    # Apply appropriate theme based on mode
    if (is_regional) {
      plot_obj <- plot_obj + ez_manhattan_theme(y_axis_style = y_axis_style)
    } else {
      plot_obj <- plot_obj + ez_theme()
    }
  } else {
    # Standard single-track Manhattan plot without grouping
    plot_obj <- ggplot2::ggplot(plotDt) +
      geom_manhattan(
        data = plotDt,
        region = region,
        chr = chr,
        bp = bp,
        p = p,
        snp = snp,
        logp = logp,
        size = size,
        color = color,
        lead_snp = lead_snp,
        r2 = r2,
        colors = colors,
        highlight_snps = highlight_snps,
        highlight_color = highlight_color,
        threshold_p = threshold_p,
        threshold_color = threshold_color,
        threshold_linetype = threshold_linetype,
        color_by = color_by,
        y_axis_label = y_axis_label,
        ...
      )

    # Apply appropriate theme based on mode
    if (is_regional) {
      plot_obj <- plot_obj + ez_manhattan_theme(y_axis_style = y_axis_style)
    } else {
      plot_obj <- plot_obj + ez_theme()
    }
  }

  return(plot_obj)
}
