#' Easy signal track visualization
#'
#' This function creates a signal track visualization from various input types.
#' It provides a flexible interface with support for grouping and multiple tracks.
#'
#' @param input A data frame, character vector of file paths, or named list of data sources
#' @param region Genomic region to display (e.g., "chr1:1000000-2000000")
#' @param track_labels Optional vector of track labels (used for character vector input)
#' @param type Type of signal visualization: "line", "area", or "heatmap" (default: "area")
#' @param color Line color (default: "steelblue"). Can be a vector for multiple colors.
#' @param fill Fill color for area plots (default: "steelblue"). Can be a vector for multiple colors.
#' @param group_var Column name for grouping data within a single data frame (default: NULL)
#' @param color_by Whether colors distinguish "group" or "track" (default: "group")
#' @param colors Color palette for groups/tracks. If NULL, uses default palette (default: NULL)
#' @param y_axis_style Y-axis style: "none", "simple", or "full" (default: "none")
#' @param y_range Y-axis range limits (default: NULL)
#' @param alpha Transparency (default: 0.5)
#' @param bin_width Width of bins in base pairs (default: NULL)
#' @param ... Additional arguments passed to geom_coverage
#' @return A ggplot2 object
#' @export
#' @importFrom ggplot2 ggplot aes scale_y_continuous coord_cartesian labs facet_wrap scale_color_manual scale_fill_manual
#' @importFrom dplyr filter mutate bind_rows
#' @examples
#' \dontrun{
#' # Single data frame with grouping
#' df <- data.frame(
#'   seqnames = "chr1", start = 1:100, end = 1:100,
#'   score = rnorm(100), sample = rep(c("A", "B"), 50)
#' )
#' ez_signal(df, "chr1:1-100", group_var = "sample")
#'
#' # Character vector of files
#' ez_signal(c("sample1.bw", "sample2.bw"), "chr1:1-100")
#'
#' # Named list with stacking
#' ez_signal(
#'   list(
#'     "ATAC-seq" = "atac.bw",
#'     "H3K27ac" = "h3k27ac.bw",
#'     "H3K4me3" = c("rep1.bw", "rep2.bw")
#'   ),
#'   "chr1:1-100"
#' )
#' }
ez_signal <- function(input, region, track_labels = NULL,
                      type = c("area", "line", "heatmap"),
                      color = "steelblue", fill = "steelblue",
                      group_var = NULL, color_by = c("group", "track"),
                      colors = NULL, y_axis_style = c("none", "simple", "full"),
                      y_range = NULL, alpha = 0.5, bin_width = NULL, ...) {
  # Validate inputs
  type <- match.arg(type)
  y_axis_style <- match.arg(y_axis_style)
  color_by <- match.arg(color_by)

  stopifnot(
    "alpha must be between 0 and 1" = alpha >= 0 && alpha <= 1,
    "region must be provided" = !missing(region),
    "bin_width must be positive integer" = is.null(bin_width) || (bin_width > 0 && is.integer(bin_width))
  )

  chr <- stringr::str_remove(stringr::str_split(region, ":")[[1]][1], "chr")

  # Default color palette function
  get_default_colors <- function(n) {
    if (n <= 9) {
      return(c("#1f4e79", "#d35400", "#27ae60", "#8e44ad", "#f1c40f", "#16a085", "#e74c3c", "#8b4513", "#5d6d7e")[1:n])
    } else {
      return(rainbow(n))
    }
  }

  # Process input using helper function
  if (is.data.frame(input)) {
    plotDt <- input
  } else {
    plotDt <- process_signal_input(input, region, track_labels)
  }

  # Determine plotting strategy
  has_track <- "track" %in% colnames(plotDt)
  has_group <- !is.null(group_var) && group_var %in% colnames(plotDt)

  # Create base plot
  if (has_group) {
    # Data has grouping - use color/fill mapping
    if (has_track) {
      # Multiple tracks with grouping
      if (color_by == "group") {
        aes_mapping <- ggplot2::aes(x = start, y = score, color = .data[[group_var]], fill = .data[[group_var]])
        color_values <- unique(plotDt[[group_var]])
        legend_name <- group_var
      } else {
        aes_mapping <- ggplot2::aes(x = start, y = score, color = track, fill = track)
        color_values <- unique(plotDt$track)
        legend_name <- "Track"
      }
    } else {
      # Single track with grouping
      aes_mapping <- ggplot2::aes(x = start, y = score, color = .data[[group_var]], fill = .data[[group_var]])
      color_values <- unique(plotDt[[group_var]])
      legend_name <- group_var
    }

    p <- ggplot2::ggplot(plotDt, aes_mapping) +
      geom_coverage(type = type, alpha = alpha, ...)

    # Apply color scales
    n_colors <- length(color_values)
    if (is.null(colors)) {
      plot_colors <- get_default_colors(n_colors)
    } else {
      plot_colors <- colors[1:n_colors]
    }
    names(plot_colors) <- color_values

    p <- p +
      ggplot2::scale_color_manual(values = plot_colors, name = legend_name) +
      ggplot2::scale_fill_manual(values = plot_colors, name = legend_name)
  } else {
    # No grouping - use single color or track-based colors
    if (has_track) {
      # Multiple tracks without grouping
      aes_mapping <- ggplot2::aes(x = start, y = score, color = track, fill = track)
      color_values <- unique(plotDt$track)
      legend_name <- "Track"

      p <- ggplot2::ggplot(plotDt, aes_mapping) +
        geom_coverage(type = type, alpha = alpha, ...)

      # Apply color scales
      n_colors <- length(color_values)
      if (is.null(colors)) {
        plot_colors <- get_default_colors(n_colors)
      } else {
        plot_colors <- colors[1:n_colors]
      }
      names(plot_colors) <- color_values

      p <- p +
        ggplot2::scale_color_manual(values = plot_colors, name = legend_name) +
        ggplot2::scale_fill_manual(values = plot_colors, name = legend_name)
    } else {
      # Single track without grouping
      p <- ggplot2::ggplot(plotDt, ggplot2::aes(x = start, y = score)) +
        geom_coverage(type = type, color = color, fill = fill, alpha = alpha, ...)
    }
  }

  # Add faceting if multiple tracks
  if (has_track) {
    p <- p + ggplot2::facet_wrap(~track, ncol = 1, scales = "free_y")
  }

  # Apply binning if requested
  if (!is.null(bin_width)) {
    p <- p + stat_bin_signal(binwidth = bin_width)
  }

  # Apply the appropriate theme and scale
  p <- p +
    ez_signal_theme(y_axis_style = y_axis_style) +
    scale_x_genome_region(region) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::coord_cartesian(ylim = y_range) +
    ggplot2::labs(x = paste0("Chr", chr))

  return(p)
}

#' Easy peak track visualization
#'
#' This function creates a peak track visualization from a BED file or data frame.
#' It is a wrapper around geom_feature that provides a simpler interface.
#'
#' @param input A BED file path or data frame with peak data
#' @param region Genomic region to display (e.g., "chr1:1000000-2000000")
#' @param color Border color of the peaks (default: "black")
#' @param fill Fill color of the peaks (default: "gray70")
#' @param alpha Transparency (default: 0.7)
#' @param height Height of the peaks (default: 0.8)
#' @param use_score Use the score column for fill color (default: FALSE)
#' @param ... Additional arguments passed to geom_feature
#' @return A ggplot2 object
#' @export
#' @importFrom ggplot2 ggplot aes scale_fill_gradient
#' @examples
#' \dontrun{
#' track <- ez_feature("peaks.bed", "chr1:1000000-2000000", use_score = TRUE)
#' }
ez_feature <- function(input, region, color = "black", fill = "gray70",
                       alpha = 0.7, height = 0.8, use_score = FALSE, ...) {
  # Check if data is a file path or data frame
  if (is.character(input) && length(input) == 1) {
    # It's a file path, use peak_track
    return(peak_track(input, region,
      color = color, fill = fill,
      alpha = alpha, height = height, use_score = use_score, ...
    ))
  } else if (is.data.frame(input)) {
    # It's a data frame, create the plot directly
    if (use_score && "score" %in% colnames(input)) {
      p <- ggplot2::ggplot(input) +
        geom_feature(ggplot2::aes(xmin = start, xmax = end, fill = score),
          color = color, alpha = alpha, height = height, ...
        ) +
        ggplot2::scale_fill_gradient(low = "white", high = fill)
    } else {
      p <- ggplot2::ggplot(input) +
        geom_feature(ggplot2::aes(xmin = start, xmax = end),
          color = color, fill = fill, alpha = alpha, height = height, ...
        )
    }

    # Apply the appropriate theme and scale
    p <- p +
      scale_x_genome_region(region) +
      ggplot2::ylim(0, 1) + # Fixed y-axis for features
      ggplot2::theme_void()

    return(p)
  } else {
    stop("Data must be a file path or data frame")
  }
}

#' Easy Manhattan plot visualization
#'
#' This function creates a Manhattan plot for GWAS results.
#' It is a wrapper around geom_manhattan that provides a simpler interface.
#'
#' @param data A data frame with GWAS results. Must include columns for chromosome, base pair position, and p-value. Optional column for SNP identifier and R-squared values.
#' @param chr Name of the chromosome column in `data` (default: "CHR").
#' @param bp Name of the base pair position column in `data` (default: "BP").
#' @param p Name of the p-value column in `data` (default: "P").
#' @param snp Name of the SNP identifier column in `data` (default: "SNP").
#' @param logp Logical. If `TRUE` (default), -log10(p-value) is used for the y-axis. If `FALSE`, raw p-value is used.
#' @param size Point size (default: 0.5).
#' @param lead.snp Vector of leading SNP identifiers to highlight.
#' @param r2 Vector of R-squared values for linkage disequilibrium (LD) coloring.
#' @param colors Vector of colors to use for alternating chromosome colors (default: c("grey", "skyblue")).
#' @param highlight_snps Data frame of SNPs to highlight, with columns `CHR`, `BP`, and `P`.
#' @param highlight_color Color for highlighted SNPs (default: "purple").
#' @param threshold_p A numeric value for the p-value threshold to draw a horizontal line (e.g., 5e-8).
#' @param threshold_color Color for the threshold line (default: "red").
#' @param threshold_linetype Linetype for the threshold line (default: 2).
#' @param colorBy Character string indicating how points should be colored.
#' Options are "chr" (default, alternating chromosome colors in colors) or "r2" (based on R-squared values).
#' @param y_axis_label Label for the y-axis (default: `expression(paste("-log"[10], "(P)"))`).
#' @param ... Additional arguments passed to `geom_manhattan()`.
#' @return A `ggplot2` object.
#' @export
ez_manhattan <- function(
  data,
  chr = "CHR", bp = "BP", p = "P", snp = "SNP", logp = TRUE, size = 0.5,
  lead.snp = NULL, r2 = NULL, colors = c("grey", "skyblue"),
  highlight_snps = NULL, highlight_color = "purple",
  threshold_p = NULL, threshold_color = "red", threshold_linetype = 2,
  colorBy = "chr", y_axis_label = expression(paste("-log"[10], "(P)")), ...
) {
  if (!is.data.frame(data)) {
    stop("Input 'data' must be a data.frame.")
  }

  # Create the plot with minimal aesthetics
  p <- ggplot2::ggplot(data) +
    geom_manhattan(
      data = data,
      chr = chr, bp = bp, p = p, snp = snp, logp = logp,
      size = size,
      lead.snp = lead.snp,
      r2 = r2,
      colors = colors,
      highlight_snps = highlight_snps,
      highlight_color = highlight_color,
      threshold_p = threshold_p,
      threshold_color = threshold_color,
      threshold_linetype = threshold_linetype,
      colorBy = colorBy,
      y_axis_label = y_axis_label,
      ...
    ) +
    ez_theme()

  return(p)
}

#' Easy gene track visualization
#'
#' This function creates a gene track visualization from a GTF/GFF file, TxDb object, or data frame.
#' It is a wrapper around geom_gene that provides a simpler interface.
#'
#' @param data A GTF/GFF file path, TxDb object, or data frame with gene data
#' @param region Genomic region to display (e.g., "chr1:1000000-2000000")
#' @param exon_height Height of exons (default: 0.75)
#' @param intron_width Width of introns (default: 0.4)
#' @param exon_color Color of exon borders (default: "black")
#' @param exon_fill Fill color of exons (default: "gray50")
#' @param intron_color Color of introns (default: "gray50")
#' @param gene_id Column name for gene ID (default: "gene_id")
#' @param gene_name Column name for gene name (default: "gene_name")
#' @param ... Additional arguments passed to geom_gene
#' @return A ggplot2 object
#' @export
#' @importFrom ggplot2 ggplot aes
#' @importFrom methods is
#' @examples
#' \dontrun{
#' # Using a GTF file
#' track1 <- ez_gene("genes.gtf", "chr1:1000000-2000000")
#'
#' # Using a TxDb object
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' track2 <- ez_gene(txdb, "chr1:1000000-2000000")
#' }
ez_gene <- function(data, region, exon_height = 0.75, intron_width = 0.4,
                    exon_color = "black", exon_fill = "gray50", intron_color = "gray50",
                    gene_id = "gene_id", gene_name = "gene_name", ...) {
  # Parse the region
  region_gr <- parse_region(region)

  # Process data based on input type
  if (is.character(data) && length(data) == 1) {
    # GTF/GFF file path
    gene_gr <- rtracklayer::import(data, which = region_gr)
    gene_data <- process_gene_data(gene_gr, gene_id = gene_id, gene_name = gene_name)
  } else if (methods::is(data, "TxDb")) {
    # TxDb object
    if (!requireNamespace("GenomicFeatures", quietly = TRUE)) {
      stop("Package 'GenomicFeatures' is required for TxDb support. Install it with: BiocManager::install('GenomicFeatures')")
    }
    gene_data <- extract_txdb_data(data, region_gr)
  } else if (is.data.frame(data)) {
    # Data frame - use as-is
    gene_data <- data
  } else {
    stop("Data must be a file path, TxDb object, or data frame")
  }

  # Create the plot
  p <- ggplot2::ggplot(gene_data) +
    geom_gene(
      ggplot2::aes(xstart = xstart, xend = xend, y = strand, type = type),
      exon_height = exon_height, intron_width = intron_width,
      exon_color = exon_color, exon_fill = exon_fill,
      intron_color = intron_color, ...
    )

  # Apply theme and scale
  p <- p +
    ggplot2::scale_y_discrete(
      expand = c(0.1, 0.1),
      drop = FALSE # Keep all levels even if not present
    ) +
    scale_x_genome_region(region) +
    ez_gene_theme()

  return(p)
}

#' Easy interaction track visualization
#'
#' This function creates an interaction track visualization from a BEDPE file or data frame.
#' It is a wrapper around geom_arc that provides a simpler interface.
#'
#' @param data A BEDPE file path or data frame with interaction data
#' @param region Genomic region to display (e.g., "chr1:1000000-2000000")
#' @param curvature Amount of curvature (default: 0.5)
#' @param color Color of the arcs (default: "gray50")
#' @param size Size of the arcs (default: 0.5)
#' @param alpha Transparency (default: 0.7)
#' @param use_score Use the score column for color (default: FALSE)
#' @param ... Additional arguments passed to geom_arc
#' @return A ggplot2 object
#' @export
#' @importFrom ggplot2 ggplot aes scale_color_gradient
#' @examples
#' \dontrun{
#' track <- ez_arc("interactions.bedpe", "chr1:1000000-2000000", use_score = TRUE)
#' }
ez_arc <- function(data, region, curvature = 0.5, color = "gray50",
                   size = 0.5, alpha = 0.7, use_score = FALSE, ...) {
  # Check if data is a file path or data frame
  if (is.character(data) && length(data) == 1) {
    # It's a file path, use interaction_track
    return(interaction_track(data, region,
      curvature = curvature, color = color,
      size = size, alpha = alpha, use_score = use_score, ...
    ))
  } else if (is.data.frame(data)) {
    # It's a data frame, create the plot directly
    if (use_score && "score" %in% colnames(data)) {
      p <- ggplot2::ggplot(data) +
        geom_arc(ggplot2::aes(x = start1, y = 0, xend = start2, yend = 0, color = score),
          curvature = curvature, size = size, alpha = alpha, ...
        ) +
        ggplot2::scale_color_gradient(low = "blue", high = "red")
    } else {
      p <- ggplot2::ggplot(data) +
        geom_arc(ggplot2::aes(x = start1, y = 0, xend = start2, yend = 0),
          curvature = curvature, color = color, size = size, alpha = alpha, ...
        )
    }

    # Apply the appropriate theme and scale
    p <- p + ez_feature_theme() +
      scale_x_genome_region(region) +
      ggplot2::ylim(-0.5, 0.5) # Fixed y-axis for arcs

    return(p)
  } else {
    stop("Data must be a file path or data frame")
  }
}

#' Easy Hi-C track visualization
#'
#' This function creates a Hi-C track visualization from a contact matrix file or data frame.
#' It is a wrapper around geom_hic that provides a simpler interface.
#'
#' @param data A contact matrix file path or data frame with Hi-C data
#' @param region Genomic region to display (e.g., "chr1:1000000-2000000")
#' @param resolution Resolution of the Hi-C data in base pairs (default: 10000)
#' @param log_transform Apply log transformation to the values (default: TRUE)
#' @param low Color for low values (default: "white")
#' @param high Color for high values (default: "red")
#' @param ... Additional arguments passed to geom_hic
#' @return A ggplot2 object
#' @export
#' @importFrom ggplot2 ggplot aes coord_fixed
#' @examples
#' \dontrun{
#' track <- ez_hic("contacts.matrix", "chr1:1000000-2000000", resolution = 10000)
#' }
ez_hic <- function(data, region, resolution = 10000, log_transform = TRUE,
                   low = "white", high = "red", ...) {
  # Check if data is a file path or data frame
  if (is.character(data) && length(data) == 1) {
    # It's a file path, use hic_track
    return(hic_track(data, region,
      resolution = resolution, log_transform = log_transform,
      low = low, high = high, ...
    ))
  } else if (is.data.frame(data)) {
    # It's a data frame, create the plot directly
    # Handle different column names for Hi-C data
    if ("bin1" %in% colnames(data) && "bin2" %in% colnames(data)) {
      # Convert bin coordinates to genomic positions
      region_gr <- parse_region(region)
      start_pos <- GenomicRanges::start(region_gr)
      bin_size <- resolution
      data$pos1 <- start_pos + (data$bin1 - 1) * bin_size
      data$pos2 <- start_pos + (data$bin2 - 1) * bin_size
    }

    p <- ggplot2::ggplot(data, ggplot2::aes(x = pos1, y = pos2, fill = count)) +
      geom_hic(low = low, high = high, ...) +
      ggplot2::coord_fixed() # Ensure the plot is square

    # Apply the appropriate theme and scale
    p <- p + ez_theme() +
      scale_x_genome_region(region) +
      scale_x_genome_region(region) # Same scale for y-axis

    return(p)
  } else {
    stop("Data must be a file path or data frame")
  }
}
