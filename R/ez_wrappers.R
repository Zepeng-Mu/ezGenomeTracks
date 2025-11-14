#' Easy coverage track visualization
#'
#' This function creates a coverage track visualization from various input types.
#' It provides a flexible interface with support for grouping and multiple tracks.
#'
#' @param input A data frame, character vector of file paths, or named list of data sources
#' @param region Genomic region to display (e.g., "chr1:1000000-2000000")
#' @param track_labels Optional vector of track labels (used for character vector input)
#' @param type Type of signal visualization: "line", "area", or "heatmap" (default: "area")
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
#' ez_coverage(df, "chr1:1-100", group_var = "sample")
#'
#' # Character vector of files
#' ez_coverage(c("sample1.bw", "sample2.bw"), "chr1:1-100")
#'
#' # Named list with stacking
#' ez_coverage(
#'   list(
#'     "ATAC-seq" = "atac.bw",
#'     "H3K27ac" = "h3k27ac.bw",
#'     "H3K4me3" = c("rep1.bw", "rep2.bw")
#'   ),
#'   "chr1:1-100"
#' )
#' }
ez_coverage <- function(
  input,
  region,
  track_labels = NULL,
  type = c("area", "line", "heatmap"),
  fill = "steelblue",
  group_var = NULL,
  color_by = c("group", "track"),
  colors = NULL,
  y_axis_style = c("none", "simple", "full"),
  y_range = NULL,
  alpha = 0.5,
  bin_width = NULL,
  ...
) {
  # Validate inputs
  type <- match.arg(type)
  y_axis_style <- match.arg(y_axis_style)
  color_by <- match.arg(color_by)

  stopifnot(
    "alpha must be between 0 and 1" = alpha >= 0 && alpha <= 1,
    "region must be provided" = !missing(region),
    "bin_width must be positive integer" = is.null(bin_width) ||
      (bin_width > 0 && is.integer(bin_width))
  )

  chr <- stringr::str_remove(stringr::str_split(region, ":")[[1]][1], "chr")

  # Default color palette function
  get_default_colors <- function(n) {
    if (n <= 9) {
      return(c(
        "#1f4e79",
        "#d35400",
        "#27ae60",
        "#8e44ad",
        "#f1c40f",
        "#16a085",
        "#e74c3c",
        "#8b4513",
        "#5d6d7e"
      )[1:n])
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
        aes_mapping <- ggplot2::aes(
          xmin = start,
          xmax = end,
          ymin = 0,
          ymax = score,
          fill = .data[[group_var]]
        )
        color_values <- unique(plotDt[[group_var]])
        legend_name <- group_var
      } else {
        aes_mapping <- ggplot2::aes(
          xmin = start,
          xmax = end,
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
        xmin = start,
        xmax = end,
        ymin = 0,
        ymax = score,
        fill = .data[[group_var]]
      )
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
      ggplot2::scale_fill_manual(values = plot_colors, name = legend_name)
  } else {
    # No grouping - use single color or track-based colors
    if (has_track) {
      # Multiple tracks without grouping
      aes_mapping <- ggplot2::aes(
        xmin = start,
        xmax = end,
        ymin = 0,
        ymax = score,
        fill = track
      )
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
        ggplot2::scale_fill_manual(values = plot_colors, name = legend_name)
    } else {
      # Single track without grouping
      p <- ggplot2::ggplot(plotDt, ggplot2::aes(xmin = start, xmax = end, ymin = 0, ymax = score)) +
        geom_coverage(
          type = type,
          fill = fill,
          alpha = alpha,
          ...
        )
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

#' @export
#' @rdname ez_coverage
ez_signal <- ez_coverage

#' Easy peak track visualization
#'
#' This function creates a peak track visualization from a BED file or data frame.
#' It is a wrapper around geom_feature that provides a simpler interface.
#'
#' @param input A BED file path or data frame with peak data
#' @param region Genomic region to display (e.g., "chr1:1000000-2000000")
#' Create a feature track from genomic regions
#'
#' @description
#' This function creates a feature track visualization from genomic regions,
#' such as peaks or other genomic annotations. It can read directly from BED files
#' or work with data frames containing genomic coordinates.
#'
#' @param input Either a file path to a BED file or a data frame containing
#'   genomic coordinates with columns for chromosome, start, and end positions.
#' @param region Genomic region to display in the format "chr:start-end".
#'   Example: "chr1:1000000-2000000"
#' @param color Border color of the features. Default: "black"
#' @param fill Fill color of the features. When `use_score = TRUE`, this will be
#'   used as the high value in the color gradient. Default: "gray70"
#' @param alpha Transparency level of the features (0 = transparent, 1 = opaque).
#'   Default: 0.7
#' @param height Height of the feature rectangles (0 to 1). Default: 0.8
#' @param use_score Logical indicating whether to use the 'score' column for
#'   fill color. If TRUE, a gradient from white to the specified fill color will
#'   be used. Default: FALSE
#' @param ... Additional arguments passed to `geom_feature()`
#'
#' @return A ggplot2 object representing the feature track.
#'
#' @details
#' The function automatically handles both file paths and data frames as input.
#' When a file path is provided, it reads the BED file and creates the track.
#' When a data frame is provided, it should contain at least 'chrom', 'start',
#' and 'end' columns. If 'score' column is present and `use_score = TRUE`,
#' features will be colored by their score values.
#'
#' @export
#' @importFrom ggplot2 ggplot aes scale_fill_gradient scale_fill_identity
#'   theme_void ylim
#' @importFrom methods is
#'
#' @examples
#' \dontrun{
#' # From a BED file with score-based coloring
#' track1 <- ez_feature(
#'   "peaks.bed",
#'   "chr1:1000000-2000000",
#'   fill = "blue",
#'   use_score = TRUE
#' )
#'
#' # From a data frame with uniform coloring
#' features <- data.frame(
#'   chrom = c("chr1", "chr1", "chr1"),
#'   start = c(1000, 3000, 5000),
#'   end = c(2000, 4000, 6000),
#'   name = c("peak1", "peak2", "peak3"),
#'   score = c(10, 30, 50)
#' )
#' track2 <- ez_feature(
#'   features,
#'   "chr1:1-10000",
#'   fill = "darkgreen",
#'   alpha = 0.8
#' )
#'
#' # Combine with other tracks using aplot
#' ez_plot(list(
#'   "Features" = track1,
#'   "Genes" = track2
#' ), "chr1:1-10000")
#' }
ez_feature <- function(
  input,
  region,
  color = "black",
  fill = "gray70",
  alpha = 0.7,
  height = 0.8,
  use_score = FALSE,
  ...
) {
  # Check if data is a file path or data frame
  if (is.character(input) && length(input) == 1) {
    # It's a file path, use peak_track
    return(peak_track(
      input,
      region,
      color = color,
      fill = fill,
      alpha = alpha,
      height = height,
      use_score = use_score,
      ...
    ))
  } else if (is.data.frame(input)) {
    # It's a data frame, create the plot directly
    if (use_score && "score" %in% colnames(input)) {
      p <- ggplot2::ggplot(input) +
        geom_feature(
          ggplot2::aes(xmin = start, xmax = end, fill = score),
          color = color,
          alpha = alpha,
          height = height,
          ...
        ) +
        ggplot2::scale_fill_gradient(low = "white", high = fill)
    } else {
      p <- ggplot2::ggplot(input) +
        geom_feature(
          ggplot2::aes(xmin = start, xmax = end),
          color = color,
          fill = fill,
          alpha = alpha,
          height = height,
          ...
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
#' @description
#' This function creates a Manhattan plot from GWAS (Genome-Wide Association Study)
#' data, which is a standard way to visualize p-values across the genome.
#'
#' @param data A data frame containing GWAS results with columns for chromosome,
#'   position, p-values, and optionally SNP names.
#' @param chr Character string specifying the column name for chromosome numbers.
#'   Default: "CHR".
#' @param bp Character string specifying the column name for base pair positions.
#'   Default: "BP".
#' @param p Character string specifying the column name for p-values.
#'   Default: "P".
#' @param snp Character string specifying the column name for SNP identifiers.
#'   Default: "SNP".
#' @param logp Logical indicating whether to plot -log10(p-values).
#'   Default: TRUE.
#' @param size Numeric value for point size in the plot.
#'   Default: 0.5.
#' @param lead.snp Character string of SNP ID to highlight as the lead variant.
#'   Default: NULL.
#' @param r2 Numeric vector of r² values for coloring points by linkage
#'   disequilibrium (LD) with lead variant. Must be same length as number of
#'   rows in data. Default: NULL.
#' @param colors Character vector of length 2 specifying colors for
#'   alternating chromosomes. Default: c("grey", "skyblue").
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
#' @param colorBy Character string specifying the variable to use for coloring points.
#'   Must be a column name in the data. Default: "chr".
#' @param y_axis_label Label for the y-axis. Default: `expression(paste("-log"[10], "(P)"))`.
#' @param ... Additional arguments passed to `geom_manhattan()`.
#'
#' @return A ggplot2 object containing the Manhattan plot.
#'
#' @details
#' The function creates a Manhattan plot with chromosomes on the x-axis and
#' -log10(p-values) on the y-axis. Points are colored by chromosome and can be
#' highlighted based on significance or LD with lead variants.
#'
#' @export
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual scale_x_continuous
#'   scale_y_continuous geom_hline labs theme_minimal theme element_blank
#' @importFrom dplyr mutate arrange group_by ungroup
#' @importFrom rlang .data
#'
#' @examples
#' \dontrun{
#' # Basic Manhattan plot
#' data(gwas_data)  # Example GWAS data
#' ez_manhattan(
#'   gwas_data,
#'   chr = "CHR",
#'   bp = "BP",
#'   p = "P",
#'   snp = "SNP",
#'   colors = c("dodgerblue", "darkblue")
#' )
#'
#' # With highlighted lead SNP and threshold line
#' ez_manhattan(
#'   gwas_data,
#'   chr = "CHR",
#'   bp = "BP",
#'   p = "P",
#'   snp = "SNP",
#'   lead.snp = "rs123456",
#'   highlight_color = "red",
#'   threshold_p = 5e-8,
#'   threshold_color = "red"
#' )
#' }
ez_manhattan <- function(
  data,
  chr = "CHR",
  bp = "BP",
  p = "P",
  snp = "SNP",
  logp = TRUE,
  size = 0.5,
  lead.snp = NULL,
  r2 = NULL,
  colors = c("grey", "skyblue"),
  highlight_snps = NULL,
  highlight_color = "purple",
  threshold_p = NULL,
  threshold_color = "red",
  threshold_linetype = 2,
  colorBy = "chr",
  y_axis_label = expression(paste("-log"[10], "(P)")),
  ...
) {
  if (!is.data.frame(data)) {
    stop("Input 'data' must be a data.frame.")
  }

  # Create the plot with minimal aesthetics
  p <- ggplot2::ggplot(data) +
    geom_manhattan(
      data = data,
      chr = chr,
      bp = bp,
      p = p,
      snp = snp,
      logp = logp,
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
#' @param color Color for both exons and introns (default: "gray50")
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
#' Create a gene track from genomic annotations
#'
#' @description
#' This function creates a gene track visualization from genomic annotations,
#' supporting various input formats including GTF/GFF files, TxDb objects, and
#' data frames. It automatically handles gene structure visualization with
#' exons, introns, and strand information.
#'
#' @param data Input data source, which can be:
#'   - A file path to a GTF/GFF file
#'   - A TxDb object from the GenomicFeatures package
#'   - A data frame with gene annotation data
#' @param region Genomic region to display in the format "chr:start-end".
#'   Example: "chr1:1000000-2000000"
#' @param exon_height Relative height of exons (0 to 1). Default: 0.75
#' @param intron_width Line width for introns. Default: 0.4
#' @param exon_color Border color for exons. Default: "gray50"
#' @param exon_fill Fill color for exons. Default: "gray50"
#' @param intron_color Color for intron lines. Default: "gray50"
#' @param gene_id Column name for gene identifiers. Default: "gene_id"
#' @param gene_name Column name for gene symbols/names. Default: "gene_name"
#' @param label Column name to use for text labels. If NULL (default), no labels
#'   are displayed. Set to a column name (e.g., "gene_name") to show labels.
#' @param label_size Size of text labels. Default: 3
#' @param label_vjust Vertical adjustment for labels. Negative values place labels
#'   above genes. Default: -1.5. Automatically adjusted based on exon_height.
#' @param label_color Color of text labels. Default: "black"
#' @param ... Additional arguments passed to `geom_gene()`
#'
#' @return A ggplot2 object representing the gene track.
#'
#' @details
#' The function automatically processes different input types:
#' - For GTF/GFF files: Uses rtracklayer to import and process the data
#' - For TxDb objects: Extracts gene models using GenomicFeatures
#' - For data frames: Expects columns for chromosome, start, end, strand, and type
#'
#' The visualization includes:
#' - Exons as filled rectangles
#' - Introns as connecting lines
#' - Strand information with arrowheads
#' - Automatic y-axis separation by strand
#'
#' @export
#' @importFrom GenomicFeatures exonsBy transcriptsBy
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
#' @importFrom rtracklayer import
#' @importFrom methods is
#'
#' @examples
#' \dontrun{
#' # From a GTF file
#' track1 <- ez_gene("genes.gtf", "chr1:1000000-2000000")
#'
#' # From a TxDb object
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' track2 <- ez_gene(txdb, "chr1:1000000-2000000")
#'
#' # With custom styling
#' track3 <- ez_gene("genes.gtf", "chr1:1000000-2000000",
#'   exon_fill = "steelblue",
#'   exon_color = "navy",
#'   intron_color = "darkblue",
#'   intron_width = 0.6
#' )
#'
#' # With gene name labels
#' track4 <- ez_gene("genes.gtf", "chr1:1000000-2000000",
#'   label = "gene_name",
#'   label_size = 3
#' )
#' }
ez_gene <- function(
  data,
  region,
  exon_height = 0.75,
  intron_width = 0.4,
  exon_color = "gray50",
  exon_fill = "gray50",
  intron_color = "gray50",
  gene_id = "gene_id",
  gene_name = "gene_name",
  label = NULL,
  label_size = 3,
  label_vjust = -2,
  label_color = "black",
  ...
) {
  # Parse the region
  region_gr <- parse_region(region)

  # Extract region limits for clipping
  region_limits <- c(
    GenomicRanges::start(region_gr),
    GenomicRanges::end(region_gr)
  )

  # Process data based on input type
  if (is.character(data) && length(data) == 1) {
    # GTF/GFF file path
    gene_gr <- rtracklayer::import(data, which = region_gr)
    gene_data <- process_gene_data(
      gene_gr,
      gene_id = gene_id,
      gene_name = gene_name
    )
  } else if (methods::is(data, "TxDb")) {
    # TxDb object
    if (!requireNamespace("GenomicFeatures", quietly = TRUE)) {
      stop(
        "Package 'GenomicFeatures' is required for TxDb support. Install it with: BiocManager::install('GenomicFeatures')"
      )
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
      exon_height = exon_height,
      intron_width = intron_width,
      exon_color = exon_color,
      exon_fill = exon_fill,
      intron_color = intron_color,
      clip_to_region = region_limits,
      ...
    )

  # Add labels if requested
  if (!is.null(label) && label %in% names(gene_data)) {
    # Get unique gene positions for labels (use gene type, not exons)
    label_data <- gene_data[gene_data$type == "gene", ]

    # Remove duplicates based on gene_id to avoid duplicate labels
    if (gene_id %in% names(label_data)) {
      label_data <- label_data[!duplicated(label_data[[gene_id]]), ]
    }

    # Calculate label positions (middle of gene)
    label_data$label_x <- (label_data$xstart + label_data$xend) / 2

    # Use strand directly (it's already a factor that matches the y-axis)
    label_data$label_y <- label_data$strand

    # Calculate effective vjust accounting for exon height
    # vjust < 0 places text above the point; we need more negative values for taller exons
    # Adjust by the exon height to ensure labels clear the top of exons
    effective_vjust <- label_vjust - exon_height

    p <- p +
      ggplot2::geom_text(
        data = label_data,
        ggplot2::aes(x = .data$label_x, y = .data$label_y, label = .data[[label]]),
        size = label_size,
        vjust = effective_vjust,
        color = label_color,
        inherit.aes = FALSE
      )
  }

  # Apply theme and scale
  p <- p +
    ggplot2::scale_y_discrete(
      expand = c(0.1, 0.1),
      drop = FALSE # Keep all levels even if not present
    ) +
    scale_x_genome_region(region) +
    ez_gene_theme() +
    ggplot2::labs(x = NULL)  # Remove x-axis title

  return(p)
}

#' Easy interaction track visualization
#'
#' This function creates an interaction track visualization from a BEDPE file or data frame.
#' It is a wrapper around geom_arc that provides a simpler interface.
#'
#' @param data A BEDPE file path or data frame with interaction data
#' Create an arc track for genomic interactions
#'
#' @description
#' This function creates an arc track visualization for genomic interactions,
#' such as those from Hi-C, ChIA-PET, or other interaction assays. It can
#' read directly from BEDPE files or work with data frames containing
#' interaction data.
#'
#' @param data Either:
#'   - A file path to a BEDPE file containing interaction data
#'   - A data frame with columns for interaction coordinates and scores
#' @param region Genomic region to display in the format "chr:start-end".
#'   Example: "chr1:1000000-2000000"
#' @param curvature Numeric value controlling the arc curvature (0-1).
#'   Higher values create more pronounced curves. Default: 0.5
#' @param color Color of the arcs. This is used when `use_score = FALSE`.
#'   Default: "gray50"
#' @param size Line width of the arcs. Default: 0.5
#' @param alpha Transparency level of the arcs (0 = transparent, 1 = opaque).
#'   Default: 0.7
#' @param use_score Logical indicating whether to use the 'score' column for
#'   arc coloring. If TRUE, a color gradient will be applied based on the
#'   interaction scores. Default: FALSE
#' @param ... Additional arguments passed to `geom_arc()`
#'
#' @return A ggplot2 object representing the arc track.
#'
#' @details
#' The function automatically handles different input types:
#' - For BEDPE files: Uses `interaction_track` to read and process the data
#' - For data frames: Expects columns for interaction coordinates (chr1, start1, end1, chr2, start2, end2)
#'   and optionally a 'score' column if `use_score = TRUE`
#'
#' The visualization includes:
#' - Arcs connecting interaction anchors
#' - Optional score-based coloring
#' - Automatic scaling to fit the specified genomic region
#'
#' @export
#' @importFrom ggplot2 ggplot aes scale_color_gradient
#' @importFrom methods is
#'
#' @examples
#' \dontrun{
#' # From a BEDPE file with score-based coloring
#' track1 <- ez_arc(
#'   "interactions.bedpe",
#'   "chr1:1000000-2000000",
#'   use_score = TRUE,
#'   high = "red",
#'   low = "blue"
#' )
#'
#' # From a data frame with uniform coloring
#' interactions <- data.frame(
#'   chr1 = c("chr1", "chr1", "chr1"),
#'   start1 = c(1000, 3000, 5000),
#'   end1 = c(2000, 4000, 6000),
#'   chr2 = c("chr1", "chr1", "chr1"),
#'   start2 = c(8000, 7000, 9000),
#'   end2 = c(9000, 8000, 10000),
#'   score = c(5, 10, 15)
#' )
#' track2 <- ez_arc(
#'   interactions,
#'   "chr1:1-15000",
#'   color = "darkblue",
#'   size = 1,
#'   alpha = 0.8
#' )
#' }
ez_arc <- function(
  data,
  region,
  curvature = 0.5,
  color = "gray50",
  size = 0.5,
  alpha = 0.7,
  use_score = FALSE,
  ...
) {
  # Check if data is a file path or data frame
  if (is.character(data) && length(data) == 1) {
    # It's a file path, use interaction_track
    return(interaction_track(
      data,
      region,
      curvature = curvature,
      color = color,
      size = size,
      alpha = alpha,
      use_score = use_score,
      ...
    ))
  } else if (is.data.frame(data)) {
    # It's a data frame, create the plot directly
    if (use_score && "score" %in% colnames(data)) {
      p <- ggplot2::ggplot(data) +
        geom_arc(
          ggplot2::aes(
            x = start1,
            y = 0,
            xend = start2,
            yend = 0,
            color = score
          ),
          curvature = curvature,
          size = size,
          alpha = alpha,
          ...
        ) +
        ggplot2::scale_color_gradient(low = "blue", high = "red")
    } else {
      p <- ggplot2::ggplot(data) +
        geom_arc(
          ggplot2::aes(x = start1, y = 0, xend = start2, yend = 0),
          curvature = curvature,
          color = color,
          size = size,
          alpha = alpha,
          ...
        )
    }

    # Apply the appropriate theme and scale
    p <- p +
      ez_feature_theme() +
      scale_x_genome_region(region) +
      ggplot2::ylim(-0.5, 0.5) # Fixed y-axis for arcs

    return(p)
    return(hic_track(
      data,
      region,
      resolution = resolution,
      log_transform = log_transform,
      low = low,
      high = high,
      ...
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
    p <- p +
      ez_theme() +
      scale_x_genome_region(region) +
      scale_x_genome_region(region) # Same scale for y-axis

    return(p)
  } else {
    stop("Data must be a file path or data frame")
  }
}

# Declare global variables used in aes() mappings to avoid R CMD check notes
utils::globalVariables(c(
  "type", "label_x", "label_y", "start1", "start2", "score",
  "resolution", "log_transform", "low", "high", "pos1", "pos2", "count",
  "track", "start", "end"
))
