#' Easy signal track visualization
#'
#' This function creates a signal track visualization from a bigWig file or data frame.
#' It is a wrapper around geom_signal that provides a simpler interface.
#'
#' @param data A bigWig file path or data frame with signal data
#' @param region Genomic region to display (e.g., "chr1:1000000-2000000")
#' @param type Type of signal visualization: "line", "area", or "heatmap" (default: "area")
#' @param color Line color (default: "steelblue")
#' @param fill Fill color for area plots (default: "steelblue")
#' @param alpha Transparency (default: 0.5)
#' @param binwidth Width of bins in base pairs (default: NULL)
#' @param ... Additional arguments passed to geom_signal
#' @return A ggplot2 object
#' @export
#' @importFrom ggplot2 ggplot aes
#' @examples
#' \dontrun{
#' track <- ez_signal("signal.bw", "chr1:1000000-2000000")
#' }
ez_signal <- function(data, region, type = "area", color = "steelblue",
                      fill = "steelblue", yrange = NULL, alpha = 0.5, binwidth = NULL, ...) {
  # Check if data is a file path or data frame
  if (is.character(data) && length(data) == 1) {
    # It's a file path, use signal_track
    return(signal_track(data, region, type = type, color = color,
                        fill = fill, alpha = alpha, binwidth = binwidth, ...))
  } else if (is.data.frame(data)) {
    # It's a data frame, create the plot directly
    p <- ggplot2::ggplot(data, ggplot2::aes(x = start, y = score)) +
      geom_signal(type = type, color = color, fill = fill, alpha = alpha, ...)

    # Apply binning if requested
    if (!is.null(binwidth)) {
      p <- p + stat_bin_signal(binwidth = binwidth)
    }

    # Apply the appropriate theme and scale
    p <- p + ez_signal_theme() +
      scale_x_genome_region(region) +
      scale_y_continuous(expand = c(0, 0), limits = yrange)

    return(p)
  } else {
    stop("Data must be a file path or data frame")
  }
}

#' Easy peak track visualization
#'
#' This function creates a peak track visualization from a BED file or data frame.
#' It is a wrapper around geom_peak that provides a simpler interface.
#'
#' @param data A BED file path or data frame with peak data
#' @param region Genomic region to display (e.g., "chr1:1000000-2000000")
#' @param color Border color of the peaks (default: "black")
#' @param fill Fill color of the peaks (default: "gray70")
#' @param alpha Transparency (default: 0.7)
#' @param height Height of the peaks (default: 0.8)
#' @param use_score Use the score column for fill color (default: FALSE)
#' @param ... Additional arguments passed to geom_peak
#' @return A ggplot2 object
#' @export
#' @importFrom ggplot2 ggplot aes scale_fill_gradient
#' @examples
#' \dontrun{
#' track <- ez_peak("peaks.bed", "chr1:1000000-2000000", use_score = TRUE)
#' }
ez_peak <- function(data, region, color = "black", fill = "gray70",
                    alpha = 0.7, height = 0.8, use_score = FALSE, ...) {
  # Check if data is a file path or data frame
  if (is.character(data) && length(data) == 1) {
    # It's a file path, use peak_track
    return(peak_track(data, region, color = color, fill = fill,
                      alpha = alpha, height = height, use_score = use_score, ...))
  } else if (is.data.frame(data)) {
    # It's a data frame, create the plot directly
    if (use_score && "score" %in% colnames(data)) {
      p <- ggplot2::ggplot(data) +
        geom_peak(ggplot2::aes(xmin = start, xmax = end, fill = score),
                  color = color, alpha = alpha, height = height, ...) +
        ggplot2::scale_fill_gradient(low = "white", high = fill)
    } else {
      p <- ggplot2::ggplot(data) +
        geom_peak(ggplot2::aes(xmin = start, xmax = end),
                  color = color, fill = fill, alpha = alpha, height = height, ...)
    }

    # Apply the appropriate theme and scale
    p <- p + ez_peak_theme() +
      scale_x_genome_region(region) +
      ggplot2::ylim(0, 1)  # Fixed y-axis for peaks

    return(p)
  } else {
    stop("Data must be a file path or data frame")
  }
}

#' Easy gene track visualization
#'
#' This function creates a gene track visualization from a GTF/GFF file or data frame.
#' It is a wrapper around geom_gene that provides a simpler interface.
#'
#' @param data A GTF/GFF file path or data frame with gene data
#' @param region Genomic region to display (e.g., "chr1:1000000-2000000")
#' @param exon_height Height of exons (default: 0.75)
#' @param intron_height Height of introns (default: 0.4)
#' @param exon_color Color of exon borders (default: "black")
#' @param exon_fill Fill color of exons (default: "gray50")
#' @param intron_color Color of introns (default: "gray50")
#' @param gene_id Column name for gene ID (default: "gene_id")
#' @param gene_name Column name for gene name (default: "gene_name")
#' @param ... Additional arguments passed to geom_gene
#' @return A ggplot2 object
#' @export
#' @importFrom ggplot2 ggplot aes
#' @examples
#' \dontrun{
#' track <- ez_gene("genes.gtf", "chr1:1000000-2000000")
#' }
ez_gene <- function(data, region, exon_height = 0.75, intron_height = 0.4,
                    exon_color = "black", exon_fill = "gray50", intron_color = "gray50",
                    gene_id = "gene_id", gene_name = "gene_name", ...) {
  # Check if data is a file path or data frame
  if (is.character(data) && length(data) == 1) {
    # It's a file path, use gene_track
    return(gene_track(data, region, exon_height = exon_height, intron_height = intron_height,
                      exon_color = exon_color, exon_fill = exon_fill, intron_color = intron_color,
                      gene_id = gene_id, gene_name = gene_name, ...))
  } else if (is.data.frame(data)) {
    # It's a data frame, create the plot directly
    p <- ggplot2::ggplot(data) +
      geom_gene(ggplot2::aes(xstart = xstart, xend = xend, y = y, strand = strand),
                exon_height = exon_height, intron_height = intron_height,
                exon_color = exon_color, exon_fill = exon_fill,
                intron_color = intron_color, ...)

    # Apply the appropriate theme and scale
    p <- p + ez_gene_theme() +
      scale_x_genome_region(region)

    return(p)
  } else {
    stop("Data must be a file path or data frame")
  }
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
    return(interaction_track(data, region, curvature = curvature, color = color,
                            size = size, alpha = alpha, use_score = use_score, ...))
  } else if (is.data.frame(data)) {
    # It's a data frame, create the plot directly
    if (use_score && "score" %in% colnames(data)) {
      p <- ggplot2::ggplot(data) +
        geom_arc(ggplot2::aes(x = start1, y = 0, xend = start2, yend = 0, color = score),
                 curvature = curvature, size = size, alpha = alpha, ...) +
        ggplot2::scale_color_gradient(low = "blue", high = "red")
    } else {
      p <- ggplot2::ggplot(data) +
        geom_arc(ggplot2::aes(x = start1, y = 0, xend = start2, yend = 0),
                 curvature = curvature, color = color, size = size, alpha = alpha, ...)
    }

    # Apply the appropriate theme and scale
    p <- p + ez_peak_theme() +
      scale_x_genome_region(region) +
      ggplot2::ylim(-0.5, 0.5)  # Fixed y-axis for arcs

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
    return(hic_track(data, region, resolution = resolution, log_transform = log_transform,
                     low = low, high = high, ...))
  } else if (is.data.frame(data)) {
    # It's a data frame, create the plot directly
    p <- ggplot2::ggplot(data, ggplot2::aes(x = pos1, y = pos2, fill = count)) +
      geom_hic(low = low, high = high, ...) +
      ggplot2::coord_fixed()  # Ensure the plot is square

    # Apply the appropriate theme and scale
    p <- p + ez_theme() +
      scale_x_genome_region(region) +
      scale_x_genome_region(region)  # Same scale for y-axis

    return(p)
  } else {
    stop("Data must be a file path or data frame")
  }
}