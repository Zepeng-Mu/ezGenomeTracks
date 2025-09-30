#' Easy signal track visualization
#'
#' This function creates a signal track visualization from a bigWig file, data frame, or list of data sources.
#' It is a wrapper around geom_signal that provides a simpler interface with support for grouping and multiple tracks.
#'
#' @param input A bigWig file path, data frame with signal data, or list of data frames/bigWig files
#' @param region Genomic region to display (e.g., "chr1:1000000-2000000")
#' @param type Type of signal visualization: "line", "area", or "heatmap" (default: "area")
#' @param color Line color (default: "steelblue"). Can be a vector for multiple colors.
#' @param fill Fill color for area plots (default: "steelblue"). Can be a vector for multiple colors.
#' @param stack Whether to stack multiple tracks or combine into one (default: TRUE)
#' @param group Column name for grouping data within a single data frame (default: NULL)
#' @param colors Color palette for groups/tracks. If NULL, uses default palette (default: NULL)
#' @param y_axis_style Y-axis style: "none", "simple", or "full" (default: "none")
#' @param y_range Y-axis range limits (default: NULL)
#' @param alpha Transparency (default: 0.5)
#' @param bin_width Width of bins in base pairs (default: NULL)
#' @param ... Additional arguments passed to geom_signal
#' @return A ggplot2 object
#' @export
#' @importFrom ggplot2 ggplot aes scale_y_continuous coord_cartesian labs facet_wrap scale_color_manual scale_fill_manual
#' @examples
#' \dontrun{
#' # Single track
#' track1 <- ez_signal("signal.bw", "chr1:1000000-2000000")
#' 
#' # Grouped data frame
#' track2 <- ez_signal(df, "chr1:1000000-2000000", group = "sample", stack = FALSE)
#' 
#' # Multiple data frames
#' track3 <- ez_signal(list(sample1 = df1, sample2 = df2), "chr1:1000000-2000000")
#' }
ez_signal <- function(input, region, type = c("area", "line", "heatmap"),
                      color = "steelblue", fill = "steelblue", stack = TRUE,
                      group = NULL, colors = NULL,
                      y_axis_style = c("none", "simple", "full"),
                      y_range = NULL, alpha = 0.5, bin_width = NULL, ...) {
  # Validate inputs
  type <- match.arg(type)
  y_axis_style <- match.arg(y_axis_style)

  stopifnot(
    "alpha must be between 0 and 1" = alpha >= 0 && alpha <= 1,
    "region must be provided" = !missing(region),
    "bin_width must be positive integer" = is.null(bin_width) || bin_width > 0 && is.integer(bin_width),
    "stack must be logical" = is.logical(stack)
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

  # Handle single element input
  if (length(input) == 1) {
    if (is(input, "data.frame")) {
      # Single data frame
      if (!all(c("start", "score") %in% colnames(input))) {
        stop("Data frame must contain 'start' and 'score' columns")
      }
      plotDt <- input
    } else if (is(input, "character")) {
      # Single file path
      if (!file.exists(input)) stop("File does not exist: ", input)
      # TODO read bigWig file
    }
  }

  # Handle list input

  if (is(input, "character")) {
    if (length(input) == 1) {
      # Single file path
      if (!file.exists(input)) stop("File does not exist: ", input)
      # TODO: Implement bigWig file reading
      stop("BigWig file support not yet implemented. Please use data frames.")
    } else {
      # Multiple file paths
      # TODO: Implement multiple bigWig file reading and plot as overlapping tracks
    }
  }
  
  # Handle list of data frames or file paths
  if (is(input, "list")) {
    if (is.null(names(input))) {
      names(input) <- paste0("Track ", seq_along(input))
    }
    
    # Process each element in the list
    plot_data_list <- list()
    for (i in seq_along(input)) {
      track_name <- names(input)[i]
      track_data <- input[[i]]
      
      if (is.character(track_data)) {
        if (length(track_data) == 1) {
          # Single file path
          if (!file.exists(track_data)) stop("File does not exist: ", track_data)
          # TODO: Implement bigWig file reading
          stop("BigWig file support not yet implemented. Please use data frames.")
        } else {
          # Vector of file paths - create grouped data
          # TODO: Implement multiple bigWig file reading
          stop("Multiple bigWig file support not yet implemented. Please use data frames.")
        }
      } else if (is.data.frame(track_data)) {
        # Validate required columns
        if (!all(c("start", "score") %in% colnames(track_data))) {
          stop("Data frame must contain 'start' and 'score' columns")
        }
        
        # Add track identifier
        track_data$track <- track_name
        
        # Handle grouping within this track
        if (!is.null(group) && group %in% colnames(track_data)) {
          track_data$group_var <- track_data[[group]]
        } else {
          track_data$group_var <- "default"
        }
        
        plot_data_list[[i]] <- track_data
      } else {
        stop("List elements must be file paths or data frames")
      }
    }
    
    # Combine all data
    plotDt <- do.call(rbind, plot_data_list)
    
    # Create the plot
    if (!is.null(group) && group %in% colnames(plotDt)) {
      # Multiple tracks with grouping
      unique_groups <- unique(plotDt$group_var)
      unique_tracks <- unique(plotDt$track)
      n_colors <- length(unique_groups)
      
      if (is.null(colors)) {
        colors <- get_default_colors(n_colors)
        names(colors) <- unique_groups
      }
      
      if (stack) {
        # Stacked tracks with grouped signals within each
        p <- ggplot2::ggplot(plotDt, ggplot2::aes(x = start, y = score, color = group_var, fill = group_var)) +
          geom_signal(type = type, alpha = alpha, ...) +
          ggplot2::facet_wrap(~ track, ncol = 1, scales = "free_y") +
          ggplot2::scale_color_manual(values = colors, name = group) +
          ggplot2::scale_fill_manual(values = colors, name = group)
      } else {
        # Overlapping tracks and groups
        p <- ggplot2::ggplot(plotDt, ggplot2::aes(x = start, y = score, color = interaction(group_var, track), fill = interaction(group_var, track))) +
          geom_signal(type = type, alpha = alpha, ...)
        
        # Create combined color palette
        combined_groups <- unique(interaction(plotDt$group_var, plotDt$track))
        n_combined <- length(combined_groups)
        if (is.null(colors)) {
          combined_colors <- get_default_colors(n_combined)
        } else {
          combined_colors <- rep(colors, length(unique_tracks))
        }
        names(combined_colors) <- combined_groups
        
        p <- p + 
          ggplot2::scale_color_manual(values = combined_colors, name = "Track.Group") +
          ggplot2::scale_fill_manual(values = combined_colors, name = "Track.Group")
      }
    } else {
      # Multiple tracks without grouping
      unique_tracks <- unique(plotDt$track)
      n_colors <- length(unique_tracks)
      
      if (is.null(colors)) {
        colors <- get_default_colors(n_colors)
        names(colors) <- unique_tracks
      }
      
      if (stack) {
        # Stacked tracks
        p <- ggplot2::ggplot(plotDt, ggplot2::aes(x = start, y = score)) +
          geom_signal(type = type, color = color, fill = fill, alpha = alpha, ...) +
          ggplot2::facet_wrap(~ track, ncol = 1, scales = "free_y")
      } else {
        # Overlapping tracks
        p <- ggplot2::ggplot(plotDt, ggplot2::aes(x = start, y = score, color = track, fill = track)) +
          geom_signal(type = type, alpha = alpha, ...) +
          ggplot2::scale_color_manual(values = colors, name = "Track") +
          ggplot2::scale_fill_manual(values = colors, name = "Track")
      }
    }
    
  } else if (is(input, "data.frame")) {
    # Single data frame
    # Validate required columns
    if (!all(c("start", "score") %in% colnames(input))) {
      stop("Data frame must contain 'start' and 'score' columns")
    }
    
    plotDt <- input
    
    if (!is.null(group) && group %in% colnames(plotDt)) {
      # Grouped data frame
      unique_groups <- unique(plotDt[[group]])
      n_colors <- length(unique_groups)
      
      if (is.null(colors)) {
        colors <- get_default_colors(n_colors)
        names(colors) <- unique_groups
      }
      
      if (stack) {
        # Stacked groups using facet_wrap
        p <- ggplot2::ggplot(plotDt, ggplot2::aes(x = start, y = score)) +
          geom_signal(type = type, color = color, fill = fill, alpha = alpha, ...) +
          ggplot2::facet_wrap(as.formula(paste("~", group)), ncol = 1, scales = "free_y")
      } else {
        # Overlapping groups
        aes_mapping <- ggplot2::aes(x = start, y = score)
        aes_mapping$color <- as.symbol(group)
        aes_mapping$fill <- as.symbol(group)
        
        p <- ggplot2::ggplot(plotDt, aes_mapping) +
          geom_signal(type = type, alpha = alpha, ...) +
          ggplot2::scale_color_manual(values = colors, name = group) +
          ggplot2::scale_fill_manual(values = colors, name = group)
      }
    } else {
      # Single track without grouping
      p <- ggplot2::ggplot(plotDt, ggplot2::aes(x = start, y = score)) +
        geom_signal(type = type, color = color, fill = fill, alpha = alpha, ...)
    }
  } else {
    stop("Input must be a file path, data frame, or list of data frames/file paths")
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

#' @importFrom ezGenomeTracks geom_manhattan ez_manhattan

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
ez_gene <- function(data, region, exon_height = 0.75, intron_height = 0.4,
                    exon_color = "black", exon_fill = "gray50", intron_color = "gray50",
                    gene_id = "gene_id", gene_name = "gene_name", ...) {
  # Check if data is a file path, TxDb object, or data frame
  if ((is.character(data) && length(data) == 1) || methods::is(data, "TxDb")) {
    # It's a file path or TxDb object, use gene_track
    return(gene_track(data, region,
      exon_height = exon_height, intron_height = intron_height,
      exon_color = exon_color, exon_fill = exon_fill, intron_color = intron_color,
      gene_id = gene_id, gene_name = gene_name, ...
    ))
  } else if (is.data.frame(data)) {
    # It's a data frame, create the plot directly
    p <- ggplot2::ggplot(data) +
      geom_gene(ggplot2::aes(xstart = xstart, xend = xend, y = y, strand = strand),
        exon_height = exon_height, intron_height = intron_height,
        exon_color = exon_color, exon_fill = exon_fill,
        intron_color = intron_color, ...
      )

    # Apply the appropriate theme and scale
    p <- p + ez_gene_theme() +
      scale_x_genome_region(region)

    return(p)
  } else {
    stop("Data must be a file path, TxDb object, or data frame")
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
