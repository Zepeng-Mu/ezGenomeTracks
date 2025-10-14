#' Geom for gene model tracks
#'
#' This function creates a geom for gene model tracks, displaying genes with exons,
#' introns, and optional directional arrows. It is designed to work with gene
#' annotation data from GTF/GFF files after lightweight preprocessing.
#'
#' Expected columns in the data:
#' - `xstart`, `xend`, `type` (required)
#' - optional: `exon_start`, `exon_end` for long-format exon rows
#' - optional: `strand` for arrow direction ("+" or "-")
#'
#' The `type` column should contain "gene" for gene body lines or "exon" for exon rectangles.
#'
#' Color mapping: The `color` aesthetic is used for both exons (as fill) and introns (as color).
#'
#' Strand separation: When `strand` is specified in aesthetics, genes are displayed in two stacked tracks:
#' positive strand genes on top, negative strand genes on bottom. The `strand_spacing` parameter controls
#' the vertical gap between tracks.
#'
#' @inheritParams ggplot2::layer
#' @param exon_height Height of exons (default: 0.75)
#' @param intron_width Line width of gene body (default: 0.4)
#' @param arrow_length Length of directional arrows in inches (default: 0)
#' @param arrow_type Type of arrow head (default: "open")
#' @param exon_color Color of exon borders (default: "black", overridden by color mapping)
#' @param exon_fill Fill color of exons (default: "gray50", overridden by color mapping)
#' @param intron_color Color of intron/gene body (default: "gray50", overridden by color mapping)
#' @param strand_spacing Vertical spacing between positive and negative strand tracks (default: 0.2)
#' @param na.rm If `TRUE`, silently drop `NA` values.
#' @return A ggplot2 layer that can be added to a plot.
#' @export
#' @importFrom ggplot2 GeomSegment GeomRect layer aes ggproto Geom arrow unit
#' @examples
#' \dontrun{
#' library(ggplot2)
#' p <- ggplot(gene_data) +
#'   geom_gene(aes(xstart = xstart, xend = xend, y = y, strand = strand, type = type, color = gene_name))
#' }
geom_gene <- function(mapping = NULL, data = NULL, stat = "identity",
                      position = "identity", ..., exon_height = 0.75,
                      intron_width = 0.4, arrow_length = 0,
                      arrow_type = "open", exon_color = "black",
                      exon_fill = "gray50", intron_color = "gray50",
                      strand_spacing = 0.2, na.rm = TRUE, show.legend = NA, inherit.aes = TRUE) {
  # Provide defaults so users don't have to map y; draw at fixed vertical band
  default_aes <- ggplot2::aes(xstart = .data$xstart, xend = .data$xend, type = .data$type)
  if (is.null(mapping)) {
    mapping <- default_aes
  } else {
    mapping <- utils::modifyList(default_aes, as.list(mapping))
    mapping <- do.call(ggplot2::aes, mapping)
  }

  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomGene,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      exon_height = exon_height,
      intron_width = intron_width,
      arrow_length = arrow_length,
      arrow_type = arrow_type,
      exon_color = exon_color,
      exon_fill = exon_fill,
      intron_color = intron_color,
      strand_spacing = strand_spacing,
      na.rm = na.rm,
      ...
    )
  )
}

#' @rdname geom_gene
#' @format NULL
#' @usage NULL
GeomGene <- ggplot2::ggproto("GeomGene", Geom,
  required_aes = c("xstart", "xend", "type"),
  optional_aes = c("strand"),
  setup_params = function(data, params) {
    if (!is.null(params$arrow_length) && params$arrow_length > 0) {
      params$arrow <- ggplot2::arrow(
        type = params$arrow_type,
        length = ggplot2::unit(params$arrow_length, "inches")
      )
    } else {
      params$arrow <- NULL
    }
    params
  },
  draw_panel = function(data, panel_params, coord,
                        exon_height = 0.75, intron_width = 0.4,
                        arrow = NULL, exon_color = "black",
                        exon_fill = "gray50", intron_color = "gray50",
                        strand_spacing = 0.2, na.rm = FALSE) {
    # Separate data by type
    exon_data <- data[data$type == "exon", ]
    gene_data <- data[data$type == "gene", ]

    # Calculate y-offset based on strand for track separation
    if ("strand" %in% names(data)) {
      # Positive strand on top, negative strand on bottom
      exon_data$y_offset <- ifelse(exon_data$strand == "+",
                                   exon_height + strand_spacing,
                                   0)
      gene_data$y_offset <- ifelse(gene_data$strand == "+",
                                   exon_height + strand_spacing,
                                   0)
    } else {
      exon_data$y_offset <- 0
      gene_data$y_offset <- 0
    }

    grobs <- list()

    # Draw exons (rectangles) for exon type
    if (nrow(exon_data) > 0) {
      if (all(c("exon_start", "exon_end") %in% names(exon_data))) {
        exons <- transform(exon_data,
          xmin = exon_start,
          xmax = exon_end,
          ymin = 0 + y_offset,
          ymax = exon_height + y_offset
        )
        # Exons use fill from mapping (color aesthetic), no border color
        exons$colour <- NA
        # Use mapped color as fill, fallback to exon_fill parameter
        if ("colour" %in% names(exon_data)) {
          exons$fill <- exon_data$colour
        } else if (!is.null(exon_fill)) {
          exons$fill <- exon_fill
        }
        grobs[[length(grobs) + 1]] <- ggplot2::GeomRect$draw_panel(exons, panel_params, coord)
      } else {
        # If no exon_start/exon_end, use xstart/xend for exon boundaries
        exons <- transform(exon_data,
          xmin = xstart,
          xmax = xend,
          ymin = 0 + y_offset,
          ymax = exon_height + y_offset
        )
        # Exons use fill from mapping (color aesthetic), no border color
        exons$colour <- NA
        # Use mapped color as fill, fallback to exon_fill parameter
        if ("colour" %in% names(exon_data)) {
          exons$fill <- exon_data$colour
        } else if (!is.null(exon_fill)) {
          exons$fill <- exon_fill
        }
        grobs[[length(grobs) + 1]] <- ggplot2::GeomRect$draw_panel(exons, panel_params, coord)
      }
    }

    # Draw gene body (line) for gene type
    if (nrow(gene_data) > 0) {
      y_center <- exon_height / 2
      body_data <- transform(gene_data,
        x = xstart,
        xend = xend,
        y = y_center + y_offset,
        yend = y_center + y_offset
      )
      # Introns use color from mapping, fallback to intron_color parameter
      if ("colour" %in% names(gene_data)) {
        body_data$colour <- gene_data$colour
      } else if (!is.null(intron_color)) {
        body_data$colour <- intron_color
      }
      if (!"linewidth" %in% names(body_data)) body_data$linewidth <- intron_width

      grobs[[length(grobs) + 1]] <- ggplot2::GeomSegment$draw_panel(
        body_data, panel_params, coord,
        arrow = arrow
      )
    }

    do.call(grid::grobTree, grobs)
  },
  default_aes = aes(
    colour = "gray50", linewidth = 0.4, linetype = 1,
    alpha = 1
  )
)

#' Process gene annotation data for visualization
#'
#' This function processes gene annotation data from GTF/GFF files for visualization
#' with geom_gene. It extracts gene, transcript, and exon information and formats it
#' for use with ggplot2.
#'
#' @param gr A GRanges object with gene annotation data
#' @param gene_id Column name for gene ID (default: "gene_id")
#' @param gene_name Column name for gene name (default: "gene_name")
#' @param transcript_id Column name for transcript ID (default: "transcript_id")
#' @param type Column name for feature type (default: "type")
#' @return A data frame with gene, transcript, and exon information
#' @export
#' @importFrom GenomicRanges strand
#' @importFrom S4Vectors mcols
#' @examples
#' \dontrun{
#' library(rtracklayer)
#' gr <- import("genes.gtf")
#' gene_data <- process_gene_data(gr)
#' }
process_gene_data <- function(gr, gene_id = "gene_id", gene_name = "gene_name",
                              transcript_id = "transcript_id", type = "type") {
  # Accept either GRanges or a data.frame already in long format
  if (methods::is(gr, "GRanges")) {
    gr_df <- granges_to_df(gr)
    if (!"strand" %in% names(gr_df)) {
      gr_df$strand <- as.character(GenomicRanges::strand(gr))
    }
  } else if (is.data.frame(gr)) {
    gr_df <- gr
  } else {
    stop("Input must be a GRanges or a data.frame")
  }

  # Normalize strand
  if (!"strand" %in% names(gr_df)) gr_df$strand <- NA_character_
  gr_df$strand[is.na(gr_df$strand)] <- "*"

  # Ensure required id/type columns are present
  if (!all(c(gene_id, type) %in% colnames(gr_df))) {
    stop("Required columns not found in the data")
  }
  if (!gene_name %in% colnames(gr_df)) {
    gr_df[[gene_name]] <- gr_df[[gene_id]]
  }

  # Prefer explicit gene features; otherwise aggregate from transcript/exon
  prefer_genes <- any(gr_df[[type]] == "gene")

  if (prefer_genes) {
    genes <- gr_df[gr_df[[type]] == "gene", c(gene_id, gene_name, "strand", "start", "end")]
  } else {
    aggregate_gene_ranges <- function(df) {
      by_keys <- list(
        gene = df[[gene_id]],
        name = df[[gene_name]],
        strand = df[["strand"]]
      )
      start_min <- tapply(df$start, by_keys, min, na.rm = TRUE)
      end_max <- tapply(df$end, by_keys, max, na.rm = TRUE)
      idx <- which(!is.na(start_min) & !is.na(end_max), arr.ind = TRUE)
      if (length(idx) == 0) {
        return(data.frame())
      }
      dims <- dimnames(start_min)
      out <- data.frame(
        gene = dims[[1]][idx[, 1]],
        name = dims[[2]][idx[, 2]],
        strand = dims[[3]][idx[, 3]],
        start = as.numeric(start_min[idx]),
        end = as.numeric(end_max[idx]),
        stringsAsFactors = FALSE
      )
      names(out)[1:2] <- c(gene_id, gene_name)
      out
    }
    if (any(gr_df[[type]] == "transcript")) {
      genes <- aggregate_gene_ranges(gr_df[gr_df[[type]] == "transcript", ])
    } else if (any(gr_df[[type]] == "exon")) {
      genes <- aggregate_gene_ranges(gr_df[gr_df[[type]] == "exon", ])
    } else {
      stop("No gene, transcript, or exon features found in the data")
    }
  }

  if (nrow(genes) == 0) stop("No gene features could be constructed from the input")

  # Add type column to genes
  genes$type <- "gene"

  # Long-format exons if present
  has_exon <- any(gr_df[[type]] == "exon") && all(c("start", "end") %in% names(gr_df))
  if (has_exon) {
    exons <- gr_df[gr_df[[type]] == "exon", c(gene_id, gene_name, "strand", "start", "end")]
    names(exons)[names(exons) == "start"] <- "exon_start"
    names(exons)[names(exons) == "end"] <- "exon_end"
    exons$type <- "exon"
    result <- merge(genes, exons, by = c(gene_id, gene_name, "strand"), all.x = TRUE)
    result$xstart <- result$start
    result$xend <- result$end
    result$y <- result[[gene_name]]
    return(result)
  } else {
    genes$xstart <- genes$start
    genes$xend <- genes$end
    genes$y <- genes[[gene_name]]
    return(genes)
  }
}

#' Extract gene data from TxDb object
#'
#' This function extracts gene and exon information from a TxDb object
#' for a specific genomic region.
#'
#' @param txdb A TxDb object (e.g., TxDb.Hsapiens.UCSC.hg19.knownGene)
#' @param region_gr A GRanges object specifying the genomic region
#' @return A GRanges object with gene and exon information
#' @importFrom IRanges subsetByOverlaps findOverlaps
#' @examples
#' \dontrun{
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' region_gr <- parse_region("chr1:1000000-2000000")
#' gene_data <- extract_txdb_data(txdb, region_gr)
#' }
extract_txdb_data <- function(txdb, region_gr) {
  # Check if GenomicFeatures is available
  if (!requireNamespace("GenomicFeatures", quietly = TRUE)) {
    stop("Package 'GenomicFeatures' is required for TxDb support. Install it with: BiocManager::install('GenomicFeatures')")
  }

  # Extract genes in the region
  all_genes <- GenomicFeatures::genes(txdb)
  region_genes <- subsetByOverlaps(all_genes, region_gr)

  # Extract exons by gene
  all_exons <- GenomicFeatures::exonsBy(txdb, by = "gene")

  # Filter exons for genes in the region
  gene_ids <- names(region_genes)
  region_exons <- all_exons[gene_ids]

  # Convert to a format similar to GTF data
  gene_list <- list()
  exon_list <- list()

  for (i in seq_along(region_genes)) {
    gene_id <- names(region_genes)[i]
    gene_gr <- region_genes[i]

    # Add gene information
    gene_gr$type <- "gene"
    gene_gr$gene_id <- gene_id
    gene_gr$gene_name <- gene_id # Use gene_id as gene_name for TxDb
    gene_list[[i]] <- gene_gr

    # Add exon information if available
    if (gene_id %in% names(region_exons)) {
      exons_gr <- region_exons[[gene_id]]
      exons_gr$type <- "exon"
      exons_gr$gene_id <- gene_id
      exons_gr$gene_name <- gene_id
      exon_list[[length(exon_list) + 1]] <- exons_gr
    }
  }

  # Combine genes and exons
  result_list <- c(gene_list, exon_list)
  if (length(result_list) > 0) {
    result_gr <- do.call(c, result_list)
  } else {
    # Return empty GRanges if no genes found
    result_gr <- GenomicRanges::GRanges()
  }

  return(result_gr)
}

#' Create a gene track from a GTF/GFF file or TxDb object
#'
#' This function creates a gene track from either a GTF/GFF file or a TxDb object.
#' It imports the data for a specific region and creates a ggplot2 layer for visualization.
#'
#' @param source Either a path to a GTF/GFF file (character) or a TxDb object
#' @param region Genomic region to display (e.g., "chr1:1000000-2000000")
#' @param exon_height Height of exons (default: 0.75)
#' @param intron_width Width of introns (default: 0.4)
#' @param exon_color Color of exon borders (default: "black")
#' @param exon_fill Fill color of exons (default: "gray50")
#' @param intron_color Color of introns (default: "gray50")
#' @param strand_spacing Vertical spacing between positive and negative strand tracks (default: 0.2)
#' @param gene_id Column name for gene ID (default: "gene_id")
#' @param gene_name Column name for gene name (default: "gene_name")
#' @param ... Additional arguments passed to geom_gene
#' @return A ggplot2 layer
#' @export
#' @importFrom ggplot2 ggplot aes
#' @importFrom methods is
#' @examples
#' \dontrun{
#' # Using a GTF file
#' p1 <- gene_track("genes.gtf", "chr1:1000000-2000000")
#'
#' # Using a TxDb object
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' p2 <- gene_track(txdb, "chr1:1000000-2000000")
#' }
gene_track <- function(source, region, exon_height = 0.75, intron_width = 0.4,
                       exon_color = "black", exon_fill = "gray50", intron_color = "gray50",
                       strand_spacing = 0.2, gene_id = "gene_id", gene_name = "gene_name", ...) {
  # Parse the region
  region_gr <- parse_region(region)

  # Determine if source is a file path or TxDb object
  if (is.character(source)) {
    # Import from file
    gene_gr <- rtracklayer::import(source, which = region_gr)
  } else if (methods::is(source, "TxDb")) {
    # Check if GenomicFeatures is available for TxDb support
    if (!requireNamespace("GenomicFeatures", quietly = TRUE)) {
      stop("Package 'GenomicFeatures' is required for TxDb support. Install it with: BiocManager::install('GenomicFeatures')")
    }
    # Extract from TxDb object
    gene_gr <- extract_txdb_data(source, region_gr)
  } else {
    stop("Source must be either a file path (character) or a TxDb object")
  }

  # Process the gene data
  gene_data <- process_gene_data(gene_gr, gene_id = gene_id, gene_name = gene_name)

  # Create the plot
  p <- ggplot2::ggplot(gene_data) +
    geom_gene(ggplot2::aes(xstart = .data$xstart, xend = .data$xend, strand = .data$strand),
      exon_height = exon_height, intron_width = intron_width,
      exon_color = exon_color, exon_fill = exon_fill,
      intron_color = intron_color, strand_spacing = strand_spacing, ...
    )

  # Apply the appropriate theme and scale
  p <- p + ez_gene_theme() +
    scale_x_genome_region(region)

  return(p)
}

# globals used in examples/aes mappings within this file
utils::globalVariables(c(".data", "xstart", "xend", "y", "strand"))
