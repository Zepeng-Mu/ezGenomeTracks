#' Geom for gene model tracks
#'
#' This function creates a geom for gene model tracks, displaying genes with exons,
#' introns, and directional arrows. It is designed to work with gene annotation data
#' from GTF/GFF files.
#'
#' @inheritParams ggplot2::geom_segment
#' @param exon_height Height of exons (default: 0.75)
#' @param intron_height Height of introns (default: 0.4)
#' @param arrow_length Length of directional arrows (default: 0.1)
#' @param arrow_type Type of arrow (default: "open")
#' @param exon_color Color of exon borders (default: "black")
#' @param exon_fill Fill color of exons (default: "gray50")
#' @param intron_color Color of introns (default: "gray50")
#' @param text_size Size of gene labels (default: 3)
#' @return A ggplot2 layer
#' @export
#' @importFrom ggplot2 geom_rect geom_segment geom_text arrow unit aes
#' @examples
#' \dontrun{
#' library(ggplot2)
#' p <- ggplot(gene_data) + geom_gene(aes(xstart = start, xend = end, y = gene_name, strand = strand))
#' }
geom_gene <- function(mapping = NULL, data = NULL, stat = "identity",
                       position = "identity", ..., exon_height = 0.75,
                       intron_height = 0.4, arrow_length = 0.1,
                       arrow_type = "open", exon_color = "black",
                       exon_fill = "gray50", intron_color = "gray50",
                       text_size = 3, show.legend = NA, inherit.aes = TRUE) {

  # This is a composite geom that will return a list of layers
  structure(
    list(
      # Function to create the layers when the plot is built
      layer_function = function(self, plot) {
        # Extract the data and mapping from the plot
        data <- plot$data
        mapping <- plot$mapping

        # Check if we have the required aesthetics
        required_aes <- c("xstart", "xend", "y", "strand")
        missing_aes <- setdiff(required_aes, names(mapping))
        if (length(missing_aes) > 0) {
          stop("Missing required aesthetics: ", paste(missing_aes, collapse = ", "))
        }

        # Create layers for introns (gene body)
        intron_layer <- ggplot2::geom_segment(
          data = data,
          mapping = ggplot2::aes(
            x = .data[[as.character(mapping$xstart)]],
            xend = .data[[as.character(mapping$xend)]],
            y = .data[[as.character(mapping$y)]],
            yend = .data[[as.character(mapping$y)]]
          ),
          size = intron_height,
          color = intron_color,
          arrow = ggplot2::arrow(
            type = arrow_type,
            length = ggplot2::unit(arrow_length, "inches")
          )
        )

        # Create layers for exons if available
        if ("exon_start" %in% colnames(data) && "exon_end" %in% colnames(data)) {
          exon_layer <- ggplot2::geom_rect(
            data = data,
            mapping = ggplot2::aes(
              xmin = exon_start,
              xmax = exon_end,
              ymin = .data[[as.character(mapping$y)]] - exon_height/2,
              ymax = .data[[as.character(mapping$y)]] + exon_height/2
            ),
            fill = exon_fill,
            color = exon_color
          )

          return(list(intron_layer, exon_layer))
        } else {
          # If no exon information, just return the intron layer
          return(list(intron_layer))
        }
      }
    ),
    class = c("GeomGene", "Geom", "ggproto")
  )
}

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
  # Convert GRanges to data frame
  gr_df <- granges_to_df(gr)

  # Check if required columns exist
  if (!all(c(gene_id, type) %in% colnames(gr_df))) {
    stop("Required columns not found in the data")
  }

  # Use gene_id as gene_name if gene_name is not available
  if (!gene_name %in% colnames(gr_df)) {
    gr_df[[gene_name]] <- gr_df[[gene_id]]
  }

  # Extract gene-level information
  genes <- gr_df[gr_df[[type]] == "gene", ]

  # If no genes found, try to infer gene boundaries from transcripts or exons
  if (nrow(genes) == 0) {
    if ("transcript" %in% gr_df[[type]]) {
      # Group by gene and get min start and max end
      transcripts <- gr_df[gr_df[[type]] == "transcript", ]
      genes <- aggregate(
        cbind(start, end) ~ get(gene_id) + get(gene_name) + strand,
        data = transcripts,
        FUN = function(x) c(min(x), max(x))
      )
      colnames(genes)[1:2] <- c(gene_id, gene_name)
      genes$start <- genes$start[, 1]
      genes$end <- genes$end[, 2]
    } else if ("exon" %in% gr_df[[type]]) {
      # Group by gene and get min start and max end
      exons <- gr_df[gr_df[[type]] == "exon", ]
      genes <- aggregate(
        cbind(start, end) ~ get(gene_id) + get(gene_name) + strand,
        data = exons,
        FUN = function(x) c(min(x), max(x))
      )
      colnames(genes)[1:2] <- c(gene_id, gene_name)
      genes$start <- genes$start[, 1]
      genes$end <- genes$end[, 2]
    } else {
      stop("No gene, transcript, or exon features found in the data")
    }
  }

  # Extract exon information if available
  if ("exon" %in% gr_df[[type]]) {
    exons <- gr_df[gr_df[[type]] == "exon", ]

    # Merge exons with genes
    result <- merge(genes, exons, by = c(gene_id, gene_name), suffixes = c("", "_exon"))

    # Rename columns for use with geom_gene
    result$xstart <- result$start
    result$xend <- result$end
    result$exon_start <- result$start_exon
    result$exon_end <- result$end_exon
    result$y <- result[[gene_name]]

    return(result)
  } else {
    # If no exons, just return gene information
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
#' @importFrom GenomicFeatures genes exonsBy
#' @importFrom GenomicRanges findOverlaps subsetByOverlaps
#' @importFrom IRanges subsetByOverlaps
#' @examples
#' \dontrun{
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' region_gr <- parse_region("chr1:1000000-2000000")
#' gene_data <- extract_txdb_data(txdb, region_gr)
#' }
extract_txdb_data <- function(txdb, region_gr) {
  # Extract genes in the region
  all_genes <- GenomicFeatures::genes(txdb)
  region_genes <- GenomicRanges::subsetByOverlaps(all_genes, region_gr)

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
    gene_gr$gene_name <- gene_id  # Use gene_id as gene_name for TxDb
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
#' @param intron_height Height of introns (default: 0.4)
#' @param exon_color Color of exon borders (default: "black")
#' @param exon_fill Fill color of exons (default: "gray50")
#' @param intron_color Color of introns (default: "gray50")
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
gene_track <- function(source, region, exon_height = 0.75, intron_height = 0.4,
                       exon_color = "black", exon_fill = "gray50", intron_color = "gray50",
                       gene_id = "gene_id", gene_name = "gene_name", ...) {
  # Parse the region
  region_gr <- parse_region(region)

  # Determine if source is a file path or TxDb object
  if (is.character(source)) {
    # Import from file
    gene_gr <- rtracklayer::import(source, which = region_gr)
  } else if (methods::is(source, "TxDb")) {
    # Extract from TxDb object
    gene_gr <- extract_txdb_data(source, region_gr)
  } else {
    stop("Source must be either a file path (character) or a TxDb object")
  }

  # Process the gene data
  gene_data <- process_gene_data(gene_gr, gene_id = gene_id, gene_name = gene_name)

  # Create the plot
  p <- ggplot2::ggplot(gene_data) +
    geom_gene(ggplot2::aes(xstart = xstart, xend = xend, y = y, strand = strand),
              exon_height = exon_height, intron_height = intron_height,
              exon_color = exon_color, exon_fill = exon_fill,
              intron_color = intron_color, ...)

  # Apply the appropriate theme and scale
  p <- p + ez_gene_theme() +
    scale_x_genome_region(region)

  return(p)
}
