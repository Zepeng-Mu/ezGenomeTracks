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
#' Strand separation: When `strand` is specified in aesthetics, genes are displayed on discrete tracks:
#' positive strand genes on one track, negative strand genes on another track.
#'
#' @inheritParams ggplot2::layer
#' @param exon_height Height of exons (default: 0.75)
#' @param intron_width Line width of gene body (default: 0.4)
#' @param arrow_length Length of directional arrows in inches (default: 0)
#' @param arrow_type Type of arrow head (default: "open")
#' @param exon_color Color of exon borders (default: "black", overridden by color mapping)
#' @param exon_fill Fill color of exons (default: "gray50", overridden by color mapping)
#' @param intron_color Color of intron/gene body (default: "gray50", overridden by color mapping)
#' @param na.rm If `TRUE`, silently drop `NA` values.
#' @return A ggplot2 layer that can be added to a plot.
#' @export
#' @importFrom ggplot2 GeomSegment GeomRect layer aes ggproto Geom arrow unit
#' @examples
#' \dontrun{
#' library(ggplot2)
#' p <- ggplot(gene_data) +
#'   geom_gene(aes(xstart = xstart, xend = xend, strand = strand, type = type, color = gene_name))
#' }
geom_gene <- function(mapping = NULL, data = NULL, stat = "identity",
                      position = "identity", ..., exon_height = 0.75,
                      intron_width = 0.4, arrow_length = 0,
                      arrow_type = "open", exon_color = "black",
                      exon_fill = "gray50", intron_color = "gray50",
                      na.rm = TRUE, show.legend = NA, inherit.aes = TRUE) {
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
  setup_data = function(data, params) {
    # Set y aesthetic based on strand for discrete positioning
    if ("strand" %in% names(data)) {
      # Use strand values directly as discrete y levels
      data$y <- ifelse(is.na(data$strand) | data$strand == "*", ".", as.character(data$strand))
    } else {
      # No strand provided, use single level
      data$y <- "."
    }
    data
  },
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
                        na.rm = FALSE) {
    # Separate data by type
    exon_data <- data[data$type == "exon", ]
    gene_data <- data[data$type == "gene", ]

    grobs <- list()

    # Draw exons (rectangles) for exon type
    if (nrow(exon_data) > 0) {
      if (all(c("exon_start", "exon_end") %in% names(exon_data))) {
        exons <- transform(exon_data,
          xmin = exon_start,
          xmax = exon_end,
          ymin = as.numeric(factor(y, levels = c(".", "-", "+"))) - exon_height/2,
          ymax = as.numeric(factor(y, levels = c(".", "-", "+"))) + exon_height/2
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
          ymin = as.numeric(factor(y, levels = c(".", "-", "+"))) - exon_height/2,
          ymax = as.numeric(factor(y, levels = c(".", "-", "+"))) + exon_height/2
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
      body_data <- transform(gene_data,
        x = xstart,
        xend = xend,
        y = as.numeric(factor(y, levels = c(".", "-", "+"))),
        yend = as.numeric(factor(y, levels = c(".", "-", "+")))
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
    return(result)
  } else {
    genes$xstart <- genes$start
    genes$xend <- genes$end
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
#' @return A data frame with gene and exon information formatted for geom_gene
#' @importFrom GenomicFeatures genes exonsBy
#' @importFrom IRanges subsetByOverlaps findOverlaps
#' @importFrom AnnotationDbi select
#' @examples
#' \dontrun{
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
#' region_gr <- parse_region("chr1:1000000-2000000")
#' gene_data <- extract_txdb_data(txdb, region_gr)
#' }
extract_txdb_data <- function(txdb, region_gr) {
  # Check dependencies
  if (!requireNamespace("GenomicFeatures", quietly = TRUE)) {
    stop("Package 'GenomicFeatures' is required for TxDb support. Install it with: BiocManager::install('GenomicFeatures')")
  }
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
    stop("Package 'AnnotationDbi' is required for TxDb support. Install it with: BiocManager::install('AnnotationDbi')")
  }

  # 1. Extract genes in region
  all_genes <- GenomicFeatures::genes(txdb)
  region_genes <- subsetByOverlaps(all_genes, region_gr)

  if (length(region_genes) == 0) {
    return(data.frame()) # empty data frame
  }

  # 2. Get gene symbols
  gene_ids <- names(region_genes)
  gene_symbols <- tryCatch({
    gene_info <- AnnotationDbi::select(txdb, keys = gene_ids,
                                       columns = c("GENEID", "SYMBOL"),
                                       keytype = "GENEID")
    # Create lookup: gene_id -> gene_symbol
    setNames(gene_info$SYMBOL, gene_info$GENEID)
  }, error = function(e) {
    # Fallback: use gene_id if symbols not available
    setNames(gene_ids, gene_ids)
  })

  # 3. Extract exons
  all_exons <- GenomicFeatures::exonsBy(txdb, by = "gene")
  region_exons <- all_exons[gene_ids]

  # 4. Build gene data frame
  gene_df <- data.frame(
    gene_id = gene_ids,
    gene_name = gene_symbols[gene_ids],
    strand = as.character(strand(region_genes)),
    start = start(region_genes),
    end = end(region_genes),
    type = "gene",
    stringsAsFactors = FALSE
  )

  # 5. Build exon data frame
  exon_list <- lapply(gene_ids, function(gid) {
    if (gid %in% names(region_exons)) {
      exons_gr <- region_exons[[gid]]
      data.frame(
        gene_id = gid,
        gene_name = gene_symbols[gid],
        strand = as.character(strand(exons_gr)),
        start = start(exons_gr),
        end = end(exons_gr),
        type = "exon",
        stringsAsFactors = FALSE
      )
    }
  })
  exon_df <- do.call(rbind, exon_list[!sapply(exon_list, is.null)])

  # 6. Combine and add required columns for geom_gene
  if (nrow(exon_df) > 0) {
    result <- rbind(gene_df, exon_df)
  } else {
    result <- gene_df
  }

  # Add required columns for geom_gene
  result$xstart <- result$start
  result$xend <- result$end

  return(result)
}


# globals used in examples/aes mappings within this file
utils::globalVariables(c(".data", "xstart", "xend", "strand"))
