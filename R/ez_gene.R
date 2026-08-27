# ezGenomeTracks - Gene functions (split from ez_wrappers.R)
#' Easy gene track visualization
#'
#' @description
#' This function creates a gene track visualization from genomic annotations,
#' supporting various input formats including GTF/GFF files, TxDb objects, and
#' data frames. It automatically handles gene structure visualization with
#' exons, introns, and strand information.
#'
#' By default, when `y = "strand"`, genes are colored by strand: plus strand
#' uses `"darkgreen"` and minus strand uses `"orange2"`. To use uniform colors
#' instead, explicitly set `exon_color`, `exon_fill`, and `intron_color`.
#'
#' @param input Input data source, which can be:
#'   - A file path to a GTF/GFF file
#'   - A TxDb object from the GenomicFeatures package
#'   - A data frame with gene annotation data
#' @param region Genomic region to display in the format "chr:start-end".
#'   Example: "chr1:1000000-2000000". Either `region` or `gene` (with `gene_db`)
#'   must be provided.
#' @param gene Gene name/symbol to look up (e.g., "PTPRC", "TP53").
#'   When provided, the region is automatically determined from the gene coordinates.
#'   If `gene_db` is not specified, `input` is used for coordinate lookup (recommended).
#'   Either `region` or `gene` must be provided.
#' @param gene_db Optional. TxDb object, GTF/GFF file path, or data frame for gene coordinate
#'   lookup when using `gene` parameter. If NULL (default), uses `input` for lookup.
#'   Only required if you want to use a different annotation source than `input`.
#' @param org_db Optional OrgDb object for gene symbol mapping. If NULL (default),
#'   auto-detects available OrgDb packages.
#' @param extend Numeric. Amount to extend the region beyond the gene body when
#'   using `gene` parameter. Default: 0.1 (10% of gene length on each side).
#' @param extend_type How to interpret `extend`: "proportion" (relative to gene
#'   length) or "bp" (absolute base pairs). Default: "proportion".
#' @param exon_height Relative height of exons (0 to 1). Default: 0.2
#' @param intron_width Line width for introns. Default: 0.6
#' @param exon_color Border color for exons. Default: NULL (uses strand-based colors
#'   when `y = "strand"`, otherwise "gray50")
#' @param exon_fill Fill color for exons. Default: NULL (uses strand-based colors
#'   when `y = "strand"`, otherwise "gray50")
#' @param intron_color Color for intron lines. Default: NULL (uses strand-based colors
#'   when `y = "strand"`, otherwise "gray50")
#' @param gene_id Column name for gene identifiers. Default: "gene_id"
#' @param gene_name Column name for gene symbols/names. Default: "gene_name"
#' @param y Column name for the y-axis grouping variable. Default: "strand"
#' @param label Column name to use for text labels. If NULL (default), no labels
#'   are displayed. Set to a column name (e.g., "gene_name") to show labels.
#' @param label_size Size of text labels. Default: 3
#' @param label_color Color of text labels. If NULL (default), uses strand-based
#'   colors when `y = "strand"`, otherwise uses exon_fill color.
#' @param label_style Strategy for handling overlapping labels. Options:
#'   - "auto" (default): Uses ggrepel if available, otherwise check_overlap
#'   - "simple": Standard geom_text with no overlap handling
#'   - "repel": Force use of ggrepel (errors if not installed)
#'   - "none": No labels displayed
#' @param max_labels Maximum number of labels to display. NULL (default) shows all.
#'   When set, labels are filtered based on label_priority.
#' @param label_priority Priority criterion for filtering labels when max_labels is set.
#'   Options: "length" (default, prioritizes longer genes), "name" (alphabetical),
#'   or a column name in the input to sort by.
#' @param repel_args Named list of additional arguments passed to geom_text_repel()
#'   when label_style = "repel" or "auto" (with ggrepel installed). Default behavior
#'   uses horizontal-only repositioning (`direction = "x"`) with no connecting lines
#'   (`segment.color = NA`). To show connecting lines, use `list(segment.color = "gray50")`.
#'   Override direction with `list(direction = "both")` for vertical repositioning too.
#'   Other useful options: `max.overlaps`, `force`, `box.padding`, `point.padding`.
#' @param border Logical. If TRUE, adds a black border around the plot panel.
#'   Default: FALSE
#' @param label_chr Logical. If `TRUE` (default), labels the x-axis with the chromosome name
#'   (e.g., "Chr1"). Set to `FALSE` to suppress the x-axis label.
#' @param auto_import Deprecated. GTF/GFF file paths are always imported via
#'   [import_gtf()] to provide a standardized data-frame schema. This parameter
#'   is retained for backward compatibility and is ignored.
#' @param ... Additional arguments passed to `geom_gene()`. Note that `color`
#'   and `colour` arguments are ignored; use `exon_color`, `exon_fill`, and
#'   `intron_color` instead.
#'
#' @return A ggplot2 object representing the gene track.
#'
#' @details
#' The function automatically processes different input types:
#' - For GTF/GFF files: Uses [import_gtf()] to create a standardized annotation data frame
#' - For TxDb objects: Extracts gene models using GenomicFeatures
#' - For data frames: Expects columns for chromosome, start, end, strand, and type
#'
#' The visualization includes:
#' - Exons as filled rectangles
#' - Introns as connecting lines
#' - Strand information with arrowheads
#' - Automatic y-axis separation by the specified y variable
#'
#' **Label Overlap Handling:**
#'
#' The function provides flexible strategies for managing overlapping gene labels:
#' - `label_style = "auto"`: Automatically uses ggrepel if installed, otherwise
#'   applies check_overlap to hide overlapping labels
#' - `label_style = "simple"`: Standard text labels with no overlap handling
#' - `label_style = "repel"`: Uses ggrepel to reposition labels. By default, labels
#'   are repositioned **horizontally only** (`direction = "x"`) with no connecting
#'   lines (`segment.color = NA`) to maintain a clean appearance while keeping
#'   labels horizontally aligned. This can be changed via `repel_args`.
#' - `label_style = "none"`: Disables all labels
#'
#' When many genes are present, use `max_labels` to limit the number of labels shown,
#' prioritized by `label_priority` (gene length by default).
#'
#' @export
#' @importFrom ggplot2 ggplot aes scale_colour_identity scale_fill_identity
#' @importFrom GenomicFeatures exonsBy transcriptsBy
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
#' @importFrom rtracklayer import
#' @importFrom methods is
#'
#' @examples
#' # From a data frame
#' data(example_genes)
#' ez_gene(example_genes, "chr1:11869-14409")
#'
#' # Limit labels to top 5 longest genes
#' ez_gene(example_genes, "chr1:42100000-42700000", max_labels = 5)
#'
#' # Use ggrepel for smart label positioning (if installed)
#' \dontrun{
#' # Default: horizontal-only repositioning, no connecting lines
#' ez_gene(example_genes, "chr1:42100000-42700000", label_style = "repel")
#'
#' # Show connecting lines to original position
#' ez_gene(example_genes, "chr1:42100000-42700000",
#'         label_style = "repel",
#'         repel_args = list(segment.color = "gray50"))
#'
#' # Allow both horizontal and vertical repositioning
#' ez_gene(example_genes, "chr1:42100000-42700000",
#'         label_style = "repel",
#'         repel_args = list(direction = "both"))
#'
#' # Custom repel settings for denser regions
#' ez_gene(example_genes, "chr1:42100000-42700000",
#'         label_style = "repel",
#'         repel_args = list(max.overlaps = 30, force = 3))
#' }
#'
#' # Hide overlapping labels automatically
#' ez_gene(example_genes, "chr1:42100000-42700000", label_style = "auto")
#'
#' # No labels
#' ez_gene(example_genes, "chr1:11869-14409", label_style = "none")
#'
#' \dontrun{
#' # Using gene name for region lookup
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
#' ez_gene(txdb, gene = "PTPRC", gene_db = txdb)
#' }
ez_gene <- function(
  input,
  region = NULL,
  gene = NULL,
  gene_db = NULL,
  org_db = NULL,
  extend = 0.1,
  extend_type = c("proportion", "bp"),
  exon_height = 0.2,
  intron_width = 0.6,
  exon_color = NULL,
  exon_fill = NULL,
  intron_color = NULL,
  gene_id = "gene_id",
  gene_name = "gene_name",
  y = "strand",
  label = "gene_name",
  label_size = 3,
  label_color = NULL,
  label_style = c("auto", "simple", "repel", "none"),
  max_labels = NULL,
  label_priority = "length",
  repel_args = list(),
  border = FALSE,
  label_chr = TRUE,
  auto_import = TRUE,
  ...
) {
  # Consolidated GTF/GFF path: always standardize through import_gtf()
  if (is.character(input) && length(input) == 1) {
    if (grepl("\\.(gtf|gff)(\\.(gz|bz2|xz))?$", input, ignore.case = TRUE)) {
      input <- import_gtf(input, region = region, verbose = FALSE)
    }
  }

  # CONSOLIDATE DATA AND GENE_DB: If gene is provided but gene_db is NULL,
  # use input as the lookup source (if it's a data frame or can be used as one)
  if (!is.null(gene) && is.null(gene_db)) {
    if (is.data.frame(input)) {
      gene_db <- input
    } else if (is.character(input) && length(input) == 1) {
      # If input is still a character file path (non-GTF), we can't use it for gene lookup
      stop(
        "When using 'gene' parameter with a file path in 'input', ",
        "provide 'gene_db' explicitly (TxDb, GTF/GFF path, or data frame from import_gtf())"
      )
    } else if (methods::is(input, "TxDb")) {
      # If input is a TxDb, use it for both gene models and lookup
      gene_db <- input
    }
  }

  # Resolve region from either region string or gene name
  region <- .resolve_region(
    region = region,
    gene = gene,
    gene_db = gene_db,
    org_db = org_db,
    extend = extend,
    extend_type = extend_type
  )

  # Match label_style argument
  label_style <- match.arg(label_style)
  # Filter out 'color' and 'colour' from ... to prevent them from overriding

  # exon_color, exon_fill, and intron_color

  dots <- list(...)
  dots <- dots[!names(dots) %in% c("color", "colour")]

  # Determine if we should use strand-based coloring

  # This is the default when y = "strand" and no explicit colors are provided

  use_strand_colors <- y == "strand" &&
    is.null(exon_color) &&
    is.null(exon_fill) &&
    is.null(intron_color)

  # Define strand color palette
  strand_colors <- c("+" = "green4", "-" = "orange2", "Unknown" = "gray50")

  # Set default colors if not using strand-based coloring

  if (!use_strand_colors) {
    if (is.null(exon_color)) {
      exon_color <- "gray50"
    }
    if (is.null(exon_fill)) {
      exon_fill <- "gray50"
    }
    if (is.null(intron_color)) intron_color <- "gray50"
  }

  # Parse the region
  region_gr <- parse_region(region)
  chr <- gsub("^chr", "", as.character(GenomicRanges::seqnames(region_gr)), ignore.case = TRUE)

  # Extract region limits for clipping
  region_limits <- c(
    GenomicRanges::start(region_gr),
    GenomicRanges::end(region_gr)
  )

  # Process input based on input type
  if (is.character(input) && length(input) == 1) {
    # Non-GTF file path: use generic importer and convert downstream
    gene_gr <- rtracklayer::import(input, which = region_gr)
    gene_data <- process_gene_data(
      gene_gr,
      gene_id = gene_id,
      gene_name = gene_name
    )
  } else if (methods::is(input, "TxDb")) {
    # TxDb object
    if (!requireNamespace("GenomicFeatures", quietly = TRUE)) {
      stop(
        "Package 'GenomicFeatures' is required for TxDb support. Install it with: BiocManager::install('GenomicFeatures')"
      )
    }
    gene_data <- extract_txdb_data(input, region_gr)
  } else if (is.data.frame(input)) {
    # Data frame - filter to region
    gene_data <- input[
      input$seqnames == as.character(GenomicRanges::seqnames(region_gr)) &
        input$end >= region_limits[1] &
        input$start <= region_limits[2],
    ]
  } else {
    stop("input must be a file path, TxDb object, or data frame")
  }

  # Ensure gene body rows exist for each gene

  # This is critical for partial gene overlaps where the intron line must be drawn
  gene_rows <- gene_data[gene_data$type == "gene", ]
  exon_rows <- gene_data[gene_data$type == "exon", ]

  if (nrow(gene_rows) == 0 && nrow(exon_rows) > 0) {
    # No gene body rows - create them from exon data
    # Group by gene_id (or other identifier) and create gene body spanning all exons
    gene_bodies <- do.call(
      rbind,
      lapply(
        split(exon_rows, exon_rows[[gene_id]]),
        function(g) {
          # Get the min/max coordinates across all exons
          x_min <- min(g$start, na.rm = TRUE)
          x_max <- max(g$end, na.rm = TRUE)

          # Create a single gene body row
          body_row <- g[1, , drop = FALSE]
          body_row$type <- "gene"
          body_row$start <- x_min
          body_row$end <- x_max
          body_row
        }
      )
    )
    gene_data <- rbind(gene_bodies, gene_data)
  }

  # Clip gene body coordinates to region limits
  # This ensures intron lines span the visible region for partial overlaps
  gene_idx <- gene_data$type == "gene"
  if (any(gene_idx)) {
    gene_data$start[gene_idx] <- pmax(
      gene_data$start[gene_idx],
      region_limits[1]
    )
    gene_data$end[gene_idx] <- pmin(gene_data$end[gene_idx], region_limits[2])
  }

  # Create the plot with strand-based coloring or uniform colors

  if (use_strand_colors) {
    # Ensure strand column exists and has proper factor levels

    if (!"strand" %in% names(gene_data)) {
      stop(
        "Strand column not found in gene data. Cannot use strand-based coloring."
      )
    }

    # Add resolved color values directly to the data based on strand
    # This bypasses ggplot2's scale system for more reliable color application
    gene_data$strand_color <- ifelse(
      gene_data$strand == "+",
      strand_colors["+"],
      ifelse(
        gene_data$strand == "-",
        strand_colors["-"],
        strand_colors["Unknown"]
      )
    )

    # Build geom_gene call with aes mapping to the resolved color values

    geom_call <- do.call(
      geom_gene,
      c(
        list(
          mapping = ggplot2::aes(
            start = start,
            end = end,
            y = .data[[y]],
            type = type,
            colour = strand_color,
            fill = strand_color
          ),
          exon_height = exon_height,
          intron_width = intron_width,
          clip_to_region = region_limits
        ),
        dots
      )
    )

    p <- ggplot2::ggplot(gene_data) +
      geom_call +
      ggplot2::scale_colour_identity() +
      ggplot2::scale_fill_identity()
  } else {
    # Uniform coloring

    geom_call <- do.call(
      geom_gene,
      c(
        list(
          mapping = ggplot2::aes(
            start = start,
            end = end,
            y = .data[[y]],
            type = type
          ),
          exon_height = exon_height,
          intron_width = intron_width,
          exon_color = exon_color,
          exon_fill = exon_fill,
          intron_color = intron_color,
          clip_to_region = region_limits
        ),
        dots
      )
    )

    p <- ggplot2::ggplot(gene_data) +
      geom_call
  }

  # Add labels if requested
  if (label_style != "none" && !is.null(label) && label %in% names(gene_data)) {
    # Get unique gene positions for labels (use gene type, not exons)
    label_data <- gene_data[gene_data$type == "gene", ]

    # Remove duplicates based on gene_id to avoid duplicate labels
    if (gene_id %in% names(label_data)) {
      label_data <- label_data[!duplicated(label_data[[gene_id]]), ]
    }

    # Apply max_labels filtering if specified
    if (!is.null(max_labels)) {
      label_data <- filter_labels(
        label_data,
        max_labels = max_labels,
        label_priority = label_priority,
        start_col = "start",
        end_col = "end"
      )
    }

    # Only proceed if there are labels to show
    if (nrow(label_data) > 0) {
      # Calculate label positions (middle of gene)
      label_data$label_x <- (label_data$start + label_data$end) / 2

      # Convert y to numeric for offset calculation
      # The discrete y-axis maps factor levels to integers (1, 2, 3, ...)
      y_values <- label_data[[y]]
      if (is.factor(y_values)) {
        label_data$label_y_num <- as.numeric(y_values)
      } else {
        label_data$label_y_num <- as.numeric(as.factor(y_values))
      }

      # Position label directly above the top of exons
      # Exons are centered on the y-position with height = exon_height * track_height
      # Track height is 1 in the coordinate system, so top of exon is at y + exon_height/2
      # Add offset to create spacing between exon and label
      # Larger offset (0.2) for better clearance and readability
      label_data$label_y <- label_data$label_y_num + (exon_height / 2) + 0.2

      # Determine which geom to use based on label_style
      use_repel <- FALSE
      use_check_overlap <- FALSE

      if (label_style == "repel") {
        # Force ggrepel usage
        if (!requireNamespace("ggrepel", quietly = TRUE)) {
          stop(
            "Package 'ggrepel' is required for label_style = 'repel'. Install it with: install.packages('ggrepel')"
          )
        }
        use_repel <- TRUE
      } else if (label_style == "auto") {
        # Use ggrepel if available, otherwise check_overlap
        if (requireNamespace("ggrepel", quietly = TRUE)) {
          use_repel <- TRUE
        } else {
          use_check_overlap <- TRUE
        }
      } else if (label_style == "simple") {
        # No overlap handling
        use_check_overlap <- FALSE
      }

      # Prepare default repel arguments for horizontal-only repositioning
      # This keeps labels horizontally aligned (no vertical nudging)
      if (use_repel) {
        default_repel_args <- list(
          direction = "x", # Only allow horizontal movement
          segment.color = NA, # Hide connecting segments
          box.padding = 0.3,
          point.padding = 0.2
        )
        # Merge with user-provided repel_args (user args take precedence)
        repel_call_args_base <- c(default_repel_args, repel_args)
        # Remove duplicates, keeping last occurrence (user's values)
        repel_call_args_base <- repel_call_args_base[
          !duplicated(names(repel_call_args_base), fromLast = TRUE)
        ]
      }

      # Prepare base aesthetic mapping
      base_aes <- ggplot2::aes(
        x = .data$label_x,
        y = .data$label_y,
        label = .data[[label]]
      )

      # Handle label color: use explicit label_color, or strand-based colors, or exon_fill
      if (!is.null(label_color)) {
        # Explicit label color provided
        if (use_repel) {
          # Use ggrepel with horizontal-only repositioning
          repel_call_args <- c(
            list(
              data = label_data,
              mapping = base_aes,
              size = label_size,
              color = label_color,
              fontface = "italic"
            ),
            repel_call_args_base
          )
          p <- p + do.call(ggrepel::geom_text_repel, repel_call_args)
        } else {
          # Use geom_text
          p <- p +
            ggplot2::geom_text(
              data = label_data,
              mapping = base_aes,
              size = label_size,
              vjust = 0,
              color = label_color,
              fontface = "italic",
              check_overlap = use_check_overlap
            )
        }
      } else if (use_strand_colors) {
        # Use strand-based colors for labels
        label_data$strand_color <- ifelse(
          label_data$strand == "+",
          strand_colors["+"],
          ifelse(
            label_data$strand == "-",
            strand_colors["-"],
            strand_colors["Unknown"]
          )
        )

        # Add color to aesthetic mapping
        color_aes <- ggplot2::aes(
          x = .data$label_x,
          y = .data$label_y,
          label = .data[[label]],
          colour = strand_color
        )

        if (use_repel) {
          # Use ggrepel with color mapping and horizontal-only repositioning
          repel_call_args <- c(
            list(
              data = label_data,
              mapping = color_aes,
              size = label_size,
              fontface = "italic"
            ),
            repel_call_args_base
          )
          p <- p +
            do.call(ggrepel::geom_text_repel, repel_call_args) +
            ggplot2::scale_colour_identity()
        } else {
          # Use geom_text with color mapping
          p <- p +
            ggplot2::geom_text(
              data = label_data,
              mapping = color_aes,
              size = label_size,
              vjust = 0,
              fontface = "italic",
              check_overlap = use_check_overlap
            ) +
            ggplot2::scale_colour_identity()
        }
      } else {
        # Use exon_fill for label color
        if (use_repel) {
          # Use ggrepel with horizontal-only repositioning
          repel_call_args <- c(
            list(
              data = label_data,
              mapping = base_aes,
              size = label_size,
              color = exon_fill,
              fontface = "italic"
            ),
            repel_call_args_base
          )
          p <- p + do.call(ggrepel::geom_text_repel, repel_call_args)
        } else {
          # Use geom_text
          p <- p +
            ggplot2::geom_text(
              data = label_data,
              mapping = base_aes,
              size = label_size,
              vjust = 0,
              color = exon_fill,
              fontface = "italic",
              check_overlap = use_check_overlap
            )
        }
      }
    }
  }

  # Apply theme and scale
  if (use_strand_colors) {
    # Determine which strands are present and build color vector accordingly
    # Discrete y-axis orders factor levels alphabetically: "-" < "+" < "Unknown"
    present_strands <- unique(as.character(gene_data$strand))
    all_strand_order <- c("-", "+", "Unknown")
    present_ordered <- all_strand_order[all_strand_order %in% present_strands]
    axis_colors <- strand_colors[present_ordered]

    p <- p +
      ggplot2::scale_y_discrete(
        expand = c(0.1, 0.1),
        drop = TRUE # Drop unused levels so colors match
      ) +
      scale_x_genome_region(region) +
      ez_gene_theme() +
      ggplot2::theme(
        axis.text.y = ggplot2::element_text(
          size = 14,
          face = "bold",
          colour = axis_colors
        )
      ) +
      ggplot2::labs(x = if (label_chr) paste0("Chr", chr) else NULL)
  } else {
    p <- p +
      ggplot2::scale_y_discrete(
        expand = c(0.1, 0.1),
        drop = TRUE
      ) +
      scale_x_genome_region(region) +
      ez_gene_theme() +
      ggplot2::labs(x = if (label_chr) paste0("Chr", chr) else NULL)
  }

  if (border) p <- apply_border_theme(p)

  return(p)
}
