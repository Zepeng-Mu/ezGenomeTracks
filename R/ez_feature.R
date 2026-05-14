# ezGenomeTracks - Feature functions (split from ez_wrappers.R)
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
#'   Either `region` or `gene` (with `gene_db`) must be provided.
#'   Example: "chr1:1000000-2000000"
#' @param gene Gene name/symbol to look up (e.g., "PTPRC", "TP53").
#'   When provided, the region is automatically determined from the gene coordinates
#'   in `gene_db`. Either `region` or `gene` must be provided.
#' @param gene_db TxDb object for gene coordinate lookup when using `gene` parameter.
#' @param org_db Optional OrgDb object for gene symbol mapping. If NULL (default),
#'   auto-detects available OrgDb packages.
#' @param extend Numeric. Amount to extend the region beyond the gene body when
#'   using `gene` parameter. Default: 0.1 (10% of gene length on each side).
#' @param extend_type How to interpret `extend`: "proportion" (relative to gene
#'   length) or "bp" (absolute base pairs). Default: "proportion".
#' @param color Border color of the features. Default: "black"
#' @param fill Fill color of the features. When `use_score = TRUE`, this will be
#'   used as the high value in the color gradient. Default: "gray70"
#' @param alpha Transparency level of the features (0 = transparent, 1 = opaque).
#'   Default: 0.7
#' @param height Height of the feature rectangles (0 to 1). Default: 0.8
#' @param use_score Logical indicating whether to use the 'score' column for
#'   fill color. If TRUE, a gradient from white to the specified fill color will
#'   be used. Default: FALSE
#' @param border Logical. If TRUE, adds a black border around the plot panel.
#'   Default: FALSE
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
#' @importFrom ggplot2 ggplot aes scale_fill_gradient scale_fill_identity ylim
#' @importFrom methods is
#'
#' @examples
#' # From a data frame with uniform coloring
#' features <- data.frame(
#'   seqnames = c("chr1", "chr1", "chr1"),
#'   start = c(1000, 3000, 5000),
#'   end = c(2000, 4000, 6000),
#'   name = c("peak1", "peak2", "peak3"),
#'   score = c(10, 30, 50)
#' )
#' track <- ez_feature(
#'   features,
#'   "chr1:1-10000",
#'   fill = "darkgreen",
#'   alpha = 0.8
#' )
ez_feature <- function(
  input,
  region = NULL,
  gene = NULL,
  gene_db = NULL,
  org_db = NULL,
  extend = 0.1,
  extend_type = c("proportion", "bp"),
  color = "black",
  fill = "gray70",
  alpha = 0.7,
  height = 0.8,
  use_score = FALSE,
  border = FALSE,
  ...
) {
  # Resolve region from either region string or gene name
  region <- .resolve_region(
    region = region,
    gene = gene,
    gene_db = gene_db,
    org_db = org_db,
    extend = extend,
    extend_type = extend_type
  )

  # Check if data is a file path or data frame
  if (is.character(input) && length(input) == 1) {
    # It's a file path, import it as a data frame
    region_gr <- parse_region(region)
    input <- import_genomic_data(input, which = region_gr)
  }

  if (is.data.frame(input)) {
    # Filter data to the specified region
    region_gr <- parse_region(region)
    input <- input[
      input$seqnames == as.character(GenomicRanges::seqnames(region_gr)) &
      input$start < GenomicRanges::end(region_gr) &
      input$end > GenomicRanges::start(region_gr),
    ]

    # Note: geom_feature's default_aes already maps start->xmin, end->xmax
    if (use_score && "score" %in% colnames(input)) {
      p <- ggplot2::ggplot(input) +
        geom_feature(
          ggplot2::aes(fill = score),
          color = color,
          alpha = alpha,
          height = height,
          ...
        ) +
        ggplot2::scale_fill_gradient(low = "white", high = fill)
    } else {
      p <- ggplot2::ggplot(input) +
        geom_feature(
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
      ez_feature_theme()

    if (border) p <- apply_border_theme(p)

    return(p)
  } else {
    stop("input must be a file path or data frame")
  }
}
