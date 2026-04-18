# ezGenomeTracks - Hic functions (split from ez_wrappers.R)
#' Easy Hi-C track visualization
#'
#' This function creates a Hi-C contact matrix visualization from various input types.
#' It provides a high-level interface supporting both square heatmap and triangle
#' (rotated) views commonly used in genome browsers.
#'
#' @param data Input data. Can be:
#'   \itemize{
#'     \item A matrix: Dense contact matrix
#'     \item A data frame: Sparse format with (pos1, pos2, score) or (bin1, bin2, score)
#'     \item A file path: Tab-delimited matrix file
#'   }
#' @param region Genomic region to display (e.g., "chr1:1000000-2000000").
#'   Either `region` or `gene` (with `gene_db`) must be provided.
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
#' @param resolution Resolution of the Hi-C data in base pairs (default: 10000)
#' @param style Visualization style: "triangle" (default, rotated view) or "square"
#' @param palette Color palette: "cooler" (red, default), "ylgnbu", "viridis", or "bwr"
#' @param transform Scale transformation: "identity" (linear), "log10" (default), "log2", "sqrt"
#' @param limits Numeric vector of length 2 for color scale limits (default: NULL, auto)
#' @param max_distance Maximum interaction distance to show in base pairs (default: NULL, show all).
#'   Only applies to triangle style.
#' @param rasterize Logical. If TRUE and ggrastr package is available, rasterize the
#'   Hi-C layer for better performance with large matrices. Default: FALSE
#' @param rasterize_dpi Resolution for rasterization in dots per inch. Default: 300
#' @param show_diagonal Logical. If TRUE, show the diagonal (self-interactions). Default: TRUE
#' @param ... Additional arguments passed to geom_hic or geom_hic_triangle
#'
#' @return A ggplot2 object
#' @export
#' @importFrom ggplot2 ggplot aes coord_fixed coord_cartesian labs scale_y_continuous
#' @examples
#' \dontrun{
#' # Create example data
#' mat <- matrix(runif(2500), nrow = 50)
#' mat <- mat + t(mat)  # Make symmetric
#' diag(mat) <- diag(mat) * 2  # Stronger diagonal
#'
#' # Triangle view (default)
#' ez_hic(mat, "chr1:1000000-1500000", resolution = 10000)
#'
#' # Square heatmap view
#' ez_hic(mat, "chr1:1000000-1500000", resolution = 10000, style = "square")
#'
#' # With log10 transformation and custom palette
#' ez_hic(mat, "chr1:1000000-1500000",
#'   resolution = 10000,
#'   trans = "log10",
#'   palette = "ylgnbu"
#' )
#'
#' # Limit maximum distance shown
#' ez_hic(mat, "chr1:1000000-1500000",
#'   resolution = 10000,
#'   max_distance = 200000
#' )
#'
#' # With rasterization for large regions
#' ez_hic(mat, "chr1:1000000-1500000",
#'   resolution = 10000,
#'   rasterize = TRUE,
#'   rasterize_dpi = 300
#' )
#'
#' # Using gene name for region lookup
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' ez_hic(hic_data, gene = "PTPRC", gene_db = TxDb.Hsapiens.UCSC.hg38.knownGene)
#' }
ez_hic <- function(
  data,
  region = NULL,
  gene = NULL,
  gene_db = NULL,
  org_db = NULL,
  extend = 0.1,
  extend_type = c("proportion", "bp"),
  resolution = 10000,
  style = c("triangle", "square"),
  palette = c("cooler", "ylgnbu", "viridis", "bwr"),
  transform = "identity",
  limits = NULL,
  max_distance = NULL,
  rasterize = FALSE,
  rasterize_dpi = 300,
  show_diagonal = TRUE,
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

  style <- match.arg(style)
  palette <- match.arg(palette)

  # Process the input data
  upper_triangle <- (style == "triangle")
  hic_df <- process_hic_data(
    data = data,
    region = region,
    resolution = resolution,
    upper_triangle = upper_triangle
  )

  # Filter out diagonal if requested
  if (!show_diagonal && style == "triangle") {
    hic_df <- hic_df[hic_df$pos1 != hic_df$pos2, ]
  }

  # Filter by max_distance for triangle view
  if (!is.null(max_distance) && style == "triangle") {
    hic_df <- hic_df[(hic_df$pos2 - hic_df$pos1) <= max_distance, ]
  }

  # Parse region for x-axis limits
  region_gr <- parse_region(region)
  chr <- as.character(GenomicRanges::seqnames(region_gr))
  start_pos <- GenomicRanges::start(region_gr)
  end_pos <- GenomicRanges::end(region_gr)

  # Create the plot based on style
  if (style == "triangle") {
    # Triangle view: x = midpoint, y = distance
    p <- ggplot2::ggplot(
      hic_df,
      ggplot2::aes(x = .data$pos1, y = .data$pos2, fill = .data$score)
    ) +
      geom_hic_triangle(
        resolution = resolution,
        max_distance = max_distance,
        rasterize = rasterize,
        rasterize_dpi = rasterize_dpi,
        ...
      )
  } else {
    # Square heatmap view
    p <- ggplot2::ggplot(
      hic_df,
      ggplot2::aes(x = .data$pos1, y = .data$pos2, fill = .data$score)
    ) +
      geom_hic(
        resolution = resolution,
        rasterize = rasterize,
        rasterize_dpi = rasterize_dpi,
        ...
      )
  }

  # Add color scale
  p <- p + scale_fill_hic(palette = palette, trans = transform, limits = limits)

  # Add appropriate scales and theme
  if (style == "triangle") {
    # For triangle, x-axis is genomic position, y-axis is distance
    # Triangle points UPWARD with diagonal at y=0 (bottom) and apex at top
    # Calculate y-axis limit based on max_distance or region span
    y_max <- if (!is.null(max_distance)) max_distance else (end_pos - start_pos)

    p <- p +
      scale_x_genome_region(region, oob = scales::oob_keep) +
      ggplot2::scale_y_continuous(
        expand = c(0, 0),
        limits = c(0, y_max), # Positive range, apex at top
        labels = function(x) format_genomic_coord(x),
        oob = scales::oob_keep
      ) +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::labs(x = paste0("Chr", chr), y = "", fill = "Score") +
      ez_hic_theme()
  } else {
    # Square view: both axes are genomic positions
    p <- p +
      scale_x_genome_region(region, oob = scales::oob_keep) +
      ggplot2::scale_y_continuous(
        expand = c(0, 0),
        limits = c(start_pos, end_pos),
        labels = function(x) format_genomic_coord(x),
        oob = scales::oob_keep
      ) +
      ggplot2::coord_fixed(ratio = 1) +
      ggplot2::labs(
        x = paste0("Chr", chr),
        y = paste0("Chr", chr),
        fill = "Score"
      ) +
      ez_hic_theme()
  }

  return(p)
}
