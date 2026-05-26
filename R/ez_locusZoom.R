# ezGenomeTracks - LocusZoom-style regional association plots

#' Easy LocusZoom-style regional association plot
#'
#' This function creates a regional association plot (LocusZoom-style) for GWAS
#' results within a specific genomic region. It supports LD-based coloring,
#' lead SNP highlighting, and is designed to stack with other track types via
#' `vstack_plot()`.
#'
#' @description
#' Creates a regional Manhattan plot focused on a single genomic locus, commonly
#' used for fine-mapping and visualization of association signals with linkage
#' disequilibrium (LD) information. The plot uses genomic coordinate formatting
#' consistent with `ez_coverage()` and `ez_gene()`, making it suitable for
#' multi-track visualizations.
#'
#' @param input A data frame containing GWAS results with columns for chromosome,
#'   position, p-values, and optionally SNP names. Supports both GWAS-style
#'   (CHR, BP, P) and GRanges-style (seqnames, start, pvalue) column naming.
#' @param region Genomic region string (e.g., "chr1:1000000-2000000"). Data is
#'   filtered to this region and the x-axis uses coordinate-based formatting.
#'   Either `region` or `gene` must be provided.
#' @param gene Gene name/symbol to look up (e.g., "PTPRC", "TP53"). When provided,
#'   the region is automatically determined from gene coordinates in `gene_db`.
#'   Either `region` or `gene` must be provided.
#' @param gene_db TxDb object for gene coordinate lookup when using `gene` parameter.
#'   Required if `gene` is provided.
#' @param org_db Optional OrgDb object for gene symbol mapping. If NULL (default),
#'   auto-detects available OrgDb packages.
#' @param extend Numeric. Amount to extend the region beyond the gene body when
#'   using `gene` parameter. Default: 0.1 (10% of gene length on each side).
#' @param extend_type How to interpret `extend`: "proportion" (relative to gene
#'   length) or "bp" (absolute base pairs). Default: "proportion".
#' @param chr Character string specifying the column name for chromosome numbers.
#'   Default: auto-detect from "CHR", "seqnames", "chrom", etc.
#' @param bp Character string specifying the column name for base pair positions.
#'   Default: auto-detect from "BP", "start", "pos", etc.
#' @param p Character string specifying the column name for p-values.
#'   Default: auto-detect from "P", "pvalue", "p.value", etc.
#' @param snp Character string specifying the column name for SNP identifiers.
#'   Default: auto-detect from "SNP", "rsid", "variant_id", etc.
#' @param logp Logical indicating whether to plot -log10(p-values).
#'   Default: TRUE.
#' @param size Numeric value for point size in the plot. Default: 1.
#' @param color Default point color when `r2` is not provided. Default: "grey50".
#' @param lead_snp Character string or vector of SNP IDs to highlight as the
#'   lead variant(s). Highlighted with `highlight_color`. Default: NULL.
#' @param r2 Numeric vector of r² values for coloring points by linkage
#'   disequilibrium with lead variant. Must be same length as number of rows
#'   in data. When provided, points are colored using a gradient from blue
#'   (low LD) to red (high LD). Default: NULL.
#' @param colors Vector of colors for the r² gradient. Default: LocusZoom palette
#'   `c("blue3", "skyblue", "green2", "orange", "red3")`.
#' @param highlight_snps Character vector of SNP IDs to highlight, or a data frame
#'   with chr, bp, p columns. Default: NULL.
#' @param highlight_color Color for highlighting lead or specified SNPs.
#'   Default: "purple".
#' @param threshold_p Numeric p-value threshold for drawing a significance line.
#'   If NULL, no line is drawn. Default: NULL.
#' @param threshold_color Color for the significance threshold line.
#'   Default: "red".
#' @param threshold_linetype Linetype for the significance threshold line.
#'   Default: 2 (dashed).
#' @param y_axis_style Y-axis style: "none", "simple", or "full".
#'   Default: "none" (suitable for stacking).
#' @param y_axis_label Label for the y-axis.
#'   Default: `expression(paste("-log"[10], "(P)"))`.
#' @param color_by How points should be colored. Can be "r2" (use the `r2`
#'   argument for LD coloring), "none" (use a single `color`), or a column name
#'   in the data for discrete/continuous coloring. Default: "r2" if `r2` is
#'   provided, otherwise "none".
#' @param border Logical. If TRUE, adds a black border around the plot panel. Default: FALSE
#' @param ... Additional arguments passed to `geom_manhattan()`.
#'
#' @return A ggplot2 object containing the regional association plot.
#'
#' @details
#' This function is designed for visualizing association results at a single

#' genomic locus, similar to the LocusZoom web tool. Key features:
#'
#' - **LD coloring**: When `r2` is provided, points are colored by linkage
#'   disequilibrium with the lead variant, using the classic LocusZoom color
#'   scheme (blue → red gradient).
#'
#' - **Gene-based regions**: Use `gene` parameter to automatically look up
#'   gene coordinates and define the viewing region.
#'
#' - **Stackable**: Uses `scale_x_genome_region()` for x-axis formatting,
#'   allowing seamless stacking with `ez_coverage()`, `ez_gene()`, and other
#'   tracks via `vstack_plot()`.
#'
#' For genome-wide Manhattan plots across multiple chromosomes, use
#' `ez_manhattan()` instead.
#'
#' @seealso
#' \code{\link{ez_manhattan}} for genome-wide Manhattan plots,
#' \code{\link{geom_manhattan}} for the underlying geom,
#' \code{\link{vstack_plot}} for combining with other tracks
#'
#' @export
#' @importFrom ggplot2 ggplot
#' @importFrom dplyr filter
#' @importFrom rlang .data
#'
#' @examples
#' # Create example data for a region
#' set.seed(42)
#' region_data <- data.frame(
#'   CHR = rep(3, 100),
#'   BP = seq(1000, 100000, length.out = 100),
#'   P = c(runif(95, 0.01, 1), runif(5, 1e-8, 1e-4)),
#'   SNP = paste0("rs", 1:100)
#' )
#'
#' # Basic regional plot
#' ez_locusZoom(region_data, region = "chr3:1000-100000")
#'
#' # With LD coloring (simulated r2 values)
#' # In practice, r2 would come from LD calculations with lead SNP
#' r2_values <- runif(100, 0, 1)
#' ez_locusZoom(
#'   region_data,
#'   region = "chr3:1000-100000",
#'   r2 = r2_values,
#'   lead_snp = "rs96",
#'   size = 2
#' )
#'
#' \dontrun{
#' # Using gene name to define region (requires TxDb)
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' ez_locusZoom(
#'   gwas_results,
#'   gene = "TP53",
#'   gene_db = TxDb.Hsapiens.UCSC.hg38.knownGene,
#'   r2 = ld_values
#' )
#'
#' # Stack with gene track
#' p1 <- ez_locusZoom(gwas_results, region = "chr17:7500000-7700000", r2 = ld)
#' p2 <- ez_gene(txdb, region = "chr17:7500000-7700000")
#' vstack_plot(p1, p2, heights = c(2, 1))
#' }
ez_locusZoom <- function(
    input,
    region = NULL,
    gene = NULL,
    gene_db = NULL,
    org_db = NULL,
    extend = 0.1,
    extend_type = c("proportion", "bp"),
    chr = NULL,
    bp = NULL,
    p = NULL,
    snp = NULL,
    logp = TRUE,
    size = 1,
    color = "grey50",
    lead_snp = NULL,
    r2 = NULL,
    colors = NULL,
    highlight_snps = NULL,
    highlight_color = "purple",
    threshold_p = NULL,
    threshold_color = "red",
    threshold_linetype = 2,
    y_axis_style = c("none", "simple", "full"),
    y_axis_label = expression(paste("-log"[10], "(P)")),
    color_by = NULL,
    border = FALSE,
    ...) {
  # Resolve region from either region string or gene name
  region <- .resolve_region(
    region = region, gene = gene, gene_db = gene_db,
    org_db = org_db, extend = extend, extend_type = extend_type
  )

  # Validate inputs
  y_axis_style <- match.arg(y_axis_style)

  # Validate input is a data frame
  if (!is.data.frame(input)) {
    stop("input must be a data frame containing GWAS results")
  }

  plotDt <- input

  # Auto-detect column names if not provided
  cols <- auto_detect_gwas_columns(plotDt, chr, bp, p, snp)
  chr <- cols$chr
  bp <- cols$bp
  p <- cols$p
  snp <- cols$snp

  # Filter data to region
  region_gr <- parse_region(region)
  region_chr <- as.character(GenomicRanges::seqnames(region_gr))
  region_start <- GenomicRanges::start(region_gr)
  region_end <- GenomicRanges::end(region_gr)

  # Normalize chromosome names for comparison (handle chr prefix)
  region_chr_normalized <- gsub("^chr", "", region_chr, ignore.case = TRUE)

  # Store original row indices for r2 subsetting
  plotDt$.original_row <- seq_len(nrow(plotDt))

  plotDt <- plotDt |>
    dplyr::filter(
      gsub("^chr", "", as.character(.data[[chr]]), ignore.case = TRUE) ==
        region_chr_normalized,
      .data[[bp]] >= region_start,
      .data[[bp]] <= region_end
    )

  if (nrow(plotDt) == 0) {
    warning("No data points found in the specified region.")
  }

  # Subset r2 values to match filtered data
  if (!is.null(r2)) {
    if (length(r2) != nrow(input)) {
      stop("Length of r2 vector must match the number of rows in input data.")
    }
    r2 <- r2[plotDt$.original_row]
  }

  # Remove helper column
  plotDt$.original_row <- NULL

  # Determine color_by based on r2 presence (user override takes priority)
  color_by <- if (!is.null(color_by)) color_by else if (!is.null(r2)) "r2" else "none"

  # Create plot using geom_manhattan in regional mode
  plot_obj <- ggplot2::ggplot(plotDt) +
    geom_manhattan(
      data = plotDt,
      region = region,
      mode = "regional",
      chr = chr,
      bp = bp,
      p = p,
      snp = snp,
      logp = logp,
      size = size,
      color = color,
      lead_snp = lead_snp,
      r2 = r2,
      colors = colors,
      highlight_snps = highlight_snps,
      highlight_color = highlight_color,
      threshold_p = threshold_p,
      threshold_color = threshold_color,
      threshold_linetype = threshold_linetype,
      color_by = color_by,
      y_axis_label = y_axis_label,
      ...
    )

  # Apply regional theme
  plot_obj <- plot_obj + ez_manhattan_theme(y_axis_style = y_axis_style)

  if (border) plot_obj <- apply_border_theme(plot_obj)

  return(plot_obj)
}
