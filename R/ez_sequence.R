#' Easy nucleotide sequence track visualization
#'
#' Creates a nucleotide sequence track for a given genomic region, fetching the
#' DNA sequence from a BSgenome object. Each nucleotide is rendered as a colored
#' tile following UCSC Genome Browser conventions (A = green, C = blue,
#' G = gold, T = red, N = grey). Single-letter labels are drawn on the tiles
#' when the region is narrow enough to display them legibly.
#'
#' @param genome A BSgenome object providing the reference genome sequence
#'   (e.g., `BSgenome.Hsapiens.UCSC.hg38::Hsapiens`). The corresponding
#'   BSgenome data package must be installed.
#' @param region Genomic region to display, as a character string in the format
#'   `"chr:start-end"` (e.g., `"chr1:1000000-1000050"`). Either `region` or
#'   `gene` (with `gene_db`) must be provided.
#' @param gene Gene name/symbol to look up (e.g., `"PTPRC"`, `"TP53"`). When
#'   provided, the region is determined automatically from the gene coordinates
#'   in `gene_db`. Either `region` or `gene` must be provided.
#' @param gene_db TxDb object for gene coordinate lookup when using `gene`.
#' @param org_db Optional OrgDb object for gene symbol mapping. When `NULL`
#'   (default), the function auto-detects available OrgDb packages.
#' @param extend Numeric. Amount to extend the region beyond the gene body when
#'   using `gene`. Default: `0.1` (10% of gene length on each side).
#' @param extend_type How to interpret `extend`: `"proportion"` (relative to
#'   gene length) or `"bp"` (absolute base pairs). Default: `"proportion"`.
#' @param colors Optional named character vector of color overrides for
#'   individual nucleotides. Only the specified bases are overridden; others
#'   retain UCSC defaults. Example: `c(A = "purple", T = "pink")`.
#' @param style Character. Visual style for the sequence track:
#'   - `"text"` (default): bold, colored nucleotide letters with no background
#'     or border. Letter color is derived from `colors` / UCSC defaults.
#'   - `"tile"`: colored background tiles with optional letter labels on top.
#' @param show_labels Logical or `"auto"`. Whether to draw nucleotide letters.
#'   When `"auto"` (default), labels are shown only when the region width is ≤
#'   `max_label_bp`. For `style = "text"`, hiding labels (`FALSE`) renders an
#'   empty panel; consider switching to `style = "tile"` for wide regions.
#' @param max_label_bp Integer. Maximum region width (in base pairs) at which
#'   `show_labels = "auto"` will display nucleotide letters. Default: `200`.
#' @param label_size Numeric. Font size for nucleotide letters (in pt).
#'   Default: `3`.
#' @param label_color Character. Color of nucleotide letters when
#'   `style = "tile"`. Ignored for `style = "text"` (color comes from the
#'   nucleotide palette). Default: `"white"`.
#' @param tile_height Numeric (0–1). Height of each nucleotide tile as a
#'   proportion of the panel height. Only used when `style = "tile"`.
#'   Default: `0.8`.
#' @param ... Additional arguments passed to `geom_sequence()`.
#'
#' @return A ggplot2 object representing the sequence track, compatible with
#'   `vstack_plot()` for stacking with other tracks.
#'
#' @details
#' The BSgenome data package for your organism of interest must be installed
#' separately. For example:
#' ```r
#' BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' ez_sequence(Hsapiens, region = "chr1:1000000-1000050")
#' ```
#'
#' UCSC nucleotide color defaults can be inspected via `ez_sequence_palette()`.
#'
#' @export
#' @importFrom ggplot2 ggplot aes scale_fill_identity ylim coord_cartesian
#' @importFrom methods is
#'
#' @examples
#' \dontrun{
#' library(BSgenome.Hsapiens.UCSC.hg38)
#'
#' # Default style: bold colored letters, no background
#' ez_sequence(Hsapiens, region = "chr1:1000000-1000030")
#'
#' # Tile style: colored background tiles with white labels
#' ez_sequence(Hsapiens, region = "chr1:1000000-1000030", style = "tile")
#'
#' # Force no labels (renders empty panel for text style)
#' ez_sequence(Hsapiens, region = "chr1:1000000-1000030", show_labels = FALSE)
#'
#' # Override colors
#' ez_sequence(
#'   Hsapiens,
#'   region = "chr1:1000000-1000030",
#'   colors = c(A = "purple", T = "pink")
#' )
#'
#' # Stack with a coverage track
#' library(aplot)
#' cov <- ez_coverage(my_bw, region = "chr1:1000000-1000100")
#' seq <- ez_sequence(Hsapiens, region = "chr1:1000000-1000100")
#' vstack_plot(list(cov, seq), heights = c(3, 1))
#' }
ez_sequence <- function(
  genome,
  region       = NULL,
  gene         = NULL,
  gene_db      = NULL,
  org_db       = NULL,
  extend       = 0.1,
  extend_type  = c("proportion", "bp"),
  colors       = NULL,
  style        = c("text", "tile"),
  show_labels  = "auto",
  max_label_bp = 200L,
  label_size   = 3,
  label_color  = "white",
  tile_height  = 0.8,
  ...
) {
  style <- match.arg(style)
  extend_type <- match.arg(extend_type)

  # Resolve region string (from region or gene + gene_db)
  region <- .resolve_region(
    region      = region,
    gene        = gene,
    gene_db     = gene_db,
    org_db      = org_db,
    extend      = extend,
    extend_type = extend_type
  )

  region_gr <- parse_region(region)

  # Build color palette (UCSC defaults + user overrides)
  palette <- ez_sequence_palette(colors = colors)

  # Fetch sequence and build per-nucleotide data frame
  df <- process_sequence_data(genome, region_gr, palette)

  # Resolve show_labels
  if (identical(show_labels, "auto")) {
    region_width <- GenomicRanges::end(region_gr) - GenomicRanges::start(region_gr) + 1L
    show_labels <- region_width <= max_label_bp
  }

  # Build the ggplot
  p <- ggplot2::ggplot(df, ggplot2::aes(
    x     = .data$position,
    label = .data$nucleotide,
    fill  = .data$fill
  )) +
    geom_sequence(
      style       = style,
      show_labels = show_labels,
      label_size  = label_size,
      label_color = label_color,
      tile_height = tile_height,
      ...
    ) +
    ggplot2::scale_fill_identity() +
    ggplot2::ylim(0, 1) +
    scale_x_genome_region(region) +
    ggplot2::coord_cartesian(clip = "off") +
    ez_sequence_theme()

  p
}
