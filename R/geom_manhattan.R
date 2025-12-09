#' @title Create a Manhattan plot
#' @description `geom_manhattan` creates a Manhattan plot for visualizing GWAS or QTL data.
#'   The function automatically detects whether to use regional mode (single chromosome)
#'   or genome-wide mode (multiple chromosomes) based on the data.
#'
#'   In **regional mode** (single chromosome): Uses `scale_x_genome_region()` for x-axis
#'   formatting consistent with other track functions like `ez_coverage` and `ez_gene`.
#'   This mode is suitable for LocusZoom-style plots and can be stacked with other tracks.
#'
#'   In **genome-wide mode** (multiple chromosomes): Uses cumulative base pair positions
#'   with chromosome labels on the x-axis and alternating colors per chromosome.
#'
#' @param mapping Set of aesthetic mappings created by `aes()`. If not specified, the default mappings are used.
#' @param data A `data.frame` containing the data to be plotted. Must include columns for chromosome, base pair position, and p-value. Optional column for SNP identifier.
#' @param region Optional genomic region string (e.g., "chr1:1000000-2000000") to force regional mode
#'   and set x-axis limits. When provided, data is NOT filtered (use `ez_manhattan` for filtering).
#' @param mode Plot mode: "auto" (default, detect from data), "regional", or "genome_wide".
#' @param stat The statistical transformation to apply to the data (default: "identity").
#' @param position Position adjustment, either as a string naming a position adjustment function, or the result of a call to a position adjustment function.
#' @param na.rm If `FALSE` (default), missing values are removed with a warning. If `TRUE`, missing values are silently removed.
#' @param show.legend Logical. Should this layer be displayed in the legend? `NA` for automatic, `TRUE` always, `FALSE` never.
#' @param inherit.aes If `FALSE`, overrides the default aesthetics; if `TRUE`, inherits them.
#' @param chr Name of the chromosome column in `data`. Supports both GWAS-style ("CHR") and
#'   GRanges-style ("seqnames") conventions. Default: auto-detect.
#' @param bp Name of the base pair position column in `data`. Supports both GWAS-style ("BP")
#'   and GRanges-style ("start") conventions. Default: auto-detect.
#' @param p Name of the p-value column in `data`. Supports "P", "pvalue", "p.value", etc.
#'   Default: auto-detect.
#' @param snp Name of the SNP identifier column in `data` (default: "SNP" or "snp").
#' @param logp Logical. If TRUE (default), -log10() transformation is applied to p-values.
#' @param size Point size (default: 0.5).
#' @param color Default point color for regional mode when colorBy is not "r2" (default: "grey50").
#' @param lead_snp Vector of SNP IDs to highlight (also accepts lead.snp for backward compatibility).
#' @param r2 Vector of R-squared values for linkage disequilibrium (LD) coloring. Should be in same order as data rows.
#' @param colors Vector of colors to use for alternating chromosome colors in genome-wide mode (default: c("grey", "skyblue")).
#' @param highlight_snps Data frame of SNPs to highlight, with columns matching chr, bp, and p.
#' @param highlight_color Color for highlighted SNPs (default: "purple").
#' @param highlight_shape Shape for highlighted SNPs (default: 18).
#' @param threshold_p A numeric value for the p-value threshold to draw a horizontal line (e.g., 5e-8).
#' @param threshold_color Color for the threshold line (default: "red").
#' @param threshold_linetype Linetype for the threshold line (default: 2).
#' @param colorBy Character string indicating how points should be colored. Options:
#'   - "chr": Alternating chromosome colors (genome-wide mode only, ignored in regional mode)
#'   - "r2": Continuous color based on R-squared values (requires r2 parameter)
#'   - "none": Single color specified by `color` parameter
#' @param x_axis_label X-axis label (default: NULL, auto-generated).
#' @param y_axis_label Label for the y-axis (default: expression for -log10(P)).
#' @param ... Additional arguments passed to `ggplot2::geom_point()`.
#' @return A list of `ggplot2` layers and scales.
#' @export
#' @importFrom ggplot2 layer ggproto aes GeomPoint
#' @importFrom dplyr mutate filter select arrange group_by summarise inner_join
#' @importFrom rlang .data
geom_manhattan <- function(
    mapping = NULL, data = NULL, region = NULL,
    mode = c("auto", "regional", "genome_wide"),
    stat = "identity", position = "identity",
    na.rm = FALSE, show.legend = NA, inherit.aes = TRUE,
    chr = NULL, bp = NULL, p = NULL, snp = NULL,
    logp = TRUE, size = 0.5, color = "grey50",
    lead_snp = NULL, r2 = NULL, colors = c("grey", "skyblue"),
    highlight_snps = NULL, highlight_color = "purple", highlight_shape = 18,
    threshold_p = NULL, threshold_color = "red", threshold_linetype = 2,
    colorBy = c("auto", "chr", "r2", "none"),
    x_axis_label = NULL, y_axis_label = NULL, ...
) {
  mode <- match.arg(mode)
  colorBy <- match.arg(colorBy)

  # Validate input data
  if (is.null(data)) stop("Data cannot be NULL.")


  # Auto-detect column names with support for both GWAS and GRanges conventions
  detect_column <- function(data, candidates, param_name) {
    for (col in candidates) {
      if (col %in% colnames(data)) return(col)
    }
    stop(paste0("Could not find ", param_name, " column. Expected one of: ",
                paste(candidates, collapse = ", ")))
  }

  # Set column names with auto-detection
  if (is.null(chr)) {
    chr <- detect_column(data, c("CHR", "chr", "seqnames", "chrom", "chromosome"), "chromosome")
  }
  if (is.null(bp)) {
    bp <- detect_column(data, c("BP", "bp", "start", "pos", "position", "POS"), "position")
  }
  if (is.null(p)) {
    p <- detect_column(data, c("P", "p", "pvalue", "p.value", "pval", "P.value"), "p-value")
  }
  if (is.null(snp)) {
    # SNP column is optional
    snp_candidates <- c("SNP", "snp", "rsid", "id", "variant_id", "marker")
    for (col in snp_candidates) {
      if (col %in% colnames(data)) {
        snp <- col
        break
      }
    }
  }

  # Validate required columns exist
  if (!chr %in% colnames(data)) {
    stop(paste0("Chromosome column '", chr, "' not found in data."))
  }
  if (!bp %in% colnames(data)) {
    stop(paste0("Position column '", bp, "' not found in data."))
  }
  if (!p %in% colnames(data)) {
    stop(paste0("P-value column '", p, "' not found in data."))
  }

  # Prepare data with standardized column names
  # First, add r2 values to original data BEFORE arranging to preserve correspondence
  if (!is.null(r2)) {
    if (length(r2) != nrow(data)) {
      stop("Length of r2 vector must match the number of rows in data.")
    }
    data$r2_value <- r2
  }

  # Build list of columns to select
  optional_cols <- "r2_value"
  if (!is.null(snp) && snp %in% colnames(data)) {
    optional_cols <- c(optional_cols, snp)
  }

  plot_data <- data |>
    dplyr::select(
      CHR = .data[[chr]],
      BP = .data[[bp]],
      P = .data[[p]],
      dplyr::any_of(optional_cols)
    )

  # Rename SNP column to standardized name if it exists and isn't already named SNP
  if (!is.null(snp) && snp %in% colnames(plot_data) && snp != "SNP") {
    names(plot_data)[names(plot_data) == snp] <- "SNP"
  }

  plot_data <- plot_data |> dplyr::arrange(.data$CHR, .data$BP)

  # Apply log transformation
  if (logp) {
    plot_data <- plot_data |> dplyr::mutate(logp = -log10(.data$P))
  } else {
    plot_data <- plot_data |> dplyr::mutate(logp = .data$P)
  }

  # Determine plot mode

  n_chr <- length(unique(plot_data$CHR))
  if (mode == "auto") {
    is_regional <- n_chr == 1 || !is.null(region)
  } else {
    is_regional <- mode == "regional"
  }

  # Handle colorBy for regional mode

  if (colorBy == "auto") {
    if (is_regional) {
      colorBy <- if (!is.null(r2)) "r2" else "none"
    } else {
      colorBy <- "chr"
    }
  }

  # Warn if chr coloring requested in regional mode

  if (is_regional && colorBy == "chr" && n_chr == 1) {
    message("colorBy='chr' is ignored in regional mode (single chromosome). Using 'none'.")
    colorBy <- "none"
  }

  # Process based on mode
  if (!is_regional) {
    # GENOME-WIDE MODE: Calculate cumulative positions
    # https://danielroelfs.com/posts/how-i-create-manhattan-plots-using-ggplot
    data_cum <- plot_data |>
      dplyr::group_by(.data$CHR) |>
      dplyr::summarise(max_BP = max(.data$BP)) |>
      dplyr::mutate(BP_add = dplyr::lag(cumsum(.data$max_BP), default = 0)) |>
      dplyr::select("CHR", "BP_add")

    plot_data <- plot_data |>
      dplyr::inner_join(data_cum, by = "CHR") |>
      dplyr::mutate(BP = .data$BP + .data$BP_add)

    # Calculate chromosome midpoints for x-axis labels
    axis_df <- plot_data |>
      dplyr::group_by(.data$CHR) |>
      dplyr::summarize(center = mean(.data$BP))

    # Assign colors based on chromosome
    plot_data$color_group <- factor(plot_data$CHR %% 2 + 1, levels = c(1, 2))
  }

  # Default mapping
 default_mapping <- aes(x = .data$BP, y = .data$logp)

  # Create mapping based on colorBy option
  if (!is.null(mapping) && inherits(mapping, "uneval")) {
    user_mapping <- mapping
  } else {
    user_mapping <- aes()
  }

  # Apply color mapping based on colorBy option
  if (colorBy == "chr" && !is_regional) {
    color_mapping <- aes(color = .data$color_group)
  } else if (colorBy == "r2") {
    if (is.null(plot_data$r2_value)) {
      stop("r2 values must be provided in 'data' or via the 'r2' parameter when colorBy is 'r2'.")
    }
    color_mapping <- aes(color = .data$r2_value)
  } else {
    # colorBy == "none" or regional without r2
    color_mapping <- NULL
  }

  # Combine mappings
  current_mapping <- default_mapping

  # Add color mapping if not overridden by user
  if (is.null(user_mapping$colour) && is.null(user_mapping$color) && !is.null(color_mapping)) {
    current_mapping$colour <- color_mapping$colour
  }

  # Add user mappings, overriding defaults where specified
  for (aes_name in names(user_mapping)) {
    current_mapping[[aes_name]] <- user_mapping[[aes_name]]
  }

  # Create a list to hold all layers
  layer_list <- list()

  # Main points layer
  point_params <- list(na.rm = na.rm, size = size, ...)
  if (colorBy == "none" || (is_regional && is.null(color_mapping))) {
    point_params$color <- color
  }

  layer_list$points <- ggplot2::layer(
    data = plot_data,
    mapping = current_mapping,
    stat = stat,
    geom = GeomPoint,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = point_params
  )

  # Add color scales based on colorBy
  if (colorBy == "chr" && !is_regional) {
    layer_list$color_scale <- ggplot2::scale_color_manual(values = colors, guide = "none")
  } else if (colorBy == "r2") {
    layer_list$color_scale <- ggplot2::scale_color_gradientn(
      colors = c("blue2", "skyblue", "green", "orange", "red2"),
      values = scales::rescale(c(0, 0.2, 0.4, 0.6, 0.8, 1)),
      limits = c(0, 1),
      na.value = "grey50",
      name = expression(r^2)
    )
  }

  # Add axis formatting based on mode
  if (is_regional) {
    # REGIONAL MODE: Use scale_x_genome_region for consistent styling
    if (!is.null(region)) {
      layer_list$x_scale <- scale_x_genome_region(region)
    } else {
      # No explicit region provided, use data range
      layer_list$x_scale <- scale_x_genome(
        limits = c(min(plot_data$BP), max(plot_data$BP))
      )
    }

    # X-axis label for regional mode
    if (is.null(x_axis_label)) {
      chr_name <- unique(plot_data$CHR)[1]
      # Remove "chr" prefix if present for cleaner display
      chr_display <- gsub("^chr", "", as.character(chr_name), ignore.case = TRUE)
      x_axis_label <- paste0("Chr", chr_display)
    }
  } else {
    # GENOME-WIDE MODE: Chromosome-labeled x-axis
    layer_list$x_scale <- ggplot2::scale_x_continuous(
      labels = axis_df$CHR,
      breaks = axis_df$center,
      expand = c(0.01, 0)
    )

    if (is.null(x_axis_label)) {
      x_axis_label <- "Chromosome"
    }
  }

  layer_list$x_label <- ggplot2::xlab(x_axis_label)

  # Y-axis scale and label
  layer_list$y_scale <- ggplot2::scale_y_continuous(
    expand = c(0.02, 0),
    limits = c(0, max(plot_data$logp) * 1.05)
  )

  if (is.null(y_axis_label)) {
    y_axis_label <- expression(paste("-log"[10], "(P)"))
  }
  layer_list$y_label <- ggplot2::ylab(y_axis_label)

  # Add threshold line if specified
  if (!is.null(threshold_p)) {
    threshold_y <- if (logp) -log10(threshold_p) else threshold_p
    layer_list$threshold <- ggplot2::geom_hline(
      yintercept = threshold_y,
      color = threshold_color,
      linetype = threshold_linetype
    )
  }

  # Add highlighted points (lead SNPs)
  if (!is.null(lead_snp) && !is.null(plot_data$SNP)) {
    highlight_data <- plot_data |> dplyr::filter(.data$SNP %in% lead_snp)
    if (nrow(highlight_data) > 0) {
      layer_list$lead_snps <- ggplot2::geom_point(
        data = highlight_data,
        mapping = aes(x = .data$BP, y = .data$logp),
        color = highlight_color,
        shape = highlight_shape,
        size = size * 2
      )
    }
  }

  # Add additional highlight SNPs from data frame
  if (!is.null(highlight_snps)) {
    # Detect columns in highlight_snps using same logic
    hl_chr <- detect_column(highlight_snps, c("CHR", "chr", "seqnames", "chrom"), "chromosome")
    hl_bp <- detect_column(highlight_snps, c("BP", "bp", "start", "pos", "position"), "position")
    hl_p <- detect_column(highlight_snps, c("P", "p", "pvalue", "p.value", "pval"), "p-value")

    hl_data <- highlight_snps |>
      dplyr::select(CHR = .data[[hl_chr]], BP = .data[[hl_bp]], P = .data[[hl_p]])

    if (logp) {
      hl_data <- hl_data |> dplyr::mutate(logp = -log10(.data$P))
    } else {
      hl_data <- hl_data |> dplyr::mutate(logp = .data$P)
    }

    # For genome-wide mode, need to adjust BP positions
    if (!is_regional && exists("data_cum")) {
      hl_data <- hl_data |>
        dplyr::inner_join(data_cum, by = "CHR") |>
        dplyr::mutate(BP = .data$BP + .data$BP_add)
    }

    if (nrow(hl_data) > 0) {
      layer_list$highlight_snps <- ggplot2::geom_point(
        data = hl_data,
        mapping = aes(x = .data$BP, y = .data$logp),
        color = highlight_color,
        shape = highlight_shape,
        size = size * 3
      )
    }
  }

  return(layer_list)
}
