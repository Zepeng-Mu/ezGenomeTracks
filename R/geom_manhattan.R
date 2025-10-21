#' @title Create a Manhattan plot
#' @description `geom_manhattan` creates a Manhattan plot, which can be used to visualize both GWAS and QTL.
#' @param mapping Set of aesthetic mappings created by `aes()`. If not specified, the default mappings are used.
#' @param data A `data.frame` containing the data to be plotted. Must include columns for chromosome, base pair position, and p-value. Optional column for SNP identifier.
#' @param stat The statistical transformation to apply to the data (default: "identity").
#' @param position Position adjustment, either as a string naming a position adjustment function, or the result of a call to a position adjustment function.
#' @param na.rm If `FALSE` (default), missing values are removed with a warning. If `TRUE`, missing values are silently removed.
#' @param show.legend Logical. Should this layer be displayed in the legend? `NA` for automatic, `TRUE` always, `FALSE` never.
#' @param inherit.aes If `FALSE`, overrides the default aesthetics; if `TRUE`, inherits them.
#' @param chr Name of the chromosome column in `data` (default: "CHR").
#' @param bp Name of the base pair position column in `data` (default: "BP").
#' @param p Name of the p-value column in `data` (default: "P").
#' @param snp Name of the SNP identifier column in `data` (default: "SNP").
#' @param logp Logical. If -log10() is used for p.
#' @param size Point size (default: 0.5).
#' @param lead.snp Vector of SNP to highlight.
#' @param r2 Vector of R-squared values for linkage disequilibrium (LD) coloring. Should be in same order as SNP.
#' @param colors Vector of colors to use for alternating chromosome colors (default: c("grey", "skyblue")).
#' @param highlight_snps Data frame of SNPs to highlight, with columns `CHR`, `BP`, and `P`.
#' @param highlight_color Color for highlighted SNPs (default: "purple").
#' @param highlight_shape Shape for highlighted SNPs (default: 18).
#' @param threshold_p A numeric value for the p-value threshold to draw a horizontal line (e.g., 5e-8).
#' @param threshold_color Color for the threshold line (default: "red").
#' @param threshold_linetype Linetype for the threshold line (default: 2).
#' @param colorBy Character string indicating how points should be colored. Options are "chr" (default, alternating chromosome colors) or "r2" (continuous color based on R-squared values).
#' @param x_axis_label X-axis label (default: NULL).
#' @param y_axis_label_prefix Y-axis label prefix added before `expression(paste("-log"[10], "(P)"))`). Default: NULL.
#' @param y_axis_label Label for the y-axis, overwrites y_axis_label_prefix. Default: NULL.
#' @param ... Additional arguments passed to `ggplot2::geom_point()`.
#' @return A `ggplot2` layer.
#' @export
#' @importFrom ggplot2 layer ggproto aes GeomPoint
#' @importFrom dplyr mutate filter
#' @importFrom rlang .data
geom_manhattan <- function(
    mapping = NULL, data = NULL, stat = "identity", position = "identity",
    na.rm = FALSE, show.legend = NA, inherit.aes = TRUE,
    chr = "CHR", bp = "BP", p = "P", snp = "SNP", logp = TRUE, size = 0.5,
    lead.snp = NULL, r2 = NULL, colors = c("grey", "skyblue"),
    highlight_snps = NULL, highlight_color = "purple", highlight_shape = 18,
    threshold_p = NULL, threshold_color = "red", threshold_linetype = 2,
    colorBy = "chr",
    x_axis_label = NULL, y_axis_label_prefix = NULL, y_axis_label = NULL, ...
) {
  # Validate input data
  if (is.null(data)) stop("Data cannot be NULL.")
  if (!all(c(chr, bp, p) %in% colnames(data))) {
    stop(paste0("Data must contain columns: ", paste(c(chr, bp, p), collapse = ", ")))
  }

  # Prepare data
  plot_data <- data %>%
    dplyr::select(CHR = .data[[chr]], BP = .data[[bp]], P = .data[[p]]) %>%
    dplyr::arrange(CHR, BP)

  if (!is.null(data[[snp]])) {
    plot_data$SNP <- data[[snp]]
  }

  if (!is.null(r2)) {
    if (length(r2) != nrow(data)) {
      stop("Length of r2 vector must match the number of rows in data.")
    }
    plot_data$r2_value <- r2
  }

  if (logp) {
    plot_data <- plot_data %>% dplyr::mutate(logp = -log10(.data$P))
  } else {
    plot_data <- plot_data %>% dplyr::mutate(logp = .data$P)
  }

  if (length(unique(plot_data$CHR)) > 1) {
    # https://danielroelfs.com/posts/how-i-create-manhattan-plots-using-ggplot
    data_cum <- plot_data %>%
      dplyr::group_by(CHR) %>%
      dplyr::summarise(max_BP = max(BP)) %>%
      dplyr::mutate(BP_add = dplyr::lag(cumsum(max_BP), default = 0)) %>%
      dplyr::select(CHR, BP_add)

    plot_data <- plot_data %>%
      dplyr::inner_join(data_cum, by = "CHR") %>%
      dplyr::mutate(BP = BP + BP_add)

    # Calculate chromosome midpoints for x-axis labels
    axis_df <- plot_data %>%
      dplyr::group_by(CHR) %>%
      dplyr::summarize(center = mean(BP))
  } else {
    axis_df <- plot_data %>%
      dplyr::mutate(center = .data$BP)
  }

  # Assign colors based on chromosome
  plot_data$color_group <- factor(plot_data$CHR %% 2 + 1, levels = c(1, 2))

  # Default mapping
  default_mapping <- aes(x = .data$BP, y = .data$logp)

  # Create mapping based on colorBy option without using modifyList
  if (!is.null(mapping) && inherits(mapping, "uneval")) {
    # Extract user-provided mappings
    user_mapping <- mapping
  } else {
    user_mapping <- aes()
  }

  # Apply color mapping based on colorBy option
  if (colorBy == "chr") {
    color_mapping <- aes(color = .data$color_group)
  } else if (colorBy == "r2") {
    if (is.null(plot_data$r2_value)) {
      stop("r2 values must be provided in 'data' or via the 'r2' parameter when colorBy is 'r2'.")
    }
    color_mapping <- aes(color = .data$r2_value)
  } else {
    stop("Invalid 'colorBy' option. Choose 'chr' or 'r2'.")
  }

  # Combine mappings without using modifyList
  # Start with default mapping
  current_mapping <- default_mapping

  # Add color mapping if not overridden by user
  if (is.null(user_mapping$colour) && is.null(user_mapping$color)) {
    current_mapping$colour <- color_mapping$colour
  }

  # Add user mappings, overriding defaults where specified
  for (aes_name in names(user_mapping)) {
    current_mapping[[aes_name]] <- user_mapping[[aes_name]]
  }

  # Y-axis
  if (!is.null(y_axis_label)) {
    y_axis_label <- y_axis_label
  } else if (is.null(y_axis_label) & !is.null(y_axis_label_prefix)) {
    y_axis_label <- paste0(y_axis_label_prefix, expression(paste("-log"[10], "(P)")))
  }

  # Create a list to hold all layers
  layer_list <- list()

  # Main points layer
  layer_list$points <- ggplot2::layer(
    data = plot_data,
    mapping = current_mapping,
    stat = stat,
    geom = GeomPoint,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, size = size, ...)
  )

  # Add color scales based on colorBy
  if (colorBy == "chr") {
    layer_list$color_scale <- ggplot2::scale_color_manual(values = colors, guide = "none")
  } else if (colorBy == "r2") {
    layer_list$color_scale <- ggplot2::scale_color_gradientn(
      colors = c("blue2", "skyblue", "green", "orange", "red2"),
      values = c(0, 0.2, 0.4, 0.6, 0.8, 1),
      na.value = "grey50",
      name = expression(r^2)
    )
  }

  # Add axis formatting
  layer_list$x_scale <- ggplot2::scale_x_continuous(
    label = axis_df$CHR,
    breaks = axis_df$center,
    expand = c(0, 0)
  )

  layer_list$y_scale <- ggplot2::scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, max(plot_data$logp) * 1.05)
  )

  if (!is.null(x_axis_label)) {
    layer_list$x_label <- ggplot2::xlab(x_axis_label)
  } else if (is.null(x_axis_label) & length(unique(plot_data$CHR)) == 1) {
    layer_list$x_label <- ggplot2::xlab(stringr::str_glue("Chr{unique(plot_data$CHR)}"))
  } else if (is.null(x_axis_label) & length(unique(plot_data$CHR)) > 1) {
    layer_list$x_label <- ggplot2::xlab(stringr::str_glue("Chr"))
  }

  layer_list$y_label <- ggplot2::ylab(y_axis_label)

  # Add threshold line if specified
  if (!is.null(threshold_p)) {
    layer_list$threshold <- ggplot2::geom_hline(
      yintercept = -log10(threshold_p),
      color = threshold_color,
      linetype = threshold_linetype
    )
  }

  # Add highlighted points in a more efficient way
  if (!is.null(lead.snp) && !is.null(plot_data$SNP)) {
    highlight_data <- plot_data %>% dplyr::filter(.data$SNP %in% lead.snp)
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

  # Add additional highlight SNPs
  if (!is.null(highlight_snps)) {
    highlight_data <- highlight_snps %>%
      dplyr::select(CHR = .data[[chr]], BP = .data[[bp]], P = .data[[p]]) %>%
      dplyr::mutate(logp = -log10(.data$P)) %>%
      dplyr::inner_join(plot_data %>% dplyr::select(CHR, BP), by = c("CHR", "BP"))

    if (nrow(highlight_data) > 0) {
      layer_list$highlight_snps <- ggplot2::geom_point(
        data = highlight_data,
        mapping = aes(x = .data$BP, y = .data$logp),
        color = highlight_color,
        shape = highlight_shape,
        size = size * 3
      )
    }
  }

  return(layer_list)
}
