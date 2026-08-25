# ezGenomeTracks - Helper functions (split from helpers.R)

`%||%` <- function(x, y) if (is.null(x)) y else x
detect_column <- function(data, candidates, param_name, required = TRUE) {
  for (col in candidates) {
    if (col %in% colnames(data)) return(col)
  }
  if (required) {
    stop(paste0(
      "Could not find ",
      param_name,
      " column. Expected one of: ",
      paste(candidates, collapse = ", ")
    ))
  }
  return(NULL)
}

# Internal helper: Validate that required columns exist in a data frame
# Throws an error with informative message if any required column is missing.
#
# @param data A data frame to validate
# @param required_cols Character vector of required column names
# @param data_name Display name of the data (used in error message, e.g. "Data frame", "Coverage data")
validate_required_columns <- function(data, required_cols, data_name = "Data frame") {
  if (!all(required_cols %in% colnames(data))) {
    stop(
      data_name,
      " must contain columns: ",
      paste(required_cols, collapse = ", ")
    )
  }
}

# Internal helper: Apply consistent facet positioning to ggplot objects
# Handles both "left" and "top" positioning with appropriate theme settings.
#
# @param plot A ggplot object to modify
# @param facet_label_position Either "left" or "top"
# @param strip_placement Either "inside" or "outside" (default: "inside")
# @return Modified ggplot object with facet theming applied
apply_facet_position <- function(plot, facet_label_position, strip_placement = "inside") {
  if (facet_label_position == "left") {
    plot <- plot +
      ggplot2::facet_wrap(
        ~track,
        ncol = 1,
        scales = "free_y",
        strip.position = "left"
      ) +
      ggplot2::theme(
        strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 1),
        strip.placement = strip_placement,
        panel.spacing.y = ggplot2::unit(0, "pt")
      )
  } else {
    plot <- plot +
      ggplot2::facet_wrap(~track, ncol = 1, scales = "free_y") +
      ggplot2::theme(panel.spacing.y = ggplot2::unit(0, "pt"))
  }
  return(plot)
}

# Internal helper: Apply a consistent black border to plot panel
# Used consistently across wrapper functions for uniform visual appearance.
#
# @param plot A ggplot object to modify
# @return Modified ggplot object with black border applied
apply_border_theme <- function(plot) {
  plot +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.5
      )
    )
}

# Internal helper: Resolve plot colors for multi-track/grouped data
# When a single color is provided but multiple groups exist, falls back to
# the default palette. Otherwise recycles colors via rep_len().
#
# @param colors Character vector of user-provided colors
# @param color_values Character vector of unique group/track names
# @return Named character vector of colors
resolve_plot_colors <- function(colors, color_values) {
  n_colors <- length(color_values)
  if (length(colors) == 1 && n_colors > 1) {
    plot_colors <- ez_default_palette(n_colors)
  } else {
    plot_colors <- rep_len(colors, n_colors)
  }
  names(plot_colors) <- color_values
  plot_colors
}

# Internal helper: Auto-detect GWAS column names from common conventions
# Checks for standard column names used in GWAS summary statistics.
#
# @param data Data frame to inspect
# @param chr Column name for chromosome (auto-detected if NULL)
# @param bp Column name for base pair position (auto-detected if NULL)
# @param p Column name for p-value (auto-detected if NULL)
# @param snp Column name for SNP ID (auto-detected if NULL, optional)
# @return Named list with chr, bp, p, snp column names
auto_detect_gwas_columns <- function(data, chr = NULL, bp = NULL, p = NULL, snp = NULL) {
  if (is.null(chr)) {
    chr <- detect_column(
      data, c("CHR", "chr", "seqnames", "chrom", "chromosome"), "chromosome"
    )
  }
  if (is.null(bp)) {
    bp <- detect_column(
      data, c("BP", "bp", "start", "pos", "position", "POS"), "position"
    )
  }
  if (is.null(p)) {
    p <- detect_column(
      data, c("P", "p", "pvalue", "p.value", "pval", "P.value"), "p-value"
    )
  }
  if (is.null(snp)) {
    snp <- detect_column(
      data, c("SNP", "snp", "rsid", "id", "variant_id", "marker"), "SNP",
      required = FALSE
    )
  }
  list(chr = chr, bp = bp, p = p, snp = snp)
}

# Internal helper: Assign default track names to a named list of inputs
# Used by process_*_input() functions to ensure consistent track naming.
#
# @param input A list of data (may or may not have names)
# @param track_labels Optional character vector of track labels
# @return The input list with names assigned if missing
ensure_track_names <- function(input, track_labels = NULL) {
  if (is.null(names(input)) && is.null(track_labels)) {
    names(input) <- paste0("Track ", seq_along(input))
  } else if (is.null(names(input)) && !is.null(track_labels)) {
    names(input) <- track_labels
  }
  input
}

#' Extract signal data for a single input element
