# Internal color utility functions
#
# These functions provide a unified color system for ezGenomeTracks.
# When the vibeColors package is installed, its palettes are used.
# Otherwise, hardcoded fallback colors are used.

# Hardcoded fallback palette (colorblind-friendly discrete palette)
.ez_fallback_discrete <- c(

  "#0E69B7", # Denim blue
  "#D00E12", # Red
  "#00B4A6", # Teal
  "#F0E441", # Yellow
  "#57258C", # Purple
  "#FE5400", # Orange
  "#57B4E9", # Light blue
  "#FFAF01", # Amber
  "#CB79A7", # Mauve
  "#AFAFAF", # Silver
  "#000000"  # Black
)

# Hardcoded fallback single-track colors
.ez_fallback_single <- list(
  coverage = "steelblue",
  feature  = "#05b1d3",
  gene     = "gray50",
  link     = "gray50",
  sashimi  = "purple3",
  manhattan = "grey50"
)

#' Get default discrete color palette
#'
#' Returns `n` discrete colors for multi-track or multi-group plots.
#' Uses vibeColors CBsafe1 palette if available, otherwise falls back to
#' a built-in colorblind-friendly palette.
#'
#' @param n Number of colors to return. If NULL, returns the full palette.
#' @return A character vector of hex color codes.
#' @keywords internal
ez_default_palette <- function(n = NULL) {
  colors <- if (requireNamespace("vibeColors", quietly = TRUE)) {
    vibeColors::vibe_palette("CBsafe1")
  } else {
    .ez_fallback_discrete
  }

  if (is.null(n)) return(unname(colors))

  if (n <= length(colors)) {
    return(unname(colors[seq_len(n)]))
  }

  # Interpolate if more colors needed than available
  if (requireNamespace("vibeColors", quietly = TRUE)) {
    unname(vibeColors::vibe_palette("CBsafe1", n = n))
  } else {
    grDevices::colorRampPalette(colors)(n)
  }
}

#' Get default single-track color
#'
#' Returns the default color for a single track of the given type.
#' Uses vibeColors named colors if available, otherwise falls back to
#' built-in defaults.
#'
#' @param type Character. One of "coverage", "feature", "gene", "link",
#'   "sashimi", "manhattan".
#' @return A single hex color code string.
#' @keywords internal
ez_default_single_color <- function(
  type = c("coverage", "feature", "gene", "link", "sashimi", "manhattan")
) {
  type <- match.arg(type)

  if (requireNamespace("vibeColors", quietly = TRUE)) {
    switch(type,
      coverage  = vibeColors::vibe_color("Azure500"),
      feature   = vibeColors::vibe_color("Turquoise400"),
      sashimi   = vibeColors::vibe_color("Plum500"),
      .ez_fallback_single[[type]]
    )
  } else {
    .ez_fallback_single[[type]]
  }
}

#' Get Hi-C sequential palette colors
#'
#' Returns colors for Hi-C sequential palettes (e.g., cooler).
#' Uses vibeColors Crimson palette if available, otherwise falls back.
#'
#' @param name Character. Palette name: "cooler" or "ylgnbu".
#' @param n Integer. Number of colors to generate.
#' @return A character vector of hex color codes.
#' @keywords internal
ez_hic_palette <- function(name = c("cooler", "ylgnbu"), n = 256) {
  name <- match.arg(name)

  if (name == "cooler" && requireNamespace("vibeColors", quietly = TRUE)) {
    return(vibeColors::vibe_palette("Crimson", n = n))
  }

  # Fallback to hardcoded palettes
  NULL
}

#' Get Hi-C diverging palette colors
#'
#' Returns colors for the Hi-C blue-white-red diverging palette.
#' Uses vibeColors AzureCrimson palette if available.
#'
#' @param n Integer. Number of colors to generate.
#' @return A character vector of hex color codes, or NULL to use fallback.
#' @keywords internal
ez_hic_diverging_palette <- function(n = 256) {
  if (requireNamespace("vibeColors", quietly = TRUE)) {
    return(vibeColors::vibe_palette("AzureCrimson", n = n))
  }
  NULL
}

# UCSC Genome Browser nucleotide color conventions (A, C, G, T, N)
.ucsc_nucleotide_colors <- c(
  A = "#00A800",
  C = "#0000CC",
  G = "#CC9900",
  T = "#CC0000",
  N = "#999999"
)

#' UCSC nucleotide color palette
#'
#' Returns the UCSC Genome Browser nucleotide color conventions as a named
#' character vector. Colors follow standard bioinformatics conventions:
#' A (green), C (blue), G (gold), T (red), N (grey).
#'
#' @param colors Optional named character vector of color overrides (e.g.,
#'   `c(A = "purple")`). Unspecified bases retain their UCSC defaults.
#' @return A named character vector with hex color codes for A, C, G, T, N.
#' @export
#' @examples
#' # Default UCSC colors
#' ez_sequence_palette()
#'
#' # Override A to purple
#' ez_sequence_palette(colors = c(A = "purple"))
ez_sequence_palette <- function(colors = NULL) {
  palette <- .ucsc_nucleotide_colors
  if (!is.null(colors)) {
    valid <- names(colors) %in% names(palette)
    if (!all(valid)) {
      warning(
        "Ignoring unknown nucleotide names in 'colors': ",
        paste(names(colors)[!valid], collapse = ", ")
      )
    }
    palette[names(colors)[valid]] <- colors[valid]
  }
  palette
}
