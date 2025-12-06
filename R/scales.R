#' Scale for genomic coordinates on the x-axis
#'
#' This function creates a continuous scale for genomic coordinates on the x-axis.
#' It formats the axis labels in a genomic coordinate style (e.g., 1Mb, 500kb).
#'
#' @param ... Additional arguments passed to scale_x_continuous
#' @param unit_suffix Suffix to use for the unit (default: "b" for base pairs)
#' @param breaks Breaks for the axis (default: waiver())
#' @param labels Labels for the axis (default: waiver())
#' @return A ggplot2 scale object
#' @export
#' @importFrom ggplot2 scale_x_continuous waiver
#' @importFrom scales label_number
#' @examples
#' \dontrun{
#' library(ggplot2)
#' p <- ggplot(data, aes(x = start)) + geom_point() + scale_x_genome()
#' }
scale_x_genome <- function(
  ...,
  unit_suffix = "b",
  breaks = waiver(),
  labels = waiver()
) {
  if (is.function(labels)) {
    # Use provided label function
    label_fun <- labels
  } else if (inherits(labels, "waiver")) {
    # Create a function to format genomic coordinates
    label_fun <- function(x) {
      scales::label_number(
        scale_cut = scales::cut_short_scale(),
        unit = unit_suffix
      )(x)
    }
  } else {
    # Use provided labels
    label_fun <- labels
  }

  ggplot2::scale_x_continuous(
    ...,
    breaks = breaks,
    labels = label_fun,
    expand = c(0, 0)
  )
}

#' Hi-C color palettes
#'
#' Predefined color palettes commonly used for Hi-C contact matrix visualization.
#'
#' @param n Number of colors to generate (default: 256)
#' @return A vector of colors
#' @name hic_palettes
#' @export
#' @examples
#' \dontrun{
#' # Get 100 colors from the cooler palette
#' cols <- hic_palette_cooler(100)
#' }
hic_palette_cooler <- function(n = 256) {
  grDevices::colorRampPalette(c(
    "white",
    "#fee5d9",
    "#fcae91",
    "#fb6a4a",
    "#de2d26",
    "#a50f15"
  ))(n)
}

#' @rdname hic_palettes
#' @export
hic_palette_ylgnbu <- function(n = 256) {
  grDevices::colorRampPalette(c(
    "#ffffd9",
    "#edf8b1",
    "#c7e9b4",
    "#7fcdbb",
    "#41b6c4",
    "#1d91c0",
    "#225ea8",
    "#0c2c84"
  ))(n)
}

#' @rdname hic_palettes
#' @export
hic_palette_viridis <- function(n = 256) {
  if (requireNamespace("viridisLite", quietly = TRUE)) {
    viridisLite::viridis(n)
  } else {
    grDevices::hcl.colors(n, "viridis")
  }
}

#' @rdname hic_palettes
#' @export
hic_palette_bwr <- function(n = 256) {
  # Blue-white-red diverging palette for log2 fold changes
  grDevices::colorRampPalette(c(
    "#2166ac",
    "#67a9cf",
    "#d1e5f0",
    "white",
    "#fddbc7",
    "#ef8a62",
    "#b2182b"
  ))(n)
}

#' Scale for Hi-C contact matrix fill colors
#'
#' This function creates a continuous fill scale optimized for Hi-C contact matrix
#' visualization. It provides common color palettes and supports log transformations.
#'
#' @param palette Character string specifying the color palette. Options:
#'   "cooler" (red gradient, default), "ylgnbu" (yellow-green-blue),
#'   "viridis", or "bwr" (blue-white-red for diverging data)
#' @param trans Transformation for the scale. Options: "identity" (linear, default),
#'   "log10", "log2", "sqrt". Can also pass a transformation object.
#' @param limits Numeric vector of length 2 giving the scale limits. Values outside
#'   this range will be squished to the nearest limit. Default: NULL (auto)
#' @param na.value Color for missing values (default: "grey50")
#' @param midpoint For diverging palettes ("bwr"), the midpoint value (default: 0)
#' @param oob Function to handle out-of-bounds values. Default: scales::squish
#' @param ... Additional arguments passed to ggplot2::scale_fill_gradientn or
#'   scale_fill_gradient2
#'
#' @return A ggplot2 scale object
#' @export
#' @importFrom ggplot2 scale_fill_gradientn scale_fill_gradient2
#' @importFrom scales squish
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' # Basic usage with default cooler palette
#' ggplot(hic_data, aes(x = x, y = y, fill = score)) +
#'   geom_hic() +
#'   scale_fill_hic()
#'
#' # Log10 transformation
#' ggplot(hic_data, aes(x = x, y = y, fill = score)) +
#'   geom_hic() +
#'   scale_fill_hic(trans = "log10")
#'
#' # Different palette with limits
#' ggplot(hic_data, aes(x = x, y = y, fill = score)) +
#'   geom_hic() +
#'   scale_fill_hic(palette = "ylgnbu", limits = c(0, 100))
#'
#' # Diverging palette for comparison data
#' ggplot(diff_data, aes(x = x, y = y, fill = log2fc)) +
#'   geom_hic() +
#'   scale_fill_hic(palette = "bwr", midpoint = 0)
#' }
scale_fill_hic <- function(
  palette = c("cooler", "ylgnbu", "viridis", "bwr"),
  trans = "identity",
  limits = NULL,
  na.value = "grey50",
  midpoint = 0,
  oob = scales::squish,
  ...
) {
  palette <- match.arg(palette)

  # Get colors based on palette choice
  colors <- switch(
    palette,
    "cooler" = hic_palette_cooler(256),
    "ylgnbu" = hic_palette_ylgnbu(256),
    "viridis" = hic_palette_viridis(256),
    "bwr" = hic_palette_bwr(256)
  )

  # Handle transformation
  if (is.character(trans)) {
    trans_obj <- switch(
      trans,
      "identity" = scales::identity_trans(),
      "log10" = scales::log10_trans(),
      "log2" = scales::log2_trans(),
      "sqrt" = scales::sqrt_trans(),
      scales::identity_trans() # fallback
    )
  } else {
    trans_obj <- trans
  }

  # For diverging palette (bwr), use gradient2

  if (palette == "bwr") {
    ggplot2::scale_fill_gradient2(
      low = colors[1],
      mid = "white",
      high = colors[length(colors)],
      midpoint = midpoint,
      limits = limits,
      oob = oob,
      na.value = na.value,
      trans = trans_obj,
      ...
    )
  } else {
    ggplot2::scale_fill_gradientn(
      colours = colors,
      limits = limits,
      oob = oob,
      na.value = na.value,
      trans = trans_obj,
      ...
    )
  }
}

#' Y-axis scale for strand tracks
#'
#' Provides y-axis labels for strand tracks produced by geom_gene.
#' Now uses discrete scale to match the factor levels created by geom_gene.
#' The default labels are "- strand" and "+ strand" as created by geom_gene.
#'
#' @param labels Character vector of labels to display. If NULL, uses the factor levels from geom_gene.
#' @param ... Additional arguments passed to ggplot2::scale_y_discrete.
#' @return A ggplot2 scale object for the y-axis.
#' @export
#' @importFrom ggplot2 scale_y_discrete
#' @examples
#' \dontrun{
#' library(ggplot2)
#' p <- ggplot(df) + geom_gene(...) + scale_y_strand()
#' # Or customize labels:
#' p <- ggplot(df) + geom_gene(...) + scale_y_strand(labels = c("-", "+", "?"))
#' }
scale_y_strand <- function(labels = NULL, ...) {
  ggplot2::scale_y_discrete(labels = labels, ...)
}

#' Format genomic coordinates
#'
#' This function formats genomic coordinates in a human-readable way (e.g., 1Mb, 500kb).
#'
#' @param x Numeric vector of genomic coordinates
#' @param unit_suffix Suffix to use for the unit (default: "b" for base pairs)
#' @return Character vector of formatted coordinates
#' @export
#' @importFrom scales label_number
#' @examples
#' \dontrun{
#' format_genomic_coord(c(1000, 1000000, 1500000))
#' # [1] "1kb" "1Mb" "1.5Mb"
#' }
format_genomic_coord <- function(x, unit_suffix = "b") {
  scales::label_number(
    scale_cut = scales::cut_short_scale(),
    unit = unit_suffix
  )(x)
}

#' Create a genomic position scale for a specific chromosome region
#'
#' This function creates a continuous scale for a specific chromosome region.
#' It sets the limits to the start and end positions of the region and formats
#' the axis labels in a genomic coordinate style.
#'
#' @param region A string specifying a genomic region (e.g., "chr1:1000000-2000000") or a GRanges object
#' @param ... Additional arguments passed to scale_x_genome
#' @return A ggplot2 scale object
#' @export
#' @importFrom GenomicRanges start end
#' @importFrom methods is
#' @examples
#' \dontrun{
#' library(ggplot2)
#' p <- ggplot(data, aes(x = start)) + geom_point() +
#'      scale_x_genome_region("chr1:1000000-2000000")
#' }
scale_x_genome_region <- function(region, ...) {
  # Parse region if it's a string
  if (is.character(region)) {
    region_gr <- parse_region(region)
  } else if (methods::is(region, "GRanges")) {
    region_gr <- region
  } else {
    stop("Region must be a character string or GRanges object")
  }

  # Extract start and end positions
  start_pos <- GenomicRanges::start(region_gr)
  end_pos <- GenomicRanges::end(region_gr)

  # Create scale with limits set to the region
  scale_x_genome(..., limits = c(start_pos, end_pos))
}
