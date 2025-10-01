#' Geom for genomic feature tracks
#'
#' This function creates a geom for genomic feature tracks, such as ChIP-seq features,
#' ATAC-seq features, or any interval-based features. It displays the features as rectangles.
#'
#' @inheritParams ggplot2::layer
#' @param mapping Aesthetic mappings (default: NULL). xmin, xmax are required.
#' @param data Dataset (default: NULL)
#' @param stat Statistic to use (default: "identity")
#' @param position Position adjustment (default: "identity")
#' @param height Height of the features (default: 0.8)
#' @param color Border color of the features (default: "#05b1d3")
#' @param fill Fill color of the features (default: "#05b1d3")
#' @param alpha Transparency (default: 0.7)
#' @param na.rm If `TRUE`, silently removes `NA` values.
#' @param ... Other arguments passed on to \code{\link[ggplot2]{layer}}. These are
#'   often aesthetics, used to set an aesthetic to a fixed value.
#' @return A ggplot2 layer
#' @export
#' @importFrom ggplot2 GeomRect layer aes ggproto Geom
#' @examples
#' \dontrun{
#' library(ggplot2)
#' p <- ggplot(feature_data) +
#'   geom_feature(aes(xmin = start, xmax = end, fill = score))
#' }
geom_feature <- function(mapping = NULL, data = NULL, stat = "identity",
                         position = "identity", ..., height = 0.8,
                         color = "#05b1d3", fill = "#05b1d3", alpha = 0.7,
                         na.rm = TRUE, show.legend = NA, inherit.aes = TRUE) {
  # Default mapping for features
  default_aes <- aes(xmin = .data$start, xmax = .data$end, ymin = 0, ymax = height)

  if (is.null(mapping)) {
    mapping <- default_aes
  } else {
    mapping <- utils::modifyList(default_aes, as.list(mapping))
    mapping <- do.call(aes, mapping)
  }

  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomFeature,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      height = height,
      color = color,
      fill = fill,
      alpha = alpha,
      na.rm = na.rm,
      ...
    )
  )
}

#' @rdname geom_feature
#' @format NULL
#' @usage NULL
GeomFeature <- ggproto("GeomFeature", GeomRect,
  required_aes = c("xmin", "xmax"),
  setup_data = function(data, params) {
    # If ymin and ymax are not provided, set them based on height
    if (!all(c("ymin", "ymax") %in% names(data))) {
      data$ymin <- 0
      data$ymax <- params$height
    }
    data
  },
  default_aes = aes(
    colour = "#05b1d3", fill = "#05b1d3",
    linewidth = 0.5, linetype = 1, alpha = 0.7
  )
)
