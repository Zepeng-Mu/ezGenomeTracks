#' Geom for genomic feature tracks
#'
#' This function creates a geom for genomic feature tracks, such as ChIP-seq features,
#' ATAC-seq features, or any interval-based features. It displays the features as rectangles.
#'
#' @inheritParams ggplot2::geom_rect
#' @param mapping Aesthetic mappings (default: NULL). xmin, xmax are required.
#' @param data Dataset (default: NULL)
#' @param stat Statistic to use (default: "identity")
#' @param position Position adjustment (default: "identity")
#' @param height Height of the features (default: 0.8)
#' @param color Border color of the features (default: "#05b1d3")
#' @param fill Fill color of the features (default: "#05b1d3")
#' @param alpha Transparency (default: 0.7)
#' @return A ggplot2 layer
#' @export
#' @importFrom ggplot2 geom_rect aes
#' @examples
#' \dontrun{
#' library(ggplot2)
#' p <- ggplot(feature_data) +
#'   geom_feature(aes(xmin = start, xmax = end, fill = score))
#' }
geom_feature <- function(mapping = NULL, data = NULL, stat = "identity",
                         position = "identity", ..., height = 0.8,
                         color = "#05b1d3", fill = "#05b1d3", alpha = 0.7,
                         show.legend = NA, inherit.aes = TRUE) {
  # Default mapping for features
  if (is.null(mapping)) {
    mapping <- ggplot2::aes(xmin = start, xmax = end, ymin = 0, ymax = height)
  } else {
    # If mapping is provided but doesn't include y positions, add them
    if (!all(c("ymin", "ymax") %in% names(mapping))) {
      mapping$ymin <- 0
      mapping$ymax <- height
    }
  }

  # Create the rectangle geom
  ggplot2::geom_rect(
    mapping = mapping, data = data, stat = stat,
    position = position, color = color, fill = fill, alpha = alpha, ...,
    show.legend = show.legend, inherit.aes = inherit.aes
  )
}
