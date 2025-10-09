#' Geom for genomic interaction arcs
#'
#' This function creates a geom for genomic interaction arcs, such as chromatin loops,
#' eQTL links, or sashimi plots. It displays the interactions as bezier curves.
#'
#' @inheritParams ggplot2::layer
#' @param curvature Amount of curvature (default: 0.5)
#' @param arrow_length Length of directional arrows (default: 0)
#' @param arrow_type Type of arrow (default: "closed")
#' @param na.rm If `TRUE`, silently drop `NA` values.
#' @param ... Additional arguments passed to [ggplot2::layer()], e.g.
#'   `color = "black"`, `linewidth = 0.8`, or `alpha = 0.6`.
#'
#' @return A ggplot2 layer that can be added to a plot.
#' @export
#' @importFrom ggplot2 GeomCurve layer aes ggproto Geom arrow unit
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' df <- data.frame(
#'   start1 = c(1000, 2000, 3000),
#'   end1 = c(1100, 2100, 3100),
#'   start2 = c(5000, 6000, 7000),
#'   end2 = c(5100, 6100, 7100),
#'   score = c(0.8, 0.6, 0.9)
#' )
#'
#' # Basic arc plot
#' ggplot(df) +
#'   geom_arc(aes(x = start1, y = 0, xend = start2, yend = 0))
#'
#' # With curvature and arrows
#' ggplot(df) +
#'   geom_arc(aes(x = start1, y = 0, xend = start2, yend = 0, color = score),
#'     curvature = 0.8, arrow_length = 0.1
#'   )
#' }
geom_arc <- function(mapping = NULL, data = NULL, stat = "identity",
                     position = "identity", ..., curvature = 0.5,
                     arrow_length = 0, arrow_type = "closed",
                     na.rm = TRUE, show.legend = NA, inherit.aes = TRUE) {
  # Default mapping for arcs
  default_aes <- aes(x = .data$start1, y = 0, xend = .data$start2, yend = 0)

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
    geom = GeomArc,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      curvature = curvature,
      arrow_length = arrow_length,
      arrow_type = arrow_type,
      na.rm = na.rm,
      ...
    )
  )
}

#' @rdname geom_arc
#' @format NULL
#' @usage NULL
GeomArc <- ggproto("GeomArc", Geom,
  required_aes = c("x", "y", "xend", "yend"),
  setup_params = function(data, params) {
    # Set up arrow if requested
    if (params$arrow_length > 0) {
      params$arrow <- ggplot2::arrow(
        type = params$arrow_type,
        length = ggplot2::unit(params$arrow_length, "inches")
      )
    } else {
      params$arrow <- NULL
    }
    params
  },
  draw_panel = function(data, panel_params, coord, curvature = 0.5,
                        arrow_length = 0, arrow_type = "closed", na.rm = FALSE) {
    # Use GeomCurve for the actual drawing
    GeomCurve$draw_panel(data, panel_params, coord, curvature = curvature, na.rm = na.rm)
  },
  default_aes = aes(
    colour = "gray50", linewidth = 0.5, linetype = 1,
    alpha = 0.7
  )
)

#' Process interaction data for visualization
#'
#' This function processes interaction data for visualization with geom_arc.
#' It formats the data for use with ggplot2.
#'
#' @param gr A GRanges object with interaction data
#' @param anchor1 Column name for the first anchor (default: "anchor1")
#' @param anchor2 Column name for the second anchor (default: "anchor2")
#' @param score Column name for interaction score (default: "score")
#' @return A data frame with interaction information
#' @export
#' @importFrom GenomicRanges start end
#' @importFrom S4Vectors mcols
#' @examples
#' \dontrun{
#' library(rtracklayer)
#' gr <- import("interactions.bedpe")
#' interaction_data <- process_interaction_data(gr)
#' }
process_interaction_data <- function(gr, anchor1 = "anchor1", anchor2 = "anchor2",
                                     score = "score") {
  # Check if it's a GRanges object
  if (methods::is(gr, "GRanges")) {
    # Convert GRanges to data frame
    gr_df <- granges_to_df(gr)

    # Check if required columns exist
    if (!all(c(anchor1, anchor2) %in% colnames(gr_df))) {
      stop("Required columns not found in the data")
    }

    # Extract interaction information
    result <- data.frame(
      start1 = gr_df[[paste0(anchor1, "_start")]],
      end1 = gr_df[[paste0(anchor1, "_end")]],
      start2 = gr_df[[paste0(anchor2, "_start")]],
      end2 = gr_df[[paste0(anchor2, "_end")]]
    )

    # Add score if available
    if (score %in% colnames(gr_df)) {
      result$score <- gr_df[[score]]
    }

    return(result)
  } else if (methods::is(gr, "GInteractions")) {
    # Handle GInteractions objects
    anchor1_gr <- gr@anchor1
    anchor2_gr <- gr@anchor2

    result <- data.frame(
      start1 = GenomicRanges::start(anchor1_gr),
      end1 = GenomicRanges::end(anchor1_gr),
      start2 = GenomicRanges::start(anchor2_gr),
      end2 = GenomicRanges::end(anchor2_gr)
    )

    # Add score if available
    if (score %in% colnames(S4Vectors::mcols(gr))) {
      result$score <- S4Vectors::mcols(gr)[[score]]
    }

    return(result)
  } else {
    stop("Input must be a GRanges or GInteractions object")
  }
}

#' Create an interaction track from a BEDPE file
#'
#' This function creates an interaction track from a BEDPE file. It imports the data
#' for a specific region and creates a ggplot2 layer for visualization.
#'
#' @param file Path to the BEDPE file
#' @param region Genomic region to display (e.g., "chr1:1000000-2000000")
#' @param curvature Amount of curvature (default: 0.5)
#' @param color Color of the arcs (default: "gray50")
#' @param size Size of the arcs (default: 0.5)
#' @param alpha Transparency (default: 0.7)
#' @param use_score Use the score column for color (default: FALSE)
#' @param ... Additional arguments passed to geom_arc
#' @return A ggplot2 layer
#' @export
#' @importFrom ggplot2 ggplot aes scale_color_gradient
#' @examples
#' \dontrun{
#' p <- interaction_track("interactions.bedpe", "chr1:1000000-2000000", use_score = TRUE)
#' }
interaction_track <- function(file, region, curvature = 0.5, color = "gray50",
                             size = 0.5, alpha = 0.7, use_score = FALSE, ...) {
  # Parse the region
  region_gr <- parse_region(region)

  # Import the data
  interaction_gr <- rtracklayer::import(file, which = region_gr)

  # Process the interaction data
  interaction_data <- process_interaction_data(interaction_gr)

  # Create the plot
  if (use_score && "score" %in% colnames(interaction_data)) {
    p <- ggplot2::ggplot(interaction_data) +
      geom_arc(ggplot2::aes(x = start1, y = 0, xend = start2, yend = 0, color = score),
               curvature = curvature, size = size, alpha = alpha, ...) +
      ggplot2::scale_color_gradient(low = "blue", high = "red")
  } else {
    p <- ggplot2::ggplot(interaction_data) +
      geom_arc(ggplot2::aes(x = start1, y = 0, xend = start2, yend = 0),
               curvature = curvature, color = color, size = size, alpha = alpha, ...)
  }

  # Apply the appropriate theme and scale
  p <- p + ez_peak_theme() +
    scale_x_genome_region(region) +
    ggplot2::ylim(-0.5, 0.5)  # Fixed y-axis for arcs

  return(p)
}