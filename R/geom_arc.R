#' Geom for genomic interaction arcs
#'
#' This function creates a geom for genomic interaction arcs, such as chromatin loops,
#' eQTL links, or sashimi plots. It displays the interactions as bezier curves.
#'
#' @inheritParams ggplot2::geom_curve
#' @param curvature Amount of curvature (default: 0.5)
#' @param arrow_length Length of directional arrows (default: 0)
#' @param arrow_type Type of arrow (default: "closed")
#' @param color Color of the arcs (default: "gray50")
#' @param size Size of the arcs (default: 0.5)
#' @param alpha Transparency (default: 0.7)
#' @return A ggplot2 layer
#' @export
#' @importFrom ggplot2 geom_curve arrow unit aes
#' @examples
#' \dontrun{
#' library(ggplot2)
#' p <- ggplot(interaction_data) + 
#'      geom_arc(aes(x = start1, y = 0, xend = start2, yend = 0, color = score))
#' }
geom_arc <- function(mapping = NULL, data = NULL, stat = "identity",
                      position = "identity", ..., curvature = 0.5,
                      arrow_length = 0, arrow_type = "closed",
                      color = "gray50", size = 0.5, alpha = 0.7,
                      show.legend = NA, inherit.aes = TRUE) {
  
  # Set up arrow if requested
  if (arrow_length > 0) {
    arrow_obj <- ggplot2::arrow(
      type = arrow_type,
      length = ggplot2::unit(arrow_length, "inches")
    )
  } else {
    arrow_obj <- NULL
  }
  
  # Create the curve geom
  ggplot2::geom_curve(
    mapping = mapping, data = data, stat = stat,
    position = position, curvature = curvature, arrow = arrow_obj,
    color = color, size = size, alpha = alpha, ...,
    show.legend = show.legend, inherit.aes = inherit.aes
  )
}

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