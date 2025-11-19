#' Geom for genomic interaction links
#'
#' This function creates a geom for genomic interaction links, such as chromatin loops,
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
#' @importFrom methods is
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' df <- data.frame(
#'     start1 = c(1000, 2000, 3000),
#'     end1 = c(1100, 2100, 3100),
#'     start2 = c(5000, 6000, 7000),
#'     end2 = c(5100, 6100, 7100),
#'     score = c(0.8, 0.6, 0.9)
#' )
#'
#' # Basic link plot
#' ggplot(df) +
#'     geom_link(aes(x = start1, y = 0, xend = start2, yend = 0))
#'
#' # With curvature and arrows
#' ggplot(df) +
#'     geom_link(aes(x = start1, y = 0, xend = start2, yend = 0, color = score),
#'         curvature = 0.8, arrow_length = 0.1
#'     )
#' }
geom_link <- function(mapping = NULL, data = NULL, stat = "identity",
                      position = "identity", ..., curvature = 0.5,
                      arrow_length = 0, arrow_type = "closed",
                      na.rm = TRUE, show.legend = NA, inherit.aes = TRUE) {
    # Handle GInteractions objects directly in data
    if (!is.null(data) && methods::is(data, "GInteractions")) {
        data <- process_interaction_data(data)
    }

    # Default mapping for links if not provided
    # We can't easily guess column names if mapping is NULL and data is a generic DF,
    # but if it came from process_interaction_data, it has start1/start2.
    # However, the user requirement says input should be x, x_end, y, y_end.
    # We will leave mapping as NULL if not provided, but we can provide a hint or default
    # if the data looks like it has standard columns.

    # If data came from our helper, it has start1, start2.
    if (is.null(mapping) && !is.null(data) && all(c("start1", "start2") %in% names(data))) {
        default_aes <- aes(x = .data$start1, y = 0, xend = .data$start2, yend = 0)
        mapping <- default_aes
    }

    layer(
        data = data,
        mapping = mapping,
        stat = stat,
        geom = GeomLink,
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

#' @rdname geom_link
#' @format NULL
#' @usage NULL
GeomLink <- ggproto("GeomLink", Geom,
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
#' This function processes interaction data for visualization with geom_link.
#' It formats the data for use with ggplot2.
#'
#' @param gr A GRanges or GInteractions object with interaction data
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
        # We need to check if GInteractions class is available/loaded, but methods::is should handle it safely if object is passed

        # Accessors for GInteractions might need package loading or direct slot access if S4
        # Usually anchor1(gr) and anchor2(gr) return GRanges
        # But to be safe and avoid extra deps if not loaded, we try to use S4 accessors if available or slots

        if (requireNamespace("InteractionSet", quietly = TRUE)) {
            a1 <- InteractionSet::anchors(gr, type = "first")
            a2 <- InteractionSet::anchors(gr, type = "second")

            result <- data.frame(
                start1 = GenomicRanges::start(a1),
                end1 = GenomicRanges::end(a1),
                start2 = GenomicRanges::start(a2),
                end2 = GenomicRanges::end(a2)
            )

            if (score %in% colnames(S4Vectors::mcols(gr))) {
                result$score <- S4Vectors::mcols(gr)[[score]]
            }
            return(result)
        } else {
            # Fallback if InteractionSet is not loaded but object is GInteractions (unlikely but possible)
            # Try direct slot access if possible, or error
            stop("InteractionSet package is required to process GInteractions objects")
        }
    } else {
        stop("Input must be a GRanges or GInteractions object")
  }
}

#' Create an interaction track from a file
#'
#' This function reads interaction data from a file and creates a track.
#'
#' @param file Path to the interaction file (e.g., BEDPE)
#' @param region Genomic region to display
#' @param ... Additional arguments passed to geom_link
#' @return A ggplot2 object
#' @export
#' @importFrom rtracklayer import
#' @importFrom ggplot2 ggplot aes
interaction_track <- function(file, region, ...) {
  # Parse region
  region_gr <- parse_region(region)
  
  # Import data
  # Note: rtracklayer::import might not support BEDPE natively with 'which' for all formats
  # But for standard BED/BEDPE it might work or we import all and filter
  # For simplicity, we try to import with 'which' if supported, otherwise import all
  
  tryCatch({
    gr <- rtracklayer::import(file, which = region_gr)
  }, error = function(e) {
    # Fallback: import all and filter
    gr <- rtracklayer::import(file)
    # TODO: Filter by region manually if needed, but process_interaction_data might handle it
    # or we rely on the plot limits.
    # For interactions, filtering is tricky (both anchors in region? or any?)
    # For now, we assume the user provides a file that rtracklayer handles or small enough.
    # Ideally we should filter here.
    gr <- IRanges::subsetByOverlaps(gr, region_gr)
  })
  
  # Process data
  df <- process_interaction_data(gr)
  
  # Create plot using ez_link logic (but we can't call ez_link because ez_link calls us)
  # So we replicate the plotting logic or make ez_link handle the plotting and this function just return data?
  # But ez_link expects this function to return a plot.
  
  # Let's make interaction_track return the plot.
  # We can reuse ez_link by passing the dataframe!
  # But ez_link calls interaction_track if input is char.
  # If we pass df to ez_link, it will use the dataframe logic.
  
  return(ez_link(df, region, ...))
}
```
