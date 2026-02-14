#' Geom for Hi-C contact matrix visualization (square heatmap)
#'
#' This function creates a geom for Hi-C contact matrix visualization. It displays
#' the contact matrix as a square heatmap where x and y represent genomic positions
#' and fill represents contact frequency.
#'
#' @inheritParams ggplot2::layer
#' @param resolution Resolution of Hi-C bins in base pairs. Used to calculate tile width/height.
#'   If NULL, attempts to infer from data (default: NULL)
#' @param na.rm If `TRUE`, silently remove missing values.
#' @param rasterize Logical. If TRUE and ggrastr package is available, rasterize
#'   the layer for better performance with large matrices. Default: FALSE
#' @param rasterize_dpi Resolution for rasterization in dots per inch. Default: 300
#' @param ... Additional arguments passed to [ggplot2::layer()]
#'
#' @return A ggplot2 layer that can be added to a plot.
#' @export
#' @importFrom ggplot2 layer aes ggproto Geom
#' @importFrom grid rectGrob gpar gTree gList
#' @importFrom methods is
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' # Create example Hi-C data
#' hic_df <- data.frame(
#'   x = rep(seq(1000000, 1500000, by = 10000), each = 51),
#'   y = rep(seq(1000000, 1500000, by = 10000), times = 51),
#'   score = runif(51 * 51)
#' )
#' ggplot(hic_df, aes(x = x, y = y, fill = score)) +
#'   geom_hic(resolution = 10000) +
#'   scale_fill_hic()
#'
#' # With rasterization for large matrices
#' ggplot(hic_df, aes(x = x, y = y, fill = score)) +
#'   geom_hic(resolution = 10000, rasterize = TRUE, rasterize_dpi = 300) +
#'   scale_fill_hic()
#' }
geom_hic <- function(
  mapping = NULL,
  data = NULL,
  stat = "identity",
  position = "identity",
  ...,
  resolution = NULL,
  na.rm = FALSE,
  rasterize = FALSE,
  rasterize_dpi = 300,
  show.legend = NA,
  inherit.aes = TRUE
) {
  layer_obj <- layer_obj <- ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomHic,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      resolution = resolution,
      na.rm = na.rm,
      ...
    )
  )

  # Apply rasterization if requested
  if (rasterize) {
    if (!requireNamespace("ggrastr", quietly = TRUE)) {
      warning(
        "Package 'ggrastr' is required for rasterization. ",
        "Install with install.packages('ggrastr')"
      )
      return(layer_obj)
    }
    layer_obj <- ggrastr::rasterise(layer_obj, dpi = rasterize_dpi)
  }

  layer_obj
}

#' @rdname geom_hic
#' @format NULL
#' @usage NULL
#' @export
GeomHic <- ggproto(
  "GeomHic",
  ggplot2::Geom,
  required_aes = c("x", "y", "fill"),
  default_aes = ggplot2::aes(
    alpha = 1,
    colour = NA
  ),
  draw_key = ggplot2::draw_key_polygon,
  setup_data = function(data, params) {
    # Calculate tile width/height from resolution or infer from data
    if (!is.null(params$resolution)) {
      data$width <- params$resolution
      data$height <- params$resolution
    } else {
      # Try to infer resolution from data
      x_vals <- sort(unique(data$x))
      if (length(x_vals) > 1) {
        resolution <- min(diff(x_vals))
        data$width <- resolution
        data$height <- resolution
      } else {
        # Default to 1 if can't infer

        data$width <- 1
        data$height <- 1
      }
    }
    data
  },
  draw_panel = function(
    data,
    panel_params,
    coord,
    resolution = NULL,
    na.rm = FALSE
  ) {
    if (nrow(data) == 0) {
      return(grid::nullGrob())
    }

    # Create rectangles for each Hi-C bin
    # Calculate corners of each tile
    data$xmin <- data$x - data$width / 2
    data$xmax <- data$x + data$width / 2
    data$ymin <- data$y - data$height / 2
    data$ymax <- data$y + data$height / 2

    # Transform coordinates
    coords <- coord$transform(data, panel_params)

    # Draw rectangles
    grid::rectGrob(
      x = coords$x,
      y = coords$y,
      width = coords$xmax - coords$xmin,
      height = coords$ymax - coords$ymin,
      default.units = "native",
      just = "center",
      gp = grid::gpar(
        col = coords$colour,
        fill = coords$fill,
        alpha = coords$alpha
      )
    )
  }
)

#' Geom for Hi-C contact matrix visualization (triangle/rotated view)
#'
#' This function creates a geom for Hi-C contact matrix visualization in a triangle
#' (rotated 45-degree) view. This is the standard view for genome browser tracks where
#' the x-axis shows genomic position and the y-axis shows interaction distance.
#'
#' @inheritParams ggplot2::layer
#' @param resolution Resolution of Hi-C bins in base pairs. Required for proper diamond sizing.
#'   If NULL, attempts to infer from data (default: NULL)
#' @param max_distance Maximum interaction distance to display in base pairs (default: NULL, show all)
#' @param na.rm If `TRUE`, silently remove missing values.
#' @param rasterize Logical. If TRUE and ggrastr package is available, rasterize
#'   the layer for better performance with large matrices. Default: FALSE
#' @param rasterize_dpi Resolution for rasterization in dots per inch. Default: 300
#' @param ... Additional arguments passed to [ggplot2::layer()]
#'
#' @details
#' The triangle view transforms Hi-C data so that:
#' \itemize{
#'   \item x-axis represents the midpoint of interacting bins: (pos1 + pos2) / 2
#'   \item y-axis represents the interaction distance: pos2 - pos1
#' }
#' The triangle points UPWARD with the diagonal (self-interactions) at y=0 at the
#' bottom, and the apex at the top showing maximum interaction distance. Each contact
#' is rendered as a diamond (rotated square). Use with `coord_fixed(ratio = 1)` for
#' a proper half-square aspect ratio where the triangle forms 45-degree edges.
#'
#' @return A ggplot2 layer that can be added to a plot.
#' @export
#' @importFrom ggplot2 layer aes ggproto Geom
#' @importFrom grid polygonGrob gpar gTree gList
#' @importFrom methods is
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' # Create example Hi-C data
#' hic_df <- expand.grid(
#'   pos1 = seq(1000000, 1500000, by = 10000),
#'   pos2 = seq(1000000, 1500000, by = 10000)
#' )
#' hic_df <- hic_df[hic_df$pos2 >= hic_df$pos1, ]  # Upper triangle only
#' hic_df$score <- runif(nrow(hic_df))
#'
#' ggplot(hic_df, aes(x = pos1, y = pos2, fill = score)) +
#'   geom_hic_triangle(resolution = 10000) +
#'   scale_fill_hic() +
#'   coord_fixed(ratio = 1)
#'
#' # With rasterization for large matrices
#' ggplot(hic_df, aes(x = pos1, y = pos2, fill = score)) +
#'   geom_hic_triangle(resolution = 10000, rasterize = TRUE) +
#'   scale_fill_hic() +
#'   coord_fixed(ratio = 1)
#' }
geom_hic_triangle <- function(
  mapping = NULL,
  data = NULL,
  stat = "identity",
  position = "identity",
  ...,
  resolution = NULL,
  max_distance = NULL,
  na.rm = FALSE,
  rasterize = FALSE,
  rasterize_dpi = 300,
  show.legend = NA,
  inherit.aes = TRUE
) {
  layer_obj <- layer_obj <- ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomHicTriangle,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      resolution = resolution,
      max_distance = max_distance,
      na.rm = na.rm,
      ...
    )
  )

  # Apply rasterization if requested
  if (rasterize) {
    if (!requireNamespace("ggrastr", quietly = TRUE)) {
      warning(
        "Package 'ggrastr' is required for rasterization. ",
        "Install with install.packages('ggrastr')"
      )
      return(layer_obj)
    }
    layer_obj <- ggrastr::rasterise(layer_obj, dpi = rasterize_dpi)
  }

  layer_obj
}

#' @rdname geom_hic_triangle
#' @format NULL
#' @usage NULL
#' @export
GeomHicTriangle <- ggproto(
  "GeomHicTriangle",
  ggplot2::Geom,
  required_aes = c("x", "y", "fill"),
  default_aes = ggplot2::aes(
    alpha = 1,
    colour = NA
  ),
  draw_key = ggplot2::draw_key_polygon,
  setup_data = function(data, params) {
    # Ensure we have pos1 and pos2 (or x/y representing bin positions)
    # The user should provide x = pos1, y = pos2 where pos1 <= pos2

    # Calculate resolution if not provided
    if (!is.null(params$resolution)) {
      resolution <- params$resolution
    } else {
      x_vals <- sort(unique(data$x))
      if (length(x_vals) > 1) {
        resolution <- min(diff(x_vals))
      } else {
        resolution <- 10000 # Default fallback
        warning("Could not infer resolution, using default 10000bp")
      }
    }
    data$resolution <- resolution

    # Filter by max_distance if specified
    if (!is.null(params$max_distance)) {
      distance <- data$y - data$x
      data <- data[distance <= params$max_distance, ]
    }

    # Transform to triangle coordinates
    # x_new = midpoint = (pos1 + pos2) / 2
    # y_new = distance = pos2 - pos1 (positive, triangle points UPWARD)
    data$x_orig <- data$x
    data$y_orig <- data$y
    data$x <- (data$x_orig + data$y_orig) / 2
    data$y <- data$y_orig - data$x_orig # Positive distance, apex at top

    data
  },
  draw_panel = function(
    data,
    panel_params,
    coord,
    resolution = NULL,
    max_distance = NULL,
    na.rm = FALSE
  ) {
    if (nrow(data) == 0) {
      return(grid::nullGrob())
    }

    resolution <- data$resolution[1]

    # Create diamond polygons for each Hi-C contact
    # Diamond vertices (in data coordinates, before transform):
    # Each diamond has width = resolution (in x) and height = resolution (in y distance)
    # Center at (x, y) where x is midpoint and y is distance
    #
    # Vertices of diamond:
    #   top:    (x, y + resolution/2)
    #   right:  (x + resolution/2, y)
    #   bottom: (x, y - resolution/2)  -- but y >= 0, so we handle edge cases

    #   left:   (x - resolution/2, y)

    # Build polygon data for all diamonds
    n <- nrow(data)

    # Diamond half-size
    half_res <- resolution / 2

    # For each point, create 4 vertices
    # Note: For points on the diagonal (distance = 0), we draw triangles instead
    # Triangle points UPWARD (y is positive distance, apex at top)
    grobs <- lapply(seq_len(n), function(i) {
      row <- data[i, ]
      cx <- row$x
      cy <- row$y

      if (cy == 0) {
        # On diagonal: draw triangle (top half of diamond, pointing up)
        xs <- c(cx - half_res, cx, cx + half_res)
        ys <- c(cy, cy + half_res, cy)
      } else {
        # Off diagonal: draw full diamond (y is positive)
        xs <- c(cx, cx + half_res, cx, cx - half_res)
        ys <- c(cy + half_res, cy, max(0, cy - half_res), cy)
      }

      # Create temporary data frame for coordinate transformation
      df <- data.frame(x = xs, y = ys)
      coords_transformed <- coord$transform(df, panel_params)

      list(
        x = coords_transformed$x,
        y = coords_transformed$y,
        fill = row$fill,
        alpha = row$alpha,
        colour = row$colour
      )
    })

    # Combine into polygon grobs
    poly_grobs <- lapply(grobs, function(g) {
      grid::polygonGrob(
        x = g$x,
        y = g$y,
        default.units = "native",
        gp = grid::gpar(
          fill = g$fill,
          col = g$colour,
          alpha = g$alpha
        )
      )
    })

    grid::gTree(children = do.call(grid::gList, poly_grobs))
  }
)

#' Process Hi-C data for visualization
#'
#' This function processes Hi-C data from various input formats for visualization
#' with geom_hic or geom_hic_triangle. It handles dense matrices, sparse data frames,
#' and file paths.
#'
#' @param data Input data. Can be:
#'   \itemize{
#'     \item A matrix: Dense contact matrix where rows and columns represent bins
#'     \item A data frame with columns (bin1, bin2, score) or (pos1, pos2, score): Sparse format
#'     \item A file path: Tab-delimited matrix file with row/column headers
#'   }
#' @param region Genomic region to display (e.g., "chr1:1000000-2000000"). Required for
#'   file input, optional for data frame/matrix if coordinates are already genomic.
#' @param resolution Resolution of the Hi-C data in base pairs (default: 10000).
#'   Used to convert bin indices to genomic coordinates for matrix input.
#' @param upper_triangle Logical. If TRUE, only return upper triangle (pos1 <= pos2).
#'   Useful for triangle visualization. Default: FALSE
#' @param symmetric Logical. If TRUE and data is a matrix, assume it's symmetric and

#'   extract upper triangle. Default: TRUE
#'
#' @return A data frame with columns:
#'   \itemize{
#'     \item pos1: Genomic position of first bin
#'     \item pos2: Genomic position of second bin
#'     \item score: Contact frequency/count
#'   }
#' @export
#' @importFrom GenomicRanges seqnames start end
#' @examples
#' \dontrun{
#' # From a matrix
#' mat <- matrix(runif(100), nrow = 10)
#' hic_df <- process_hic_data(mat, "chr1:1000000-1100000", resolution = 10000)
#'
#' # From a sparse data frame
#' sparse_df <- data.frame(pos1 = c(1e6, 1e6), pos2 = c(1e6, 1.01e6), score = c(100, 50))
#' hic_df <- process_hic_data(sparse_df)
#'
#' # From a file
#' hic_df <- process_hic_data("contacts.matrix", "chr1:1000000-2000000")
#' }
process_hic_data <- function(
  data,
  region = NULL,
  resolution = 10000,
  upper_triangle = FALSE,
  symmetric = TRUE
) {
  # Case 1: File path
  if (is.character(data) && length(data) == 1 && file.exists(data)) {
    if (is.null(region)) {
      stop("region is required when reading from a file")
    }

    # Parse the region
    region_gr <- parse_region(region)
    start_pos <- GenomicRanges::start(region_gr)
    end_pos <- GenomicRanges::end(region_gr)

    # Calculate bin indices
    start_bin <- floor(start_pos / resolution)
    end_bin <- ceiling(end_pos / resolution)
    n_bins <- end_bin - start_bin + 1

    # Read the matrix
    hic_matrix <- as.matrix(utils::read.table(
      data,
      header = TRUE,
      row.names = 1
    ))

    # Extract the region of interest
    bin_indices <- seq_len(n_bins)
    if (n_bins > nrow(hic_matrix) || n_bins > ncol(hic_matrix)) {
      stop("Region is outside the bounds of the Hi-C matrix")
    }
    hic_subset <- hic_matrix[bin_indices, bin_indices, drop = FALSE]

    # Convert to data frame
    hic_df <- expand.grid(bin1 = seq_len(n_bins), bin2 = seq_len(n_bins))
    hic_df$score <- as.vector(hic_subset)
    hic_df$pos1 <- (start_bin + hic_df$bin1 - 1) * resolution + resolution / 2
    hic_df$pos2 <- (start_bin + hic_df$bin2 - 1) * resolution + resolution / 2
    hic_df <- hic_df[, c("pos1", "pos2", "score")]
  } else if (is.matrix(data)) {
    # Case 2: Matrix
    n_bins <- nrow(data)

    # Parse region if provided for coordinate conversion
    if (!is.null(region)) {
      region_gr <- parse_region(region)
      start_pos <- GenomicRanges::start(region_gr)
    } else {
      start_pos <- 0
    }

    start_bin <- floor(start_pos / resolution)

    # Convert matrix to data frame
    hic_df <- expand.grid(bin1 = seq_len(n_bins), bin2 = seq_len(n_bins))
    hic_df$score <- as.vector(data)
    hic_df$pos1 <- (start_bin + hic_df$bin1 - 1) * resolution + resolution / 2
    hic_df$pos2 <- (start_bin + hic_df$bin2 - 1) * resolution + resolution / 2
    hic_df <- hic_df[, c("pos1", "pos2", "score")]
  } else if (is.data.frame(data)) {
    # Case 3: Data frame (sparse format)
    # Check which columns are present
    if (all(c("pos1", "pos2") %in% names(data))) {
      # Already has genomic coordinates
      hic_df <- data
      if (!"score" %in% names(hic_df) && "count" %in% names(hic_df)) {
        hic_df$score <- hic_df$count
      }

      # Filter to region if specified
      if (!is.null(region)) {
        region_gr <- parse_region(region)
        start_pos <- GenomicRanges::start(region_gr)
        end_pos <- GenomicRanges::end(region_gr)
        hic_df <- hic_df[
          hic_df$pos1 >= start_pos &
            hic_df$pos1 <= end_pos &
            hic_df$pos2 >= start_pos &
            hic_df$pos2 <= end_pos,
        ]
      }
    } else if (all(c("bin1", "bin2") %in% names(data))) {
      # Has bin indices, convert to genomic coordinates
      if (!is.null(region)) {
        region_gr <- parse_region(region)
        start_pos <- GenomicRanges::start(region_gr)
        start_bin <- floor(start_pos / resolution)
      } else {
        start_bin <- 0
      }

      hic_df <- data
      hic_df$pos1 <- (start_bin + hic_df$bin1 - 1) * resolution + resolution / 2
      hic_df$pos2 <- (start_bin + hic_df$bin2 - 1) * resolution + resolution / 2

      if (!"score" %in% names(hic_df) && "count" %in% names(hic_df)) {
        hic_df$score <- hic_df$count
      }
    } else {
      stop("Data frame must have columns (pos1, pos2) or (bin1, bin2)")
    }

    # Ensure we have required columns
    if (!"score" %in% names(hic_df)) {
      stop("Data frame must have a 'score' or 'count' column")
    }

    hic_df <- hic_df[, c("pos1", "pos2", "score")]
  } else {
    stop("data must be a file path, matrix, or data frame")
  }

  # Complete the grid: fill missing bin pairs with score = 0
  # This ensures all bins in the plotting region are represented
  all_positions <- sort(unique(c(hic_df$pos1, hic_df$pos2)))
  complete_grid <- expand.grid(pos1 = all_positions, pos2 = all_positions)
  hic_df <- merge(complete_grid, hic_df, by = c("pos1", "pos2"), all.x = TRUE)
  hic_df$score[is.na(hic_df$score)] <- 0

  # Filter to upper triangle if requested
  if (upper_triangle) {
    hic_df <- hic_df[hic_df$pos1 <= hic_df$pos2, ]
  }

  return(hic_df)
}
