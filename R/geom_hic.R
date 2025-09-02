#' Geom for Hi-C contact matrix visualization
#'
#' This function creates a geom for Hi-C contact matrix visualization. It displays
#' the contact matrix as a heatmap.
#'
#' @inheritParams ggplot2::geom_tile
#' @param low Color for low values (default: "white")
#' @param high Color for high values (default: "red")
#' @param midpoint Midpoint for diverging color scales (default: NULL)
#' @param transform Function to transform the values (default: NULL)
#' @param alpha Transparency (default: 1)
#' @return A ggplot2 layer
#' @export
#' @importFrom ggplot2 geom_tile scale_fill_gradient scale_fill_gradient2 aes
#' @examples
#' \dontrun{
#' library(ggplot2)
#' p <- ggplot(hic_data, aes(x = bin1, y = bin2, fill = count)) + geom_hic()
#' }
geom_hic <- function(mapping = NULL, data = NULL, stat = "identity",
                      position = "identity", ..., low = "white", high = "red",
                      midpoint = NULL, transform = NULL, alpha = 1,
                      show.legend = NA, inherit.aes = TRUE) {
  
  # Create a list to hold the layers
  layers <- list()
  
  # Add the tile geom
  layers[[1]] <- ggplot2::geom_tile(
    mapping = mapping, data = data, stat = stat,
    position = position, alpha = alpha, ...,
    show.legend = show.legend, inherit.aes = inherit.aes
  )
  
  # Add the appropriate scale
  if (is.null(midpoint)) {
    layers[[2]] <- ggplot2::scale_fill_gradient(low = low, high = high, transform = transform)
  } else {
    layers[[2]] <- ggplot2::scale_fill_gradient2(low = low, mid = "white", high = high,
                                                 midpoint = midpoint, transform = transform)
  }
  
  # Return the layers
  return(layers)
}

#' Process Hi-C data for visualization
#'
#' This function processes Hi-C data for visualization with geom_hic.
#' It formats the data for use with ggplot2.
#'
#' @param file Path to the Hi-C contact matrix file
#' @param region Genomic region to display (e.g., "chr1:1000000-2000000")
#' @param resolution Resolution of the Hi-C data in base pairs (default: 10000)
#' @param log_transform Apply log transformation to the values (default: TRUE)
#' @param pseudo_count Pseudo-count to add before log transformation (default: 1)
#' @return A data frame with Hi-C contact information
#' @export
#' @examples
#' \dontrun{
#' hic_data <- process_hic_data("contacts.matrix", "chr1:1000000-2000000", resolution = 10000)
#' }
process_hic_data <- function(file, region, resolution = 10000, log_transform = TRUE,
                             pseudo_count = 1) {
  # Parse the region
  region_gr <- parse_region(region)
  chr <- as.character(GenomicRanges::seqnames(region_gr))
  start_pos <- GenomicRanges::start(region_gr)
  end_pos <- GenomicRanges::end(region_gr)
  
  # Calculate bin indices
  start_bin <- floor(start_pos / resolution)
  end_bin <- ceiling(end_pos / resolution)
  n_bins <- end_bin - start_bin + 1
  
  # Read the Hi-C matrix
  # This is a simplified version that assumes a specific format
  # In practice, you would need to handle different Hi-C file formats
  if (file.exists(file)) {
    # Read the matrix (assuming a tab-delimited format with row and column headers)
    hic_matrix <- as.matrix(read.table(file, header = TRUE, row.names = 1))
    
    # Extract the region of interest
    bin_indices <- (start_bin:end_bin) - min(start_bin) + 1
    if (max(bin_indices) <= nrow(hic_matrix) && max(bin_indices) <= ncol(hic_matrix)) {
      hic_subset <- hic_matrix[bin_indices, bin_indices]
    } else {
      stop("Region is outside the bounds of the Hi-C matrix")
    }
  } else {
    # If the file doesn't exist, create a dummy matrix for demonstration
    warning("File not found. Creating a dummy Hi-C matrix for demonstration.")
    hic_subset <- matrix(runif(n_bins^2), nrow = n_bins, ncol = n_bins)
    hic_subset <- hic_subset * t(hic_subset)  # Make it symmetric
  }
  
  # Convert the matrix to a data frame for ggplot2
  hic_df <- expand.grid(bin1 = 1:n_bins, bin2 = 1:n_bins)
  hic_df$count <- as.vector(hic_subset)
  
  # Convert bin indices to genomic coordinates
  hic_df$pos1 <- (start_bin + hic_df$bin1 - 1) * resolution
  hic_df$pos2 <- (start_bin + hic_df$bin2 - 1) * resolution
  
  # Apply log transformation if requested
  if (log_transform) {
    hic_df$count <- log10(hic_df$count + pseudo_count)
  }
  
  return(hic_df)
}

#' Create a Hi-C track from a contact matrix file
#'
#' This function creates a Hi-C track from a contact matrix file. It processes the data
#' for a specific region and creates a ggplot2 layer for visualization.
#'
#' @param file Path to the Hi-C contact matrix file
#' @param region Genomic region to display (e.g., "chr1:1000000-2000000")
#' @param resolution Resolution of the Hi-C data in base pairs (default: 10000)
#' @param log_transform Apply log transformation to the values (default: TRUE)
#' @param low Color for low values (default: "white")
#' @param high Color for high values (default: "red")
#' @param ... Additional arguments passed to geom_hic
#' @return A ggplot2 layer
#' @export
#' @importFrom ggplot2 ggplot aes coord_fixed
#' @examples
#' \dontrun{
#' p <- hic_track("contacts.matrix", "chr1:1000000-2000000", resolution = 10000)
#' }
hic_track <- function(file, region, resolution = 10000, log_transform = TRUE,
                       low = "white", high = "red", ...) {
  # Process the Hi-C data
  hic_data <- process_hic_data(file, region, resolution, log_transform)
  
  # Create the plot
  p <- ggplot2::ggplot(hic_data, ggplot2::aes(x = pos1, y = pos2, fill = count)) +
    geom_hic(low = low, high = high, ...) +
    ggplot2::coord_fixed()  # Ensure the plot is square
  
  # Apply the appropriate theme and scale
  p <- p + ez_theme() +
    scale_x_genome_region(region) +
    scale_x_genome_region(region)  # Same scale for y-axis
  
  return(p)
}