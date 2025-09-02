# Tests for scales.R functions

test_that("format_genomic_coord formats genomic coordinates correctly", {
  # Test various coordinate values
  expect_equal(format_genomic_coord(0), "0")
  expect_equal(format_genomic_coord(10), "10")
  expect_equal(format_genomic_coord(1000), "1kb")
  expect_equal(format_genomic_coord(1500), "1.5kb")
  expect_equal(format_genomic_coord(10000), "10kb")
  expect_equal(format_genomic_coord(1000000), "1Mb")
  expect_equal(format_genomic_coord(1500000), "1.5Mb")
  expect_equal(format_genomic_coord(1000000000), "1Gb")
  
  # Test with custom units
  expect_equal(format_genomic_coord(1000, units = c("b", "kb", "Mb", "Gb")), "1kb")
  expect_equal(format_genomic_coord(1000000, units = c("b", "kb", "Mb", "Gb")), "1Mb")
  
  # Test with custom breaks
  expect_equal(format_genomic_coord(500, breaks = c(0, 100, 1000, 1000000)), "500")
  expect_equal(format_genomic_coord(5000, breaks = c(0, 100, 1000, 1000000)), "5kb")
})

test_that("scale_x_genome creates a proper scale", {
  # Create a simple ggplot object
  p <- ggplot2::ggplot(data.frame(x = c(1, 1000000), y = c(1, 2)), ggplot2::aes(x = x, y = y))
  
  # Add scale_x_genome
  p_with_scale <- p + scale_x_genome()
  
  # Extract the scales from the plot
  scales <- p_with_scale$scales$scales
  x_scale <- scales[[1]]
  
  # Check that it's a continuous scale
  expect_s3_class(x_scale, "ScaleContinuousPosition")
  
  # Check that it has the correct name
  expect_equal(x_scale$name, "Genomic Position")
})

test_that("scale_x_genome_region creates a proper scale for a specific region", {
  # Create a simple ggplot object
  p <- ggplot2::ggplot(data.frame(x = c(1000, 2000), y = c(1, 2)), ggplot2::aes(x = x, y = y))
  
  # Add scale_x_genome_region
  p_with_scale <- p + scale_x_genome_region("chr1:1000-2000")
  
  # Extract the scales from the plot
  scales <- p_with_scale$scales$scales
  x_scale <- scales[[1]]
  
  # Check that it's a continuous scale
  expect_s3_class(x_scale, "ScaleContinuousPosition")
  
  # Check that it has the correct name
  expect_equal(x_scale$name, "chr1:1,000-2,000")
  
  # Check that it has the correct limits
  expect_equal(x_scale$limits, c(1000, 2000))
})