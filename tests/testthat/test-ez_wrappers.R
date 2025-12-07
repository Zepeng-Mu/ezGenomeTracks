# Tests for ez_wrappers.R functions

test_that("ez_coverage handles various input types correctly", {
  data(example_signal)

  # Test with data frame
  p1 <- ez_coverage(example_signal, region = "chr1:1000000-2000000")
  expect_s3_class(p1, "ggplot")
  expect_true(length(p1$layers) > 0)

  # Test different types
  p2 <- ez_coverage(example_signal, region = "chr1:1000000-2000000", type = "line")
  expect_s3_class(p2, "ggplot")

  # Test parameter validation
  expect_error(ez_coverage(example_signal, region = "chr1:1000000-2000000", alpha = 2))
  expect_error(ez_coverage(example_signal, region = "chr1:1000000-2000000", type = "invalid"))

  # Test with missing required columns
  bad_data <- example_signal[, -which(names(example_signal) == "score")]
  expect_error(ez_coverage(bad_data, region = "chr1:1000000-2000000"))
})

test_that("ez_feature creates a feature track", {
  # Load example data
  data(example_peaks)

  # Create a feature track
  p <- ez_feature(data = example_peaks, region = "chr1:1000000-2000000")

  # Check that it returns a ggplot object
  expect_s3_class(p, "ggplot")

  # Check that it has at least one layer
  expect_true(length(p$layers) > 0)
})

test_that("ez_gene creates a gene track from data frame", {
  # Load example data
  data(example_genes)

  # Create a gene track
  p <- ez_gene(data = example_genes, region = "chr1:1000000-2000000")

  # Check that it returns a ggplot object
  expect_s3_class(p, "ggplot")

  # Check that it has at least one layer
  expect_true(length(p$layers) > 0)
})

test_that("ez_gene handles TxDb objects", {
  # Skip if TxDb.Hsapiens.UCSC.hg19.knownGene is not available
  skip_if_not_installed("TxDb.Hsapiens.UCSC.hg19.knownGene")

  # Load TxDb package
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

  # Test that ez_gene accepts TxDb objects
  expect_no_error({
    p <- ez_gene(data = txdb, region = "chr1:1000000-2000000")
  })

  # Create a gene track from TxDb
  p <- ez_gene(data = txdb, region = "chr1:1000000-2000000")

  # Check that it returns a ggplot object
  expect_s3_class(p, "ggplot")

  # Check that it has at least one layer
  expect_true(length(p$layers) > 0)
})

test_that("ez_link creates an interaction track", {
  # Load example data
  data(example_interactions)

  # Create an interaction track
  p <- ez_link(data = example_interactions, region = "chr1:1000000-1100000")

  # Check that it returns a ggplot object
  expect_s3_class(p, "ggplot")

  # Check that it has at least one layer
  expect_true(length(p$layers) > 0)
})

test_that("ez_hic creates a Hi-C track", {
  # Load example data
  data(example_hic)

  # Create a Hi-C track
  p <- ez_hic(data = example_hic, region = "chr1:1000000-1100000")

  # Check that it returns a ggplot object
  expect_s3_class(p, "ggplot")

  # Check that it has at least one layer
  expect_true(length(p$layers) > 0)
})

test_that("ez_sashimi creates a sashimi plot", {
  # Create test coverage data
  coverage_data <- data.frame(
    seqnames = "chr1",
    start = seq(1000, 5000, by = 50),
    end = seq(1050, 5050, by = 50),
    score = c(runif(20, 5, 10), rep(0, 20), runif(20, 5, 10), rep(0, 20), runif(20, 5, 10))
  )

  # Create test junction data
  junction_data <- data.frame(
    seqnames = "chr1",
    start = c(1500, 2500, 3500),
    end = c(2000, 3000, 4000),
    score = c(25, 50, 15)
  )

  # Basic sashimi plot
  p1 <- ez_sashimi(coverage_data, junction_data, region = "chr1:1000-5000")
  expect_s3_class(p1, "ggplot")
  expect_true(length(p1$layers) > 0)

  # Test with different junction directions
  p2 <- ez_sashimi(
    coverage_data, junction_data, region = "chr1:1000-5000",
    junction_direction = "up"
  )
  expect_s3_class(p2, "ggplot")

  p3 <- ez_sashimi(
    coverage_data, junction_data, region = "chr1:1000-5000",
    junction_direction = "down"
  )
  expect_s3_class(p3, "ggplot")

  # Test with score transformations
  p4 <- ez_sashimi(
    coverage_data, junction_data, region = "chr1:1000-5000",
    score_transform = "log10"
  )
  expect_s3_class(p4, "ggplot")

  p5 <- ez_sashimi(
    coverage_data, junction_data, region = "chr1:1000-5000",
    score_transform = "sqrt"
  )
  expect_s3_class(p5, "ggplot")

  # Test with labels disabled
  p6 <- ez_sashimi(
    coverage_data, junction_data, region = "chr1:1000-5000",
    show_labels = FALSE
  )
  expect_s3_class(p6, "ggplot")

  # Test with custom colors and styling
  p7 <- ez_sashimi(
    coverage_data, junction_data, region = "chr1:1000-5000",
    coverage_fill = "darkblue",
    junction_color = "red",
    linewidth_range = c(1, 5)
  )
  expect_s3_class(p7, "ggplot")

  # Test with y-axis styles
  p8 <- ez_sashimi(
    coverage_data, junction_data, region = "chr1:1000-5000",
    y_axis_style = "full"
  )
  expect_s3_class(p8, "ggplot")

  # Test parameter validation
  expect_error(ez_sashimi(
    coverage_data, junction_data, region = "chr1:1000-5000",
    alpha = 2
  ))
  expect_error(ez_sashimi(
    coverage_data, junction_data, region = "chr1:1000-5000",
    junction_direction = "invalid"
  ))
})

test_that("ez_sashimi handles empty junction data", {
  # Create test coverage data
  coverage_data <- data.frame(
    seqnames = "chr1",
    start = seq(1000, 2000, by = 50),
    end = seq(1050, 2050, by = 50),
    score = runif(20, 5, 10)
  )

  # Empty junction data
  junction_data <- data.frame(
    seqnames = character(0),
    start = numeric(0),
    end = numeric(0),
    score = numeric(0)
  )

  # Should work with just coverage
  p <- ez_sashimi(coverage_data, junction_data, region = "chr1:1000-2000")
  expect_s3_class(p, "ggplot")
})

test_that("process_sashimi_data processes inputs correctly", {
  # Create test data
  coverage_data <- data.frame(
    seqnames = "chr1",
    start = seq(1000, 2000, by = 100),
    end = seq(1100, 2100, by = 100),
    score = c(5, 8, 10, 12, 15, 12, 10, 8, 5, 3, 2)
  )

  junction_data <- data.frame(
    seqnames = "chr1",
    start = c(1100, 1500, 1800),
    end = c(1400, 1700, 1950),
    score = c(10, 20, 5)
  )

  # Test alternate direction
  result <- process_sashimi_data(
    coverage_data, junction_data, "chr1:1000-2000",
    junction_direction = "alternate"
  )

  expect_type(result, "list")
  expect_true("coverage_df" %in% names(result))
  expect_true("junction_df" %in% names(result))
  expect_true("arc_direction" %in% colnames(result$junction_df))
  expect_true("y_start" %in% colnames(result$junction_df))
  expect_true("y_end" %in% colnames(result$junction_df))

  # Check alternating directions
  expect_equal(result$junction_df$arc_direction, c("up", "down", "up"))

  # Test up direction
  result_up <- process_sashimi_data(
    coverage_data, junction_data, "chr1:1000-2000",
    junction_direction = "up"
  )
  expect_true(all(result_up$junction_df$arc_direction == "up"))

  # Test down direction
  result_down <- process_sashimi_data(
    coverage_data, junction_data, "chr1:1000-2000",
    junction_direction = "down"
  )
  expect_true(all(result_down$junction_df$arc_direction == "down"))
  expect_true(all(result_down$junction_df$y_start == 0))
  expect_true(all(result_down$junction_df$y_end == 0))
})
