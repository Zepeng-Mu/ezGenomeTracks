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

test_that("ez_arc creates an interaction track", {
  # Load example data
  data(example_interactions)

  # Create an interaction track
  p <- ez_arc(data = example_interactions, region = "chr1:1000000-1100000")

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
