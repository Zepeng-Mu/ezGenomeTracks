# Tests for ez_wrappers.R functions

test_that("ez_signal creates a signal track", {
  # Load example data
  data(example_signal)
  
  # Create a signal track
  p <- ez_signal(data = example_signal)
  
  # Check that it returns a ggplot object
  expect_s3_class(p, "ggplot")
  
  # Check that it has at least one layer
  expect_true(length(p$layers) > 0)
})

test_that("ez_peak creates a peak track", {
  # Load example data
  data(example_peaks)
  
  # Create a peak track
  p <- ez_peak(data = example_peaks)
  
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
  p <- ez_arc(data = example_interactions)
  
  # Check that it returns a ggplot object
  expect_s3_class(p, "ggplot")
  
  # Check that it has at least one layer
  expect_true(length(p$layers) > 0)
})

test_that("ez_hic creates a Hi-C track", {
  # Load example data
  data(example_hic)
  
  # Create a Hi-C track
  p <- ez_hic(data = example_hic)
  
  # Check that it returns a ggplot object
  expect_s3_class(p, "ggplot")
  
  # Check that it has at least one layer
  expect_true(length(p$layers) > 0)
})