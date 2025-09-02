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

test_that("ez_gene creates a gene track", {
  # Load example data
  data(example_genes)
  
  # Create a gene track
  p <- ez_gene(data = example_genes)
  
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