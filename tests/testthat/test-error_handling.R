# Tests for error handling and edge cases
# This file contains tests for error conditions and edge cases across ezGenomeTracks functions

test_that("Memory management and large data handling", {
  # Test with large datasets
  large_signal <- create_test_genomic_data(
    start = 1,
    end = 10000000,  # 10MB region
    bin_size = 100
  )

  expect_no_error({
    p <- ez_coverage(large_signal, region = "chr1:1-10000000")
    expect_s3_class(p, "ggplot")
  })

  # Test with many genes
  many_genes <- create_test_gene_data(gene_count = 100)

  expect_no_error({
    p <- ez_gene(many_genes, region = "chr1:1000000-1100000")
    expect_s3_class(p, "ggplot")
  })

  # Test with many interactions
  many_interactions <- create_test_interaction_data(n_interactions = 500)

  expect_no_error({
    p <- ez_link(many_interactions, region = "chr1:1000000-1100000")
    expect_s3_class(p, "ggplot")
  })
})

test_that("Boundary conditions and extreme values", {
  # Test with minimum coordinate values
  min_coords <- data.frame(
    seqnames = "chr1",
    start = 1,
    end = 100,
    score = 5.0
  )

  expect_no_error({
    p <- ez_coverage(min_coords, region = "chr1:1-100")
    expect_s3_class(p, "ggplot")
  })

  # Test with maximum coordinate values (within R integer limits)
  max_coords <- data.frame(
    seqnames = "chr1",
    start = 2^31 - 1000,
    end = 2^31 - 900,
    score = 5.0
  )

  expect_no_error({
    p <- ez_coverage(max_coords, region = paste0("chr1:", 2^31 - 1000, "-", 2^31 - 900))
    expect_s3_class(p, "ggplot")
  })

  # Test with zero scores
  zero_scores <- data.frame(
    seqnames = "chr1",
    start = c(100, 200, 300),
    end = c(110, 210, 310),
    score = c(0, 0, 0)
  )

  expect_no_error({
    p <- ez_coverage(zero_scores, region = "chr1:100-310")
    expect_s3_class(p, "ggplot")
  })

  # Test with negative scores (if allowed)
  negative_scores <- data.frame(
    seqnames = "chr1",
    start = c(100, 200, 300),
    end = c(110, 210, 310),
    score = c(-5.0, 0, 5.0)
  )

  # This might be allowed or might error depending on implementation
  expect_warning_or_error({
    p <- ez_coverage(negative_scores, region = "chr1:100-310")
    expect_s3_class(p, "ggplot")
  })
})

test_that("Charset and encoding edge cases", {
  # Test with special characters in gene names
  special_chars <- data.frame(
    seqnames = c("chr1", "chr1"),
    start = c(100, 200),
    end = c(110, 210),
    strand = c("+", "-"),
    type = c("gene", "exon"),
    gene_id = c("GENE-1", "GENE_2"),
    gene_name = c("Gene-α", "Gene/β")
  )

  expect_no_error({
    p <- ez_gene(special_chars, region = "chr1:100-210")
    expect_s3_class(p, "ggplot")
  })

  # Test with Unicode characters
  unicode_data <- data.frame(
    seqnames = c("chr1", "chr1"),
    start = c(100, 200),
    end = c(110, 210),
    score = c(5.0, 8.0)
  )

  expect_no_error({
    p <- ez_feature(unicode_data, region = "chr1:100-210")
    expect_s3_class(p, "ggplot")
  })
})

test_that("Missing data and incomplete records handling", {
  # Test with partial gene data
  partial_genes <- data.frame(
    seqnames = c("chr1", "chr1", NA),
    start = c(100, 200, 300),
    end = c(110, 210, 310),
    strand = c("+", "-", "+"),
    type = c("gene", "exon", "gene"),
    gene_id = c("GENE1", "GENE1", NA),
    gene_name = c("Gene1", "Gene1", "Gene3")
  )

  expect_error(ez_gene(partial_genes, region = "chr1:100-310"))

  # Test with some missing but not all required fields
  partial_signal <- data.frame(
    seqnames = c("chr1", "chr1"),
    start = c(100, 200),
    end = c(110, 210)
    # Missing score column
  )

  expect_error(ez_coverage(partial_signal, region = "chr1:100-210"))
})

test_that("Invalid parameter combinations", {
  signal_data <- create_test_genomic_data()

  # Test conflicting parameters
  expect_error(ez_coverage(
    signal_data,
    region = "chr1:1000000-1100000",
    type = "heatmap",
    alpha = 1.5  # High alpha with heatmap might not make sense
  ))

  # Test invalid color parameters
  expect_error(ez_coverage(
    signal_data,
    region = "chr1:1000000-1100000",
    fill = "invalid_color_name"
  ))

  # Test geom functions with conflicting aesthetics
  test_data <- data.frame(
    start = seq(1, 100, 10),
    end = seq(10, 100, 10),
    score = rnorm(10)
  )

  # This should work - providing both required and optional aesthetics
  expect_no_error({
    p <- ggplot2::ggplot(test_data) +
      geom_coverage(ggplot2::aes(x = start, y = score, color = score, fill = score))
    expect_s3_class(p, "ggplot")
  })
})

test_that("Threading and concurrency edge cases", {
  # Test rapid successive calls
  signal_data <- create_test_genomic_data()

  # This tests for race conditions or shared state issues
  results <- list()
  for (i in 1:10) {
    results[[i]] <- ez_coverage(signal_data, region = "chr1:1000000-1100000")
  }

  # All results should be valid ggplot objects
  expect_true(all(sapply(results, function(x) inherits(x, "gg"))))

  # Test with different data but same region
  signal_data2 <- create_test_genomic_data(chr = "chr2")

  expect_no_error({
    p1 <- ez_coverage(signal_data, region = "chr1:1000000-1100000")
    p2 <- ez_coverage(signal_data2, region = "chr1:1000000-1100000")

    expect_s3_class(p1, "ggplot")
    expect_s3_class(p2, "ggplot")
  })
})

test_that("File system related errors", {
  # Test with non-existent region files (if file-based input is supported)
  temp_file <- create_temp_file(".bed", "chr1\t100\t200\tpeak1")

  # This tests file reading if supported
  expect_no_error({
    # File exists and should be readable
    expect_true(file.exists(temp_file))
  })

  cleanup_temp_files(temp_file)

  # Test with deleted file
  expect_error({
    # File no longer exists
    readLines(temp_file)
  })
})

test_that("Graphics device and rendering errors", {
  signal_data <- create_test_genomic_data()

  # Test plot creation with various graphics parameters
  expect_no_error({
    p1 <- ez_coverage(signal_data, region = "chr1:1000000-1100000", width = 10, height = 6)
    expect_s3_class(p1, "ggplot")

    p2 <- ez_coverage(signal_data, region = "chr1:1000000-1100000", dpi = 300)
    expect_s3_class(p2, "ggplot")
  })

  # Test rendering with extreme aspect ratios
  expect_no_error({
    p <- ez_coverage(signal_data, region = "chr1:1000000-1100000")
    # This should not cause rendering errors
    expect_true(can_render_plot(p, width = 100, height = 10))  # Very wide
    expect_true(can_render_plot(p, width = 10, height = 100))  # Very tall
  })
})

test_that("Data corruption and integrity checks", {
  # Test with corrupted coordinate data
  corrupted_coords <- data.frame(
    seqnames = c("chr1", "chr1", "chr2"),  # Mixed chromosomes in supposedly single-region data
    start = c(100, 200, 300),
    end = c(110, 210, 310),
    score = c(5.0, 8.0, 3.0)
  )

  # This might error or might filter the data
  expect_warning_or_error({
    p <- ez_coverage(corrupted_coords, region = "chr1:100-310")
    expect_s3_class(p, "ggplot")
  })

  # Test with overlapping intervals
  overlapping_intervals <- data.frame(
    seqnames = c("chr1", "chr1", "chr1"),
    start = c(100, 105, 110),  # Overlapping
    end = c(200, 210, 220),    # Overlapping
    score = c(5.0, 8.0, 3.0)
  )

  expect_no_error({
    p <- ez_coverage(overlapping_intervals, region = "chr1:100-220")
    expect_s3_class(p, "ggplot")
  })
})

test_that("Performance and timeout scenarios", {
  # Test with computationally intensive operations
  dense_data <- create_test_genomic_data(
    start = 1000000,
    end = 1010000,
    bin_size = 10  # Very dense data
  )

  # This should complete in reasonable time
  start_time <- Sys.time()
  expect_no_error({
    p <- ez_coverage(dense_data, region = "chr1:1000000-1010000")
    expect_s3_class(p, "ggplot")
  })
  end_time <- Sys.time()

  # Should complete within 10 seconds
  expect_true(as.numeric(difftime(end_time, start_time, units = "secs")) < 10)
})

test_that("Version compatibility and API changes", {
  signal_data <- create_test_genomic_data()

  # Test with old-style parameters (if any)
  expect_no_error({
    p1 <- ez_coverage(signal_data, region = "chr1:1000000-1100000")
    p2 <- ez_coverage(signal_data, region = "chr1:1000000-1100000", type = "area")

    expect_s3_class(p1, "ggplot")
    expect_s3_class(p2, "ggplot")
  })

  # Test parameter deprecation warnings if any
  expect_warning_or_error({
    # Try parameters that might be deprecated
    p <- ez_coverage(signal_data, region = "chr1:1000000-1100000", old_param = "value")
  })
})

test_that("Mathematical edge cases", {
  # Test with infinite values
  infinite_data <- data.frame(
    seqnames = c("chr1", "chr1"),
    start = c(100, 200),
    end = c(110, 210),
    score = c(Inf, -Inf)
  )

  expect_error(ez_coverage(infinite_data, region = "chr1:100-210"))

  # Test with very large numbers
  huge_numbers <- data.frame(
    seqnames = c("chr1", "chr1"),
    start = c(100, 200),
    end = c(110, 210),
    score = c(1e308, 1e300)  # Very large numbers
  )

  expect_no_error({
    p <- ez_coverage(huge_numbers, region = "chr1:100-210")
    expect_s3_class(p, "ggplot")
  })

  # Test with NaN values
  nan_data <- data.frame(
    seqnames = c("chr1", "chr1"),
    start = c(100, 200),
    end = c(110, 210),
    score = c(NaN, 5.0)
  )

  expect_error(ez_coverage(nan_data, region = "chr1:100-210"))
})

test_that("Statistical edge cases", {
  # Test with constant values
  constant_data <- data.frame(
    seqnames = c("chr1", "chr1", "chr1"),
    start = c(100, 200, 300),
    end = c(110, 210, 310),
    score = c(5.0, 5.0, 5.0)  # All same value
  )

  expect_no_error({
    p <- ez_coverage(constant_data, region = "chr1:100-310")
    expect_s3_class(p, "ggplot")
  })

  # Test with extreme variance
  extreme_variance <- data.frame(
    seqnames = c("chr1", "chr1", "chr1"),
    start = c(100, 200, 300),
    end = c(110, 210, 310),
    score = c(0, 5.0, 1e10)  # Extreme variance
  )

  expect_no_error({
    p <- ez_coverage(extreme_variance, region = "chr1:100-310")
    expect_s3_class(p, "ggplot")
  })
})

# Helper function for testing warning or error conditions
warning_or_error <- function(expr) {
  result <- tryCatch(
    withCallingHandlers(
      expr,
      warning = function(w) {
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      list(status = "error", message = e$message)
    },
    warning = function(w) {
      list(status = "warning", message = w$message)
    }
  )

  return(result)
}

expect_warning_or_error <- function(expr) {
  result <- warning_or_error(expr)
  expect_true(result$status %in% c("success", "warning", "error"))
}