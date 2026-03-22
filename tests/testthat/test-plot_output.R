# Tests for plot output validation
# This file contains tests for validating the structural integrity of plot outputs

test_that("ez_coverage plots have correct structure", {
  signal_data <- create_test_genomic_data()

  # Test different coverage plot types
  coverage_types <- c("area", "line", "heatmap")

  for (plot_type in coverage_types) {
    p <- ez_coverage(signal_data, region = "chr1:1000000-1100000", type = plot_type)

    # Basic structural validation
    expect_true(is_valid_ggplot(p))
    expect_true(has_expected_layers(p, c("Geom")))
    expect_true(can_render_plot(p))

    # Check for genomic scales
    expect_true(has_expected_scales(p, c("ScaleXContinuous")))
  }

  # Test with different parameters
  p_custom <- ez_coverage(
    signal_data,
    region = "chr1:1000000-1100000",
    fill = "blue",
    alpha = 0.7,
    color = "darkblue"
  )

  expect_true(is_valid_ggplot(p_custom))
  expect_true(can_render_plot(p_custom))
})

test_that("ez_feature plots have correct structure", {
  feature_data <- data.frame(
    seqnames = c("chr1", "chr1", "chr1"),
    start = c(1000, 2000, 3000),
    end = c(1100, 2100, 3100),
    score = c(0.5, 0.8, 0.3),
    name = c("Peak1", "Peak2", "Peak3")
  )

  p <- ez_feature(feature_data, region = "chr1:1000-3200")

  # Basic structural validation
  expect_true(is_valid_ggplot(p))
  expect_true(has_expected_layers(p, c("GeomRect", "GeomText")))
  expect_true(can_render_plot(p))

  # Check that data is preserved
  expect_equal(nrow(p$data), nrow(feature_data))
})

test_that("ez_gene plots have correct structure", {
  gene_data <- create_test_gene_data()

  p <- ez_gene(gene_data, region = "chr1:1000000-1100000")

  # Basic structural validation
  expect_true(is_valid_ggplot(p))
  expect_true(has_expected_layers(p, c("GeomSegment", "GeomRect")))
  expect_true(can_render_plot(p))

  # Check for strand-based y positioning if applicable
  expect_true(has_expected_scales(p, c("ScaleXContinuous")))
})

test_that("ez_link plots have correct structure", {
  interaction_data <- create_test_interaction_data()

  p <- ez_link(interaction_data, region = "chr1:1000000-1100000")

  # Basic structural validation
  expect_true(is_valid_ggplot(p))
  expect_true(has_expected_layers(p, c("GeomSegment", "GeomCurve")))
  expect_true(can_render_plot(p))

  # Check that interaction data is preserved
  expect_true(nrow(p$data) > 0)
})

test_that("ez_hic plots have correct structure", {
  hic_data <- create_test_hic_data()

  p <- ez_hic(hic_data, region = "chr1:1000000-1100000")

  # Basic structural validation
  expect_true(is_valid_ggplot(p))
  expect_true(has_expected_layers(p, c("GeomTile", "GeomRect")))
  expect_true(can_render_plot(p))

  # Check for Hi-C specific scales
  expect_true(has_expected_scales(p, c("ScaleXContinuous", "ScaleYContinuous")))
})

test_that("ez_sashimi plots have correct structure", {
  coverage_data <- create_test_genomic_data()
  junction_data <- data.frame(
    seqnames = "chr1",
    start = c(1050000, 1070000),
    end = c(1060000, 1080000),
    score = c(25, 50)
  )

  p <- ez_sashimi(coverage_data, junction_data, region = "chr1:1000000-1100000")

  # Basic structural validation
  expect_true(is_valid_ggplot(p))
  expect_true(has_expected_layers(p, c("GeomArea", "GeomCurve", "GeomText")))
  expect_true(can_render_plot(p))

  # Check that both coverage and junction data are present
  expect_true(nrow(p$data) >= nrow(coverage_data))
})

test_that("vstack_plot stacking works correctly", {
  # Create multiple tracks
  signal_data <- create_test_genomic_data()
  feature_data <- data.frame(
    seqnames = "chr1",
    start = 1050000,
    end = 1060000,
    score = 0.8
  )
  gene_data <- create_test_gene_data()

  track1 <- ez_coverage(signal_data, region = "chr1:1000000-1100000")
  track2 <- ez_feature(feature_data, region = "chr1:1000000-1100000")
  track3 <- ez_gene(gene_data, region = "chr1:1000000-1100000")

  # Stack tracks
  stacked_plot <- vstack_plot(track1, track2, track3)

  # Basic structural validation
  expect_true(is_valid_ggplot(stacked_plot))
  expect_true(can_render_plot(stacked_plot))

  # Should have layers from all tracks
  expect_true(length(stacked_plot$layers) >= length(track1$layers) +
              length(track2$layers) + length(track3$layers))
})

test_that("add_* annotation functions work correctly", {
  # Create base track
  signal_data <- create_test_genomic_data()
  base_track <- ez_coverage(signal_data, region = "chr1:1000000-1100000")

  # Test add_vline
  track_with_vline <- add_vline(base_track, position = 1050000)
  expect_true(is_valid_ggplot(track_with_vline))
  expect_true(has_expected_layers(track_with_vline, c("GeomVline")))
  expect_true(can_render_plot(track_with_vline))

  # Test add_hline
  track_with_hline <- add_hline(base_track, y = 10)
  expect_true(is_valid_ggplot(track_with_hline))
  expect_true(has_expected_layers(track_with_hline, c("GeomHline")))
  expect_true(can_render_plot(track_with_hline))

  # Test add_rect
  track_with_rect <- add_rect(base_track, xmin = 1050000, xmax = 1060000)
  expect_true(is_valid_ggplot(track_with_rect))
  expect_true(has_expected_layers(track_with_rect, c("GeomRect")))
  expect_true(can_render_plot(track_with_rect))

  # Test add_text
  track_with_text <- add_text(base_track, x = 1050000, y = 10, label = "Test Region")
  expect_true(is_valid_ggplot(track_with_text))
  expect_true(has_expected_layers(track_with_text, c("GeomText")))
  expect_true(can_render_plot(track_with_text))

  # Test add_highlight
  highlight_data <- data.frame(
    start = 1050000,
    end = 1060000,
    score = 0.8
  )
  track_with_highlight <- add_highlight(base_track, highlight_data)
  expect_true(is_valid_ggplot(track_with_highlight))
  expect_true(can_render_plot(track_with_highlight))
})

test_that("geom_* functions produce valid plots", {
  # Test geom_coverage
  coverage_data <- data.frame(
    start = seq(100, 1000, 50),
    end = seq(150, 1050, 50),
    score = runif(18, 0, 10)
  )

  p1 <- ggplot2::ggplot(coverage_data) + geom_coverage()
  expect_true(is_valid_ggplot(p1))
  expect_true(can_render_plot(p1))

  # Test geom_feature
  feature_data <- data.frame(
    start = c(200, 400, 600),
    end = c(250, 450, 650),
    score = c(0.5, 0.8, 0.3)
  )

  p2 <- ggplot2::ggplot(feature_data) + geom_feature()
  expect_true(is_valid_ggplot(p2))
  expect_true(can_render_plot(p2))

  # Test geom_gene
  gene_data <- data.frame(
    start = c(100, 200),
    end = c(300, 400),
    type = c("gene", "exon"),
    strand = c("+", "-"),
    gene_name = c("Gene1", "Gene1")
  )

  p3 <- ggplot2::ggplot(gene_data) + geom_gene()
  expect_true(is_valid_ggplot(p3))
  expect_true(can_render_plot(p3))

  # Test geom_manhattan
  manhattan_data <- data.frame(
    CHR = rep(1:3, each = 10),
    BP = rep(seq(1000, 10000, 1000), 3),
    P = runif(30, 1e-8, 1),
    SNP = paste0("rs", 1:30)
  )

  p4 <- ggplot2::ggplot(manhattan_data) +
    geom_manhattan(ggplot2::aes(CHR = CHR, BP = BP, P = P))
  expect_true(is_valid_ggplot(p4))
  expect_true(can_render_plot(p4))

  # Test geom_arc
  arc_data <- data.frame(
    start1 = c(100, 200),
    start2 = c(300, 400),
    score = c(0.5, 0.8)
  )

  p5 <- ggplot2::ggplot(arc_data) +
    geom_arc(ggplot2::aes(x = start1, y = 0, xend = start2, yend = 0))
  expect_true(is_valid_ggplot(p5))
  expect_true(can_render_plot(p5))

  # Test geom_hic
  hic_data <- data.frame(
    pos1 = rep(seq(1000, 5000, 1000), 3),
    pos2 = rep(seq(1000, 5000, 1000), each = 3),
    count = runif(9, 0, 100)
  )

  p6 <- ggplot2::ggplot(hic_data) +
    geom_hic(ggplot2::aes(x = pos1, y = pos2, fill = count))
  expect_true(is_valid_ggplot(p6))
  expect_true(can_render_plot(p6))
})

test_that("Plot themes and styling are applied correctly", {
  signal_data <- create_test_genomic_data()

  # Test custom themes
  p1 <- ez_coverage(signal_data, region = "chr1:1000000-1100000")
  expect_true(!is.null(p1$theme))

  # Test theme functions
  p2 <- ez_coverage_theme() + ez_coverage(signal_data, region = "chr1:1000000-1100000")
  expect_true(is_valid_ggplot(p2))
  expect_true(can_render_plot(p2))

  # Test with different palettes
  if (requireNamespace("RColorBrewer", quietly = TRUE)) {
    p3 <- ez_coverage(
      signal_data,
      region = "chr1:1000000-1100000",
      palette = "Blues"
    )
    expect_true(is_valid_ggplot(p3))
    expect_true(can_render_plot(p3))
  }
})

test_that("Multi-sample and multi-condition plots work", {
  # Create multi-sample signal data
  multi_sample_data <- create_test_signal_data(
    n_samples = 3,
    n_conditions = 2
  )

  # This should work with faceting or grouping
  expect_no_error({
    p <- ez_coverage(multi_sample_data, region = "chr1:1000000-1100000")
    expect_true(is_valid_ggplot(p))
    expect_true(can_render_plot(p))
  })

  # Test faceting if supported
  expect_no_error({
    p_faceted <- ez_coverage(
      multi_sample_data,
      region = "chr1:1000000-1100000",
      facet_by = "sample"
    )
    expect_true(is_valid_ggplot(p_faceted))
    expect_true(can_render_plot(p_faceted))
  })
})

test_that("Coordinate transformation and scaling works", {
  # Test with different coordinate ranges
  large_region_data <- data.frame(
    seqnames = "chr1",
    start = c(1e8, 1e8 + 1000, 1e8 + 2000),
    end = c(1e8 + 100, 1e8 + 1100, 1e8 + 2100),
    score = c(5.0, 8.0, 3.0)
  )

  p_large <- ez_coverage(
    large_region_data,
    region = paste0("chr1:", 1e8, "-", 1e8 + 2100)
  )

  expect_true(is_valid_ggplot(p_large))
  expect_true(can_render_plot(p_large))

  # Test with small regions
  small_region_data <- data.frame(
    seqnames = "chr1",
    start = c(1000, 1005, 1010),
    end = c(1002, 1007, 1012),
    score = c(5.0, 8.0, 3.0)
  )

  p_small <- ez_coverage(small_region_data, region = "chr1:1000-1015")
  expect_true(is_valid_ggplot(p_small))
  expect_true(can_render_plot(p_small))
})

test_that("Plot output formats and dimensions", {
  signal_data <- create_test_genomic_data()
  p <- ez_coverage(signal_data, region = "chr1:1000000-1100000")

  # Test different output formats
  expect_true(can_render_plot(p, width = 800, height = 600))   # Standard
  expect_true(can_render_plot(p, width = 1200, height = 800))  # Large
  expect_true(can_render_plot(p, width = 400, height = 300))   # Small
  expect_true(can_render_plot(p, width = 1600, height = 900))  # HD

  # Test extreme aspect ratios
  expect_true(can_render_plot(p, width = 1000, height = 100))   # Very wide
  expect_true(can_render_plot(p, width = 100, height = 1000))   # Very tall
})

test_that("Data preservation in plots", {
  # Create test data with known properties
  test_signal <- data.frame(
    seqnames = c("chr1", "chr1", "chr1"),
    start = c(1000, 2000, 3000),
    end = c(1100, 2100, 3100),
    score = c(5.5, 8.2, 3.1)
  )

  p <- ez_coverage(test_signal, region = "chr1:1000-3200")

  # Check that original data is preserved or transformed correctly
  expect_true(nrow(p$data) >= nrow(test_signal))

  # Check that coordinate columns are present
  expect_true(all(c("start", "end") %in% colnames(p$data)))

  # Check that score information is preserved
  if ("score" %in% colnames(p$data)) {
    expect_true(length(unique(p$data$score)) >= 1)
  }
})

test_that("Edge case plot structures", {
  # Plot with single data point
  single_point <- data.frame(
    seqnames = "chr1",
    start = 1000,
    end = 1100,
    score = 5.0
  )

  p_single <- ez_coverage(single_point, region = "chr1:1000-1100")
  expect_true(is_valid_ggplot(p_single))
  expect_true(can_render_plot(p_single))

  # Plot with very sparse data
  sparse_data <- data.frame(
    seqnames = c("chr1", "chr1", "chr1"),
    start = c(1000, 50000, 100000),
    end = c(1100, 50100, 100100),
    score = c(5.0, 8.0, 3.0)
  )

  p_sparse <- ez_coverage(sparse_data, region = "chr1:1000-100100")
  expect_true(is_valid_ggplot(p_sparse))
  expect_true(can_render_plot(p_sparse))

  # Plot with very dense data
  dense_data <- data.frame(
    seqnames = rep("chr1", 1000),
    start = seq(1000, 2000, length.out = 1000),
    end = seq(1001, 2001, length.out = 1000),
    score = runif(1000, 0, 10)
  )

  p_dense <- ez_coverage(dense_data, region = "chr1:1000-2001")
  expect_true(is_valid_ggplot(p_dense))
  expect_true(can_render_plot(p_dense))
})

test_that("Plot layering and composition", {
  signal_data <- create_test_genomic_data()
  base_track <- ez_coverage(signal_data, region = "chr1:1000000-1100000")

  # Add multiple annotations
  track_annotated <- base_track %>%
    add_vline(position = 1050000) %>%
    add_hline(y = 10) %>%
    add_rect(xmin = 1055000, xmax = 1065000, alpha = 0.3) %>%
    add_text(x = 1060000, y = 15, label = "Important Region")

  expect_true(is_valid_ggplot(track_annotated))
  expect_true(can_render_plot(track_annotated))

  # Should have layers from base + annotations
  expect_true(length(track_annotated$layers) > length(base_track$layers))
})

test_that("Performance validation for plot generation", {
  # Test that plots are generated within reasonable time
  large_data <- create_test_genomic_data(
    start = 1000000,
    end = 2000000,
    bin_size = 100
  )

  start_time <- Sys.time()
  p <- ez_coverage(large_data, region = "chr1:1000000-2000000")
  end_time <- Sys.time()

  expect_true(is_valid_ggplot(p))
  expect_true(can_render_plot(p))

  # Should complete within 5 seconds
  expect_true(as.numeric(difftime(end_time, start_time, units = "secs")) < 5)
})