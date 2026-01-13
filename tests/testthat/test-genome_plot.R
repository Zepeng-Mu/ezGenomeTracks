# Tests for vstack_plot.R functions

test_that("vstack_plot stacks tracks correctly", {
  # Create simple tracks
  track1 <- ggplot2::ggplot() + ggplot2::geom_blank()
  track2 <- ggplot2::ggplot() + ggplot2::geom_blank()

  # Stack tracks
  p <- vstack_plot(track1, track2)

  # Check that it returns a ggplot object
  expect_s3_class(p, "ggplot")
  expect_true(is_valid_ggplot(p))
  expect_true(can_render_plot(p))

  # Check that it has the correct number of tracks
  expect_true(length(p$layers) >= length(track1$layers) + length(track2$layers))
})

test_that("vstack_plot handles multiple tracks with different data", {
  # Create tracks with different types of genomic data
  signal_data <- create_test_genomic_data()
  gene_data <- create_test_gene_data()

  track1 <- ez_coverage(signal_data, region = "chr1:1000000-1100000")
  track2 <- ez_gene(gene_data, region = "chr1:1000000-1100000")

  # Stack tracks
  stacked_plot <- vstack_plot(track1, track2)

  expect_s3_class(stacked_plot, "ggplot")
  expect_true(is_valid_ggplot(stacked_plot))
  expect_true(can_render_plot(stacked_plot))

  # Should preserve data from both tracks
  expect_true(nrow(stacked_plot$data) >= nrow(track1$data))
})

test_that("add_vline adds vertical lines with various parameters", {
  # Create a track with real data
  signal_data <- create_test_genomic_data()
  track <- ez_coverage(signal_data, region = "chr1:1000000-1100000")

  # Test basic vertical line
  track_with_vline <- add_vline(track, position = 1050000)
  expect_s3_class(track_with_vline, "ggplot")
  expect_true(has_expected_layers(track_with_vline, c("GeomVline")))
  expect_true(can_render_plot(track_with_vline))

  # Test with custom styling
  track_with_styled_vline <- add_vline(
    track,
    position = 1050000,
    color = "red",
    linetype = "dashed",
    size = 2
  )
  expect_s3_class(track_with_styled_vline, "ggplot")
  expect_true(can_render_plot(track_with_styled_vline))

  # Test with multiple lines
  track_with_multiple_vlines <- track %>%
    add_vline(position = 1040000) %>%
    add_vline(position = 1060000) %>%
    add_vline(position = 1080000)

  expect_s3_class(track_with_multiple_vlines, "ggplot")
  expect_true(can_render_plot(track_with_multiple_vlines))

  # Count vline layers
  vline_count <- sum(sapply(track_with_multiple_vlines$layers, function(l) inherits(l$geom, "GeomVline")))
  expect_true(vline_count >= 3)
})

test_that("add_hline adds horizontal lines with various parameters", {
  # Create a track with real data
  signal_data <- create_test_genomic_data()
  track <- ez_coverage(signal_data, region = "chr1:1000000-1100000")

  # Test basic horizontal line
  track_with_hline <- add_hline(track, y = 10)
  expect_s3_class(track_with_hline, "ggplot")
  expect_true(has_expected_layers(track_with_hline, c("GeomHline")))
  expect_true(can_render_plot(track_with_hline))

  # Test with custom styling
  track_with_styled_hline <- add_hline(
    track,
    y = 10,
    color = "blue",
    linetype = "dotted",
    size = 1.5
  )
  expect_s3_class(track_with_styled_hline, "ggplot")
  expect_true(can_render_plot(track_with_styled_hline))

  # Test with multiple lines
  track_with_multiple_hlines <- track %>%
    add_hline(y = 5) %>%
    add_hline(y = 10) %>%
    add_hline(y = 15)

  expect_s3_class(track_with_multiple_hlines, "ggplot")
  expect_true(can_render_plot(track_with_multiple_hlines))

  # Count hline layers
  hline_count <- sum(sapply(track_with_multiple_hlines$layers, function(l) inherits(l$geom, "GeomHline")))
  expect_true(hline_count >= 3)
})

test_that("add_rect adds rectangles with various parameters", {
  # Create a track with real data
  signal_data <- create_test_genomic_data()
  track <- ez_coverage(signal_data, region = "chr1:1000000-1100000")

  # Test basic rectangle
  track_with_rect <- add_rect(track, xmin = 1050000, xmax = 1060000)
  expect_s3_class(track_with_rect, "ggplot")
  expect_true(has_expected_layers(track_with_rect, c("GeomRect")))
  expect_true(can_render_plot(track_with_rect))

  # Test with custom styling
  track_with_styled_rect <- add_rect(
    track,
    xmin = 1050000,
    xmax = 1060000,
    fill = "lightblue",
    color = "darkblue",
    alpha = 0.5
  )
  expect_s3_class(track_with_styled_rect, "ggplot")
  expect_true(can_render_plot(track_with_styled_rect))

  # Test with y-range
  track_with_y_rect <- add_rect(
    track,
    xmin = 1050000,
    xmax = 1060000,
    ymin = 5,
    ymax = 15
  )
  expect_s3_class(track_with_y_rect, "ggplot")
  expect_true(can_render_plot(track_with_y_rect))

  # Test with multiple rectangles
  track_with_multiple_rects <- track %>%
    add_rect(xmin = 1040000, xmax = 1050000, alpha = 0.3) %>%
    add_rect(xmin = 1060000, xmax = 1070000, alpha = 0.3) %>%
    add_rect(xmin = 1080000, xmax = 1090000, alpha = 0.3)

  expect_s3_class(track_with_multiple_rects, "ggplot")
  expect_true(can_render_plot(track_with_multiple_rects))

  # Count rect layers
  rect_count <- sum(sapply(track_with_multiple_rects$layers, function(l) inherits(l$geom, "GeomRect")))
  expect_true(rect_count >= 3)
})

test_that("add_text adds text annotations with various parameters", {
  # Create a track with real data
  signal_data <- create_test_genomic_data()
  track <- ez_coverage(signal_data, region = "chr1:1000000-1100000")

  # Test basic text annotation
  track_with_text <- add_text(track, x = 1050000, y = 10, label = "Peak Region")
  expect_s3_class(track_with_text, "ggplot")
  expect_true(has_expected_layers(track_with_text, c("GeomText")))
  expect_true(can_render_plot(track_with_text))

  # Test with custom styling
  track_with_styled_text <- add_text(
    track,
    x = 1050000,
    y = 10,
    label = "Important",
    color = "red",
    size = 12,
    fontface = "bold"
  )
  expect_s3_class(track_with_styled_text, "ggplot")
  expect_true(can_render_plot(track_with_styled_text))

  # Test with special characters
  track_with_special_chars <- add_text(
    track,
    x = 1050000,
    y = 10,
    label = "αβγ δε",
    color = "darkgreen"
  )
  expect_s3_class(track_with_special_chars, "ggplot")
  expect_true(can_render_plot(track_with_special_chars))

  # Test with multiple text annotations
  track_with_multiple_texts <- track %>%
    add_text(x = 1040000, y = 5, label = "Start") %>%
    add_text(x = 1060000, y = 10, label = "Middle") %>%
    add_text(x = 1080000, y = 15, label = "End")

  expect_s3_class(track_with_multiple_texts, "ggplot")
  expect_true(can_render_plot(track_with_multiple_texts))

  # Count text layers
  text_count <- sum(sapply(track_with_multiple_texts$layers, function(l) inherits(l$geom, "GeomText")))
  expect_true(text_count >= 3)
})

test_that("add_highlight adds highlight regions with various parameters", {
  # Create a track with real data
  signal_data <- create_test_genomic_data()
  track <- ez_coverage(signal_data, region = "chr1:1000000-1100000")

  # Test basic highlight with data frame
  highlight_data <- data.frame(
    start = 1050000,
    end = 1060000,
    score = 0.8
  )

  track_with_highlight <- add_highlight(track, highlight_data)
  expect_s3_class(track_with_highlight, "ggplot")
  expect_true(can_render_plot(track_with_highlight))

  # Test with custom styling
  track_with_styled_highlight <- add_highlight(
    track,
    highlight_data,
    fill = "yellow",
    alpha = 0.4,
    color = "orange"
  )
  expect_s3_class(track_with_styled_highlight, "ggplot")
  expect_true(can_render_plot(track_with_styled_highlight))

  # Test with multiple highlight regions
  multiple_highlights <- data.frame(
    start = c(1040000, 1060000, 1080000),
    end = c(1045000, 1065000, 1085000),
    score = c(0.6, 0.8, 0.7)
  )

  track_with_multiple_highlights <- add_highlight(track, multiple_highlights)
  expect_s3_class(track_with_multiple_highlights, "ggplot")
  expect_true(can_render_plot(track_with_multiple_highlights))

  # Verify data preservation
  expect_true(nrow(track_with_highlight$data) >= nrow(highlight_data))
})

test_that("annotation functions handle edge cases correctly", {
  # Create a track with real data
  signal_data <- create_test_genomic_data()
  track <- ez_coverage(signal_data, region = "chr1:1000000-1100000")

  # Test with coordinates outside the region
  expect_no_error({
    track_with_outside_vline <- add_vline(track, position = 500000)  # Outside region
    expect_s3_class(track_with_outside_vline, "ggplot")
  })

  # Test with extreme coordinates
  expect_no_error({
    track_with_extreme_coords <- add_rect(
      track,
      xmin = 1,
      xmax = 2^31 - 1
    )
    expect_s3_class(track_with_extreme_coords, "ggplot")
  })

  # Test with empty labels
  expect_no_error({
    track_with_empty_text <- add_text(track, x = 1050000, y = 10, label = "")
    expect_s3_class(track_with_empty_text, "ggplot")
  })

  # Test with NA values
  expect_error(add_vline(track, position = NA))
  expect_error(add_text(track, x = 1050000, y = 10, label = NA))
})

test_that("combined annotation workflow", {
  # Create a track with real data
  signal_data <- create_test_genomic_data()
  track <- ez_coverage(signal_data, region = "chr1:1000000-1100000")

  # Create a fully annotated track
  annotated_track <- track %>%
    # Add highlight region
    add_highlight(
      data.frame(
        start = 1050000,
        end = 1060000,
        score = 0.8
      ),
      fill = "lightblue",
      alpha = 0.3
    ) %>%
    # Add reference lines
    add_vline(position = 1055000, color = "red", linetype = "dashed") %>%
    add_hline(y = 10, color = "green", linetype = "dotted") %>%
    # Add annotations
    add_text(
      x = 1055000,
      y = 15,
      label = "Highlighted Region",
      color = "darkblue",
      fontface = "bold",
      size = 14
    ) %>%
    # Add boundary rectangles
    add_rect(
      xmin = 1000000,
      xmax = 1010000,
      ymin = 0,
      ymax = 20,
      alpha = 0.1,
      fill = "gray"
    )

  # Verify the annotated track
  expect_s3_class(annotated_track, "ggplot")
  expect_true(is_valid_ggplot(annotated_track))
  expect_true(can_render_plot(annotated_track))

  # Verify all annotation types are present
  layer_classes <- sapply(annotated_track$layers, function(l) class(l$geom)[1])
  expect_true("GeomRect" %in% layer_classes)
  expect_true("GeomVline" %in% layer_classes)
  expect_true("GeomHline" %in% layer_classes)
  expect_true("GeomText" %in% layer_classes)

  # Should have more layers than the original track
  expect_true(length(annotated_track$layers) > length(track$layers))
})

test_that("annotation functions preserve plot functionality", {
  # Create real genomic tracks
  signal_data <- create_test_genomic_data()
  gene_data <- create_test_gene_data()

  base_track1 <- ez_coverage(signal_data, region = "chr1:1000000-1100000")
  base_track2 <- ez_gene(gene_data, region = "chr1:1000000-1100000")

  # Annotate both tracks
  annotated_track1 <- base_track1 %>%
    add_highlight(
      data.frame(start = 1050000, end = 1060000, score = 0.8)
    ) %>%
    add_vline(position = 1055000)

  annotated_track2 <- base_track2 %>%
    add_text(x = 1050000, y = 0.5, label = "Gene Region")

  # Stack annotated tracks
  final_plot <- vstack_plot(annotated_track1, annotated_track2)

  # Verify final plot functionality
  expect_s3_class(final_plot, "ggplot")
  expect_true(is_valid_ggplot(final_plot))
  expect_true(can_render_plot(final_plot))

  # Verify data is preserved
  expect_true(nrow(final_plot$data) > 0)
})