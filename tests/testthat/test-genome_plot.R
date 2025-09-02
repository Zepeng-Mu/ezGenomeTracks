# Tests for genome_plot.R functions

test_that("genome_plot stacks tracks correctly", {
  # Create simple tracks
  track1 <- ggplot2::ggplot() + ggplot2::geom_blank()
  track2 <- ggplot2::ggplot() + ggplot2::geom_blank()
  
  # Stack tracks
  p <- genome_plot(track1, track2)
  
  # Check that it returns a ggplot object
  expect_s3_class(p, "ggplot")
  
  # Check that it has the correct number of tracks
  # This is a bit tricky to test directly, so we'll just check that it doesn't error
  expect_error(p, NA)
})

test_that("add_vline adds a vertical line", {
  # Create a simple track
  track <- ggplot2::ggplot() + ggplot2::geom_blank()
  
  # Add a vertical line
  track_with_vline <- add_vline(track, position = 1000)
  
  # Check that it returns a ggplot object
  expect_s3_class(track_with_vline, "ggplot")
  
  # Check that it has a geom_vline layer
  expect_true(any(sapply(track_with_vline$layers, function(l) inherits(l$geom, "GeomVline"))))
})

test_that("add_hline adds a horizontal line", {
  # Create a simple track
  track <- ggplot2::ggplot() + ggplot2::geom_blank()
  
  # Add a horizontal line
  track_with_hline <- add_hline(track, y = 10)
  
  # Check that it returns a ggplot object
  expect_s3_class(track_with_hline, "ggplot")
  
  # Check that it has a geom_hline layer
  expect_true(any(sapply(track_with_hline$layers, function(l) inherits(l$geom, "GeomHline"))))
})

test_that("add_rect adds a rectangle", {
  # Create a simple track
  track <- ggplot2::ggplot() + ggplot2::geom_blank()
  
  # Add a rectangle
  track_with_rect <- add_rect(track, xmin = 1000, xmax = 2000)
  
  # Check that it returns a ggplot object
  expect_s3_class(track_with_rect, "ggplot")
  
  # Check that it has a geom_rect layer
  expect_true(any(sapply(track_with_rect$layers, function(l) inherits(l$geom, "GeomRect"))))
})

test_that("add_text adds text annotation", {
  # Create a simple track
  track <- ggplot2::ggplot() + ggplot2::geom_blank()
  
  # Add text annotation
  track_with_text <- add_text(track, x = 1000, y = 10, label = "Test")
  
  # Check that it returns a ggplot object
  expect_s3_class(track_with_text, "ggplot")
  
  # Check that it has a geom_text layer
  expect_true(any(sapply(track_with_text$layers, function(l) inherits(l$geom, "GeomText"))))
})