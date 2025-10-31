# Tests for theme functions

test_that("ez_feature_theme returns a ggplot2 theme object", {
  # Get the theme
  theme <- ez_feature_theme()

  # Check that it's a theme object
  expect_s3_class(theme, "theme")

  # Check some key theme elements
  expect_true("axis.text.y" %in% names(theme))
  expect_true("axis.title.y" %in% names(theme))
  expect_true("panel.background" %in% names(theme))
})

test_that("ez_theme returns a ggplot2 theme object", {
  # Get the theme
  theme <- ez_theme()

  # Check that it's a theme object
  expect_s3_class(theme, "theme")

  # Check some key theme elements
  expect_true("axis.line.x" %in% names(theme))
  expect_true("panel.grid" %in% names(theme))
  expect_true("panel.background" %in% names(theme))
})

test_that("ez_signal_theme returns a ggplot2 theme object", {
  # Get the theme
  theme <- ez_signal_theme()

  # Check that it's a theme object
  expect_s3_class(theme, "theme")

  # Check some key theme elements
  expect_true("axis.text.y" %in% names(theme))
  expect_true("axis.title.y" %in% names(theme))
})

test_that("ez_gene_theme returns a ggplot2 theme object", {
  # Get the theme
  theme <- ez_gene_theme()

  # Check that it's a theme object
  expect_s3_class(theme, "theme")

  # Check some key theme elements
  expect_true("axis.text.y" %in% names(theme))
  expect_true("axis.title.y" %in% names(theme))
})

test_that("ez_signal_theme handles y_axis_style parameter", {
  # Test different y_axis_style options
  theme_none <- ez_signal_theme(y_axis_style = "none")
  theme_simple <- ez_signal_theme(y_axis_style = "simple")
  theme_full <- ez_signal_theme(y_axis_style = "full")

  # All should be theme objects
  expect_s3_class(theme_none, "theme")
  expect_s3_class(theme_simple, "theme")
  expect_s3_class(theme_full, "theme")

  # They should be different (though this is hard to test directly)
  # We just verify they don't error and return theme objects
})

test_that("themes can be applied to plots", {
  # Create a simple test plot
  test_data <- data.frame(x = 1:10, y = 1:10)
  p <- ggplot2::ggplot(test_data, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point()

  # Test applying different themes
  p1 <- p + ez_theme()
  expect_s3_class(p1, "ggplot")

  p2 <- p + ez_signal_theme()
  expect_s3_class(p2, "ggplot")

  p3 <- p + ez_gene_theme()
  expect_s3_class(p3, "ggplot")

  p4 <- p + ez_feature_theme()
  expect_s3_class(p4, "ggplot")

  # All should still be valid ggplot objects
  expect_true(length(p1$layers) > 0)
  expect_true(length(p2$layers) > 0)
  expect_true(length(p3$layers) > 0)
  expect_true(length(p4$layers) > 0)
})
