# Tests for theme.R functions

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

test_that("ez_peak_theme returns a ggplot2 theme object", {
  # Get the theme
  theme <- ez_peak_theme()
  
  # Check that it's a theme object
  expect_s3_class(theme, "theme")
  
  # Check some key theme elements
  expect_true("axis.text.y" %in% names(theme))
  expect_true("axis.title.y" %in% names(theme))
})