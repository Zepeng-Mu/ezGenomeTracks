test_that("ez_default_palette returns correct number of colors", {
  # Full palette
  full <- ez_default_palette()
  expect_true(length(full) >= 9)
  expect_true(all(grepl("^#", full)))

  # Subset
  p5 <- ez_default_palette(5)
  expect_length(p5, 5)

  # Interpolated (more than palette size)
  p20 <- ez_default_palette(20)
  expect_length(p20, 20)

  # Single
  p1 <- ez_default_palette(1)
  expect_length(p1, 1)
})

test_that("ez_default_single_color returns valid colors for all types", {
  types <- c("coverage", "feature", "gene", "link", "sashimi", "manhattan")
  for (type in types) {
    color <- ez_default_single_color(type)
    expect_true(is.character(color), info = paste("Type:", type))
    expect_length(color, 1)
  }
})

test_that("ez_hic_palette returns colors or NULL", {
  result <- ez_hic_palette("cooler", 10)
  if (!is.null(result)) {
    expect_length(result, 10)
    expect_true(all(grepl("^#", result)))
  }
})

test_that("ez_hic_diverging_palette returns colors or NULL", {
  result <- ez_hic_diverging_palette(10)
  if (!is.null(result)) {
    expect_length(result, 10)
    expect_true(all(grepl("^#", result)))
  }
})

test_that("fallback colors work when vibeColors is not available", {
  # Test that the fallback palette is accessible
  expect_true(length(.ez_fallback_discrete) >= 9)
  expect_true(all(grepl("^#", .ez_fallback_discrete)))

  # Test that all single-color types have fallbacks
  types <- c("coverage", "feature", "gene", "link", "sashimi", "manhattan")
  for (type in types) {
    expect_true(!is.null(.ez_fallback_single[[type]]), info = paste("Type:", type))
  }
})
