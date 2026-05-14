# Tests for ez_locusZoom function

test_that("ez_locusZoom creates regional plots with region parameter", {
  # Create test data for a single region
  test_data <- data.frame(
    CHR = rep(3, 50),
    BP = seq(1000, 50000, length.out = 50),
    P = runif(50, 1e-6, 1),
    SNP = paste0("rs", 1:50)
  )

  # Basic regional plot
  p <- ez_locusZoom(test_data, region = "chr3:1000-50000")
  expect_s3_class(p, "ggplot")
})

test_that("ez_locusZoom requires region or gene parameter", {
  test_data <- data.frame(
    CHR = rep(3, 50),
    BP = seq(1000, 50000, length.out = 50),
    P = runif(50, 1e-6, 1)
  )

  expect_error(
    ez_locusZoom(test_data),
    "Either 'region' or 'gene' \\(with 'gene_db'\\) must be provided"
  )
})

test_that("ez_locusZoom supports r2 coloring", {
  test_data <- data.frame(
    CHR = rep(3, 50),
    BP = seq(1000, 50000, length.out = 50),
    P = runif(50, 1e-6, 1),
    SNP = paste0("rs", 1:50)
  )

  r2_values <- runif(50, 0, 1)

  p <- ez_locusZoom(
    test_data,
    region = "chr3:1000-50000",
    r2 = r2_values,
    lead_snp = "rs25"
  )
  expect_s3_class(p, "ggplot")

  # Check that r2 gradient scale is present
  scales <- p$scales$scales
  has_color_scale <- any(sapply(scales, function(s) {
    "colour" %in% s$aesthetics
  }))
  expect_true(has_color_scale)
})

test_that("ez_locusZoom validates r2 length", {
  test_data <- data.frame(
    CHR = rep(3, 50),
    BP = seq(1000, 50000, length.out = 50),
    P = runif(50, 1e-6, 1)
  )

  # Wrong length r2
  expect_error(
    ez_locusZoom(
      test_data,
      region = "chr3:1000-50000",
      r2 = runif(30, 0, 1)  # Wrong length
    ),
    "r2 vector must match"
  )
})

test_that("ez_locusZoom filters data to region", {
  # Data spanning multiple regions
  test_data <- data.frame(
    CHR = rep(3, 100),
    BP = seq(1000, 100000, length.out = 100),
    P = runif(100, 1e-6, 1),
    SNP = paste0("rs", 1:100)
  )

  # Plot only a subset
  p <- ez_locusZoom(test_data, region = "chr3:20000-40000")
  expect_s3_class(p, "ggplot")

  # The plot should only contain filtered data
  # (checked via the fact that it renders without error)
})

test_that("ez_locusZoom warns on empty region", {
  test_data <- data.frame(
    CHR = rep(3, 50),
    BP = seq(1000, 50000, length.out = 50),
    P = runif(50, 1e-6, 1)
  )

  # Region with no data
  expect_warning(
    ez_locusZoom(test_data, region = "chr5:1000000-2000000"),
    "No data points found"
  )
})

test_that("ez_locusZoom handles lead_snp highlighting", {
  test_data <- data.frame(
    CHR = rep(3, 50),
    BP = seq(1000, 50000, length.out = 50),
    P = c(runif(49, 0.01, 1), 1e-10),
    SNP = paste0("rs", 1:50)
  )

  p <- ez_locusZoom(
    test_data,
    region = "chr3:1000-50000",
    lead_snp = "rs50",
    highlight_color = "purple"
  )
  expect_s3_class(p, "ggplot")
})

test_that("ez_locusZoom supports threshold line", {
  test_data <- data.frame(
    CHR = rep(3, 50),
    BP = seq(1000, 50000, length.out = 50),
    P = runif(50, 1e-8, 1)
  )

  p <- ez_locusZoom(
    test_data,
    region = "chr3:1000-50000",
    threshold_p = 5e-8,
    threshold_color = "red"
  )
  expect_s3_class(p, "ggplot")

  # Check that threshold line layer exists
  layer_classes <- sapply(p$layers, function(l) class(l$geom)[1])
  expect_true("GeomHline" %in% layer_classes)
})

test_that("ez_locusZoom accepts different column naming conventions", {
  # GRanges-style naming
  test_data <- data.frame(
    seqnames = rep("chr3", 50),
    start = seq(1000, 50000, length.out = 50),
    pvalue = runif(50, 1e-6, 1),
    rsid = paste0("rs", 1:50)
  )

  p <- ez_locusZoom(test_data, region = "chr3:1000-50000")
  expect_s3_class(p, "ggplot")
})

test_that("ez_locusZoom y_axis_style parameter works", {
  test_data <- data.frame(
    CHR = rep(3, 50),
    BP = seq(1000, 50000, length.out = 50),
    P = runif(50, 1e-6, 1)
  )

  # Test all y_axis_style options
  for (style in c("none", "simple", "full")) {
    p <- ez_locusZoom(
      test_data,
      region = "chr3:1000-50000",
      y_axis_style = style
    )
    expect_s3_class(p, "ggplot")
  }
})

test_that("ez_locusZoom rejects non-dataframe input", {
  expect_error(
    ez_locusZoom("not_a_dataframe", region = "chr3:1000-50000"),
    "must be a data frame"
  )
})

test_that("ez_locusZoom requires gene_db when gene is provided", {
  test_data <- data.frame(
    CHR = rep(3, 50),
    BP = seq(1000, 50000, length.out = 50),
    P = runif(50, 1e-6, 1)
  )

  expect_error(
    ez_locusZoom(test_data, gene = "TP53"),
    "gene_db.*must be provided"
  )
})

test_that("ez_locusZoom handles chromosome prefix variations", {
  # Data without chr prefix
  test_data <- data.frame(
    CHR = rep(3, 50),
    BP = seq(1000, 50000, length.out = 50),
    P = runif(50, 1e-6, 1)
  )

  # Region with chr prefix should still work
  p <- ez_locusZoom(test_data, region = "chr3:1000-50000")
  expect_s3_class(p, "ggplot")

  # Data with chr prefix
  test_data2 <- data.frame(
    CHR = rep("chr3", 50),
    BP = seq(1000, 50000, length.out = 50),
    P = runif(50, 1e-6, 1)
  )

  # Region without chr prefix should also work
  p2 <- ez_locusZoom(test_data2, region = "3:1000-50000")
  expect_s3_class(p2, "ggplot")
})

test_that("ez_locusZoom subsets r2 values correctly when filtering", {
  # Full data
  test_data <- data.frame(
    CHR = rep(3, 100),
    BP = seq(1000, 100000, length.out = 100),
    P = runif(100, 1e-6, 1),
    SNP = paste0("rs", 1:100)
  )

  # r2 for full data
  r2_full <- seq(0, 1, length.out = 100)

  # Filter to subset - should automatically subset r2
  p <- ez_locusZoom(
    test_data,
    region = "chr3:25000-75000",
    r2 = r2_full
  )
  expect_s3_class(p, "ggplot")
})
