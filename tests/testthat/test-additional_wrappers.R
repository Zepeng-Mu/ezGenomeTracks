# Tests for additional ez wrapper functions

test_that("ez_manhattan creates Manhattan plots", {
  # Create test GWAS data
  test_data <- data.frame(
    CHR = rep(1:3, each = 10),
    BP = rep(seq(1000, 10000, 1000), 3),
    P = runif(30, 1e-8, 1),
    SNP = paste0("rs", 1:30)
  )

  # Test basic Manhattan plot
  p1 <- ez_manhattan(test_data)
  expect_s3_class(p1, "ggplot")
  expect_true(length(p1$layers) > 0)

  # Test with custom parameters
  p2 <- ez_manhattan(test_data, chr = "CHR", bp = "BP", p = "P",
                     colors = c("red", "blue"), size = 1)
  expect_s3_class(p2, "ggplot")
  expect_true(length(p2$layers) > 0)

  # Test with highlighting
  highlight_data <- test_data[test_data$P < 1e-5, ]
  p3 <- ez_manhattan(test_data, highlight_snps = highlight_data)
  expect_s3_class(p3, "ggplot")
  expect_true(length(p3$layers) > 0)

  # Test with threshold line
  p4 <- ez_manhattan(test_data, threshold_p = 5e-8)
  expect_s3_class(p4, "ggplot")
  expect_true(length(p4$layers) > 0)

  # Test with R-squared coloring
  test_data$R2 <- runif(30, 0, 1)
  p5 <- ez_manhattan(test_data, r2 = test_data$R2, colorBy = "r2")
  expect_s3_class(p5, "ggplot")
  expect_true(length(p5$layers) > 0)

  # Test error handling
  expect_error(ez_manhattan("not_a_dataframe"))
})

test_that("ez_feature creates feature tracks from data frames", {
  # Create test peak data
  test_data <- data.frame(
    seqnames = "chr1",
    start = c(1000, 2000, 3000),
    end = c(1500, 2500, 3500),
    score = c(0.5, 0.8, 0.3)
  )

  # Test basic feature plot
  p1 <- ez_feature(test_data, "chr1:1000-4000")
  expect_s3_class(p1, "ggplot")
  expect_true(length(p1$layers) > 0)

  # Test with score coloring
  p2 <- ez_feature(test_data, "chr1:1000-4000", use_score = TRUE)
  expect_s3_class(p2, "ggplot")
  expect_true(length(p2$layers) > 0)

  # Test with custom parameters
  p3 <- ez_feature(test_data, "chr1:1000-4000",
                   color = "red", fill = "blue", alpha = 0.8, height = 0.5)
  expect_s3_class(p3, "ggplot")
  expect_true(length(p3$layers) > 0)

  # Test error handling
  expect_error(ez_feature("not_a_dataframe", "chr1:1000-4000"))
})

test_that("ez_arc creates interaction tracks from data frames", {
  # Create test interaction data
  test_data <- data.frame(
    seqnames1 = "chr1",
    start1 = c(1000, 2000, 3000),
    end1 = c(1100, 2100, 3100),
    seqnames2 = "chr1",
    start2 = c(5000, 6000, 7000),
    end2 = c(5100, 6100, 7100),
    score = c(0.5, 0.8, 0.3)
  )

  # Test basic arc plot
  p1 <- ez_arc(test_data, "chr1:1000-8000")
  expect_s3_class(p1, "ggplot")
  expect_true(length(p1$layers) > 0)

  # Test with score coloring
  p2 <- ez_arc(test_data, "chr1:1000-8000", use_score = TRUE)
  expect_s3_class(p2, "ggplot")
  expect_true(length(p2$layers) > 0)

  # Test with custom parameters
  p3 <- ez_arc(test_data, "chr1:1000-8000",
               curvature = 0.3, color = "red", size = 1, alpha = 0.8)
  expect_s3_class(p3, "ggplot")
  expect_true(length(p3$layers) > 0)

  # Test error handling
  expect_error(ez_arc("not_a_dataframe", "chr1:1000-8000"))
})

test_that("ez_hic creates Hi-C tracks from data frames", {
  # Create test Hi-C data
  test_data <- data.frame(
    pos1 = rep(seq(1000, 5000, 1000), 5),
    pos2 = rep(seq(1000, 5000, 1000), each = 5),
    count = runif(25, 0, 100)
  )

  # Test basic Hi-C plot
  p1 <- ez_hic(test_data, "chr1:1000-5000")
  expect_s3_class(p1, "ggplot")
  expect_true(length(p1$layers) > 0)

  # Test with custom parameters
  p2 <- ez_hic(test_data, "chr1:1000-5000",
               resolution = 5000, log_transform = FALSE,
               low = "blue", high = "red")
  expect_s3_class(p2, "ggplot")
  expect_true(length(p2$layers) > 0)

  # Test error handling
  expect_error(ez_hic("not_a_dataframe", "chr1:1000-5000"))
})

test_that("ez_coverage handles complex input scenarios", {
  # Test with grouping
  test_data <- data.frame(
    seqnames = "chr1",
    start = rep(seq(1000, 2000, 100), 2),
    end = rep(seq(1100, 2100, 100), 2),
    score = rnorm(22),
    sample = rep(c("Sample1", "Sample2"), each = 11)
  )

  # Test with grouping
  p1 <- ez_coverage(test_data, "chr1:1000-2000", group_var = "sample")
  expect_s3_class(p1, "ggplot")
  expect_true(length(p1$layers) > 0)

  # Test with multiple tracks (list input)
  track_list <- list(
    "Track1" = test_data[test_data$sample == "Sample1", ],
    "Track2" = test_data[test_data$sample == "Sample2", ]
  )

  p2 <- ez_coverage(track_list, "chr1:1000-2000")
  expect_s3_class(p2, "ggplot")
  expect_true(length(p2$layers) > 0)

  # Test with different visualization types
  p3 <- ez_coverage(test_data, "chr1:1000-2000", type = "line")
  expect_s3_class(p3, "ggplot")
  expect_true(length(p3$layers) > 0)

  p4 <- ez_coverage(test_data, "chr1:1000-2000", type = "heatmap")
  expect_s3_class(p4, "ggplot")
  expect_true(length(p4$layers) > 0)

  # Test with binning
  p5 <- ez_coverage(test_data, "chr1:1000-2000", bin_width = 50L)
  expect_s3_class(p5, "ggplot")
  expect_true(length(p5$layers) > 0)

  # Test with custom colors
  p6 <- ez_coverage(test_data, "chr1:1000-2000", group_var = "sample",
                  colors = c("red", "blue"))
  expect_s3_class(p6, "ggplot")
  expect_true(length(p6$layers) > 0)

  # Test error handling
  expect_error(ez_coverage(test_data, "chr1:1000-2000", alpha = 2)) # Invalid alpha
  expect_error(ez_coverage(test_data, "chr1:1000-2000", type = "invalid")) # Invalid type
})

test_that("ez_gene handles different input types", {
  # Test with data frame
  test_data <- data.frame(
    seqnames = "chr1",
    start = c(1000, 1200, 1400),
    end = c(1100, 1300, 1500),
    strand = c("+", "+", "-"),
    gene_id = c("gene1", "gene1", "gene2"),
    gene_name = c("Gene1", "Gene1", "Gene2"),
    type = c("gene", "exon", "gene")
  )

  p1 <- ez_gene(test_data, "chr1:1000-2000")
  expect_s3_class(p1, "ggplot")
  expect_true(length(p1$layers) > 0)

  # Test with custom parameters
  p2 <- ez_gene(test_data, "chr1:1000-2000",
                exon_height = 0.5, intron_width = 0.2,
                exon_color = "red", exon_fill = "blue")
  expect_s3_class(p2, "ggplot")
  expect_true(length(p2$layers) > 0)

  # Test error handling
  expect_error(ez_gene("not_a_dataframe", "chr1:1000-2000"))
})
