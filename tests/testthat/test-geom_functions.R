# Tests for geom functions

test_that("geom_coverage creates different visualization types", {
  # Create test data
  test_data <- data.frame(
    start = seq(1, 100, 10),
    end = seq(10, 100, 10),
    score = rnorm(10)
  )

  # Test area plot (default)
  p1 <- ggplot2::ggplot(test_data) + geom_coverage()
  expect_s3_class(p1, "ggplot")
  expect_true(length(p1$layers) > 0)

  # Test line plot
  p2 <- ggplot2::ggplot(test_data) + geom_coverage(type = "line")
  expect_s3_class(p2, "ggplot")
  expect_true(length(p2$layers) > 0)

  # Test heatmap
  p3 <- ggplot2::ggplot(test_data) + geom_coverage(type = "heatmap")
  expect_s3_class(p3, "ggplot")
  expect_true(length(p3$layers) > 0)

  # Test with custom mapping
  p4 <- ggplot2::ggplot(test_data) +
    geom_coverage(ggplot2::aes(x = start, y = score, color = score))
  expect_s3_class(p4, "ggplot")
  expect_true(length(p4$layers) > 0)
})

test_that("geom_feature creates feature tracks", {
  # Create test data
  test_data <- data.frame(
    start = c(10, 30, 50),
    end = c(20, 40, 60),
    score = c(0.5, 0.8, 0.3)
  )

  # Test basic feature plot
  p1 <- ggplot2::ggplot(test_data) + geom_feature()
  expect_s3_class(p1, "ggplot")
  expect_true(length(p1$layers) > 0)

  # Test with custom mapping
  p2 <- ggplot2::ggplot(test_data) +
    geom_feature(ggplot2::aes(xmin = start, xmax = end, fill = score))
  expect_s3_class(p2, "ggplot")
  expect_true(length(p2$layers) > 0)

  # Test with custom parameters
  p3 <- ggplot2::ggplot(test_data) +
    geom_feature(height = 0.5, color = "red", fill = "blue", alpha = 0.8)
  expect_s3_class(p3, "ggplot")
  expect_true(length(p3$layers) > 0)
})

test_that("geom_gene creates gene tracks", {
  # Create test gene data
  test_data <- data.frame(
    xstart = c(10, 20, 30),
    xend = c(15, 25, 35),
    type = c("gene", "exon", "gene"),
    strand = c("+", "+", "-"),
    gene_name = c("Gene1", "Gene1", "Gene2")
  )

  # Test basic gene plot
  p1 <- ggplot2::ggplot(test_data) + geom_gene()
  expect_s3_class(p1, "ggplot")
  expect_true(length(p1$layers) > 0)

  # Test with strand separation
  p2 <- ggplot2::ggplot(test_data) +
    geom_gene(ggplot2::aes(xstart = xstart, xend = xend, type = type, strand = strand))
  expect_s3_class(p2, "ggplot")
  expect_true(length(p2$layers) > 0)

  # Test with custom parameters
  p3 <- ggplot2::ggplot(test_data) +
    geom_gene(exon_height = 0.5, intron_width = 0.2)
  expect_s3_class(p3, "ggplot")
  expect_true(length(p3$layers) > 0)
})

test_that("geom_manhattan creates Manhattan plots", {
  # Create test GWAS data
  test_data <- data.frame(
    CHR = rep(1:3, each = 10),
    BP = rep(seq(1000, 10000, 1000), 3),
    P = runif(30, 1e-8, 1),
    SNP = paste0("rs", 1:30)
  )

  # Test basic Manhattan plot
  p1 <- ggplot2::ggplot(test_data) +
    geom_manhattan(ggplot2::aes(CHR = CHR, BP = BP, P = P))
  expect_s3_class(p1, "ggplot")
  expect_true(length(p1$layers) > 0)

  # Test with highlighting
  highlight_data <- test_data[test_data$P < 1e-5, ]
  p2 <- ggplot2::ggplot(test_data) +
    geom_manhattan(ggplot2::aes(CHR = CHR, BP = BP, P = P),
                   highlight_snps = highlight_data)
  expect_s3_class(p2, "ggplot")
  expect_true(length(p2$layers) > 0)

  # Test with threshold line
  p3 <- ggplot2::ggplot(test_data) +
    geom_manhattan(ggplot2::aes(CHR = CHR, BP = BP, P = P),
                   threshold_p = 5e-8)
  expect_s3_class(p3, "ggplot")
  expect_true(length(p3$layers) > 0)
})

test_that("geom_arc creates interaction tracks", {
  # Create test interaction data
  test_data <- data.frame(
    start1 = c(10, 20, 30),
    start2 = c(50, 60, 70),
    score = c(0.5, 0.8, 0.3)
  )

  # Test basic arc plot
  p1 <- ggplot2::ggplot(test_data) +
    geom_arc(ggplot2::aes(x = start1, y = 0, xend = start2, yend = 0))
  expect_s3_class(p1, "ggplot")
  expect_true(length(p1$layers) > 0)

  # Test with score coloring
  p2 <- ggplot2::ggplot(test_data) +
    geom_arc(ggplot2::aes(x = start1, y = 0, xend = start2, yend = 0, color = score))
  expect_s3_class(p2, "ggplot")
  expect_true(length(p2$layers) > 0)

  # Test with custom parameters
  p3 <- ggplot2::ggplot(test_data) +
    geom_arc(ggplot2::aes(x = start1, y = 0, xend = start2, yend = 0),
             curvature = 0.3, size = 1, alpha = 0.8)
  expect_s3_class(p3, "ggplot")
  expect_true(length(p3$layers) > 0)
})

test_that("geom_hic creates Hi-C tracks", {
  # Create test Hi-C data
  test_data <- data.frame(
    pos1 = rep(seq(1000, 5000, 1000), 5),
    pos2 = rep(seq(1000, 5000, 1000), each = 5),
    count = runif(25, 0, 100)
  )

  # Test basic Hi-C plot
  p1 <- ggplot2::ggplot(test_data) +
    geom_hic(ggplot2::aes(x = pos1, y = pos2, fill = count))
  expect_s3_class(p1, "ggplot")
  expect_true(length(p1$layers) > 0)

  # Test with custom colors
  p2 <- ggplot2::ggplot(test_data) +
    geom_hic(ggplot2::aes(x = pos1, y = pos2, fill = count),
             low = "blue", high = "red")
  expect_s3_class(p2, "ggplot")
  expect_true(length(p2$layers) > 0)
})

test_that("stat_bin_signal bins signal data", {
  # Create test signal data
  test_data <- data.frame(
    start = seq(1, 1000, 1),
    score = rnorm(1000)
  )

  # Test with binwidth
  p1 <- ggplot2::ggplot(test_data) +
    stat_bin_signal(ggplot2::aes(x = start, y = score), binwidth = 50)
  expect_s3_class(p1, "ggplot")
  expect_true(length(p1$layers) > 0)

  # Test with number of bins
  p2 <- ggplot2::ggplot(test_data) +
    stat_bin_signal(ggplot2::aes(x = start, y = score), bins = 20)
  expect_s3_class(p2, "ggplot")
  expect_true(length(p2$layers) > 0)
})
