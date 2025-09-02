# Tests for granges_helpers.R functions

test_that("granges_to_df converts GRanges to data frame", {
  # Create a simple GRanges object
  gr <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1", "chr2"),
    ranges = IRanges::IRanges(start = c(1, 100, 200), end = c(50, 150, 250)),
    strand = c("+", "-", "*"),
    score = c(10, 20, 30),
    name = c("A", "B", "C")
  )
  
  # Convert to data frame
  df <- granges_to_df(gr)
  
  # Check structure
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 3)
  expect_true(all(c("seqnames", "start", "end", "strand", "score", "name") %in% colnames(df)))
  
  # Check content
  expect_equal(df$seqnames, c("chr1", "chr1", "chr2"))
  expect_equal(df$start, c(1, 100, 200))
  expect_equal(df$end, c(50, 150, 250))
  expect_equal(df$strand, c("+", "-", "*"))
  expect_equal(df$score, c(10, 20, 30))
  expect_equal(df$name, c("A", "B", "C"))
})

test_that("df_to_granges converts data frame to GRanges", {
  # Create a simple data frame
  df <- data.frame(
    seqnames = c("chr1", "chr1", "chr2"),
    start = c(1, 100, 200),
    end = c(50, 150, 250),
    strand = c("+", "-", "*"),
    score = c(10, 20, 30),
    name = c("A", "B", "C"),
    stringsAsFactors = FALSE
  )
  
  # Convert to GRanges
  gr <- df_to_granges(df)
  
  # Check structure
  expect_s4_class(gr, "GRanges")
  expect_equal(length(gr), 3)
  expect_true(all(c("score", "name") %in% names(GenomicRanges::mcols(gr))))
  
  # Check content
  expect_equal(as.character(GenomicRanges::seqnames(gr)), c("chr1", "chr1", "chr2"))
  expect_equal(GenomicRanges::start(gr), c(1, 100, 200))
  expect_equal(GenomicRanges::end(gr), c(50, 150, 250))
  expect_equal(as.character(GenomicRanges::strand(gr)), c("+", "-", "*"))
  expect_equal(gr$score, c(10, 20, 30))
  expect_equal(gr$name, c("A", "B", "C"))
})

test_that("parse_region parses genomic region string", {
  # Test valid region strings
  gr1 <- parse_region("chr1:1000-2000")
  expect_s4_class(gr1, "GRanges")
  expect_equal(as.character(GenomicRanges::seqnames(gr1)), "chr1")
  expect_equal(GenomicRanges::start(gr1), 1000)
  expect_equal(GenomicRanges::end(gr1), 2000)
  
  gr2 <- parse_region("chrX:1,000,000-2,000,000")
  expect_s4_class(gr2, "GRanges")
  expect_equal(as.character(GenomicRanges::seqnames(gr2)), "chrX")
  expect_equal(GenomicRanges::start(gr2), 1000000)
  expect_equal(GenomicRanges::end(gr2), 2000000)
  
  # Test invalid region strings
  expect_error(parse_region("invalid"))
  expect_error(parse_region("chr1:1000"))
  expect_error(parse_region("chr1-1000-2000"))
})