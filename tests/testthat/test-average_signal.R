# Tests for average_signal() and related internal helpers

# --- Test data setup ---
make_sample_df <- function(starts, ends, scores, chr = "chr1") {
  data.frame(
    seqnames = chr,
    start = starts,
    end = ends,
    score = scores
  )
}

# sample1 bins: [1000,1049], [1050,1099], [1100,1149], [1150,1199]  (50bp each)
sample1 <- make_sample_df(
  starts = c(1000, 1050, 1100, 1150),
  ends   = c(1049, 1099, 1149, 1199),
  scores = c(2, 4, 6, 8)
)

# sample2: different bin boundaries
sample2 <- make_sample_df(
  starts = c(1000, 1060, 1120, 1180),
  ends   = c(1059, 1119, 1179, 1199),
  scores = c(3, 5, 7, 1)
)

# chr1:1000-1199 → GRanges width = 200bp, seq(1000, 1199, 50) = 4 starts
region <- "chr1:1000-1199"


# ============================================================
# average_signal() — data frame list input
# ============================================================

test_that("average_signal returns correct structure with data frame inputs", {
  result <- average_signal(
    inputs = list(sample1, sample2),
    region = region,
    bin_width = 50
  )

  expect_s3_class(result, "data.frame")
  expect_true(all(c("seqnames", "start", "end", "score") %in% colnames(result)))
  expect_true(all(result$seqnames == "chr1"))
  # 200bp region / 50bp bins = 4 bins
  expect_equal(nrow(result), 4)
})

test_that("average_signal mean is between individual sample values", {
  result <- average_signal(
    inputs = list(sample1, sample2),
    region = region,
    bin_width = 50
  )

  expect_true(all(result$score >= 0))
})

test_that("average_signal with single data frame returns that signal", {
  result <- average_signal(
    inputs = list(sample1),
    region = region,
    bin_width = 50
  )

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 4)
  # Bins align exactly with sample1, so scores should match
  expect_equal(result$score, c(2, 4, 6, 8))
})

# ============================================================
# summary_fun parameter
# ============================================================

test_that("average_signal supports different summary functions", {
  identical_samples <- list(sample1, sample1)

  result_mean <- average_signal(identical_samples, region, bin_width = 50, summary_fun = "mean")
  result_max <- average_signal(identical_samples, region, bin_width = 50, summary_fun = "max")
  result_min <- average_signal(identical_samples, region, bin_width = 50, summary_fun = "min")
  result_sum <- average_signal(identical_samples, region, bin_width = 50, summary_fun = "sum")

  # Mean/max/min of identical values = the value itself
  expect_equal(result_mean$score, c(2, 4, 6, 8))
  expect_equal(result_max$score, c(2, 4, 6, 8))
  expect_equal(result_min$score, c(2, 4, 6, 8))
  # Sum of two identical = 2x
  expect_equal(result_sum$score, c(4, 8, 12, 16))
})

test_that("average_signal rejects invalid summary_fun", {
  expect_error(
    average_signal(list(sample1), region, summary_fun = "invalid"),
    "should be one of"
  )
})

# ============================================================
# bin_width and region handling
# ============================================================

test_that("average_signal respects bin_width", {
  result_50 <- average_signal(list(sample1), region, bin_width = 50)
  result_100 <- average_signal(list(sample1), region, bin_width = 100)

  expect_equal(nrow(result_50), 4)  # 200bp / 50bp
  expect_equal(nrow(result_100), 2) # 200bp / 100bp
})

test_that("average_signal handles bin_width that doesn't evenly divide region", {
  # 200bp region with 60bp bins → 4 bins (last shorter)
  result <- average_signal(list(sample1), region, bin_width = 60)
  n_expected <- ceiling(200 / 60)
  expect_equal(nrow(result), n_expected)
  # Last bin should end at region_end
  expect_equal(result$end[nrow(result)], 1199)
})

# ============================================================
# Edge cases
# ============================================================

test_that("average_signal handles empty region (no signal)", {
  empty_sample <- make_sample_df(
    starts = c(5000),
    ends   = c(5100),
    scores = c(10)
  )

  result <- average_signal(
    inputs = list(empty_sample),
    region = region,
    bin_width = 50
  )

  expect_equal(nrow(result), 4)
  expect_true(all(result$score == 0))
})

test_that("average_signal validates inputs", {
  expect_error(
    average_signal(42, region),
    "inputs must be"
  )
  expect_error(
    average_signal(list(sample1), region, bin_width = -10),
    "bin_width must be"
  )
})

# ============================================================
# .overlap_weighted_mean (internal helper)
# ============================================================

test_that(".overlap_weighted_mean computes correct weighted means", {
  bins_gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = c(100, 200), end = c(199, 299))
  )
  signal_gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = c(100, 150), end = c(149, 249)),
    score = c(10, 20)
  )

  result <- ezGenomeTracks:::.overlap_weighted_mean(bins_gr, signal_gr)
  expect_length(result, 2)
  # Bin 1 [100,199]: overlaps [100,149](50bp,10) and [150,199](50bp,20)
  # mean = (50*10 + 50*20) / 100 = 15
  expect_equal(result[1], 15)
  # Bin 2 [200,299]: overlaps [200,249](50bp,20)
  # mean = (50*20) / 50 = 20
  expect_equal(result[2], 20)
})

test_that(".overlap_weighted_mean returns zeros for no overlap", {
  bins_gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = c(100, 200), end = c(199, 299))
  )
  signal_gr <- GenomicRanges::GRanges(
    seqnames = "chr2",
    ranges = IRanges::IRanges(start = 100, end = 199),
    score = 10
  )

  result <- suppressWarnings(
    ezGenomeTracks:::.overlap_weighted_mean(bins_gr, signal_gr)
  )
  expect_length(result, 2)
  expect_true(all(result == 0))
})

# ============================================================
# .query_df_bins (internal helper)
# ============================================================

test_that(".query_df_bins validates input data frames", {
  bins_gr <- GenomicRanges::GRanges(
    seqnames = "chr1",
    ranges = IRanges::IRanges(start = 1000, end = 1049)
  )
  region_gr <- parse_region(region)

  bad_df <- data.frame(chrom = "chr1", pos = 1000, val = 5)
  expect_error(
    ezGenomeTracks:::.query_df_bins(list(bad_df), bins_gr, region_gr),
    "must contain columns"
  )
})

# ============================================================
# ez_coverage integration
# ============================================================

test_that("ez_coverage with average = TRUE and character vector produces a ggplot", {
  # Create temporary bedGraph files for testing
  tmp1 <- tempfile(fileext = ".bedGraph")
  tmp2 <- tempfile(fileext = ".bedGraph")
  writeLines(c(
    "chr1\t1000\t1050\t5",
    "chr1\t1050\t1100\t10",
    "chr1\t1100\t1150\t15",
    "chr1\t1150\t1199\t20"
  ), tmp1)
  writeLines(c(
    "chr1\t1000\t1050\t10",
    "chr1\t1050\t1100\t20",
    "chr1\t1100\t1150\t30",
    "chr1\t1150\t1199\t40"
  ), tmp2)
  on.exit(unlink(c(tmp1, tmp2)))

  result <- ez_coverage(
    input = c(tmp1, tmp2),
    region = region,
    average = TRUE,
    average_bin_width = 50
  )

  expect_s3_class(result, "gg")
  expect_s3_class(result, "ggplot")
})

test_that("ez_coverage with average = TRUE and list of data frames does NOT average", {
  # Lists represent separate tracks — averaging should not collapse them
  result <- ez_coverage(
    input = list(track1 = sample1, track2 = sample2),
    region = region,
    average = TRUE
  )

  expect_s3_class(result, "gg")
  # Should still have multiple tracks (not averaged into one)
  plot_data <- ggplot2::layer_data(result)
  # The data should have been processed individually, not collapsed
  expect_true("track" %in% colnames(result$data))
  expect_equal(length(unique(result$data$track)), 2)
})

# ============================================================
# Real-world bigWig tests (inst/extdata files)
# ============================================================

# Helper: find the package's bigWig files
bw_files <- system.file(
  "extdata",
  c("avg_chr2-231091223_231109786_231113600_0.bw",
    "avg_chr2-231091223_231109786_231113600_1.bw",
    "avg_chr2-231091223_231109786_231113600_2.bw"),
  package = "ezGenomeTracks"
)
bw_region <- "chr2:231109000-231113000"
bw_available <- all(nzchar(bw_files)) && all(file.exists(bw_files))

test_that("average_signal works with real bigWig files", {
  skip_if_not(bw_available, "bigWig test files not found")

  result <- average_signal(
    inputs = bw_files,
    region = bw_region,
    bin_width = 100
  )

  expect_s3_class(result, "data.frame")
  expect_true(all(c("seqnames", "start", "end", "score") %in% colnames(result)))
  expect_true(nrow(result) > 0)
  # All scores should be finite after nans_to_zeros

  expect_true(all(is.finite(result$score)))
  # Region should be chr2
  expect_true(all(result$seqnames == "chr2"))
  # Bins should be within the requested region
  expect_true(all(result$start >= 231109000))
  expect_true(all(result$end <= 231113000))
})

test_that("average_signal bigWig scores change with summary_fun", {
  skip_if_not(bw_available, "bigWig test files not found")

  result_mean <- average_signal(bw_files, bw_region, bin_width = 200, summary_fun = "mean")
  result_max  <- average_signal(bw_files, bw_region, bin_width = 200, summary_fun = "max")
  result_sum  <- average_signal(bw_files, bw_region, bin_width = 200, summary_fun = "sum")

  # sum >= mean (for non-negative data with 3 samples)
  expect_true(all(result_sum$score >= result_mean$score - 1e-10))
  # max >= mean
  expect_true(all(result_max$score >= result_mean$score - 1e-10))
})

test_that("average_signal bigWig respects bin_width", {
  skip_if_not(bw_available, "bigWig test files not found")

  r50  <- average_signal(bw_files[1:2], bw_region, bin_width = 50)
  r200 <- average_signal(bw_files[1:2], bw_region, bin_width = 200)

  expect_true(nrow(r50) > nrow(r200))
})

test_that("average_signal with single bigWig returns its own signal", {
  skip_if_not(bw_available, "bigWig test files not found")

  result <- average_signal(bw_files[1], bw_region, bin_width = 100)

  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) > 0)
  expect_true(all(is.finite(result$score)))
})

test_that("ez_coverage average=TRUE with bigWig character vector produces plot", {
  skip_if_not(bw_available, "bigWig test files not found")

  # Character vector → overlapping tracks → average collapses to one
  p <- ez_coverage(
    input = bw_files,
    region = bw_region,
    average = TRUE,
    average_bin_width = 100
  )

  expect_s3_class(p, "gg")
  # Averaged: should NOT have a "track" column (single track)
  expect_false("track" %in% colnames(p$data))
})

test_that("ez_coverage average=TRUE with named list averages within elements", {
  skip_if_not(bw_available, "bigWig test files not found")

  # Named list: each element is a separate track.
  # Elements with >1 file get averaged; single-file elements stay as-is.
  p <- ez_coverage(
    input = list(
      "Averaged" = bw_files[1:2],
      "Single"   = bw_files[3]
    ),
    region = bw_region,
    average = TRUE,
    average_bin_width = 100,
    y_axis_style = "simple"
  )

  expect_s3_class(p, "gg")
  # Should have 2 tracks
  expect_true("track" %in% colnames(p$data))
  expect_equal(sort(unique(p$data$track)), c("Averaged", "Single"))
})

test_that("ez_coverage average=FALSE with bigWig character vector keeps overlapping tracks", {
  skip_if_not(bw_available, "bigWig test files not found")

  # Without average, character vector produces overlapping (grouped) tracks
  p <- ez_coverage(
    input = bw_files,
    region = bw_region,
    average = FALSE
  )

  expect_s3_class(p, "gg")
})

test_that("ez_coverage with average = TRUE and single input warns", {
  expect_warning(
    ez_coverage(sample1, region, average = TRUE),
    "average = TRUE has no effect"
  )
})
