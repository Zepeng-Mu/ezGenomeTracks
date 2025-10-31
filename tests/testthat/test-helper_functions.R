# Tests for helper functions

test_that("import_genomic_data imports files correctly", {
  # Test with a simple BED file (we'll create one temporarily)
  bed_content <- "chr1\t1000\t2000\tpeak1\t100\t+\nchr1\t3000\t4000\tpeak2\t200\t-"
  temp_bed <- tempfile(fileext = ".bed")
  writeLines(bed_content, temp_bed)

  # Test import
  result <- import_genomic_data(temp_bed)
  expect_s3_class(result, "data.frame")
  expect_true(all(c("seqnames", "start", "end") %in% colnames(result)))
  expect_equal(nrow(result), 2)

  # Clean up
  unlink(temp_bed)
})

test_that("get_single_signal extracts data correctly", {
  # Create test data frame
  test_data <- data.frame(
    seqnames = c("chr1", "chr1", "chr2"),
    start = c(1000, 1500, 2000),
    end = c(1200, 1700, 2200),
    score = c(0.5, 0.8, 0.3)
  )

  # Test with data frame input
  result1 <- get_single_signal(test_data, "chr1:1000-2000", name = "test_track")
  expect_s3_class(result1, "data.frame")
  expect_true("name" %in% colnames(result1))
  expect_equal(result1$name[1], "test_track")
  expect_equal(nrow(result1), 2) # Only chr1 features

  # Test with file input (create temporary file)
  bed_content <- "chr1\t1000\t2000\tpeak1\t100\t+\nchr1\t3000\t4000\tpeak2\t200\t-"
  temp_bed <- tempfile(fileext = ".bed")
  writeLines(bed_content, temp_bed)

  result2 <- get_single_signal(temp_bed, "chr1:1000-2000", name = "file_track")
  expect_s3_class(result2, "data.frame")
  expect_true("name" %in% colnames(result2))
  expect_equal(result2$name[1], "file_track")

  # Clean up
  unlink(temp_bed)
})

test_that("process_signal_input handles different input types", {
  # Test data frame input
  test_df <- data.frame(
    seqnames = "chr1",
    start = seq(1000, 2000, 100),
    end = seq(1100, 2100, 100),
    score = rnorm(11)
  )

  result1 <- process_signal_input(test_df, "chr1:1000-2000")
  expect_s3_class(result1, "data.frame")
  expect_true(all(c("seqnames", "start", "end", "score") %in% colnames(result1)))

  # Test character vector input (create temporary files)
  bed_content1 <- "chr1\t1000\t2000\tpeak1\t100\t+\nchr1\t3000\t4000\tpeak2\t200\t-"
  bed_content2 <- "chr1\t1500\t2500\tpeak3\t150\t+\nchr1\t3500\t4500\tpeak4\t250\t-"

  temp_bed1 <- tempfile(fileext = ".bed")
  temp_bed2 <- tempfile(fileext = ".bed")
  writeLines(bed_content1, temp_bed1)
  writeLines(bed_content2, temp_bed2)

  result2 <- process_signal_input(c(temp_bed1, temp_bed2), "chr1:1000-4000",
                                 track_labels = c("Sample1", "Sample2"))
  expect_s3_class(result2, "data.frame")
  expect_true("track" %in% colnames(result2))
  expect_true("group" %in% colnames(result2))

  # Test list input
  result3 <- process_signal_input(
    list("Track1" = test_df, "Track2" = c(temp_bed1, temp_bed2)),
    "chr1:1000-4000"
  )
  expect_s3_class(result3, "data.frame")
  expect_true("track" %in% colnames(result3))

  # Clean up
  unlink(c(temp_bed1, temp_bed2))
})

test_that("plot_signal_df creates signal plots with grouping", {
  # Create test data with grouping
  test_data <- data.frame(
    start = rep(seq(1000, 2000, 100), 2),
    end = rep(seq(1100, 2100, 100), 2),
    score = rnorm(22),
    sample = rep(c("Sample1", "Sample2"), each = 11),
    condition = rep(c("Control", "Treatment"), each = 11)
  )

  # Test basic plot
  p1 <- plot_signal_df(test_data, "chr1:1000-2000")
  expect_s3_class(p1, "ggplot")
  expect_true(length(p1$layers) > 0)

  # Test with track separation
  p2 <- plot_signal_df(test_data, "chr1:1000-2000", track_by = "sample")
  expect_s3_class(p2, "ggplot")
  expect_true(length(p2$layers) > 0)

  # Test with grouping
  p3 <- plot_signal_df(test_data, "chr1:1000-2000",
                       track_by = "sample", group_by = "condition")
  expect_s3_class(p3, "ggplot")
  expect_true(length(p3$layers) > 0)

  # Test with custom colors
  p4 <- plot_signal_df(test_data, "chr1:1000-2000",
                       track_by = "sample", group_by = "condition",
                       colors = c("red", "blue"))
  expect_s3_class(p4, "ggplot")
  expect_true(length(p4$layers) > 0)

  # Test with binning
  p5 <- plot_signal_df(test_data, "chr1:1000-2000", binwidth = 50)
  expect_s3_class(p5, "ggplot")
  expect_true(length(p5$layers) > 0)
})

test_that("process_gene_data processes gene annotation correctly", {
  # Create test GRanges data
  gr <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr1", "chr1"),
    ranges = IRanges::IRanges(start = c(1000, 1200, 1400), end = c(1100, 1300, 1500)),
    strand = c("+", "+", "-"),
    gene_id = c("gene1", "gene1", "gene2"),
    gene_name = c("Gene1", "Gene1", "Gene2"),
    type = c("gene", "exon", "gene")
  )

  # Test processing
  result <- process_gene_data(gr)
  expect_s3_class(result, "data.frame")
  expect_true(all(c("gene_id", "gene_name", "type", "xstart", "xend", "y") %in% colnames(result)))

  # Test with data frame input
  gr_df <- granges_to_df(gr)
  result2 <- process_gene_data(gr_df)
  expect_s3_class(result2, "data.frame")
  expect_true(all(c("gene_id", "gene_name", "type", "xstart", "xend", "y") %in% colnames(result2)))
})

test_that("extract_txdb_data extracts data from TxDb objects", {
  # Skip if TxDb.Hsapiens.UCSC.hg19.knownGene is not available
  skip_if_not_installed("TxDb.Hsapiens.UCSC.hg19.knownGene")

  # Load TxDb package
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

  # Test extraction
  region_gr <- parse_region("chr1:1000000-2000000")
  result <- extract_txdb_data(txdb, region_gr)

  expect_s3_class(result, "data.frame")
  expect_true(all(c("gene_id", "gene_name", "type", "xstart", "xend", "y") %in% colnames(result)))
})

test_that("error handling works correctly", {
  # Test invalid input types
  expect_error(import_genomic_data("nonexistent_file.bed"))
  expect_error(get_single_signal("invalid_input", "chr1:1000-2000"))
  expect_error(process_signal_input("invalid_input", "chr1:1000-2000"))

  # Test missing required columns
  bad_data <- data.frame(start = 1:10, end = 11:20) # Missing seqnames and score
  expect_error(process_signal_input(bad_data, "chr1:1000-2000"))

  # Test invalid region format
  expect_error(parse_region("invalid_region"))
  expect_error(parse_region("chr1:1000")) # Missing end
})
