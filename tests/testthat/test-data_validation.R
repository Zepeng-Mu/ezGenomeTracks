# Tests for comprehensive data validation
# This file contains tests for input validation across all ezGenomeTracks functions

test_that("Region string validation works correctly", {
  # Valid region strings
  expect_true(is_valid_region_string("chr1:1000-2000"))
  expect_true(is_valid_region_string("chrX:50000-60000"))
  expect_true(is_valid_region_string("chr2_frag:100-200"))

  # Invalid region strings
  expect_false(is_valid_region_string("invalid_format"))
  expect_false(is_valid_region_string("chr1:2000-1000"))  # end < start
  expect_false(is_valid_region_string("chr1:0-1000"))     # start <= 0
  expect_false(is_valid_region_string("chr1:1000"))       # missing end
  expect_false(is_valid_region_string("1000-2000"))       # missing chr
  expect_false(is_valid_region_string("chr1:abc-def"))    # non-numeric coords
  expect_false(is_valid_region_string(""))                # empty string
  expect_false(is_valid_region_string(NA))                # NA value
  expect_false(is_valid_region_string(NULL))              # NULL value
})

test_that("Genomic data frame validation works for signal data", {
  # Valid signal data
  valid_signal <- data.frame(
    seqnames = c("chr1", "chr1"),
    start = c(1000, 1100),
    end = c(1100, 1200),
    score = c(5.5, 8.2)
  )
  expect_true(is_valid_genomic_df(valid_signal, c("seqnames", "start", "end", "score")))

  # Missing required columns
  invalid_missing_cols <- data.frame(
    seqnames = c("chr1", "chr1"),
    start = c(1000, 1100)
  )
  expect_false(is_valid_genomic_df(invalid_missing_cols, c("seqnames", "start", "end", "score")))

  # Invalid data types
  invalid_types <- data.frame(
    seqnames = c("chr1", "chr1"),
    start = c("1000", "1100"),  # should be numeric
    end = c(1100, 1200),
    score = c(5.5, 8.2)
  )
  expect_false(is_valid_genomic_df(invalid_types, c("seqnames", "start", "end", "score")))

  # Invalid coordinate ranges
  invalid_coords <- data.frame(
    seqnames = c("chr1", "chr1"),
    start = c(1100, 1000),  # start >= end
    end = c(1000, 1100),
    score = c(5.5, 8.2)
  )
  expect_false(is_valid_genomic_df(invalid_coords, c("seqnames", "start", "end", "score")))

  # Empty data frame
  empty_df <- data.frame(
    seqnames = character(0),
    start = integer(0),
    end = integer(0),
    score = numeric(0)
  )
  expect_false(is_valid_genomic_df(empty_df, c("seqnames", "start", "end", "score")))
})

test_that("Genomic data frame validation works for gene data", {
  # Valid gene data
  valid_genes <- data.frame(
    seqnames = c("chr1", "chr1"),
    start = c(1000, 2000),
    end = c(1500, 2500),
    strand = c("+", "-"),
    type = c("gene", "exon"),
    gene_id = c("GENE1", "GENE1"),
    gene_name = c("Gene1", "Gene1")
  )
  expect_true(is_valid_genomic_df(valid_genes, c("seqnames", "start", "end", "strand", "type")))

  # Invalid strand values
  invalid_strand <- data.frame(
    seqnames = c("chr1", "chr1"),
    start = c(1000, 2000),
    end = c(1500, 2500),
    strand = c("invalid", "+"),  # invalid strand value
    type = c("gene", "exon")
  )
  expect_false(is_valid_genomic_df(invalid_strand, c("seqnames", "start", "end", "strand", "type")))
})

test_that("Genomic data frame validation works for interaction data", {
  # Valid interaction data
  valid_interactions <- data.frame(
    seqnames1 = c("chr1", "chr1"),
    start1 = c(1000, 2000),
    end1 = c(1100, 2100),
    seqnames2 = c("chr1", "chr1"),
    start2 = c(1500, 2500),
    end2 = c(1600, 2600),
    score = c(5.5, 8.2)
  )
  expect_true(is_valid_genomic_df(valid_interactions, c("seqnames1", "start1", "end1", "seqnames2", "start2", "end2", "score")))

  # Invalid interaction coordinates
  invalid_interaction <- data.frame(
    seqnames1 = c("chr1", "chr1"),
    start1 = c(2000, 1000),  # start1 >= end1
    end1 = c(1000, 1100),
    seqnames2 = c("chr1", "chr1"),
    start2 = c(1500, 2500),
    end2 = c(1600, 2600),
    score = c(5.5, 8.2)
  )
  expect_false(is_valid_genomic_df(invalid_interaction, c("seqnames1", "start1", "end1", "seqnames2", "start2", "end2", "score")))
})

test_that("ez_* wrapper functions validate inputs correctly", {
  # Test signal track validation
  signal_data <- create_test_genomic_data()

  # Valid inputs should work
  expect_no_error({
    p1 <- ez_coverage(signal_data, region = "chr1:1000000-1100000")
    expect_s3_class(p1, "ggplot")
  })

  # Invalid region format
  expect_error(ez_coverage(signal_data, region = "invalid_region"))
  expect_error(ez_coverage(signal_data, region = "chr1:2000-1000"))

  # Missing required columns
  bad_signal <- signal_data[, -which(names(signal_data) == "score")]
  expect_error(ez_coverage(bad_signal, region = "chr1:1000000-1100000"))

  # Invalid alpha values
  expect_error(ez_coverage(signal_data, region = "chr1:1000000-1100000", alpha = -0.5))
  expect_error(ez_coverage(signal_data, region = "chr1:1000000-1100000", alpha = 2))

  # Invalid type parameter
  expect_error(ez_coverage(signal_data, region = "chr1:1000000-1100000", type = "invalid_type"))
})

test_that("ez_feature validates feature data correctly", {
  feature_data <- data.frame(
    seqnames = c("chr1", "chr1"),
    start = c(1000, 2000),
    end = c(1100, 2100),
    score = c(0.5, 0.8)
  )

  # Valid inputs should work
  expect_no_error({
    p1 <- ez_feature(feature_data, region = "chr1:1000-2500")
    expect_s3_class(p1, "ggplot")
  })

  # Missing required columns
  bad_feature <- feature_data[, -which(names(feature_data) == "score")]
  expect_error(ez_feature(bad_feature, region = "chr1:1000-2500"))

  # Invalid region
  expect_error(ez_feature(feature_data, region = "invalid"))
})

test_that("ez_gene validates gene data correctly", {
  gene_data <- create_test_gene_data()

  # Valid inputs should work
  expect_no_error({
    p1 <- ez_gene(gene_data, region = "chr1:1000000-1100000")
    expect_s3_class(p1, "ggplot")
  })

  # Missing required columns
  bad_gene <- gene_data[, -which(names(gene_data) == "strand")]
  expect_error(ez_gene(bad_gene, region = "chr1:1000000-1100000"))

  # Invalid strand values
  bad_strand_gene <- gene_data
  bad_strand_gene$strand <- c("invalid", rep("+", nrow(bad_strand_gene) - 1))
  expect_error(ez_gene(bad_strand_gene, region = "chr1:1000000-1100000"))
})

test_that("ez_link validates interaction data correctly", {
  interaction_data <- create_test_interaction_data()

  # Valid inputs should work
  expect_no_error({
    p1 <- ez_link(interaction_data, region = "chr1:1000000-1100000")
    expect_s3_class(p1, "ggplot")
  })

  # Missing required columns
  bad_interaction <- interaction_data[, -which(names(interaction_data) == "score")]
  expect_error(ez_link(bad_interaction, region = "chr1:1000000-1100000"))

  # Invalid interaction coordinates
  bad_coords <- interaction_data
  bad_coords$start1 <- bad_coords$end1  # start1 >= end1
  expect_error(ez_link(bad_coords, region = "chr1:1000000-1100000"))
})

test_that("ez_hic validates Hi-C data correctly", {
  hic_data <- create_test_hic_data()

  # Valid inputs should work
  expect_no_error({
    p1 <- ez_hic(hic_data, region = "chr1:1000000-1100000")
    expect_s3_class(p1, "ggplot")
  })

  # Missing required columns
  bad_hic <- hic_data[, -which(names(hic_data) == "count")]
  expect_error(ez_hic(bad_hic, region = "chr1:1000000-1100000"))

  # Negative counts
  bad_counts <- hic_data
  bad_counts$count <- -1
  expect_error(ez_hic(bad_counts, region = "chr1:1000000-1100000"))
})

test_that("ez_sashimi validates sashimi plot inputs correctly", {
  coverage_data <- create_test_genomic_data()
  junction_data <- data.frame(
    seqnames = "chr1",
    start = c(1500, 2500),
    end = c(2000, 3000),
    score = c(25, 50)
  )

  # Valid inputs should work
  expect_no_error({
    p1 <- ez_sashimi(coverage_data, junction_data, region = "chr1:1000000-1100000")
    expect_s3_class(p1, "ggplot")
  })

  # Invalid junction direction
  expect_error(ez_sashimi(coverage_data, junction_data, region = "chr1:1000000-1100000",
                         junction_direction = "invalid"))

  # Invalid score transform
  expect_error(ez_sashimi(coverage_data, junction_data, region = "chr1:1000000-1100000",
                         score_transform = "invalid"))

  # Invalid alpha values
  expect_error(ez_sashimi(coverage_data, junction_data, region = "chr1:1000000-1100000",
                         alpha = -0.5))

  # Missing junction data columns
  bad_junction <- junction_data[, -which(names(junction_data) == "score")]
  expect_error(ez_sashimi(coverage_data, bad_junction, region = "chr1:1000000-1100000"))
})

test_that("geom_* functions validate their inputs correctly", {
  # Test geom_coverage validation
  test_data <- data.frame(
    start = seq(1, 100, 10),
    end = seq(10, 100, 10),
    score = rnorm(10)
  )

  # Valid inputs should work
  expect_no_error({
    p1 <- ggplot2::ggplot(test_data) + geom_coverage()
    expect_s3_class(p1, "ggplot")
  })

  # Missing required aesthetics
  expect_error(ggplot2::ggplot(test_data) + geom_coverage(ggplot2::aes(x = start)))

  # Invalid type parameter
  expect_error(ggplot2::ggplot(test_data) + geom_coverage(type = "invalid"))
})

test_that("Helper functions validate parameters correctly", {
  # Test create_temp_file validation
  expect_error(create_temp_file("invalid extension"))
  expect_error(create_temp_file(NULL))

  # Test is_valid_ggplot validation
  expect_false(is_valid_ggplot(NULL))
  expect_false(is_valid_ggplot("not_a_plot"))
  expect_false(is_valid_ggplot(data.frame()))

  # Test are_equivalent_dfs validation
  df1 <- data.frame(a = 1:3, b = c("x", "y", "z"))
  df2 <- data.frame(a = 1:3, b = c("x", "y", "z"))
  df3 <- data.frame(a = 1:4, b = c("w", "x", "y", "z"))  # Different dimensions

  expect_true(are_equivalent_dfs(df1, df2))
  expect_false(are_equivalent_dfs(df1, df3))
  expect_false(are_equivalent_dfs(df1, NULL))
  expect_false(are_equivalent_dfs(NULL, df2))
})

test_that("File I/O validation works correctly", {
  # Test with temporary files
  temp_files <- create_temp_files(c(".bed", ".gtf"), list(
    "chr1\t100\t200\tpeak1",
    "chr1\t.\tgene\t100\t200\t.\t+\t.\tgene_id \"GENE1\";"
  ))

  expect_true(all(file.exists(unlist(temp_files))))

  # Clean up
  cleanup_temp_files(temp_files)
  expect_false(any(file.exists(unlist(temp_files))))

  # Test with invalid file extensions
  expect_error(create_temp_file("invalid_extension"))

  # Test with mismatched lengths
  expect_error(create_temp_files(c(".bed", ".gtf"), list("content1")))
})

test_that("Edge cases for genomic coordinates are handled correctly", {
  # Test very large coordinates
  large_coords <- data.frame(
    seqnames = "chr1",
    start = 1e9,
    end = 1e9 + 1000,
    score = 5.0
  )
  expect_no_error({
    p <- ez_coverage(large_coords, region = "chr1:1000000000-1000001100")
    expect_s3_class(p, "ggplot")
  })

  # Test zero/negative coordinates
  invalid_coords <- data.frame(
    seqnames = "chr1",
    start = 0,
    end = 100,
    score = 5.0
  )
  expect_error(ez_coverage(invalid_coords, region = "chr1:0-100"))

  # Test start > end
  reversed_coords <- data.frame(
    seqnames = "chr1",
    start = 200,
    end = 100,
    score = 5.0
  )
  expect_error(ez_coverage(reversed_coords, region = "chr1:100-200"))
})

test_that("NA and NULL values are handled correctly", {
  # Test data with NA values
  na_data <- data.frame(
    seqnames = c("chr1", NA, "chr1"),
    start = c(100, 200, 300),
    end = c(110, 210, 310),
    score = c(5.0, NA, 8.0)
  )

  expect_error(ez_coverage(na_data, region = "chr1:100-310"))

  # Test NULL data
  expect_error(ez_coverage(NULL, region = "chr1:100-200"))

  # Test empty character vectors
  empty_data <- data.frame(
    seqnames = character(0),
    start = integer(0),
    end = integer(0),
    score = numeric(0)
  )
  expect_error(ez_coverage(empty_data, region = "chr1:100-200"))
})

test_that("Data type coercion and validation works", {
  # Test automatic type coercion
  mixed_types <- data.frame(
    seqnames = c("chr1", "chr1"),  # Will be coerced to character
    start = c("100", "200"),       # Should be numeric but is character
    end = c("110", "210"),         # Should be numeric but is character
    score = c(5.0, 8.0)             # Numeric
  )

  # This should fail due to type mismatch
  expect_error(ez_coverage(mixed_types, region = "chr1:100-210"))

  # Test factor vs character
  factor_data <- data.frame(
    seqnames = factor(c("chr1", "chr1")),
    start = c(100, 200),
    end = c(110, 210),
    score = c(5.0, 8.0)
  )

  # Should work as factors are often converted to characters
  expect_no_error({
    p <- ez_coverage(factor_data, region = "chr1:100-210")
    expect_s3_class(p, "ggplot")
  })
})