test_that("validate_gtf_df validates data frames", {
  # Create a valid GTF-like data frame
  valid_df <- data.frame(
    seqnames = "chr1",
    start = 1000,
    end = 2000,
    strand = "+",
    type = "gene",
    gene_id = "GENE1",
    stringsAsFactors = FALSE
  )

  # Should pass validation
  expect_true(validate_gtf_df(valid_df))

  # Create invalid data frame (missing required columns)
  invalid_df <- data.frame(
    seqnames = "chr1",
    start = 1000,
    end = 2000
  )

  expect_error(validate_gtf_df(invalid_df), "required columns")
})

test_that("extract_gtf_attribute extracts attribute values", {
  test_attrs <- data.frame(
    gene_id = c("GENE1", "GENE2", NA),
    gene_name = c("GeneA", NA, "GeneC"),
    stringsAsFactors = FALSE
  )

  # Extract with defaults
  result <- .extract_gtf_attribute(test_attrs, "gene_id", default = "UNKNOWN")
  expect_equal(result[1], "GENE1")
  expect_equal(result[3], "UNKNOWN")

  # Extract non-existent column
  result <- .extract_gtf_attribute(test_attrs, "missing", default = "DEFAULT")
  expect_equal(length(result), 3)
  expect_true(all(result == "DEFAULT"))
})

test_that("import_gtf function exists and is exported", {
  # Verify the function exists
  expect_true(exists("import_gtf"))
  expect_true("import_gtf" %in% getNamespaceExports("ezGenomeTracks"))
})

test_that("ez_gene auto_import parameter exists", {
  # Check that ez_gene has the auto_import parameter
  params <- names(formals(ez_gene))
  expect_true("auto_import" %in% params)
})
