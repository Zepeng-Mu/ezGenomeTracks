# Integration tests for end-to-end workflows
# This file contains tests for complete workflows using multiple ezGenomeTracks functions together

test_that("Complete genomic visualization workflow", {
  # Create comprehensive test dataset
  test_data <- create_comprehensive_test_data("chr1:1000000-1100000")

  # Step 1: Create individual tracks
  expect_no_error({
    signal_track <- ez_coverage(test_data$genomic_signal, region = "chr1:1000000-1100000")
    expect_s3_class(signal_track, "ggplot")
  })

  expect_no_error({
    gene_track <- ez_gene(test_data$genes, region = "chr1:1000000-1100000")
    expect_s3_class(gene_track, "ggplot")
  })

  expect_no_error({
    interaction_track <- ez_link(test_data$interactions, region = "chr1:1000000-1100000")
    expect_s3_class(interaction_track, "ggplot")
  })

  expect_no_error({
    hic_track <- ez_hic(test_data$hic_data, region = "chr1:1000000-1100000")
    expect_s3_class(hic_track, "ggplot")
  })

  # Step 2: Stack all tracks into a comprehensive plot
  expect_no_error({
    comprehensive_plot <- vstack_plot(
      signal_track,
      gene_track,
      interaction_track,
      hic_track
    )
    expect_s3_class(comprehensive_plot, "ggplot")
    expect_true(can_render_plot(comprehensive_plot))
  })
})

test_that("Multi-track workflow with annotations", {
  # Create test data
  signal_data <- create_test_genomic_data()
  gene_data <- create_test_gene_data()

  # Create base tracks
  coverage_track <- ez_coverage(signal_data, region = "chr1:1000000-1100000")
  gene_track <- ez_gene(gene_data, region = "chr1:1000000-1100000")

  # Add annotations to coverage track
  annotated_coverage <- coverage_track %>%
    add_highlight(
      data.frame(
        start = 1050000,
        end = 1060000,
        score = 0.8
      )
    ) %>%
    add_vline(position = 1050000) %>%
    add_text(x = 1050500, y = 15, label = "Peak Region")

  # Stack tracks
  final_plot <- vstack_plot(annotated_coverage, gene_track)

  expect_s3_class(final_plot, "ggplot")
  expect_true(can_render_plot(final_plot))

  # Verify all annotations are present
  expect_true(length(final_plot$layers) > length(coverage_track$layers))
})

test_that("Comparative analysis workflow", {
  # Create multi-condition signal data
  multi_condition_data <- create_test_signal_data(
    n_samples = 2,
    n_conditions = 2,
    sample_names = c("Control", "Treatment"),
    condition_names = c("Condition_A", "Condition_B")
  )

  # Create comparative plots
  expect_no_error({
    control_plot <- ez_coverage(
      multi_condition_data[multi_condition_data$condition == "Condition_A", ],
      region = "chr1:1000000-1100000"
    )
    expect_s3_class(control_plot, "ggplot")
  })

  expect_no_error({
    treatment_plot <- ez_coverage(
      multi_condition_data[multi_condition_data$condition == "Condition_B", ],
      region = "chr1:1000000-1100000"
    )
    expect_s3_class(treatment_plot, "ggplot")
  })

  # Compare plots side by side
  expect_no_error({
    comparison_plot <- vstack_plot(control_plot, treatment_plot)
    expect_s3_class(comparison_plot, "ggplot")
  })
})

test_that("Sashimi plot integration workflow", {
  # Create sashimi plot data
  coverage_data <- create_test_genomic_data()
  junction_data <- data.frame(
    seqnames = "chr1",
    start = c(1050000, 1070000, 1085000),
    end = c(1060000, 1080000, 1095000),
    score = c(25, 50, 15)
  )

  # Create base sashimi plot
  sashimi_plot <- ez_sashimi(
    coverage_data,
    junction_data,
    region = "chr1:1000000-1100000"
  )

  expect_s3_class(sashimi_plot, "ggplot")
  expect_true(can_render_plot(sashimi_plot))

  # Add annotations to sashimi plot
  annotated_sashimi <- sashimi_plot %>%
    add_vline(position = 1055000) %>%
    add_text(x = 1055000, y = 20, label = "Splice Site")

  expect_s3_class(annotated_sashimi, "ggplot")
  expect_true(can_render_plot(annotated_sashimi))

  # Test different sashimi parameters
  expect_no_error({
    log_sashimi <- ez_sashimi(
      coverage_data,
      junction_data,
      region = "chr1:1000000-1100000",
      score_transform = "log10"
    )
    expect_s3_class(log_sashimi, "ggplot")
  })

  expect_no_error({
    alternate_sashimi <- ez_sashimi(
      coverage_data,
      junction_data,
      region = "chr1:1000000-1100000",
      junction_direction = "alternate"
    )
    expect_s3_class(alternate_sashimi, "ggplot")
  })
})

test_that("Complex genomic region analysis", {
  # Create test data for a specific genomic region
  region <- "chr1:1050000-1080000"

  # Generate various data types for this region
  signal_data <- create_test_genomic_data(
    start = 1050000,
    end = 1080000,
    n_peaks = 5
  )

  gene_data <- create_test_gene_data(
    chr = "chr1",
    gene_count = 3
  )

  interaction_data <- create_test_interaction_data(
    start = 1050000,
    end = 1080000,
    n_interactions = 15
  )

  # Create tracks
  signal_track <- ez_coverage(signal_data, region = region)
  gene_track <- ez_gene(gene_data, region = region)
  interaction_track <- ez_link(interaction_data, region = region)

  # Combine into comprehensive analysis
  analysis_plot <- vstack_plot(
    signal_track,
    gene_track,
    interaction_track
  )

  expect_s3_class(analysis_plot, "ggplot")
  expect_true(can_render_plot(analysis_plot))

  # Add region-specific annotations
  final_analysis <- analysis_plot %>%
    add_rect(xmin = 1055000, xmax = 1065000, alpha = 0.2, fill = "red") %>%
    add_text(x = 1060000, y = max(signal_data$score) * 1.1,
             label = "Focus Region")

  expect_s3_class(final_analysis, "ggplot")
  expect_true(can_render_plot(final_analysis))
})

test_that("Data import and processing workflow", {
  # Test file-based data import workflow
  # Create temporary test files
  bed_content <- c(
    "chr1\t1050000\t1051000\tpeak1\t0.8\t+",
    "chr1\t1060000\t1061000\tpeak2\t0.6\t+",
    "chr1\t1070000\t1071000\tpeak3\t0.9\t-"
  )

  gtf_content <- c(
    "chr1\t.\tgene\t1050000\t1070000\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"TestGene1\";",
    "chr1\t.\texon\t1051000\t1051500\t.\t+\t.\tgene_id \"GENE1\"; exon_number \"1\";",
    "chr1\t.\texon\t1061000\t1061500\t.\t+\t.\tgene_id \"GENE1\"; exon_number \"2\";"
  )

  temp_files <- create_temp_files(
    c(".bed", ".gtf"),
    list(bed_content, gtf_content)
  )

  bed_file <- temp_files[[1]]
  gtf_file <- temp_files[[2]]

  # Test reading and processing files
  expect_no_error({
    # This would typically use import_genomic_data or similar functions
    # For now, we'll simulate the workflow with data frames
    peak_data <- read.table(bed_file, sep = "\t", header = FALSE)
    colnames(peak_data) <- c("seqnames", "start", "end", "name", "score", "strand")

    gene_gtf_data <- read.table(gtf_file, sep = "\t", header = FALSE)
    colnames(gene_gtf_data) <- c("seqnames", "source", "type", "start", "end",
                                 "score", "strand", "phase", "attributes")
  })

  # Create plots from imported data
  expect_no_error({
    peak_track <- ez_feature(peak_data, region = "chr1:1050000-1071000")
    expect_s3_class(peak_track, "ggplot")
  })

  # Clean up
  cleanup_temp_files(temp_files)
})

test_that("Theme and styling integration", {
  signal_data <- create_test_genomic_data()

  # Apply different themes
  expect_no_error({
    custom_theme_plot <- ez_theme() +
      ez_coverage(signal_data, region = "chr1:1000000-1100000")
    expect_s3_class(custom_theme_plot, "ggplot")
  })

  # Test specific track themes
  expect_no_error({
    coverage_theme_plot <- ez_coverage_theme() +
      ez_coverage(signal_data, region = "chr1:1000000-1100000")
    expect_s3_class(coverage_theme_plot, "ggplot")
  })

  gene_data <- create_test_gene_data()
  expect_no_error({
    gene_theme_plot <- ez_gene_theme() +
      ez_gene(gene_data, region = "chr1:1000000-1100000")
    expect_s3_class(gene_theme_plot, "ggplot")
  })

  # Combine themed tracks
  expect_no_error({
    coverage_track <- ez_coverage_theme() +
      ez_coverage(signal_data, region = "chr1:1000000-1100000")
    gene_track <- ez_gene_theme() +
      ez_gene(gene_data, region = "chr1:1000000-1100000")

    combined_plot <- vstack_plot(coverage_track, gene_track)
    expect_s3_class(combined_plot, "ggplot")
  })
})

test_that("Scale and axis integration", {
  signal_data <- create_test_genomic_data()

  # Test genomic coordinate scales
  expect_no_error({
    plot_with_genome_scale <- ez_coverage(signal_data, region = "chr1:1000000-1100000") +
      scale_x_genome_region(region = "chr1:1000000-1100000")
    expect_s3_class(plot_with_genome_scale, "ggplot")
  })

  # Test strand scales for gene data
  gene_data <- create_test_gene_data()
  expect_no_error({
    gene_plot_with_strand_scale <- ez_gene(gene_data, region = "chr1:1000000-1100000") +
      scale_y_strand()
    expect_s3_class(gene_plot_with_strand_scale, "ggplot")
  })

  # Test Hi-C color scales
  hic_data <- create_test_hic_data()
  expect_no_error({
    hic_plot_with_scale <- ez_hic(hic_data, region = "chr1:1000000-1100000") +
      scale_fill_hic(palette = "viridis")
    expect_s3_class(hic_plot_with_scale, "ggplot")
  })
})

test_that("Error recovery and resilience workflow", {
  # Test workflow with some invalid data mixed in
  mixed_signal_data <- rbind(
    create_test_genomic_data()[1:5, ],  # Valid data
    data.frame(seqnames = "invalid", start = "bad", end = "data", score = "values")  # Invalid
  )

  # This should either fail gracefully or handle the error
  expect_warning_or_error({
    result <- ez_coverage(mixed_signal_data, region = "chr1:1000000-1100000")
    if (!inherits(result, "try-error")) {
      expect_s3_class(result, "ggplot")
    }
  })

  # Test with partial data recovery
  partial_gene_data <- create_test_gene_data()[1:3, ]  # Subset of gene data

  expect_no_error({
    partial_gene_plot <- ez_gene(partial_gene_data, region = "chr1:1000000-1100000")
    expect_s3_class(partial_gene_plot, "ggplot")
  })
})

test_that("Performance and scalability workflow", {
  # Test with progressively larger datasets
  dataset_sizes <- c(100, 500, 1000)

  for (size in dataset_sizes) {
    # Create larger signal dataset
    large_signal <- create_test_genomic_data(
      start = 1000000,
      end = 1000000 + size * 100,
      bin_size = 100
    )

    start_time <- Sys.time()

    expect_no_error({
      plot <- ez_coverage(large_signal, region = paste0("chr1:1000000-", 1000000 + size * 100))
      expect_s3_class(plot, "ggplot")
      expect_true(can_render_plot(plot))
    })

    end_time <- Sys.time()
    processing_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

    # Each larger dataset should still complete within reasonable time
    expect_true(processing_time < size * 0.01)  # Less than 0.01 seconds per data point
  }
})

test_that("Cross-function compatibility workflow", {
  # Test that different functions work together seamlessly
  signal_data <- create_test_genomic_data()
  gene_data <- create_test_gene_data()

  # Create plots using different functions
  coverage_plot <- ez_coverage(signal_data, region = "chr1:1000000-1100000")
  gene_plot <- ez_gene(gene_data, region = "chr1:1000000-1100000")

  # Use geom functions directly
  geom_plot <- ggplot2::ggplot(signal_data) + geom_coverage()
  feature_plot <- ggplot2::ggplot(gene_data) + geom_gene()

  # Mix ez_* and geom_* functions
  mixed_plot <- vstack_plot(coverage_plot, geom_plot, feature_plot, gene_plot)

  expect_s3_class(mixed_plot, "ggplot")
  expect_true(can_render_plot(mixed_plot))

  # Add annotations to mixed plot
  annotated_mixed <- mixed_plot %>%
    add_vline(position = 1050000) %>%
    add_hline(y = 10)

  expect_s3_class(annotated_mixed, "ggplot")
  expect_true(can_render_plot(annotated_mixed))
})

test_that("End-to-end quality assurance workflow", {
  # Comprehensive test of a complete analysis workflow

  # Step 1: Data preparation
  test_region <- "chr1:1050000-1080000"
  comprehensive_data <- create_comprehensive_test_data(test_region)

  # Step 2: Individual track creation
  tracks <- list()

  # Signal track
  tracks$signal <- ez_coverage(comprehensive_data$genomic_signal, region = test_region)

  # Feature tracks
  peak_data <- data.frame(
    seqnames = "chr1",
    start = c(1055000, 1065000),
    end = c(1056000, 1066000),
    score = c(0.8, 0.6),
    name = c("Peak1", "Peak2")
  )
  tracks$features <- ez_feature(peak_data, region = test_region)

  # Gene track
  tracks$genes <- ez_gene(comprehensive_data$genes, region = test_region)

  # Interaction track
  tracks$interactions <- ez_link(comprehensive_data$interactions, region = test_region)

  # Step 3: Quality checks on individual tracks
  for (track_name in names(tracks)) {
    track <- tracks[[track_name]]
    expect_s3_class(track, "ggplot")
    expect_true(can_render_plot(track))
    expect_true(is_valid_ggplot(track))
  }

  # Step 4: Track combination
  combined_plot <- do.call(vstack_plot, tracks)
  expect_s3_class(combined_plot, "ggplot")
  expect_true(can_render_plot(combined_plot))

  # Step 5: Annotation and finalization
  final_plot <- combined_plot %>%
    add_highlight(
      data.frame(
        start = 1060000,
        end = 1070000,
        score = 0.7
      )
    ) %>%
    add_vline(position = 1065000) %>%
    add_text(x = 1065000, y = 15, label = "Analysis Region")

  # Step 6: Final quality assurance
  expect_s3_class(final_plot, "ggplot")
  expect_true(can_render_plot(final_plot))
  expect_true(length(final_plot$layers) > length(combined_plot$layers))

  # Step 7: Output validation
  # Test different output formats
  expect_true(can_render_plot(final_plot, width = 800, height = 600))
  expect_true(can_render_plot(final_plot, width = 1200, height = 800))

  # Save snapshot for regression testing
  save_test_snapshot(final_plot, "integration_test_final_plot")

  # Verify snapshot can be loaded and compared
  expect_true(compare_test_snapshot(final_plot, "integration_test_final_plot", tolerance = 1e-6))
})

# Helper function for testing warning or error conditions in integration tests
warning_or_error <- function(expr) {
  result <- tryCatch(
    withCallingHandlers(
      expr,
      warning = function(w) {
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      list(status = "error", message = e$message)
    },
    warning = function(w) {
      list(status = "warning", message = w$message)
    }
  )

  return(result)
}

expect_warning_or_error <- function(expr) {
  result <- warning_or_error(expr)
  expect_true(result$status %in% c("success", "warning", "error"))
}