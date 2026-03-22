#!/usr/bin/env Rscript

#' Coverage Analysis Script for ezGenomeTracks Package
#'
#' This script performs code coverage analysis using the covr package.
#' It runs tests, generates coverage reports, and checks coverage thresholds.
#'
#' Usage: Rscript tests/test-coverage.R [options]
#'
#' Options:
#'   --threshold=NUM    Minimum coverage threshold (default: 80)
#'   --output=DIR      Output directory for reports (default: coverage)
#'   --quiet           Suppress output except errors
#'   --help            Show this help message

# Load optparse if available, otherwise use basic argument parsing
if (requireNamespace("optparse", quietly = TRUE)) {
  library(optparse)
  use_optparse <- TRUE
} else {
  use_optparse <- FALSE
}

# Parse command line arguments
if (use_optparse) {
  option_list <- list(
    make_option(c("--threshold"), type = "numeric", default = 80,
                help = "Minimum coverage threshold percentage [default: %default]"),
    make_option(c("--output"), type = "character", default = "coverage",
                help = "Output directory for coverage reports [default: %default]"),
    make_option(c("--quiet"), action = "store_true", default = FALSE,
                help = "Suppress output except errors [default: %default]"),
    make_option(c("--help"), action = "store_true", default = FALSE,
                help = "Show this help message and exit")
  )

  opt_parser <- OptionParser(option_list = option_list)
  opt <- parse_args(opt_parser)

  # Show help if requested
  if (opt$help) {
    print_help(opt_parser)
    quit(status = 0)
  }
} else {
  # Basic argument parsing without optparse
  args <- commandArgs(trailingOnly = TRUE)

  opt <- list(
    threshold = 80,
    output = "coverage",
    quiet = FALSE,
    help = FALSE
  )

  # Parse simple arguments
  for (i in seq_along(args)) {
    if (args[i] == "--help" || args[i] == "-h") {
      opt$help <- TRUE
    } else if (args[i] == "--quiet") {
      opt$quiet <- TRUE
    } else if (grepl("^--threshold=", args[i])) {
      opt$threshold <- as.numeric(sub("^--threshold=", "", args[i]))
    } else if (grepl("^--output=", args[i])) {
      opt$output <- sub("^--output=", "", args[i])
    }
  }

  if (opt$help) {
    cat("Usage: Rscript tests/test-coverage.R [options]\n")
    cat("Options:\n")
    cat("  --threshold=NUM    Minimum coverage threshold (default: 80)\n")
    cat("  --output=DIR       Output directory for reports (default: coverage)\n")
    cat("  --quiet            Suppress output except errors\n")
    cat("  --help             Show this help message\n")
    quit(status = 0)
  }
}

# Create output directory if it doesn't exist
if (!dir.exists(opt$output)) {
  dir.create(opt$output, recursive = TRUE)
}

# Check if covr is available
if (!requireNamespace("covr", quietly = TRUE)) {
  stop("covr package is required but not installed. Please install with: install.packages('covr')")
}

# Load required packages
suppressPackageStartupMessages({
  library(covr)
  library(testthat)
})

# Function to generate coverage report
generate_coverage_report <- function(coverage_data, output_dir, package_name = "ezGenomeTracks") {

  # Create HTML report
  html_file <- file.path(output_dir, "coverage_report.html")
  report(coverage_data, html_file)

  # Create text summary
  text_file <- file.path(output_dir, "coverage_summary.txt")
  capture.output({
    cat("Coverage Analysis Report for", package_name, "\n")
    cat("================================================\n\n")
    print(coverage_data)
  }, file = text_file)

  # Create coverage data CSV
  csv_file <- file.path(output_dir, "coverage_data.csv")
  coverage_df <- as.data.frame(coverage_data)
  write.csv(coverage_df, csv_file, row.names = FALSE)

  # Return summary statistics
  list(
    total_lines = nrow(coverage_df),
    covered_lines = sum(coverage_df$coverage > 0),
    percent_covered = round(mean(coverage_df$coverage) * 100, 2),
    files = list(
      html = html_file,
      text = text_file,
      csv = csv_file
    )
  )
}

# Function to check coverage threshold
check_coverage_threshold <- function(coverage_data, threshold) {
  coverage_percent <- mean(coverage_data$coverage) * 100

  if (coverage_percent >= threshold) {
    return(list(
      passed = TRUE,
      coverage = coverage_percent,
      message = sprintf("Coverage %.2f%% meets threshold of %.0f%%", coverage_percent, threshold)
    ))
  } else {
    return(list(
      passed = FALSE,
      coverage = coverage_percent,
      message = sprintf("Coverage %.2f%% below threshold of %.0f%%", coverage_percent, threshold)
    ))
  }
}

# Main coverage analysis
main <- function() {

  if (!opt$quiet) {
    cat("Starting coverage analysis for ezGenomeTracks package...\n")
    cat(sprintf("Coverage threshold: %.0f%%\n", opt$threshold))
    cat(sprintf("Output directory: %s\n\n", opt$output))
  }

  # Run package tests and calculate coverage
  if (!opt$quiet) {
    cat("Running tests and calculating coverage...\n")
  }

  # Calculate coverage using different methods for robustness
  coverage_result <- tryCatch({
    # Try package coverage first
    package_coverage()
  }, error = function(e) {
    if (!opt$quiet) {
      cat("Package coverage failed, trying directory coverage...\n")
    }
    # Fallback to directory coverage
    tryCatch({
      coverage <- covr::package_coverage(quiet = TRUE)
      return(coverage)
    }, error = function(e2) {
      if (!opt$quiet) {
        cat("Trying local directory coverage...\n")
      }
      # Last resort: coverage of current directory
      covr::coverage_file()
    })
  })

  if (is.null(coverage_result) || length(coverage_result) == 0) {
    stop("Failed to calculate coverage. Please ensure tests are runnable.")
  }

  # Generate reports
  if (!opt$quiet) {
    cat("Generating coverage reports...\n")
  }

  report_summary <- generate_coverage_report(coverage_result, opt$output)

  # Check threshold
  threshold_check <- check_coverage_threshold(coverage_result, opt$threshold)

  # Print summary
  if (!opt$quiet) {
    cat("\nCoverage Analysis Summary:\n")
    cat("=========================\n")
    cat(sprintf("Total lines analyzed: %d\n", report_summary$total_lines))
    cat(sprintf("Lines with coverage: %d\n", report_summary$covered_lines))
    cat(sprintf("Overall coverage: %.2f%%\n", report_summary$percent_covered))
    cat(sprintf("Threshold check: %s\n", threshold_check$message))
    cat(sprintf("\nReports saved to: %s\n", opt$output))
    cat(sprintf("  HTML report: %s\n", report_summary$files$html))
    cat(sprintf("  Text summary: %s\n", report_summary$files$text))
    cat(sprintf("  CSV data: %s\n", report_summary$files$csv))
  }

  # Return status for CI integration
  if (threshold_check$passed) {
    if (!opt$quiet) {
      cat("\n✓ Coverage analysis PASSED\n")
    }
    return(0)
  } else {
    if (!opt$quiet) {
      cat("\n✗ Coverage analysis FAILED\n")
      cat("Consider adding more tests to improve coverage.\n")
    }
    return(1)
  }
}

# Run main function and handle errors
tryCatch({
  exit_code <- main()
  quit(status = exit_code)
}, error = function(e) {
  if (!opt$quiet) {
    cat("ERROR:", conditionMessage(e), "\n")
    cat("Coverage analysis failed.\n")
  }
  quit(status = 1)
})