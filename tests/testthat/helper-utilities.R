#' @title Testing Utilities for ezGenomeTracks Package
#' @description Helper functions for common testing patterns and utilities
#' @name test-utilities
#' @keywords internal testing

#' Temporary File Management Utilities
#' @description Functions to manage temporary files for testing

#' Create a temporary file with specific extension
#' @param extension File extension including dot (e.g., ".bed", ".gtf")
#' @param content Optional content to write to file
#' @return Path to temporary file
#' @export
#' @examples
#' temp_file <- create_temp_file(".bed", "chr1\t100\t200\tpeak1")
create_temp_file <- function(extension = ".txt", content = NULL) {
  temp_file <- tempfile(pattern = "test_", fileext = extension)

  if (!is.null(content)) {
    writeLines(content, temp_file)
  }

  return(temp_file)
}

#' Create multiple temporary files
#' @param extensions Vector of file extensions
#' @param contents Optional list of content for each file
#' @return List of temporary file paths
#' @export
#' @examples
#' temp_files <- create_temp_files(c(".bed", ".gtf"),
#'                                list("chr1\t100\t200", "chr1\t.\tgene\t100\t200"))
create_temp_files <- function(extensions, contents = NULL) {
  if (length(extensions) != length(contents) && !is.null(contents)) {
    stop("Length of extensions and contents must match")
  }

  temp_files <- vector("list", length(extensions))

  for (i in seq_along(extensions)) {
    if (is.null(contents)) {
      temp_files[[i]] <- create_temp_file(extensions[i])
    } else {
      temp_files[[i]] <- create_temp_file(extensions[i], contents[[i]])
    }
  }

  return(temp_files)
}

#' Clean up temporary files
#' @param file_paths Vector of file paths to remove
#' @export
#' @examples
#' temp_files <- create_temp_files(c(".bed"))
#' cleanup_temp_files(temp_files)
cleanup_temp_files <- function(file_paths) {
  invisible(file.remove(unlist(file_paths)))
}

#' Plot Testing Utilities
#' @description Functions to test ggplot objects and plots

#' Check if object is a valid ggplot object
#' @param plot_obj Object to check
#' @return logical indicating if object is a ggplot
#' @export
#' @examples
#' p <- ggplot2::ggplot() + ggplot2::geom_blank()
#' expect_true(is_valid_ggplot(p))
is_valid_ggplot <- function(plot_obj) {
  inherits(plot_obj, "gg") && !is.null(plot_obj$layers) && !is.null(plot_obj$data)
}

#' Check if ggplot has expected layers
#' @param plot_obj ggplot object
#' @param layer_classes Expected layer classes (character vector)
#' @return logical indicating if plot has expected layers
#' @export
#' @examples
#' p <- ggplot2::ggplot() + ggplot2::geom_line()
#' expect_true(has_expected_layers(p, "GeomLine"))
has_expected_layers <- function(plot_obj, layer_classes) {
  if (!is_valid_ggplot(plot_obj)) {
    return(FALSE)
  }

  actual_classes <- sapply(plot_obj$layers, function(layer) class(layer$geom)[1])
  return(all(layer_classes %in% actual_classes))
}

#' Check if ggplot has expected scales
#' @param plot_obj ggplot object
#' @param scale_types Expected scale types (e.g., "ScaleXContinuous", "ScaleYContinuous")
#' @return logical indicating if plot has expected scales
#' @export
#' @examples
#' p <- ggplot2::ggplot() + ggplot2::scale_x_continuous()
#' expect_true(has_expected_scales(p, "ScaleXContinuous"))
has_expected_scales <- function(plot_obj, scale_types) {
  if (!is_valid_ggplot(plot_obj)) {
    return(FALSE)
  }

  actual_classes <- sapply(plot_obj$scales$scales, function(scale) class(scale)[1])
  return(all(scale_types %in% actual_classes))
}

#' Test ggplot rendering without errors
#' @param plot_obj ggplot object
#' @param width Plot width (default: 800)
#' @param height Plot height (default: 600)
#' @return logical indicating if plot renders without error
#' @export
#' @examples
#' p <- ggplot2::ggplot() + ggplot2::geom_blank()
#' expect_true(can_render_plot(p))
can_render_plot <- function(plot_obj, width = 800, height = 600) {
  if (!is_valid_ggplot(plot_obj)) {
    return(FALSE)
  }

  # Try to render the plot
  result <- tryCatch({
    # Create a temporary device
    temp_file <- tempfile(fileext = ".png")
    grDevices::png(temp_file, width = width, height = height)
    print(plot_obj)
    grDevices::dev.off()

    # Check if file was created and has content
    file.exists(temp_file) && file.info(temp_file)$size > 0
  }, error = function(e) {
    FALSE
  }, warning = function(w) {
    # Warnings are acceptable, only fail on errors
    TRUE
  })

  return(result)
}

#' Data Validation Helpers
#' @description Functions to validate genomic data structures

#' Validate genomic data frame structure
#' @param df Data frame to validate
#' @param required_cols Required column names
#' @param check_types Whether to check data types
#' @return logical indicating if data frame is valid
#' @export
#' @examples
#' test_df <- data.frame(seqnames = "chr1", start = 100, end = 200)
#' expect_true(is_valid_genomic_df(test_df, c("seqnames", "start", "end")))
is_valid_genomic_df <- function(df, required_cols, check_types = TRUE) {
  # Check if it's a data frame
  if (!is.data.frame(df)) {
    return(FALSE)
  }

  # Check required columns
  if (!all(required_cols %in% colnames(df))) {
    return(FALSE)
  }

  if (check_types) {
    # Check numeric columns are actually numeric
    numeric_cols <- intersect(required_cols, c("start", "end", "score", "count"))
    for (col in numeric_cols) {
      if (!is.numeric(df[[col]])) {
        return(FALSE)
      }
    }

    # Check character columns are actually character
    char_cols <- intersect(required_cols, c("seqnames", "strand", "type", "gene_id", "gene_name"))
    for (col in char_cols) {
      if (!is.character(df[[col]])) {
        return(FALSE)
      }
    }
  }

  # Check for reasonable values
  if ("start" %in% required_cols && "end" %in% required_cols) {
    if (any(df$start >= df$end, na.rm = TRUE)) {
      return(FALSE)
    }
  }

  return(TRUE)
}

#' Validate genomic region string format
#' @param region_string Region string (e.g., "chr1:1000-2000")
#' @return logical indicating if region format is valid
#' @export
#' @examples
#' expect_true(is_valid_region_string("chr1:1000-2000"))
#' expect_false(is_valid_region_string("invalid_format"))
is_valid_region_string <- function(region_string) {
  # Check basic format
  if (!grepl("^[^:]+:[0-9]+-[0-9]+$", region_string)) {
    return(FALSE)
  }

  # Parse and validate coordinates
  parts <- strsplit(region_string, ":")[[1]]
  coords <- strsplit(parts[2], "-")[[1]]
  start <- as.integer(coords[1])
  end <- as.integer(coords[2])

  # Check coordinates are positive and end > start
  return(start > 0 && end > start)
}

#' Compare data frames with tolerance
#' @param df1 First data frame
#' @param df2 Second data frame
#' @param tolerance Numeric tolerance for floating point comparison
#' @param check_order Whether to check row order
#' @return logical indicating if data frames are equivalent
#' @export
#' @examples
#' df1 <- data.frame(a = c(1.0, 2.0), b = c("x", "y"))
#' df2 <- data.frame(a = c(1.001, 2.001), b = c("x", "y"))
#' expect_true(are_equivalent_dfs(df1, df2, tolerance = 0.01))
are_equivalent_dfs <- function(df1, df2, tolerance = 1e-10, check_order = FALSE) {
  # Check dimensions
  if (ncol(df1) != ncol(df2) || nrow(df1) != nrow(df2)) {
    return(FALSE)
  }

  # Check column names
  if (!identical(colnames(df1), colnames(df2))) {
    return(FALSE)
  }

  # Compare each column
  for (col in colnames(df1)) {
    if (is.numeric(df1[[col]]) && is.numeric(df2[[col]])) {
      # Numeric comparison with tolerance
      diff <- abs(df1[[col]] - df2[[col]])
      if (any(diff > tolerance, na.rm = TRUE)) {
        return(FALSE)
      }
    } else {
      # Character/factor comparison
      if (!identical(df1[[col]], df2[[col]])) {
        return(FALSE)
      }
    }
  }

  # Check row order if requested
  if (check_order) {
    return(identical(df1, df2))
  }

  return(TRUE)
}

#' Common Testing Patterns
#' @description Pre-defined test patterns for common scenarios

#' Test function with various input types
#' @param func Function to test
#' @param test_cases List of test cases with inputs and expected outputs
#' @param error_cases List of inputs expected to cause errors
#' @return List of test results
#' @export
#' @examples
#' test_cases <- list(
#'   list(input = "valid_input", expected = "expected_output"),
#'   list(input = "another_valid", expected = "another_expected")
#' )
#' error_cases <- list("invalid_input", NULL)
#' results <- test_function_variations(my_function, test_cases, error_cases)
test_function_variations <- function(func, test_cases, error_cases = NULL) {
  results <- list(success = 0, failures = list(), errors = list())

  # Test valid cases
  for (i in seq_along(test_cases)) {
    test_case <- test_cases[[i]]
    result <- tryCatch({
      actual <- do.call(func, test_case$input)
      expected <- test_case$expected

      if (isTRUE(all.equal(actual, expected))) {
        results$success <- results$success + 1
        list(status = "success", case = i)
      } else {
        results$failures[[length(results$failures) + 1]] <- list(
          status = "failure",
          case = i,
          input = test_case$input,
          expected = expected,
          actual = actual
        )
        list(status = "failure", case = i)
      }
    }, error = function(e) {
      results$errors[[length(results$errors) + 1]] <- list(
        status = "error",
        case = i,
        input = test_case$input,
        error = e$message
      )
      list(status = "error", case = i)
    })
  }

  # Test error cases
  if (!is.null(error_cases)) {
    for (i in seq_along(error_cases)) {
      error_input <- error_cases[[i]]
      result <- tryCatch({
        do.call(func, error_input)
        # If we get here, no error was thrown when expected
        results$failures[[length(results$failures) + 1]] <- list(
          status = "no_error",
          case = i,
          input = error_input
        )
        list(status = "no_error", case = i)
      }, error = function(e) {
        # Error was expected and thrown
        results$success <- results$success + 1
        list(status = "expected_error", case = i)
      })
    }
  }

  return(results)
}

#' Mock external dependencies for testing
#' @description Functions to mock external package dependencies
#' @param package Package name to mock
#' @param functions Named list of function names to mock
#' @param mock_functions Named list of mock function definitions
#' @export
#' @examples
#' mock_functions <- list(
#'   read.table = function(...) data.frame(V1 = "mock")
#' )
#' with_mocked_functions("utils", c("read.table"), mock_functions, {
#'   result <- read.table("fake_file.txt")
#'   expect_equal(result$V1, "mock")
#' })
with_mocked_functions <- function(package, functions, mock_functions, expr) {
  if (!requireNamespace("mockr", quietly = TRUE)) {
    stop("mockr package required for mocking functions")
  }

  mockr::with_mock(
    setNames(mock_functions, paste0(package, "::", functions)),
    eval(expr)
  )
}

#' Helper to skip tests when dependencies are not available
#' @param package Package name to check
#' @param message Optional message to display when skipping
#' @export
#' @examples
#' skip_if_missing_dependency("BiocManager")
skip_if_missing_dependency <- function(package, message = NULL) {
  if (!requireNamespace(package, quietly = TRUE)) {
    skip_msg <- if (is.null(message)) {
      sprintf("Package '%s' not available, skipping test", package)
    } else {
      message
    }
    testthat::skip(skip_msg)
  }
}

#' Helper to create reproducible test snapshots
#' @param data Data to snapshot
#' @param filename Snapshot filename (without extension)
#' @param path Path to save snapshots (default: "tests/testthat/_snaps")
#' @export
#' @examples
#' test_data <- create_test_genomic_data()
#' save_test_snapshot(test_data, "genomic_data")
save_test_snapshot <- function(data, filename, path = "tests/testthat/_snaps") {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }

  snapshot_file <- file.path(path, paste0(filename, ".rds"))
  saveRDS(data, snapshot_file)
}

#' Helper to load and compare test snapshots
#' @param data Current data
#' @param filename Snapshot filename (without extension)
#' @param path Path to load snapshots from (default: "tests/testthat/_snaps")
#' @param tolerance Numeric tolerance for comparison
#' @export
#' @examples
#' test_data <- create_test_genomic_data()
#' expect_true(compare_test_snapshot(test_data, "genomic_data"))
compare_test_snapshot <- function(data, filename, path = "tests/testthat/_snaps", tolerance = 1e-10) {
  snapshot_file <- file.path(path, paste0(filename, ".rds"))

  if (!file.exists(snapshot_file)) {
    warning(sprintf("Snapshot file %s does not exist", snapshot_file))
    return(FALSE)
  }

  saved_data <- readRDS(snapshot_file)

  if (is.data.frame(data) && is.data.frame(saved_data)) {
    return(are_equivalent_dfs(data, saved_data, tolerance = tolerance))
  } else {
    return(isTRUE(all.equal(data, saved_data, tolerance = tolerance)))
  }
}