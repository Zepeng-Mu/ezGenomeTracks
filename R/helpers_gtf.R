# ezGenomeTracks - GTF parsing helper functions
# Internal utilities for import_gtf and related functions

#' Extract GTF attribute values
#'
#' Internal helper to safely extract attribute values from a GTF attributes data frame.
#' When rtracklayer imports a GTF, it extracts common attributes into mcols automatically.
#' This function retrieves these values with a sensible fallback.
#'
#' @param attributes_df Data frame of attributes extracted by rtracklayer
#' @param col_name Attribute name to extract (e.g., "gene_id", "gene_name")
#' @param default Default value if attribute not found (recycled to match length)
#'
#' @return Character vector of extracted values
#'
#' @keywords internal
#' @noRd
.extract_gtf_attribute <- function(attributes_df, col_name, default = NA_character_) {
  if (col_name %in% names(attributes_df)) {
    vals <- attributes_df[[col_name]]
    # Replace NA/empty with default
    if (length(default) == 1) {
      vals[is.na(vals) | vals == ""] <- default
    } else {
      mask <- is.na(vals) | vals == ""
      vals[mask] <- default[mask]
    }
    return(as.character(vals))
  }

  # Column not found; return defaults
  return(rep_len(as.character(default), nrow(attributes_df)))
}

#' Validate GTF data frame format
#'
#' Check that a data frame has required columns for use with ez_* functions.
#' Raises informative errors if columns are missing or have incorrect types.
#'
#' @param df Data frame to validate
#' @param required_cols Character vector of required column names
#'
#' @return TRUE invisibly if valid, otherwise stops with error message
#'
#' @keywords internal
#' @noRd
validate_gtf_df <- function(
  df,
  required_cols = c("seqnames", "start", "end", "strand", "type", "gene_id")
) {
  if (!is.data.frame(df)) {
    stop("Input must be a data frame. Got: ", class(df)[1])
  }

  # Check for required columns
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(
      "GTF data frame missing required columns:\n  ",
      paste(missing_cols, collapse = ", "),
      "\nAvailable columns: ", paste(names(df), collapse = ", ")
    )
  }

  # Check data types and attempt coercion if needed
  if (!is.numeric(df$start) || !is.numeric(df$end)) {
    stop("Columns 'start' and 'end' must be numeric")
  }

  if (!is.character(df$strand)) {
    df$strand <- as.character(df$strand)
  }

  if (!is.character(df$type)) {
    df$type <- as.character(df$type)
  }

  return(TRUE)
}

#' Clean and standardize GTF data frame
#'
#' Ensures a GTF data frame has consistent formatting and required columns
#' for use with ez_* functions. Performs normalization like strand conversion
#' and type coercion.
#'
#' @param df GTF data frame (e.g., output from import_gtf)
#' @param gene_id_col Name of gene_id column (default: "gene_id")
#' @param gene_name_col Name of gene_name column (default: "gene_name")
#'
#' @return Cleaned data frame with standardized columns
#'
#' @keywords internal
#' @noRd
.standardize_gtf_df <- function(
  df,
  gene_id_col = "gene_id",
  gene_name_col = "gene_name"
) {
  # Validate structure
  validate_gtf_df(df)

  # Ensure strand values are standard: "+", "-", or "*"
  df$strand <- as.character(df$strand)
  valid_strands <- c("+", "-", "*")
  df$strand[!df$strand %in% valid_strands] <- "*"

  # Ensure type column has consistent casing
  df$type <- tolower(df$type)

  # Ensure gene_name defaults to gene_id if missing
  if (is.null(df[[gene_name_col]]) || all(is.na(df[[gene_name_col]]))) {
    df[[gene_name_col]] <- df[[gene_id_col]]
  }

  return(df)
}
