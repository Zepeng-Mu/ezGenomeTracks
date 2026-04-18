# ezGenomeTracks - Helper functions (split from helpers.R)
#' Process Manhattan plot input into standardized data frame
#'
#' This function converts any input type (data.frame or named list) into a
#' standardized data frame with consistent columns for Manhattan plot visualization.
#' Supports both GWAS-style (CHR, BP, P) and GRanges-style (seqnames, start, pvalue)
#' column naming conventions with auto-detection.
#'
#' @param input Either a data frame or named list of data frames
#' @param chr Column name for chromosome. If NULL, auto-detects from common names
#'   (CHR, chr, seqnames, chrom, chromosome). Default: NULL
#' @param bp Column name for base pair position. If NULL, auto-detects from common names
#'   (BP, bp, start, pos, position, POS). Default: NULL
#' @param p Column name for p-value. If NULL, auto-detects from common names
#'   (P, p, pvalue, p.value, pval). Default: NULL
#' @param snp Column name for SNP identifier. If NULL, auto-detects from common names
#'   (SNP, snp, rsid, id, variant_id, marker). Default: NULL (optional column)
#' @param track_labels Optional vector of track labels (used for unnamed list input)
#' @return A data frame with standardized columns and optionally track column
#' @export
#' @importFrom dplyr bind_rows mutate
#' @examples
#' # Data frame input with GWAS-style columns
#' df <- data.frame(CHR = 1, BP = 1:100, P = runif(100), SNP = paste0("rs", 1:100))
#' process_manhattan_input(df)
#'
#' # Data frame input with GRanges-style columns
#' df2 <- data.frame(seqnames = "chr1", start = 1:100, pvalue = runif(100))
#' process_manhattan_input(df2)
#'
#' # List input
#' data_list <- list("GWAS1" = df1, "GWAS2" = df2)
#' process_manhattan_input(data_list)
process_manhattan_input <- function(
  input,
  chr = NULL,
  bp = NULL,
  p = NULL,
  snp = NULL,
  track_labels = NULL
) {
  # Helper function to auto-detect column names
  detect_column <- function(data, candidates, param_name, required = TRUE) {
    for (col in candidates) {
      if (col %in% colnames(data)) return(col)
    }
    if (required) {
      stop(paste0(
        "Could not find ",
        param_name,
        " column. Expected one of: ",
        paste(candidates, collapse = ", ")
      ))
    }
    return(NULL)
  }

  # Define candidate column names
  chr_candidates <- c("CHR", "chr", "seqnames", "chrom", "chromosome")
  bp_candidates <- c("BP", "bp", "start", "pos", "position", "POS")
  p_candidates <- c("P", "p", "pvalue", "p.value", "pval", "P.value")
  snp_candidates <- c("SNP", "snp", "rsid", "id", "variant_id", "marker")

  if (is.data.frame(input)) {
    # Case 1: Data frame input
    # Auto-detect column names if not provided
    if (is.null(chr)) {
      chr <- detect_column(input, chr_candidates, "chromosome")
    }
    if (is.null(bp)) {
      bp <- detect_column(input, bp_candidates, "position")
    }
    if (is.null(p)) {
      p <- detect_column(input, p_candidates, "p-value")
    }
    if (is.null(snp)) {
      snp <- detect_column(input, snp_candidates, "SNP", required = FALSE)
    }

    return(input)
  } else if (is.list(input)) {
    # Case 2: List input (named or unnamed)
    if (is.null(names(input)) && is.null(track_labels)) {
      names(input) <- paste0("Track ", seq_along(input))
    } else if (is.null(names(input)) && !is.null(track_labels)) {
      names(input) <- track_labels
    }

    track_data_list <- list()
    for (i in seq_along(input)) {
      track_name <- names(input)[i]
      track_element <- input[[i]]

      if (is.data.frame(track_element)) {
        # Auto-detect column names for each track element if not provided
        local_chr <- if (is.null(chr)) {
          detect_column(track_element, chr_candidates, "chromosome")
        } else {
          chr
        }
        local_bp <- if (is.null(bp)) {
          detect_column(track_element, bp_candidates, "position")
        } else {
          bp
        }
        local_p <- if (is.null(p)) {
          detect_column(track_element, p_candidates, "p-value")
        } else {
          p
        }

        # Add track column
        track_element$track <- track_name
        track_data_list[[i]] <- track_element
      } else {
        stop("List elements must be data frames")
      }
    }

    names(track_data_list) <- names(input)
    return(dplyr::bind_rows(track_data_list))
  } else {
    stop("Input must be a data frame or named list of data frames")
  }
}
