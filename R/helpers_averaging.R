# ezGenomeTracks - Helper functions (split from helpers.R)
#' Average signal across multiple samples on a common binned grid
#'
#' Creates a uniform grid of bins across the specified region, computes signal
#' scores per bin for each input sample, and then summarizes (e.g. mean) across
#' all samples. This solves the problem of different bigWig/bedGraph files having
#' different internal bin boundaries, which makes direct cross-sample arithmetic
#' impossible without a shared coordinate system.
#'
#' For bigWig inputs, signals are queried efficiently via
#' `rtracklayer::summary()` which computes per-bin statistics directly from the
#' indexed file. For data frame inputs, an overlap-based weighted average is
#' calculated against the uniform bins.
#'
#' @param inputs A character vector of file paths (bigWig, bedGraph, etc.) or a
#'   list of data frames, each with columns `seqnames`, `start`, `end`, `score`.
#' @param region A genomic region string ("chr:start-end") or a GRanges object.
#' @param bin_width Bin size in base pairs (default: 50).
#' @param summary_fun Summary function to apply across samples. One of
#'   `"mean"`, `"median"`, `"max"`, `"min"`, `"sum"` (default: `"mean"`).
#' @param nans_to_zeros Logical. Convert NaN/NA values to zero before
#'   summarizing (default: TRUE). Recommended for bigWig files that may have
#'   missing data in some regions.
#' @return A data frame with columns `seqnames`, `start`, `end`, `score`
#'   representing the summarized signal on the uniform grid.
#' @export
#' @importFrom GenomicRanges GRanges start end tileGenome
#' @importFrom IRanges IRanges
#' @importFrom rtracklayer BigWigFile
#' @importFrom methods is
#' @examples
#' \dontrun{
#' # Average two bigWig files
#' avg_df <- average_signal(
#'   c("sample1.bw", "sample2.bw"),
#'   region = "chr1:1000000-2000000",
#'   bin_width = 100
#' )
#'
#' # Use with ez_coverage for plotting
#' ez_coverage(avg_df, "chr1:1000000-2000000")
#'
#' # Average data frames
#' avg_df <- average_signal(
#'   list(df1, df2, df3),
#'   region = "chr1:1000000-2000000",
#'   summary_fun = "median"
#' )
#' }
average_signal <- function(
  inputs,
  region,
  bin_width = 50,
  summary_fun = c("mean", "median", "max", "min", "sum"),
  nans_to_zeros = TRUE
) {
  summary_fun <- match.arg(summary_fun)
  stopifnot(
    "bin_width must be a positive number" = is.numeric(bin_width) && bin_width > 0,
    "inputs must be a character vector of file paths or a list of data frames" =
      is.character(inputs) || is.list(inputs)
  )

  # Parse region
  region_gr <- parse_region(region)
  region_chr <- as.character(GenomicRanges::seqnames(region_gr))
  region_start <- GenomicRanges::start(region_gr)
  region_end <- GenomicRanges::end(region_gr)


  # Create uniform bins across the region
  bin_starts <- seq(region_start, region_end, by = bin_width)
  bin_ends <- pmin(bin_starts + bin_width - 1, region_end)
  # Drop degenerate last bin if it starts beyond the region

  keep <- bin_starts <= region_end
  bin_starts <- bin_starts[keep]
  bin_ends <- bin_ends[keep]

  bins_gr <- GenomicRanges::GRanges(
    seqnames = region_chr,
    ranges = IRanges::IRanges(start = bin_starts, end = bin_ends)
  )
  n_bins <- length(bins_gr)

  # Collect score matrix: rows = bins, cols = samples
  if (is.character(inputs)) {
    score_matrix <- .query_bigwig_bins(inputs, bins_gr, nans_to_zeros)
  } else if (is.list(inputs)) {
    score_matrix <- .query_df_bins(inputs, bins_gr, region_gr, nans_to_zeros)
  }

  # Apply summary function across samples
  fun <- switch(summary_fun,
    mean   = rowMeans,
    median = function(x, ...) apply(x, 1, stats::median, ...),
    max    = function(x, ...) apply(x, 1, max, ...),
    min    = function(x, ...) apply(x, 1, min, ...),
    sum    = rowSums
  )
  summarized_scores <- fun(score_matrix, na.rm = TRUE)

  # Build output data frame
  data.frame(
    seqnames = region_chr,
    start = bin_starts,
    end = bin_ends,
    score = summarized_scores
  )
}

#' Query bigWig files at uniform bins
#'
#' Internal helper: for each bigWig file, uses `rtracklayer::summary()` to get
#' the mean score in each bin. Returns an n_bins x n_samples matrix.
#'
#' @param files Character vector of bigWig file paths.
#' @param bins_gr GRanges of uniform bins.
#' @param nans_to_zeros Replace NaNs with 0.
#' @return A numeric matrix (n_bins x n_files).
#' @keywords internal
.query_bigwig_bins <- function(files, bins_gr, nans_to_zeros = TRUE) {
  n_bins <- length(bins_gr)
  n_files <- length(files)
  score_matrix <- matrix(NA_real_, nrow = n_bins, ncol = n_files)

  for (i in seq_along(files)) {
    file <- files[i]
    if (!file.exists(file)) {
      stop(sprintf("File not found: %s", file))
    }

    # Use rtracklayer::summary for BigWig files (efficient C-level query)
    ext <- tolower(tools::file_ext(file))
    if (ext %in% c("bw", "bigwig")) {
      bw <- rtracklayer::BigWigFile(file)
      # summary() returns a CompressedGRangesList; extract numeric scores
      bin_summary <- rtracklayer::summary(bw, bins_gr, type = "mean")
      scores <- as.numeric(unlist(bin_summary)$score)
    } else {
      # For non-bigWig formats (bedGraph, BED, etc.), import and overlap
      imported <- import_genomic_data(file, which = bins_gr)
      if (nrow(imported) == 0) {
        scores <- rep(0, n_bins)
      } else {
        imported_gr <- df_to_granges(imported)
        if ("score" %in% colnames(S4Vectors::mcols(imported_gr))) {
          scores <- .overlap_weighted_mean(bins_gr, imported_gr)
        } else {
          warning(sprintf("File %s has no 'score' column. Using zeros.", file))
          scores <- rep(0, n_bins)
        }
      }
    }

    if (nans_to_zeros && any(is.na(scores) | is.nan(scores))) {
      scores[is.na(scores) | is.nan(scores)] <- 0
    }
    score_matrix[, i] <- scores
  }

  score_matrix
}

#' Compute overlap-weighted mean for data frame inputs on uniform bins
#'
#' Internal helper: for each data frame, computes the overlap-weighted mean
#' score onto uniform bins. Returns an n_bins x n_samples matrix.
#'
#' @param dfs List of data frames with seqnames, start, end, score.
#' @param bins_gr GRanges of uniform bins.
#' @param region_gr GRanges of the full region (for filtering).
#' @param nans_to_zeros Replace NaNs with 0.
#' @return A numeric matrix (n_bins x n_dfs).
#' @keywords internal
.query_df_bins <- function(dfs, bins_gr, region_gr, nans_to_zeros = TRUE) {
  n_bins <- length(bins_gr)
  n_dfs <- length(dfs)
  score_matrix <- matrix(NA_real_, nrow = n_bins, ncol = n_dfs)

  region_chr <- as.character(GenomicRanges::seqnames(region_gr))

  for (i in seq_along(dfs)) {
    df <- dfs[[i]]
    if (!is.data.frame(df)) {
      stop("All elements of inputs list must be data frames")
    }
    required_cols <- c("seqnames", "start", "end", "score")
    if (!all(required_cols %in% colnames(df))) {
      stop("Data frames must contain columns: ",
           paste(required_cols, collapse = ", "))
    }

    # Filter to region
    df_filtered <- df[
      df$seqnames == region_chr &
      df$end >= GenomicRanges::start(region_gr) &
      df$start <= GenomicRanges::end(region_gr),
    ]

    if (nrow(df_filtered) == 0) {
      scores <- rep(0, n_bins)
    } else {
      signal_gr <- df_to_granges(df_filtered)
      scores <- .overlap_weighted_mean(bins_gr, signal_gr)
    }

    if (nans_to_zeros && any(is.na(scores) | is.nan(scores))) {
      scores[is.na(scores) | is.nan(scores)] <- 0
    }
    score_matrix[, i] <- scores
  }

  score_matrix
}

#' Overlap-weighted mean of signal onto bins
#'
#' For each bin, finds overlapping signal ranges and computes a weighted mean
#' score where weights are the number of bases each signal range overlaps
#' the bin.
#'
#' @param bins_gr GRanges of uniform bins.
#' @param signal_gr GRanges of signal ranges with a `score` metadata column.
#' @return Numeric vector of length `length(bins_gr)`.
#' @keywords internal
#' @importFrom GenomicRanges findOverlaps
#' @importFrom IRanges pintersect
#' @importFrom S4Vectors queryHits subjectHits
.overlap_weighted_mean <- function(bins_gr, signal_gr) {
  n_bins <- length(bins_gr)
  scores <- rep(0, n_bins)

  hits <- GenomicRanges::findOverlaps(bins_gr, signal_gr)
  if (length(hits) == 0) return(scores)

  q_idx <- S4Vectors::queryHits(hits)
  s_idx <- S4Vectors::subjectHits(hits)

  # Compute overlap widths
  overlap_ranges <- IRanges::pintersect(
    bins_gr[q_idx],
    signal_gr[s_idx]
  )
  overlap_widths <- GenomicRanges::width(overlap_ranges)

  # Weighted mean per bin
  signal_scores <- signal_gr$score[s_idx]
  weighted_scores <- overlap_widths * signal_scores

  # Aggregate by bin index
  sum_weighted <- tapply(weighted_scores, q_idx, sum)
  sum_widths <- tapply(overlap_widths, q_idx, sum)

  bin_indices <- as.integer(names(sum_weighted))
  scores[bin_indices] <- sum_weighted / sum_widths

  scores
}
