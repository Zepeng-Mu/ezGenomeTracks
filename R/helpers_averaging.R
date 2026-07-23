# ezGenomeTracks - Helper functions for averaging coverage across samples

#' Average coverage across multiple samples on a common binned grid
#'
#' Tiles the specified region into `n_bins` equal-width bins, queries each input
#' at those bins, then summarizes (e.g. mean) across all samples. Because every
#' sample is queried on the **same** bin grid, no re-alignment is needed.
#'
#' For BigWig inputs, megadepth is used directly via [.query_megadepth_bins()].
#' For data frame inputs, an overlap-based weighted average is calculated.
#'
#' @param inputs A character vector of BigWig file paths, or a list of data
#'   frames each with columns `seqnames`, `start`, `end`, `score`.
#' @param region A genomic region string ("chr:start-end") or a GRanges object.
#' @param n_bins Number of equal-width bins to tile the region with (default: 2000).
#' @param summary_fun Summary function applied across samples per bin. One of
#'   `"mean"`, `"median"`, `"max"`, `"min"`, `"sum"` (default: `"mean"`).
#' @param nans_to_zeros Logical. Replace NA/NaN per-bin scores with 0 before
#'   summarizing (default: TRUE).
#' @param verbose Logical. If `TRUE`, print megadepth progress output when
#'   querying BigWig files. If `FALSE` (default), suppress routine messages.
#' @return A data frame with columns `seqnames`, `start`, `end`, `score`.
#' @export
#' @importFrom GenomicRanges GRanges start end seqnames
#' @importFrom IRanges IRanges
#' @importFrom methods is
#' @examples
#' \dontrun{
#' # Average two BigWig files
#' avg_df <- average_coverage(
#'   c("sample1.bw", "sample2.bw"),
#'   region = "chr1:1000000-2000000",
#'   n_bins = 1000
#' )
#'
#' # Use with ez_coverage for plotting
#' ez_coverage(avg_df, "chr1:1000000-2000000")
#'
#' # Average data frames
#' avg_df <- average_coverage(
#'   list(df1, df2, df3),
#'   region = "chr1:1000000-2000000",
#'   summary_fun = "median"
#' )
#' }
average_coverage <- function(
  inputs,
  region,
  n_bins = 2000L,
  summary_fun = c("mean", "median", "max", "min", "sum"),
  nans_to_zeros = TRUE,
  verbose = FALSE
) {
  summary_fun <- match.arg(summary_fun)
  stopifnot(
    "n_bins must be a positive integer" =
      is.numeric(n_bins) && length(n_bins) == 1 && n_bins >= 1,
    "inputs must be a character vector of file paths or a list of data frames" =
      is.character(inputs) || is.list(inputs)
  )

  # Parse region
  region_gr    <- parse_region(region)
  region_chr   <- as.character(GenomicRanges::seqnames(region_gr))
  region_start <- GenomicRanges::start(region_gr)
  region_end   <- GenomicRanges::end(region_gr)

  # Build the shared bin grid from n_bins
  region_width <- region_end - region_start + 1L
  bin_width    <- max(1L, as.integer(ceiling(region_width / as.integer(n_bins))))
  bin_starts   <- seq.int(region_start, region_end, by = bin_width)
  bin_ends     <- pmin(bin_starts + bin_width - 1L, region_end)

  bins_gr <- GenomicRanges::GRanges(
    seqnames = region_chr,
    ranges   = IRanges::IRanges(start = bin_starts, end = bin_ends)
  )

  # Collect score matrix: rows = bins, cols = samples
  if (is.character(inputs)) {
    score_matrix <- .query_bigwig_bins(
      inputs,
      bins_gr,
      nans_to_zeros = nans_to_zeros,
      verbose = verbose
    )
  } else {
    score_matrix <- .query_df_bins(inputs, bins_gr, region_gr, nans_to_zeros)
  }

  # Apply summary function across samples (columns)
  fun <- switch(summary_fun,
    mean   = rowMeans,
    median = function(x, ...) apply(x, 1, stats::median, ...),
    max    = function(x, ...) apply(x, 1, max, ...),
    min    = function(x, ...) apply(x, 1, min, ...),
    sum    = rowSums
  )

  data.frame(
    seqnames = region_chr,
    start    = bin_starts,
    end      = bin_ends,
    score    = fun(score_matrix, na.rm = TRUE)
  )
}

#' Query BigWig files at shared uniform bins via megadepth
#'
#' Internal helper used by [average_coverage()]. Queries each BigWig file at
#' exactly the same pre-computed `bins_gr`, so no bin-alignment is needed when
#' stacking scores across samples.
#'
#' @param files Character vector of file paths.
#' @param bins_gr Pre-computed GRanges of uniform bins (shared across all files).
#' @param nans_to_zeros Replace NA/NaN with 0.
#' @param verbose Logical. If `TRUE`, print megadepth progress output.
#' @return A numeric matrix (n_bins × n_files).
#' @keywords internal
.query_bigwig_bins <- function(files, bins_gr, nans_to_zeros = TRUE, verbose = FALSE) {
  n_bins  <- length(bins_gr)
  n_files <- length(files)
  score_matrix <- matrix(0, nrow = n_bins, ncol = n_files)

  for (i in seq_along(files)) {
    file <- files[i]
    ext  <- tolower(tools::file_ext(file))
    is_bigwig <- ext %in% c("bw", "bigwig") || grepl("^https?://", file)

    if (is_bigwig) {
      # Use megadepth with the shared bins — fast and already unified
      score_matrix[, i] <- .query_megadepth_bins(
        file,
        bins_gr,
        op = "mean",
        nans_to_zeros = nans_to_zeros,
        verbose = verbose
      )
    } else {
      # Non-BigWig (bedGraph, BED, etc.): import full region then overlap-weight
      if (!grepl("^https?://", file) && !file.exists(file)) {
        stop(sprintf("File not found: %s", file))
      }
      region_hull <- GenomicRanges::GRanges(
        seqnames = as.character(GenomicRanges::seqnames(bins_gr)[1]),
        ranges   = IRanges::IRanges(
          start = min(GenomicRanges::start(bins_gr)),
          end   = max(GenomicRanges::end(bins_gr))
        )
      )
      imported    <- import_genomic_data(file, which = region_hull)
      if (nrow(imported) == 0) {
        scores <- rep(0, n_bins)
      } else {
        imported_gr <- df_to_granges(imported)
        if ("score" %in% names(S4Vectors::mcols(imported_gr))) {
          scores <- .overlap_weighted_mean(bins_gr, imported_gr)
        } else {
          warning(sprintf("File %s has no 'score' column. Using zeros.", file))
          scores <- rep(0, n_bins)
        }
      }
      if (nans_to_zeros) scores[is.na(scores) | is.nan(scores)] <- 0
      score_matrix[, i] <- scores
    }
  }

  score_matrix
}

#' Compute overlap-weighted mean for data frame inputs on uniform bins
#'
#' @param dfs List of data frames with seqnames, start, end, score.
#' @param bins_gr GRanges of uniform bins.
#' @param region_gr GRanges of the full region (for filtering).
#' @param nans_to_zeros Replace NaNs with 0.
#' @return A numeric matrix (n_bins x n_dfs).
#' @keywords internal
.query_df_bins <- function(dfs, bins_gr, region_gr, nans_to_zeros = TRUE) {
  n_bins <- length(bins_gr)
  n_dfs  <- length(dfs)
  score_matrix <- matrix(NA_real_, nrow = n_bins, ncol = n_dfs)

  region_chr <- as.character(GenomicRanges::seqnames(region_gr))

  for (i in seq_along(dfs)) {
    df <- dfs[[i]]
    if (!is.data.frame(df)) stop("All elements of inputs list must be data frames")
    required_cols <- c("seqnames", "start", "end", "score")
    if (!all(required_cols %in% colnames(df))) {
      stop("Data frames must contain columns: ", paste(required_cols, collapse = ", "))
    }

    df_filtered <- df[
      df$seqnames == region_chr &
      df$end   >= GenomicRanges::start(region_gr) &
      df$start <= GenomicRanges::end(region_gr),
    ]

    if (nrow(df_filtered) == 0) {
      scores <- rep(0, n_bins)
    } else {
      signal_gr <- df_to_granges(df_filtered)
      scores    <- .overlap_weighted_mean(bins_gr, signal_gr)
    }

    if (nans_to_zeros) scores[is.na(scores) | is.nan(scores)] <- 0
    score_matrix[, i] <- scores
  }

  score_matrix
}

#' Overlap-weighted mean of signal onto bins
#'
#' @param bins_gr GRanges of uniform bins.
#' @param signal_gr GRanges of signal ranges with a `score` metadata column.
#' @return Numeric vector of length `length(bins_gr)`.
#' @keywords internal
#' @importFrom GenomicRanges findOverlaps width
#' @importFrom IRanges pintersect
#' @importFrom S4Vectors queryHits subjectHits
.overlap_weighted_mean <- function(bins_gr, signal_gr) {
  n_bins <- length(bins_gr)
  scores <- rep(0, n_bins)

  hits <- GenomicRanges::findOverlaps(bins_gr, signal_gr)
  if (length(hits) == 0) return(scores)

  q_idx <- S4Vectors::queryHits(hits)
  s_idx <- S4Vectors::subjectHits(hits)

  overlap_widths  <- GenomicRanges::width(IRanges::pintersect(bins_gr[q_idx], signal_gr[s_idx]))
  signal_scores   <- signal_gr$score[s_idx]
  weighted_scores <- overlap_widths * signal_scores

  sum_weighted <- tapply(weighted_scores, q_idx, sum)
  sum_widths   <- tapply(overlap_widths,  q_idx, sum)

  bin_indices <- as.integer(names(sum_weighted))
  scores[bin_indices] <- sum_weighted / sum_widths

  scores
}
