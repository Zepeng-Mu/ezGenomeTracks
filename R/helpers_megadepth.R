# ezGenomeTracks - Megadepth helpers for fast BigWig processing

# ---------------------------------------------------------------------------
# Low-level shared helper
# ---------------------------------------------------------------------------

#' Query a BigWig file at pre-defined GRanges bins via megadepth
#'
#' Writes `bins_gr` directly to a temporary BED file (no re-tiling) and calls
#' `megadepth::get_coverage()`. Results are aligned back to `bins_gr` by start
#' position, so the returned vector always has length `length(bins_gr)`.
#' Bins with no BigWig coverage receive a score of 0.
#'
#' @param file Path to BigWig file (local or remote URL).
#' @param bins_gr GRanges of pre-computed bins to query.
#' @param op Summarization operation: "mean" (default), "sum", "min", "max".
#' @param nans_to_zeros Replace NA/NaN results with 0 (default TRUE).
#' @return Numeric vector of length `length(bins_gr)`.
#' @keywords internal
#' @importFrom GenomicRanges seqnames start end
#' @importFrom S4Vectors mcols
.query_megadepth_bins <- function(file, bins_gr, op = "mean", nans_to_zeros = TRUE) {
  op <- match.arg(op, c("sum", "mean", "min", "max"))

  temp_bed   <- tempfile(fileext = ".bed")
  out_prefix <- tempfile("megadepth_cov_")
  on.exit({
    unlink(temp_bed)
    unlink(Sys.glob(paste0(out_prefix, "*")))
  }, add = TRUE)

  # Write bins directly to BED (0-based start coordinates)
  bed_data <- data.frame(
    chr   = as.character(GenomicRanges::seqnames(bins_gr)),
    start = GenomicRanges::start(bins_gr) - 1L,
    end   = GenomicRanges::end(bins_gr)
  )
  utils::write.table(
    bed_data,
    file      = temp_bed,
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

  result_gr <- megadepth::get_coverage(
    bigwig_file = file,
    op          = op,
    annotation  = temp_bed,
    prefix      = out_prefix
  )

  # Align result back to bins_gr by 1-based start position.
  # megadepth may return fewer rows (zero-coverage bins are absent) or in
  # BigWig chromosome order rather than BED order.
  n <- length(bins_gr)
  scores <- rep(0, n)

  if (length(result_gr) > 0 && "score" %in% names(S4Vectors::mcols(result_gr))) {
    idx <- match(GenomicRanges::start(result_gr), GenomicRanges::start(bins_gr))
    valid <- !is.na(idx)
    scores[idx[valid]] <- result_gr$score[valid]
  }

  if (nans_to_zeros) scores[is.na(scores) | is.nan(scores)] <- 0

  scores
}

# ---------------------------------------------------------------------------
# Public helper used by import_genomic_data()
# ---------------------------------------------------------------------------

#' Import BigWig coverage using megadepth
#'
#' Fast region-based BigWig queries using the megadepth command-line tool.
#' The region is tiled into `n_bins` equal-width bins; each bin receives the
#' summary score computed by megadepth.
#'
#' @param file Path to BigWig file (local or remote URL).
#' @param which GRanges specifying the genomic region to query.
#' @param op Summary operation: "mean" (default), "sum", "min", "max".
#' @param n_bins Number of bins used to tile `which`. Default: 2000.
#' @return GRanges with a `score` metadata column.
#' @noRd
#' @importFrom GenomicRanges reduce GRanges width start end seqnames
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors mcols
#' @keywords internal
import_bigwig_megadepth <- function(file, which, op = "mean", n_bins = 2000L) {
  op <- match.arg(op, c("sum", "mean", "min", "max"))
  stopifnot(
    "which must be a non-empty GRanges object" =
      methods::is(which, "GRanges") && length(which) > 0,
    "n_bins must be a positive integer" =
      is.numeric(n_bins) && length(n_bins) == 1 && n_bins >= 1
  )

  # Tile the (reduced) region into n_bins equal-width bins
  which_reduced <- GenomicRanges::reduce(which)
  total_width   <- sum(GenomicRanges::width(which_reduced))
  bin_width     <- max(1L, as.integer(ceiling(total_width / as.integer(n_bins))))

  bins_gr <- do.call(c, lapply(seq_along(which_reduced), function(i) {
    chr_i   <- as.character(GenomicRanges::seqnames(which_reduced)[i])
    start_i <- GenomicRanges::start(which_reduced)[i]
    end_i   <- GenomicRanges::end(which_reduced)[i]
    starts  <- seq.int(start_i, end_i, by = bin_width)
    ends    <- pmin(starts + bin_width - 1L, end_i)
    GenomicRanges::GRanges(
      seqnames = chr_i,
      ranges   = IRanges::IRanges(start = starts, end = ends)
    )
  }))

  scores <- .query_megadepth_bins(file, bins_gr, op = op)

  S4Vectors::mcols(bins_gr)$score <- scores
  bins_gr
}

