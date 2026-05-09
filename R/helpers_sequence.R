# ezGenomeTracks - Sequence helper functions

#' Fetch and process nucleotide sequence data for a genomic region
#'
#' Internal helper that retrieves the DNA sequence for a given region from a
#' BSgenome object and returns a per-nucleotide data frame ready for plotting.
#'
#' @param genome A BSgenome object (e.g., `BSgenome.Hsapiens.UCSC.hg38::Hsapiens`).
#' @param region_gr A GRanges object specifying the region to fetch.
#' @param palette A named character vector mapping nucleotides to colors (A, C, G, T, N).
#'
#' @return A data frame with columns:
#'   - `position`: Genomic coordinate of each nucleotide.
#'   - `nucleotide`: Single-character nucleotide (A, C, G, T, or N).
#'   - `fill`: Hex color string for the nucleotide.
#' @noRd
process_sequence_data <- function(genome, region_gr, palette) {
  if (!requireNamespace("BSgenome", quietly = TRUE)) {
    stop(
      "The 'BSgenome' package is required to use ez_sequence(). ",
      "Install it with: BiocManager::install('BSgenome')"
    )
  }

  if (!methods::is(genome, "BSgenome")) {
    stop(
      "'genome' must be a BSgenome object (e.g., BSgenome.Hsapiens.UCSC.hg38::Hsapiens). ",
      "Install the appropriate BSgenome data package and pass the genome object directly."
    )
  }

  # Fetch sequence for the region
  seq <- BSgenome::getSeq(genome, region_gr)

  # Convert DNAStringSet → character → per-nucleotide vector
  seq_char <- as.character(seq[[1]])
  nucs <- strsplit(seq_char, "")[[1]]

  # Genomic positions (1-based, matching GRanges start)
  region_start <- GenomicRanges::start(region_gr)
  positions <- seq(region_start, by = 1L, length.out = length(nucs))

  # Map nucleotides to UCSC colors; unknown characters get the N color
  fill_colors <- palette[nucs]
  fill_colors[is.na(fill_colors)] <- palette["N"]

  data.frame(
    position  = positions,
    nucleotide = nucs,
    fill      = unname(fill_colors),
    stringsAsFactors = FALSE
  )
}
