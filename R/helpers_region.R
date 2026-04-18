# ezGenomeTracks - Helper functions (split from helpers.R)
#' Parse a genomic region string into a GRanges object
#'
#' This function parses a genomic region string in the format "chr:start-end" into a GRanges object.
#'
#' @param region A string specifying a genomic region (e.g., "chr1:1000000-2000000")
#' @return A GRanges object representing the specified region
#' @export
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @examples
#' \dontrun{
#' region_gr <- parse_region("chr1:1000000-2000000")
#' }
parse_region <- function(region) {
  if (!is.character(region) || length(region) != 1) {
    stop("Region must be a single character string")
  }

  # Remove commas (thousand separators) from the region string
  region <- gsub(",", "", region)

  # Parse the region string using regular expression
  # Match pattern: chromosome followed by any separator (: _ -) then start and end positions
  matches <- regexpr("^(.+?)[:_-](\\d+)[:_-](\\d+)$", region, perl = TRUE)
  if (matches == -1) {
    stop(
      "Region must be in the format 'chr*start*end' where * can be `:`, `_`, or `-`"
    )
  }

  # Extract the matched groups
  parsed <- regmatches(region, matches)
  groups <- stringr::str_match(parsed, "^(.+?)[:_-](\\d+)[:_-](\\d+)$")

  chr <- groups[, 2]
  start <- as.numeric(groups[, 3])
  end <- as.numeric(groups[, 4])

  if (is.na(start) || is.na(end)) {
    stop("Start and end positions must be numeric")
  }

  # Create GRanges object
  gr <- GenomicRanges::GRanges(
    seqnames = chr,
    ranges = IRanges::IRanges(start = start, end = end)
  )

  return(gr)
}
