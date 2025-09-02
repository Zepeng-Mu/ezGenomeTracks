#' Example signal data
#'
#' A small example dataset containing signal values for a genomic region.
#' This dataset is in a format similar to what would be extracted from a bigWig file.
#'
#' @format A data frame with 1000 rows and 3 variables:
#' \describe{
#'   \item{seqnames}{Chromosome name}
#'   \item{start}{Start position}
#'   \item{end}{End position}
#'   \item{score}{Signal value}
#' }
#' @source Simulated data
"example_signal" <- NULL

#' Example peak data
#'
#' A small example dataset containing genomic peaks, similar to what would be
#' found in a BED file.
#'
#' @format A data frame with 50 rows and 4 variables:
#' \describe{
#'   \item{seqnames}{Chromosome name}
#'   \item{start}{Start position}
#'   \item{end}{End position}
#'   \item{score}{Peak score}
#' }
#' @source Simulated data
"example_peaks" <- NULL

#' Example gene data
#'
#' A small example dataset containing gene annotations, similar to what would be
#' extracted from a GTF/GFF file.
#'
#' @format A data frame with 20 rows and 7 variables:
#' \describe{
#'   \item{seqnames}{Chromosome name}
#'   \item{start}{Start position}
#'   \item{end}{End position}
#'   \item{strand}{Strand (+ or -)}
#'   \item{type}{Feature type (gene, exon, CDS, etc.)}
#'   \item{gene_id}{Gene identifier}
#'   \item{gene_name}{Gene name}
#' }
#' @source Simulated data
"example_genes" <- NULL

#' Example interaction data
#'
#' A small example dataset containing genomic interactions, similar to what would be
#' found in a BEDPE file.
#'
#' @format A data frame with 30 rows and 7 variables:
#' \describe{
#'   \item{seqnames1}{Chromosome name for the first anchor}
#'   \item{start1}{Start position for the first anchor}
#'   \item{end1}{End position for the first anchor}
#'   \item{seqnames2}{Chromosome name for the second anchor}
#'   \item{start2}{Start position for the second anchor}
#'   \item{end2}{End position for the second anchor}
#'   \item{score}{Interaction score}
#' }
#' @source Simulated data
"example_interactions" <- NULL

#' Example Hi-C data
#'
#' A small example dataset containing Hi-C contact matrix data.
#'
#' @format A data frame with 400 rows and 3 variables:
#' \describe{
#'   \item{bin1}{Bin index for the first position}
#'   \item{bin2}{Bin index for the second position}
#'   \item{count}{Contact count}
#' }
#' @source Simulated data
"example_hic" <- NULL