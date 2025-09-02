# Script to create example datasets for the ezGenomeTracks package

library(GenomicRanges)
library(IRanges)
library(usethis)

# Set seed for reproducibility
set.seed(42)

# Create example signal data (similar to bigWig)
chrom <- "chr1"
start <- seq(1000000, 1099900, by = 100)
end <- start + 100
score <- c(
  # Create some peaks and valleys for visualization
  runif(250, 0, 3),
  runif(100, 2, 8),  # A medium peak
  runif(150, 0, 2),
  runif(100, 5, 15), # A high peak
  runif(150, 0, 3),
  runif(100, 3, 10), # Another medium peak
  runif(150, 0, 4)
)

example_signal <- data.frame(
  seqnames = chrom,
  start = start,
  end = end,
  score = score
)

# Create example peak data (similar to BED)
peak_count <- 50
peak_width <- sample(500:2000, peak_count, replace = TRUE)
peak_start <- sample(seq(1000000, 1100000 - max(peak_width)), peak_count)
peak_end <- peak_start + peak_width
peak_score <- runif(peak_count, 1, 10)

example_peaks <- data.frame(
  seqnames = chrom,
  start = peak_start,
  end = peak_end,
  score = peak_score
)

# Create example gene data (similar to GTF/GFF)
gene_count <- 5
gene_width <- sample(5000:15000, gene_count, replace = TRUE)
gene_start <- seq(1010000, 1090000, length.out = gene_count)
gene_end <- gene_start + gene_width
gene_strand <- sample(c("+", "-"), gene_count, replace = TRUE)
gene_id <- paste0("GENE", 1:gene_count)
gene_name <- paste0("Gene", 1:gene_count)

# Create a data frame to hold all gene features
example_genes <- data.frame()

# For each gene, create gene, exons, and CDS features
for (i in 1:gene_count) {
  # Gene feature
  gene_row <- data.frame(
    seqnames = chrom,
    start = gene_start[i],
    end = gene_end[i],
    strand = gene_strand[i],
    type = "gene",
    gene_id = gene_id[i],
    gene_name = gene_name[i]
  )
  
  # Create 2-5 exons per gene
  exon_count <- sample(2:5, 1)
  exon_width <- sample(100:500, exon_count, replace = TRUE)
  
  # Distribute exons along the gene
  gene_length <- gene_end[i] - gene_start[i]
  exon_positions <- sort(sample(seq(0, gene_length - max(exon_width)), exon_count))
  exon_start <- gene_start[i] + exon_positions
  exon_end <- exon_start + exon_width
  
  # Create exon features
  exon_rows <- data.frame(
    seqnames = chrom,
    start = exon_start,
    end = exon_end,
    strand = gene_strand[i],
    type = "exon",
    gene_id = gene_id[i],
    gene_name = gene_name[i]
  )
  
  # Combine gene and exon features
  example_genes <- rbind(example_genes, gene_row, exon_rows)
}

# Create example interaction data (similar to BEDPE)
interaction_count <- 30

# Create anchors with some clustering
anchor1_start <- sample(seq(1000000, 1050000), interaction_count, replace = TRUE)
anchor1_end <- anchor1_start + sample(100:500, interaction_count, replace = TRUE)

anchor2_start <- sample(seq(1050000, 1100000), interaction_count, replace = TRUE)
anchor2_end <- anchor2_start + sample(100:500, interaction_count, replace = TRUE)

interaction_score <- runif(interaction_count, 1, 10)

example_interactions <- data.frame(
  seqnames1 = chrom,
  start1 = anchor1_start,
  end1 = anchor1_end,
  seqnames2 = chrom,
  start2 = anchor2_start,
  end2 = anchor2_end,
  score = interaction_score
)

# Create example Hi-C data
# Define bins for a small region
bin_size <- 5000
bins <- seq(1000000, 1100000, by = bin_size)
bin_count <- length(bins) - 1

# Create a sparse matrix representation
bin1 <- c()
bin2 <- c()
count <- c()

# Add some structured signal (diagonal-heavy with some off-diagonal interactions)
for (i in 1:bin_count) {
  for (j in i:bin_count) {
    # Higher counts for nearby bins (diagonal), lower for distant bins
    if (i == j) {
      # Diagonal elements have highest counts
      bin1 <- c(bin1, i)
      bin2 <- c(bin2, j)
      count <- c(count, sample(50:100, 1))
    } else if (abs(i - j) <= 3) {
      # Near-diagonal elements have medium counts
      bin1 <- c(bin1, i)
      bin2 <- c(bin2, j)
      count <- c(count, sample(10:50, 1))
    } else if (runif(1) < 0.3) {
      # Some random off-diagonal interactions
      bin1 <- c(bin1, i)
      bin2 <- c(bin2, j)
      count <- c(count, sample(1:20, 1))
    }
  }
}

example_hic <- data.frame(
  bin1 = bin1,
  bin2 = bin2,
  count = count
)

# Save the example datasets using usethis
usethis::use_data(example_signal, overwrite = TRUE)
usethis::use_data(example_peaks, overwrite = TRUE)
usethis::use_data(example_genes, overwrite = TRUE)
usethis::use_data(example_interactions, overwrite = TRUE)
usethis::use_data(example_hic, overwrite = TRUE)

# Also save some example files in inst/extdata for demonstration

# Save example bigWig-like file
write.table(
  example_signal,
  file = "inst/extdata/example_signal.bedgraph",
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

# Save example BED-like file
write.table(
  example_peaks[, c("seqnames", "start", "end", "score")],
  file = "inst/extdata/example_peaks.bed",
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

# Save example GTF-like file
gtf_data <- example_genes
gtf_data$source <- "ezGenomeTracks"
gtf_data$score <- "."
gtf_data$phase <- "."
gtf_data$attributes <- paste0(
  "gene_id \"", gtf_data$gene_id, "\"; ",
  "gene_name \"", gtf_data$gene_name, "\";"
)

write.table(
  gtf_data[, c("seqnames", "source", "type", "start", "end", "score", "strand", "phase", "attributes")],
  file = "inst/extdata/example_genes.gtf",
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

# Save example BEDPE-like file
write.table(
  example_interactions,
  file = "inst/extdata/example_interactions.bedpe",
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

# Save example Hi-C matrix
write.table(
  example_hic,
  file = "inst/extdata/example_hic.matrix",
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

# Create a README file for the example data
readme_content <- "# Example Data for ezGenomeTracks

This directory contains example data files for demonstrating the ezGenomeTracks package:

- example_signal.bedgraph: A bedGraph-like file containing signal values
- example_peaks.bed: A BED-like file containing genomic peaks
- example_genes.gtf: A GTF-like file containing gene annotations
- example_interactions.bedpe: A BEDPE-like file containing genomic interactions
- example_hic.matrix: A matrix file containing Hi-C contact data

These files are used in the package vignettes and examples.
"

writeLines(readme_content, "inst/extdata/README.md")

cat("Example datasets created and saved successfully.\n")