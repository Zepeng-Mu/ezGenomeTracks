#' @title Test Data Generation Utilities for ezGenomeTracks
#' @description Helper functions to generate synthetic genomic data for testing purposes
#' @name test-data-generators
#' @keywords internal testing

# Set seed for reproducibility in tests
set.seed(123)

#' Generate synthetic genomic signal data
#' @description Creates test signal track data similar to bigWig/BedGraph format
#' @param chr Chromosome name (default: "chr1")
#' @param start Start position (default: 1000000)
#' @param end End position (default: 1100000)
#' @param bin_size Size of each bin (default: 100)
#' @param n_peaks Number of peak regions to add (default: 3)
#' @param peak_height Range for peak heights (default: c(5, 15))
#' @param background_range Range for background signal (default: c(0, 3))
#' @return data.frame with columns: seqnames, start, end, score
#' @export
#' @examples
#' signal_data <- create_test_genomic_data()
#' head(signal_data)
create_test_genomic_data <- function(chr = "chr1",
                                    start = 1000000,
                                    end = 1100000,
                                    bin_size = 100,
                                    n_peaks = 3,
                                    peak_height = c(5, 15),
                                    background_range = c(0, 3)) {

  # Generate basic bins
  bins <- seq(start, end, by = bin_size)
  positions <- bins[1:floor((end - start) / bin_size)]

  # Create background signal
  score <- runif(length(positions), background_range[1], background_range[2])

  # Add peaks at random positions
  for (i in 1:n_peaks) {
    peak_start_idx <- sample(1:(length(positions) - 5), 1)
    peak_length <- sample(3:8, 1)
    peak_height_val <- runif(1, peak_height[1], peak_height[2])

    end_idx <- min(peak_start_idx + peak_length - 1, length(positions))
    score[peak_start_idx:end_idx] <- peak_height_val +
      rnorm(end_idx - peak_start_idx + 1, 0, peak_height_val * 0.2)
  }

  # Create data frame
  result <- data.frame(
    seqnames = chr,
    start = positions,
    end = positions + bin_size,
    score = pmax(0, score)  # Ensure non-negative scores
  )

  return(result)
}

#' Generate synthetic gene annotation data
#' @description Creates test gene annotation data similar to GTF/GFF format
#' @param chr Chromosome name (default: "chr1")
#' @param gene_count Number of genes to create (default: 5)
#' @param gene_length_range Range for gene lengths (default: c(5000, 20000))
#' @param exon_per_gene_range Range for number of exons per gene (default: c(2, 6))
#' @return data.frame with columns: seqnames, start, end, strand, type, gene_id, gene_name
#' @export
#' @examples
#' gene_data <- create_test_gene_data()
#' head(gene_data)
create_test_gene_data <- function(chr = "chr1",
                                 gene_count = 5,
                                 gene_length_range = c(5000, 20000),
                                 exon_per_gene_range = c(2, 6)) {

  # Calculate positions for genes
  start_pos <- seq(1000000, 1100000, length.out = gene_count + 1)[1:gene_count]
  gene_lengths <- sample(seq(gene_length_range[1], gene_length_range[2]), gene_count, replace = TRUE)
  end_pos <- start_pos + gene_lengths

  # Gene attributes
  gene_ids <- paste0("GENE", seq_len(gene_count))
  gene_names <- paste0("TestGene", seq_len(gene_count))
  strands <- sample(c("+", "-"), gene_count, replace = TRUE)

  result <- data.frame(
    seqnames = character(),
    start = integer(),
    end = integer(),
    strand = character(),
    type = character(),
    gene_id = character(),
    gene_name = character(),
    stringsAsFactors = FALSE
  )

  # Create genes and exons
  for (i in seq_len(gene_count)) {
    # Gene feature
    gene_row <- data.frame(
      seqnames = chr,
      start = start_pos[i],
      end = end_pos[i],
      strand = strands[i],
      type = "gene",
      gene_id = gene_ids[i],
      gene_name = gene_names[i],
      stringsAsFactors = FALSE
    )

    result <- rbind(result, gene_row)

    # Exon features
    n_exons <- sample(exon_per_gene_range[1]:exon_per_gene_range[2], 1)
    exon_lengths <- sample(seq(100, 800), n_exons, replace = TRUE)

    # Distribute exons along gene
    gene_length <- end_pos[i] - start_pos[i]
    exon_positions <- sort(sample(seq(0, gene_length - max(exon_lengths)), n_exons))

    # Add exons for this gene
    for (j in seq_len(n_exons)) {
      exon_start <- start_pos[i] + exon_positions[j]
      exon_end <- exon_start + exon_lengths[j]

      exon_row <- data.frame(
        seqnames = chr,
        start = exon_start,
        end = exon_end,
        strand = strands[i],
        type = "exon",
        gene_id = gene_ids[i],
        gene_name = gene_names[i],
        stringsAsFactors = FALSE
      )

      result <- rbind(result, exon_row)
    }
  }

  return(result)
}

#' Generate synthetic signal track data
#' @description Creates test signal track data with multiple samples/conditions
#' @param chr Chromosome name (default: "chr1")
#' @param start Start position (default: 1000000)
#' @param end End position (default: 1100000)
#' @param n_samples Number of samples (default: 2)
#' @param n_conditions Number of conditions (default: 2)
#' @param sample_names Names for samples (default: NULL, auto-generated)
#' @param condition_names Names for conditions (default: NULL, auto-generated)
#' @return data.frame with columns: seqnames, start, end, score, sample, condition
#' @export
#' @examples
#' signal_data <- create_test_signal_data()
#' head(signal_data)
create_test_signal_data <- function(chr = "chr1",
                                   start = 1000000,
                                   end = 1100000,
                                   n_samples = 2,
                                   n_conditions = 2,
                                   sample_names = NULL,
                                   condition_names = NULL) {

  # Set default names if not provided
  if (is.null(sample_names)) {
    sample_names <- paste0("Sample", seq_len(n_samples))
  }
  if (is.null(condition_names)) {
    condition_names <- paste0("Condition", seq_len(n_conditions))
  }

  # Generate base signal for each sample-condition combination
  result <- data.frame(
    seqnames = character(),
    start = integer(),
    end = integer(),
    score = numeric(),
    sample = character(),
    condition = character(),
    stringsAsFactors = FALSE
  )

  bin_size <- 100
  positions <- seq(start, end, by = bin_size)[1:floor((end - start) / bin_size)]

  for (sample_idx in seq_len(n_samples)) {
    for (condition_idx in seq_len(n_conditions)) {
      # Generate signal with some variation per sample/condition
      base_signal <- runif(length(positions), 0, 3)

      # Add some peaks
      n_peaks <- sample(2:5, 1)
      for (peak in seq_len(n_peaks)) {
        peak_pos <- sample(length(positions), 1)
        peak_height <- runif(1, 5, 15)
        peak_width <- sample(2:5, 1)

        end_idx <- min(peak_pos + peak_width - 1, length(positions))
        base_signal[peak_pos:end_idx] <- peak_height +
          rnorm(end_idx - peak_pos + 1, 0, peak_height * 0.2)
      }

      # Create data for this sample-condition
      sample_data <- data.frame(
        seqnames = chr,
        start = positions,
        end = positions + bin_size,
        score = pmax(0, base_signal),
        sample = sample_names[sample_idx],
        condition = condition_names[condition_idx],
        stringsAsFactors = FALSE
      )

      result <- rbind(result, sample_data)
    }
  }

  return(result)
}

#' Generate synthetic interaction data
#' @description Creates test genomic interaction data similar to BEDPE format
#' @param chr Chromosome name (default: "chr1")
#' @param start Start position (default: 1000000)
#' @param end End position (default: 1100000)
#' @param n_interactions Number of interactions to create (default: 20)
#' @param cluster_interactions Whether to cluster interactions (default: TRUE)
#' @param interaction_range Range for interaction scores (default: c(1, 10))
#' @return data.frame with columns: seqnames1, start1, end1, seqnames2, start2, end2, score
#' @export
#' @examples
#' interaction_data <- create_test_interaction_data()
#' head(interaction_data)
create_test_interaction_data <- function(chr = "chr1",
                                        start = 1000000,
                                        end = 1100000,
                                        n_interactions = 20,
                                        cluster_interactions = TRUE,
                                        interaction_range = c(1, 10)) {

  result <- data.frame(
    seqnames1 = character(),
    start1 = integer(),
    end1 = integer(),
    seqnames2 = character(),
    start2 = integer(),
    end2 = integer(),
    score = numeric(),
    stringsAsFactors = FALSE
  )

  if (cluster_interactions) {
    # Create clustered interactions
    n_clusters <- sample(2:4, 1)
    cluster_centers <- seq(start, end, length.out = n_clusters + 1)[1:n_clusters]

    interactions_per_cluster <- round(n_interactions / n_clusters)

    for (cluster_idx in seq_len(n_clusters)) {
      cluster_center <- cluster_centers[cluster_idx]

      for (i in seq_len(interactions_per_cluster)) {
        # Create anchors around cluster center
        anchor1_center <- cluster_center + rnorm(1, 0, 10000)
        anchor2_center <- cluster_center + rnorm(1, 20000, 15000)

        # Ensure anchors are within bounds
        anchor1_center <- max(start + 1000, min(anchor1_center, end - 5000))
        anchor2_center <- max(start + 5000, min(anchor2_center, end - 1000))

        # Create anchor regions
        anchor1_width <- sample(100:800, 1)
        anchor2_width <- sample(100:800, 1)

        anchor1_start <- anchor1_center - anchor1_width / 2
        anchor1_end <- anchor1_center + anchor1_width / 2
        anchor2_start <- anchor2_center - anchor2_width / 2
        anchor2_end <- anchor2_center + anchor2_width / 2

        # Ensure proper ordering
        if (anchor1_start > anchor2_start) {
          temp <- anchor1_start
          anchor1_start <- anchor2_start
          anchor2_start <- temp
          temp <- anchor1_end
          anchor1_end <- anchor2_end
          anchor2_end <- temp
        }

        interaction_score <- runif(1, interaction_range[1], interaction_range[2])

        interaction_row <- data.frame(
          seqnames1 = chr,
          start1 = as.integer(anchor1_start),
          end1 = as.integer(anchor1_end),
          seqnames2 = chr,
          start2 = as.integer(anchor2_start),
          end2 = as.integer(anchor2_end),
          score = interaction_score,
          stringsAsFactors = FALSE
        )

        result <- rbind(result, interaction_row)
      }
    }
  } else {
    # Random interactions
    for (i in seq_len(n_interactions)) {
      # Create random anchors
      anchor1_start <- sample(seq(start, end - 5000), 1)
      anchor1_end <- anchor1_start + sample(100:800, 1)
      anchor2_start <- sample(seq(anchor1_end + 1000, end - 1000), 1)
      anchor2_end <- anchor2_start + sample(100:800, 1)

      interaction_score <- runif(1, interaction_range[1], interaction_range[2])

      interaction_row <- data.frame(
        seqnames1 = chr,
        start1 = anchor1_start,
        end1 = anchor1_end,
        seqnames2 = chr,
        start2 = anchor2_start,
        end2 = anchor2_end,
        score = interaction_score,
        stringsAsFactors = FALSE
      )

      result <- rbind(result, interaction_row)
    }
  }

  return(result)
}

#' Generate synthetic Hi-C data
#' @description Creates test Hi-C contact matrix data
#' @param chr Chromosome name (default: "chr1")
#' @param start Start position (default: 1000000)
#' @param end End position (default: 1100000)
#' @param bin_size Size of Hi-C bins (default: 5000)
#' @param enrichment_factor Factor for diagonal enrichment (default: 10)
#' @param distance_decay Exponential decay factor for distance (default: 0.001)
#' @return data.frame with columns: bin1, bin2, count
#' @export
#' @examples
#' hic_data <- create_test_hic_data()
#' head(hic_data)
create_test_hic_data <- function(chr = "chr1",
                                start = 1000000,
                                end = 1100000,
                                bin_size = 5000,
                                enrichment_factor = 10,
                                distance_decay = 0.001) {

  # Create bins
  bins <- seq(start, end, by = bin_size)
  bin_starts <- bins[1:(length(bins) - 1)]
  bin_ends <- bins[2:length(bins)]
  n_bins <- length(bin_starts)

  result <- data.frame(
    bin1 = integer(),
    bin2 = integer(),
    count = numeric(),
    stringsAsFactors = FALSE
  )

  # Generate contacts with distance-dependent decay
  for (i in seq_len(n_bins)) {
    for (j in seq_len(n_bins)) {
      if (j >= i) {  # Upper triangle only
        # Calculate genomic distance
        distance <- abs(bin_starts[j] - bin_starts[i])

        # Base count with distance decay
        base_count <- exp(-distance_decay * distance)

        # Add diagonal enrichment
        if (i == j) {
          base_count <- base_count * enrichment_factor
        } else if (abs(i - j) <= 2) {
          base_count <- base_count * (enrichment_factor * 0.5)
        }

        # Add some noise
        final_count <- max(0, rpois(1, base_count))

        if (final_count > 0) {
          contact_row <- data.frame(
            bin1 = i,
            bin2 = j,
            count = final_count,
            stringsAsFactors = FALSE
          )

          result <- rbind(result, contact_row)
        }
      }
    }
  }

  return(result)
}

#' Create comprehensive test dataset
#' @description Generates a complete set of test data for integration testing
#' @param region Genomic region string (default: "chr1:1000000-1100000")
#' @param include_all Whether to include all data types (default: TRUE)
#' @return List containing all test datasets
#' @export
#' @examples
#' test_data <- create_comprehensive_test_data()
#' names(test_data)
create_comprehensive_test_data <- function(region = "chr1:1000000-1100000",
                                          include_all = TRUE) {

  # Parse region
  region_parts <- strsplit(region, ":")[[1]]
  chr <- region_parts[1]
  range_parts <- strsplit(region_parts[2], "-")[[1]]
  start <- as.integer(range_parts[1])
  end <- as.integer(range_parts[2])

  test_data <- list(
    genomic_signal = create_test_genomic_data(chr, start, end),
    genes = create_test_gene_data(chr, gene_count = 3),
    signal_tracks = create_test_signal_data(chr, start, end, n_samples = 2),
    interactions = create_test_interaction_data(chr, start, end),
    hic_data = create_test_hic_data(chr, start, end)
  )

  return(test_data)
}