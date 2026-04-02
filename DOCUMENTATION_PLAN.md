# ezGenomeTracks Documentation Plan

This document outlines the plan for creating a complete pkgdown website with comprehensive vignettes for every track type and feature.

## Summary

### Current State
- **Existing vignettes**: get_started, coverage, gene, sashimi, manhattan
- **Missing vignettes**: feature (peaks), link (arcs), hic, composing tracks, showcase, landing page

### Proposed Structure

```
Home (index.Rmd)                    <- NEW: Independent landing page
├── Get Started                      <- Existing: get_started.Rmd
├── Tracks (dropdown)
│   ├── Coverage Track               <- Existing: tracks_coverage.Rmd (needs advanced section)
│   ├── Gene Track                   <- Existing: tracks_gene.Rmd (needs expansion)
│   ├── Feature/Peak Track           <- NEW: tracks_feature.Rmd
│   ├── Link/Arc Track               <- NEW: tracks_link.Rmd
│   ├── Sashimi Track                <- Existing: tracks_sashimi.Rmd (needs advanced section)
│   ├── Manhattan Track              <- Existing: tracks_manhattan.Rmd (needs advanced section)
│   └── Hi-C Track                   <- NEW: tracks_hic.Rmd
├── Examples (dropdown)
│   ├── Composing Tracks             <- NEW: composing_tracks.Rmd
│   └── Showcase: plotgardenerData   <- NEW: showcase_plotgardener.Rmd
└── Reference                        <- Auto-generated from roxygen
```

---

## 1. Landing Page (`index.Rmd`)

**Location**: `vignettes/articles/index.Rmd` (or configure pkgdown to use a custom home)

**Purpose**: A visually appealing, independent landing page (different from README.md)

### Content Structure

```markdown
# ezGenomeTracks

**Minimal, ggplot2-based genome browser tracks for R**

## Hero Section
- Package logo/banner (optional)
- One-liner tagline
- Installation code snippet
- Quick example with output image

## Why ezGenomeTracks?
- Comparison with other tools (pyGenomeTracks, Gviz, ggbio)
- Key differentiators:
  - Pure ggplot2 foundation
  - ez_* wrappers for quick plots
  - geom_* layers for customization
  - Seamless vertical stacking

## At a Glance
Gallery grid showing 6 track types with small example images:
- Coverage | Gene | Feature
- Link/Arc | Sashimi | Hi-C | Manhattan

## Quick Start
```r
library(ezGenomeTracks)
region <- "chr1:100000-200000"

# Create individual tracks
cov <- ez_coverage("signal.bw", region)
genes <- ez_gene(txdb, region)

# Stack them
vstack_plot(cov, genes, region = region)
```

## Learn More
- Links to vignettes
- Links to reference documentation
```

---

## 2. Track Vignettes

### Standard Structure for Each Track Vignette

Every track vignette should follow this template:

```markdown
---
title: "[Track Type] Track"
---

# Introduction
- What this track type is for
- Common use cases
- Input data types supported

# Basic Usage
## Loading Example Data
## Creating a Single Track
## Customizing Colors and Appearance

# Multiple Tracks
## Stacked Tracks (faceted)
## Overlapping Tracks (same panel)

# Advanced Usage
## Using geom_* directly
## Custom themes and scales
## Integration with other ggplot2 layers

# Tips and Best Practices
- Performance considerations
- Common pitfalls
```

---

### 2.1 Coverage Track (`tracks_coverage.Rmd`) - ENHANCE

**Status**: Exists, needs advanced section

**Add**:
- [ ] Section on `average = TRUE` for averaging multiple bigWigs
- [ ] Section on `type = "heatmap"` option
- [ ] Advanced: Using `geom_coverage()` with custom aesthetics
- [ ] Advanced: Adding annotations (vlines, rectangles)
- [ ] Tips: Memory considerations for large bigWigs

---

### 2.2 Gene Track (`tracks_gene.Rmd`) - EXPAND

**Status**: Exists but minimal

**Add**:
- [ ] Basic: From GTF/GFF files
- [ ] Basic: From data frames
- [ ] Customization: Exon/intron colors
- [ ] Customization: Strand-based coloring
- [ ] Labels: `label`, `label_style`, `max_labels`, `label_priority`
- [ ] Labels: Using ggrepel for smart positioning
- [ ] Advanced: Using `geom_gene()` directly
- [ ] Advanced: Custom y-axis grouping (not just strand)

---

### 2.3 Feature/Peak Track (`tracks_feature.Rmd`) - NEW

**Status**: Does not exist

**Content**:

```markdown
---
title: "Feature/Peak Track"
---

# Introduction
- Visualizing genomic intervals (peaks, enhancers, promoters)
- Input: BED files, data frames

# Basic Usage
## From a BED file
```r
ez_feature("peaks.bed", region = "chr1:1000000-2000000")
```

## From a data frame
```r
peaks_df <- data.frame(
  seqnames = "chr1",
  start = c(1000, 3000, 5000),
  end = c(2000, 4000, 6000),
  score = c(10, 30, 50)
)
ez_feature(peaks_df, region, fill = "darkgreen")
```

# Customization
## Score-based coloring
```r
ez_feature(peaks_df, region, use_score = TRUE, fill = "red")
```

## Adjusting height and transparency
```r
ez_feature(peaks_df, region, height = 0.6, alpha = 0.5)
```

# Advanced Usage
## Using geom_feature()
```r
ggplot(peaks_df) +
  geom_feature(aes(fill = score), height = 0.8) +
  scale_fill_gradient(low = "white", high = "red") +
  scale_x_genome_region(region) +
  theme_void()
```

# Stacking with Other Tracks
```r
vstack_plot(
  ez_coverage(bw_file, region),
  ez_feature(peaks_df, region),
  ez_gene(txdb, region),
  region = region
)
```
```

---

### 2.4 Link/Arc Track (`tracks_link.Rmd`) - NEW

**Status**: Does not exist

**Content**:

```markdown
---
title: "Link/Arc Track"
---

# Introduction
- Visualizing genomic interactions (chromatin loops, enhancer-promoter)
- Input: BEDPE format, data frames, GInteractions

# Basic Usage
## Creating interaction data
```r
interactions_df <- data.frame(
  start1 = c(1000, 5000),
  end1 = c(1500, 5500),
  start2 = c(3000, 8000),
  end2 = c(3500, 8500),
  score = c(10, 20)
)
ez_link(interactions_df, region = "chr1:1-10000")
```

# Customization
## Arc direction
```r
ez_link(interactions_df, region, direction = "up")
ez_link(interactions_df, region, direction = "down")
```

## Score-based coloring
```r
ez_link(interactions_df, region, use_score = TRUE, colors = "purple")
```

## Curvature and height
```r
ez_link(interactions_df, region, curvature = 0.7, height_factor = 0.2)
```

# Multiple Tracks
## Grouped interactions
```r
ez_link(interactions_df, region, group_var = "cell_type")
```

## Stacked tracks (faceted)
```r
ez_link(list("K562" = k562_loops, "GM12878" = gm_loops), region)
```

# Advanced Usage
## Using geom_link() directly
```r
ggplot(interactions_df) +
  geom_link(
    aes(start1 = start1, end1 = end1, start2 = start2, end2 = end2, color = score),
    direction = "down",
    height_factor = 0.15
  ) +
  scale_color_viridis_c() +
  scale_x_genome_region(region)
```
```

---

### 2.5 Sashimi Track (`tracks_sashimi.Rmd`) - ENHANCE

**Status**: Exists, needs advanced section

**Add**:
- [ ] Explanation of input format (coverage + junction files)
- [ ] Customization: `junction_direction = "alternate"/"up"/"down"`
- [ ] Customization: `score_transform` for wide-range scores
- [ ] Advanced: Junction filtering by score
- [ ] Advanced: Custom label formatting
- [ ] Tips: Preparing junction files from BAM

---

### 2.6 Manhattan Track (`tracks_manhattan.Rmd`) - ENHANCE

**Status**: Exists, needs expansion

**Add**:
- [ ] Genome-wide mode: chromosome coloring options
- [ ] Regional mode: LocusZoom-style with r2
- [ ] Highlighting SNPs: `highlight_snps`, `lead_snp`
- [ ] Significance threshold lines
- [ ] Multiple traits: grouped Manhattan
- [ ] Advanced: Using `geom_manhattan()` directly
- [ ] Tips: Column name auto-detection

---

### 2.7 Hi-C Track (`tracks_hic.Rmd`) - NEW

**Status**: Does not exist

**Content**:

```markdown
---
title: "Hi-C Track"
---

# Introduction
- Visualizing 3D chromatin organization
- Input: Contact matrices, sparse data frames

# Basic Usage
## From a matrix
```r
mat <- matrix(runif(2500), nrow = 50)
mat <- mat + t(mat)  # Make symmetric
ez_hic(mat, region = "chr1:1000000-1500000", resolution = 10000)
```

## From a sparse data frame
```r
hic_df <- data.frame(
  pos1 = c(100000, 100000, 110000),
  pos2 = c(110000, 120000, 120000),
  score = c(50, 30, 45)
)
ez_hic(hic_df, region, resolution = 10000)
```

# Visualization Styles
## Triangle view (default)
```r
ez_hic(mat, region, style = "triangle")
```

## Square heatmap view
```r
ez_hic(mat, region, style = "square")
```

# Customization
## Color palettes
```r
ez_hic(mat, region, palette = "cooler")   # Red (default)
ez_hic(mat, region, palette = "ylgnbu")   # Yellow-Green-Blue
ez_hic(mat, region, palette = "viridis")  # Viridis
ez_hic(mat, region, palette = "bwr")      # Blue-White-Red
```

## Scale transformation
```r
ez_hic(mat, region, transform = "log10")
ez_hic(mat, region, transform = "sqrt")
```

## Maximum distance filter
```r
ez_hic(mat, region, max_distance = 200000)
```

# Advanced Usage
## Rasterization for large matrices
```r
ez_hic(mat, region, rasterize = TRUE, rasterize_dpi = 300)
```

## Using geom_hic_triangle() directly
```r
ggplot(hic_df) +
  geom_hic_triangle(aes(x = pos1, y = pos2, fill = score), resolution = 10000) +
  scale_fill_hic(palette = "cooler", trans = "log10") +
  scale_x_genome_region(region)
```
```

---

## 3. Composing Tracks (`composing_tracks.Rmd`) - NEW

**Purpose**: How to combine multiple tracks into publication-quality figures

**Content**:

```markdown
---
title: "Composing Tracks"
---

# Introduction
- The philosophy of modular track composition
- When to use vstack_plot vs hstack_plot

# Vertical Stacking with vstack_plot()
## Basic stacking
```r
cov <- ez_coverage(bw_file, region)
genes <- ez_gene(txdb, region)
peaks <- ez_feature(peaks_df, region)

vstack_plot(cov, genes, peaks, region = region)
```

## Controlling track heights
```r
vstack_plot(cov, genes, peaks, region = region, heights = c(2, 1, 0.5))
```

# Horizontal Stacking with hstack_plot()
## Comparing regions
```r
region1 <- "chr1:1000000-2000000"
region2 <- "chr2:5000000-6000000"

p1 <- ez_coverage(bw, region1)
p2 <- ez_coverage(bw, region2)

hstack_plot(p1, p2, widths = c(1, 1))
```

# Adding Annotations
## Vertical lines
```r
track <- ez_coverage(bw, region)
track <- add_vline(track, position = 1500000, color = "red")
```

## Highlight regions
```r
track <- add_rect(track, xmin = 1200000, xmax = 1400000, fill = "yellow")
```

## Text labels
```r
track <- add_text(track, x = 1300000, y = 50, label = "Peak")
```

# Real-World Examples
## Gene locus browser
## QTL visualization
## Epigenomic landscape
```

---

## 4. Showcase: plotgardenerData (`showcase_plotgardener.Rmd`) - NEW

**Purpose**: Demonstrate ezGenomeTracks with realistic, publication-quality data from Bioconductor

**Why plotgardenerData?**
- Provides real experimental data (not synthetic)
- Includes multiple data types that showcase all track types
- Users can reproduce examples exactly
- Establishes comparison with plotgardener package

### Available Datasets from plotgardenerData

| Dataset | Type | Cell Line | Description |
|---------|------|-----------|-------------|
| `IMR90_HiC_10kb` | Hi-C | IMR90 | 10kb resolution contact matrix |
| `GM12878_HiC_10kb` | Hi-C | GM12878 | 10kb resolution contact matrix |
| `IMR90_DNAloops_pairs` | Loops/BEDPE | IMR90 | Chromatin loop calls |
| `IMR90_ChIP_H3K27ac_signal` | Signal | IMR90 | H3K27ac ChIP-seq |
| `GM12878_ChIP_H3K27ac_signal` | Signal | GM12878 | H3K27ac ChIP-seq |
| `GM12878_ChIP_CTCF_signal` | Signal | GM12878 | CTCF ChIP-seq |
| `IMR90_ChIP_CTCF_signal` | Signal | IMR90 | CTCF ChIP-seq |
| `hg19_insulin_GWAS` | GWAS | - | Insulin-related GWAS results |

**Files in extdata:**
- `test.bw` - Example bigWig file
- `test_chr22.hic` - Example .hic file

**Common Region**: `chr21:28000000-30300000` (hg19)

### Content Structure

```markdown
---
title: "Showcase: Multi-omic Visualization with plotgardenerData"
---

# Introduction
- About plotgardenerData (Bioconductor package)
- What we'll build: A publication-ready multi-omic browser view
- Prerequisites: installing plotgardenerData

# Setup

```r
# Install plotgardenerData from Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("plotgardenerData")

# Load packages
library(ezGenomeTracks)
library(plotgardenerData)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

# Define common region
region <- "chr21:28000000-30300000"
```

# Example 1: Comparing Cell Types (Hi-C + Signal)

Visualize chromatin architecture differences between IMR90 and GM12878.

```r
# Load Hi-C data
data("IMR90_HiC_10kb")
data("GM12878_HiC_10kb")

# Create Hi-C tracks
hic_imr90 <- ez_hic(IMR90_HiC_10kb, region, resolution = 10000,
                    style = "triangle", palette = "cooler")
hic_gm12878 <- ez_hic(GM12878_HiC_10kb, region, resolution = 10000,
                      style = "triangle", palette = "cooler")

# Load signal data
data("IMR90_ChIP_H3K27ac_signal")
data("GM12878_ChIP_H3K27ac_signal")

# Create signal tracks
sig_imr90 <- ez_coverage(IMR90_ChIP_H3K27ac_signal, region,
                         colors = "orange2", y_axis_style = "simple")
sig_gm12878 <- ez_coverage(GM12878_ChIP_H3K27ac_signal, region,
                           colors = "purple2", y_axis_style = "simple")

# Gene track
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
genes <- ez_gene(txdb, region)

# Stack: IMR90 view
vstack_plot(
  hic_imr90,
  sig_imr90,
  genes,
  region = region,
  heights = c(2, 1, 0.5)
)
```

# Example 2: Chromatin Loops with Signal Context

Overlay loop calls on Hi-C and correlate with H3K27ac signal.

```r
# Load loop data
data("IMR90_DNAloops_pairs")

# Create tracks
hic_track <- ez_hic(IMR90_HiC_10kb, region, resolution = 10000)
loop_track <- ez_link(IMR90_DNAloops_pairs, region,
                      direction = "down", colors = "blue")
signal_track <- ez_coverage(IMR90_ChIP_H3K27ac_signal, region)
gene_track <- ez_gene(txdb, region)

# Composed view
vstack_plot(
  hic_track,
  loop_track,
  signal_track,
  gene_track,
  region = region,
  heights = c(2, 0.8, 1, 0.5)
)
```

# Example 3: GWAS Integration

Combine genetic association data with epigenomic context.

```r
# Load GWAS data
data("hg19_insulin_GWAS")

# Create Manhattan track (regional view)
gwas_region <- "chr21:28500000-29500000"
manhattan_track <- ez_manhattan(hg19_insulin_GWAS, region = gwas_region,
                                 y_axis_style = "full")

# Signal context
h3k27ac_track <- ez_coverage(IMR90_ChIP_H3K27ac_signal, gwas_region)
gene_track <- ez_gene(txdb, gwas_region)

# Integrated QTL-style view
vstack_plot(
  manhattan_track,
  h3k27ac_track,
  gene_track,
  region = gwas_region,
  heights = c(1.5, 1, 0.5)
)
```

# Example 4: Side-by-Side Cell Type Comparison

Use hstack_plot to compare IMR90 vs GM12878 at the same locus.

```r
# IMR90 stack
imr90_view <- vstack_plot(

  ez_hic(IMR90_HiC_10kb, region, resolution = 10000),
  ez_coverage(IMR90_ChIP_H3K27ac_signal, region, colors = "orange"),
  ez_gene(txdb, region),
  region = region,
  heights = c(2, 1, 0.5)
)

# GM12878 stack
gm12878_view <- vstack_plot(
  ez_hic(GM12878_HiC_10kb, region, resolution = 10000),
  ez_coverage(GM12878_ChIP_H3K27ac_signal, region, colors = "purple"),
  ez_gene(txdb, region),
  region = region,
  heights = c(2, 1, 0.5)
)

# Side by side
hstack_plot(imr90_view, gm12878_view, widths = c(1, 1))
```

# Example 5: Complete Multi-omic Browser

Build a comprehensive view with all track types.

```r
region <- "chr21:28500000-29500000"

# All tracks
p_hic <- ez_hic(IMR90_HiC_10kb, region, resolution = 10000,
                style = "triangle", max_distance = 500000)
p_loops <- ez_link(IMR90_DNAloops_pairs, region, direction = "down")
p_h3k27ac <- ez_coverage(IMR90_ChIP_H3K27ac_signal, region,
                          colors = "orange2", y_axis_style = "minmax")
p_ctcf <- ez_coverage(IMR90_ChIP_CTCF_signal, region,
                       colors = "green4", y_axis_style = "minmax")
p_genes <- ez_gene(txdb, region, label = "gene_name")

# Final composition
vstack_plot(
  p_hic,
  p_loops,
  p_h3k27ac,
  p_ctcf,
  p_genes,
  region = region,
  heights = c(2, 0.6, 1, 1, 0.8)
)
```

# Tips for Working with plotgardenerData

1. **Assembly**: All data uses hg19 - ensure your TxDb matches
2. **Resolution**: Hi-C data is at 10kb resolution
3. **Memory**: Hi-C matrices can be large; use `max_distance` to limit
4. **Regions**: Start with chr21:28000000-30300000 for best coverage
```

---

## 5. Updated `_pkgdown.yml`

```yaml
destination: docs

url: https://zepeng-mu.github.io/ezGenomeTracks/

home:
  title: ezGenomeTracks
  description: A ggplot2-based R package for creating publication-quality genomic track visualizations

template:
  bootstrap: 5
  bootswatch: journal

articles:
- title: Getting Started
  navbar: ~
  contents:
  - get_started

- title: Track Types
  desc: Documentation for each track type
  contents:
  - tracks_coverage
  - tracks_gene
  - tracks_feature
  - tracks_link
  - tracks_sashimi
  - tracks_manhattan
  - tracks_hic

- title: Advanced Topics
  desc: Composing and customizing track visualizations
  contents:
  - composing_tracks

- title: Examples & Showcases
  desc: Real-world examples with realistic data
  contents:
  - showcase_plotgardener

reference:
- title: "Easy Wrappers"
  desc: >-
    High-level functions for quick track creation
  contents:
  - ez_coverage
  - ez_gene
  - ez_feature
  - ez_link
  - ez_sashimi
  - ez_manhattan
  - ez_hic

- title: "Track Composition"
  desc: "Functions for combining multiple tracks"
  contents:
  - vstack_plot
  - hstack_plot

- title: "Annotation Helpers"
  desc: "Functions for adding annotations to tracks"
  contents:
  - add_vline
  - add_hline
  - add_rect
  - add_text

- title: "Geom Functions"
  desc: "Low-level geometric objects for genomic data"
  contents:
  - geom_coverage
  - geom_feature
  - geom_gene
  - geom_link
  - geom_manhattan
  - geom_hic
  - geom_hic_triangle

- title: "Data Processing Functions"
  desc: "Functions for importing and processing track data"
  contents:
  - process_signal_input
  - process_gene_data
  - process_interaction_data
  - process_interaction_input
  - process_sashimi_data
  - process_sashimi_input
  - process_manhattan_input
  - process_hic_data
  - import_genomic_data
  - average_signal
  - get_single_signal

- title: "Scales and Themes"
  desc: "Functions for customizing track appearance"
  contents:
  - scale_x_genome
  - scale_x_genome_region
  - scale_y_strand
  - scale_fill_hic
  - hic_palettes
  - ez_theme
  - ez_coverage_theme
  - ez_feature_theme
  - ez_gene_theme
  - ez_sashimi_theme
  - ez_manhattan_theme
  - ez_hic_theme

- title: "Helper Functions"
  desc: "Utility functions for data conversion and calculation"
  contents:
  - granges_to_df
  - df_to_granges
  - parse_region
  - format_genomic_coord
  - extract_txdb_data
  - calculate_link_ylim
  - stat_bin_signal
  - plot_signal_df

- title: "Color Palettes"
  desc: "Color utilities for genomic visualizations"
  contents:
  - ez_default_palette
  - ez_default_single_color
  - ez_hic_palette
  - ez_hic_diverging_palette

- title: "Example Datasets"
  desc: "Example datasets for demonstration"
  contents:
  - example_signal
  - example_peaks
  - example_genes
  - example_interactions
  - example_junctions
  - example_hic

navbar:
  structure:
    left: [intro, tracks, examples, reference, news]
    right: [github]
  components:
    intro:
      text: Get Started
      href: articles/get_started.html
    tracks:
      text: Tracks
      menu:
      - text: Coverage Track
        href: articles/tracks_coverage.html
      - text: Gene Track
        href: articles/tracks_gene.html
      - text: Feature/Peak Track
        href: articles/tracks_feature.html
      - text: Link/Arc Track
        href: articles/tracks_link.html
      - text: "---"
      - text: Sashimi Track
        href: articles/tracks_sashimi.html
      - text: Manhattan Track
        href: articles/tracks_manhattan.html
      - text: Hi-C Track
        href: articles/tracks_hic.html
    examples:
      text: Examples
      menu:
      - text: Composing Tracks
        href: articles/composing_tracks.html
      - text: "---"
      - text: "Showcase: plotgardenerData"
        href: articles/showcase_plotgardener.html
    github:
      icon: fab fa-github fa-lg
      href: https://github.com/Zepeng-Mu/ezGenomeTracks
```

---

## 6. Work Checklist

### Phase 1: New Vignettes (Priority: High)
- [x] Create `vignettes/articles/tracks_feature.Rmd`
- [x] Create `vignettes/articles/tracks_link.Rmd`
- [x] Create `vignettes/articles/tracks_hic.Rmd`
- [x] Create `vignettes/articles/composing_tracks.Rmd`
- [x] Create `vignettes/articles/showcase_plotgardener.Rmd`

### Phase 2: Enhance Existing Vignettes (Priority: Medium)
- [x] Expand `tracks_gene.Rmd` with label options, custom colors
- [x] Add advanced section to `tracks_coverage.Rmd` (averaging, heatmap)
- [x] Add advanced section to `tracks_sashimi.Rmd`
- [x] Add advanced section to `tracks_manhattan.Rmd`

### Phase 3: Landing Page (Priority: Medium)
- [x] Create custom `index.Rmd` for pkgdown home
- [x] Generate gallery images for each track type (quick example + multi-sample)
- [x] Update `_pkgdown.yml` to use custom home

### Phase 4: Polish (Priority: Low)
- [ ] Ensure example data available for all vignettes
- [ ] Add plotgardenerData to Suggests in DESCRIPTION
- [ ] Add gallery/screenshots section
- [ ] Review and update function documentation
- [ ] Test all vignette code chunks

---

## Notes

1. **Example Data**: Ensure `inst/extdata/` contains sample files for all track types, or use `data()` objects from the package.

2. **plotgardenerData**: The showcase vignette requires `plotgardenerData` from Bioconductor. Add to `Suggests` in DESCRIPTION and document installation in the vignette.

3. **Consistent Styling**: All vignettes should use the same region where possible to show how tracks stack together. For plotgardenerData examples, use `chr21:28000000-30300000` (hg19).

4. **Real-World Examples**: The showcase page demonstrates complete workflows (comparing cell types, integrating GWAS with epigenomics, multi-omic browser).

5. **Performance Notes**: Add tips for handling large files (bigWig, Hi-C matrices) in each relevant vignette.

6. **Comparison with plotgardener**: The showcase implicitly demonstrates that ezGenomeTracks can reproduce similar visualizations to plotgardener using pure ggplot2 syntax, highlighting the different design philosophies.
