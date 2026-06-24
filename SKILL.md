---
name: ezgenometracks
description: >
  Helper for writing R code with the ezGenomeTracks package — a ggplot2-based genome browser track library.
  Use this skill whenever the user is working with ezGenomeTracks, asks how to plot genomic data in R
  (coverage, genes, peaks, Hi-C, interactions, Manhattan, sashimi), wants to stack genome tracks, or
  asks about any of its functions (ez_coverage, ez_gene, ez_feature, ez_link, ez_hic, ez_manhattan,
  ez_sashimi, vstack_plot). Also trigger when the user mentions bigWig, bedGraph, GTF, BEDPE, or
  TxDb objects in the context of plotting in R. Do NOT wait for the user to explicitly say
  "use the skill" — if the task involves genome track visualization in R, use this skill proactively.
---

# ezGenomeTracks Helper

`ezGenomeTracks` is a ggplot2-based R package for genome browser-style track visualization.
It provides `ez_*()` wrappers (quick, high-level) and `geom_*()` layers (full ggplot2 control).

## Installation

```r
# install.packages("pak")
pak::pak("Zepeng-Mu/ezGenomeTracks")
```

---

## Core Concepts

### 1. Region string
Every function needs a region: `"chr1:1000000-2000000"`. All tracks in the same figure must use the **same** region string so their x-axes align.

### 2. Input types
All `ez_*()` functions accept flexible input:
- A **data frame** with genomic columns
- A **file path** (BED, bedGraph, bigWig, GTF, BEDPE …)
- A **named list** of the above — creates multiple faceted sub-tracks

### 3. Required column names
Follow GRanges-style naming:
| Column | Meaning |
|---|---|
| `seqnames` | chromosome, e.g. `"chr1"` |
| `start` | start position (1-based) |
| `end` | end position |
| `strand` | `"+"`, `"-"`, or `"Unknown"` — must be a **factor** with levels `c("minus","plus","Unknown")` for gene tracks |
| `score` | numeric signal value (coverage, peak score, etc.) |

### 4. All tracks are ggplot2 objects
`ez_*()` returns a `ggplot` object — you can add ggplot2 layers (`+ theme(...)`, `+ labs(...)`) on top of any track.

---

## Track Types & `ez_*` Wrappers

### Coverage / Signal Track — `ez_coverage()`

Plot a signal over intervals (bigWig, bedGraph, or a data frame with `seqnames/start/end/score`).

```r
library(ezGenomeTracks)
region <- "chr1:100000-101000"

# From a data frame
cov_df <- data.frame(
  seqnames = "chr1",
  start = seq(100000, 100990, by = 10),
  end   = seq(100009, 100999, by = 10),
  score = rnorm(100, 10, 2)
)
track_cov <- ez_coverage(cov_df, region = region)

# From a file path
track_cov <- ez_coverage("path/to/signal.bw", region = region)

# Multiple tracks via named list
track_cov <- ez_coverage(
  list(Sample1 = "s1.bw", Sample2 = "s2.bw"),
  region = region,
  colors = c("steelblue", "tomato")
)
```

**Key parameters:**
- `type`: `"area"` (default), `"line"`, or `"heatmap"`
- `colors`: single color or vector of colors for multiple tracks/groups
- `group_var`: column name to color by within a single data frame
- `y_axis_style`: `"none"` (default) | `"simple"` | `"minmax"` | `"full"`
- `alpha`: transparency (default `0.5`)
- `average = TRUE`: average multiple files into one track
- `y_range`: fix y-axis limits, e.g. `c(0, 50)`

---

### Gene / Transcript Track — `ez_gene()`

Plot gene models from a GTF file, TxDb object, or a data frame.

```r
# From a data frame (must have seqnames, start, end, strand, type, gene_id, transcript_id)
genes_df <- data.frame(
  seqnames      = "chr1",
  start         = c(100050, 100200),
  end           = c(100150, 100350),
  strand        = factor(c("+", "+"), levels = c("-", "+", "Unknown")),
  type          = "exon",
  gene_id       = "GENE1",
  transcript_id = "TX1",
  gene_name     = "GENE1"
)
track_gene <- ez_gene(genes_df, region = region)

# From a GTF file (auto-imports with region filtering)
track_gene <- ez_gene("genes.gtf", region = region)

# From a TxDb object
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
track_gene <- ez_gene(TxDb.Hsapiens.UCSC.hg38.knownGene, region = region)

# Look up region by gene name instead of coordinates
track_gene <- ez_gene(txdb, gene = "PTPRC", gene_db = txdb)
```

**Key parameters:**
- `label`: column to label genes with (default `"gene_name"`); set `label = NULL` to hide
- `label_style`: `"auto"` | `"repel"` | `"simple"` | `"none"`
- `max_labels`: limit number of labels shown
- `exon_fill`, `exon_color`, `intron_color`: explicit colors (bypasses strand-based coloring)
- Default colors: plus strand = `"green4"`, minus strand = `"orange2"`

**Strand factor levels pitfall:** If you build the data frame manually, use:
```r
strand = factor(strand_vec, levels = c("-", "+", "Unknown"))
```

---

### Feature / Peak Track — `ez_feature()`

Plot genomic intervals as rectangles (e.g. ChIP-seq peaks, regulatory elements).

```r
features_df <- data.frame(
  seqnames = "chr1",
  start    = c(100120, 100300),
  end      = c(100180, 100360),
  score    = c(12, 8)
)
track_feat <- ez_feature(features_df, region = region)

# Color by score (gradient)
track_feat <- ez_feature(features_df, region = region,
                         use_score = TRUE, fill = "darkred")

# From a BED file
track_feat <- ez_feature("peaks.bed", region = region)
```

**Key parameters:**
- `fill`: fill color (or high end of gradient when `use_score = TRUE`)
- `color`: border color (default `"black"`)
- `alpha`: transparency (default `0.7`)
- `height`: rectangle height 0–1 (default `0.8`)
- `use_score`: color-code by `score` column (default `FALSE`)

---

### Interaction / Arc Track — `ez_link()`

Draw arcs between pairs of genomic intervals (enhancer–promoter contacts, etc.).

Input data frame needs columns: `chr1`, `start1`, `end1`, `chr2`, `start2`, `end2`, optionally `score`.

```r
interactions_df <- data.frame(
  chr1 = "chr1", start1 = 100120, end1 = 100160,
  chr2 = "chr1", start2 = 100300, end2 = 100340,
  score = 5
)
track_arc <- ez_link(interactions_df, region = region)

# Score-based coloring
track_arc <- ez_link(interactions_df, region = region,
                     use_score = TRUE, colors = "steelblue")
```

**Key parameters:**
- `direction`: `"down"` (default) | `"up"`
- `curvature`: arc curvature 0–1 (default `0.5`)
- `size`: line width (default `0.5`)
- `use_score`: color arcs by score
- `colors`: color or vector of colors

---

### Manhattan Track — `ez_manhattan()`

Plot association statistics (e.g. GWAS -log10 p-values).

```r
manhattan_df <- data.frame(
  seqnames = rep("chr1", 50),
  start    = seq(100000, 100000 + 49 * 200, by = 200),
  end      = seq(100001, 100001 + 49 * 200, by = 200),
  score    = rexp(50, rate = 0.3)
)
track_manh <- ez_manhattan(manhattan_df, region = region)
```

---

### Hi-C Contact Matrix — `ez_hic()`

Plot binned pairwise contact frequencies as a triangular heatmap (default) or square matrix.

Input: sparse data frame with `seqnames`, `start1`, `start2`, `score`; or a dense matrix; or a file path.

```r
hic_df <- tidyr::expand_grid(
  bin1 = seq(100000, 100000 + 4 * 200, by = 200),
  bin2 = seq(100000, 100000 + 4 * 200, by = 200)
) |>
  dplyr::filter(bin2 >= bin1) |>
  dplyr::mutate(score = runif(dplyr::n()), seqnames = "chr1") |>
  dplyr::rename(start1 = bin1, start2 = bin2)

track_hic <- ez_hic(hic_df, region = region, resolution = 200)

# Square view
track_hic <- ez_hic(hic_df, region = region, resolution = 200, style = "square")
```

**Key parameters:**
- `resolution`: bin size in bp (default `10000`)
- `style`: `"triangle"` (default) | `"square"`
- `palette`: `"cooler"` (default) | `"ylgnbu"` | `"viridis"` | `"bwr"`
- `transform`: `"identity"` | `"log10"` (default) | `"log2"` | `"sqrt"`
- `max_distance`: trim triangle to a max interaction distance (bp)

---

### Sashimi Plot — `ez_sashimi()`

RNA-seq coverage + splice junction arcs in one track.

```r
coverage <- data.frame(
  seqnames = "chr1", start = 1:100, end = 2:101,
  score = c(runif(30, 5, 10), rep(0, 20), runif(50, 5, 10))
)
junctions <- data.frame(
  seqnames = "chr1", start = 30, end = 50, score = 15
)
track_sashimi <- ez_sashimi(coverage, junctions, "chr1:1-100")

# Multi-sample
track_sashimi <- ez_sashimi(
  coverage_data = list(S1 = cov1, S2 = cov2),
  junction_data = list(S1 = junc1, S2 = junc2),
  region = region,
  colors = c("purple", "orange")
)
```

**Key parameters:**
- `junction_direction`: `"alternate"` (default) | `"up"` | `"down"`
- `show_labels`: show arc read-count labels (default `TRUE`)
- `score_transform`: `"identity"` | `"log10"` | `"sqrt"` (scale arc thickness)

---

## Composing Multiple Tracks

### Vertical stacking — `vstack_plot()`

All tracks must use the **same region**. Pass them as a list or as `...` arguments.

```r
combined <- vstack_plot(
  list(
    coverage   = track_cov,
    genes      = track_gene,
    peaks      = track_feat,
    arcs       = track_arc
  ),
  region  = region,
  heights = c(2, 1, 1, 1)   # one value per additional track (length = n_tracks - 1)
)
combined
```

> **Important:** `heights` must have length **n_tracks − 1**, not n_tracks.

### Horizontal stacking — `hstack_plot()`

Side-by-side regions (no shared x-axis required).

```r
combined_h <- hstack_plot(track_cov, track_gene, widths = c(2, 1))
```

---

## Annotation Helpers

Add lines, rectangles, and text to any track:

```r
track <- add_vline(track_cov, position = 100500, color = "red", linetype = "dashed")
track <- add_hline(track_cov, y = 15, color = "blue")
track <- add_rect(track_cov, xmin = 100200, xmax = 100400, fill = "yellow", alpha = 0.2)
track <- add_text(track_cov, x = 100350, y = 20, label = "peak")
```

---

## Using Gene Name Instead of Coordinates

Any `ez_*()` function accepts `gene =` + `gene_db =` to auto-resolve the region:

```r
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

track_cov  <- ez_coverage(signal_df, gene = "PTPRC", gene_db = txdb)
track_gene <- ez_gene(txdb, gene = "PTPRC", gene_db = txdb)

# Extend the window by 20%
track_cov <- ez_coverage(signal_df, gene = "PTPRC", gene_db = txdb, extend = 0.2)
# Or by exact bp
track_cov <- ez_coverage(signal_df, gene = "PTPRC", gene_db = txdb,
                          extend = 5000, extend_type = "bp")
```

---

## Common Patterns

### Quick multi-sample coverage panel

```r
region <- "chr6:29900000-30100000"

tracks <- ez_coverage(
  list(
    "Naive B"  = "naive_b.bw",
    "Memory B" = "memory_b.bw",
    "Plasma"   = "plasma.bw"
  ),
  region = region,
  colors = c("#4C72B0", "#DD8452", "#55A868"),
  y_axis_style = "simple"
)

genes <- ez_gene("hg38_genes.gtf", region = region, label = "gene_name")

vstack_plot(list(tracks, genes), region = region, heights = c(3, 1))
```

### Fix y-axis across samples for fair comparison

```r
t1 <- ez_coverage(df1, region = region, y_range = c(0, 50))
t2 <- ez_coverage(df2, region = region, y_range = c(0, 50))
vstack_plot(list(t1, t2), region = region, heights = c(1))
```

### Average biological replicates

```r
track_avg <- ez_coverage(
  c("rep1.bw", "rep2.bw", "rep3.bw"),
  region = region,
  average = TRUE,
  summary_fun = "mean"
)
```

---

## Tips & Common Pitfalls

- **Same region for all tracks** — `vstack_plot()` relies on this; mismatched regions will produce misaligned x-axes.
- **`heights` length** — must be `n_tracks − 1`. For 3 tracks pass `heights = c(2, 1)`.
- **Strand factor levels** — when building gene data frames manually: `factor(strand_col, levels = c("-", "+", "Unknown"))`.
- **Column names are exact** — `seqnames` not `chr`, `start` not `chromStart`. Wrappers look for these names specifically.
- **Large files** — always supply `region =` to enable region-filtered import; loading the full genome into R is slow.
- **ggplot2 composability** — every `ez_*()` result is a plain ggplot object. Add `+ theme(...)`, `+ ggtitle(...)`, etc. freely.
