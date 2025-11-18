# ezGenomeTracks AI Coding Guide

## Project Overview
`ezGenomeTracks` is a ggplot2-based R package for creating genome browser-style visualizations. It provides custom geoms and wrapper functions for plotting genomic data types (signal tracks, peaks, gene models, interactions, Hi-C) with a common genomic coordinate system, stacked using the `aplot` package.

## Architecture & Design Patterns

### Three-Tier API Design
1. **Low-level geoms** (`geom_*`): Pure ggplot2 geoms for direct data visualization
   - `geom_gene()`, `geom_coverage()`, `geom_feature()`, `geom_arc()`, `geom_hic()`, `geom_manhattan()`
   - Located in `R/geom_*.R` files
   - Use ggproto objects (e.g., `GeomGene`, `GeomCoverage`)

2. **Easy wrappers** (`ez_*`): High-level convenience functions that combine geoms with themes/scales
   - `ez_coverage()`, `ez_gene()`, `ez_feature()`, `ez_arc()`, `ez_manhattan()`
   - All in `R/ez_wrappers.R`
   - Handle multiple input types: file paths, data frames, named lists
   - Apply appropriate themes automatically

3. **Track composition** (`genome_plot()`): Stacks multiple tracks with aligned x-axes
   - Uses `aplot::plot_list()` for vertical stacking
   - Located in `R/genome_plot.R`

### Data Flow Pattern
```
Input (file/df/list) → process_*_data() → ez_*() → geom_*() → theme/scales → genome_plot()
```

**Critical helpers in `R/helpers.R`:**
- `parse_region()`: Converts "chr1:1000-2000" strings to GRanges objects
- `process_signal_input()`: Normalizes diverse input types (df/files/lists) into standardized format
- `granges_to_df()`/`df_to_granges()`: Bidirectional GRanges ↔ data frame conversion
- `import_genomic_data()`: Wraps rtracklayer::import for file reading

### Column Name Conventions
- **Coordinates**: `seqnames`, `start`, `end`, `strand` (GRanges-compatible)
- **Signal data**: `score` for numeric values
- **Gene data**: `gene_id`, `gene_name`, `transcript_id`, `type` ("gene"/"exon")
- **Interactions**: `chr1`, `start1`, `end1`, `chr2`, `start2`, `end2`
- **Required aesthetics**: Use `.data$` pronoun for NSE (e.g., `.data$start`)

## R Package Development Workflows

### Building & Testing
```bash
# Standard R package workflow (use R console or devtools)
devtools::load_all()           # Load package for testing
devtools::document()           # Update NAMESPACE and .Rd files
devtools::check()              # Run R CMD check
devtools::test()               # Run testthat tests
pkgdown::build_site()          # Build documentation site
```

### Documentation Standards
- **Roxygen2**: All exported functions must have `@export` tags
- Use `@importFrom` for specific imports (avoid `@import` for entire packages)
- Include `@examples` with `\dontrun{}` for functions requiring external data
- Document ggproto objects with `@rdname`, `@format NULL`, `@usage NULL`

### Testing Conventions
- Tests in `tests/testthat/test-*.R` follow pattern: `test-<source_file>.R`
- Use `testthat` edition 3 (`Config/testthat/edition: 3`)
- Mock/example data in `data/` (`.rda` files) and `inst/extdata/` (raw files)

## Critical Integration Points

### Bioconductor Dependencies
- **GenomicRanges/IRanges**: Core data structures; use `methods::is()` for type checking
- **rtracklayer**: File I/O; wrap imports in `which = region_gr` for region filtering
- **GenomicFeatures**: TxDb support (optional, check with `requireNamespace()`)

### ggplot2 Extension Mechanics
- **Custom geoms**: Inherit from `ggplot2::Geom`, implement `draw_panel()` and `setup_data()`
- **Aesthetics**: Define `required_aes`, `optional_aes`, and `default_aes` in ggproto
- **Discrete y-axes**: Use `ggplot2:::compute_data_size()` for proper height calculation (see `GeomGene`)

### aplot Stacking
- `genome_plot()` passes tracks to `aplot::plot_list()` with `guides = "collect"`
- All tracks MUST share the same x-axis limits (via `scale_x_genome_region()`)
- Use `heights` parameter for relative track sizing

## Project-Specific Conventions

### Function Naming
- `ez_*`: User-facing wrappers (e.g., `ez_coverage()`, `ez_gene()`)
- `geom_*`: ggplot2 geom layers
- `process_*_data()`: Data transformation utilities
- `add_*`: Plot annotations (`add_vline()`, `add_hline()`, `add_rect()`, `add_text()`)

### Theme System
- **Base**: `ez_theme()` - Minimal theme removing gridlines/backgrounds
- **Track-specific**: `ez_coverage_theme()`, `ez_gene_theme()`, `ez_feature_theme()`
- **Y-axis control**: `y_axis_style = c("none", "simple", "full")` parameter pattern

### Genomic Region Handling
- Always use `parse_region()` to convert strings to GRanges
- Support flexible separators: "chr1:1000-2000", "chr1_1000-2000", "chr1:1000:2000"
- Apply `scale_x_genome_region()` to ensure consistent coordinate scales

### Example Data Pattern
- Raw files in `inst/extdata/` (GTF, BED, bedgraph, etc.)
- Processed `.rda` objects in `data/` (loaded via `data()`)
- Creation script in `data-raw/create_example_data.R`

## Common Pitfalls

1. **Coordinate mismatch**: Tracks won't align if x-scale limits differ; always apply `scale_x_genome_region()` with the same region
2. **Missing strand levels**: `geom_gene()` expects factor levels `c("minus", "plus", "Unknown")` with labels `c("-", "+", "Unknown")`
3. **GRanges vs data.frame**: Use `methods::is(obj, "GRanges")` to type-check; convert with `granges_to_df()`
4. **File imports**: Always pass `which = region_gr` to `rtracklayer::import()` for large files
5. **Namespace conflicts**: Use explicit `ggplot2::` prefixes in ggproto methods (not standard functions)

## Key Files Reference
- **Core API**: `R/ez_wrappers.R`, `R/genome_plot.R`
- **Geoms**: `R/geom_gene.R`, `R/geom_coverage.R`, `R/geom_*.R`
- **Utilities**: `R/helpers.R`, `R/scales.R`, `R/theme.R`
- **Package metadata**: `DESCRIPTION`, `NAMESPACE` (auto-generated)
- **Vignettes**: `vignettes/introduction.Rmd`, `vignettes/advanced_usage.Rmd`
