# ezGenomeTracks 0.0.15

## Changes

* Consolidated GTF/GFF loading in `ez_gene()`:
  - GTF/GFF file paths are now standardized through `import_gtf()` to produce a consistent data-frame schema.
  - Retained direct `rtracklayer::import()` handling only for non-GTF generic file paths.
  - Updated `ez_gene()` documentation to reflect the canonical `import_gtf()` pathway and deprecate `auto_import` behavior.

# ezGenomeTracks (development version)

## Major Changes
* BigWig file processing now uses `megadepth` for significantly improved performance
  - Region-based BigWig queries are 2-50x faster depending on the number of regions
  - Particularly beneficial for remote BigWig files (HTTP/HTTPS URLs)
  - Built-in parallelization support for processing multiple BigWig files
  - All existing APIs remain unchanged; upgrade is transparent to users

* v0.0.1 release.
