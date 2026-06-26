# ezGenomeTracks (development version)

## Major Changes
* BigWig file processing now uses `megadepth` for significantly improved performance
  - Region-based BigWig queries are 2-50x faster depending on the number of regions
  - Particularly beneficial for remote BigWig files (HTTP/HTTPS URLs)
  - Built-in parallelization support for processing multiple BigWig files
  - All existing APIs remain unchanged; upgrade is transparent to users

* v0.0.1 release.
