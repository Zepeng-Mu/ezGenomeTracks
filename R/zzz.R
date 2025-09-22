# Package startup functions

#' @importFrom utils packageVersion
.onAttach <- function(libname, pkgname) {
  # Get package version
  version <- utils::packageVersion("ezGenomeTracks")

  # Create welcome message
  msg <- c(
    paste0("ezGenomeTracks v", version),
    "Easy and flexible genomic track visualization",
    "Use citation('ezGenomeTracks') to see how to cite this package",
    "For documentation and examples, visit: https://github.com/zmu/ezGenomeTracks"
  )

  # Display message when package is attached
  packageStartupMessage(paste(msg, collapse = "\n"))
}

#' @importFrom utils packageVersion
.onLoad <- function(libname, pkgname) {
  # Any package initialization code can go here
  # This function runs before .onAttach
  invisible()
}
