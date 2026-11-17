#' @details
#' For developers: adding a new track type is a matter of defining its
#' `info_<key>()`, `settings_<key>()` and `render_<key>()` functions -- see
#' [viz-module-contract] for the convention and [viz-module-discovery] for how
#' modules are found.
#'
#' @keywords internal
#' @import shiny
#' @import plotgardener
#' @import GenomicRanges
#' @import IRanges
#' @importFrom GenomicFeatures cdsBy exonsBy transcripts genes
#' @importFrom AnnotationDbi loadDb select keys
#' @importFrom rtracklayer import
#' @importFrom GenomeInfoDb seqlevels seqlengths
#' @importFrom S4Vectors metadata
#' @importFrom jsonlite fromJSON toJSON
#' @importFrom svglite svglite
#' @importFrom grDevices dev.cur dev.off
#' @importFrom methods as
#' @importFrom stats aggregate setNames
#' @importFrom utils head modifyList read.delim tail
"_PACKAGE"

# Where every cache lives (genome files, protein domains, UniProt canonical
# isoforms, the transcripts shown per gene): the open project's cache/ folder,
# carried as settings$cache_dir. There is no per-user fallback, so deleting the
# project folder removes everything PlotTrellis wrote.
.pt_cache_root <- function(cache_dir) {
  if (is.null(cache_dir) || !nzchar(cache_dir))
    stop(
      "PlotTrellis keeps its caches in the project folder, and none was given ",
      "(open a project with pt_project() or run_plottrellis()).",
      call. = FALSE
    )
  cache_dir
}

# Writes one cache file, through a temporary file renamed into place so a reader
# never sees half a file. A cache only saves time: one that can't be written (a
# folder that can't be created, or is read-only) is one warning naming the file,
# never an error, so the caller still returns its result.
.pt_cache_save <- function(value, path) {
  tmp <- NULL
  ok <- tryCatch({
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    tmp <- tempfile(tmpdir = dirname(path))
    suppressWarnings(saveRDS(value, tmp))
    suppressWarnings(file.rename(tmp, path))
  }, error = function(e) FALSE)
  if (!isTRUE(ok)) {
    if (!is.null(tmp)) unlink(tmp)
    warning("PlotTrellis: could not write the cache file ", path,
            "; carrying on without caching this result.",
            call. = FALSE)
  }
  invisible(isTRUE(ok))
}
