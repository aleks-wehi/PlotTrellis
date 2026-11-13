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
#' @importFrom sortable rank_list bucket_list add_rank_list update_rank_list
#' @importFrom tools R_user_dir
#' @importFrom BiocManager install
#' @importFrom R.utils isGzipped
#' @importFrom grDevices dev.cur dev.off
#' @importFrom methods as
#' @importFrom stats aggregate setNames
#' @importFrom utils head modifyList read.delim tail
"_PACKAGE"

# Where the gene-level caches (protein domains, UniProt canonical isoforms, the
# transcripts shown per gene) live: the open project's cache/ folder, carried as
# settings$cache_dir, or the per-user cache folder when there is no project
# (cache_dir NULL). Computed per call: the downloaded R script gets the project's
# folder from its settings, and the per-user folder of the machine that wrote it is
# never baked in.
.pt_cache_root <- function(cache_dir) cache_dir %||% tools::R_user_dir("PlotTrellis", which = "cache")

# Writes one cache file, through a temporary file renamed into place so a reader never
# sees half a file. A cache only saves time: one that can't be written (a folder that
# can't be created, or is read-only) is one warning naming the file, never an error, so
# the caller still returns its result.
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
    warning("PlotTrellis: could not write the cache file ", path, "; carrying on without caching this result.",
            call. = FALSE)
  }
  invisible(isTRUE(ok))
}
