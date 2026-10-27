#' @details
#' For developers: adding a new viz type is a matter of defining a
#' `settings_<key>()`/`render_<key>()` pair -- see [viz-module-contract] for
#' the convention and [viz-module-discovery] for how modules are found.
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
#' @importFrom grDevices dev.cur dev.off
#' @importFrom methods as
#' @importFrom stats aggregate setNames
#' @importFrom utils head modifyList read.delim tail
"_PACKAGE"

# Cache locations (bugs 53, 57). These are top-level constants, not computed
# per call, because the downloaded R script embeds them as literals. So they
# are (re)set here: per user at load time (not at install time, which would
# bake in the build machine's path), then to the project's own cache/ folder
# whenever a project is opened (.pt_project()), so a project keeps all its
# caches in its folder.
.pt_set_cache_root <- function(root = NULL) {
  ns <- environment(.pt_set_cache_root)
  gene_root <- root %||% tools::R_user_dir("PlotTrellis", which = "cache")
  vals <- list(
    GENOME_CACHE_DIR = if (is.null(root)) NULL else file.path(root, "genome"),
    PROTEIN_DOMAINS_CACHE_DIR = file.path(gene_root, "protein_domains"),
    UNIPROT_CANONICAL_CACHE_PATH = file.path(gene_root, "uniprot_canonical_cache.rds"),
    GENE_DISPLAY_TRANSCRIPTS_CACHE_PATH = file.path(gene_root, "gene_display_transcripts_cache.rds")
  )
  for (n in names(vals)) {
    locked <- bindingIsLocked(n, ns)
    if (locked) unlockBinding(n, ns)
    assign(n, vals[[n]], envir = ns)
    if (locked) lockBinding(n, ns)
  }
}

.onLoad <- function(libname, pkgname) .pt_set_cache_root(NULL)
