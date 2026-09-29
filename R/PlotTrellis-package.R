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
#' @importFrom grDevices dev.cur dev.off
#' @importFrom methods as
#' @importFrom stats aggregate setNames
#' @importFrom utils head modifyList read.delim tail
"_PACKAGE"
