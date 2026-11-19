# The bundled example project. Its data and README are in
# inst/extdata/example/, rebuilt by data-raw/build-example-data.R.

#' Example project: a ready-made samples.txt and its data
#'
#' The path of PlotTrellis's small example project's `samples.txt`, ready for
#' [parse_samples_txt()]; the same file as `pt_example("mammary")`. The data is
#' mouse mammary basal cells and luminal progenitors (Milevskiy et al.
#' 2023, *Cell Genomics*, GEO GSE227750) at five genes, Eya2, Lrrfip1, Acot7,
#' Dzip1 and Cxcl14 (mm39, GENCODE vM39), with every data type: RNA-seq
#' coverage, splice junctions and transcript expression, ATAC and H3K36me3
#' coverage and peaks, IgG and Hi-C. Search for Eya2 in the composer: basal
#' cells mostly use its isoform from an upstream promoter, luminal progenitors a
#' shorter one.
#'
#' The example needs the mouse OrgDb package, `org.Mm.eg.db`, which PlotTrellis
#' doesn't install for you: `BiocManager::install("org.Mm.eg.db")`. For your own
#' data, install your species' package instead (e.g. `org.Hs.eg.db` for human).
#'
#' Its data is read where it is installed: a project keeps everything it writes
#' (settings, caches) in its own folder. Open the returned `samples.txt` to see
#' how a project is described; a copy of it is a good starting point for your
#' own.
#'
#' @return The path of the example's `samples.txt`.
#' @seealso [parse_samples_txt()], [run_plottrellis()]
#' @examples
#' samples_txt <- plottrellis_example()
#' writeLines(readLines(samples_txt))
#' # Needs org.Mm.eg.db: BiocManager::install("org.Mm.eg.db")
#' if (interactive()) {
#'   project <- tempfile()
#'   parse_samples_txt(samples_txt, project_dir = project)
#'   run_plottrellis(project_dir = project)
#' }
#' @export
plottrellis_example <- function() pt_example("mammary")
