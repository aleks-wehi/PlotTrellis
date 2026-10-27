# The bundled example project (phase-I.2 step 23). Its data and README are in
# inst/extdata/example/, rebuilt by data-raw/build-example-data.R.

#' Example project: a ready-made samples.txt and its data
#'
#' Copies PlotTrellis's small example project into a new folder and returns the path of its
#' `samples.txt`, ready for [parse_samples_txt()]. The data is mouse mammary basal cells and
#' luminal progenitors (Milevskiy et al. 2023, *Cell Genomics*, GEO GSE227750): ATAC-seq
#' coverage and Salmon transcript quantification for each, over the Ehf, Elf5 and Cat genes
#' (mm10, GENCODE vM25). Search for Elf5 in the composer: it is almost silent in basal cells
#' and highly expressed in luminal progenitors.
#'
#' The example needs the mouse OrgDb package, `org.Mm.eg.db`, which PlotTrellis doesn't
#' install for you: `BiocManager::install("org.Mm.eg.db")`. For your own data, install your
#' species' package instead (e.g. `org.Hs.eg.db` for human).
#'
#' Each call makes a fresh copy (in a new folder under [tempdir()], so it is deleted when R
#' exits) because PlotTrellis writes genome caches next to the GTF on first use, and the
#' installed package folder may not be writable. Open the returned `samples.txt` to see how a
#' project is described; it is a good starting point for your own.
#'
#' @return The path of the copied `samples.txt`.
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
plottrellis_example <- function() {
  src <- system.file("extdata", "example", package = "PlotTrellis", mustWork = TRUE)
  dest <- tempfile("plottrellis-example-")
  dir.create(dest)
  # copy.mode = FALSE: an installed package's files may be read-only, and the copy is
  # meant to be edited and to receive genome caches.
  ok <- file.copy(list.files(src, full.names = TRUE), dest, recursive = TRUE, copy.mode = FALSE)
  if (!all(ok)) stop("Couldn't copy the example data to '", dest, "'.", call. = FALSE)
  file.path(dest, "samples.txt")
}
