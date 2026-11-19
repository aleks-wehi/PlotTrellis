# The bundled example project (plottrellis_example(),
# inst/extdata/example/). Uses only the installed example data, never the
# data/ fixtures, so these run on R-universe too.

# helper-gz-snapshots.R only unpacks the gzipped golden images when the data/
# fixtures exist; this file's one has to be there without them, so unpack it here.
local({
  gz <- test_path("_snaps_gz", "example-data", "example-eya2.svg.gz")
  svg <- test_path("_snaps", "example-data", "example-eya2.svg")
  if (file.exists(gz)) {
    dir.create(dirname(svg), recursive = TRUE, showWarnings = FALSE)
    con <- gzfile(gz, "rt")
    writeLines(readLines(con, warn = FALSE), svg)
    close(con)
  }
})

# Every file of the installed example with its md5, to show nothing writes there.
example_fingerprint <- function() {
  src <- system.file("extdata", "example", package = "PlotTrellis")
  files <- list.files(src, recursive = TRUE, all.files = TRUE)
  setNames(unname(tools::md5sum(file.path(src, files))), files)
}

test_that("plottrellis_example() gives the installed example, not a copy in tempdir()", {
  a <- plottrellis_example()
  expect_identical(a, pt_example("mammary"))
  expect_same_path(dirname(a), system.file("extdata", "example", package = "PlotTrellis"))
  expect_false(startsWith(normalizePath(a), normalizePath(tempdir())))
})

test_that("the example samples.txt imports with no errors or warnings", {
  txt <- plottrellis_example()
  dir <- normalizePath(dirname(txt), winslash = "/")
  res <- import_samples_txt(txt, defaults = .pt_default_settings())
  expect_length(res$errors, 0)
  expect_length(res$warnings, 0)

  skip_if_not_installed("org.Mm.eg.db")
  project <- tempfile("pt-example-project-")
  expect_no_warning(settings <- parse_samples_txt(txt, project_dir = project))
  expect_true(file.exists(file.path(project, "project_settings.json")))
  expect_same_path(settings$data_dir, dir)
  expect_same_path(settings$gtf_path, file.path(dir, "annotation.gtf.gz"))
  expect_identical(settings$orgdb, "org.Mm.eg.db")
  expect_identical(settings$genome_label, "mm39")
  expect_identical(vapply(settings$samples, `[[`, "", "condition"), c("Basal", "LP"))
  paths <- function(s) c(sprintf(c("rna/%s_RNA.bw", "rna/%s_junctions.bed", "rna/%s_salmon", "atac/%s_ATAC.bw",
                                   "atac/%s_ATAC_peaks.narrowPeak", "cuttag/%s_H3K36me3.bw",
                                   "cuttag/%s_H3K36me3_peaks.broadPeak", "cuttag/%s_IgG.bw", "hic/%s_HiC.mcool"), s))
  expect_identical(lapply(settings$samples, function(s) vapply(s$data, `[[`, "", "path")),
                   list(paths("Basal"), paths("LP")))
})

test_that("the example's Eya2 figure renders where it is installed, offline, writing only to the project; golden image", {
  skip_if_not_installed("org.Mm.eg.db")
  before <- example_fingerprint()
  txt <- plottrellis_example()
  project <- tempfile("pt-example-project-")
  # The shipped chromosome sizes must be used, not fetched from UCSC.
  local_mocked_bindings(build_chrom_sizes = function(...) stop("tried to download chromosome sizes"))
  # txdbmaker's two notices for any GENCODE GTF are expected; anything else still warns.
  proj <- withCallingHandlers(pt_project(project, txt), warning = function(w) {
    if (grepl("\"phase\" metadata column|genome version information", conditionMessage(w)))
      invokeRestart("muffleWarning")
  })
  cache_dir <- file.path(normalizePath(project), "cache")
  # Every genome file is in the project's cache, the chromosome sizes copied from the shipped ones.
  genome <- proj$genome
  for (f in unlist(genome[c("txdb_path", "tx_name_lookup_path", "tx_meta_lookup_path", "chrom_sizes_path")]))
    expect_true(startsWith(normalizePath(f), normalizePath(file.path(cache_dir, "genome"))), info = f)
  expect_identical(readRDS(genome$chrom_sizes_path),
                   readRDS(file.path(dirname(txt), "annotation_chrom_sizes.rds")))

  # Eya2's promoter switch: the transcripts, ATAC and transcript expression of both samples.
  fig <- pt_figure(proj, "Eya2", c("transcripts", "ATAC", "transcript_expression"))
  draw <- function() render_all(fig$entries, fig$region_params, fig$settings)

  # Rendered and checked outside vdiffr first: expect_doppelganger() skips when
  # NOT_CRAN is unset, and these checks must run everywhere.
  svg <- tempfile(fileext = ".svg")
  svglite::svglite(svg)
  draw()
  grDevices::dev.off()
  text <- paste(readLines(svg, warn = FALSE), collapse = "\n")
  for (label in c("Eya2-201", "Eya2-202", "Basal", "LP", "Transcript expression"))
    expect_match(text, label, fixed = TRUE)
  expect_identical(example_fingerprint(), before)

  skip_if_not_installed("vdiffr")
  vdiffr::expect_doppelganger("example-eya2", draw)
})
