# Phase I.2 step 23: the bundled example project (plottrellis_example(),
# inst/extdata/example/). Uses only the installed example data, never the
# data/ fixtures, so these run on R-universe too.

# helper-gz-snapshots.R only unpacks the gzipped golden images when the data/
# fixtures exist; this file's one has to be there without them, so unpack it here.
local({
  gz <- test_path("_snaps_gz", "example-data", "example-elf5.svg.gz")
  svg <- test_path("_snaps", "example-data", "example-elf5.svg")
  if (file.exists(gz)) {
    dir.create(dirname(svg), recursive = TRUE, showWarnings = FALSE)
    con <- gzfile(gz, "rt")
    writeLines(readLines(con, warn = FALSE), svg)
    close(con)
  }
})

# No project open (an earlier test file may have opened one): caches go beside the GTF.
.pt_set_cache_root(NULL)

# Every file of the installed example with its md5, to show nothing writes there.
example_fingerprint <- function() {
  src <- system.file("extdata", "example", package = "PlotTrellis")
  files <- list.files(src, recursive = TRUE, all.files = TRUE)
  setNames(unname(tools::md5sum(file.path(src, files))), files)
}

test_that("plottrellis_example() copies the example to a new folder each call", {
  a <- plottrellis_example()
  b <- plottrellis_example()
  expect_false(identical(dirname(a), dirname(b)))
  expect_true(startsWith(normalizePath(a), normalizePath(tempdir())))
  expect_setequal(list.files(dirname(a), recursive = TRUE), names(example_fingerprint()))
  expect_true(file.access(dirname(a), 2) == 0)
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
  expect_same_path(settings$gtf_path, file.path(dir, "gencode.vM25.example.gtf.gz"))
  expect_identical(settings$orgdb, "org.Mm.eg.db")
  expect_identical(settings$genome_label, "mm10")
  expect_identical(vapply(settings$samples, `[[`, "", "condition"), c("Basal", "LP"))
  expect_identical(
    lapply(settings$samples, function(s) vapply(s$data, `[[`, "", "path")),
    list(c("atac/Basal_ATAC.bw", "quant/Basal_rep1"), c("atac/LP_ATAC.bw", "quant/LP_rep1")))
})

test_that("the example's Elf5 figure renders from a fresh copy, offline; golden image", {
  skip_if_not_installed("org.Mm.eg.db")
  before <- example_fingerprint()
  txt <- plottrellis_example()
  settings <- parse_samples_txt(txt, project_dir = tempfile("pt-example-project-"))
  # The shipped chromosome sizes must be used, not fetched from UCSC.
  local_mocked_bindings(build_chrom_sizes = function(...) stop("tried to download chromosome sizes"))
  # txdbmaker's two notices for any GENCODE GTF are expected; anything else still warns.
  genome <- withCallingHandlers(
    ensure_genome_ready(settings$gtf_path, settings$orgdb, settings$genome_label),
    warning = function(w) {
      if (grepl("\"phase\" metadata column|genome version information", conditionMessage(w)))
        invokeRestart("muffleWarning")
    })
  expect_identical(normalizePath(dirname(genome$txdb_path)), normalizePath(dirname(txt)))

  txdb <- AnnotationDbi::loadDb(genome$txdb_path)
  coords <- gene_to_coords("Elf5", txdb, flank_pct = settings$gene_flank_pct / 100, orgdb_pkg = genome$orgdb)
  region <- list(chrom = coords$chrom, start = coords$start, end = coords$end, gene = "Elf5")
  # The fields the composer's RENDER_CFG adds to the project settings.
  cfg <- modifyList(settings, list(
    txdb_path = genome$txdb_path, orgdb_pkg = genome$orgdb, genome_label = genome$genome_label,
    tx_name_lookup_path = genome$tx_name_lookup_path, tx_meta_lookup_path = genome$tx_meta_lookup_path,
    title = "", output_path = tempfile(fileext = ".svg")))
  rows <- load_samples_from_settings(settings$samples, settings$data_dir)
  names(rows) <- vapply(rows, `[[`, "", "display")
  specs <- list(
    single_track("transcripts", NULL),
    single_track("bw", rows$Basal), single_track("bw", rows$LP),
    single_track("tpm_usage", rows$Basal), single_track("tpm_usage", rows$LP))

  # Rendered and checked outside vdiffr first: expect_doppelganger() skips when
  # NOT_CRAN is unset, and these checks must run everywhere.
  svg <- tempfile(fileext = ".svg")
  svglite::svglite(svg)
  render_all(specs, region, cfg)
  grDevices::dev.off()
  text <- paste(readLines(svg, warn = FALSE), collapse = "\n")
  for (label in c("Elf5-201", "Elf5-207", "Basal", "LP", "Transcript counts"))
    expect_match(text, label, fixed = TRUE)
  expect_identical(example_fingerprint(), before)

  skip_if_not_installed("vdiffr")
  vdiffr::expect_doppelganger("example-elf5", test_fig(specs, region, cfg))
})
