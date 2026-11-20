# Project isolation (W28): everything PlotTrellis writes is inside the project folder, so
# deleting it removes it all. Nothing beside the GTF or in data folders, no per-user cache.

# A one-gene GTF in a folder of its own.
isolation_gtf <- function(dir = tempfile("pt-gtf-")) {
  dir.create(dir, showWarnings = FALSE)
  path <- file.path(dir, "mini.gtf")
  g <- 'gene_id "ENSMUSG00000000001.1"; gene_name "Mini";'
  t <- paste(g, 'transcript_id "ENSMUST00000000001.1"; transcript_name "Mini-201"; transcript_type "protein_coding";')
  row <- function(type, s, e, attr) paste("chr1", "mini", type, s, e, ".", "+", ".", attr, sep = "\t")
  writeLines(c(row("gene", 2001, 8000, g), row("transcript", 2001, 8000, t),
               row("exon", 2001, 3000, paste(t, 'exon_number 1;')),
               row("exon", 7001, 8000, paste(t, 'exon_number 2;'))), path)
  path
}

ready_in <- function(gtf_path, cache_dir) {
  local_mocked_bindings(build_chrom_sizes = function(genome_label, out_path) saveRDS(c(chr1 = 100000L), out_path))
  suppressMessages(suppressWarnings(
    ensure_genome_ready(gtf_path, "org.Mm.eg.db", "synthetic", cache_dir = cache_dir)))
}

test_that("there is no cache without a project folder", {
  expect_error(.pt_cache_root(NULL), "project folder")
  expect_error(.pt_cache_root(""), "project folder")
  expect_error(.genome_cache_path(isolation_gtf(), "_txdb.sqlite", NULL), "project folder")
})

test_that("genome cache files go to <project>/cache/genome/<hash>/, by GTF path, size and time", {
  gtf <- isolation_gtf()
  cache <- tempfile("pt-cache-")
  p <- .genome_cache_path(gtf, "_txdb.sqlite", cache)
  expect_identical(basename(p), "mini_txdb.sqlite")
  expect_match(basename(dirname(p)), "^[0-9a-f]{12}$")
  expect_same_path(dirname(dirname(p)), file.path(cache, "genome"))
  expect_identical(.genome_cache_path(gtf, "_txdb.sqlite", cache), p)
  # A GTF of the same name elsewhere, or this one changed, gets its own folder.
  expect_false(identical(.genome_cache_path(isolation_gtf(), "_txdb.sqlite", cache), p))
  # R's display options don't change it.
  old <- options(digits = 3, OutDec = ",")
  same <- .genome_cache_path(gtf, "_txdb.sqlite", cache)
  options(old)
  expect_identical(same, p)
  cat("# edited\n", file = gtf, append = TRUE)
  expect_false(identical(.genome_cache_path(gtf, "_txdb.sqlite", cache), p))
})

test_that("building a genome writes only inside the project's cache, and ignores caches beside the GTF", {
  skip_if_not_installed("org.Mm.eg.db")
  skip_if_not_installed("txdbmaker")
  gtf <- isolation_gtf()
  # Left beside the GTF by an older version: not a TxDb, so reading it would fail.
  writeLines("not a TxDb", file.path(dirname(gtf), "mini_txdb.sqlite"))
  before <- tools::md5sum(list.files(dirname(gtf), full.names = TRUE))
  cache <- tempfile("pt-cache-")
  g <- ready_in(gtf, cache)
  expect_identical(tools::md5sum(list.files(dirname(gtf), full.names = TRUE)), before)
  files <- unlist(g[c("txdb_path", "tx_name_lookup_path", "tx_meta_lookup_path", "chrom_sizes_path")])
  expect_true(all(file.exists(files)))
  expect_true(all(startsWith(normalizePath(files), normalizePath(file.path(cache, "genome")))))
  expect_length(AnnotationDbi::keys(AnnotationDbi::loadDb(g$txdb_path), "GENEID"), 1)
})

test_that("Edit project settings shows the project folder and keeps Clear caches and Reset under a closed Advanced", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  proj <- normalizePath(tempfile("pt-isolation-"), mustWork = FALSE)
  dir.create(proj)
  on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  seed_project_genome(proj)
  json <- file.path(proj, "project_settings.json")
  jsonlite::write_json(list(data_dir = FIXTURES$mouse$data_dir, samples = FIXTURES$mouse$samples,
                            gtf_path = GENOMES$mouse$gtf_path, orgdb = GENOMES$mouse$orgdb,
                            genome_label = GENOMES$mouse$genome_label), json, auto_unbox = TRUE)
  app <- suppressMessages(.pt_composer_app(json, .pt_read_settings(json, .pt_default_settings())))

  # The page has Edit project settings (its button sends open_settings); Clear caches and Reset
  # are only in that dialog.
  ui <- environment(app$httpHandler)$ui
  html <- htmltools::renderTags(if (is.function(ui)) ui(NULL) else ui)$html
  expect_match(html, 'id="editBtn"', fixed = TRUE)
  expect_match(html, "setInputValue('open_settings'", fixed = TRUE)
  expect_false(grepl("clear_caches_btn|reset_project_btn", html))

  shown <- list()
  local_mocked_bindings(showModal = function(ui, ...) shown[[length(shown) + 1]] <<- as.character(ui))
  scheduled <- list()
  local_mocked_bindings(later = function(func, ...) {
    scheduled[[length(scheduled) + 1]] <<- func
    invisible(NULL)
  }, .package = "later")
  # Something outside the project's cache folder, in the project: kept.
  writeLines("keep", file.path(proj, "keep.txt"))
  writeLines("cached", file.path(proj, "cache", "gene_display_transcripts_cache.rds"))
  stamp <- function() file.info(list.files(dirname(GENOMES$mouse$gtf_path), full.names = TRUE))[c("size", "mtime")]
  gtf_dir_before <- stamp()
  shiny::testServer(app, {
    session$flushReact()
    session$setInputs(open_settings = 1)
    session$setInputs(confirm_clear_caches = 1)
  })
  dialog <- shown[[1]]
  expect_match(dialog, proj, fixed = TRUE)
  expect_match(dialog, "<details id=\"settings_advanced\">", fixed = TRUE)   # closed: no `open`
  expect_match(dialog, "<summary class=\"btn btn-default btn-sm\">Advanced</summary>", fixed = TRUE)  # a button (bug 192)
  expect_match(dialog, "clear_caches_btn", fixed = TRUE)
  expect_match(dialog, "reset_project_btn", fixed = TRUE)
  # Clear all caches removed the project's cache folder and nothing else.
  expect_false(dir.exists(file.path(proj, "cache")))
  expect_true(file.exists(file.path(proj, "keep.txt")))
  expect_true(file.exists(json))
  expect_identical(stamp(), gtf_dir_before)
  expect_true(all(file.exists(unlist(GENOMES$mouse[c("txdb_path", "tx_name_lookup_path", "tx_meta_lookup_path")]))))
  stops <- unlist(lapply(scheduled, function(f) paste(deparse(body(f)), collapse = " ")))
  expect_identical(grep("stopApp", stops, value = TRUE), "stopApp(returnValue = \"restart\")")
})
