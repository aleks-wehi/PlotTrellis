# After Clear all caches: an open project's figure says how to rebuild the genome files, and
# the app's dialog says when files could not be deleted.

test_that("a figure whose genome files were deleted says to re-run pt_project()", {
  expect_error(.pt_annotation(list(txdb_path = file.path(tempfile(), "genome.sqlite"))),
               "genome files are gone .*Run pt_project\\(\\) on the project again")
})

test_that("Clear all caches says what it could not delete", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  proj <- normalizePath(tempfile("pt-clear-"), mustWork = FALSE)
  dir.create(proj)
  stuck <- file.path(proj, "cache", "stuck")
  on.exit({ Sys.chmod(stuck, "0755"); unlink(proj, recursive = TRUE) }, add = TRUE)
  seed_project_genome(proj)
  json <- file.path(proj, "project_settings.json")
  jsonlite::write_json(list(data_dir = FIXTURES$mouse$data_dir, samples = FIXTURES$mouse$samples,
                            gtf_path = GENOMES$mouse$gtf_path, orgdb = GENOMES$mouse$orgdb,
                            genome_label = GENOMES$mouse$genome_label), json, auto_unbox = TRUE)
  app <- suppressMessages(.pt_composer_app(json, .pt_read_settings(json, .pt_default_settings())))
  dir.create(stuck, recursive = TRUE)
  writeLines("x", file.path(stuck, "held.rds"))
  Sys.chmod(stuck, "0555")
  skip_if(file.access(stuck, 2) == 0, "can't make a read-only folder here (running as root?)")
  shown <- list()
  local_mocked_bindings(showModal = function(ui, ...) shown[[length(shown) + 1]] <<- as.character(ui))
  local_mocked_bindings(later = function(...) invisible(NULL), .package = "later")
  shiny::testServer(app, {
    session$flushReact()
    session$setInputs(open_settings = 1)
    session$setInputs(confirm_clear_caches = 1)
  })
  dialog <- shown[[length(shown)]]
  expect_match(dialog, "Caches partly cleared", fixed = TRUE)
  expect_match(dialog, "stuck/held.rds", fixed = TRUE)
  expect_false(grepl("All caches have been deleted", dialog, fixed = TRUE))
})

test_that("closing the open genome databases lets their files go", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  # Meaningful where tempdir() is on a network drive (as on the HPC this was found on).
  dir <- tempfile("pt-close-")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  db <- file.path(dir, "genome.sqlite")
  file.copy(GENOMES$mouse$txdb_path, db)
  .pt_annotation(list(txdb_path = db))
  expect_true(db %in% ls(.pt_txdb_memo))
  .pt_close_txdbs()
  expect_false(db %in% ls(.pt_txdb_memo))
  unlink(db)
  expect_length(list.files(dir, all.files = TRUE, no.. = TRUE), 0)  # no .nfs stand-in left
})
