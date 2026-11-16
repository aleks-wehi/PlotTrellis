# Reset project configuration in a folder that can't be written: the settings file stays, so
# the app must say nothing was reset and not restart.

test_that("Reset in a read-only project folder says what could not be deleted and does not restart", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  proj <- tempfile("pt-reset-ro-"); dir.create(proj)
  on.exit({ Sys.chmod(proj, "0755"); unlink(proj, recursive = TRUE) }, add = TRUE)
  seed_project_genome(proj)
  json <- file.path(proj, "project_settings.json")
  writeLines(jsonlite::toJSON(list(data_dir = FIXTURES$mouse$data_dir, samples = FIXTURES$mouse$samples,
                                   gtf_path = GENOMES$mouse$gtf_path, orgdb = GENOMES$mouse$orgdb,
                                   genome_label = GENOMES$mouse$genome_label), auto_unbox = TRUE), json)
  app <- suppressMessages(.pt_composer_app(json, .pt_read_settings(json, .pt_default_settings())))
  Sys.chmod(proj, "0555")
  skip_if(file.access(proj, 2) == 0, "can't make a read-only folder here (running as root?)")
  shown <- list(); restarted <- FALSE
  local_mocked_bindings(showModal = function(ui, ...) shown[[length(shown) + 1]] <<- ui)
  local_mocked_bindings(later = function(...) restarted <<- TRUE, .package = "later")
  # Twice: the observer has ignoreInit = TRUE, and testServer() counts the first value as init.
  shiny::testServer(app, { session$setInputs(confirm_reset_project = 1); session$setInputs(confirm_reset_project = 2) })
  expect_true(file.exists(json))
  expect_false(restarted)
  html <- as.character(shown[[length(shown)]])
  expect_match(html, "Nothing was reset")
  expect_match(html, "Couldn't delete\\s+\\S*project_settings.json")
  expect_match(html, "(Permission denied)", fixed = TRUE)
})
