# A hand-edited project_settings.json with a blank genome label is refused when the project is
# opened (pt_project() and the composer app), with the same message the wizard's Save gives,
# instead of failing at the first launch of a new GTF with a console error. Uses the bundled
# example, so needs no fixture data.
test_that("opening a project with a blank genome label stops with 'Genome label is required'", {
  proj <- tempfile("pt-blank-label-")
  on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  # Setup only: building the example's TxDb warns about its GTF, unrelated to this test.
  suppressWarnings(suppressMessages(pt_project(proj, pt_example())))
  json <- file.path(proj, "project_settings.json")
  settings <- jsonlite::fromJSON(json, simplifyDataFrame = FALSE)
  expect_identical(settings$genome_label, "mm10")  # the label the edits below blank out

  for (blank in list("", "  ", NULL)) {
    edited <- settings
    edited["genome_label"] <- list(blank)  # NULL is kept as JSON null
    writeLines(jsonlite::toJSON(edited, auto_unbox = TRUE, null = "null", digits = NA), json)
    expect_error(pt_project(proj), "Genome label is required: the UCSC name of the genome build", fixed = TRUE)
    expect_error(.pt_composer_app(json, .pt_read_settings(json, .pt_default_settings())),
                 "Genome label is required", fixed = TRUE)
  }
  # The message says where to fix it. Matched on the project folder's name, not the full path:
  # the message prints the normalised path, which differs from tempfile()'s where the temp dir
  # is a symlink (macOS /var -> /private/var).
  expect_error(pt_project(proj), file.path(basename(proj), "project_settings.json"), fixed = TRUE)

  # The unedited project still opens.
  writeLines(jsonlite::toJSON(settings, auto_unbox = TRUE, null = "null", digits = NA), json)
  expect_no_error(suppressWarnings(suppressMessages(pt_project(proj))))
})
