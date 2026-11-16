# Bug 191: after a Save that failed, the wizard's "Not saved yet" list re-checks as the fields
# change: a fixed problem leaves the list, and the list goes once everything passes, without
# another Save. Drives the real wizard_server() headlessly (first-time wizard, no settings yet).

test_that("the Not saved yet list follows the fields after a failed Save", {
  proj <- tempfile("pt-saveerr-")
  dir.create(file.path(proj, "data", "atac"), recursive = TRUE)
  file.create(file.path(proj, "data", "atac", "a.bw"))
  on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  txt <- tempfile(fileext = ".txt")
  on.exit(unlink(txt), add = TRUE)
  writeLines(c("condition\tcolour\tbw:ATAC", "Ba\t#1B9E77\tatac/a.bw"), txt)
  json <- file.path(proj, "project_settings.json")
  defaults <- .pt_default_settings()
  shown <- function(out) if (is.null(out)) "" else paste(out$html, collapse = "")

  eval(bquote(shiny::testServer(wizard_server(.(json), .pt_read_settings(.(json), .(defaults)), .(defaults)), {
    session$setInputs(wiz_mode_import = 1)
    session$setInputs(wiz_data_root = "", wiz_gtf_path = "", wiz_orgdb = "org.Mm.eg.db", wiz_genome_label = "mm10")
    expect_identical(.(shown)(output$wiz_completion_error), "")  # nothing before a Save
    session$setInputs(wiz_save = 1)
    html <- .(shown)(output$wiz_completion_error)
    expect_match(html, "Not saved yet", fixed = TRUE)
    expect_match(html, "Data root is required.", fixed = TRUE)
    expect_match(html, "Add at least one sample before saving.", fixed = TRUE)
    expect_match(html, "GTF path is required.", fixed = TRUE)

    # Import fills the samples; the folder is typed: those two leave the list, the GTF stays.
    session$setInputs(wiz_import_path = .(txt), wiz_import_load = 1)
    session$setInputs(wiz_data_root = .(file.path(proj, "data")))
    html <- .(shown)(output$wiz_completion_error)
    expect_match(html, "GTF path is required.", fixed = TRUE)
    expect_no_match(html, "Data root is required.", fixed = TRUE)
    expect_no_match(html, "Add at least one sample", fixed = TRUE)

    # The last one fixed: the list goes, with no new Save, and nothing was written.
    session$setInputs(wiz_gtf_path = "/x/genes.gtf.gz")
    expect_identical(.(shown)(output$wiz_completion_error), "")
    expect_false(file.exists(.(json)))
  })))
})

test_that("Save re-checks what changed outside the fields, and Back clears the list", {
  proj <- tempfile("pt-saveerr-")
  dir.create(proj)
  on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  txt <- tempfile(fileext = ".txt")
  on.exit(unlink(txt), add = TRUE)
  writeLines(c("condition\tcolour\tbw:ATAC", "Ba\t#1B9E77\tatac/a.bw"), txt)
  json <- file.path(proj, "project_settings.json")
  data <- file.path(proj, "data")
  defaults <- .pt_default_settings()
  shown <- function(out) if (is.null(out)) "" else paste(out$html, collapse = "")

  eval(bquote(shiny::testServer(wizard_server(.(json), .pt_read_settings(.(json), .(defaults)), .(defaults)), {
    session$setInputs(wiz_mode_import = 1)
    session$setInputs(wiz_import_path = .(txt), wiz_import_load = 1)
    session$setInputs(wiz_data_root = .(data), wiz_gtf_path = "/x/genes.gtf.gz",
                      wiz_orgdb = "org.Mm.eg.db", wiz_genome_label = "mm10")
    session$setInputs(wiz_save = 1)  # the data folder doesn't exist yet
    expect_match(.(shown)(output$wiz_completion_error), "Data root directory does not exist.", fixed = TRUE)
    expect_false(file.exists(.(json)))

    # Back to the start, then Import again: no list before a new Save.
    session$setInputs(wiz_mode_back = 1)
    session$setInputs(wiz_mode_import = 2)
    expect_identical(.(shown)(output$wiz_completion_error), "")

    # The folder made outside the app (no field changed): the next Save goes through.
    dir.create(.(data))
    session$setInputs(wiz_save = 2)
    expect_true(file.exists(.(json)))
    expect_identical(.(shown)(output$wiz_completion_error), "")
  })))
})
