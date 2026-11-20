# The wizard's Save (.pt_wizard_save_server(), R/wizard-save.R) on its own: plain reactives for
# the fields, a temp project folder, no genome. The "saved" dialog is caught, not shown, and
# the delayed stopApp() is caught, not run.

new_save_project <- function() {
  proj <- tempfile("pt-save-")
  dir.create(file.path(proj, "data", "atac"), recursive = TRUE)
  file.create(file.path(proj, "data", "atac", "a.bw"))
  proj
}

ba_sample <- list(list(condition = "Ba", colour = "#1B9E77", data = list(list(type = "bw", name = "ATAC", path = "atac/a.bw"))))

# A server made of nothing but the Save block for the project in `proj`. The genome fields
# are inputs (Save reads them there); the samples, folder and collected settings are plain
# reactives.
save_app <- function(proj, samples = ba_sample) {
  function(input, output, session) {
    data_root   <- reactiveVal(file.path(proj, "data"))
    wiz_samples <- reactiveVal(samples)
    save_failed <- reactiveVal(FALSE)
    write_error <- reactiveVal(NULL)
    flat <- function() c(list(gtf_path = wiz_clean_optional_text(input$wiz_gtf_path),
                              orgdb = wiz_clean_optional_text(input$wiz_orgdb),
                              genome_label = wiz_clean_optional_text(input$wiz_genome_label)),
                         .pt_collect_display_settings(list(), "", .pt_default_settings()))
    res <- .pt_wizard_save_server(input, output, file.path(proj, "project_settings.json"),
                                  file.path(proj, "project_samples.txt"), data_root, wiz_samples, flat,
                                  save_failed, write_error)
  }
}

# Catches the dialogs shown and the stopApp() calls scheduled, for the rest of the calling test.
catch_saved <- function(env = parent.frame()) {
  got <- new.env(); got$shown <- character(0); got$stops <- character(0)
  local_mocked_bindings(showModal = function(ui, ...) got$shown <- c(got$shown, as.character(ui)), .env = env)
  local_mocked_bindings(later = function(func, ...) {
    got$stops <- c(got$stops, grep("stopApp", deparse(body(func)), value = TRUE))
    invisible(NULL)
  }, .package = "later", .env = env)
  got
}

shown_error <- function(out) if (is.null(out)) "" else paste(as.character(out$html), collapse = "")

test_that("a Save with problems writes nothing; its list then follows the fields", {
  proj <- new_save_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  got <- catch_saved()
  shiny::testServer(save_app(proj), {
    session$setInputs(wiz_gtf_path = "", wiz_orgdb = "org.Mm.eg.db", wiz_genome_label = "mm10")
    expect_identical(shown_error(output$wiz_completion_error), "")  # nothing before a Save
    session$setInputs(wiz_save = 1)
    expect_true(save_failed())
    html <- shown_error(output$wiz_completion_error)
    expect_match(html, "Not saved yet. Please fix:", fixed = TRUE)
    expect_match(html, "GTF path is required.", fixed = TRUE)
    expect_identical(res$save_problems(), "GTF path is required.")
    wiz_samples(list())
    session$flushReact()
    expect_match(shown_error(output$wiz_completion_error), "Add at least one sample before saving.", fixed = TRUE)
    wiz_samples(ba_sample)
    session$setInputs(wiz_gtf_path = "/x/genes.gtf.gz")
    expect_identical(shown_error(output$wiz_completion_error), "")
  })
  expect_false(file.exists(file.path(proj, "project_settings.json")))
  expect_length(got$shown, 0)
  expect_length(got$stops, 0)
})

test_that("a Save writes the settings and the samples copy, then restarts into the app", {
  proj <- new_save_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  got <- catch_saved()
  expect_message(shiny::testServer(save_app(proj), {
    session$setInputs(wiz_gtf_path = "/x/genes.gtf.gz", wiz_orgdb = "org.Mm.eg.db", wiz_genome_label = "mm10")
    session$setInputs(wiz_save = 1)
    expect_false(save_failed())
    expect_null(write_error())
  }), "settings saved", fixed = TRUE)
  json <- jsonlite::fromJSON(file.path(proj, "project_settings.json"), simplifyDataFrame = FALSE)
  expect_identical(json$data_dir, file.path(proj, "data"))
  expect_identical(json$samples[[1]]$condition, "Ba")
  expect_identical(json$gtf_path, "/x/genes.gtf.gz")
  expect_identical(json$genome_label, "mm10")
  expect_null(json$start_view)
  expect_true(file.exists(file.path(proj, "project_samples.txt")))
  expect_match(got$shown[[1]], "Settings saved", fixed = TRUE)
  expect_match(got$shown[[1]], file.path(proj, "project_samples.txt"), fixed = TRUE)
  expect_identical(trimws(got$stops), 'stopApp(returnValue = "restart")')
})

test_that("a re-Save keeps the project's samples.txt source", {
  proj <- new_save_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  writeLines('{"samples_txt_source": "/data/samples.txt"}', file.path(proj, "project_settings.json"))
  catch_saved()
  suppressMessages(shiny::testServer(save_app(proj), {
    session$setInputs(wiz_gtf_path = "/x/genes.gtf.gz", wiz_orgdb = "org.Mm.eg.db", wiz_genome_label = "mm10")
    session$setInputs(wiz_save = 1)
  }))
  json <- jsonlite::fromJSON(file.path(proj, "project_settings.json"), simplifyDataFrame = FALSE)
  expect_identical(json$samples_txt_source, "/data/samples.txt")
  expect_identical(json$samples[[1]]$condition, "Ba")
})

test_that("a folder that can't be written gives a message, cleared by the next Save", {
  proj <- new_save_project()
  on.exit({ Sys.chmod(proj, "0755"); unlink(proj, recursive = TRUE) }, add = TRUE)
  Sys.chmod(proj, "0555")
  skip_if(file.access(proj, 2) == 0, "can't make a read-only folder here (running as root?)")
  got <- catch_saved()
  shiny::testServer(save_app(proj), {
    session$setInputs(wiz_gtf_path = "/x/genes.gtf.gz", wiz_orgdb = "org.Mm.eg.db", wiz_genome_label = "mm10")
    session$setInputs(wiz_save = 1)
    expect_match(write_error(), "Couldn't write .*project_settings.json.*Permission denied")
    expect_match(shown_error(output$wiz_completion_error), "Not saved:", fixed = TRUE)
    # The folder fixed, a field emptied: the next Save lists that field, not the old write error.
    Sys.chmod(proj, "0755")
    session$setInputs(wiz_gtf_path = "")
    session$setInputs(wiz_save = 2)
    expect_null(write_error())
    expect_match(shown_error(output$wiz_completion_error), "Not saved yet. Please fix", fixed = TRUE)
    expect_no_match(shown_error(output$wiz_completion_error), "Permission denied", fixed = TRUE)
    session$setInputs(wiz_gtf_path = "/x/genes.gtf.gz")
    suppressMessages(session$setInputs(wiz_save = 3))
    expect_null(write_error())
    expect_identical(shown_error(output$wiz_completion_error), "")
  })
  expect_length(got$shown, 1)
  expect_true(file.exists(file.path(proj, "project_settings.json")))
})
