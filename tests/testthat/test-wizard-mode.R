# The wizard's landing screen and mode switch (.pt_wizard_mode_server(), R/wizard-mode.R) on
# their own: plain reactiveVals, no genome, no project folder.

# A server made of nothing but the mode block, starting in mode `start`. A failed Save is
# pending (save_failed TRUE, a write error) so Back's clearing of it shows.
mode_app <- function(start = NULL, is_edit = FALSE, has_mirror = FALSE, mirror_txt = "") {
  function(input, output, session) {
    wiz_mode    <- reactiveVal(start)
    save_failed <- reactiveVal(TRUE)
    write_error <- reactiveVal("could not write")
    .pt_wizard_mode_server(input, output, wiz_mode, save_failed, write_error, .pt_default_settings(),
                           is_edit, has_mirror, mirror_txt)
  }
}

body_html <- function(out) paste(as.character(out$html), collapse = "")

test_that("the landing screen offers the two ways in and no form", {
  shiny::testServer(mode_app(), {
    html <- body_html(output$wiz_body)
    expect_match(html, 'id="wiz_mode_manual"', fixed = TRUE)
    expect_match(html, 'id="wiz_mode_import"', fixed = TRUE)
    expect_match(html, "What a project needs", fixed = TRUE)
    expect_no_match(html, 'id="wiz_gtf_path"', fixed = TRUE)
    expect_no_match(html, 'id="wiz_mode_back"', fixed = TRUE)
  })
})

test_that("Set up by hand shows the form with a Back link and no import box", {
  shiny::testServer(mode_app(), {
    session$setInputs(wiz_mode_manual = 1)
    expect_identical(wiz_mode(), "manual")
    html <- body_html(output$wiz_body)
    expect_match(html, 'id="wiz_gtf_path"', fixed = TRUE)
    expect_match(html, 'id="wiz_mode_back"', fixed = TRUE)
    expect_match(html, "Save and open the app", fixed = TRUE)
    expect_no_match(html, 'id="wiz_import_path"', fixed = TRUE)
    expect_no_match(html, 'id="wiz_close"', fixed = TRUE)
  })
})

test_that("Import shows the import box, empty, above the form", {
  shiny::testServer(mode_app(), {
    session$setInputs(wiz_mode_import = 1)
    expect_identical(wiz_mode(), "import")
    html <- body_html(output$wiz_body)
    expect_match(html, 'id="wiz_import_path"[^>]*value=""')
    expect_match(html, "The file fills in the three sections below", fixed = TRUE)
    expect_match(html, 'id="wiz_gtf_path"', fixed = TRUE)
  })
})

test_that("a project_samples.txt found at start prefills the import box with its path", {
  mirror <- file.path(tempdir(), "proj", "project_samples.txt")
  shiny::testServer(mode_app("import", has_mirror = TRUE, mirror_txt = mirror), {
    html <- body_html(output$wiz_body)
    expect_match(html, sprintf('id="wiz_import_path"[^>]*value="%s"', mirror))
    expect_match(html, "Click Import to restore it", fixed = TRUE)
  })
})

test_that("Back returns to the landing screen and clears a failed Save's messages", {
  shiny::testServer(mode_app(), {
    session$flushReact()  # Back ignores its initial value
    session$setInputs(wiz_mode_manual = 1)
    expect_true(save_failed())
    session$setInputs(wiz_mode_back = 1)
    expect_null(wiz_mode())
    expect_false(save_failed())
    expect_null(write_error())
    expect_match(body_html(output$wiz_body), 'id="wiz_mode_manual"', fixed = TRUE)
  })
})

test_that("edit mode shows the form with Save and return and Close, and no Back link", {
  shiny::testServer(mode_app("edit", is_edit = TRUE), {
    html <- body_html(output$wiz_body)
    expect_match(html, "Save and return to the app", fixed = TRUE)
    expect_match(html, 'id="wiz_close"', fixed = TRUE)
    expect_no_match(html, 'id="wiz_mode_back"', fixed = TRUE)
    expect_no_match(html, 'id="wiz_import_path"', fixed = TRUE)
  })
})
