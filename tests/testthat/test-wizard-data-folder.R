# The wizard's Data folder field (.pt_wizard_data_folder_server(), R/wizard-data-folder.R) on
# its own: the folder it returns and the status line under the field.

data_folder_app <- function(input, output, session) {
  res <- .pt_wizard_data_folder_server(input, output)
}

status_html <- function(out) if (is.null(out)) "" else paste(as.character(out$html), collapse = "")

test_that("a blank or all-space folder is NULL and shows no status line", {
  shiny::testServer(data_folder_app, {
    session$setInputs(wiz_data_root = "")
    expect_null(res$data_root())
    expect_identical(status_html(output$wiz_data_root_status), "")
    session$setInputs(wiz_data_root = "   ")
    expect_null(res$data_root())
  })
})

test_that("an existing folder is returned trimmed and says it was found", {
  dir <- tempfile("pt-data-"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  shiny::testServer(data_folder_app, {
    session$setInputs(wiz_data_root = paste0("  ", dir, " "))
    expect_identical(res$data_root(), dir)
    expect_match(status_html(output$wiz_data_root_status), "wiz-ok.*Folder found")
  })
})

test_that("a folder that doesn't exist is still returned, with an error line", {
  dir <- file.path(tempdir(), "pt-no-such-folder")
  shiny::testServer(data_folder_app, {
    session$setInputs(wiz_data_root = dir)
    expect_identical(res$data_root(), dir)
    expect_match(status_html(output$wiz_data_root_status), "wiz-error.*No folder at this path")
  })
})
