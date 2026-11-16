# Building the composer's UI must not warn that a navbar child isn't a nav panel: the page's
# <head> tags are passed as the navbar's header, not as loose children (bslib warns once per
# loose child). Uses the bundled example, so needs no fixture data.
test_that("the composer UI is built without bslib's 'Navigation containers' warning", {
  proj <- tempfile("pt-nav-")
  on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  # Setup only: building the example's TxDb warns about its GTF, unrelated to the UI.
  suppressWarnings(suppressMessages(pt_project(proj, pt_example())))
  json <- file.path(proj, "project_settings.json")
  warns <- character(0)
  app <- withCallingHandlers(
    suppressMessages(.pt_composer_app(json, .pt_read_settings(json, .pt_default_settings()))),
    warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
  expect_false(any(grepl("Navigation containers", warns, fixed = TRUE)))
  # The head content moved into `header` must still end up in <head>.
  ui <- environment(app$httpHandler)$ui
  head_html <- htmltools::renderTags(if (is.function(ui)) ui(NULL) else ui)$head
  expect_true(grepl("TRACK_TYPE_MAP", head_html, fixed = TRUE))
  expect_true(grepl("--bio-font-base", head_html, fixed = TRUE))
})
