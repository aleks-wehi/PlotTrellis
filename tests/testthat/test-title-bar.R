# W91: the setup wizard has the composer's title bar, from the same builder.
test_that("the wizard page has the composer's title bar, with its stylesheet in the head", {
  page <- lapply(htmltools::renderTags(wizard_ui(file.path(tempdir(), "no-such-settings.json")))[c("head", "html")],
                 as.character)
  expect_match(page[["html"]], '<nav class="navbar navbar-default navbar-static-top"', fixed = TRUE)
  expect_match(page[["html"]], 'class="plottrellis-version"', fixed = TRUE)
  expect_match(page[["html"]], paste0("v", utils::packageVersion("PlotTrellis")), fixed = TRUE)
  expect_match(page[["head"]], ".navbar.navbar-default", fixed = TRUE)  # the theme, in the head
  expect_false(grepl("\\.navbar\\.navbar-default\\s*\\{", page[["html"]]))
  # Full width: the bar comes before the page's (width-capped) container, not inside it.
  expect_lt(regexpr("<nav class=", page[["html"]], fixed = TRUE),
            regexpr('class="container-fluid wiz-page"', page[["html"]], fixed = TRUE))
  # The composer's navbarPage title is the same builder: checked in the source, which only exists
  # in a source tree (not under R CMD check or an installed package).
  app_src <- testthat::test_path("..", "..", "R", "app-ui.R")
  skip_if_not(file.exists(app_src), "package source not available")
  expect_true(any(grepl("title = .pt_brand()", readLines(app_src, warn = FALSE), fixed = TRUE)))
})

# ---- added since ---------------------------------------------------------------------------

# W91: the setup wizard's title bar.
test_that("the wizard page has its title bar, with its stylesheet in the head", {
  page <- lapply(htmltools::renderTags(wizard_ui(file.path(tempdir(), "no-such-settings.json")))[c("head", "html")],
                 as.character)
  expect_match(page[["html"]], '<nav class="navbar navbar-default navbar-static-top"', fixed = TRUE)
  expect_match(page[["html"]], 'class="plottrellis-version"', fixed = TRUE)
  expect_match(page[["html"]], paste0("v", utils::packageVersion("PlotTrellis")), fixed = TRUE)
  expect_match(page[["head"]], ".navbar.navbar-default", fixed = TRUE)  # the theme, in the head
  expect_false(grepl("\\.navbar\\.navbar-default\\s*\\{", page[["html"]]))
  # Full width: the bar comes before the page's (width-capped) container, not inside it.
  expect_lt(regexpr("<nav class=", page[["html"]], fixed = TRUE),
            regexpr('class="container-fluid wiz-page"', page[["html"]], fixed = TRUE))
})
