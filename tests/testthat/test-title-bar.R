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
