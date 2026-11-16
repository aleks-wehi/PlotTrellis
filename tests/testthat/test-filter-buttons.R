# The composer's track filter offers only the track types the project has: the bundled
# mammary example has no BED, Hi-C or junction tracks, so no buttons for them.
test_that("the composer shows filter buttons only for the project's track types", {
  proj <- tempfile("pt-filters-")
  on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  suppressWarnings(suppressMessages(pt_project(proj, pt_example("mammary"))))
  json <- file.path(proj, "project_settings.json")
  app <- suppressWarnings(suppressMessages(.pt_composer_app(json, .pt_read_settings(json, .pt_default_settings()))))
  ui <- environment(app$httpHandler)$ui
  html <- as.character(htmltools::renderTags(if (is.function(ui)) ui(NULL) else ui)$html)
  btns <- regmatches(html, gregexpr('filter-btn[^>]*data-filter="[^"]+"', html))[[1]]
  filters <- sub('.*data-filter="([^"]+)"', "\\1", btns)
  expect_false(any(c("bed", "hic", "junctions") %in% filters))
  expect_true(all(c("all", "bw", "annotation") %in% filters))
})
