# A hand-edited project_settings.json with a display setting of the wrong type or length:
# each falls back to its default with one warning, instead of stopping the figure.

test_that("unusable display settings fall back to their defaults with one warning", {
  json <- tempfile(fileext = ".json")
  writeLines('{"page_w": "wide", "fontsize": [], "show_title": "yes", "atac_height": [1, 2],
               "zrange_cap": "abc", "hic_height": {}, "gene_flank_pct": 20, "font": "Arial"}', json)
  expect_warning(cfg <- .pt_read_settings(json, .pt_default_settings()),
                 "for fontsize, page_w, show_title, zrange_cap, atac_height; using the defaults")
  d <- .pt_display_defaults()
  expect_identical(cfg[c("page_w", "fontsize", "show_title", "atac_height")], d[c("page_w", "fontsize", "show_title", "atac_height")])
  expect_true("zrange_cap" %in% names(cfg))
  expect_null(cfg$zrange_cap)
  expect_equal(cfg$gene_flank_pct, 20)  # good values are kept
  expect_equal(cfg$font, "Arial")
  # A clean file gives no warning, nor do boxes the app saved empty ("NA", null).
  writeLines('{"page_w": 7, "fontsize": 10, "zrange_cap": null, "atac_height": "NA", "gene_flank_pct": null}', json)
  expect_no_warning(cfg <- .pt_read_settings(json, .pt_default_settings()))
  expect_equal(cfg$page_w, 7)
  expect_equal(cfg$atac_height, d$atac_height)
  expect_equal(cfg$gene_flank_pct, d$gene_flank_pct)
})

test_that("a page_w of 0 or below takes the default with the same warning", {
  json <- tempfile(fileext = ".json")
  for (w in c(0, -2)) {
    writeLines(sprintf('{"page_w": %s}', w), json)
    expect_warning(cfg <- .pt_read_settings(json, .pt_default_settings()), "for page_w; using the default", info = paste("page_w", w))
    expect_equal(cfg$page_w, .pt_display_defaults()$page_w, info = paste("page_w", w))
  }
  # A narrow page still draws (its narrow tracks say so themselves): kept.
  writeLines('{"page_w": 2.5}', json)
  expect_no_warning(cfg <- .pt_read_settings(json, .pt_default_settings()))
  expect_equal(cfg$page_w, 2.5)
})

test_that("a figure still draws from a project with a text page_w and an empty fontsize", {
  skip_if_not_installed("org.Hs.eg.db")
  dir <- tempfile("pt-bad-display-")
  suppressWarnings(pt_project(dir, pt_example("human")))
  json <- file.path(dir, "project_settings.json")
  s <- jsonlite::fromJSON(json, simplifyDataFrame = FALSE)
  s$page_w <- "wide"; s$fontsize <- list()
  writeLines(jsonlite::toJSON(s, auto_unbox = TRUE, null = "null", digits = NA), json)
  expect_warning(proj <- pt_project(dir), "for fontsize, page_w; using the defaults")
  fig <- pt_figure(proj, "PSME2", "genes")
  png <- pt_save(fig, tempfile(fileext = ".png"), dpi = 50)
  expect_true(file.exists(png))
  unlink(dir, recursive = TRUE)
})
