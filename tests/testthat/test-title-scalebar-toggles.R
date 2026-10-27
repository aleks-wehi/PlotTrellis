# Feature 40 (Phase I.3 step 28): show_title / show_scale_bar switch off the
# figure title and the first coverage track's scale bar. Both default to on,
# so a settings list without them (test_settings(), or a script calling
# render_all() directly) must give today's figure.
test_that("show_title and show_scale_bar remove the title and the scale bar", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  rows          <- test_rows("mouse")
  region_params <- test_region("mouse")
  spec          <- single_track("bw", rows[["Ba"]])
  base          <- test_settings("mouse")
  base$title    <- "Toggle test title"

  svg_of <- function(settings) {
    f <- tempfile(fileext = ".svg")
    on.exit(unlink(f), add = TRUE)
    render_all(list(spec), region_params, settings,
               device_fn = function(w, h) svglite::svglite(f, width = w, height = h))
    grDevices::dev.off()
    paste(readLines(f), collapse = "\n")
  }
  page_h <- function(svg) as.numeric(sub(".*height='([0-9.]+)pt'.*", "\\1", svg))
  scale_label <- format_bp_generic(nice_bp(region_params$end - region_params$start))

  on_svg <- svg_of(base)
  expect_identical(on_svg, svg_of(modifyList(base, list(show_title = TRUE, show_scale_bar = TRUE))))
  expect_match(on_svg, "Toggle test title", fixed = TRUE)
  expect_match(on_svg, paste0(">", scale_label, "<"), fixed = TRUE)

  no_title <- svg_of(modifyList(base, list(show_title = FALSE)))
  expect_no_match(no_title, "Toggle test title", fixed = TRUE)
  expect_match(no_title, paste0(">", scale_label, "<"), fixed = TRUE)
  expect_equal(page_h(on_svg) - page_h(no_title), 0.4 * 72, tolerance = 0.5)

  no_bar <- svg_of(modifyList(base, list(show_scale_bar = FALSE)))
  expect_match(no_bar, "Toggle test title", fixed = TRUE)
  expect_no_match(no_bar, paste0(">", scale_label, "<"), fixed = TRUE)
})
