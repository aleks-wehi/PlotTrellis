# Phase I.3 step 24 (BUGS.md #51): render_all()'s fallback diff colours must
# match the app's (#CC0000 / #0033CC, as in DEFAULT_SETTINGS and
# test_settings()), so a caller whose settings leave them out (e.g. a script
# calling render_all() directly) gets the same figure as the app.
# Compares against test_settings()'s copy of the values (the app's own are
# in .pt_default_settings(), R/app.R).
test_that("diff colours fall back to the app's defaults when unset", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  rows          <- test_rows("mouse")
  region_params <- test_region("mouse")
  spec          <- pair_track("bw_diff_stacked", rows[["Ba"]], rows[["LP"]])
  explicit      <- test_settings("mouse")
  unset         <- explicit
  unset$diff_colour_hi <- NULL
  unset$diff_colour_lo <- NULL

  svg_of <- function(settings) {
    f <- tempfile(fileext = ".svg")
    svglite::svglite(f, width = 9, height = 4)
    on.exit(unlink(f), add = TRUE)
    dev <- grDevices::dev.cur()
    on.exit(if (dev %in% grDevices::dev.list()) grDevices::dev.off(dev), add = TRUE, after = FALSE)
    render_all(list(spec), region_params, settings)
    grDevices::dev.off(dev)
    paste(readLines(f), collapse = "\n")
  }
  fallback <- svg_of(unset)
  # Guards against a vacuous pass: the fixture must actually draw both colours.
  expect_match(fallback, "#CC0000", fixed = TRUE)
  expect_match(fallback, "#0033CC", fixed = TRUE)
  expect_identical(fallback, svg_of(explicit))
})
