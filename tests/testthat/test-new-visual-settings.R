# Phase I.3 step 28: settings added for testers' asks. Each defaults to the
# value that used to be hard-coded, so a settings list without it gives
# today's figure.

test_that("junction_height sets the junction track's height, 0.85 by default", {
  # No fixture data needed: an unreadable path gives an empty track, which
  # still gets its height.
  entry <- list(viz = "sashimi", name = "Junctions", path = tempfile(fileext = ".bedpe"))
  region_params <- list(chrom = "chr1", start = 1e6, end = 1.01e6)
  height_with <- function(settings) suppressWarnings(settings_sashimi(list(), entry, settings, 7.2, region_params))$entry$height
  expect_identical(height_with(list()), 0.85)
  expect_identical(height_with(list(junction_height = 1.4)), 1.4)
  # An explicit per-track height still wins, as for every other track type.
  entry$height <- 0.5
  expect_identical(height_with(list(junction_height = 1.4)), 0.5)
})

test_that("diff_shared_alpha sets the stacked diff's shared-layer tint, 0.15 by default", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  rows          <- test_rows("mouse")
  region_params <- test_region("mouse")
  spec          <- pair_track("bw_diff_stacked", rows[["Ba"]], rows[["LP"]])
  base          <- test_settings("mouse")

  svg_of <- function(settings) {
    f <- tempfile(fileext = ".svg")
    on.exit(unlink(f), add = TRUE)
    render_all(list(spec), region_params, settings,
               device_fn = function(w, h) svglite::svglite(f, width = w, height = h))
    grDevices::dev.off()
    paste(readLines(f), collapse = "\n")
  }
  # Same blend as render_bw_diff_stacked()'s faint_solid(), towards white.
  tint <- function(alpha) {
    m <- alpha * grDevices::col2rgb(spec$control_colour) / 255 + (1 - alpha)
    grDevices::rgb(m[1], m[2], m[3])
  }
  default_svg <- svg_of(base)
  expect_identical(default_svg, svg_of(modifyList(base, list(diff_shared_alpha = 0.15))))
  expect_match(default_svg, tint(0.15), fixed = TRUE)
  strong <- svg_of(modifyList(base, list(diff_shared_alpha = 0.6)))
  expect_match(strong, tint(0.6), fixed = TRUE)
  expect_no_match(strong, tint(0.15), fixed = TRUE)
})
