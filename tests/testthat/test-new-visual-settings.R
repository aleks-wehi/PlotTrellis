# Display settings testers asked for. Each defaults to the value that was
# hard-coded before the setting existed, so a settings list without it gives
# the same figure.

test_that("junction_height sets the height of a junction track's arcs, 0.85 by default, under a 0.3 in heading band", {
  # No fixture data needed: an empty file gives an empty track, which still gets its height.
  junction_file <- tempfile(fileext = ".bedpe")
  file.create(junction_file)
  on.exit(unlink(junction_file), add = TRUE)
  entry <- list(viz = "junctions", name = "Junctions", path = junction_file)
  region_params <- list(chrom = "chr1", start = 1e6, end = 1.01e6)
  # Called directly, a module's settings function gets the settings render_all() would give it.
  height_with <- function(settings)
    suppressWarnings(settings_junctions(list(), entry, .pt_fill_settings(settings), 7.2, region_params))$entry$height
  expect_identical(height_with(list()), 0.3 + 0.85)
  expect_identical(height_with(list(junction_height = 1.4)), 0.3 + 1.4)
  # A mirrored comparison is the heading band plus exactly two junction_heights.
  mirrored <- function(settings, entry = list(viz = "junctions_mirrored_diff"))
    settings_junctions_mirrored_diff(list(), entry, .pt_fill_settings(settings), 7.2, region_params)$entry$height
  expect_identical(mirrored(list()), 0.3 + 2 * 0.85)
  expect_identical(mirrored(list(junction_height = 1.4)), 0.3 + 2 * 1.4)
  # An explicit per-track height still wins (it then covers the heading band and the arcs).
  entry$height <- 0.5
  expect_identical(height_with(list(junction_height = 1.4)), 0.5)
  expect_identical(mirrored(list(junction_height = 1.4), list(viz = "junctions_mirrored_diff", height = 0.9)), 0.9)
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
