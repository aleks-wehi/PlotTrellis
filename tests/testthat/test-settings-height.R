# Phase D.3 step 4a: non-default settings$atac_height/hic_height/hic_style,
# asserting the rendered track height actually reflects the override (not
# just that a render happens to succeed either way).
#
# render_all() computes entry$height internally, during its settings-fold
# loop (pipeline/_render_all.R calls each item's settings_<viz>()) -- there
# is no pre-render height on a caller-built entry to assert on directly. So
# each case below calls the real settings_bw()/settings_hic() function
# itself (expected_settings_height(), below) to get the computed height, the
# same call render_all() makes internally, then separately builds/renders
# the actual figure to get a doppelganger baseline too.
#
# Each case asserts the NUMERIC height differs from default first (the
# meaningful-difference signal, survives any future spacing/layout
# regeneration) and adds one doppelganger baseline to catch a visual
# regression too.

expected_settings_height <- function(viz, entry, settings, region_params, page_w = 9) {
  track_w <- page_w - 0.8 - 1.0
  fn <- switch(viz, bw = settings_bw, hic = settings_hic)
  fn(list(), entry, settings, track_w, region_params)$entry$height
}

# app.R:build_entries()'s hic_height() formula, reproduced here since app.R
# itself (Shiny reactive code) isn't sourced by this harness -- pipeline/viz/
# hic.R's settings_hic() uses this exact same formula as its own default.
expected_hic_height <- function(page_w = 9, hic_style = "rectangle") {
  track_w <- page_w - 0.8 - 1.0
  if (hic_style == "square") track_w else round(track_w / 2, 2)
}

test_that("settings$atac_height override changes bw track height (mouse dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  rows          <- test_rows("mouse")
  region_params <- test_region("mouse")

  settings_default <- test_settings("mouse")
  entry_default <- single_track("bw", rows[["Ba"]])
  height_default <- expected_settings_height("bw", entry_default, settings_default, region_params)
  expect_identical(height_default, 0.85)

  settings_override <- test_settings("mouse", overrides = list(atac_height = 1.6))
  entry_override <- single_track("bw", rows[["Ba"]])
  height_override <- expected_settings_height("bw", entry_override, settings_override, region_params)
  expect_identical(height_override, 1.6)
  expect_false(identical(height_default, height_override))

  fig <- test_fig(list(entry_override), region_params, settings_override,
                          perf_id = "atac-height-override-mouse", module = "bw")
  vdiffr::expect_doppelganger("atac-height-override-mouse", fig)
})

test_that("settings$hic_height override changes hic track height (human dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  rows          <- test_rows("human")
  region_params <- test_region("human", flank_pct = 25)

  settings_default <- test_settings("human")
  entry_default <- single_track("hic", rows[["GM12878"]], settings = settings_default)
  height_default <- expected_settings_height("hic", entry_default, settings_default, region_params)
  expect_identical(height_default, expected_hic_height(9, "rectangle"))

  settings_override <- test_settings("human", overrides = list(hic_height = 3.0))
  entry_override <- single_track("hic", rows[["GM12878"]], settings = settings_override)
  height_override <- expected_settings_height("hic", entry_override, settings_override, region_params)
  expect_identical(height_override, 3.0)
  expect_false(identical(height_default, height_override))

  fig <- test_fig(list(entry_override), region_params, settings_override,
                          perf_id = "hic-height-override-human", module = "hic")
  vdiffr::expect_doppelganger("hic-height-override-human", fig)
})

test_that("settings$hic_style square vs rectangle changes hic track height (mouse dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  rows          <- test_rows("mouse")
  region_params <- test_region("mouse")

  settings_rect   <- test_settings("mouse")
  settings_square <- test_settings("mouse", overrides = list(hic_style = "square"))
  entry_rect   <- single_track("hic", rows[["Ba"]], settings = settings_rect)
  entry_square <- single_track("hic", rows[["Ba"]], settings = settings_square)

  height_rect   <- expected_settings_height("hic", entry_rect, settings_rect, region_params)
  height_square <- expected_settings_height("hic", entry_square, settings_square, region_params)

  expect_false(identical(height_rect, height_square))
  expect_identical(height_square, expected_hic_height(9, "square"))
  expect_identical(entry_square$hic_style, "square")

  # No doppelganger here: BUGS.md #23 (plotHicSquare()'s width==height check
  # tripping on a ~9e-16 float mismatch between two independently-recomputed
  # copies of track_w) was a property of the OLD codegen, which recomputed
  # track_w once as a literal embedded in the assembled script and again at
  # render_code() call time. pipeline/_render_all.R computes track_w ONCE and
  # threads the same value through both settings_hic() and ctx$track_w, so
  # the original root cause may no longer apply -- not reverified here
  # (out of this migration's scope); the numeric assertions above are what
  # this test guards regardless.
})
