test_that("smoke: single bw track renders (mouse dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  rows     <- test_rows("mouse")
  region   <- test_region("mouse")
  settings <- test_settings("mouse")
  spec     <- single_track("bw", rows[["Ba"]])
  fig      <- test_fig(list(spec), region, settings,
                               perf_id = "smoke-bw-mouse-ba", module = "bw")
  vdiffr::expect_doppelganger("smoke-bw-mouse-ba", fig)
})
