# Phase I.1 (2026-10-27): render_all() isolates each track's errors. A failing
# settings_/render_ call becomes a red in-figure message naming the track and
# the reason, and the other tracks still render. Asserts on the SVG text
# rather than golden images, so no new snapshot files.

render_svg_text <- function(items, region, settings) {
  tmp <- tempfile(fileext = ".svg")
  on.exit(unlink(tmp), add = TRUE)
  render_all(items, region, settings,
             device_fn = function(w, h) svglite::svglite(tmp, width = w, height = h))
  grDevices::dev.off()
  paste(readLines(tmp, warn = FALSE), collapse = "\n")
}

test_that("a settings_ error fails only its own track, with the reason shown", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  rows <- test_rows("mouse")
  local_mocked_bindings(settings_sashimi = function(...) stop("boom in settings"))
  items <- list(single_track("bw", rows[["Ba"]]), single_track("sashimi", rows[["Ba"]]))
  expect_message(svg <- render_svg_text(items, test_region("mouse"), test_settings("mouse")),
                 "sashimi \\(.*\\) track failed: boom in settings")
  expect_match(svg, "track failed: boom in settings", fixed = TRUE)
  expect_match(svg, "<polygon", fixed = TRUE)  # the bw signal track still drew
})

test_that("a render_ error fails only its own track, with the reason shown", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  rows <- test_rows("mouse")
  local_mocked_bindings(render_bw = function(...) stop("boom in render"))
  items <- list(single_track("bw", rows[["Ba"]]), single_track("sashimi", rows[["Ba"]]))
  expect_message(svg <- render_svg_text(items, test_region("mouse"), test_settings("mouse")),
                 "bw \\(.*\\) track failed: boom in render")
  expect_match(svg, "track failed: boom in render", fixed = TRUE)
  expect_match(svg, ">No data<", fixed = TRUE)  # the sashimi track (no junctions here) still drew
})

test_that("the error label survives a message with invalid bytes, and covers partial drawing", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  rows <- test_rows("mouse")
  local_mocked_bindings(render_bw = function(entry, dict, ctx) {
    plotText("PARTIAL", x = ctx$track_x, y = ctx$y, just = c("left", "top"))
    stop("bad \xff byte")
  })
  svg <- suppressMessages(render_svg_text(list(single_track("bw", rows[["Ba"]])), test_region("mouse"), test_settings("mouse")))
  expect_match(svg, "track failed: bad ? byte", fixed = TRUE)
  # a white blanking rect is drawn after (so on top of) the partial text
  partial_at <- regexpr(">PARTIAL<", svg, fixed = TRUE)
  expect_gt(partial_at, 0)
  expect_true(any(gregexpr("fill: #FFFFFF", svg, fixed = TRUE)[[1]] > partial_at))
})

test_that("sashimi heading accepts a non-integer junction score (BUGS.md #20)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  rows   <- test_rows("mouse")
  region <- test_region("mouse")
  mid    <- round((region$start + region$end) / 2)
  # A FitHiChIP-style BEDPE score, fed in as the pooled junction data.
  local_mocked_bindings(load_junctions = function(...)
    data.frame(intron_start = mid - 5000L, intron_end = mid + 5000L, nSupport = 3.457920544))
  svg <- render_svg_text(list(single_track("sashimi", rows[["Ba"]])), region, test_settings("mouse"))
  expect_match(svg, "max_h=3.46 reads", fixed = TRUE)
  expect_false(grepl("track failed", svg, fixed = TRUE))
})
