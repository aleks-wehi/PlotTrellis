# Phase I.3 step 28 (stage B): per-setting effect tests for the Hi-C settings
# (hic_norm, zrange_cap, hic_style triangle vs rectangle). Each test renders
# or loads the same track at two values of one setting and asserts a specific
# difference. Needs the mouse fixture (Ba .hic around Elf5, a ~45Kb window at
# 10Kb resolution, so each render is well under a second).
#
# zrange_cap reaches render_all() as settings$zrange = c(0, cap) (the
# composer renames it, R/app.R settings_base); that rename is the composer's
# wiring and is not covered here.

hic_fixture <- function() {
  list(rows = test_rows("mouse"), region = test_region("mouse"), settings = test_settings("mouse"))
}

hic_svg <- function(entry, region, settings) {
  f <- tempfile(fileext = ".svg")
  on.exit(unlink(f), add = TRUE)
  dev <- NULL
  on.exit(if (!is.null(dev) && dev %in% grDevices::dev.list()) grDevices::dev.off(dev), add = TRUE, after = FALSE)
  utils::capture.output(render_all(list(entry), region, settings, device_fn = function(w, h) {
    svglite::svglite(f, width = w, height = h)
    dev <<- grDevices::dev.cur()
  }))
  grDevices::dev.off(dev)
  dev <- NULL
  paste(readLines(f), collapse = "\n")
}

svg_labels <- function(svg) {
  sub(".*>([^<]*)</text>$", "\\1", regmatches(svg, gregexpr("<text[^>]*>[^<]*</text>", svg))[[1]])
}

# plotgardener draws each Hi-C pixel as one <polygon>; returns their fills and x coordinates.
hic_pixels <- function(svg) {
  polys <- regmatches(svg, gregexpr("<polygon[^>]*>", svg))[[1]]
  fills <- ifelse(grepl("fill: #", polys, fixed = TRUE), sub(".*fill: (#[0-9A-Fa-f]+).*", "\\1", polys), NA_character_)
  pts <- sub(".*points='([^']*)'.*", "\\1", polys)
  xs <- lapply(strsplit(trimws(pts), " "), function(p) as.numeric(sub(",.*", "", p)))
  list(fill = fills, min_x = vapply(xs, min, numeric(1)), max_x = vapply(xs, max, numeric(1)))
}

test_that("hic_norm changes the loaded contact matrix (mouse dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  fx <- hic_fixture()
  entry <- single_track("hic", fx$rows[["Ba"]], settings = fx$settings)
  kr   <- load_hic(entry, fx$region, fx$settings)$data$matrix
  none <- load_hic(entry, fx$region, modifyList(fx$settings, list(hic_norm = "NONE")))$data$matrix
  expect_gt(nrow(kr), 0)
  expect_identical(kr[, c("x", "y")], none[, c("x", "y")])
  # NONE is raw read counts (integers); KR-balanced values are not.
  expect_true(all(none$counts == round(none$counts)))
  expect_false(isTRUE(all.equal(kr$counts, none$counts)))
})

test_that("hic_norm changes the drawn heatmap at a fixed colour scale (mouse dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  fx <- hic_fixture()
  entry <- single_track("hic", fx$rows[["Ba"]], settings = fx$settings)
  # Fixed zrange so only plotgardener's own re-read (ctx$hic_norm) can change the pixels.
  fixed <- modifyList(fx$settings, list(zrange = c(0, 20)))
  kr   <- hic_pixels(hic_svg(entry, fx$region, fixed))$fill
  none <- hic_pixels(hic_svg(entry, fx$region, modifyList(fixed, list(hic_norm = "NONE"))))$fill
  expect_gt(sum(!is.na(kr)), 0)
  expect_false(identical(kr, none))
  # Auto colour scale: the legend maximum follows the normalised matrix (KR max 48, raw max 74).
  auto_kr   <- svg_labels(hic_svg(entry, fx$region, fx$settings))
  auto_none <- svg_labels(hic_svg(entry, fx$region, modifyList(fx$settings, list(hic_norm = "NONE"))))
  expect_true("48" %in% auto_kr)
  expect_true("74" %in% auto_none)
})

test_that("zrange_cap (settings$zrange) sets the heatmap colour-scale maximum (mouse dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  fx <- hic_fixture()
  entry <- single_track("hic", fx$rows[["Ba"]], settings = fx$settings)
  auto   <- hic_svg(entry, fx$region, fx$settings)
  capped <- hic_svg(entry, fx$region, modifyList(fx$settings, list(zrange = c(0, 7))))
  expect_true("48" %in% svg_labels(auto))   # data max, today's default
  expect_false("7" %in% svg_labels(auto))
  expect_true("7" %in% svg_labels(capped))
  expect_false("48" %in% svg_labels(capped))
  # More pixels saturate at the cap colour.
  expect_false(identical(hic_pixels(auto)$fill, hic_pixels(capped)$fill))
})

test_that("hic_style triangle draws a different heatmap shape from rectangle (mouse dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  fx <- hic_fixture()
  rect <- single_track("hic", fx$rows[["Ba"]], settings = fx$settings)
  tri  <- single_track("hic", fx$rows[["Ba"]], settings = modifyList(fx$settings, list(hic_style = "triangle")))
  expect_identical(tri$hic_style, "triangle")
  track_lo <- 0.8 * 72
  track_hi <- (0.8 + fx$settings$page_w - 1.8) * 72
  px_rect <- hic_pixels(hic_svg(rect, fx$region, fx$settings))
  px_tri  <- hic_pixels(hic_svg(tri,  fx$region, fx$settings))
  expect_gt(length(px_tri$fill), 0)
  # A rectangle fills its corners with pixels from outside the window (drawn past
  # the track edges, then clipped); a triangle only draws in-window contacts.
  expect_true(any(px_rect$min_x < track_lo - 1))
  expect_true(all(px_tri$min_x >= track_lo - 1 & px_tri$max_x <= track_hi + 1))
})
