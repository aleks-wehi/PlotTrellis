# Every SVG PlotTrellis writes is tidied for vector editors
# (.pt_tidy_svg(): no invisible rects, no textLength width lock, a sans-serif
# fallback, xml:space='preserve'), and the junction/Hi-C legends are vector rects,
# not raster images. Needs the mouse fixture.

svg_lines <- function(items, region, settings) {
  f <- tempfile(fileext = ".svg")
  on.exit(unlink(f), add = TRUE)
  utils::capture.output(render_all(items, region, settings,
                                   device_fn = function(w, h) svglite::svglite(f, width = w, height = h)))
  grDevices::dev.off()
  readLines(f)
}
rect_styles <- function(s) ifelse(grepl("^<rect .*style='", s), sub("^<rect .*style='([^']*)'.*$", "\\1", s), NA)

test_that(".pt_tidy_svg() removes only invisible rects and fixes text and fonts (mouse dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  rows <- test_rows("mouse")
  raw <- svg_lines(list(single_track("bw", rows[["Ba"]]), single_track("transcripts", rows[["Ba"]])),
                   test_region("mouse"), test_settings("mouse"))
  f <- tempfile(fileext = ".svg")
  on.exit(unlink(f), add = TRUE)
  writeLines(raw, f)
  .pt_tidy_svg(f)
  .pt_tidy_svg(f)  # idempotent
  s <- readLines(f)

  invisible <- function(x) { st <- rect_styles(x); !is.na(st) & grepl("stroke: none", st) & !grepl("fill:", st) }
  expect_true(any(grepl("^<rect width='100%'", raw)) && any(invisible(raw)))  # the raw file has them
  expect_false(any(grepl("^<rect width='100%'", s)))
  expect_false(any(invisible(s)))
  # Everything else is kept line for line: visible rects, clip-path rects, all marks and text.
  expect_identical(sum(grepl("^<rect", s)), sum(grepl("^<rect", raw)) - sum(invisible(raw)) - 1L)
  expect_identical(sum(grepl("<text", s)), sum(grepl("<text", raw)))
  expect_identical(sum(grepl("<polygon", s)), sum(grepl("<polygon", raw)))
  expect_true(any(grepl("textLength=", raw)))
  expect_false(any(grepl("textLength=|lengthAdjust=", s)))
  fonts <- unlist(regmatches(s, gregexpr("font-family: [^;']*", s)))
  expect_gt(length(fonts), 0)
  expect_true(all(fonts == 'font-family: "DejaVu Sans", sans-serif'))
  expect_identical(sum(grepl("xml:space='preserve'", s, fixed = TRUE)), 1L)
  expect_true(grepl("^<svg xml:space='preserve' ", s[grepl("^<svg ", s)]))
})

test_that("pt_save() writes a tidied SVG", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  local_mocked_bindings(render_all = function(items, region_params, settings, device_fn) {
    device_fn(2, 1)
    grid::grid.rect(gp = grid::gpar(col = NA, fill = NA))
    grid::grid.text("label", gp = grid::gpar(fontfamily = "DejaVu Sans"))
  })
  fig <- structure(list(entries = list(), region_params = list(), settings = list()), class = "pt_figure")
  f <- tempfile(fileext = ".svg")
  on.exit(unlink(f), add = TRUE)
  pt_save(fig, f)
  s <- readLines(f)
  expect_true(any(grepl("xml:space='preserve'", s, fixed = TRUE)))
  expect_false(any(grepl("^<rect width='100%'|textLength=", s)))
})

test_that("junction opacity and Hi-C colour legends are vector, not raster (mouse dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  rows <- test_rows("mouse")
  region <- test_region("mouse")
  settings <- test_settings("mouse")
  jmd <- svg_lines(list(pair_track("junctions_mirrored_diff", rows[["Ba"]], rows[["P6"]])), region, settings)
  hic <- svg_lines(list(single_track("hic", rows[["Ba"]])), region, settings)
  for (s in list(jmd, hic)) {
    expect_false(any(grepl("<image", s, fixed = TRUE)))
    expect_gte(sum(grepl("^<rect .*fill: #", s)), 100)  # the 100-step colour bar
  }
  expect_true(any(grepl(">Opacity<", jmd, fixed = TRUE)))
})

test_that(".pt_tidy_svg() edits only opening tags, can keep the width lock, and skips a missing file", {
  f <- tempfile(fileext = ".svg")
  on.exit(unlink(f), add = TRUE)
  writeLines(c("<svg xmlns='http://www.w3.org/2000/svg' width='10pt' height='10pt'>",
               "<text x='1' y='2' style='font-family: \"DejaVu Sans\";' textLength='5px' lengthAdjust='spacingAndGlyphs'>a font-family: Evil; textLength='1' lengthAdjust='x' b</text>",
               "</svg>"), f)
  g <- tempfile(fileext = ".svg")
  on.exit(unlink(g), add = TRUE)
  file.copy(f, g)
  .pt_tidy_svg(f)
  s <- readLines(f)
  expect_identical(s[2], paste0("<text x='1' y='2' style='font-family: \"DejaVu Sans\", sans-serif;'>",
                                "a font-family: Evil; textLength='1' lengthAdjust='x' b</text>"))
  .pt_tidy_svg(g, unlock_text = FALSE)
  expect_true(grepl("textLength='5px' lengthAdjust='spacingAndGlyphs'>", readLines(g)[2], fixed = TRUE))
  expect_silent(.pt_tidy_svg(tempfile(fileext = ".svg")))
})

test_that("pt_save() closes only its own device and leaves the caller's current one current", {
  local_mocked_bindings(render_all = function(items, region_params, settings, device_fn) {
    device_fn(2, 1)
    grid::grid.text("label")
  })
  fig <- structure(list(entries = list(), region_params = list(), settings = list()), class = "pt_figure")
  pdfs <- replicate(3, tempfile(fileext = ".pdf"))
  devs <- vapply(pdfs, function(p) { grDevices::pdf(p); grDevices::dev.cur() }, integer(1))
  on.exit({ for (d in intersect(devs, grDevices::dev.list())) grDevices::dev.off(d); unlink(pdfs) }, add = TRUE)
  grDevices::dev.set(devs[2])
  f <- tempfile(fileext = ".svg")
  on.exit(unlink(f), add = TRUE)
  pt_save(fig, f)
  expect_identical(unname(grDevices::dev.list()[grDevices::dev.list() %in% devs]), unname(devs))
  expect_identical(unname(grDevices::dev.cur()), unname(devs[2]))
})
