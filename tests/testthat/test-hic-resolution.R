# Hi-C auto-resolution-tier fixtures, small-span and large-span per
# Hi-C-having dataset (both mouse and human have real per-condition .hic
# files, see helper-fixtures.R's FIXTURES table). Asserts the auto-resolution
# logic (load_hic() in data-loader-hic.R, "auto" hic_res_mode) actually picks the
# expected tier for each span, not just that a render happens to succeed
# either way.
#
# Mechanism: calls load_hic() directly, the
# same resolution-selection logic render_all() invokes internally via
# settings_hic() -- avoids a full render just to read back which tier was
# picked. `avail` is recomputed with the same two lines load_hic() itself
# uses (it doesn't expose the candidate list on the entry).
#
# The mouse case doubles as a check of the mouse .hic files' chrom-naming
# convention against the UCSC-style query chrom.

hic_tier_check <- function(dataset, sample, chrom, start, end) {
  rows          <- test_rows(dataset)
  region_params <- list(chrom = chrom, start = start, end = end, gene = NA_character_)
  settings      <- test_settings(dataset)
  entry         <- single_track("hic", rows[[sample]], settings = settings)
  entry         <- load_hic(entry, region_params, settings)
  is_cool <- tolower(tools::file_ext(entry$path)) %in% c("cool", "mcool")
  avail <- sort(if (is_cool) plotgardener::readCoolBpResolutions(entry$path) else strawr::readHicBpResolutions(entry$path))
  list(res = entry$data$resolution, avail = avail)
}

test_that("hic auto-resolution picks the finest tier for a small (gene-scale) span (mouse dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  region_params <- test_region("mouse")  # Elf5 default, ~45kb
  result <- hic_tier_check("mouse", "Ba", region_params$chrom, region_params$start, region_params$end)
  expect_identical(result$res, min(result$avail))

  rows     <- test_rows("mouse")
  settings <- test_settings("mouse")
  spec     <- single_track("hic", rows[["Ba"]], settings = settings)
  fig      <- test_fig(list(spec), region_params, settings,
                               perf_id = "hic-tier-small-mouse", module = "hic")
  vdiffr::expect_doppelganger("hic-tier-small-mouse", fig)
})

test_that("hic auto-resolution picks a coarser tier for a large (30Mb) span (mouse dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  region_params <- list(chrom = "chr2", start = 80000000, end = 110000000, gene = NA_character_)
  result <- hic_tier_check("mouse", "Ba", region_params$chrom, region_params$start, region_params$end)
  expect_gt(result$res, min(result$avail))
  expect_identical(result$res, 50000)  # checked against the real fixture

  rows     <- test_rows("mouse")
  settings <- test_settings("mouse")
  spec     <- single_track("hic", rows[["Ba"]], settings = settings)
  fig      <- test_fig(list(spec), region_params, settings,
                               perf_id = "hic-tier-large-mouse", module = "hic")
  vdiffr::expect_doppelganger("hic-tier-large-mouse", fig)
})

# Not a true gene-scale (~140kb) span like the mouse small-span case above --
# the human .hic fixture's declared finest resolution (1000bp) has ZERO
# usable KR-normalized signal anywhere in the file (confirmed by direct
# strawr scan, not guessed), so the default gene-scale window always falls
# back to that empty tier and draws an empty track, regardless of which small
# window is chosen.
# A 3Mb window is still an order of magnitude smaller than the 30Mb
# large-span case below and demonstrates real tier differentiation (5000bp
# vs 50000bp) without hitting the broken 1000bp floor -- multiple candidate
# resolutions satisfy HIC_MAX_BINS_PER_VIEW at this span (1000 AND 5000bp
# both do), so "auto" naturally picks the coarser of those two (5000bp)
# instead of ever touching the unusable finest tier.
test_that("hic auto-resolution picks a finer (not the largest available) tier for a moderate (3Mb) span (human dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  region_params <- list(chrom = "chr1", start = 198600000, end = 201600000, gene = NA_character_)
  result <- hic_tier_check("human", "GM12878", region_params$chrom, region_params$start, region_params$end)
  expect_lt(result$res, max(result$avail))
  expect_identical(result$res, 5000)  # checked against the real fixture

  rows     <- test_rows("human")
  settings <- test_settings("human")
  spec     <- single_track("hic", rows[["GM12878"]], settings = settings)
  fig      <- test_fig(list(spec), region_params, settings,
                               perf_id = "hic-tier-small-human", module = "hic")
  vdiffr::expect_doppelganger("hic-tier-small-human", fig)
})

test_that("hic auto-resolution picks a coarser tier for a large (30Mb) span (human dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  region_params <- list(chrom = "chr1", start = 180000000, end = 210000000, gene = NA_character_)
  result <- hic_tier_check("human", "GM12878", region_params$chrom, region_params$start, region_params$end)
  expect_gt(result$res, min(result$avail))
  expect_identical(result$res, 50000)  # checked against the real fixture

  rows     <- test_rows("human")
  settings <- test_settings("human")
  spec     <- single_track("hic", rows[["GM12878"]], settings = settings)
  fig      <- test_fig(list(spec), region_params, settings,
                               perf_id = "hic-tier-large-human", module = "hic")
  vdiffr::expect_doppelganger("hic-tier-large-human", fig)
})

# ---- Fixed bin size (hic_res_mode = "fixed", hic_resolution) ----------------------------------
# An exact bin size, or an error naming the closest ones (never another size drawn instead).

hic_fixed_settings <- function(dataset, size) {
  test_settings(dataset, overrides = list(hic_res_mode = "fixed", hic_resolution = size))
}

test_that("fixed mode uses exactly the asked-for bin size, whatever auto would pick (mouse dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  region <- list(chrom = "chr2", start = 95000000, end = 105000000, gene = NA_character_)   # auto picks 20000 here
  auto <- hic_tier_check("mouse", "Ba", region$chrom, region$start, region$end)
  want <- min(auto$avail)                                                                   # not what auto picks
  expect_false(identical(auto$res, want))
  settings <- hic_fixed_settings("mouse", want)
  entry <- load_hic(single_track("hic", test_rows("mouse")[["Ba"]], settings = settings), region, settings)
  expect_identical(entry$data$resolution, want)
  expect_false(isTRUE(entry$data$failed))
  # Too fine to draw across a 30 Mb view (3,000 bins; plotgardener draws at most 1,000): the finest
  # size that fits, which the colour range is then taken from (bug 137).
  wide <- list(chrom = "chr2", start = 80000000, end = 110000000, gene = NA_character_)
  entry <- load_hic(single_track("hic", test_rows("mouse")[["Ba"]], settings = settings), wide, settings)
  expect_identical(entry$data$resolution, 50000)
})

test_that("fixed mode with a bin size the file lacks stops with the closest sizes the file has (mouse dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  region <- test_region("mouse")
  avail <- hic_tier_check("mouse", "Ba", region$chrom, region$start, region$end)$avail
  want <- 12345
  expect_false(want %in% avail)
  near <- c(10000, 20000, 50000)   # the mouse Ba file's three bin sizes nearest 12,345 (checked against the fixture)
  expect_true(all(near %in% avail))
  settings <- hic_fixed_settings("mouse", want)
  spec <- single_track("hic", test_rows("mouse")[["Ba"]], settings = settings)
  expect_error(load_hic(spec, region, settings),
               paste0("no 12,345 bp bins; closest: ", paste(format(near, big.mark = ",", trim = TRUE), collapse = ", "), " bp"),
               fixed = TRUE)
  # Blank bin size in fixed mode: an error too, not a silent default.
  expect_error(load_hic(spec, region, hic_fixed_settings("mouse", NULL)), "Hi-C Resolution is Fixed but no Bin size (bp) is set.", fixed = TRUE)
})

test_that("a Hi-C track with a missing bin size shows the error in the figure and the other tracks still draw", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  region <- test_region("mouse")
  settings <- hic_fixed_settings("mouse", 12345)
  rows <- test_rows("mouse")
  specs <- list(single_track("hic", rows[["Ba"]], settings = settings), single_track("bw", rows[["Ba"]], settings = settings))
  tmp <- tempfile(fileext = ".svg")
  on.exit(unlink(tmp), add = TRUE)
  suppressMessages(suppressWarnings(
    render_all(specs, region, settings, device_fn = function(w, h) svglite::svglite(tmp, width = w, height = h))))
  grDevices::dev.off()
  svg <- paste(readLines(tmp, warn = FALSE), collapse = "\n")
  expect_match(svg, "track failed: this file has no 12,345 bp bins", fixed = TRUE)
  # Exactly one track failed, and the bigWig track beside it drew.
  expect_false(grepl("No data", svg, fixed = TRUE))
  expect_identical(lengths(regmatches(svg, gregexpr("track failed", svg, fixed = TRUE))), 1L)
})

test_that("the fixed bin size reaches a downloaded script", {
  skip_if_not(PT_TEST_MOUSE_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  settings <- hic_fixed_settings("mouse", 25000)
  spec <- single_track("hic", test_rows("mouse")[["Ba"]], settings = settings)
  script <- assemble_download_script(list(spec), test_region("mouse"), settings)
  expect_match(script, "hic_res_mode = \"fixed\"", fixed = TRUE)
  expect_match(script, "hic_resolution = 25000", fixed = TRUE)
})

test_that("hic_resolution is a NULL-default number that the samples.txt importer reads as a number", {
  expect_null(.pt_default_settings()$hic_resolution)
  expect_identical(clean_optional_numeric(NA), NULL)
  coerced <- .pt_coerce_auto_overrides(list(hic_res_mode = "fixed", hic_resolution = "10000"), .pt_default_settings())
  expect_identical(coerced$overrides$hic_resolution, 10000)
  expect_identical(coerced$overrides$hic_res_mode, "fixed")
})


# A view whose bins all lack a value under the normalisation (the human ChIA-PET fixture has no KR
# values at its finest, 1,000 bp, bins, which auto picks for an 800 kb view) says so in the
# track instead of drawing an empty panel.
test_that("an all-NaN Hi-C view names the normalisation and what to change; NONE draws (human dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  region <- list(chrom = "chr1", start = 198300000, end = 199100000, gene = NA_character_)
  settings <- test_settings("human")
  spec <- single_track("hic", test_rows("human")[["GM12878"]], settings = settings)
  msg <- "no KR values at 1,000 bp bins in view: in Figure settings, Hi-C, set Normalisation (hic_norm) to NONE, or Resolution Fixed, Bin size 5000."
  expect_error(load_hic(spec, region, settings), msg, fixed = TRUE)
  none <- load_hic(spec, region, test_settings("human", list(hic_norm = "NONE")))
  expect_true(any(is.finite(none$data$matrix$counts)))
  # 1,000 bp bins would be over 1,000 across this view, which plotgardener can't draw: the next
  # size, 5,000 bp, which has values (no error).
  wide <- list(chrom = "chr1", start = 198300000, end = 199500000, gene = NA_character_)
  expect_identical(load_hic(spec, wide, settings)$data$resolution, 5000)

  tmp <- tempfile(fileext = ".svg")
  on.exit(unlink(tmp), add = TRUE)
  rows <- test_rows("human")
  specs <- list(spec, single_track("bw", rows[["GM12878"]], settings = settings))
  suppressMessages(suppressWarnings(
    render_all(specs, region, settings, device_fn = function(w, h) svglite::svglite(tmp, width = w, height = h))))
  grDevices::dev.off()
  svg <- paste(readLines(tmp, warn = FALSE), collapse = "\n")
  # The whole message, not cut in the middle (render_all() cuts one over 140 characters). It is
  # wider than the track, so it is drawn as several lines (bug 77): read them as one.
  text <- paste(regmatches(svg, gregexpr("(?<=>)[^<]*(?=</text>)", svg, perl = TRUE))[[1]], collapse = " ")
  expect_match(text, paste("track failed:", gsub("&", "&amp;", msg, fixed = TRUE)), fixed = TRUE)
  expect_identical(lengths(regmatches(svg, gregexpr("track failed", svg, fixed = TRUE))), 1L)
})

# The figure's SVG, drawn with render_all(), as one string.
hic_svg <- function(specs, region, settings) {
  tmp <- tempfile(fileext = ".svg")
  on.exit(unlink(tmp), add = TRUE)
  suppressMessages(suppressWarnings(
    render_all(specs, region, settings, device_fn = function(w, h) svglite::svglite(tmp, width = w, height = h))))
  grDevices::dev.off()
  paste(readLines(tmp, warn = FALSE), collapse = "\n")
}

test_that("a Hi-C view whose maximum is below 0.00005 draws instead of failing (mouse dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  region <- test_region("mouse")
  settings <- test_settings("mouse")
  spec <- single_track("hic", test_rows("mouse")[["Ba"]], settings = settings)
  # The same view with every count scaled down so the colour scale's top (the 95th percentile) is
  # 0.00001: 4-decimal rounding used to give a colour range of c(0, 0), which plotgardener can't draw.
  real_load <- load_hic
  local_mocked_bindings(load_hic = function(entry, region_params, settings, pad = 0) {
    entry <- real_load(entry, region_params, settings, pad)
    m <- entry$data$matrix
    # Scaled by the 95th percentile of the contacts drawn (the track width of a page_w 9 figure).
    drawn <- .pt_hic_drawn_counts(m, "rectangle", region_params$start, region_params$end,
                                  entry$data$resolution, entry$height / (9 - 0.8 - 1.0))
    entry$data$matrix$counts <- m$counts / stats::quantile(drawn, 0.95, na.rm = TRUE, names = FALSE) * 1e-5
    entry
  })
  svg <- hic_svg(list(spec), region, settings)
  expect_false(grepl("track failed", svg, fixed = TRUE))
  expect_match(svg, ">0.00001<", fixed = TRUE)   # the colour legend's top
})

test_that("a Hi-C track that fails keeps its usual height (mouse dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  region <- test_region("mouse")
  row <- test_rows("mouse")[["Ba"]]
  height <- function(size) {
    settings <- hic_fixed_settings("mouse", size)
    svg <- hic_svg(list(single_track("hic", row, settings = settings)), region, settings)
    list(failed = grepl("track failed", svg, fixed = TRUE),
         height = regmatches(svg, regexpr("height='[0-9.]+pt'", svg)))
  }
  ok <- height(10000)
  failed <- height(12345)   # this file has no 12,345 bp bins
  expect_false(ok$failed)
  expect_true(failed$failed)
  expect_identical(failed$height, ok$height)
})

# A real .cool from the test data (1 Mb bins, hg38, chromosomes named "1", "2", ..., 455 contigs,
# the last 25 without contacts).
cool_spec <- function(settings) {
  spec <- single_track("hic", test_rows("human")[["GM12878"]], settings = settings)
  spec$path <- file.path(DATA_ROOT, "hic_test", "hg38_LEUK_HEK_1Mb.cool")
  spec
}

test_that("a .cool contig without contacts at the end of the file's list says so in the track (human dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  skip_if_not(file.exists(file.path(DATA_ROOT, "hic_test", "hg38_LEUK_HEK_1Mb.cool")))
  settings <- test_settings("human", list(hic_norm = "NONE"))
  region <- list(chrom = "UN_KI270394V1", start = 1L, end = 970L, gene = NA_character_)
  expect_error(load_hic(cool_spec(settings), region, settings), "this file has no contacts on UN_KI270394V1.", fixed = TRUE)
})

test_that("a Hi-C colour cap of 0 says so in the track, keeping its height (mouse dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  region <- test_region("mouse")
  row <- test_rows("mouse")[["Ba"]]
  draw <- function(cap) {
    settings <- test_settings("mouse", list(zrange_cap = cap))
    hic_svg(list(single_track("hic", row, settings = settings)), region, settings)
  }
  ok <- draw(100)
  zero <- draw(0)
  expect_false(grepl("track failed", ok, fixed = TRUE))
  expect_match(zero, "hic (Ba) track failed: the Hi-C colour cap must be above 0 (it is 0).", fixed = TRUE)
  height <- function(svg) regmatches(svg, regexpr("height='[0-9.]+pt'", svg))
  expect_identical(height(zero), height(ok))
})
