# Hi-C auto-resolution-tier fixtures, small-span and large-span per
# Hi-C-having dataset (both mouse and human have real per-condition .hic
# files, see helper-fixtures.R's FIXTURES table). Asserts the auto-resolution
# logic (load_hic() in loader-hic.R, "auto" hic_res_mode) actually picks the
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
  region <- list(chrom = "chr2", start = 80000000, end = 110000000, gene = NA_character_)   # auto picks 50000 here
  auto <- hic_tier_check("mouse", "Ba", region$chrom, region$start, region$end)
  want <- min(auto$avail)                                                                   # not what auto picks
  expect_false(identical(auto$res, want))
  settings <- hic_fixed_settings("mouse", want)
  entry <- load_hic(single_track("hic", test_rows("mouse")[["Ba"]], settings = settings), region, settings)
  expect_identical(entry$data$resolution, want)
  expect_false(isTRUE(entry$data$failed))
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
  expect_error(load_hic(spec, region, hic_fixed_settings("mouse", NULL)), "no bin size is set", fixed = TRUE)
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
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
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

