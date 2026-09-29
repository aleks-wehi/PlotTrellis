# Phase D.3 step 4b: Hi-C auto-resolution-tier fixtures, small-span and
# large-span per Hi-C-having dataset (both mouse and human now have real
# per-condition .hic files, see helper-fixtures.R's FIXTURES table).
# Asserts Phase A.2's auto-resolution logic (pipeline/loaders/load_hic.R's
# load_hic(), "auto" hic_res_mode) actually picks the expected tier for each
# span, not just that a render happens to succeed either way.
#
# Mechanism: calls pipeline/loaders/load_hic.R's load_hic() directly, the
# same resolution-selection logic render_all() invokes internally via
# settings_hic() -- avoids a full render just to read back which tier was
# picked. `avail` is recomputed with the same two lines load_hic() itself
# uses (it doesn't expose the candidate list on the entry).
#
# The mouse case doubles as the first real-pipeline verification of the
# mouse .hic files' chrom-naming convention (data/mouse_mammary_test/
# README.md flagged this as unverified) -- confirmed working live
# 2026-09-01, no BUGS.md #12 chrom-naming failure.

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
  expect_identical(result$res, 50000)  # confirmed live 2026-09-01 against the real fixture

  rows     <- test_rows("mouse")
  settings <- test_settings("mouse")
  spec     <- single_track("hic", rows[["Ba"]], settings = settings)
  fig      <- test_fig(list(spec), region_params, settings,
                               perf_id = "hic-tier-large-mouse", module = "hic")
  vdiffr::expect_doppelganger("hic-tier-large-mouse", fig)
})

# Not a true gene-scale (~140kb) span like the mouse small-span case above --
# see BUGS.md #23's severity-raised note: the human .hic fixture's declared
# finest resolution (1000bp) has ZERO usable KR-normalized signal anywhere
# in the file (confirmed by direct strawr scan, not guessed), so the
# default gene-scale window always falls back to that broken tier and
# crashes with bug 14's error, regardless of which small window is chosen.
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
  expect_identical(result$res, 5000)  # confirmed live 2026-09-01 against the real fixture

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
  expect_identical(result$res, 50000)  # confirmed live 2026-09-01 against the real fixture

  rows     <- test_rows("human")
  settings <- test_settings("human")
  spec     <- single_track("hic", rows[["GM12878"]], settings = settings)
  fig      <- test_fig(list(spec), region_params, settings,
                               perf_id = "hic-tier-large-human", module = "hic")
  vdiffr::expect_doppelganger("hic-tier-large-human", fig)
})
