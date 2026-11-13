# Cross-module composition fixtures, on top of test-modules.R's per-module grid.
#
# junctions/bed are excluded from every MOUSE composition below (same
# reasons as test-modules.R: no BED-interval file, and the mouse
# "Junctions" data are Hi-C loop calls, not splice junctions) -- both are
# fully available on the human dataset instead.
#
# Every region_params binding is fresh per test_that() -- plain style
# hygiene (render_all() has no global-state side effects).

# ---- Cross-module fixture: non-scale-bar module first ---------------------
# Mixes several track types in one figure with the scale-bar module not
# first. "is_first" (the scale bar's position-sensitive draw, render_all()'s
# ctx$is_first) is track-order-sensitive -- see the "large multi-track
# composition 2" fixtures below, which deliberately put a scale-bar
# module (scale_bar in its info_<key>()) first for exactly that reason (neither bug57 fixture does, by
# design, since they're testing the no-scale-bar-first case instead).

test_that("cross-module fixture with no scale-bar module first renders correctly (mouse dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  rows          <- test_rows("mouse")
  region_params <- test_region("mouse")
  settings      <- test_settings("mouse")
  specs <- list(
    single_track("genes", rows[["Ba"]]),
    pair_track("junctions_mirrored_diff", rows[["Ba"]], rows[["P6"]]),
    single_track("bw", rows[["Ba"]])
  )
  fig <- test_fig(specs, region_params, settings, perf_id = "bug57-mouse", module = "composition")
  vdiffr::expect_doppelganger("bug57-mouse", fig)
})

test_that("cross-module fixture with no scale-bar module first renders correctly (human dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  rows          <- test_rows("human")
  region_params <- test_region("human")
  settings      <- test_settings("human")
  specs <- list(
    single_track("genes", rows[["GM12878"]]),
    pair_track("bw_diff_stacked", rows[["GM12878"]], rows[["K562"]]),
    single_track("junctions", rows[["GM12878"]])
  )
  fig <- test_fig(specs, region_params, settings, perf_id = "bug57-human", module = "composition")
  vdiffr::expect_doppelganger("bug57-human", fig)
})

# Render-determinism: the mixed cross-module fixture is the one most likely
# to be order/state-sensitive, so prove vdiffr's "one baseline is
# enough" assumption actually holds for it: render it twice in one session,
# diff the raw SVG bytes.
test_that("cross-module fixture renders deterministically (mouse dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  rows          <- test_rows("mouse")
  region_params <- test_region("mouse")
  settings      <- test_settings("mouse")
  specs <- list(
    single_track("genes", rows[["Ba"]]),
    pair_track("junctions_mirrored_diff", rows[["Ba"]], rows[["P6"]]),
    single_track("bw", rows[["Ba"]])
  )
  fig <- test_fig(specs, region_params, settings)
  f1 <- tempfile(fileext = ".svg"); f2 <- tempfile(fileext = ".svg")
  svglite::svglite(f1); fig(); dev.off()
  svglite::svglite(f2); fig(); dev.off()
  expect_identical(unname(tools::md5sum(f1)), unname(tools::md5sum(f2)))
})

# ---- No-"genes"-track fixture (project-wide, one is enough) ----------------
# render_all() always appends genome_axis as its own track, whether or not
# a "genes" track is present. A regression check that nothing in the render
# path silently assumes a "genes" track exists somewhere in the figure.
test_that("multi-track fixture with no genes track renders correctly (mouse dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  rows          <- test_rows("mouse")
  region_params <- test_region("mouse")
  settings      <- test_settings("mouse")
  specs <- list(
    single_track("transcripts", rows[["Ba"]]),
    single_track("bw", rows[["Ba"]]),
    single_track("protein_domains", rows[["Ba"]]),
    single_track("tpm_usage", rows[["Ba"]])
  )
  fig <- test_fig(specs, region_params, settings, perf_id = "no-genes-track-mouse", module = "composition")
  vdiffr::expect_doppelganger("no-genes-track-mouse", fig)
})

# ---- Larger multi-track compositions, order deliberately varied -----------
# Not exhaustive permutation testing -- just enough real variation (which
# module leads, which scale-bar module is where) to exercise
# order-dependent failure modes beyond the cross-module fixture above.

test_that("large multi-track composition 1 renders correctly (mouse dataset, scale-bar module in middle)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  rows          <- test_rows("mouse")
  region_params <- test_region("mouse")
  settings      <- test_settings("mouse")
  specs <- list(
    single_track("transcripts", rows[["Ba"]]),
    single_track("protein_domains", rows[["Ba"]]),
    single_track("bw", rows[["Ba"]]),
    single_track("tpm_usage", rows[["Ba"]]),
    single_track("transcripts_annotated", rows[["Ba"]])
  )
  fig <- test_fig(specs, region_params, settings, perf_id = "multitrack1-mouse", module = "composition")
  vdiffr::expect_doppelganger("multitrack1-mouse", fig)
})

test_that("large multi-track composition 2 renders correctly (mouse dataset, scale-bar module first)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  rows          <- test_rows("mouse")
  region_params <- test_region("mouse")
  settings      <- test_settings("mouse")
  specs <- list(
    pair_track("bw_diff_stacked", rows[["Ba"]], rows[["LP"]]),
    single_track("genes", rows[["Ba"]]),
    single_track("transcripts_all", rows[["Ba"]]),
    single_track("hic", rows[["Ba"]]),
    pair_track("junctions_mirrored_diff", rows[["Ba"]], rows[["P6"]])
  )
  fig <- test_fig(specs, region_params, settings, perf_id = "multitrack2-mouse", module = "composition")
  vdiffr::expect_doppelganger("multitrack2-mouse", fig)
})

test_that("large multi-track composition 1 renders correctly (human dataset, scale-bar module in middle)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  rows          <- test_rows("human")
  region_params <- test_region("human")
  settings      <- test_settings("human")
  specs <- list(
    single_track("transcripts", rows[["GM12878"]]),
    single_track("protein_domains", rows[["GM12878"]]),
    single_track("bw", rows[["GM12878"]]),
    single_track("tpm_usage", rows[["GM12878"]]),
    single_track("bed", rows[["GM12878"]]),
    single_track("junctions", rows[["GM12878"]])
  )
  fig <- test_fig(specs, region_params, settings, perf_id = "multitrack1-human", module = "composition")
  vdiffr::expect_doppelganger("multitrack1-human", fig)
})

# flank_pct=25 here too (see test-modules.R's hic/human comment) -- this
# composition includes hic, and all tracks in one figure share one region.
test_that("large multi-track composition 2 renders correctly (human dataset, different scale-bar module last)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  rows          <- test_rows("human")
  region_params <- test_region("human", flank_pct = 25)
  settings      <- test_settings("human")
  specs <- list(
    single_track("hic", rows[["GM12878"]]),
    single_track("genes", rows[["GM12878"]]),
    single_track("transcripts_all", rows[["GM12878"]]),
    single_track("transcripts_annotated", rows[["GM12878"]]),
    pair_track("bw_diff_stacked", rows[["GM12878"]], rows[["K562"]])
  )
  fig <- test_fig(specs, region_params, settings, perf_id = "multitrack2-human", module = "composition")
  vdiffr::expect_doppelganger("multitrack2-human", fig)
})
