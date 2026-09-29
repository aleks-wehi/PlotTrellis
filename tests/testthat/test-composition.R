# Phase D.3 step 3: cross-module composition fixtures (phase-D.3.md step
# 15's bullets 1-3), on top of Step 2's per-module grid.
#
# sashimi/bed are excluded from every MOUSE composition below (same
# reasons as test-modules.R: no BED-interval file at all, and every real
# loop-derived "Junctions" region hits BUGS.md #19/#20) -- both are fully
# available on the human dataset instead.
#
# Every region_params binding is fresh per test_that() -- plain style
# hygiene, kept from the pre-cutover suite (render_all() has no
# global-state side effects, so this is no longer a correctness necessity).

# ---- Bug 5/7 cross-module fixture: non-SCALE_BAR_TYPES module first --------
# Historical rationale (pre-cutover, viz/*.R + viz_comparisons/*.R
# architecture): per-module-isolated tests alone couldn't catch a clobber
# between the two families' duplicated copies of bioeda_mixed_text()/the
# scale bar, since whichever copy's setup_code() was sourced last won for
# the WHOLE figure. That specific failure mode is structurally impossible
# now -- pipeline/_boilerplate.R is the single shared source of
# mixed_text()/draw_scale_bar() every viz file uses, no duplicated copies
# left to clobber each other. Kept as a general cross-module composition
# fixture (still real coverage: mixing many track types in one figure,
# scale-bar attachment on a non-first track) even though its original
# clobber-catching purpose no longer applies. "is_first" (the scale bar's
# position-sensitive draw, pipeline/_render_all.R's ctx$is_first) IS still
# genuinely track_specs-position-sensitive -- see the "large multi-track
# composition 2" fixtures below, which deliberately put a SCALE_BAR_TYPES
# module first for exactly that reason (neither bug57 fixture does, by
# design, since they're testing the no-scale-bar-first case instead).

test_that("cross-module bug 5/7 fixture renders correctly (mouse dataset)", {
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

test_that("cross-module bug 5/7 fixture renders correctly (human dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  rows          <- test_rows("human")
  region_params <- test_region("human")
  settings      <- test_settings("human")
  specs <- list(
    single_track("genes", rows[["GM12878"]]),
    pair_track("bw_diff_stacked", rows[["GM12878"]], rows[["K562"]]),
    single_track("sashimi", rows[["GM12878"]])
  )
  fig <- test_fig(specs, region_params, settings, perf_id = "bug57-human", module = "composition")
  vdiffr::expect_doppelganger("bug57-human", fig)
})

# Extra check 1 (phase log): render-determinism -- this is the fixture
# most likely to be order/state-sensitive (it deliberately exercises the
# still-unfixed duplicated-helper risk), so prove vdiffr's "one baseline is
# enough" assumption actually holds for it: render it twice in one session,
# diff the raw SVG bytes.
test_that("bug 5/7 fixture renders deterministically (mouse dataset)", {
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
# Historical rationale (pre-cutover): the OLD codegen's compute_axis_index()
# had a distinct fallback (attach the axis to the LAST track) when no
# "genes"-type track was present. render_all() no longer has that branch at
# all -- genome_axis is now always its own appended track (pipeline/
# _render_all.R), regardless of whether "genes" is present. Kept anyway as a
# general regression check that nothing in the pipeline (title generation's
# types-minus-"genes" filter, GENE_REQUIRED_TYPES handling, etc.) silently
# assumes a "genes" track exists somewhere in the figure.
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
# module leads, which SCALE_BAR_TYPES module is where) to exercise
# order-dependent failure modes beyond the single bug5/7 fixture above.

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
    single_track("sashimi", rows[["GM12878"]])
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
