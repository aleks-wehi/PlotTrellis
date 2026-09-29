# Phase D.3 step 5: "settings-actually-wired QC" (phase-D.3.md step 15's
# last bullet) -- confirms a setting doesn't just reach the UI, it changes
# the render. Restricted to settings with a live consumer today:
# control_colour/treatment_colour (both comparison vizzes), diff_colour_lo/
# diff_colour_hi (pipeline/viz/bw_diff_stacked.R specifically),
# junction_opacity_cap (pipeline/viz/junctions_mirrored_diff.R specifically),
# and hic_res_mode/hic_max_bins_per_view (single-sample pipeline/viz/hic.R
# only -- no 2-sample hic comparison viz exists).
#
# CORRECTION (Step 5 adversarial review finding 1, 2026-09-01): this file
# originally also excluded diff_colour_lo/diff_colour_hi/
# junction_opacity_cap as "dead settings" per BUGS.md #16 -- that bug's
# claim was WRONG for these three (confirmed live: bw_diff_stacked.R DOES
# consume DIFF_COLOUR_HI/DIFF_COLOUR_LO, junctions_mirrored_diff.R DOES
# consume JUNCTION_OPACITY_CAP -- see BUGS.md #16's own correction note).
# Only diff_log2fc_cap is genuinely dead (zero code references anywhere)
# and stays excluded.
#
# Each case asserts a SPECIFIC, semantic signal (a colour hex string
# present in the assembled script text; a resolved bin/resolution value) --
# never a whole-file diff or an exact coordinate -- so a future pure
# spacing/layout change can't break these assertions (Step 1's
# "regeneration-friendly" design decision).

test_that("control_colour/treatment_colour are wired through bw_diff_stacked (mouse dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  rows          <- test_rows("mouse")
  region_params <- test_region("mouse")
  settings      <- test_settings("mouse")
  spec <- pair_track("bw_diff_stacked", rows[["Ba"]], rows[["LP"]],
                             extra = list(control_colour = "#123456", treatment_colour = "#ABCDEF"))
  expect_identical(spec$control_colour, "#123456")
  expect_identical(spec$treatment_colour, "#ABCDEF")

  script <- assemble_download_script(list(spec), region_params, settings)
  expect_true(grepl("123456", script, ignore.case = TRUE))
  expect_true(grepl("ABCDEF", script, ignore.case = TRUE))

  fig <- test_fig(list(spec), region_params, settings,
                          perf_id = "qc-bw-diff-stacked-colour-mouse", module = "bw_diff_stacked")
  vdiffr::expect_doppelganger("qc-bw-diff-stacked-colour-mouse", fig)
})

test_that("control_colour/treatment_colour are wired through junctions_mirrored_diff (human dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  rows          <- test_rows("human")
  region_params <- test_region("human")
  settings      <- test_settings("human")
  spec <- pair_track("junctions_mirrored_diff", rows[["GM12878"]], rows[["K562"]],
                             extra = list(control_colour = "#123456", treatment_colour = "#ABCDEF"))
  expect_identical(spec$control_colour, "#123456")
  expect_identical(spec$treatment_colour, "#ABCDEF")

  script <- assemble_download_script(list(spec), region_params, settings)
  expect_true(grepl("123456", script, ignore.case = TRUE))
  expect_true(grepl("ABCDEF", script, ignore.case = TRUE))

  fig <- test_fig(list(spec), region_params, settings,
                          perf_id = "qc-junctions-mirrored-diff-colour-human", module = "junctions_mirrored_diff")
  vdiffr::expect_doppelganger("qc-junctions-mirrored-diff-colour-human", fig)
})

test_that("diff_colour_lo/diff_colour_hi are wired through bw_diff_stacked (mouse dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  rows          <- test_rows("mouse")
  region_params <- test_region("mouse")
  settings      <- test_settings("mouse", overrides = list(diff_colour_lo = "#111111", diff_colour_hi = "#222222"))
  spec          <- pair_track("bw_diff_stacked", rows[["Ba"]], rows[["LP"]])

  script <- assemble_download_script(list(spec), region_params, settings)
  expect_true(grepl("111111", script, ignore.case = TRUE))
  expect_true(grepl("222222", script, ignore.case = TRUE))

  fig <- test_fig(list(spec), region_params, settings,
                          perf_id = "qc-bw-diff-stacked-diffcolour-mouse", module = "bw_diff_stacked")
  vdiffr::expect_doppelganger("qc-bw-diff-stacked-diffcolour-mouse", fig)
})

test_that("junction_opacity_cap is wired through junctions_mirrored_diff (human dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  rows          <- test_rows("human")
  region_params <- test_region("human")
  settings      <- test_settings("human", overrides = list(junction_opacity_cap = 0.5))
  spec          <- pair_track("junctions_mirrored_diff", rows[["GM12878"]], rows[["K562"]])

  # settings gets deparse()d verbatim into the assembled download script's
  # data block (utils/assemble.R) -- render_junctions_mirrored_diff() itself
  # only reads ctx$junction_opacity_cap (render_all()'s per-track context,
  # set from settings$junction_opacity_cap), so the literal only appears
  # once, at the settings list's own declaration.
  script <- assemble_download_script(list(spec), region_params, settings)
  expect_true(grepl("junction_opacity_cap = 0.5", script, fixed = TRUE))

  fig <- test_fig(list(spec), region_params, settings,
                          perf_id = "qc-junctions-mirrored-diff-opacitycap-human", module = "junctions_mirrored_diff")
  vdiffr::expect_doppelganger("qc-junctions-mirrored-diff-opacitycap-human", fig)
})

# hic_res_mode/hic_max_bins_per_view: reuse test-hic-resolution.R's
# load_hic()-direct-call mechanism (see that file's own header comment).
hic_setting_check <- function(dataset, sample, chrom, start, end, settings) {
  rows <- test_rows(dataset)
  region_params <- list(chrom = chrom, start = start, end = end, gene = NA_character_)
  entry <- single_track("hic", rows[[sample]], settings = settings)
  entry <- load_hic(entry, region_params, settings)
  entry$data$resolution
}

# No doppelganger for the two cases below (deliberate, not an oversight):
# Step 4 already committed a 108MB baseline SVG for this identical
# large-span mouse window (test-hic-resolution.R's large-span case), and a
# second/third 100MB+ binary just to re-prove settings threading would be
# pure repo bloat. What IS proven here is that the setting reaches
# load_hic()'s own resolution-selection logic (the resolution differs).
# CAVEAT (Step 5 adversarial review finding 2, confirmed live): for the
# "finest" case specifically, plotgardener's own pixel-count safety valve
# ("Attempting to plot too many Hi-C pixels. Adjusting to a resolution of
# 50000 BP.") silently re-coarsens the ACTUAL displayed matrix back to the
# same 50000bp grid as "auto" on this large a span -- so the resolution
# load_hic() picked differing is real and correctly asserted, but
# don't read that as proof the two renders look visibly different at this
# span; they don't. hic_res_mode's effect IS visibly real at a span where
# plotgardener's own coarsening doesn't kick in (any span small enough
# that "finest" stays within plotgardener's own pixel budget) -- not
# reproduced here to avoid a 3rd giant baseline.
test_that("hic_res_mode is wired through: auto vs finest pick different resolutions (mouse dataset, large span)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  settings_auto   <- test_settings("mouse")  # hic_res_mode="auto" (default)
  settings_finest <- test_settings("mouse", overrides = list(hic_res_mode = "finest"))
  res_auto   <- hic_setting_check("mouse", "Ba", "chr2", 80000000, 110000000, settings_auto)
  res_finest <- hic_setting_check("mouse", "Ba", "chr2", 80000000, 110000000, settings_finest)
  expect_identical(res_auto, 50000)    # confirmed live, test-hic-resolution.R's large-span mouse case
  expect_identical(res_finest, 10000)  # finest available resolution for this fixture
  expect_false(identical(res_auto, res_finest))
})

test_that("hic_max_bins_per_view is wired through: lowering the threshold picks a coarser tier (mouse dataset, large span)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  settings_default <- test_settings("mouse")  # hic_max_bins_per_view=500 (default)
  settings_lower   <- test_settings("mouse", overrides = list(hic_max_bins_per_view = 50))
  res_default <- hic_setting_check("mouse", "Ba", "chr2", 80000000, 110000000, settings_default)
  res_lower   <- hic_setting_check("mouse", "Ba", "chr2", 80000000, 110000000, settings_lower)
  expect_identical(res_default, 50000)
  expect_identical(res_lower, 200000)  # confirmed live: lowering the bin-count floor lets a coarser resolution qualify
  expect_false(identical(res_default, res_lower))
})
