# Phase D.4 step 6/7 originally wrote these as testthat::expect_error() TDD
# tripwires: no "no data"/"failed" track appearance existed yet (that was
# Phase E.4's job), so a crash was the correctly-asserted expected behavior
# at the time -- see git history for that original version if needed.
#
# Phase E.4 step 12 fixed each underlying crash (viz/hic.R, viz/sashimi.R,
# viz/bw.R now catch their own failure and render the shared "No data"
# placeholder instead of throwing), so these are now real
# vdiffr::expect_doppelganger() golden-image tests against that placeholder
# look, per phase-D.4.md step 6 and this phase's own kickoff ambiguity 1
# (solutions/phase-logs/phase-E.4.md).
#
# A fourth case is added here (not one of phase-D.4.md's original three):
# a GENE-required module's own is.na(GENE) guard, previously an explicit
# stop() -- Phase E.4's kickoff audit found this gap (its own Verification
# block names "a deliberately... gene-less region" as a required check,
# but D.4's three cases didn't cover it) and its kickoff question 3
# confirmed adding a case for it. transcripts_all.R is used as the
# representative case (simplest of the six GENE-required modules); the
# other five (transcripts.R, transcripts_annotated.R, protein_domains.R,
# tpm_usage.R, junctions_mirrored_diff_tpm_scaled.R) share the identical
# is.na(GENE) guard pattern and were live-verified during implementation
# (see solutions/phase-logs/phase-E.4.md) but aren't each given their own
# golden image here, to avoid five near-duplicate snapshots of the same
# placeholder text.
#
# One case per category, each a real, live-verified crash-turned-render
# (probed with module load R/4.5.3, NOT_CRAN=true, directly against the D.2
# fixtures before being written here) -- not a hypothetical or a full
# module x dataset grid, see kickoff ambiguity 2 (phase-D.4.md).

test_that("out-of-bounds region: hic renders 'No data' instead of crashing (BUGS.md #14, human)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  rows_h    <- test_rows("human")
  settings_h <- test_settings("human")
  # Genuinely past the end of any real chr1 (hg38 chr1 is ~249Mb) -- not a
  # sparse-but-in-range locus like BUGS.md #14's original PTPRC finding,
  # a coordinate that doesn't exist on the chromosome at all.
  oob_region <- list(chrom = "chr1", start = 990000000, end = 990010000, gene = "PTPRC")
  spec <- single_track("hic", rows_h[["GM12878"]])
  fig  <- test_fig(list(spec), oob_region, settings_h)
  vdiffr::expect_doppelganger("hic-oob-human-gm12878", fig)
})

test_that("empty GRanges: sashimi renders 'No data' instead of crashing (BUGS.md #19, mouse)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  rows_m    <- test_rows("mouse")
  region_m  <- test_region("mouse")
  settings_m <- test_settings("mouse")
  # Mouse's "Junctions" slot is repurposed OmniC/FitHiChIP loop data (BEDPE) --
  # the Ba/Elf5 window has zero overlapping loop calls, hitting
  # pool_junctions()'s zero-row data.frame() crash directly (same case
  # test-modules.R's sashimi/mouse skip() already documents as blocked).
  spec <- single_track("sashimi", rows_m[["Ba"]])
  fig  <- test_fig(list(spec), region_m, settings_m)
  vdiffr::expect_doppelganger("sashimi-empty-mouse-ba", fig)
})

test_that("missing data type: bw renders 'No data' instead of crashing (BUGS.md #24, mouse)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  rows_m    <- test_rows("mouse")
  region_m  <- test_region("mouse")
  settings_m <- test_settings("mouse")
  # Simulates a project_settings.json referencing a bigWig that's since been
  # moved/deleted -- nothing stops render_all() itself from being handed a
  # bad path. `path` is the real entry field name (pipeline/viz/bw.R).
  spec <- single_track("bw", rows_m[["Ba"]])
  spec$path <- file.path(FIXTURES$mouse$data_dir, "atac", "DOES_NOT_EXIST.bw")
  fig <- test_fig(list(spec), region_m, settings_m)
  vdiffr::expect_doppelganger("bw-missing-file-mouse-ba", fig)
})

test_that("gene-less region: transcripts_all renders 'No data' instead of crashing (Phase E.4 audit finding, coords mode with no gene, human)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  rows_h     <- test_rows("human")
  settings_h <- test_settings("human")
  # Coords mode with no gene resolved -- exactly the state app.R's
  # GENE_REQUIRED_TYPES_JS red-box guard (Phase E.2/E.3) exists to prevent a
  # biologist from reaching in the UI; this test bypasses that guard
  # directly to exercise the module's own defense-in-depth path.
  no_gene_region <- list(chrom = "chr1", start = 1000000, end = 1010000, gene = NA_character_)
  spec <- single_track("transcripts_all", rows_h[["GM12878"]])
  fig  <- test_fig(list(spec), no_gene_region, settings_h)
  vdiffr::expect_doppelganger("transcripts-all-gene-less-human-gm12878", fig)
})
