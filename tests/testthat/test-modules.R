# Phase D.3 step 2: one correct-render golden image per surviving module,
# per fixture dataset (mouse/human) -- the phase's own "every surviving
# module has at least one correct-render test" bar. Later steps in this
# file's siblings (cross-module composition, non-default settings, Hi-C
# resolution tiers, settings-wired QC) build on top of this grid, they
# don't replace it.
#
# Every region_params/rows/settings binding below is fresh inside its own
# test_that() (never shared across blocks) -- plain style hygiene, kept from
# the pre-cutover suite; render_all() takes region_params/settings as
# explicit arguments and has no global-state side effects, so cross-block
# reuse would no longer be a correctness risk the way it once was.
#
# viz/bed.R gap: mouse_mammary_test has no standalone BED-interval file
# (D.2 only produced bigWig ATAC coverage for mouse, no peak calls) -- only
# the human dataset's ATAC IDR-peaks BED exercises this module. Logged as
# an open, known fixture gap (not papered over with fabricated mouse data);
# revisit if a mouse peaks BED ever lands.
#
# hic/human uses a wider region (flank_pct=25, i.e. a 2500% flank
# MULTIPLIER, not a literal percent -- don't "fix" this down to 0.25) --
# ~4Mb around PTPRC, instead of the dataset default. The human .hic
# fixture is CTCF ChIA-PET data,
# much sparser than a real Hi-C matrix, and the default gene-scale window
# (~150kb) has no bins with any signal, hitting the pre-existing empty-
# matrix crash logged as BUGS.md #14 ("Cannot add heatmap legend to a plot
# without a colour palette"). Confirmed live 2026-09-01: flank_pct=25 is
# the first value of {2, 10, 25, 50} that actually renders. Out of D.3's
# surgical scope to fix #14 itself; this just picks a region where the
# real bug doesn't fire, same principle as sashimi's mouse case below.
single_cases <- list(
  list(module = "bw",                    dataset = "mouse", sample = "Ba"),
  list(module = "bw",                    dataset = "human", sample = "GM12878"),
  list(module = "genes",                 dataset = "mouse", sample = "Ba"),
  list(module = "genes",                 dataset = "human", sample = "GM12878"),
  list(module = "hic",                   dataset = "mouse", sample = "Ba"),
  list(module = "hic",                   dataset = "human", sample = "GM12878", flank_pct = 25),
  list(module = "protein_domains",       dataset = "mouse", sample = "Ba"),
  list(module = "protein_domains",       dataset = "human", sample = "GM12878"),
  list(module = "sashimi",               dataset = "human", sample = "GM12878"),
  list(module = "tpm_usage",             dataset = "mouse", sample = "Ba"),
  list(module = "tpm_usage",             dataset = "human", sample = "GM12878"),
  list(module = "transcripts",           dataset = "mouse", sample = "Ba"),
  list(module = "transcripts",           dataset = "human", sample = "GM12878"),
  list(module = "transcripts_all",       dataset = "mouse", sample = "Ba"),
  list(module = "transcripts_all",       dataset = "human", sample = "GM12878"),
  list(module = "transcripts_annotated", dataset = "mouse", sample = "Ba"),
  list(module = "transcripts_annotated", dataset = "human", sample = "GM12878"),
  list(module = "bed",                   dataset = "human", sample = "GM12878")
)

for (cs in single_cases) {
  local({
    module    <- cs$module
    dataset   <- cs$dataset
    sample    <- cs$sample
    flank_pct <- if (is.null(cs$flank_pct)) 0.15 else cs$flank_pct
    title     <- paste0(gsub("_", "-", module), "-", dataset, "-", tolower(sample))
    test_that(paste0(module, " renders correctly (", dataset, " dataset)"), {
      skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
      rows          <- test_rows(dataset)
      region_params <- test_region(dataset, flank_pct = flank_pct)
      settings      <- test_settings(dataset)
      spec          <- single_track(module, rows[[sample]])
      fig           <- test_fig(list(spec), region_params, settings, perf_id = title, module = module)
      vdiffr::expect_doppelganger(title, fig)
    })
  })
}

# sashimi's mouse case is a known, logged, pre-existing gap, not a
# golden-image case: the mouse dataset's only "Junctions"-slot data is
# repurposed OmniC/FitHiChIP Hi-C loop calls (BEDPE), not real splice
# junctions. Every window with a real loop call hits BUGS.md #20 (crashes
# formatting a non-integer nSupport with "%d"); an empty window hits #19
# (zero-row data.frame() crash). Both are pre-existing bugs, out of this
# phase's surgical scope. Human's real BED12 junction data already gives
# sashimi its required correct-render coverage above.
test_that("sashimi mouse dataset: known-blocked by pre-existing bugs 19/20, not a golden-image case", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  skip(paste(
    "mouse fixture's Junctions slot is repurposed HiChIP loop data with",
    "continuous (non-integer) scores -- any region with a real loop hits",
    "BUGS.md #20's heading-format crash, an empty region hits #19's",
    "zero-row crash. Pre-existing bugs, not fixed by this phase."
  ))
})

# ---- Pair (viz_comparisons/*.R) modules -------------------------------------
# junctions_mirrored_diff needs both control/treatment to actually carry
# junction data -- mouse's LP has none (see helper-fixtures.R), so its pair
# uses Ba (control) vs P6 (treatment) instead of the Ba/LP pairing used for
# bw_diff_stacked (bw coverage is universal across mouse conditions).
#
# junctions_mirrored_diff_tpm_scaled needs both sides to carry junction AND
# quantdir data (its per-side gene-TPM lookup requires quantdir). No mouse
# pair satisfies both: Ba has junction+quantdir but LP/P6/P7 each have only
# one of the two (see helper-fixtures.R's FIXTURES$mouse). Human-only
# coverage below, same principle as the bed/sashimi mouse gaps above --
# a real fixture gap, not papered over.

pair_cases <- list(
  list(module = "bw_diff_stacked",         dataset = "mouse", control = "Ba",      treatment = "LP"),
  list(module = "bw_diff_stacked",         dataset = "human", control = "GM12878", treatment = "K562"),
  list(module = "junctions_mirrored_diff", dataset = "mouse", control = "Ba",      treatment = "P6"),
  list(module = "junctions_mirrored_diff", dataset = "human", control = "GM12878", treatment = "K562"),
  list(module = "junctions_mirrored_diff_tpm_scaled", dataset = "human", control = "GM12878", treatment = "K562")
)

for (cs in pair_cases) {
  local({
    module    <- cs$module
    dataset   <- cs$dataset
    control   <- cs$control
    treatment <- cs$treatment
    title     <- paste0(gsub("_", "-", module), "-", dataset, "-pair")
    test_that(paste0(module, " pair renders correctly (", dataset, " dataset)"), {
      skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
      rows          <- test_rows(dataset)
      region_params <- test_region(dataset)
      settings      <- test_settings(dataset)
      spec          <- pair_track(module, rows[[control]], rows[[treatment]])
      fig           <- test_fig(list(spec), region_params, settings, perf_id = title, module = module)
      vdiffr::expect_doppelganger(title, fig)
    })
  })
}
