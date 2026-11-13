# One correct-render golden image per module, per fixture dataset
# (mouse/human), so every module has at least one correct-render test. The
# sibling files (cross-module composition, non-default settings, Hi-C
# resolution tiers, settings-wired QC) build on top of this grid, they
# don't replace it.
#
# Every region_params/rows/settings binding below is fresh inside its own
# test_that() (never shared across blocks) -- plain style hygiene;
# render_all() takes region_params/settings as explicit arguments and has
# no global-state side effects.
#
# bed gap: mouse_mammary_test has no standalone BED-interval file (only
# bigWig ATAC coverage for mouse, no peak calls) -- only the human
# dataset's ATAC IDR-peaks BED exercises this module. A known fixture gap
# (not papered over with fabricated mouse data); revisit if a mouse peaks
# BED ever lands.
#
# hic/human uses a wider region (flank_pct=25, i.e. a 2500% flank
# MULTIPLIER, not a literal percent -- don't "fix" this down to 0.25) --
# ~4Mb around PTPRC, instead of the dataset default. The human .hic
# fixture is CTCF ChIA-PET data,
# much sparser than a real Hi-C matrix, and the default gene-scale window
# (~150kb) has no bins with any signal, so it would only draw an empty track
# (test-soft-fail.R covers that). flank_pct=25 is the first value of
# {2, 10, 25, 50} that draws a matrix.
single_cases <- list(
  list(module = "bw",                    dataset = "mouse", sample = "Ba"),
  list(module = "bw",                    dataset = "human", sample = "GM12878"),
  list(module = "genes",                 dataset = "mouse", sample = "Ba"),
  list(module = "genes",                 dataset = "human", sample = "GM12878"),
  list(module = "hic",                   dataset = "mouse", sample = "Ba"),
  list(module = "hic",                   dataset = "human", sample = "GM12878", flank_pct = 25),
  list(module = "protein_domains",       dataset = "mouse", sample = "Ba"),
  list(module = "protein_domains",       dataset = "human", sample = "GM12878"),
  list(module = "junctions",             dataset = "human", sample = "GM12878"),
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

# junctions has no mouse case: the mouse dataset's "Junctions" slot holds Hi-C loop
# calls (BEDPE), not splice junctions. Human's BED12 junctions cover it above.

# ---- Comparison (two-sample) modules ---------------------------------------
# junctions_mirrored_diff needs both control/treatment to actually carry
# junction data -- mouse's LP has none (see helper-fixtures.R), so its pair
# uses Ba (control) vs P6 (treatment) instead of the Ba/LP pairing used for
# bw_diff_stacked (bw coverage is universal across mouse conditions).
#
# junctions_mirrored_diff_tpm_scaled needs both sides to carry junction AND
# quantdir data (its per-side gene-TPM lookup requires quantdir). No mouse
# pair satisfies both: Ba has junction+quantdir but LP/P6/P7 each have only
# one of the two (see helper-fixtures.R's FIXTURES$mouse). Human-only
# coverage below, same principle as the bed/junctions mouse gaps above --
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
