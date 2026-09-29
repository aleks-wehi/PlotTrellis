# Phase D.4 step 8: writes _perf_baseline.csv the first time this suite runs
# clean with every correct-render fixture instrumented (test_fig()'s
# perf_id/module params, see helper-fixtures.R). Named to sort alphabetically
# last among test-*.R files (testthat::test_dir() runs them in filename
# order), so every fixture above has already recorded into PERF_TIMES by the
# time this file runs.
#
# Never overwrites an existing baseline -- once _perf_baseline.csv exists,
# every later run just asserts against it (2x the recorded seconds, Phase D
# Q3) inline at each fixture's own render, inside record_perf() itself.
# A deliberate re-baseline is a manual `rm` + re-run, not something a normal
# test run does automatically.
#
# The write guard below is "file absent" plus a completeness check against
# the vdiffr _snaps/ SVG count (an independent, already-committed source of
# truth for how many golden-image fixtures exist) -- NOT a check that every
# OTHER test in this run passed. A render that errors skips its own
# record_perf() call (system.time() never completes), so a partial/
# crashed run undercounts against _snaps/ and fails loudly here instead of
# silently committing an incomplete baseline. Found by Step 8's adversarial
# review (finding 1, phase-D.4.md); a vdiffr *mismatch* (not a crash) still
# records its timing before the doppelganger comparison fails, so this is a
# completeness check, not a full pass/fail gate -- treat a fresh
# _perf_baseline.csv as trustworthy only when this test run was otherwise green.
test_that("performance baseline is captured / up to date", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  recorded <- PERF_TIMES$records
  skip_if(length(recorded) == 0, "no fixtures recorded a perf_id this run")
  df <- do.call(rbind, lapply(recorded, function(r) data.frame(
    fixture = r$fixture, module = r$module, seconds = round(r$seconds, 3), stringsAsFactors = FALSE
  )))

  if (is.null(PERF_BASELINE)) {
    snaps_dir <- testthat::test_path("_snaps")
    n_snaps <- length(list.files(snaps_dir, pattern = "\\.svg$", recursive = TRUE))
    if (nrow(df) == n_snaps) {
      write.csv(df, PERF_BASELINE_PATH, row.names = FALSE)
      succeed(sprintf("Wrote new _perf_baseline.csv with %d fixtures (first clean run).", nrow(df)))
    } else {
      fail(sprintf(
        "recorded %d fixture(s) this run but _snaps/ has %d golden-image SVGs -- a render likely crashed or was skipped; refusing to write an incomplete _perf_baseline.csv",
        nrow(df), n_snaps
      ))
    }
  } else {
    missing <- setdiff(df$fixture, PERF_BASELINE$fixture)
    expect_true(
      length(missing) == 0,
      info = sprintf(
        "%d fixture(s) rendered this run have no recorded baseline: %s -- remove _perf_baseline.csv and re-run once to recapture, or add rows manually.",
        length(missing), paste(missing, collapse = ", ")
      )
    )
  }
})
