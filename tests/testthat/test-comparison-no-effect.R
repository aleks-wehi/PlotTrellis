# pt_figure(comparison =) with no track that needs a sample used to be ignored silently; it
# now warns, but still draws the figure. Uses the bundled mammary example (no fixture data).
test_that("pt_figure(comparison =) warns when no track needs a sample, and not otherwise", {
  dir <- tempfile("pt-cmp-noeffect-")
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  proj <- suppressWarnings(suppressMessages(pt_project(dir, pt_example())))
  samples <- pt_samples(proj)$name

  expect_warning(fig <- pt_figure(proj, "Eya2", c("genes", "transcripts"), comparison = rev(samples)),
                 "comparison has no effect")
  # Same figure as without it.
  expect_identical(fig$keys, pt_figure(proj, "Eya2", c("genes", "transcripts"))$keys)
  # By track number, as in the bug report.
  expect_warning(pt_figure(proj, "Eya2", c(1, 2), comparison = c(2, 1)), "comparison has no effect")

  # A track that needs a sample uses it: no warning.
  expect_no_warning(pt_figure(proj, "Eya2", c("ATAC", "genes"), comparison = samples))
})
