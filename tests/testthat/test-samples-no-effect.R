# pt_figure(samples =) with no track that needs a sample used to be ignored silently; it now
# warns, like comparison = does, but still draws the figure. Uses the bundled mammary example.
test_that("pt_figure(samples =) warns when no track needs a sample, and not otherwise", {
  dir <- tempfile("pt-smp-noeffect-")
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  proj <- suppressWarnings(suppressMessages(pt_project(dir, pt_example())))

  expect_warning(fig <- pt_figure(proj, "Elf5", c("genes", "transcripts"), samples = "LP"),
                 "samples has no effect")
  expect_identical(fig$keys, pt_figure(proj, "Elf5", c("genes", "transcripts"))$keys)
  # One entry per track (NA for tracks without a sample) is just as ineffective.
  expect_warning(pt_figure(proj, "Elf5", c("genes", "transcripts"), samples = c(NA, NA)), "samples has no effect")
  # By track and sample number.
  expect_warning(pt_figure(proj, "Elf5", c(1, 4), samples = 2), "samples has no effect")

  # A track that needs a sample uses it: no warning.
  expect_no_warning(pt_figure(proj, "Elf5", c("atac", "genes"), samples = "LP"))
  expect_no_warning(pt_figure(proj, "Elf5", c("atac", "genes"), samples = c("LP", NA)))
})
