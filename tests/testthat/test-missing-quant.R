# One sample's missing quant.sf must not fail the project's transcript tracks: the named
# transcripts come from the readable samples, only that sample's own expression row shows the
# message naming its file, and a run with an unreadable file is not cached.

without_quant <- function(sample) {
  sample$data <- lapply(sample$data, function(d) {
    if (identical(d$type, "quantdir")) d$path <- "quant/not_there"
    d
  })
  sample
}

test_that("gene_display_transcripts() pools only the readable samples and does not cache the run (mouse)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  fx <- FIXTURES$mouse
  anno <- .pt_annotation(test_settings("mouse"))
  cache_path <- tempfile(fileext = ".rds")
  on.exit(unlink(cache_path), add = TRUE)

  one_gone <- fx$samples
  one_gone[[2]] <- without_quant(one_gone[[2]])
  expect_message(res <- gene_display_transcripts(fx$gene, anno, one_gone, fx$data_dir, cache_path = cache_path),
                 "named transcripts leave out LP \\(file not found: quant.sf")
  from_first <- gene_display_transcripts(fx$gene, anno, fx$samples[1], fx$data_dir,
                                         cache_path = tempfile(fileext = ".rds"))
  expect_identical(res$named, from_first$named)
  expect_false(file.exists(cache_path))

  # Every quant file missing: the ordinal fallback, still no error and still not cached.
  all_gone <- lapply(fx$samples, without_quant)
  res <- suppressMessages(gene_display_transcripts(fx$gene, anno, all_gone, fx$data_dir, cache_path = cache_path))
  expect_gt(length(res$named), 0)
  expect_false(file.exists(cache_path))
})

test_that("a missing quant.sf fails only that sample's expression row; the transcript tracks still draw", {
  # The copy sits under a long folder, so the drawn message is cut, as on a deep project path.
  # Three levels: long enough to be cut, short enough for Windows' 260-character paths.
  top <- tempfile("pt-long-")
  proj_dir <- tempfile("pt-proj-")
  copy <- dirname(example_copy("mammary"))
  on.exit(unlink(c(proj_dir, top, copy), recursive = TRUE), add = TRUE)
  long <- do.call(file.path, as.list(c(top, rep("a-long-data-folder-name", 3))))
  dir.create(long, recursive = TRUE)
  stopifnot(file.rename(copy, file.path(long, "example")))
  txt <- file.path(long, "example", "samples.txt")
  local_offline_chrom_sizes()
  unlink(file.path(dirname(txt), "quant", "LP_rep1", "quant.sf"))
  proj <- suppressMessages(suppressWarnings(pt_project(proj_dir, txt)))
  fig <- pt_figure(proj, "Elf5", c("transcripts", "transcript_expression"))
  out <- tempfile(fileext = ".svg")
  suppressMessages(pt_save(fig, out))
  svg <- paste(readLines(out, warn = FALSE), collapse = "\n")
  # One failure in the whole figure (the LP row), naming the file; the Basal row and the
  # transcripts track drew. Before the fix the transcripts track failed too.
  expect_identical(lengths(regmatches(svg, gregexpr("failed", svg, fixed = TRUE))), 1L)
  expect_match(svg, "quant.sf", fixed = TRUE)
  # Only LP_rep1's quant.sf was deleted, so the one failed row is LP's, and the drawn message
  # names its folder even when the data folder's path is long.
  expect_match(svg, "LP_rep1", fixed = TRUE)
  expect_match(svg, "Transcripts (filtered)", fixed = TRUE)
  expect_match(svg, "Elf5-207", fixed = TRUE)   # named from Basal, the readable sample
})
