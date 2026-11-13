# Golden-image tests for a track with nothing to draw, against the fixtures: an empty view
# draws the track as usual with nothing in it (the hic and junctions cases), and a missing file
# fails the track in red (the bw case, read off the SVG text: the message holds a path).
# test-empty-view.R covers every track type on synthetic files.
#
# The gene-less case uses transcripts_all as the representative of the
# modules that need a gene (transcripts, transcripts_annotated,
# protein_domains, tpm_usage and junctions_mirrored_diff_tpm_scaled share the
# same no-gene check, .pt_no_gene()); they aren't each given their own
# golden image here, to avoid five near-duplicate snapshots of the same
# placeholder text.
#
# One case per category, each a real crash-turned-render against the
# fixtures -- not a hypothetical or a full module x dataset grid.

test_that("out-of-bounds region: hic draws its empty track instead of crashing (human)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  rows_h    <- test_rows("human")
  settings_h <- test_settings("human")
  # Genuinely past the end of any real chr1 (hg38 chr1 is ~249Mb) -- not a
  # sparse-but-in-range locus like the PTPRC case below,
  # a coordinate that doesn't exist on the chromosome at all.
  oob_region <- list(chrom = "chr1", start = 990000000, end = 990010000, gene = "PTPRC")
  spec <- single_track("hic", rows_h[["GM12878"]])
  fig  <- test_fig(list(spec), oob_region, settings_h)
  # plotgardener's own "No data found in region." warning is expected for an empty view.
  withCallingHandlers(vdiffr::expect_doppelganger("hic-oob-human-gm12878", fig),
    warning = function(w) if (grepl("No data found in region", conditionMessage(w), fixed = TRUE)) invokeRestart("muffleWarning"))
})

test_that("sparse in-range region: hic draws its empty track with no -Inf warning (human)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  # The PTPRC window on the GM12878 ChIA-PET .hic is on the chromosome and
  # in range, but auto resolution picks 1000 bp and `straw` returns a few rows whose counts are all
  # NaN, so there is no colour palette for the legend. Not the 0-row out-of-bounds case above.
  region <- list(chrom = "chr1", start = 143800000, end = 144200000, gene = "PTPRC")
  spec   <- single_track("hic", test_rows("human")[["GM12878"]])
  # Guard: still rows, all without a usable count (else this would silently become the 0-row case).
  m <- load_hic(spec, region, test_settings("human"))$data$matrix
  expect_gt(nrow(m), 0)
  expect_true(all(is.na(m$counts)))
  tmp <- tempfile(fileext = ".svg")
  on.exit(unlink(tmp), add = TRUE)
  # plotgardener's own "No data found in region." warning is expected here and muffled; the point is that
  # settings_hic() adds no "returning -Inf" warning of its own for an all-NaN matrix.
  expect_no_warning(
    withCallingHandlers(
      render_all(list(spec), region, test_settings("human"),
                 device_fn = function(w, h) svglite::svglite(tmp, width = w, height = h)),
      warning = function(w) if (grepl("No data found in region", conditionMessage(w), fixed = TRUE)) invokeRestart("muffleWarning")))
  grDevices::dev.off()
  svg <- paste(readLines(tmp, warn = FALSE), collapse = "\n")
  # The heading as usual, no colour legend (no counts anywhere), no message.
  expect_match(svg, ">hic</text>", fixed = TRUE)
  expect_false(grepl("No data", svg, fixed = TRUE))
  expect_false(grepl("track failed", svg, fixed = TRUE))
})

test_that("empty GRanges: junctions draws its empty track instead of crashing (mouse)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  rows_m    <- test_rows("mouse")
  region_m  <- test_region("mouse")
  settings_m <- test_settings("mouse")
  # Mouse's "Junctions" slot is repurposed OmniC/FitHiChIP loop data (BEDPE) --
  # the Ba/Elf5 window has zero overlapping loop calls, so pool_junctions()
  # returns a zero-row data.frame.
  spec <- single_track("junctions", rows_m[["Ba"]])
  fig  <- test_fig(list(spec), region_m, settings_m)
  vdiffr::expect_doppelganger("junctions-empty-mouse-ba", fig)
})

test_that("missing file: bw fails its track in red, the figure still draws (mouse)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  rows_m    <- test_rows("mouse")
  region_m  <- test_region("mouse")
  settings_m <- test_settings("mouse")
  # Simulates a project_settings.json referencing a bigWig that's since been
  # moved/deleted -- nothing stops render_all() itself from being handed a
  # bad path. `path` is the real entry field name (viz-bw.R).
  spec <- single_track("bw", rows_m[["Ba"]])
  spec$path <- file.path(FIXTURES$mouse$data_dir, "atac", "DOES_NOT_EXIST.bw")
  tmp <- tempfile(fileext = ".svg")
  on.exit(unlink(tmp), add = TRUE)
  expect_message(render_all(list(spec, single_track("genes", rows_m[["Ba"]])), region_m, settings_m,
                            device_fn = function(w, h) svglite::svglite(tmp, width = w, height = h)),
                 "bw \\(Ba\\) track failed: file not found: DOES_NOT_EXIST\\.bw")
  grDevices::dev.off()
  svg <- paste(readLines(tmp, warn = FALSE), collapse = "\n")
  expect_match(svg, "fill: #CC0000;[^>]*>bw \\(Ba\\) track failed: file not found: DOES_NOT_EXIST\\.bw")
  expect_match(svg, ">Elf5</text>", fixed = TRUE)  # the gene track still drew
})

test_that("gene-less region: transcripts_all renders its message instead of crashing (coords mode with no gene, human)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  rows_h     <- test_rows("human")
  settings_h <- test_settings("human")
  # Coords mode with no gene resolved -- exactly the state the app's
  # red-box guard on gene-only tracks exists to prevent a
  # biologist from reaching in the UI; this test bypasses that guard
  # directly to exercise the module's own defense-in-depth path.
  no_gene_region <- list(chrom = "chr1", start = 1000000, end = 1010000, gene = NA_character_)
  spec <- single_track("transcripts_all", rows_h[["GM12878"]])
  fig  <- test_fig(list(spec), no_gene_region, settings_h)
  vdiffr::expect_doppelganger("transcripts-all-gene-less-human-gm12878", fig)
})

test_that("a chromosome the .hic lists but holds no contacts on: hic draws its empty track, not a failure (mouse)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  spec <- single_track("hic", test_rows("mouse")[["Ba"]])
  region <- list(chrom = "chrM", start = 1L, end = 16000L, gene = NA_character_)
  # Guard: straw itself fails on this chromosome (no chrM_chrM block), the case under test.
  expect_error(strawr::straw("NONE", spec$path, "chrM:1:16000", "chrM:1:16000", "BP", load_hic(spec, region, test_settings("mouse"))$data$resolution),
               "chr_chr map")
  tmp <- tempfile(fileext = ".svg")
  on.exit(unlink(tmp), add = TRUE)
  suppressWarnings(render_all(list(spec), region, test_settings("mouse"),
                              device_fn = function(w, h) svglite::svglite(tmp, width = w, height = h)))
  grDevices::dev.off()
  svg <- paste(readLines(tmp, warn = FALSE), collapse = "\n")
  expect_match(svg, ">hic</text>", fixed = TRUE)
  expect_false(grepl("track failed|No data|Chromosome not found", svg))
})
