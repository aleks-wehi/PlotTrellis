# .pt_no_gene(): the one test of "this view has no gene" every gene track and the Gene info
# panel use.

test_that(".pt_no_gene() is TRUE for no value, several values, NA and a blank string", {
  expect_true(.pt_no_gene(NULL))
  expect_true(.pt_no_gene(character(0)))
  expect_true(.pt_no_gene(NA))
  expect_true(.pt_no_gene(NA_character_))
  expect_true(.pt_no_gene(""))
  expect_true(.pt_no_gene("  "))
  expect_true(.pt_no_gene(c("a", "b")))
  expect_false(.pt_no_gene("Elf5"))
})

test_that("a blank gene draws what no gene draws, on every track that needs one (human dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  rows <- test_rows("human")
  entries <- list(single_track("transcripts_all", rows[["GM12878"]]), single_track("transcripts", rows[["GM12878"]]),
                  single_track("transcripts_annotated", rows[["GM12878"]]),
                  single_track("tpm_usage", rows[["GM12878"]]), single_track("tpm_usage", rows[["K562"]]),
                  single_track("protein_domains", rows[["GM12878"]]),
                  pair_track("junctions_mirrored_diff_tpm_scaled", rows[["GM12878"]], rows[["K562"]]))
  svg_for <- function(gene) {
    region <- utils::modifyList(test_region("human"), list(gene = gene))
    f <- tempfile(fileext = ".svg")
    suppressWarnings(suppressMessages(render_all(entries, region, test_settings("human"),
                                                 device_fn = function(w, h) svglite::svglite(f, width = w, height = h))))
    grDevices::dev.off()
    readLines(f, warn = FALSE)
  }
  none <- svg_for(NA_character_)
  expect_false(any(grepl("track failed", none, fixed = TRUE)))
  expect_identical(svg_for(""), none)
  expect_identical(svg_for("  "), none)
})
