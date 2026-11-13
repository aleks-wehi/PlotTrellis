# AnnotationDbi prints a "'select()' returned 1:1 mapping between keys and columns"
# note on every select() call, several times per figure. Every select() call is
# wrapped in suppressMessages().

test_that("rendering gene and transcript tracks prints no select() mapping notes (mouse dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  rows <- test_rows("mouse")
  specs <- list(single_track("genes", rows[["Ba"]]), single_track("transcripts", rows[["Ba"]]),
                single_track("transcripts_annotated", rows[["Ba"]]), single_track("tpm_usage", rows[["Ba"]]))
  f <- tempfile(fileext = ".svg")
  on.exit(unlink(f), add = TRUE)
  # render_all() unwrapped on purpose (it prints plotgardener's grob lines here): .pt_quiet()
  # also muffles select() notes as a backstop, so wrapping this call would pass even without
  # the per-call suppressMessages().
  expect_no_message(
    render_all(specs, test_region("mouse"), test_settings("mouse"),
               device_fn = function(w, h) svglite::svglite(f, width = w, height = h)),
    message = "select\\(\\)' returned")
  grDevices::dev.off()
})

test_that("gene search and the Gene info transcript table print no select() mapping notes (mouse dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  fx <- FIXTURES[["mouse"]]
  genome <- GENOMES[[fx$genome]]
  expect_no_message(gene_to_coords(fx$gene, TXDBS[[fx$genome]], orgdb_pkg = genome$orgdb),
                    message = "select\\(\\)' returned")
  expect_no_message(transcript_annotations(fx$gene, TXDBS[[fx$genome]], genome$orgdb,
                                           genome$tx_name_lookup_path, genome$tx_meta_lookup_path),
                    message = "select\\(\\)' returned")
  # The protein track's transcript lookup (offline; the UniProt fetch is a separate step).
  expect_no_message(.pd_gene_transcript_tags(fx$gene, genome$tx_meta_lookup_path, TXDBS[[fx$genome]], genome$orgdb),
                    message = "select\\(\\)' returned")
})

test_that(".pt_quiet() hides plotgardener grob lines and select() notes, nothing else", {
  # The app's figure render and pt_save() both go through it.
  shown <- character()
  build_msg <- .pt_txdb_build_message(tempfile("a", fileext = ".gtf"))  # the real first-build line
  withCallingHandlers(.pt_quiet({
    message("text[text2]"); message("segments[segments2]"); message("signal[signal1_h]")
    message("'select()' returned 1:1 mapping between keys and columns")
    message("PlotTrellis: genes track failed: BOOM genes [x]")
    message("PlotTrellis: signal track failed: x[i]")  # a bracket pair inside a message
    message("genes[genes1] extra")                     # a grob-like start, then more text
    message(build_msg)
  }), message = function(m) {
    shown <<- c(shown, trimws(conditionMessage(m)))
    invokeRestart("muffleMessage")
  })
  expect_identical(shown, c("PlotTrellis: genes track failed: BOOM genes [x]",
                            "PlotTrellis: signal track failed: x[i]", "genes[genes1] extra",
                            build_msg))
})

test_that("the first-build line gives a time estimate only for a large GTF", {
  small <- tempfile(fileext = ".gtf"); writeLines("x", small)
  big   <- tempfile(fileext = ".gtf"); writeBin(raw(1.1e7), big)
  on.exit(unlink(c(small, big)))
  expect_identical(.pt_txdb_build_message(small),
                   paste0("Building TxDb from ", small, " (this only happens once per GTF file)..."))
  expect_identical(.pt_txdb_build_message(big),
                   paste0("Building TxDb from ", big,
                          " (this only happens once per GTF file, and takes about 1-3 minutes, longer on a slow disk)..."))
  expect_identical(.pt_txdb_build_message("/no/such.gtf"),
                   "Building TxDb from /no/such.gtf (this only happens once per GTF file)...")
})
