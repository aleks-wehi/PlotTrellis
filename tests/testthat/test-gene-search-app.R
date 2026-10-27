# Phase I.1 step 4a review: the real app server's gene observer. testServer()
# doesn't echo updateSelectizeInput() back, so each echo the browser would
# send after a resolved ID or case fix is set by hand.

test_that("a typed Ensembl ID keeps its own region through the symbol echo", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  proj <- tempfile("pt-gene-app-"); dir.create(proj)
  on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  json <- file.path(proj, "project_settings.json")
  jsonlite::write_json(list(data_dir = FIXTURES$mouse$data_dir, samples = FIXTURES$mouse$samples,
                            gtf_path = GENOMES$mouse$gtf_path, orgdb = GENOMES$mouse$orgdb,
                            genome_label = GENOMES$mouse$genome_label), json, auto_unbox = TRUE)
  app <- .pt_composer_app(json, .pt_read_settings(json, .pt_default_settings()))
  shiny::testServer(app, {
    # Amy2a1 has 3 Ensembl IDs; ENSMUSG00000093931 is also named Amy2a3
    want <- GenomicFeatures::genes(TXDBS$mouse, filter = list(gene_id = "ENSMUSG00000093931"))
    session$setInputs(gene_name = "ENSMUSG00000093931")
    by_id <- state$gene_region
    expect_equal(state$gene_symbol, "Amy2a1")
    expect_true(by_id$start < GenomicRanges::start(want) && by_id$end > GenomicRanges::end(want))
    expect_match(state$gene_note, "also named Amy2a3")
    session$setInputs(gene_name = "Amy2a1")  # the echo
    expect_identical(state$gene_region, by_id)
    session$setInputs(gene_name = "pax5")
    expect_equal(state$gene_symbol, "Pax5")
    expect_null(state$gene_note)
    session$setInputs(gene_name = "Pax5")    # the echo
    expect_equal(state$gene_symbol, "Pax5")
    session$setInputs(gene_name = "")
    expect_null(state$gene_symbol)
  })
})
