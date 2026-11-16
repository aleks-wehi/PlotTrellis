# The real app server's gene observer. testServer()
# doesn't echo updateSelectizeInput() back, so each echo the browser would
# send after a resolved ID or case fix is set by hand.

test_that("a typed Ensembl ID keeps its own region through the symbol echo", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  proj <- tempfile("pt-gene-app-"); dir.create(proj)
  seed_project_genome(proj)
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

test_that("clearing the gene box clears the region, so Generate can't draw the old gene (bug 101)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  proj <- tempfile("pt-gene-app-"); dir.create(proj)
  seed_project_genome(proj)
  on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  json <- file.path(proj, "project_settings.json")
  jsonlite::write_json(list(data_dir = FIXTURES$mouse$data_dir, samples = FIXTURES$mouse$samples,
                            gtf_path = GENOMES$mouse$gtf_path, orgdb = GENOMES$mouse$orgdb,
                            genome_label = GENOMES$mouse$genome_label), json, auto_unbox = TRUE)
  app <- .pt_composer_app(json, .pt_read_settings(json, .pt_default_settings()))
  drawn <- 0L
  local_mocked_bindings(render_all = function(...) { drawn <<- drawn + 1L; invisible(NULL) })
  shiny::testServer(app, {
    session$setInputs(view_mode = "compose", region_mode = "gene", gene_name = "Elf5")
    expect_false(is.null(state$gene_region))
    session$setInputs(gene_name = "Nofoundgene")
    expect_false(is.null(state$gene_error))
    session$setInputs(gene_name = "")
    expect_null(state$gene_error)  # a "not found" for text that is gone goes too
    session$setInputs(gene_name = "Elf5")
    expect_false(is.null(state$gene_region))
    session$setInputs(gene_name = "")
    expect_null(state$gene_region)
    expect_null(state$gene_symbol)
    expect_null(state$gene_error)
    expect_null(resolved_region())
    session$setInputs(figure_tracks = "genes", generate = 1)
    expect_identical(drawn, 0L)
  })
})

test_that("an empty or non-numeric Gene flank % never ends the session: the default is used (bug 193)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  default_region <- NULL
  for (saved in list(15, NULL, NA, "abc")) {  # "abc": a hand-edited project file
    proj <- tempfile("pt-flank-"); dir.create(proj)
    seed_project_genome(proj)
    json <- file.path(proj, "project_settings.json")
    s <- list(data_dir = FIXTURES$mouse$data_dir, samples = FIXTURES$mouse$samples,
              gtf_path = GENOMES$mouse$gtf_path, orgdb = GENOMES$mouse$orgdb,
              genome_label = GENOMES$mouse$genome_label)
    s["gene_flank_pct"] <- list(saved)  # a project saved with the box empty stores null
    writeLines(jsonlite::toJSON(s, auto_unbox = TRUE, null = "null", na = "null"), json)
    app <- suppressWarnings(suppressMessages(.pt_composer_app(json, .pt_read_settings(json, .pt_default_settings()))))  # "abc": an unusable-value warning
    shiny::testServer(app, {
      session$setInputs(view_mode = "compose", region_mode = "gene", gene_name = "Elf5")
      ok <- state$gene_region
      expect_false(is.null(ok), info = deparse(saved))
      if (is.null(default_region)) default_region <<- ok
      expect_identical(ok, default_region, info = deparse(saved))  # the default 15%
      for (typed in list(NA, NA_real_, "abc")) {  # NA is what an empty box sends
        session$setInputs(s_gene_flank_pct = typed)
        expect_identical(state$gene_region, default_region, info = paste(deparse(saved), deparse(typed)))
      }
    })
    unlink(proj, recursive = TRUE)
  }
})
