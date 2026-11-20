# The app's gene lookup on the mouse fixture genome: .pt_main_locus() (R/app-screen-api.R),
# directly and through the server's "locus" request with the live Gene flank % (pattern of
# test-main-server.R). The resolver itself is test-gene-search.R's.

# A project on the mouse fixture with the given extra settings, and its startup context.
gene_ctx <- function(extra = list(), proj = tempfile("pt-gene-app-")) {
  dir.create(proj)
  seed_project_genome(proj)
  json <- file.path(proj, "project_settings.json")
  s <- c(list(data_dir = FIXTURES$mouse$data_dir, samples = FIXTURES$mouse$samples,
              gtf_path = GENOMES$mouse$gtf_path, orgdb = GENOMES$mouse$orgdb,
              genome_label = GENOMES$mouse$genome_label), extra)
  writeLines(jsonlite::toJSON(s, auto_unbox = TRUE, null = "null", na = "null"), json)
  cfg <- .pt_read_settings(json, .pt_default_settings())
  list(proj = proj, json = json, ctx = suppressMessages(.pt_composer_context(json, cfg)))
}

test_that("a typed Ensembl ID opens its own gene under its first symbol", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  g <- gene_ctx()
  on.exit(unlink(g$proj, recursive = TRUE), add = TRUE)
  # ENSMUSG00000093931's symbols depend on the org.Mm.eg.db build (Amy2a1 and
  # Amy2a3 on Milton's build, Amy2a1 having 3 IDs), so they are read from its own tables, as
  # in test-gene-search.R.
  map <- merge(AnnotationDbi::toTable(org.Mm.eg.db::org.Mm.egENSEMBL),
               AnnotationDbi::toTable(org.Mm.eg.db::org.Mm.egSYMBOL), by = "gene_id")
  syms <- sort(unique(map$symbol[map$ensembl_id == "ENSMUSG00000093931"]), method = "radix")
  skip_if(length(syms) == 0, "this org.Mm.eg.db build no longer knows ENSMUSG00000093931")
  want <- GenomicFeatures::genes(TXDBS$mouse, filter = list(gene_id = "ENSMUSG00000093931"))
  by_id <- .pt_main_locus(g$ctx, "ENSMUSG00000093931")
  expect_identical(by_id$gene, syms[1])
  expect_identical(c(by_id$s, by_id$e), c(GenomicRanges::start(want), GenomicRanges::end(want)))
  expect_true(by_id$view$s < GenomicRanges::start(want) && by_id$view$e > GenomicRanges::end(want))
  # Any case: the symbol as the OrgDb spells it.
  expect_identical(.pt_main_locus(g$ctx, "pax5")$gene, "Pax5")
  expect_identical(.pt_main_locus(g$ctx, "Pax5")$gene, "Pax5")
})

test_that("an empty or non-numeric Gene flank % never ends the session: the default is used (bug 193)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  default_view <- NULL
  for (saved in list(15, NULL, NA, "abc")) {  # "abc": a hand-edited project file
    extra <- list()
    extra["gene_flank_pct"] <- list(saved)  # a project saved with the box empty stores null
    g <- suppressWarnings(gene_ctx(extra))  # "abc": an unusable-value warning
    rec <- recording_session()
    view_of <- function(id) {
      rec$session$setInputs(pt_rpc = as.character(jsonlite::toJSON(
        list(list(id = id, op = "locus", args = list(q = "Elf5"))), auto_unbox = TRUE, digits = NA)))
      r <- Filter(function(m) identical(m$type, "pt_rpc_reply") && identical(m$message$id, id), rec$custom)
      expect_length(r, 1)
      reply <- jsonlite::fromJSON(r[[1]]$message$json, simplifyVector = FALSE)
      expect_null(reply$error)
      reply$view
    }
    shiny::testServer(suppressMessages(.pt_main_app(g$ctx, g$json)), {
      session$flushReact()
      ok <- view_of("saved")
      expect_false(is.null(ok), info = deparse(saved))
      if (is.null(default_view)) default_view <<- ok
      expect_identical(ok, default_view, info = deparse(saved))  # the default 15%
      for (typed in list(NA, NA_real_, "abc")) {  # NA is what an empty box sends
        session$setInputs(s_gene_flank_pct = typed)
        expect_identical(view_of(paste("typed", deparse(typed))), default_view,
                         info = paste(deparse(saved), deparse(typed)))
      }
    }, session = rec$session)
    unlink(g$proj, recursive = TRUE)
  }
})
