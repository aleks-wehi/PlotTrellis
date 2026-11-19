# The Gene info answer, .pt_main_geneinfo() (R/app-screen-api.R), on the mouse fixture
# genome: every transcript of the gene, and the UniProt lookups behind it, cached in the
# project. test-main-api.R covers its fields on the bundled example. UniProt is never
# contacted and no disk cache is read or written: its cache, network reads and the
# protein-domain lookup are replaced in memory.

gi_ctx <- function() {
  proj <- tempfile("pt-gene-info-"); dir.create(proj)
  seed_project_genome(proj)
  json <- file.path(proj, "project_settings.json")
  jsonlite::write_json(list(data_dir = FIXTURES$mouse$data_dir, samples = FIXTURES$mouse$samples,
                            gtf_path = GENOMES$mouse$gtf_path, orgdb = GENOMES$mouse$orgdb,
                            genome_label = GENOMES$mouse$genome_label), json, auto_unbox = TRUE)
  list(proj = proj, ctx = suppressMessages(.pt_composer_context(json, .pt_read_settings(json, .pt_default_settings()))))
}

# In-memory UniProt: the canonical-isoform cache, canned reads (counted: a search hit
# Q99999, then its 99 aa record, which matches no transcript), and a counted protein-domain
# lookup that records the cache folder and CDS index it is given.
local_fake_uniprot <- function(env = parent.frame()) {
  mem <- new.env()
  mem$cache <- list(); mem$reads <- 0L; mem$pd_calls <- 0L
  local_mocked_bindings(
    .uniprot_canonical_cache_load = function(path) {
      mem$uniprot_cache_path <- path
      mem$cache
    },
    .uniprot_canonical_cache_save = function(cache, path) mem$cache <- cache,
    .uniprot_read_lines = function(url) {
      mem$reads <- mem$reads + 1L
      if (grepl("/search?", url, fixed = TRUE)) "Q99999"
      else c("ID   TEST_MOUSE              Reviewed;          99 AA.",
             "SQ   SEQUENCE   99 AA;  11111 MW;  0000000000000000 CRC64;")
    },
    gene_protein_domains = function(gene_symbol, ...) {
      mem$pd_calls <- mem$pd_calls + 1L
      mem$pd_dir <- list(...)$cache_dir
      mem$pd_got_cds <- methods::is(list(...)$cds_by_tx, "GRangesList")  # the app's CDS_BY_TX is passed on
      list(status = "no_domains", gene = gene_symbol)
    },
    .env = env)
  mem
}

test_that("Gene info lists every transcript of the gene", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  local_fake_uniprot()
  g <- gi_ctx()
  on.exit(unlink(g$proj, recursive = TRUE), add = TRUE)
  gm <- GENOMES$mouse
  n_tx <- nrow(transcript_annotations("Nfatc1", TXDBS$mouse, gm$orgdb, gm$tx_name_lookup_path, gm$tx_meta_lookup_path))
  expect_gt(n_tx, 1)
  gi <- .pt_main_geneinfo(g$ctx, "Nfatc1")
  expect_identical(gi$gene, "Nfatc1")
  names_ <- vapply(gi$transcripts, `[[`, "", "name")
  expect_length(names_, n_tx)
  expect_true("Nfatc1-201" %in% names_)
  expect_identical(anyDuplicated(names_), 0L)
})

test_that("the UniProt lookups are cached in the project and read the network once", {
  skip_if_not(PT_TEST_MOUSE_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  mem <- local_fake_uniprot()
  g <- gi_ctx()
  on.exit(unlink(g$proj, recursive = TRUE), add = TRUE)
  cold <- .pt_main_geneinfo(g$ctx, "Elf5")
  expect_equal(mem$reads, 2L)  # one search + one record, then cached
  expect_true("org.Mm.eg.db|Elf5" %in% names(mem$cache))
  expect_equal(mem$pd_calls, 1L)
  expect_true(mem$pd_got_cds)
  # Both caches are the project's own.
  expect_identical(dirname(mem$uniprot_cache_path), file.path(g$proj, "cache"))
  expect_identical(mem$pd_dir, file.path(g$proj, "cache", "protein_domains"))
  # The 99 aa record matches no transcript: no UniProt tag; no features, said as such.
  expect_true("Elf5-201" %in% vapply(cold$transcripts, `[[`, "", "name"))
  expect_false(any(grepl("UniProt", unlist(lapply(cold$transcripts, `[[`, "tags")), fixed = TRUE)))
  expect_identical(cold$protein_status, .pd_status_text("no_domains", panel = TRUE))
  expect_match(cold$protein_status, "No annotated domain/region/motif features", fixed = TRUE)

  # Warm: the same answer, with no new network read.
  expect_identical(.pt_main_geneinfo(g$ctx, "Elf5"), cold)
  expect_equal(mem$reads, 2L)
})
