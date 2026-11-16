# Deterministic protein-domain outcomes are cached per gene like a resolved one: no reviewed
# UniProt entry, no domain features, no transcript match. Only network failures and errors are
# looked up again on the next render.

pd_status_lookup <- function(cache_dir, gene = "G", cds = NULL)
  gene_protein_domains(gene, NULL, "unused_tx_meta.rds", "org.Mm.eg.db", "mm10",
                       cache_dir = cache_dir, cds_by_tx = cds)

# After `first` is cached, a second lookup must not touch UniProt at all.
expect_cached <- function(cache_dir, gene, status) {
  local_mocked_bindings(.uniprot_read_lines = function(url) stop("looked up again"))
  again <- pd_status_lookup(cache_dir, gene)
  expect_identical(again$status, status)
  expect_true(file.exists(file.path(cache_dir, "mm10", paste0(gene, ".rds"))))
}

test_that("no reviewed UniProt entry is cached", {
  cache_dir <- tempfile("pd-")
  local_mocked_bindings(.uniprot_read_lines = function(url) character(0))
  expect_identical(pd_status_lookup(cache_dir)$status, "no_uniprot_entry")
  expect_cached(cache_dir, "G", "no_uniprot_entry")
})

test_that("a record without domain features is cached", {
  cache_dir <- tempfile("pd-")
  local_mocked_bindings(.uniprot_read_lines = function(url)
    if (grepl("/search?", url, fixed = TRUE)) "Q99999" else c("ID   TEST_MOUSE   Reviewed;   99 AA.", "SQ   SEQUENCE   99 AA;"))
  expect_identical(pd_status_lookup(cache_dir)$status, "no_domains")
  expect_cached(cache_dir, "G", "no_domains")
})

test_that("a canonical isoform matching no transcript is cached", {
  cache_dir <- tempfile("pd-")
  local_mocked_bindings(
    .uniprot_read_lines = function(url) {
      if (grepl("/search?", url, fixed = TRUE)) return("Q99999")
      if (grepl("\\.txt$", url)) return(c("FT   DOMAIN          1..10", "FT                   /note=\"Test domain\""))
      c(">sp|Q99999|TEST", strrep("M", 99))
    },
    .pd_resolve_canonical_transcript = function(...) list(transcript_id = NA_character_))
  cds <- GenomicRanges::GRangesList(T1.1 = GenomicRanges::GRanges("chr1", IRanges::IRanges(101, 400), strand = "+"))
  expect_identical(pd_status_lookup(cache_dir, cds = cds)$status, "no_canonical_transcript")
  expect_cached(cache_dir, "G", "no_canonical_transcript")
})

test_that("a local error is still not cached, so the next render retries", {
  cache_dir <- tempfile("pd-")
  local_mocked_bindings(
    .uniprot_read_lines = function(url) if (grepl("/search?", url, fixed = TRUE)) "Q99999" else stop("disk full"),
    .pd_fetch_uniprot_record = function(acc) stop("cannot open the record file"))
  res <- pd_status_lookup(cache_dir)
  expect_identical(res$status, "error")
  expect_false(file.exists(file.path(cache_dir, "mm10", "G.rds")))
})
