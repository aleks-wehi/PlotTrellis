# Cache-correctness contracts for the caches under the R_user_dir() cache root.

test_that("protein_domains never persists an error-status cache entry", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "no /vast test data on this machine")

  tmp_cache <- tempfile()
  dir.create(tmp_cache)
  on.exit(unlink(tmp_cache, recursive = TRUE), add = TRUE)
  gene <- FIXTURES$mouse$gene

  # A tx_meta_lookup_path that doesn't resolve to a real file forces
  # .pd_gene_transcript_tags()'s readRDS() to fail *after* a successful
  # UniProt search/fetch, landing in gene_protein_domains()'s catch-all
  # (status = "error") rather than a deterministic outcome.
  res1 <- gene_protein_domains(gene, TXDBS$mouse, "/nonexistent/genome_tx_meta.rds",
                                GENOMES$mouse$orgdb, "mm10_error_status_test",
                                cache_dir = tmp_cache)
  # A live UniProt lookup happens before the forced readRDS() failure this
  # test is actually targeting -- if rest.uniprot.org is unreachable from
  # this machine right now, that earlier step fails first with
  # "uniprot_unreachable" instead, which isn't the scenario under test.
  skip_if(res1$status == "uniprot_unreachable", "UniProt unreachable from this machine right now")
  expect_identical(res1$status, "error")
  expect_false(file.exists(file.path(tmp_cache, "mm10_error_status_test", paste0(gene, ".rds"))))

  # Calling again must re-attempt (and re-fail) rather than serve a stale
  # cached "ok" from some earlier, unrelated write.
  res2 <- gene_protein_domains(gene, TXDBS$mouse, "/nonexistent/genome_tx_meta.rds",
                                GENOMES$mouse$orgdb, "mm10_error_status_test",
                                cache_dir = tmp_cache)
  expect_identical(res2$status, "error")
})

test_that("gene_display_transcripts cache key includes the samples list", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "no /vast test data on this machine")

  fx <- FIXTURES$mouse
  anno <- .pt_annotation(test_settings("mouse"))
  cache_path <- tempfile(fileext = ".rds")
  on.exit(unlink(cache_path), add = TRUE)
  samples_ba <- fx$samples[1]  # "Ba" condition only
  samples_lp <- fx$samples[2]  # "LP" condition only

  res_ba <- gene_display_transcripts(fx$gene, anno, samples_ba, fx$data_dir, cache_path = cache_path)
  res_lp <- gene_display_transcripts(fx$gene, anno, samples_lp, fx$data_dir, cache_path = cache_path)

  # Same gene/gtf_path/data_dir, different samples[] -- if the cache key
  # dropped samples, this second call would wrongly return the first call's
  # cached result instead of recomputing against the LP condition's own
  # quant.sf.
  expect_false(identical(res_ba$named, res_lp$named))
})

test_that("gene_display_transcripts cache key includes min_display_tpm", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "no /vast test data on this machine")

  fx <- FIXTURES$mouse
  anno <- .pt_annotation(test_settings("mouse"))
  cache_path <- tempfile(fileext = ".rds")
  on.exit(unlink(cache_path), add = TRUE)
  low <- gene_display_transcripts(fx$gene, anno, fx$samples, fx$data_dir, min_display_tpm = 1, cache_path = cache_path)
  # The first call populated the disk cache; with min_display_tpm missing from
  # the key, this call returned the same (cached) named set.
  high <- gene_display_transcripts(fx$gene, anno, fx$samples, fx$data_dir, min_display_tpm = 1e9, cache_path = cache_path)
  expect_gt(length(low$named), 0)
  expect_length(high$named, 0)
})

test_that("a cache that can't be written costs one warning naming the file, never the result (no fixture data)", {
  blocker <- tempfile()
  writeLines("", blocker)
  on.exit(unlink(blocker), add = TRUE)
  bad <- file.path(blocker, "cache")  # can't be created: its parent is a file
  warnings_of <- function(expr) {
    w <- character(0)
    withCallingHandlers(expr, warning = function(c) {
      w <<- c(w, conditionMessage(c))
      invokeRestart("muffleWarning")
    })
    w
  }

  local_mocked_bindings(gtf_transcript_structure = function(gene, anno)
    data.frame(transcript_id = "T1.1", transcript_name = "G-201", tx_id_stripped = "T1"))
  path <- file.path(bad, "gene_display_transcripts_cache.rds")
  w <- warnings_of(res <- gene_display_transcripts("G", list(gtf_path = "g.gtf"), list(), "", cache_path = path))
  expect_identical(res$named, "G-201")
  expect_length(w, 1)
  expect_match(w, path, fixed = TRUE)

  # A "no reviewed entry" answer is cached by uniprot_canonical_match().
  local_mocked_bindings(.uniprot_read_lines = function(url) character(0))
  path <- file.path(bad, "uniprot_canonical_cache.rds")
  w <- warnings_of(res <- uniprot_canonical_match("G", "ENSMUST00000000001.1", list(), orgdb = "org.Mm.eg.db",
                                                  cache_path = path))
  expect_identical(nrow(res), 1L)
  expect_length(w, 1)
  expect_match(w, path, fixed = TRUE)

  # A resolved protein-domain result is cached.
  local_mocked_bindings(
    .uniprot_read_lines = function(url) {
      if (grepl("/search?", url, fixed = TRUE)) return("Q99999")
      if (grepl("\\.txt$", url)) return(c("FT   DOMAIN          1..10", "FT                   /note=\"Test domain\""))
      c(">sp|Q99999|TEST", strrep("M", 99))
    },
    .pd_gene_transcript_tags = function(...) data.frame(transcript_id = "T1.1", ccdsid = NA_character_))
  cds <- GenomicRanges::GRangesList(T1.1 = GenomicRanges::GRanges("chr1", IRanges::IRanges(101, 400), strand = "+"))
  w <- warnings_of(res <- gene_protein_domains("G", NULL, "unused_tx_meta.rds", "org.Mm.eg.db", "mm10",
                                               cache_dir = file.path(bad, "protein_domains"), cds_by_tx = cds))
  expect_identical(res$status, "ok")
  expect_length(w, 1)
  expect_match(w, file.path(bad, "protein_domains", "mm10", "G.rds"), fixed = TRUE)
})
