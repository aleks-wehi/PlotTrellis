# Phase G.1 step 4 verify: two cache-correctness contracts that must survive
# the R_user_dir() path migration unchanged (see phase-G.1.md step 4).

test_that("protein_domains never persists an error-status cache entry", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "no /vast test data on this machine")

  tmp_cache <- tempfile()
  dir.create(tmp_cache)
  on.exit(unlink(tmp_cache, recursive = TRUE), add = TRUE)
  gene <- FIXTURES$mouse$gene

  # A gtf_path that doesn't resolve to a real tx_meta lookup forces
  # .pd_gene_transcript_tags()'s readRDS() to fail *after* a successful
  # UniProt search/fetch, landing in gene_protein_domains()'s catch-all
  # (status = "error") rather than a deterministic outcome.
  res1 <- gene_protein_domains(gene, TXDBS$mouse, "/nonexistent/genome.gtf.gz",
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
  res2 <- gene_protein_domains(gene, TXDBS$mouse, "/nonexistent/genome.gtf.gz",
                                GENOMES$mouse$orgdb, "mm10_error_status_test",
                                cache_dir = tmp_cache)
  expect_identical(res2$status, "error")
})

test_that("gene_display_transcripts cache key includes the samples list", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "no /vast test data on this machine")

  fx <- FIXTURES$mouse
  gm <- GENOMES$mouse
  .set_transcript_filter_globals(test_settings("mouse"))
  samples_ba <- fx$samples[1]  # "Ba" condition only
  samples_lp <- fx$samples[2]  # "LP" condition only

  res_ba <- gene_display_transcripts(fx$gene, gm$gtf_path, samples_ba, fx$data_dir)
  res_lp <- gene_display_transcripts(fx$gene, gm$gtf_path, samples_lp, fx$data_dir)

  # Same gene/gtf_path/data_dir, different samples[] -- if the cache key
  # dropped samples, this second call would wrongly return the first call's
  # cached result instead of recomputing against the LP condition's own
  # quant.sf.
  expect_false(identical(res_ba$named, res_lp$named))
})

test_that("gene_display_transcripts cache key includes MIN_DISPLAY_TPM (bug 61)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "no /vast test data on this machine")

  fx <- FIXTURES$mouse
  gm <- GENOMES$mouse
  s <- test_settings("mouse")
  on.exit(.set_transcript_filter_globals(s), add = TRUE)
  .set_transcript_filter_globals(utils::modifyList(s, list(min_display_tpm = 1)))
  low <- gene_display_transcripts(fx$gene, gm$gtf_path, fx$samples, fx$data_dir)
  # The first call populated the disk cache; with MIN_DISPLAY_TPM missing from
  # the key, this call returned the same (cached) named set.
  .set_transcript_filter_globals(utils::modifyList(s, list(min_display_tpm = 1e9)))
  high <- gene_display_transcripts(fx$gene, gm$gtf_path, fx$samples, fx$data_dir)
  expect_gt(length(low$named), 0)
  expect_length(high$named, 0)
})
