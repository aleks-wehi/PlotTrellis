# The UniProt-CCDS cross-reference must match whether the transcript's CCDS id
# carries a version (GENCODE "CCDS12345.1") or not (Ensembl "CCDS12345"); UniProt's
# own CCDS ids are versioned.
test_that("uniprot_canonical_match() matches a CCDS id with or without a version suffix", {
  cache_path <- tempfile(fileext = ".rds")
  on.exit(unlink(cache_path), add = TRUE)
  entry <- list(accession = "Q99999", protein_len = 999,
                ensembl = data.frame(id = character(0), isoform = character(0)),
                ccds = data.frame(id = "CCDS12345.1", isoform = NA_character_),
                isoforms = data.frame(number = character(0), isoid = character(0), name = character(0)),
                canonical_isoform = NA_character_,
                canonical_isoform_id = "Q99999", canonical_isoform_name = NA_character_)
  cache <- list("org.Mm.eg.db|Testgene" = entry)
  attr(cache, "schema_version") <- UNIPROT_CACHE_SCHEMA_VERSION
  saveRDS(cache, cache_path)

  tx <- c("T1.1", "T2.1", "T3.1")
  res <- uniprot_canonical_match("Testgene", tx, list(), tx_ccds_id = c("CCDS12345", "CCDS12345.1", "CCDS99999"),
                                 orgdb = "org.Mm.eg.db", cache_path = cache_path)
  expect_identical(res$uniprot_source, c("ccds_xref", "ccds_xref", "length_match"))
  expect_identical(res$uniprot_match, c(TRUE, TRUE, FALSE))
})
