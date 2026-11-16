# build_tx_meta_lookup() must read the CCDS id from GENCODE's `ccdsid` and from
# Ensembl's `ccds_id` attribute (bug 178); a transcript with neither gets NA.
test_that("build_tx_meta_lookup() reads ccds_id from GENCODE (ccdsid) and Ensembl (ccds_id) GTF rows", {
  row <- function(feature, attrs) paste("chr1", "src", feature, 1, 100, ".", "+", ".", attrs, sep = "\t")
  gtf <- tempfile(fileext = ".gtf.gz")
  con <- gzfile(gtf, "w")
  writeLines(c(
    row("gene", 'gene_id "G1"; gene_name "Elf5";'),
    row("transcript", 'gene_id "G1"; transcript_id "T_GENCODE.1"; transcript_type "protein_coding"; tag "basic"; ccdsid "CCDS11111.1";'),
    row("transcript", 'gene_id "G1"; gene_version "1"; transcript_id "T_ENSEMBL"; transcript_version "1"; transcript_biotype "protein_coding"; tag "basic"; ccds_id "CCDS37873"; tag "Ensembl_canonical";'),
    row("transcript", 'gene_id "G1"; transcript_id "T_NONE"; transcript_biotype "lncRNA"; tag "basic";')
  ), con)
  close(con)
  out <- tempfile(fileext = ".rds")
  meta <- build_tx_meta_lookup(gtf, out)
  expect_equal(meta["T_GENCODE", "ccds_id"], "CCDS11111.1")
  expect_equal(meta["T_ENSEMBL", "ccds_id"], "CCDS37873")
  expect_true(is.na(meta["T_NONE", "ccds_id"]))
  expect_equal(readRDS(out)["T_ENSEMBL", "ccds_id"], "CCDS37873")
})
