# Bug 44: the Gene info table was always empty because
# transcript_annotations() read the OrgDb and lookup paths as globals that, in
# the package, exist only inside the composer. It now takes them as arguments;
# this checks it returns a real table for a known gene in a fresh session.
test_that("transcript_annotations() returns a table for a known gene, no app globals needed", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  g <- GENOMES$mouse
  ann <- transcript_annotations("Nfatc1", TXDBS$mouse, g$orgdb, g$tx_name_lookup_path, g$tx_meta_lookup_path)
  expect_s3_class(ann, "data.frame")
  expect_gt(nrow(ann), 1)
  expect_true(all(grepl("^Nfatc1-", ann$transcript_name)))
  expect_true(all(c("coding_status", "ensembl_canonical", "gencode_basic") %in% names(ann)))
})
