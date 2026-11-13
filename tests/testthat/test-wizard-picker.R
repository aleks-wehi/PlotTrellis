# The wizard's file picker (wizard_scoped_choices()) lists the files each data type's loader reads,
# gzipped bed and paired-end files included (bug 142: they load, the picker used to hide them).

test_that("the picker lists gzipped bed and paired-end files, not other types' files", {
  root <- tempfile("pt-picker-"); dir.create(file.path(root, "sub"), recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  files <- c("a.bed", "b.bed.gz", "c.narrowPeak.gz", "sub/d.bedpe.gz", "e.bed12.gz", "f.bw", "g.bw.gz",
             "h.hic", "i.bed.gz.tbi", "j.txt.gz")
  file.create(file.path(root, files))
  expect_setequal(wizard_scoped_choices(root, "bed"), c("a.bed", "b.bed.gz", "c.narrowPeak.gz"))
  expect_setequal(wizard_scoped_choices(root, "paired-end"), c("a.bed", "b.bed.gz", "sub/d.bedpe.gz", "e.bed12.gz"))
  # bigWig and Hi-C files are binary formats their readers take only as they are.
  expect_equal(wizard_scoped_choices(root, "bw"), "f.bw")
  expect_equal(wizard_scoped_choices(root, "hic"), "h.hic")
})

test_that("on the human test data the picker lists its .bed.gz peaks and .bedpe.gz loops", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  root <- file.path(DATA_ROOT, "human_encode_test")
  bed <- wizard_scoped_choices(root, "bed")
  pe  <- wizard_scoped_choices(root, "paired-end")
  expect_true(all(c("atac/GM12878_ATAC_IDRpeaks.bed.gz", "atac/K562_ATAC_IDRpeaks.bed.gz") %in% bed))
  expect_true(all(c("loops/GM12878_CTCF_ChIAPET_loops_rep1.bedpe.gz", "loops/K562_CTCF_ChIAPET_loops_rep2.bedpe.gz",
                    "junctions/GM12878_rep1.bed") %in% pe))
})
