# The wizard's reason for refusing a file added as paired-end arcs (wiz_paired_end_problem()):
# none for the formats the arc tracks read; the BED type for a peaks file or a 3-5 column file;
# "no data lines" for a file with none; and a list naming STAR SJ.out.tab for anything else.
test_that("the paired-end check accepts arc formats and explains what it refuses", {
  dir <- tempfile("pe-"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  write <- function(name, ...) { p <- file.path(dir, name); writeLines(c(...), p); p }
  sj     <- write("SJ.out.tab", "chr2\t103410100\t103412000\t1\t1\t1\t12\t0\t40")
  bedpe  <- write("loops.bedpe", "chr2\t100\t200\tchr2\t5000\t5100\tloop1\t7")
  bed12  <- write("junctions.bed", paste("chr2", 100, 900, "j1", 5, "+", 100, 900, 0, 2, "50,50", "0,750", sep = "\t"))
  peaks  <- write("LP.narrowPeak", "chr2\t100\t600\tpeak1\t500\t.\t8.2\t12.1\t9.9\t250")
  plain  <- write("intervals.bed", "chr2\t100\t600\tregion1")
  odd    <- write("eleven.txt", paste(c("chr2", 1:10), collapse = "\t"))
  empty  <- write("empty.bed", "# no junctions", "track name=x")
  for (ok in c(sj, bedpe, bed12)) expect_null(wiz_paired_end_problem(ok), label = basename(ok))
  expect_match(wiz_paired_end_problem(peaks), "looks like a BED peaks file", fixed = TRUE)
  expect_match(wiz_paired_end_problem(peaks), "add it as Intervals (BED) instead", fixed = TRUE)
  expect_match(wiz_paired_end_problem(plain), "has 4 columns (arcs need 6 or more); if it is BED intervals, add it as Intervals (BED) instead", fixed = TRUE)
  other <- wiz_paired_end_problem(odd)
  expect_match(other, "STAR SJ.out.tab", fixed = TRUE)
  expect_false(grepl("peaks", other, fixed = TRUE))
  expect_match(wiz_paired_end_problem(empty), "has no data lines", fixed = TRUE)
})
