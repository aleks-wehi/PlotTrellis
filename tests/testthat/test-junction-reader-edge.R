# Junction/arc file edge cases (bugs 145, 146, 147, 151, 157, 170, 225) and an absolute path in a
# samples.txt cell (bug 215). The wizard's check, the arc-kind check and the reader must agree.

edge_file <- function(lines, ext = ".bed") {
  f <- tempfile("synth-edge-", fileext = ext)
  writeLines(lines, f)
  f
}
bed6_lines <- c("chr1\t3000\t7000\t.\t12\t+", "chr1\t4000\t6000\t.\t2\t-")
bed6_want <- data.frame(intron_start = c(3000L, 4000L), intron_end = c(7000L, 6000L), nSupport = c(12L, 2L))
bedpe_lines <- c("chr1\t2990\t3000\tchr1\t7000\t7010\tj1\t12", "chr1\t3990\t4000\tchr1\t6000\t6010\tj2\t2")
read_sorted <- function(f) {
  d <- read_junction_file(f, "chr1", 1000L, 10000L)
  d <- d[order(d$intron_start), ]
  rownames(d) <- NULL
  d
}

test_that("header lines before the data are skipped: comment, track, browser, column names (bug 145) (no fixture data)", {
  heads <- list(comment = "# converted from SJ.out.tab",
                track = "track name=junctions description=\"TopHat junctions\"",
                browser = c("browser position chr1:1-10000", "track name=j"),
                names = "chrom\tstart\tend\tname\tscore\tstrand",
                blank_then_comment = c("", "#chrom\tstart"))
  for (h in names(heads)) {
    f <- edge_file(c(heads[[h]], bed6_lines))
    expect_true(is_paired_end_structural(f), info = h)
    expect_identical(.pt_arc_file_kind(f), "junctions", info = h)
    expect_identical(read_sorted(f), bed6_want, info = h)
  }
  # A TopHat junctions.bed: BED12 under a track line.
  tophat <- edge_file(c("track name=junctions description=\"TopHat junctions\"",
                        "chr1\t2990\t7010\tJUNC1\t12\t+\t2990\t7010\t255,0,0\t2\t10,10\t0,4010"))
  expect_true(is_paired_end_structural(tophat))
  expect_identical(.pt_arc_file_kind(tophat), "junctions")
  expect_identical(read_sorted(tophat), data.frame(intron_start = 3000L, intron_end = 7000L, nSupport = 12L))
  bedpe <- edge_file(c("chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tname\tscore", bedpe_lines), ".bedpe")
  expect_identical(read_sorted(bedpe), bed6_want)
  # A first row written as "1e+05" or "3000.0" (R's write.table, pandas floats) is data, not names.
  sci <- edge_file(c("chr1\t1e+05\t110000\tchr1\t3e+05\t310000\tl1\t5",
                     "chr1\t250000\t260000\tchr1\t4e+05\t410000\tl2\t7"), ".bedpe")
  expect_identical(nrow(read_junction_file(sci, "chr1", 1L, 1000000L)), 2L)
  flt <- edge_file(c("chr1\t3000.0\t7000.0\t.\t12\t+", bed6_lines[2]))
  expect_identical(nrow(read_junction_file(flt, "chr1", 1000L, 10000L)), 2L)
})

test_that("a '.' score counts each junction or loop once (bug 146) (no fixture data)", {
  bed6 <- edge_file(sub("\t[0-9]+\t([-+])$", "\t.\t\\1", bed6_lines))
  bedpe <- edge_file(sub("\t[0-9]+$", "\t.", bedpe_lines), ".bedpe")
  mixed <- edge_file(c(bedpe_lines[1], sub("\t[0-9]+$", "\t.", bedpe_lines[2])), ".bedpe")
  ones <- data.frame(intron_start = c(3000L, 4000L), intron_end = c(7000L, 6000L), nSupport = c(1L, 1L))
  expect_equal(read_sorted(bed6), ones)
  expect_equal(read_sorted(bedpe), ones)
  expect_equal(read_sorted(mixed)$nSupport, c(12, 1))
  expect_equal(pool_junctions(list(list(junction_path = bed6), list(junction_path = bed6)),
                              "chr1", 1000L, 10000L, 0)$nSupport, c(2, 2))
  expect_identical(.pt_arc_file_kind(bedpe), "unscored")  # no score at all: no height in the heading
  expect_identical(.pt_arc_file_kind(bed6), "junctions")  # still splice junctions: keeps its TPM-scaled track
  expect_identical(.pt_arc_file_kind(mixed), "loops")
  bad <- edge_file(sub("\t[0-9]+$", "\tlots", bedpe_lines), ".bedpe")
  expect_error(read_junction_file(bad, "chr1", 1000L, 10000L), "column 8 .*a number or '\\.'.*found 'lots'")
})

test_that("a trailing tab on every line reads the same (bug 147) (no fixture data)", {
  bed6 <- edge_file(paste0(bed6_lines, "\t"))
  bedpe <- edge_file(paste0(bedpe_lines, "\t"), ".bedpe")
  expect_true(is_paired_end_structural(bed6))
  expect_identical(.pt_arc_file_kind(bed6), "junctions")
  expect_identical(read_sorted(bed6), bed6_want)
  expect_identical(.pt_arc_file_kind(bedpe), "loops")
  expect_identical(read_sorted(bedpe), bed6_want)
})

test_that("a narrowPeak chosen as paired-end stops with a peaks message, not an empty track (bug 151) (no fixture data)", {
  np <- edge_file(c("chr1\t3000\t3500\tpeak1\t500\t.\t12.5\t-1\t3.2\t250",
                    "chr1\t6000\t6400\tpeak2\t300\t.\t8.1\t-1\t2.0\t200"), ".narrowPeak")
  expect_false(is_paired_end_structural(np))
  expect_identical(.pt_arc_file_kind(np), NA_character_)
  expect_error(read_junction_file(np, "chr1", 1000L, 10000L), "looks like a BED peaks file")
})

test_that("an empty, whitespace-only or empty .gz junction file has no junctions (bug 157) (no fixture data)", {
  blank <- edge_file(c("", "  ", "\t"))
  gz <- tempfile("synth-edge-", fileext = ".bed.gz")
  close(gzfile(gz, "w"))
  expect_gt(file.size(gz), 0)
  for (f in c(blank, gz)) {
    expect_identical(nrow(read_junction_file(f, "chr1", 1000L, 10000L)), 0L, info = f)
    expect_identical(nrow(load_junctions(f, list(chrom = "chr1", start = 1000L, end = 10000L))), 0L, info = f)
  }
})

test_that("reader messages don't repeat the file path (bug 170) (no fixture data)", {
  three <- edge_file("chr1\t100\t200")
  strand <- edge_file("chr1\t100\t200\t.\t3\tx")
  for (f in c(three, strand)) {
    msg <- tryCatch(read_junction_file(f, "chr1", 1L, 1000L), error = conditionMessage)
    expect_false(grepl(basename(f), msg, fixed = TRUE), info = msg)
  }
})

# STAR SJ.out.tab: chrom, first and last intron base (1-based), strand 0/1/2, motif, annotated,
# unique reads, multi-mapping reads, overhang.
sj_lines <- c("chr1\t3001\t7000\t1\t1\t1\t12\t3\t40",
              "chr1\t4001\t6000\t2\t2\t0\t2\t0\t25",
              "chr1\t101\t300\t0\t0\t0\t8\t1\t10",
              "chr2\t3001\t7000\t1\t1\t1\t40\t0\t50")

test_that("a STAR SJ.out.tab reads as junctions with unique-read counts (bug 225) (no fixture data)", {
  sj <- edge_file(sj_lines, ".SJ.out.tab")
  bare <- edge_file(sub("^chr", "", sj_lines), ".SJ.out.tab")
  expect_true(is_paired_end_structural(sj))
  expect_identical(.pt_arc_file_kind(sj), "junctions")
  for (f in c(sj, bare)) expect_identical(read_sorted(f), bed6_want, info = f)
  # The same introns as a junction BED6 draw the same arcs.
  expect_identical(read_sorted(sj), read_sorted(edge_file(bed6_lines)))
  # A 9-column BEDPE (strand1 in column 9) is still a BEDPE, even on chromosomes named 1 and 2.
  bedpe9 <- edge_file(c("1\t2990\t3000\t1\t7000\t7010\tj1\t12\t+", "2\t2990\t3000\t2\t7000\t7010\tj2\t5\t."), ".bedpe")
  expect_identical(.pt_arc_file_kind(bedpe9), "loops")
  # Even with a number for a name and in column 9: column 6 (end2) is not STAR's annotated 0/1.
  numeric9 <- edge_file("1\t2990\t3000\t1\t7000\t7010\t7\t12\t4", ".bedpe")
  expect_identical(read_junction_file(numeric9, "chr1", 1000L, 10000L),
                   data.frame(intron_start = 3000L, intron_end = 7000L, nSupport = 12L))
  expect_identical(read_junction_file(bedpe9, "chr1", 1000L, 10000L),
                   data.frame(intron_start = 3000L, intron_end = 7000L, nSupport = 12L))
})

test_that("an absolute path in a samples.txt cell is named as the problem, not as a missing file (bug 215)", {
  root <- tempfile("pt-root-")
  dir.create(root)
  elsewhere <- tempfile("pt-elsewhere-", fileext = ".bw")
  writeLines("x", elsewhere)
  writeLines("x", file.path(root, "rel.bw"))
  txt <- tempfile(fileext = ".txt")
  writeLines(c(paste0("# data_root: ", root), "condition\tbw:ATAC",
               paste0("A\t", elsewhere), "B\trel.bw"), txt)
  res <- import_samples_txt(txt)
  expect_length(res$errors, 0)
  expect_length(res$warnings, 1)
  expect_match(res$warnings, "Row 2 (A)", fixed = TRUE)
  expect_match(res$warnings, "absolute path", fixed = TRUE)
  expect_match(res$warnings, "relative to the data root", fixed = TRUE)
})
