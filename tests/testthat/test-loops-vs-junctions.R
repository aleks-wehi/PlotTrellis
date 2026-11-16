# Loops (BEDPE) and splice junctions share the paired-end column type. Only junctions get the
# gene-TPM-scaled comparison and a "reads" label; a loop file's score isn't reads, and a BEDPE
# without a score has no height to report. No fixture data needed (a copy of the human example).

test_that("loops get no _diff_tpm track and no reads label; junctions keep both (no fixture data)", {
  skip_if_not_installed("org.Hs.eg.db")
  local_offline_chrom_sizes()
  txt <- example_copy("human")
  dir <- dirname(txt)
  # Two loops in the PSME2 view, one file without a score column, one with.
  writeLines(c("chr14\t24143500\t24144000\tchr14\t24146000\t24146500",
               "chr14\t24144000\t24144500\tchr14\t24147000\t24147500"), file.path(dir, "unscored.bedpe"))
  writeLines(c("chr14\t24143500\t24144000\tchr14\t24146000\t24146500\tl1\t7.5",
               "chr14\t24144000\t24144500\tchr14\t24147000\t24147500\tl2\t2"), file.path(dir, "scored.bedpe"))
  lines <- readLines(txt)
  hdr <- grep("^condition\t", lines)
  lines[hdr] <- paste0(lines[hdr], "\tpaired-end:loops\tpaired-end:scored")
  for (i in hdr + 1:2) lines[i] <- paste0(lines[i], "\tunscored.bedpe\tscored.bedpe")
  writeLines(lines, txt)
  proj <- suppressWarnings(pt_project(tempfile("pt-loops-"), txt))

  tracks <- pt_tracks(proj)$track
  expect_true(all(c("junctions_diff_tpm", "loops_diff", "scored_diff") %in% tracks))
  expect_false(any(c("loops_diff_tpm", "scored_diff_tpm") %in% tracks))

  svg_file <- tempfile(fileext = ".svg")
  on.exit(unlink(svg_file), add = TRUE)
  suppressMessages(pt_save(pt_figure(proj, "PSME2", c("junctions", "loops", "scored"), samples = "GM12878"), svg_file))
  svg <- paste(readLines(svg_file, warn = FALSE), collapse = "\n")
  expect_match(svg, "max_h=[0-9]+ reads")                      # junctions: read counts
  expect_match(svg, "max_h=7.5)", fixed = TRUE)                # scored loops: the score, not reads
  expect_identical(lengths(regmatches(svg, gregexpr("max_h=", svg, fixed = TRUE))), 2L)  # unscored: none
  expect_identical(lengths(regmatches(svg, gregexpr(" reads", svg, fixed = TRUE))), 1L)
})

test_that(".pt_arc_file_kind() tells junctions, scored loops and unscored loops apart (no fixture data)", {
  f <- function(...) { p <- tempfile(fileext = ".bed"); writeLines(c(...), p); p }
  bed12 <- "chr1\t10\t500\tj1\t5\t+\t10\t500\t0\t2\t10,10\t0,480"
  expect_identical(.pt_arc_file_kind(f(paste(c("#Chr", 2:12), collapse = "\t"), bed12)), "junctions")
  # A header read_junction_file() can't take (not 12 columns): unknown, so the track stays and shows its error.
  expect_identical(.pt_arc_file_kind(f("#Chr\tStart", bed12)), NA_character_)
  expect_identical(.pt_arc_file_kind(f("chr1\t10\t500\tj1\t5\t+")), "junctions")
  expect_identical(.pt_arc_file_kind(f("chr1\t10\t20\tchr1\t500\t510")), "unscored")
  expect_identical(.pt_arc_file_kind(f("chr1\t10\t20\tchr1\t500\t510\tl1\t3.2")), "loops")
  expect_identical(.pt_arc_file_kind(f("chr1\t10\t20\tchr1\t500\t510\tl1")), "unscored")
  expect_identical(.pt_arc_file_kind(f(paste(rep("x", 11), collapse = "\t"))), NA_character_)
  expect_identical(.pt_arc_file_kind(f("")), NA_character_)
  expect_identical(.pt_arc_file_kind(tempfile()), NA_character_)
})

test_that("a 7-column BEDPE takes arc heights from a numeric 7th column (ENCODE ChIA-PET) (no fixture data)", {
  f <- function(...) { p <- tempfile(fileext = ".bedpe"); writeLines(c(...), p); p }
  counted <- f("chr1\t100\t200\tchr1\t5000\t5100\t4", "chr1\t300\t400\tchr1\t6000\t6100\t8")
  expect_identical(read_junction_file(counted, "chr1", 1, 10000)$nSupport, c(4L, 8L))
  expect_identical(.pt_arc_file_kind(counted), "loops")
  # A text 7th column is a name: every arc counts once, as before.
  named <- f("chr1\t100\t200\tchr1\t5000\t5100\tloopA", "chr1\t300\t400\tchr1\t6000\t6100\tloopB")
  expect_identical(read_junction_file(named, "chr1", 1, 10000)$nSupport, c(1L, 1L))
  expect_identical(.pt_arc_file_kind(named), "unscored")
  # An NA on any row: every arc counts once (none dropped), and the heading shows no height.
  gap <- f("chr1\t100\t200\tchr1\t5000\t5100\tNA", "chr1\t300\t400\tchr1\t6000\t6100\t8")
  expect_identical(read_junction_file(gap, "chr1", 1, 10000)$nSupport, c(1L, 1L))
  expect_identical(.pt_arc_file_kind(gap), "unscored")
  # An 8-column BEDPE's score column still wins over a numeric name.
  scored <- f("chr1\t100\t200\tchr1\t5000\t5100\t17\t2.5")
  expect_identical(read_junction_file(scored, "chr1", 1, 10000)$nSupport, 2.5)
})
