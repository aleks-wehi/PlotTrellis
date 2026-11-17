# The bed track's reader on awkward BED files: no line is lost to the file's layout, a file
# with nothing to draw is an empty track, and a file that isn't BED says why. No fixture data.

bed_lines <- function(lines, ext = ".bed") {
  f <- tempfile(fileext = ext)
  writeLines(lines, f)
  f
}
starts <- function(lines) .pt_read_bed(bed_lines(lines))$V2

test_that("a track, browser, # or blank line in the middle loses no interval", {
  expect_equal(starts(c("chr1\t100\t200", "track name=x", "chr1\t300\t400", "chr1\t500\t600")), c(100, 300, 500))
  expect_equal(starts(c("chr1\t100\t200", "browser position chr1:1-1000", "chr1\t300\t400")), c(100, 300))
  expect_equal(starts(c("chr1\t100\t200", "chr1\t300\t400", "# comment", "chr1\t500\t600")), c(100, 300, 500))
  expect_equal(starts(c("chr1\t100\t200", "", "  ", "chr1\t300\t400")), c(100, 300))
  # A tab-only line (a spreadsheet export's empty row), mid-file or last.
  expect_equal(starts(c("chr1\t100\t200", "\t\t\t", "chr1\t300\t400", " \t ")), c(100, 300))
  # Past the first lines, and past what fread() samples of a file.
  many <- sprintf("chr1\t%d\t%d", 1:3000 * 100, 1:3000 * 100 + 50)
  for (odd in c("# comment", "track name=y", "", "\t\t", "chr1 1 2", "chr1\t1\t2\tx\t1\t+\textra"))
    expect_length(starts(c(many, odd, "chr1\t7\t8")), 3001 + (odd %in% c("chr1 1 2", "chr1\t1\t2\tx\t1\t+\textra")))
})

test_that("rows of different lengths and spaces for tabs lose no interval", {
  expect_equal(starts(c("chr1\t100\t200\tp1\t5\t+", "chr1\t300\t400", "chr1\t500\t600\tp3\t5\t+")), c(100, 300, 500))
  expect_equal(starts(c("chr1\t100\t200", "chr1\t300\t400\tp2\t5\t+", "chr1\t500\t600")), c(100, 300, 500))
  expect_equal(starts(c("chr1\t100\t200", "chr1 300 400", "chr1  500\t600")), c(100, 300, 500))
  expect_equal(starts(c("chr1 100 200 p1", "chr1  300   400 p2", "chr1 500 600")), c(100, 300, 500))
  # A chromosome named like "tig..." or "b..." is not a track or browser line.
  expect_equal(.pt_read_bed(bed_lines(c("tig1\t100\t200", "bchr\t300\t400")))$V1, c("tig1", "bchr"))
})

test_that("a header row without '#', a compressed file and a numeric chromosome name read as before", {
  expect_equal(starts(c("chrom\tstart\tend", "chr1\t100\t200", "chr1\t300\t400")), c(100, 300))
  gz <- tempfile(fileext = ".bed.gz")
  con <- gzfile(gz, "w"); writeLines(c("track name=x", "chr1\t100\t200", "", "chr1\t300\t400"), con); close(con)
  expect_equal(.pt_read_bed(gz)$V2, c(100, 300))
  expect_identical(.pt_read_bed(bed_lines("1\t100\t200"))$V1, "1")
  expect_named(.pt_read_bed(bed_lines("chr1\t100\t200\tp\t5\t+")), c("V1", "V2", "V3"))
  # Line ends converted twice ("\r\r\n") after a header row; a commented-out interval.
  f <- tempfile(fileext = ".bed")
  writeBin(charToRaw("chrom\tstart\tend\r\r\nchr1\t100\t200\r\r\nchr1\t300\t400\r\r\n"), f)
  expect_equal(.pt_read_bed(f)$V2, c(100, 300))
  expect_equal(starts(c("chr1\t100\t200", "#chr1\t150\t250", "chr1\t300\t400")), c(100, 300))
})

test_that("a file with only a header row, or only blank lines, is an empty track", {
  region <- list(chrom = "chr1", start = 1, end = 1000)
  expect_true(load_bed(list(path = bed_lines("chrom\tstart\tend")), region)$data$empty)
  expect_true(load_bed(list(path = bed_lines(c("  ", "", "\t"))), region)$data$empty)
  expect_true(load_bed(list(path = bed_lines(c("track name=x", "#chrom\tstart\tend"))), region)$data$empty)
  expect_true(load_bed(list(path = bed_lines(c("# a", "   ", "\t"))), region)$data$empty)
  expect_true(load_bed(list(path = bed_lines("chr2\t100\t200")), region)$data$no_chrom)
})

test_that("a file that isn't BED stops with a plain reason", {
  f <- bed_lines(c("chr1\ta\tb", "chr1\tc\td"))
  expect_error(.pt_read_bed(f), paste0(basename(f), " is not a BED file: no line has numbers for start and end"),
               fixed = TRUE)
  f <- bed_lines(c("chrom start end", "chr1 100 x"))   # the example is past the header row
  expect_error(.pt_read_bed(f), "no line has numbers for start and end (columns 2 and 3), e.g. line 2: chr1 100 x",
               fixed = TRUE)
  f <- bed_lines(c("chr1\t100\t200", "chr1\t300"))
  expect_error(.pt_read_bed(f), paste0(basename(f), " is not a BED file: line 2 has 2 fields, and needs at least 3"),
               fixed = TRUE)
  f <- bed_lines(c("chr1\t100\t200", "", "chr1\t300\tx"))
  expect_error(.pt_read_bed(f), paste0(basename(f), " is not a BED file: line 3 has a start or end that isn't a number"),
               fixed = TRUE)
})

# render_bed() draws the rows load_bed() read (the box counts are checked in
# test-setting-effects-display.R).
test_that("the bed track gets every interval of a file with a track line and a blank line mid-file", {
  f <- bed_lines(c("chr1\t2000\t2500", "track name=x", "chr1\t5000\t5800", "", "chr1\t9000\t9500"))
  e <- load_bed(list(path = f), list(chrom = "chr1", start = 1, end = 10000))
  expect_identical(e$data$chrom, "chr1")
  expect_equal(e$data$bed$V2, c(2000, 5000, 9000))
})
