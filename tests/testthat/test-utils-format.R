# The composer's number, coordinate and label formatting (R/util-format.R), as plain
# functions: no app, no genome.

test_that(".pt_format_track_label bolds what precedes the first ' \u00B7 '", {
  expect_identical(.pt_format_track_label("atac \u00B7 Ba"), "<strong>atac</strong> \u00B7 Ba")
  # Only the first separator: a sample name can hold one too.
  expect_identical(.pt_format_track_label("a \u00B7 b \u00B7 c"), "<strong>a</strong> \u00B7 b \u00B7 c")
  # No separator: the whole label is bold (a sample-free track).
  expect_identical(.pt_format_track_label("Gene annotation"), "<strong>Gene annotation</strong>")
})

test_that(".pt_format_track_label leaves blank, NA and separator-first labels alone, and is vectorised", {
  expect_identical(.pt_format_track_label(""), "")
  expect_identical(.pt_format_track_label(NA_character_), NA_character_)
  expect_identical(.pt_format_track_label(" \u00B7 x"), " \u00B7 x")
  expect_identical(.pt_format_track_label(character(0)), character(0))
  expect_identical(.pt_format_track_label(c("x \u00B7 A", "Genes")),
                   c("<strong>x</strong> \u00B7 A", "<strong>Genes</strong>"))
})

test_that(".pt_format_track_label does not bold a name holding a bare middle dot", {
  # current behaviour, possibly a bug: a data name with a "\u00B7" not between spaces
  # (here "a\u00B7b") stops the match, so nothing of the label is bolded.
  expect_identical(.pt_format_track_label("a\u00B7b \u00B7 S"), "a\u00B7b \u00B7 S")
})

test_that(".pt_sort_chrs puts numbered chromosomes in number order, then X, Y, M, then the rest", {
  chrs <- c("chr10", "chr2", "chr1", "chrM", "chrY", "chrX", "chrUn_x", "chr1_random", "chr2_alt")
  expect_identical(.pt_sort_chrs(chrs),
                   c("chr1", "chr2", "chr10", "chrX", "chrY", "chrM", "chr1_random", "chr2_alt", "chrUn_x"))
  # X, Y, M keep that order whatever the input order, and a missing one is skipped.
  expect_identical(.pt_sort_chrs(c("chrM", "chrX", "chr1")), c("chr1", "chrX", "chrM"))
  # Contigs are kept, never filtered out.
  expect_setequal(.pt_sort_chrs(chrs), chrs)
})

test_that(".pt_sort_chrs handles an empty vector and names without 'chr'", {
  expect_identical(.pt_sort_chrs(character(0)), character(0))
  # current behaviour, possibly a bug: names without the "chr" prefix (Ensembl style) are
  # all treated as contigs and sorted as text, so "10" comes before "2".
  expect_identical(.pt_sort_chrs(c("2", "10", "1", "X", "MT")), c("1", "10", "2", "MT", "X"))
})

test_that(".pt_fmt_num rounds, then drops trailing zeros and a bare decimal point", {
  expect_identical(.pt_fmt_num(1.5), "1.5")
  expect_identical(.pt_fmt_num(2), "2")
  expect_identical(.pt_fmt_num(1.234), "1.23")
  expect_identical(.pt_fmt_num(1.23456, 3), "1.235")
  expect_identical(.pt_fmt_num(0), "0")
  expect_identical(.pt_fmt_num(10), "10")
  expect_identical(.pt_fmt_num(100), "100")
  expect_identical(.pt_fmt_num(-1.5), "-1.5")
  expect_identical(.pt_fmt_num(0.004), "0")
  expect_identical(.pt_fmt_num(c(1.5, 2)), c("1.5", "2"))
})

test_that(".pt_fmt_num with digits = 0 and with NA", {
  # current behaviour, possibly a bug: with digits = 0 there is no decimal point, so the
  # trailing-zero strip eats the integer's own zeros (10 -> "1", 100 -> "1"). No caller
  # passes digits = 0 today.
  expect_identical(.pt_fmt_num(10, 0), "1")
  expect_identical(.pt_fmt_num(100, 0), "1")
  expect_identical(.pt_fmt_num(7, 0), "7")
  # NA formats with formatC()'s leading space.
  expect_identical(.pt_fmt_num(NA), " NA")
})

test_that(".pt_format_bp picks bp, Kb or Mb at the 1,000 and 1,000,000 boundaries", {
  expect_identical(.pt_format_bp(0), "0bp")
  expect_identical(.pt_format_bp(1), "1bp")
  expect_identical(.pt_format_bp(999), "999bp")
  expect_identical(.pt_format_bp(1000), "1Kb")
  expect_identical(.pt_format_bp(1500), "1.5Kb")
  expect_identical(.pt_format_bp(1e6), "1Mb")
  expect_identical(.pt_format_bp(2.5e6), "2.5Mb")
  expect_identical(.pt_format_bp(1234567), "1.23Mb")
  expect_identical(.pt_format_bp(123456789), "123.46Mb")
  # current behaviour, possibly a bug: just under 1Mb rounds up to "1000Kb" rather than "1Mb".
  expect_identical(.pt_format_bp(999999), "1000Kb")
})

test_that(".pt_format_bp takes one number only", {
  expect_error(.pt_format_bp(c(1, 2)), "the condition has length > 1")
  expect_error(.pt_format_bp(NA), "missing value where TRUE/FALSE needed")
})

test_that(".pt_format_coord keeps three decimals and spaces the bp unit", {
  expect_identical(.pt_format_coord(0), "0 bp")
  expect_identical(.pt_format_coord(1L), "1 bp")
  expect_identical(.pt_format_coord(999), "999 bp")
  expect_identical(.pt_format_coord(1000), "1Kb")
  expect_identical(.pt_format_coord(1500), "1.5Kb")
  expect_identical(.pt_format_coord(1234.5), "1.234Kb")
  expect_identical(.pt_format_coord(999999), "999.999Kb")
  expect_identical(.pt_format_coord(1e6), "1Mb")
  expect_identical(.pt_format_coord(17490000), "17.49Mb")
  expect_identical(.pt_format_coord(1234567), "1.235Mb")
  expect_error(.pt_format_coord(NA), "missing value where TRUE/FALSE needed")
})

test_that(".pt_parse_bp reads plain integers and Kb/Mb values, in any case, as integers", {
  expect_identical(.pt_parse_bp("17000000"), 17000000L)
  expect_identical(.pt_parse_bp(" 17000000 "), 17000000L)
  expect_identical(.pt_parse_bp("0"), 0L)
  expect_identical(.pt_parse_bp("17Mb"), 17000000L)
  expect_identical(.pt_parse_bp("17mB"), 17000000L)
  expect_identical(.pt_parse_bp("17.56Mb"), 17560000L)
  expect_identical(.pt_parse_bp("15kb"), 15000L)
  expect_identical(.pt_parse_bp("150KB"), 150000L)
  expect_identical(.pt_parse_bp(".5Mb"), 500000L)
  expect_identical(.pt_parse_bp("1.5 kb"), 1500L)
  expect_identical(.pt_parse_bp("17 Mb "), 17000000L)
  # Rounded to a whole base.
  expect_identical(.pt_parse_bp("1.2345kb"), 1234L)
  expect_identical(.pt_parse_bp("0.0001kb"), 0L)
})

test_that(".pt_parse_bp gives NA for anything else", {
  for (s in c("17,000,000", "-5", "5bp", "", "   ", "abc", "1e6", "17.Mb", "1.5", "17m", "17 Gb"))
    expect_identical(.pt_parse_bp(s), NA_integer_, info = s)
  expect_identical(.pt_parse_bp(NA_character_), NA_integer_)
})

test_that(".pt_parse_bp past the integer range gives NA with R's coercion warning", {
  expect_warning(v <- .pt_parse_bp("99999999999"), "NAs introduced by coercion to integer range")
  expect_identical(v, NA_integer_)
  expect_warning(v <- .pt_parse_bp("3000Mb"), "NAs introduced by coercion to integer range")
  expect_identical(v, NA_integer_)
})

test_that(".pt_parse_bp needs one string", {
  expect_error(.pt_parse_bp(NULL), "argument is of length zero")
})
