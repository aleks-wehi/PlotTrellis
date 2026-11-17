# Bug 110: a chromosome missing from the chrom sizes (a GTF scaffold UCSC does not list) gets
# the caller's default length; `CHR_LENGTHS[[chrom]] %||% default` stopped with "subscript out
# of bounds", which ended the app session in zoom_to().

test_that(".pt_chr_length() returns the length, or the default for a chromosome it lacks", {
  lens <- c(chr1 = 195154279L, chrX = 169476592L)
  expect_identical(.pt_chr_length(lens, "chrX", 2e8), 169476592L)
  expect_identical(.pt_chr_length(lens, "GL456210.1", 2e8), 2e8)
  expect_identical(.pt_chr_length(lens, NA_character_, Inf), Inf)
  expect_identical(.pt_chr_length(lens, character(0), 200000000L), 200000000L)
})
