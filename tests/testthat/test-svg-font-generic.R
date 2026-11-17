# An exported SVG's generic font fallback follows the chosen font (.pt_with_generic_font(), used by
# .pt_tidy_svg()): a viewer without the exact font then shows a face of the same kind.
test_that("the generic fallback is serif, monospace or sans-serif by the font's name", {
  expect_identical(.pt_with_generic_font('"DejaVu Sans"'), '"DejaVu Sans", sans-serif')
  expect_identical(.pt_with_generic_font("Helvetica"), "Helvetica, sans-serif")
  expect_identical(.pt_with_generic_font('"Nimbus Roman"'), '"Nimbus Roman", serif')
  expect_identical(.pt_with_generic_font('"Times New Roman"'), '"Times New Roman", serif')
  expect_identical(.pt_with_generic_font('"DejaVu Serif"'), '"DejaVu Serif", serif')
  expect_identical(.pt_with_generic_font('"Nimbus Mono PS"'), '"Nimbus Mono PS", monospace')
  expect_identical(.pt_with_generic_font('"DejaVu Sans Mono"'), '"DejaVu Sans Mono", monospace')
  expect_identical(.pt_with_generic_font("Courier"), "Courier, monospace")
  expect_identical(.pt_with_generic_font('"Fira Code"'), '"Fira Code", monospace')
  for (m in c("Consolas", "Inconsolata", "Lucida Console", "FreeMono", "Monaco"))
    expect_identical(.pt_with_generic_font(m), paste0(m, ", monospace"))
  expect_identical(.pt_with_generic_font("FreeSerif"), "FreeSerif, serif")
  # A word inside a longer one is not the word: these are sans faces (or not monospace).
  expect_identical(.pt_with_generic_font('"Encode Sans"'), '"Encode Sans", sans-serif')
  expect_identical(.pt_with_generic_font('"Lucida Sans Unicode"'), '"Lucida Sans Unicode", sans-serif')
  expect_identical(.pt_with_generic_font('"Monotype Corsiva"'), '"Monotype Corsiva", sans-serif')
  expect_identical(.pt_with_generic_font('"PT Serif"'), '"PT Serif", serif')
  expect_identical(.pt_with_generic_font('"Source Sans Pro"'), '"Source Sans Pro", sans-serif')
  # Already ending in a generic family: kept as it is, so tidying twice adds nothing.
  expect_identical(.pt_with_generic_font('"Nimbus Roman", serif'), '"Nimbus Roman", serif')
  expect_identical(.pt_with_generic_font("serif"), "serif")
})

test_that(".pt_tidy_svg() gives each text its own font's generic family, once", {
  f <- tempfile(fileext = ".svg")
  on.exit(unlink(f), add = TRUE)
  writeLines(c("<svg xmlns='http://www.w3.org/2000/svg' width='10pt' height='10pt'>",
               "<text x='1' y='2' style='font-family: \"Nimbus Roman\";'>serif</text>",
               "<text x='1' y='4' style='font-family: \"DejaVu Sans\";'>sans</text>",
               "<text x='1' y='6' style='font-family: \"Nimbus Mono PS\";'>mono</text>",
               "</svg>"), f)
  .pt_tidy_svg(f)
  .pt_tidy_svg(f)
  s <- readLines(f)
  expect_match(s[2], "font-family: \"Nimbus Roman\", serif;", fixed = TRUE)
  expect_match(s[3], "font-family: \"DejaVu Sans\", sans-serif;", fixed = TRUE)
  expect_match(s[4], "font-family: \"Nimbus Mono PS\", monospace;", fixed = TRUE)
})
