# Bug 77: a long in-track message ("<track> track failed: <why>") was drawn at full length,
# past its track and the page. no_data_text() fits it to the track's slot.

measure <- function(fit, fontfamily = "sans") {
  g <- grid::textGrob(fit$label, gp = grid::gpar(fontsize = fit$fontsize, fontfamily = fontfamily))
  c(grid::convertWidth(grid::grobWidth(g), "inches", valueOnly = TRUE),
    grid::convertHeight(grid::grobHeight(g), "inches", valueOnly = TRUE))
}

test_that("a long track error is wrapped (and if need be shrunk) to its track's width and height", {
  svglite::svglite(tempfile(fileext = ".svg"), width = 8, height = 4)
  on.exit(grDevices::dev.off(), add = TRUE)
  label <- paste("junctions_diff_tpm (K562 vs GM12878) track failed: file not found: quant.sf",
                 "(in /a/rather/long/folder/name/for/this/sample/quant/K562_rep1_with_a_long_suffix)")
  w <- 5.18
  h <- 0.9
  # Unfitted it is wider than the track (the bug).
  expect_gt(measure(list(label = label, fontsize = 7))[1], w)

  fit <- .pt_fit_label(label, w, h, 7, "sans")
  expect_true(grepl("\n", fit$label, fixed = TRUE))
  expect_identical(gsub("\n", " ", fit$label, fixed = TRUE), label)
  expect_equal(fit$fontsize, 7)
  expect_lte(measure(fit)[1], w)
  expect_lte(measure(fit)[2], h)

  # A path with no spaces wider than the track, in a short slot: smaller text.
  fit <- .pt_fit_label(strrep("long_folder_name/", 7), 1, 0.2, 7, "sans")
  expect_lt(fit$fontsize, 7)
  expect_lte(measure(fit)[1], 1 + 1e-6)
  expect_lte(measure(fit)[2], 0.2 + 1e-6)

  # A message that fits is drawn as it is.
  expect_identical(.pt_fit_label("No gene selected for this region", w, h, 7, "sans"),
                   list(label = "No gene selected for this region", fontsize = 7))
})
