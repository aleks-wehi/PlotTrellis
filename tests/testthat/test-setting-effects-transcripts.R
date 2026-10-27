# Phase I.3 step 28 (stage B): per-setting effect tests for the transcript /
# TPM tracks (min_display_tpm through a real transcripts-track render,
# fontsize on the tpm_usage track's ggplot text). test-caching.R already
# covers min_display_tpm at the gene_display_transcripts() level (bug 61);
# the render test here also catches the setting not reaching the track.
# Needs the mouse fixture (Elf5, quant.sf).

render_svg <- function(items, region, settings) {
  f <- tempfile(fileext = ".svg")
  on.exit(unlink(f), add = TRUE)
  dev <- NULL
  on.exit(if (!is.null(dev) && dev %in% grDevices::dev.list()) grDevices::dev.off(dev), add = TRUE, after = FALSE)
  utils::capture.output(render_all(items, region, settings, device_fn = function(w, h) {
    svglite::svglite(f, width = w, height = h)
    dev <<- grDevices::dev.cur()
  }))
  grDevices::dev.off(dev)
  dev <- NULL
  paste(readLines(f), collapse = "\n")
}

svg_texts <- function(svg) regmatches(svg, gregexpr("<text[^>]*>[^<]*</text>", svg))[[1]]
text_label <- function(t) sub(".*>([^<]*)</text>$", "\\1", t)
text_size  <- function(t) as.numeric(sub(".*font-size: ([0-9.]+)px.*", "\\1", t))

test_that("min_display_tpm changes which transcripts the transcripts track draws (mouse dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  rows   <- test_rows("mouse")
  region <- test_region("mouse")
  spec   <- single_track("transcripts", rows[["Ba"]])
  base   <- test_settings("mouse")
  on.exit(.set_transcript_filter_globals(base), add = TRUE)
  # Fresh cache file per value, so a stale cache entry can't mask or fake the effect.
  local_mocked_bindings(GENE_DISPLAY_TRANSCRIPTS_CACHE_PATH = tempfile(fileext = ".rds"))
  def  <- text_label(svg_texts(render_svg(list(spec), region, modifyList(base, list(min_display_tpm = 1)))))
  local_mocked_bindings(GENE_DISPLAY_TRANSCRIPTS_CACHE_PATH = tempfile(fileext = ".rds"))
  none <- text_label(svg_texts(render_svg(list(spec), region, modifyList(base, list(min_display_tpm = 1e9)))))
  expect_true(any(grepl("^Elf5-", def)))
  expect_false("No retained transcripts for this gene in view" %in% def)
  expect_false(any(grepl("^Elf5-", none)))
  expect_true("No retained transcripts for this gene in view" %in% none)
})

test_that("fontsize scales the tpm_usage track's bar labels and legend (mouse dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  rows   <- test_rows("mouse")
  region <- test_region("mouse")
  spec   <- single_track("tpm_usage", rows[["Ba"]])
  sizes <- function(fontsize) {
    t <- svg_texts(render_svg(list(spec), region, test_settings("mouse", overrides = list(fontsize = fontsize))))
    lab <- text_label(t)
    list(bar = text_size(t[grepl("^[0-9]+%$", lab)]), legend = text_size(t[grepl("^Elf5-", lab)]))
  }
  def <- sizes(8L)
  big <- sizes(16L)
  expect_gt(length(def$bar), 0)
  expect_gt(length(def$legend), 0)
  # Default reproduces the old fixed sizes: geom_text size 2.6 (mm) and legend text 7.
  expect_equal(def$bar, rep(2.6 * 72.27 / 25.4, length(def$bar)), tolerance = 0.01)
  expect_equal(def$legend, rep(7, length(def$legend)), tolerance = 0.01)
  expect_equal(big$bar, def$bar * 2, tolerance = 0.01)
  expect_equal(big$legend, rep(15, length(big$legend)), tolerance = 0.01)
})
