# The transcript tracks size as heading + one row per
# transcript + half a gap (no fixed height, no minimum); rows are 0.16in at the
# default font size and grow with bigger fonts; the first box sits one gap below
# the heading block; "+n not shown" shares the heading's baseline. Needs the
# mouse fixture (Elf5: 8 annotated transcripts, 2 kept at the default
# min_display_tpm).

track_height <- function(viz, settings, region) {
  rows <- test_rows("mouse")
  fn <- switch(viz, transcripts = settings_transcripts, transcripts_all = settings_transcripts_all,
               transcripts_annotated = settings_transcripts_annotated)
  # The settings render_all() would pass on (.pt_fill_settings()).
  fn(list(), single_track(viz, rows[["Ba"]]), .pt_fill_settings(settings), 9 - 0.8 - 1.0, region)$entry$height
}

layout_svg <- function(items, region, settings) {
  f <- tempfile(fileext = ".svg")
  on.exit(unlink(f), add = TRUE)
  utils::capture.output(render_all(items, region, settings,
                                   device_fn = function(w, h) svglite::svglite(f, width = w, height = h)))
  grDevices::dev.off()
  paste(readLines(f), collapse = "\n")
}
layout_texts <- function(svg) {
  t <- regmatches(svg, gregexpr("<text[^>]*>[^<]*</text>", svg))[[1]]
  data.frame(label = sub(".*>([^<]*)</text>$", "\\1", t),
             y = as.numeric(sub(".* y='([0-9.]+)'.*", "\\1", t)))
}

test_that("transcript track heights are heading + 0.16in per transcript, no minimum (mouse dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  region   <- test_region("mouse")
  settings <- test_settings("mouse")
  header <- 0.02 + 7 / 72            # default fontsize 8, heading at 7pt
  pad    <- (0.16 - 2 / 25.4) / 2     # half a gap below the last box
  expect_equal(track_height("transcripts_all", settings, region), header + 8 * 0.16 + pad)
  expect_equal(track_height("transcripts", settings, region), header + 2 * 0.16 + pad)
  expect_equal(track_height("transcripts_annotated", settings, region), header + 2 * 0.16 + pad)
  # A bigger font grows the heading block and the rows, by the 16pt the labels grew.
  row24 <- 0.16 + 16 / 72
  big <- modifyList(settings, list(fontsize = 24L))
  expect_equal(track_height("transcripts", big, region),
               0.02 + 23 / 72 + 2 * row24 + (row24 - 2 / 25.4) / 2)
  # Nothing kept: one row for the message, not a fixed default height.
  expect_equal(track_height("transcripts", modifyList(settings, list(min_display_tpm = 1e9, cache_dir = tempfile("pt-cache-"))),
                            region),
               header + 0.16 + pad)
})

test_that("transcript rows are evenly spaced below the heading, '+n not shown' on its baseline (mouse dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  rows     <- test_rows("mouse")
  region   <- test_region("mouse")
  # No title: the first track then starts at y = 0.1in.
  settings <- modifyList(test_settings("mouse"), list(show_title = FALSE))
  for (viz in c("transcripts", "transcripts_annotated")) {
    svg <- layout_svg(list(single_track(viz, rows[["Ba"]])), region, settings)
    tx <- layout_texts(svg)
    heading <- tx$y[grepl("(filtered)", tx$label, fixed = TRUE)]
    expect_lt(abs(tx$y[tx$label == "+6 not shown"] - heading), 0.01)  # svglite units: 1in = 72
    expect_lt(abs(diff(tx$y[tx$label %in% c("Elf5-201", "Elf5-207")]) - 0.16 * 72), 0.01)
    # The first transcript's full-height (2mm) CDS boxes start one gap below the heading block.
    r <- regmatches(svg, gregexpr("<rect [^>]* y='[0-9.]+' [^>]*height='[0-9.]+'[^>]*>", svg))[[1]]
    r_y <- as.numeric(sub(".* y='([0-9.]+)'.*", "\\1", r))
    r_h <- as.numeric(sub(".* height='([0-9.]+)'.*", "\\1", r))
    box_top <- min(r_y[abs(r_h - 2 / 25.4 * 72) < 0.05])
    expect_lt(abs(box_top - (0.1 + 0.02 + 7 / 72 + 0.16 - 2 / 25.4) * 72), 0.02)
  }
  # Fontsize 24: rows grow to 0.16 + 16/72 in (27.5pt), so 22.5pt labels don't overlap.
  tx <- layout_texts(layout_svg(list(single_track("transcripts", rows[["Ba"]])), region,
                                modifyList(settings, list(fontsize = 24L))))
  expect_lt(abs(diff(tx$y[tx$label %in% c("Elf5-201", "Elf5-207")]) - (0.16 + 16 / 72) * 72), 0.01)
})
