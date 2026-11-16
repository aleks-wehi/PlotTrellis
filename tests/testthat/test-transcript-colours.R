# One transcript, one colour (W64): every track that colours transcripts takes its colour from
# one map per figure, so a transcript keeps its colour across Transcripts, Transcripts (filtered),
# Transcripts annotated (filtered) and Transcript expression, even when their sets differ.
# Needs the mouse (Elf5, 8 transcripts) and human (PTPRC, more than 8) fixtures.

test_that("transcript_colours() repeats the palette past its length", {
  names <- sprintf("T%02d", 1:10)
  cols <- transcript_colours(names)
  expect_equal(unname(cols[1:8]), TX_PALETTE)
  expect_equal(unname(cols[9:10]), TX_PALETTE[1:2])
})

test_that("transcript_colours() gives `first` the first palette slots", {
  cols <- transcript_colours(sprintf("T%02d", 1:10), first = c("T10", "T09"))
  expect_equal(unname(cols[c("T09", "T10")]), TX_PALETTE[1:2])
  expect_equal(unname(cols["T01"]), TX_PALETTE[3])
})

test_that(".pt_add_tx_names() keeps one copy of each name across tracks", {
  dict <- .pt_add_tx_names(list(), c("b", "a"))
  dict <- .pt_add_tx_names(dict, c("a", "c"), first = TRUE)
  expect_setequal(dict[["tx_names"]], c("a", "b", "c"))
  expect_setequal(dict[["tx_names_first"]], c("a", "c"))
})

# Renders the four transcript tracks for one sample in one figure and returns every
# (track, transcript, colour) the tracks drew, "Other" left out.
figure_tx_colours <- function(dataset, sample) {
  rows     <- test_rows(dataset)
  region   <- test_region(dataset)
  settings <- test_settings(dataset, overrides = list(min_display_tpm = 1, cache_dir = tempfile("pt-cache-")))
  tx_all <- load_transcripts(list(viz = "transcripts_all"), region, settings)$data$tx_all
  id_to_name <- setNames(tx_all$transcript_name, tx_all$transcript_id)

  seen <- list()
  record <- function(track, name, colour) seen[[length(seen) + 1]] <<- data.frame(track = track, name = name, colour = colour)
  current <- NULL
  dispatch <- .pt_dispatch_fn
  local_mocked_bindings(
    plotTranscripts = function(transcriptFilter, fill, ...) record(current, id_to_name[[transcriptFilter]], fill),
    draw_transcript_counts = function(tpm, transcript_colours, ...)
      record("tpm_usage", colnames(tpm), transcript_colours),
    # Notes which track is drawing, for the plotTranscripts() records.
    .pt_dispatch_fn = function(prefix, key) {
      fn <- dispatch(prefix, key)
      if (prefix != "render") return(fn)
      function(...) { current <<- key; fn(...) }
    }
  )
  items <- lapply(c("transcripts_all", "transcripts", "transcripts_annotated", "tpm_usage"),
                  single_track, sample_row = rows[[sample]])
  svg <- tempfile(fileext = ".svg")
  on.exit(unlink(svg), add = TRUE)
  .pt_quiet(utils::capture.output(render_all(items, region, settings,
                                             device_fn = function(w, h) svglite::svglite(svg, width = w, height = h))))
  grDevices::dev.off()
  seen <- do.call(rbind, seen)
  seen[!grepl("^Other \\(n=", seen$name), ]
}

expect_one_colour_per_transcript <- function(seen) {
  expect_setequal(unique(seen$track), c("transcripts_all", "transcripts", "transcripts_annotated", "tpm_usage"))
  # The tracks' sets differ (else this test would not tell one map from four).
  expect_lt(sum(seen$track == "transcripts"), sum(seen$track == "transcripts_all"))
  per_name <- tapply(seen$colour, seen$name, function(x) length(unique(x)))
  expect_true(all(per_name == 1), info = paste(names(per_name)[per_name > 1], collapse = ", "))
  # The filtered tracks' transcripts (8 or fewer here) never share a colour, however many the
  # all-transcripts track adds.
  named <- unique(seen$name[seen$track == "tpm_usage"])
  expect_lte(length(named), length(TX_PALETTE))
  expect_equal(length(unique(seen$colour[match(named, seen$name)])), length(named))
}

test_that("a transcript has one colour in every transcript track of a figure (mouse dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  expect_one_colour_per_transcript(figure_tx_colours("mouse", "Ba"))
})

test_that("a transcript has one colour everywhere when the figure shows more than 8 (human dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  seen <- figure_tx_colours("human", "GM12878")
  expect_gt(length(unique(seen$name)), length(TX_PALETTE))
  expect_one_colour_per_transcript(seen)
})
