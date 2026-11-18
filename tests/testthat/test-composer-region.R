# The composer's coordinates block (R/app-region.R) on its own: a server made of
# nothing but .pt_region_server(), two made-up chromosomes, no genome. The page's input
# updates are recorded by the mock session (recording_session(), helper-composer-fakes.R);
# the mock session does not apply them, so each test sets the inputs the browser would.

CHRS <- c(chr1 = 1e7, chr2 = 5e6)

region_app <- function(chr_lengths = CHRS) function(input, output, session) {
  state <- reactiveValues(start_parse_err = NULL, end_parse_err = NULL)
  rg <- .pt_region_server(input, output, session, state, chr_lengths)
}

# The values sent to one input (updateTextInput()/updateSliderInput()), in order. A slider
# value comes as text, so it is turned back into numbers.
sent_values <- function(rec, id) {
  msgs <- Filter(function(m) identical(m$id, id), rec$input)
  lapply(msgs, function(m) if (id == "region_range") as.numeric(trimws(m$message$value)) else m$message$value)
}

test_that("the chromosome select lists the chromosomes in genome order, the first selected", {
  rec <- recording_session()
  shiny::testServer(region_app(c(chr10 = 1e6, chrX = 2e6, chr2 = 3e6, chrUn_a = 10)), {
    html <- ui_html(output$chr_select_ui)
    expect_match(html, 'id="region_chr"', fixed = TRUE)
    opts <- regmatches(html, gregexpr('option value="[^"]+"', html))[[1]]
    expect_identical(sub('option value="([^"]+)"', "\\1", opts), c("chr2", "chr10", "chrX", "chrUn_a"))
    expect_match(html, '<option value="chr2" selected>', fixed = TRUE)
  }, session = rec$session)
})

test_that("the slider and the text boxes span the chosen chromosome", {
  rec <- recording_session()
  shiny::testServer(region_app(), {
    session$setInputs(region_chr = "chr2")
    expect_identical(rg$chr_len_r(), 5e6)
    slider <- ui_html(output$range_slider_ui)
    expect_match(slider, 'data-min="1" data-max="5000000" data-from="1" data-to="5000000" data-step="1000"', fixed = TRUE)
    expect_match(ui_html(output$start_num_ui), 'value="1 bp"', fixed = TRUE)
    expect_match(ui_html(output$end_num_ui), 'value="5Mb"', fixed = TRUE)
    session$setInputs(region_chr = "chr1")
    expect_match(ui_html(output$end_num_ui), 'value="10Mb"', fixed = TRUE)
  }, session = rec$session)
})

test_that("chr_len_r waits for a chromosome, and an unknown one gets the default length", {
  rec <- recording_session()
  shiny::testServer(region_app(), {
    expect_error(rg$chr_len_r(), class = "shiny.silent.error")
    session$setInputs(region_chr = "chrZ")
    # Bug 110: a chromosome the chrom sizes lack gets the 200 Mb default.
    expect_identical(rg$chr_len_r(), 200000000L)
  }, session = rec$session)
})

test_that("moving the slider writes both ends into the text boxes", {
  rec <- recording_session()
  shiny::testServer(region_app(), {
    # The slider's first value (set in the same flush as the chromosome) is not echoed: the
    # text boxes start from their own values.
    session$setInputs(region_chr = "chr1", region_range = c(1, 1e7))
    expect_length(rec$input, 0)
    session$setInputs(region_range = c(1500, 2e6))
    expect_identical(sent_values(rec, "region_start_text"), list("1.5Kb"))
    expect_identical(sent_values(rec, "region_end_text"), list("2Mb"))
    session$setInputs(region_range = c(17490000, 17567890))
    expect_identical(sent_values(rec, "region_start_text")[[2]], "17.49Mb")
    expect_identical(sent_values(rec, "region_end_text")[[2]], "17.568Mb")
  }, session = rec$session)
})

test_that("a typed start or end moves that end of the slider", {
  rec <- recording_session()
  shiny::testServer(region_app(), {
    session$setInputs(region_chr = "chr1", region_range = c(1, 1e7))
    session$setInputs(region_start_text = "1.5Mb")
    expect_identical(sent_values(rec, "region_range"), list(c(1500000, 1e7)))
    session$setInputs(region_end_text = "2500kb")
    expect_identical(sent_values(rec, "region_range")[[2]], c(1, 2500000))
    expect_null(state$start_parse_err)
    expect_null(state$end_parse_err)
  }, session = rec$session)
})

test_that("a typed value equal to the slider's sends no update", {
  rec <- recording_session()
  shiny::testServer(region_app(), {
    session$setInputs(region_chr = "chr1", region_range = c(1000, 2000))
    session$setInputs(region_start_text = "1Kb")
    session$setInputs(region_end_text = "2000")
    expect_length(sent_values(rec, "region_range"), 0)
  }, session = rec$session)
})

test_that("typed values are not checked against the chromosome or each other", {
  rec <- recording_session()
  shiny::testServer(region_app(), {
    session$setInputs(region_chr = "chr2", region_range = c(1, 2e6))
    # Past chr2's end (5Mb): sent as typed; the slider's own max holds it in the browser.
    session$setInputs(region_end_text = "20Mb")
    expect_identical(sent_values(rec, "region_range")[[1]], c(1, 2e7))
    # A start after the end: sent as typed, then region_valid() reports it.
    session$setInputs(region_start_text = "3Mb")
    expect_identical(sent_values(rec, "region_range")[[2]], c(3e6, 2e6))
  }, session = rec$session)
})

test_that("an unreadable start or end sets its exact parse error, a good one clears it", {
  rec <- recording_session()
  shiny::testServer(region_app(), {
    session$setInputs(region_chr = "chr1", region_range = c(1, 1e7))
    session$setInputs(region_start_text = "17 megabases")
    expect_identical(state$start_parse_err, "Cannot parse start. Accepts: 17000000, 17Mb, 15.5kb (case-insensitive)")
    session$setInputs(region_end_text = "-5")
    expect_identical(state$end_parse_err, "Cannot parse end. Accepts: 17000000, 17Mb, 15.5kb (case-insensitive)")
    expect_length(sent_values(rec, "region_range"), 0)
    session$setInputs(region_start_text = "2Mb")
    expect_null(state$start_parse_err)
    session$setInputs(region_end_text = "3Mb")
    expect_null(state$end_parse_err)
  }, session = rec$session)
})

test_that("a blank text box changes nothing, not even an earlier parse error", {
  rec <- recording_session()
  shiny::testServer(region_app(), {
    session$setInputs(region_chr = "chr1", region_range = c(1, 1e7))
    session$setInputs(region_start_text = "abc")
    session$setInputs(region_start_text = "   ")
    # current behaviour, possibly a bug: emptying the box leaves the old parse error showing.
    expect_identical(state$start_parse_err, "Cannot parse start. Accepts: 17000000, 17Mb, 15.5kb (case-insensitive)")
    expect_length(sent_values(rec, "region_range"), 0)
  }, session = rec$session)
})

test_that("typed text before the slider exists is ignored, even when unreadable", {
  rec <- recording_session()
  shiny::testServer(region_app(), {
    session$setInputs(region_start_text = "abc", region_end_text = "2Mb")
    expect_null(state$start_parse_err)
    expect_null(state$end_parse_err)
    expect_length(rec$input, 0)
  }, session = rec$session)
})

test_that("region_valid accepts a region inside the chromosome, up to its last base", {
  rec <- recording_session()
  shiny::testServer(region_app(), {
    expect_error(rg$region_valid(), class = "shiny.silent.error")
    session$setInputs(region_chr = "chr2", region_range = c(1, 5e6))
    expect_identical(rg$region_valid(), list(valid = TRUE, chrom = "chr2", start = 1L, end = 5000000L))
    # Fractional slider values are truncated.
    session$setInputs(region_range = c(1000.9, 2000.7))
    expect_identical(rg$region_valid()[c("start", "end")], list(start = 1000L, end = 2000L))
  }, session = rec$session)
})

test_that("region_valid's three errors: start below 1, end not after start, end past the chromosome", {
  rec <- recording_session()
  shiny::testServer(region_app(), {
    session$setInputs(region_chr = "chr2", region_range = c(0, 100))
    expect_identical(rg$region_valid(), list(valid = FALSE, msg = "Start must be \u2265 1"))
    session$setInputs(region_range = c(100, 100))
    expect_identical(rg$region_valid(), list(valid = FALSE, msg = "End must be > start"))
    session$setInputs(region_range = c(200, 100))
    expect_identical(rg$region_valid()$msg, "End must be > start")
    session$setInputs(region_range = c(1, 5e6 + 1))
    expect_identical(rg$region_valid(), list(valid = FALSE, msg = "End exceeds chromosome length"))
  }, session = rec$session)
})

test_that("region_valid takes a missing start or end as an error", {
  # Set in the first flush, so the slider-to-text observer (ignoreInit) does not run: it
  # cannot format an NA, which the browser's slider never sends.
  rec <- recording_session()
  shiny::testServer(region_app(), {
    session$setInputs(region_chr = "chr2", region_range = c(NA, 100))
    expect_identical(rg$region_valid(), list(valid = FALSE, msg = "Start must be \u2265 1"))
  }, session = rec$session)
  rec <- recording_session()
  shiny::testServer(region_app(), {
    session$setInputs(region_chr = "chr2", region_range = c(100, NA))
    expect_identical(rg$region_valid(), list(valid = FALSE, msg = "End must be > start"))
  }, session = rec$session)
})

test_that("region_valid on an unknown chromosome has no upper bound", {
  rec <- recording_session()
  shiny::testServer(region_app(), {
    session$setInputs(region_chr = "chrZ", region_range = c(1, 1e9))
    # Bug 110: a chromosome the chrom sizes lack validates against Inf.
    expect_true(rg$region_valid()$valid)
    session$setInputs(region_mode = "coords")
    expect_null(output$region_error_ui)
  }, session = rec$session)
})

test_that("the region error box shows only in Coordinates mode", {
  rec <- recording_session()
  shiny::testServer(region_app(), {
    session$setInputs(region_chr = "chr2", region_range = c(200, 100))
    expect_null(output$region_error_ui)
    session$setInputs(region_mode = "gene")
    expect_null(output$region_error_ui)
    session$setInputs(region_mode = "coords")
    expect_match(ui_html(output$region_error_ui), '<div class="alert alert-danger">End must be &gt; start</div>', fixed = TRUE)
    session$setInputs(region_range = c(1, 100))
    expect_null(output$region_error_ui)
  }, session = rec$session)
})

test_that("parse errors take the place of the region error, joined by ' | '", {
  rec <- recording_session()
  shiny::testServer(region_app(), {
    session$setInputs(region_mode = "coords", region_chr = "chr2", region_range = c(200, 100))
    session$setInputs(region_start_text = "x")
    html <- ui_html(output$region_error_ui)
    expect_match(html, '^<div class="alert alert-warning">Cannot parse start\\.')
    expect_false(grepl("End must be", html, fixed = TRUE))
    session$setInputs(region_end_text = "y")
    expect_match(ui_html(output$region_error_ui),
                 "(case-insensitive) | Cannot parse end. Accepts", fixed = TRUE)
  }, session = rec$session)
})

test_that("the region error box is empty before a chromosome and region exist", {
  rec <- recording_session()
  shiny::testServer(region_app(), {
    session$setInputs(region_mode = "coords")
    expect_null(output$region_error_ui)
  }, session = rec$session)
})
