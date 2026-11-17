# The composer's Output section (R/composer-output.R) on its own: a state holding a small
# hand-written SVG file and R script, no render.

output_app <- function() function(input, output, session) {
  state <- reactiveValues(last_region = NULL, plot_error = NULL, svg_path = NULL, rcode = NULL)
  .pt_output_server(output, state)
}

# A small SVG of the kind svglite writes, white page background first.
small_svg <- function() {
  f <- tempfile(fileext = ".svg")
  writeLines(c("<?xml version='1.0' encoding='UTF-8' ?>",
               "<svg xmlns='http://www.w3.org/2000/svg' width='648.00pt' height='100.00pt' viewBox='0 0 648.00 100.00'>",
               "<rect width='100%' height='100%' style='stroke: none; fill: #FFFFFF;'/>",
               "<text x='10' y='20' style='font-size: 8.00px; font-family: DejaVu Sans;' textLength='20.00px' lengthAdjust='spacingAndGlyphs'>Pax5</text>",
               "</svg>"), f)
  f
}

test_that("nothing is shown before a figure exists", {
  shiny::testServer(output_app(), {
    expect_error(output$zoom_ui, class = "shiny.silent.error")
    expect_error(output$render_error_ui, class = "shiny.silent.error")
    expect_error(output$svg_display_ui, class = "shiny.silent.error")
    expect_identical(output$rcode_output, "")
  })
})

test_that("the zoom buttons appear once a figure has a region", {
  shiny::testServer(output_app(), {
    state$last_region <- list(chrom = "chr1", start = 1L, end = 1000L)
    session$flushReact()
    html <- ui_html(output$zoom_ui)
    for (id in c("zoom_in_15x", "zoom_in_3x", "zoom_in_10x", "zoom_out_15x", "zoom_out_3x", "zoom_out_10x"))
      expect_match(html, sprintf('id="%s"', id), fixed = TRUE)
    expect_match(html, '<div class="zoom-group-label">Zoom in</div>', fixed = TRUE)
    expect_match(html, '<div class="zoom-group-label">Zoom out</div>', fixed = TRUE)
    expect_identical(lengths(regmatches(html, gregexpr("1.5\u00D7", html, fixed = TRUE))), 2L)
  })
})

test_that("a render error shows twice, as text and in a pre block, escaped", {
  shiny::testServer(output_app(), {
    state$plot_error <- "bad <value> & more"
    session$flushReact()
    html <- ui_html(output$render_error_ui)
    expect_match(html, '^<div class="alert alert-danger">\\s*<strong>Error: </strong>')
    expect_match(html, "</strong>\n  bad &lt;value&gt; &amp; more\n  <br/>", fixed = TRUE)
    expect_match(html, '<pre style="font-size:11px;margin-top:8px;">bad &lt;value&gt; &amp; more</pre>', fixed = TRUE)
  })
})

test_that("the figure is shown inline only while its file exists", {
  f <- small_svg(); on.exit(unlink(f), add = TRUE)
  shiny::testServer(output_app(), {
    state$svg_path <- f
    session$flushReact()
    html <- ui_html(output$svg_display_ui)
    expect_match(html, "max-width:1400px", fixed = TRUE)
    expect_match(html, "<svg style='width:100%;height:auto;' xmlns=", fixed = TRUE)
    unlink(f)
    state$svg_path <- NULL; session$flushReact()
    state$svg_path <- f; session$flushReact()
    expect_error(output$svg_display_ui, class = "shiny.silent.error")
  })
})

test_that("the R script output shows the script", {
  shiny::testServer(output_app(), {
    # assemble_download_script() gives the script as one string.
    state$rcode <- "library(PlotTrellis)\nrender_all()"
    session$flushReact()
    expect_identical(output$rcode_output, "library(PlotTrellis)\nrender_all()")
  })
})

test_that("Download SVG gives a tidied copy, named by the time, and keeps the shown file as it is", {
  f <- small_svg(); on.exit(unlink(f), add = TRUE)
  before <- readLines(f)
  shiny::testServer(output_app(), {
    state$svg_path <- f
    session$flushReact()
    got <- output$dl_svg
    expect_match(basename(got), "^plottrellis_[0-9]{8}_[0-9]{6}\\.svg$")
    svg <- readLines(got)
    # .pt_tidy_svg(): the page background goes, the text unlocks, a font fallback is added.
    expect_false(any(grepl("<rect width='100%'", svg, fixed = TRUE)))
    expect_false(any(grepl("textLength", svg, fixed = TRUE)))
    expect_true(any(grepl("font-family: DejaVu Sans, sans-serif;", svg, fixed = TRUE)))
  })
  expect_identical(readLines(f), before)
})

test_that("Download R script writes the script lines", {
  shiny::testServer(output_app(), {
    state$rcode <- c("x <- 1", "y <- 2")
    session$flushReact()
    got <- output$dl_rcode
    expect_match(basename(got), "^plottrellis_[0-9]{8}_[0-9]{6}\\.R$")
    expect_identical(readLines(got), c("x <- 1", "y <- 2"))
  })
})

test_that("the downloads stop with no figure or no script", {
  shiny::testServer(output_app(), {
    expect_error(output$dl_svg, class = "shiny.silent.error")
    expect_error(output$dl_rcode, class = "shiny.silent.error")
    state$svg_path <- file.path(tempdir(), "no-such-figure.svg")
    session$flushReact()
    expect_error(output$dl_svg, class = "shiny.silent.error")
  })
})
