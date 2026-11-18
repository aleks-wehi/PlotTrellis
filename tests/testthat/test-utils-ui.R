# The composer's shared UI pieces (R/app-utils-ui.R): the Generate button, the collapsible
# section shell and the inline SVG wrapper. Plain HTML, no app.

test_that("an enabled Generate button is an action button with the generate-btn class", {
  html <- as.character(.pt_generate_button_ui("generate", TRUE))
  expect_match(html, 'id="generate"', fixed = TRUE)
  expect_match(html, "action-button", fixed = TRUE)
  expect_match(html, "btn-success btn-lg generate-btn", fixed = TRUE)
  expect_match(html, "Generate figure \u2192", fixed = TRUE)
  expect_false(grepl("disabled", html, fixed = TRUE))
})

test_that("a disabled Generate button is a plain disabled button with no input id", {
  html <- as.character(.pt_generate_button_ui("generate", FALSE))
  expect_identical(html, paste0('<button class="btn btn-success btn-lg generate-btn" disabled="disabled">',
                                'Generate figure \u2192</button>'))
  expect_false(grepl("id=", html, fixed = TRUE))
  expect_false(grepl("action-button", html, fixed = TRUE))
})

test_that("the button id is the one asked for", {
  expect_match(as.character(.pt_generate_button_ui("other_btn", TRUE)), 'id="other_btn"', fixed = TRUE)
})

test_that("an open collapsible panel shows its body, the down arrow and its children", {
  html <- as.character(.pt_collapsible_panel("sec", "Title", p("x"), "y"))
  expect_match(html, '<div class="panel panel-default plottrellis-section ">', fixed = TRUE)
  expect_match(html, '<a href="#sec-body" data-toggle="collapse">', fixed = TRUE)
  expect_match(html, '<span class="collapse-arrow">\u25BE</span>', fixed = TRUE)
  expect_match(html, " Title", fixed = TRUE)
  expect_match(html, '<div id="sec-body" class="collapse in">', fixed = TRUE)
  expect_match(html, '<div class="panel-body">\\s*<p>x</p>\\s*y\\s*</div>')
  expect_false(grepl("collapsed-panel", html, fixed = TRUE))
  # No wrapper id or style unless asked for.
  expect_false(grepl('plottrellis-section[^>]*(id|style)=', html))
})

test_that("a collapsed panel hides its body and shows the right arrow", {
  html <- as.character(.pt_collapsible_panel("sec", "Title", collapsed = TRUE))
  expect_match(html, "plottrellis-section collapsed-panel", fixed = TRUE)
  expect_match(html, '<span class="collapse-arrow">\u25B8</span>', fixed = TRUE)
  expect_match(html, '<div id="sec-body" class="collapse ">', fixed = TRUE)
  expect_false(grepl("collapse in", html, fixed = TRUE))
  expect_match(html, '<div class="panel-body"></div>', fixed = TRUE)
})

test_that("a panel's wrapper takes the id and style given", {
  html <- as.character(.pt_collapsible_panel("output", "Output", wrapper_id = "output-section",
                                             wrapper_style = "display:none;"))
  expect_match(html, '<div id="output-section" class="panel panel-default plottrellis-section " style="display:none;">',
               fixed = TRUE)
})

test_that(".pt_inline_svg sizes the first <svg> tag only and wraps it in the capped box", {
  f <- tempfile(fileext = ".svg"); on.exit(unlink(f), add = TRUE)
  writeLines(c('<?xml version="1.0"?>', '<svg width="10" height="5"><svg x="1"/>',
               "<text>a &amp; b</text>", "</svg>"), f)
  html <- as.character(.pt_inline_svg(f))
  expect_match(html, '^<div style="width:100%; max-width:1400px; margin:0 auto; overflow-x:auto;">')
  expect_match(html, "<svg style='width:100%;height:auto;' width=\"10\" height=\"5\"><svg x=\"1\"/>", fixed = TRUE)
  expect_equal(lengths(regmatches(html, gregexpr("width:100%;height:auto;", html, fixed = TRUE))), 1)
  # The SVG is inserted as is, not escaped, with its lines joined by newlines.
  expect_match(html, "\n<text>a &amp; b</text>\n</svg></div>", fixed = TRUE)
})

test_that(".pt_inline_svg leaves a file with no <svg tag unsized, and needs the file to exist", {
  f <- tempfile(fileext = ".svg"); on.exit(unlink(f), add = TRUE)
  writeLines("<g></g>", f)
  expect_identical(as.character(.pt_inline_svg(f)),
                   '<div style="width:100%; max-width:1400px; margin:0 auto; overflow-x:auto;"><g></g></div>')
  expect_error(suppressWarnings(.pt_inline_svg(file.path(tempdir(), "no-such.svg"))), "cannot open")
})
