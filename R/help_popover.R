# Help popovers: the circled-i button beside a field label, in the setup wizard and the
# composer's Figure settings alike. Bootstrap 3's popover (bundled with shiny), opened by
# hover and by click (a click pins it, so it works on touch too; Escape or a click elsewhere
# closes it), driven by inst/www/plottrellis-help.js with delegated handlers, so an icon
# Shiny renders later needs no re-initialisation. Replaces the native `title` tooltip, which
# was hover-only, small, delayed and unusable on touch.

# The stylesheet and script, for a page's <head>; both apps include them next to
# .pt_theme_css().
.pt_help_assets <- function() {
  tagList(
    includeCSS(system.file("www", "plottrellis-help.css", package = "PlotTrellis", mustWork = TRUE)),
    includeScript(system.file("www", "plottrellis-help.js", package = "PlotTrellis", mustWork = TRUE))
  )
}

# A label with a help button after it. A <button> rather than a <span>: inside a <label>
# a click on a button does not activate the label's input (a span's would toggle a
# checkbox). The text travels in data-pt-help, read by the script when the popover opens.
wiz_help <- function(label_text, help_text) {
  tagList(label_text, tags$button(
    type = "button", class = "pt-help", `aria-label` = "Help", `data-pt-help` = help_text, "i"
  ))
}
