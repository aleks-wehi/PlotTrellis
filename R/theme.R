# The app's look (colours, type, buttons, track slots), as one <style> tag for the composer's
# and the setup wizard's page head: inst/www/plottrellis.css, read at page build.
.pt_theme_css <- function() {
  includeCSS(system.file("www", "plottrellis.css", package = "PlotTrellis", mustWork = TRUE))
}
