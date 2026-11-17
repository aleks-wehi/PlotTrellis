# The app's look (colours, type, buttons, track slots), as one <style> tag for
# the composer's and the setup wizard's page head: inst/www/plottrellis.css,
# read at page build.
.pt_theme_css <- function() {
  includeCSS(system.file("www", "plottrellis.css", package = "PlotTrellis",
                         mustWork = TRUE))
}

# The title bar's name and version: navbarPage()'s title in the composer, and
# the same in the setup wizard's .pt_title_bar().
.pt_brand <- function() {
  tags$span("PlotTrellis",
            tags$small(class = "plottrellis-version",
                       paste0("v", utils::packageVersion("PlotTrellis"))))
}

# The composer's title bar, for a page that is not a navbarPage() (the setup
# wizard): the markup navbarPage() writes, without its (hidden) tab list, so the
# stylesheet draws the two the same.
.pt_title_bar <- function() {
  tags$nav(class = "navbar navbar-default navbar-static-top",
           role = "navigation",
    div(class = "container-fluid",
      div(class = "navbar-header",
          tags$span(class = "navbar-brand", .pt_brand()))))
}
