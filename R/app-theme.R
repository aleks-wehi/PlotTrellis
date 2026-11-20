# The app's look (colours, type, buttons), as one <style> tag for the main
# screen's and the setup wizard's page head: inst/www/plottrellis.css, read at
# page build.
.pt_theme_css <- function() {
  includeCSS(system.file("www", "plottrellis.css", package = "PlotTrellis",
                         mustWork = TRUE))
}

# The name and version in the setup wizard's title bar (.pt_title_bar()).
.pt_brand <- function() {
  tags$span("PlotTrellis",
            tags$small(class = "plottrellis-version",
                       paste0("v", utils::packageVersion("PlotTrellis"))))
}

# The setup wizard's title bar: the markup navbarPage() writes, without its
# (hidden) tab list, which is what the stylesheet's .navbar rules draw.
.pt_title_bar <- function() {
  tags$nav(class = "navbar navbar-default navbar-static-top",
           role = "navigation",
    div(class = "container-fluid",
      div(class = "navbar-header",
          tags$span(class = "navbar-brand", .pt_brand()))))
}
