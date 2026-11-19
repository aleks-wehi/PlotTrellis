# The main screen's page (.pt_main_app(), R/app-screen-server.R): the markup,
# stylesheet and script of inst/www/main/, plus the one part that stays Shiny
# inputs, the Figure settings panel (its inputs feed live_cfg(), see
# R/app-settings.R). The page keeps its own state in the browser and asks
# the session for everything else (R/app-screen-server.R).

# Which screen the app opens with: "new" (this one) or "classic" (the composer
# of R/app-ui.R, kept as the fallback and as what the older tests drive).
# options(PlotTrellis.screen = ) wins over the PLOTTRELLIS_SCREEN environment
# variable, which a launcher script can set without touching R code.
.pt_screen <- function() {
  x <- getOption("PlotTrellis.screen", Sys.getenv("PLOTTRELLIS_SCREEN", "new"))
  if (identical(x, "classic")) "classic" else "new"
}

# A file of the main screen (inst/www/main/).
.pt_main_www <- function(file)
  system.file("www", "main", file, package = "PlotTrellis", mustWork = TRUE)

# The Figure settings panel: the same inputs as the classic screen's section
# (plot title, every display setting, Save as project default, Revert), in a
# panel the page's "Figure settings" button shows beside the figure. It is in
# the page from the start and only hidden, never removed: live_cfg() reads its
# inputs on every render.
.pt_main_settings_panel <- function() {
  tags$aside(id = "ptSettings", hidden = NA, `aria-label` = "Figure settings",
    div(class = "pt-set-head",
      tags$h3("Figure settings"),
      tags$button(id = "ptSettingsClose", type = "button", class = "close",
                  `aria-label` = "Close Figure settings", HTML("&times;"))
    ),
    .pt_plot_title_input(),
    p(class = "text-muted",
      "Changes apply to the next figure you generate, for this session only. Save as project default keeps them for next time."),
    uiOutput("display_settings_ui"),
    div(class = "pt-set-foot",
      actionButton("save_settings", "Save as project default",
                   class = "btn-default pt-once"),
      actionButton("revert_settings", "Revert to saved",
                   class = "btn-default pt-once"),
      uiOutput("settings_dirty_ui", inline = TRUE)
    )
  )
}

# The page. bootstrapPage(), not navbarPage(): the title bar is the page's own,
# and Bootstrap is loaded only for the Shiny parts (the settings panel's inputs
# and the project dialogs); main.css is written to sit on top of it. The script
# comes last: it runs once, after the markup it addresses.
.pt_main_ui <- function(ctx) {
  bootstrapPage(
    title = "PlotTrellis",
    tags$head(.pt_theme_css(), .pt_once_js(), .pt_help_assets(),
              includeCSS(.pt_main_www("main.css"))),
    includeHTML(.pt_main_www("main.html")),
    .pt_main_settings_panel(),
    includeScript(.pt_main_www("main.js"))
  )
}
