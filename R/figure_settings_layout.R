# One subsection of the composer's Figure settings panel (.pt_display_settings_ui()):
# a native <details>, collapsed by default. Its inputs stay in the DOM while it is
# closed, so they still bind and report their values (the panel's output is rendered
# with suspendWhenHidden = FALSE). Styled inline, like the wizard's
# `details > summary` rule, because the composer page has no such rule and a global
# one would also restyle the Gene info table's <details> rows.
.pt_settings_group <- function(title, ...) {
  tags$details(class = "pt-settings-group", style = "margin-bottom:8px;",
    tags$summary(style = "display:list-item; cursor:pointer; font-weight:600;", title),
    tags$div(style = "margin-top:8px;", ...)
  )
}

# The per-figure title field, first in the Figure settings panel. Not a saved
# setting (no s_ prefix, so Save as project default never writes it).
.pt_plot_title_input <- function() {
  fluidRow(column(4,
    textInput("plot_title", wiz_help("Title (optional):", "This figure's title. Kept for zoom and Apply to figure; not saved by Save as project default."),
              placeholder = "Leave blank for no title", width = "100%")
  ))
}
