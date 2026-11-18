# Shared UI pieces of the composer (R/app-main.R): the collapsible section shell, the
# Generate button in both states, and the inline SVG wrapper of a rendered
# figure.

# The Generate button; its "generate-btn" class puts it under the page script's
# spinner/double-click gate, which matches that class, not an id.
.pt_generate_button_ui <- function(id, enabled) {
  if (enabled)
    actionButton(id, "Generate figure \u2192",
                 class = "btn-success btn-lg generate-btn")
  else
    tags$button("Generate figure \u2192",
                class = "btn btn-success btn-lg generate-btn",
                disabled = "disabled")
}

# Shared collapsible ".plottrellis-section" panel shell; the app-wide
# "hide.bs.collapse"/"show.bs.collapse" JS handlers key off any
# ".plottrellis-section .collapse" element (not a specific id), so this same
# shell works for every section.
.pt_collapsible_panel <- function(section_id, title, ..., wrapper_id = NULL,
                                  wrapper_style = NULL, collapsed = FALSE) {
  div(id = wrapper_id, class = paste("panel panel-default plottrellis-section",
                                     if (collapsed) "collapsed-panel"),
      style = wrapper_style,
    div(class = "panel-heading",
      tags$h4(class = "panel-title",
        tags$a(href = paste0("#", section_id, "-body"),
               `data-toggle` = "collapse",
          tags$span(class = "collapse-arrow",
                    if (collapsed) "\u25B8" else "\u25BE"),
          paste0(" ", title)
        )
      )
    ),
    div(id = paste0(section_id, "-body"),
        class = paste("collapse", if (!collapsed) "in"),
      div(class = "panel-body", ...)
    )
  )
}

.pt_inline_svg <- function(svg_path) {
  lines    <- readLines(svg_path, warn = FALSE)
  svg_text <- paste(lines, collapse = "\n")
  svg_text <- base::sub("<svg ", "<svg style='width:100%;height:auto;' ",
                        svg_text)
  # Capped + centred, not full-bleed: navbarPage is fluid, so an uncapped
  # width:100% stretches the SVG to however wide the monitor is -- since
  # height:auto preserves aspect ratio, that also inflates the height,
  # pushing tracks near the bottom (genes/transcripts) well past one
  # viewport. Capping keeps the figure at a legible, comparable size
  # regardless of screen width; narrower screens still shrink it via the
  # width:100% below the cap. No max-height/overflow-y here, deliberately: a
  # capped-height scroll box makes tall figures feel clipped/hidden.
  div(style = "width:100%; max-width:1400px; margin:0 auto; overflow-x:auto;",
      HTML(svg_text))
}
