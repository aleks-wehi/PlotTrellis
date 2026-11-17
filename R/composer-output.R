# The composer's Output section (R/app.R): zoom buttons, the error box, the
# figure, the R script and the two downloads, all read from `state`.
.pt_output_server <- function(output, state) {
  # ── Output page ─────────────────────────────────────────────────────────────
  output$zoom_ui <- renderUI({
    req(!is.null(state$last_region))
    tags$div(
      style = "display:flex; justify-content:center; align-items:flex-start; gap:28px; padding:8px 0 12px;",
      tags$div(
        style = "text-align:center;",
        tags$div(class = "zoom-group-label",
                 "Zoom in"),
        actionButton("zoom_in_15x", "1.5\u00D7", class = "btn-default"),
        actionButton("zoom_in_3x",  "3\u00D7",   class = "btn-default"),
        actionButton("zoom_in_10x", "10\u00D7",  class = "btn-default")
      ),
      tags$div(
        style = "width:1px; background:var(--pt-border); align-self:stretch; margin-top:22px;"
      ),
      tags$div(
        style = "text-align:center;",
        tags$div(class = "zoom-group-label",
                 "Zoom out"),
        actionButton("zoom_out_15x", "1.5\u00D7", class = "btn-default"),
        actionButton("zoom_out_3x",  "3\u00D7",   class = "btn-default"),
        actionButton("zoom_out_10x", "10\u00D7",  class = "btn-default")
      )
    )
  })

  output$render_error_ui <- renderUI({
    req(state$plot_error)
    div(class = "alert alert-danger",
        strong("Error: "), state$plot_error,
        br(),
        tags$pre(style="font-size:11px;margin-top:8px;", state$plot_error))
  })

  output$svg_display_ui <- renderUI({
    req(state$svg_path, file.exists(state$svg_path))
    .pt_inline_svg(state$svg_path)
  })

  output$rcode_output <- renderText({
    state$rcode %||% ""
  })

  output$dl_svg <- downloadHandler(
    filename = function() paste0("plottrellis_",
                                 format(Sys.time(), "%Y%m%d_%H%M%S"), ".svg"),
    content  = function(file) {
      req(state$svg_path, file.exists(state$svg_path))
      file.copy(state$svg_path, file)
      .pt_tidy_svg(file)
    }
  )

  output$dl_rcode <- downloadHandler(
    filename = function() paste0("plottrellis_",
                                 format(Sys.time(), "%Y%m%d_%H%M%S"), ".R"),
    content  = function(file) {
      req(state$rcode)
      writeLines(state$rcode, file)
    }
  )

  invisible(NULL)
}
