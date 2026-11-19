# The setup wizard's Export these settings as a samples.txt: one block of
# wizard_server(). Returns list(wiz_collect_flat_settings =,
# wiz_display_settings =): the genome fields, the suggested genes and the
# display settings Save, Export and the dirty guard all write or compare, and
# the display settings alone (`cfg`'s with `import_visual`'s on top).
.pt_wizard_export_server <- function(input, output, data_root, wiz_samples,
                                     import_visual, cfg, DEFAULT_SETTINGS) {
  wiz_collect_flat_settings <- function() {
    c(list(
      gtf_path            = wiz_clean_optional_text(input$wiz_gtf_path),
      orgdb               = wiz_clean_optional_text(input$wiz_orgdb),
      genome_label        = wiz_clean_optional_text(input$wiz_genome_label),
      # no wizard input: the project's, or an imported samples.txt's line
      suggested_genes     = import_visual()$suggested_genes %||%
        cfg$suggested_genes %||% ""
    ), wiz_display_settings())
  }
  # The saved display settings with any imported ones on top (no wizard inputs
  # for them).
  wiz_display_settings <- function() {
    x <- cfg
    iv <- import_visual()
    x[names(iv)] <- iv
    .pt_collect_display_settings(list(), "", x)
  }

  # A separate action from Save, never triggered automatically.
  observeEvent(input$wiz_export_samples, {
    dr <- data_root()
    showModal(modalDialog(
      title = "Export settings to txt", easyClose = TRUE,
      textInput("wiz_export_path", "Write to (full path):",
                value = if (is.null(dr)) "" else file.path(dr, "samples.txt"),
                width = "100%"),
      uiOutput("wiz_export_status"),
      footer = tagList(modalButton("Cancel"),
                       actionButton("wiz_confirm_export", "Export",
                                    class = "pt-once btn-primary"))
    ))
  })

  output$wiz_export_status <- renderUI(NULL)

  observeEvent(input$wiz_confirm_export, {
    path <- trimws(input$wiz_export_path %||% "")
    if (nchar(path) == 0) {
      output$wiz_export_status <- renderUI(div(
        class = "wiz-error", "Enter a path first."
      ))
      return()
    }
    settings <- c(
      list(data_root = data_root(), samples = wiz_samples()),
      wiz_collect_flat_settings()
    )
    ok <- tryCatch({
      writeLines(export_samples_txt(settings, defaults = DEFAULT_SETTINGS),
                 path)
      TRUE
    }, error = function(e) {
      output$wiz_export_status <- renderUI(div(
        class = "wiz-error",
        sprintf("Failed to write: %s", conditionMessage(e))
      ))
      FALSE
    })
    if (ok)
      output$wiz_export_status <- renderUI(div(
        class = "wiz-ok", sprintf("Wrote %s.", path)
      ))
  })

  list(wiz_collect_flat_settings = wiz_collect_flat_settings,
       wiz_display_settings = wiz_display_settings)
}
