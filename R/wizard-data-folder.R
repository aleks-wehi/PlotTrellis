# The setup wizard's Data folder field: one block of wizard_server(). Returns
# list(data_root =): the typed folder as a reactive, NULL when blank.
.pt_wizard_data_folder_server <- function(input, output) {
  data_root <- reactive({
    dr <- trimws(input$wiz_data_root %||% "")
    if (nchar(dr) == 0) NULL else dr
  })

  output$wiz_data_root_status <- renderUI({
    dr <- data_root()
    if (is.null(dr)) return(NULL)
    if (dir.exists(dr)) div(class = "wiz-status wiz-ok", "\u2713 Folder found.")
    else div(class = "wiz-status wiz-error",
             "\u2717 No folder at this path (it must be a path on the computer running R).")
  })

  list(data_root = data_root)
}
