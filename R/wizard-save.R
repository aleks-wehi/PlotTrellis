# The setup wizard's Save: one block of wizard_server(). Checks the fields
# (wiz_save_problems()), writes `SETTINGS_JSON` and the `mirror_txt` copy beside
# it, then stops the app with "restart" so run_plottrellis()'s loop opens the
# composer. A failed Save sets `save_failed` (the "Not saved yet" list then
# follows the fields) or `write_error`. Returns list(save_problems =): what a
# Save now would refuse, as a reactive.
.pt_wizard_save_server <- function(input, output, SETTINGS_JSON, mirror_txt,
                                   data_root, wiz_samples,
                                   wiz_collect_flat_settings, save_failed,
                                   write_error) {
  save_problems <- reactive({
    input$wiz_save  # each Save click re-checks too (see the Save observer)
    wiz_save_problems(data_root(), wiz_samples(), input$wiz_gtf_path,
                      input$wiz_orgdb, input$wiz_genome_label)
  })
  output$wiz_completion_error <- renderUI({
    if (!is.null(write_error()))
      return(wiz_error_list("Not saved:", list(write_error())))
    if (!save_failed()) return(NULL)
    errs <- save_problems()
    if (length(errs) > 0) wiz_error_list("Not saved yet. Please fix:", errs)
  })

  observeEvent(input$wiz_save, {
    write_error(NULL)
    dr <- data_root()
    # Checked afresh, not from save_problems(): a folder created or an OrgDb
    # installed since changes the answer without changing any field.
    errs <- wiz_save_problems(dr, wiz_samples(), input$wiz_gtf_path,
                              input$wiz_orgdb, input$wiz_genome_label)
    if (length(errs) > 0) {
      save_failed(TRUE)
      return()
    }

    settings <- c(
      list(
        data_dir        = dr,
        samples         = wiz_samples()
      ),
      wiz_collect_flat_settings()
    )
    # Kept across edits: pt_project() uses it to tell this project's samples.txt
    # from a different one pointed at the same folder.
    if (file.exists(SETTINGS_JSON)) {
      old <- tryCatch(jsonlite::read_json(SETTINGS_JSON),
                      error = function(e) list())
      settings$samples_txt_source <- old$samples_txt_source
    }
    # Always keep a plain-text project_samples.txt beside it, the human-readable
    # record the project can be rebuilt from (mirror_txt). A folder that can't
    # be written (permissions, a full disk) says so here instead of ending the
    # session.
    txt_ok <- tryCatch(.pt_write_project_files(settings, SETTINGS_JSON),
                       error = function(e) {
      write_error(paste0("Couldn't write ", SETTINGS_JSON, " (",
                         conditionMessage(e), "). ",
                         "Check that you can write to that folder, then Save again."))
      NULL
    })
    if (is.null(txt_ok)) return()

    save_failed(FALSE)
    showModal(modalDialog(
      title = "Settings saved", easyClose = FALSE, footer = NULL,
      p("project_settings.json has been written."),
      p(if (txt_ok) paste0("A copy of your samples is in ", mirror_txt,
                           "; it can rebuild the project if the settings ever break.")
        else "Couldn't write project_samples.txt beside it (use Export samples.txt to save a copy yourself)."),
      p(strong("You can close this tab now."),
        " The app reopens with your project in a new tab from RStudio in a few seconds. ",
        "If nothing opens, go back to RStudio: the app's address is in the Console.")
    ))
    message("PlotTrellis: settings saved, launching the app.")
    later::later(function() stopApp(returnValue = "restart"), 0.2)
  })

  list(save_problems = save_problems)
}
