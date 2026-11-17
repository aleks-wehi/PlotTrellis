# The setup wizard's Import a samples.txt box: one block of wizard_server(). An
# import fills the committed samples (`wiz_samples`), the data folder and genome
# fields, and keeps the file's lines (`imported`, for its start view) and its
# display settings (`import_visual`) for Save. `DEFAULT_SETTINGS` says which
# `# key: value` lines are settings.
.pt_wizard_import_server <- function(input, output, session, wiz_samples,
                                     imported, import_visual,
                                     DEFAULT_SETTINGS) {
  output$wiz_import_status <- renderUI(NULL)

  observeEvent(input$wiz_import_load, {
    path <- trimws(input$wiz_import_path %||% "")
    if (nchar(path) == 0) {
      output$wiz_import_status <- renderUI(div(
        class = "wiz-error", "Enter a path to samples.txt first."
      ))
      return()
    }
    if (!is_real_file(path)) {
      output$wiz_import_status <- renderUI(div(
        class = "wiz-error", sprintf("No such file: %s", path)
      ))
      return()
    }
    res <- tryCatch(
      # defaults = DEFAULT_SETTINGS: lets any `# <key>: value` line in the
      # file be parsed back into `res$overrides`, applied below: the 3 genome
      # fields (gtf_path/orgdb/genome_label) fill their inputs, the display
      # settings are kept for Save. A line for a key that isn't a setting
      # (e.g. a cache path like `# txdb_path: ...`) is ignored.
      import_samples_txt(path, defaults = DEFAULT_SETTINGS),
      error = function(e) list(errors = conditionMessage(e))
    )
    if (length(res$errors) > 0) {
      output$wiz_import_status <- renderUI(
        wiz_error_list("Nothing imported. The file has these problems:",
                       res$errors)
      )
      return()
    }
    if (!is.null(res$data_root))
      updateTextInput(session, "wiz_data_root", value = res$data_root)
    wiz_samples(res$samples)
    imported(readLines(path, warn = FALSE))

    # Apply any `# <key>: value` overrides: the genome fields onto their inputs,
    # the display settings into import_visual() for Save. Numbers for the
    # NULL-default ("auto") keys are converted the way parse_samples_txt() does.
    coerced <- .pt_coerce_auto_overrides(res$overrides, DEFAULT_SETTINGS)
    ov <- coerced$overrides
    res$warnings <- c(res$layout_warning, res$warnings, coerced$ignored)
    visual_keys <- setdiff(names(DEFAULT_SETTINGS), .PT_CORE_SETTINGS)
    import_visual(ov[intersect(names(ov), visual_keys)])
    if (!is.null(ov$gtf_path))
      updateTextInput(session, "wiz_gtf_path", value = ov$gtf_path)
    if (!is.null(ov$orgdb))
      updateTextInput(session, "wiz_orgdb", value = ov$orgdb)
    if (!is.null(ov$genome_label))
      updateTextInput(session, "wiz_genome_label", value = ov$genome_label)

    status_lines <- list(div(class = "wiz-ok",
                             sprintf("Imported %d sample(s).",
                                     length(res$samples))))
    if (length(res$warnings) > 0)
      status_lines <- c(status_lines,
                        list(wiz_error_list("Imported, with warnings:",
                                            res$warnings, kind = "warn")))
    output$wiz_import_status <- renderUI(tagList(status_lines))
  })
}
