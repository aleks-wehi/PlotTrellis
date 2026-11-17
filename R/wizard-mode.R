# The setup wizard's landing screen and the switch from it to the form: one
# block of wizard_server(), a top-level function taking the session objects and
# the state it reads. `wiz_mode` is NULL on the landing screen, "import" /
# "manual" for the entry path chosen and "edit" for an existing project; every
# non-NULL mode shows wiz_shared_form() for `cfg`, below the import box in
# "import" mode (prefilled with `mirror_txt` when `has_mirror`). Back to the
# start clears a failed Save's messages (`save_failed`, `write_error`).
.pt_wizard_mode_server <- function(input, output, wiz_mode, save_failed,
                                   write_error, cfg, is_edit,
                                   has_mirror, mirror_txt) {
  observeEvent(input$wiz_mode_import, wiz_mode("import"))
  observeEvent(input$wiz_mode_manual, wiz_mode("manual"))
  observeEvent(input$wiz_mode_back,
               { wiz_mode(NULL); save_failed(FALSE); write_error(NULL) },
               ignoreInit = TRUE)

  output$wiz_body <- renderUI({
    mode <- wiz_mode()
    if (is.null(mode)) {
      need <- function(title, ...) div(class = "wiz-need", tags$b(title), ...)
      return(tagList(
        h3("What a project needs", style = "margin-bottom:12px;"),
        div(class = "wiz-needs",
          need("Your data files",
               "Coverage (bigWig), Hi-C maps, BED intervals, splice-junction or loop files, salmon quant folders. ",
               "In one folder on the computer running R."),
          need("A GTF",
               "The gene annotation from GENCODE or Ensembl, for the same genome build and release your data was quantified with. ",
               "Not RefSeq. The next screen says where to get it."),
          need("The species' annotation package",
               "org.Mm.eg.db (mouse) or org.Hs.eg.db (human), installed in R. Gene search uses it to find genes by name.")
        ),
        h3("How do you want to start?", style = "margin-bottom:12px;"),
        div(class = "wiz-choices",
          div(class = "wiz-choice",
              h3("Fill in a form"),
              p("Type the genome and data folder, then add each sample and its files from a list. Best for a handful of samples."),
              actionButton("wiz_mode_manual", "Set up by hand",
                           class = "btn-primary btn-lg")),
          div(class = "wiz-choice",
              h3("Import a samples.txt"),
              p("A tab-separated table listing your samples, written by you, a colleague or a chat AI (the vignette shows how). It prefills the form."),
              actionButton("wiz_mode_import", "Import a samples.txt",
                           class = "btn-default btn-lg"))
        )
      ))
    }
    tagList(
      if (!identical(mode, "edit"))
        actionLink("wiz_mode_back", "\u2190 Back to the start",
                   class = "wiz-back"),
      if (identical(mode, "import")) div(class = "wiz-section",
        div(class = "wiz-section-head", tags$h3("Import a samples.txt")),
        if (has_mirror)
          p(class = "wiz-section-intro",
            "A project_samples.txt saved with your last project was found. Click Import to restore it.")
        else
          p(class = "wiz-section-intro",
            "The file fills in the three sections below; check them and add anything missing before saving."),
        fluidRow(
          column(9, textInput(
            "wiz_import_path",
            "Path to the samples.txt (on the computer running R)",
            value = if (has_mirror) mirror_txt else "", width = "100%",
            placeholder = "e.g. /data/smithlab/tcell_2026/samples.txt"
          )),
          column(3, style = "margin-top: 25px;",
                 actionButton("wiz_import_load", "Import",
                              class = "pt-once btn-primary"))
        ),
        uiOutput("wiz_import_status")
      ),
      wiz_shared_form(cfg, is_edit)
    )
  })
}
