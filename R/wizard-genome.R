# The setup wizard's Genome section: one block of wizard_server(). The species
# preset fills the OrgDb and genome label fields, and the live status lines
# under the GTF, OrgDb and genome label fields say what the checks in
# wizard-checks.R find.

# One status line under a field, for a check's list(state =, text =); NULL for
# "blank".
.pt_wizard_status_line <- function(chk) {
  if (identical(chk$state, "blank")) return(NULL)
  mark <- switch(chk$state, ok = "\u2713 ", warn = "\u26a0 ", "\u2717 ")
  div(class = paste("wiz-status", switch(chk$state, ok = "wiz-ok",
                                         warn = "wiz-warn", "wiz-error")),
      paste0(mark, chk$text))
}

# The species preset and the three status lines. Returns list(gtf_check =): the
# GTF check, as a reactive.
.pt_wizard_genome_server <- function(input, output, session) {
  # Quick-fill the genome identity fields from a bundled GENOME_PRESETS
  # entry -- pure convenience, nothing is stored under the preset's id; Save
  # always writes the plain field values. GENOME_PRESETS carries no
  # gtf_path, so the GTF path field is left alone: always the user's own
  # file, and a path already entered survives a preset pick.
  observeEvent(input$wiz_genome_preset, {
    preset <- GENOME_PRESETS[[input$wiz_genome_preset]]
    if (is.null(preset)) return()
    updateTextInput(session, "wiz_orgdb", value = preset$orgdb)
    updateTextInput(session, "wiz_genome_label", value = preset$genome_label)
  }, ignoreInit = TRUE)

  # Live checks under the GTF and OrgDb fields (wiz_gtf_check(),
  # wiz_orgdb_check()); advice, never a block: Save's own checks are
  # wiz_save_problems(). The typed path is debounced, so the file is read once
  # per pause, not per keystroke.
  gtf_path_settled <- debounce(reactive(input$wiz_gtf_path), 400)
  gtf_check <- reactive(wiz_gtf_check(gtf_path_settled()))
  output$wiz_gtf_status   <- renderUI(.pt_wizard_status_line(gtf_check()))
  output$wiz_orgdb_status <- renderUI(.pt_wizard_status_line(
    wiz_orgdb_check(input$wiz_orgdb, gtf_check()$species)
  ))
  output$wiz_genome_label_status <- renderUI(.pt_wizard_status_line(
    wiz_genome_label_check(input$wiz_genome_label, gtf_check()$build,
                           gtf_check()$species)
  ))
  # The genome label's placeholder names the GTF's build (mm10 for GRCm38), not
  # a fixed mm39.
  observe({
    build <- gtf_check()$build
    hint <- if (!is.null(build) && build %in% names(WIZ_LABEL_BY_BUILD))
              WIZ_LABEL_BY_BUILD[[build]]
            else if (identical(gtf_check()$species$species, "human")) "hg38"
            else "mm39"
    updateTextInput(session, "wiz_genome_label", placeholder = hint)
  })

  list(gtf_check = gtf_check)
}
