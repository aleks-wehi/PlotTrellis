# The setup wizard's page: wizard_ui() and the form wizard_server() renders into
# it (the Genome, Data folder and Samples sections, Save and Export). The page's
# own styles are inst/www/wizard.css.

# wiz_help() (the help button beside a label) lives in app-help-popover.R.

# The messages a failed Save or import shows: one bullet each, never one joined
# line. `kind` "error" is a red alert; "warn" an amber note (an import that
# succeeded with warnings).
wiz_error_list <- function(heading, msgs, kind = "error") {
  div(class = if (identical(kind, "warn"))
                "wiz-warn wiz-warn-box" else "wiz-error wiz-error-box",
      role = if (identical(kind, "warn")) "status" else "alert",
      if (!is.null(heading)) tags$b(heading),
      tags$ul(style = "margin:4px 0 0 0; padding-left:20px;",
              lapply(msgs, tags$li)))
}

# One numbered section card of the setup form (Genome, Data folder, Samples).
wiz_section <- function(number, title, intro, ...) {
  div(class = "wiz-section",
      div(class = "wiz-section-head",
          tags$span(class = "wiz-step", number),
          tags$h3(title)),
      if (!is.null(intro)) p(class = "wiz-section-intro", intro),
      ...)
}

# The "Which GTF?" note under the GTF field: GENCODE or Ensembl, and where to
# download it (a RefSeq GTF finds no gene: its IDs are not Ensembl IDs).
wiz_gtf_guide <- function() {
  link <- function(href, text) tags$a(href = href, target = "_blank", text)
  tags$details(class = "wiz-guide",
    tags$summary("Which GTF?"),
    p(style = "margin:8px 0 0 0;",
      "Must be GENCODE or Ensembl, the comprehensive gene annotation: GENCODE ",
      link("https://www.gencodegenes.org/human/", "human"), " or ",
      link("https://www.gencodegenes.org/mouse/", "mouse"), ", or ",
      link("https://www.ensembl.org/info/data/ftp/index.html", "Ensembl"), ".")
  )
}

# Shared Genome/Data folder/Samples/Save form for every entry path (import,
# manual, edit-mode re-entry; import only prefills it); `cfg` supplies initial
# values (CFG when editing, DEFAULT_SETTINGS when fresh), `is_edit` governs what
# else differs (data-root's blank default, Save's label, the extra Close
# button).
wiz_shared_form <- function(cfg, is_edit) {
  tagList(
    wiz_section(1, "Genome",
      "The gene annotation (GTF)",
      fluidRow(
        column(5, selectizeInput("wiz_genome_preset",
                    wiz_help("Species and genome build",
                             "Fills in the OrgDb and genome label below; the GTF is always your own file."),
                    choices = c(
                      "Choose, or fill the fields below yourself" = "",
                      setNames(names(GENOME_PRESETS),
                               vapply(GENOME_PRESETS, function(g) g$label,
                                      character(1)))
                    ),
                    width = "100%"))
      ),
      textInput("wiz_gtf_path",
                wiz_help("GTF file",
                         "A GENCODE or Ensembl GTF of the genome build your reads were aligned to (RefSeq GTFs do not work)."),
                value = cfg$gtf_path %||% "", width = "100%",
                placeholder = "e.g. /data/reference/gencode.vM32.annotation.gtf.gz"),
      uiOutput("wiz_gtf_status"),
      wiz_gtf_guide(),
      fluidRow(
        column(6, textInput("wiz_orgdb",
                            wiz_help("OrgDb package",
                                     "The installed gene-name package: org.Mm.eg.db (mouse) or org.Hs.eg.db (human)."),
                            value = cfg$orgdb %||% "", width = "100%",
                            placeholder = "org.Mm.eg.db"),
               uiOutput("wiz_orgdb_status")),
        column(6, textInput("wiz_genome_label",
                            wiz_help("Genome label (UCSC name)",
                                     "The UCSC name of the GTF's genome build, e.g. mm10, mm39 or hg38."),
                            value = cfg$genome_label %||% "", width = "100%",
                            placeholder = "mm39"),
               uiOutput("wiz_genome_label_status"))
      )
    ),
    wiz_section(2, "Data folder (full path)", NULL,
      textInput("wiz_data_root", NULL,
                value = if (is_edit) cfg$data_dir else "", width = "100%",
                placeholder = "e.g. /data/smithlab/tcell_2026"),
      uiOutput("wiz_data_root_status")
    ),
    wiz_section(3, "Samples", NULL,
      actionButton("wiz_add_sample", "+ Add a sample",
                   class = "pt-once btn-primary"),
      uiOutput("wiz_samples_ui")
    ),
    uiOutput("wiz_completion_error"),
    div(class = "wiz-actions",
      if (is_edit) tagList(
        actionButton("wiz_save", "Save and return to the app",
                     class = "pt-once btn-success btn-lg"),
        actionButton("wiz_close", "Close without saving",
                     class = "pt-once btn-default btn-lg")
      ) else actionButton("wiz_save", "Save and open the app",
                          class = "pt-once btn-success btn-lg"),
      # Export only, never import, and never triggered by Save: kept apart from
      # the primary buttons and captioned so it cannot be mistaken for a way to
      # load samples.
      div(class = "wiz-export",
          actionButton("wiz_export_samples",
                       "Export these settings as a samples.txt",
                       class = "pt-once btn-link btn-sm"),
          div(class = "text-muted",
              "A text copy you can keep, share, or import next time. Export only."))
    )
  )
}

# The wizard's page for the project whose settings file is (or will be)
# `SETTINGS_JSON`.
wizard_ui <- function(SETTINGS_JSON) {
  is_edit <- file.exists(SETTINGS_JSON)
  # The composer's title bar, full width above the page (.pt_title_bar()). The
  # theme goes in the head, so the bar is never drawn without it. Dark mode as
  # the screen picks it (main.js): ?mode=dark, else the mode last chosen there.
  tagList(tags$head(.pt_theme_css(), tags$script(HTML(
    "(function () { var m = new URLSearchParams(location.search).get('mode');",
    "if (!m) { try { m = localStorage.getItem('pt-mode'); } catch (e) {} }",
    "if (m === 'dark') document.documentElement.setAttribute('data-pt-mode', 'dark'); })();"))),
    .pt_title_bar(), fluidPage(
    class = "wiz-page",
    tags$head(tags$title(if (is_edit) "PlotTrellis: project settings" else
                         "PlotTrellis: set up your project")),
    .pt_once_js(),
    .pt_help_assets(),
    includeCSS(system.file("www", "wizard.css", package = "PlotTrellis",
                           mustWork = TRUE)),
    if (is_edit) tagList(
      h1("Project settings"),
      p(class = "wiz-lead",
        "Change your project's genome, data folder or samples. Save restarts the app with the new settings; Close keeps the old ones.")
    ) else h1("Welcome to PlotTrellis"),
    uiOutput("wiz_body")
  ))
}
