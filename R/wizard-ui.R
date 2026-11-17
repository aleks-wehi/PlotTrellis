# The setup wizard's page: wizard_ui() and the form wizard_server() renders into
# it (the Genome, Data folder and Samples sections, Save and Export). The page's
# own styles are inst/www/wizard.css.

# wiz_help() (the help button beside a label) lives in help_popover.R.

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

# The "Which GTF?" guide under the GTF field: what works and why, where to get
# it, and the OrgDb package and genome label that go with each. Every claim here
# is checked against ensure_genome_ready() (Ensembl IDs required, UCSC
# chromosome names applied, lengths fetched by genome_label); a plain Ensembl
# GTF was run end to end before the Ensembl claims were written.
wiz_gtf_guide <- function() {
  row <- function(...) tags$tr(lapply(list(...), tags$td))
  tags$details(class = "wiz-guide",
    tags$summary("Which GTF? Where do I get it?"),
    tags$ul(
      tags$li(tags$b("GENCODE or Ensembl, not RefSeq. "),
              "PlotTrellis finds a gene by looking its name up in the OrgDb package, which gives the ",
              "gene's Ensembl ID (ENSG..., ENSMUSG...), then finds that ID in the GTF. GENCODE and ",
              "Ensembl GTFs use Ensembl IDs; RefSeq (NCBI) GTFs use other IDs, so no gene would be ",
              "found. GENCODE gives the most: an Ensembl GTF draws every track, but has no APPRIS labels ",
              "for the Gene info table. Mouse and human are supported."),
      tags$li(tags$b("The release your quantification used. "),
              "Transcripts in your quant.sf files are matched to the GTF by ID, so use the GTF from ",
              "the same GENCODE or Ensembl release as the transcript sequences Salmon was run on; a ",
              "transcript the GTF does not have is left out of the expression tracks. For GENCODE, ",
              "Salmon must have been run with --gencode, or no transcript name matches. A different ",
              "genome build (mm10 vs mm39) puts every gene at the wrong coordinates."),
      tags$li(tags$b("Keep it gzipped "),
              "(.gtf.gz). The first launch builds an index from it, ",
              "about 1 to 3 minutes for a whole genome, once per GTF file.")
    ),
    p(tags$b("Where to download it: "),
      tags$a(href = "https://www.gencodegenes.org/human/", target = "_blank",
             "GENCODE human"), " or ",
      tags$a(href = "https://www.gencodegenes.org/mouse/", target = "_blank",
             "GENCODE mouse"),
      ", the \"Comprehensive gene annotation\" GTF (region CHR); earlier releases under \"Releases\" on the same page. ",
      "Ensembl: ",
      tags$code("https://ftp.ensembl.org/pub/release-<N>/gtf/<species>/"),
      ", the file ",
      tags$code("<Species>.<build>.<N>.gtf.gz"), "; human GRCh37 is under ",
      tags$code("pub/grch37/"), " instead."),
    tags$table(class = "table table-condensed wiz-guide-table",
      tags$thead(tags$tr(tags$th("Species, build"), tags$th("GENCODE releases"),
                         tags$th("Ensembl releases"), tags$th("OrgDb package"),
                         tags$th("Genome label"))),
      tags$tbody(
        row("Mouse, GRCm39", "M26 and later", "103 and later",
            tags$code("org.Mm.eg.db"), tags$code("mm39")),
        row("Mouse, GRCm38", "M25 and earlier", "68 to 102",
            tags$code("org.Mm.eg.db"), tags$code("mm10")),
        row("Human, GRCh38", "20 and later", "76 and later",
            tags$code("org.Hs.eg.db"), tags$code("hg38")),
        row("Human, GRCh37", "3c to 19 (not the \"lift37\" files)",
            "55 to 75, or the GRCh37 archive", tags$code("org.Hs.eg.db"),
            tags$code("hg19"))
      )
    ),
    p(class = "text-muted", style = "margin-bottom:0;",
      "The build must be the one your reads were aligned to. GENCODE's \"lift37\" GTFs carry IDs such as ",
      tags$code("ENSG00000290825.1_1"),
      " that match no gene name: use GENCODE 19 instead. ",
      "Install the OrgDb package once in R: ",
      tags$code("BiocManager::install(\"org.Mm.eg.db\")"),
      " (mouse) or ", tags$code("BiocManager::install(\"org.Hs.eg.db\")"),
      " (human).")
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
      "The gene annotation (GTF) PlotTrellis draws genes and transcripts from, and the R package it looks gene names up in.",
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
                wiz_help("GTF file (path on the computer running R)",
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
    wiz_section(2, "Data folder",
      "The folder your data files are in (bigWig, Hi-C, BED, junction or loop files, salmon quant folders), on the computer running R. Each sample's files are then picked from it, so one folder for everything is easiest.",
      textInput("wiz_data_root",
                wiz_help("Data folder (full path)",
                         "Files in any subfolder of it are found."),
                value = if (is_edit) cfg$data_dir else "", width = "100%",
                placeholder = "e.g. /data/smithlab/tcell_2026"),
      uiOutput("wiz_data_root_status")
    ),
    wiz_section(3, "Samples",
      "A sample is one condition you want to look at (Control, KO, Basal, LP...), with its files. Add one per condition; the figure then draws them side by side or compares two of them.",
      actionButton("wiz_add_sample", "+ Add a sample",
                   class = "pt-once btn-primary"),
      uiOutput("wiz_samples_ui")
    ),
    p(class = "text-muted",
      "Fonts, track heights and other display settings are changed in the app's Figure settings panel, not here."),
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
  # theme goes in the head, so the bar is never drawn without it.
  tagList(tags$head(.pt_theme_css()), .pt_title_bar(), fluidPage(
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
    ) else tagList(
      h1("Welcome to PlotTrellis"),
      p(class = "wiz-lead",
        "Let's set up your project: where your data is, which genome annotation to use, and what your samples are. ",
        "It takes a few minutes, and everything can be changed later from the app.")
    ),
    uiOutput("wiz_body")
  ))
}
