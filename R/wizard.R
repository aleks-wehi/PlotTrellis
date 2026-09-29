# rationale: see solutions/concepts/19-relocated-code-rationale.md, "utils/wizard.R:1"

WIZARD_DEFAULT_COLOUR <- "#1f78b4"  # PALETTE[1] -- see app.R; not re-cycled per sample, a minor known simplification

WIZARD_DATA_TYPES <- list(
  bw          = list(label = "BigWig coverage (ATAC / ChIP-seq / CUT&RUN / RNA-seq)", ext = "\\.(bw|bigWig)$"),
  hic         = list(label = "Hi-C contact map",                                     ext = "\\.(hic|cool|mcool)$"),
  bed         = list(label = "BED intervals",                                       ext = "\\.(bed|narrowPeak|broadPeak|gappedPeak)$"),
  `paired-end` = list(label = "Paired-interval arcs (BED12 splice junctions or BEDPE loop calls)",
                       ext  = "\\.(bed|bed12|bedpe)$")
)

# Returns paths relative to data_root (for populating an instance's path selectize); quantdir instead lists directories containing quant.sf, not extension-matched files.
wizard_scoped_choices <- function(data_root, type) {
  if (is.null(data_root) || !nchar(data_root) || !dir.exists(data_root)) return(character(0))
  if (identical(type, "quantdir")) {
    quant_files <- list.files(data_root, pattern = "^quant\\.sf$", recursive = TRUE, full.names = TRUE)
    dirs <- unique(dirname(quant_files))
    return(substring(dirs, nchar(data_root) + 2))
  }
  ext <- WIZARD_DATA_TYPES[[type]]$ext
  files <- list.files(data_root, pattern = ext, recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  substring(files, nchar(data_root) + 2)
}

# A literal "\u00b7" not surrounded by spaces breaks format_track_label()'s bolding regex (Phase C.0 open issue #2) -- reject at entry time.
wizard_name_has_bad_dot <- function(name) {
  grepl("\u00b7", name, fixed = TRUE) && !grepl("\\s\u00b7\\s", name)
}

# Treats NULL/NA/non-finite/list alike as unset (list too, for jsonlite's NULL/{} round-trip -- see app.R's copy); duplicated from app.R's clean_optional_numeric() since this file's separate shinyApp() never reaches app.R's main-app-only copy; used for zrange_cap/hic_height.
wiz_clean_optional_numeric <- function(v) {
  if (is.null(v) || is.list(v) || is.na(v) || !is.finite(v)) NULL else v
}

# Same idea as wiz_clean_optional_numeric(): normalizes a blank input to NULL, not "", matching how an unbuilt path is represented elsewhere (CFG, GENOME_PRESETS) and avoiding a spurious is_dirty()/export delta trip against NULL.
wiz_clean_optional_text <- function(v) {
  v <- trimws(v %||% "")
  if (nchar(v) == 0) NULL else v
}

# A native browser tooltip (title attribute), not a JS/CSS framework widget -- used wherever a field's purpose isn't obvious from its name alone.
wiz_help <- function(label_text, help_text) {
  tagList(label_text, tags$span(
    "ⓘ", title = help_text, class = "wiz-help-icon",
    style = "cursor:help; color:#888; margin-left:4px; font-size:12px;"
  ))
}

# Shared Genome/Data-root/Samples/Other-settings/Save form for every entry path (import, manual, edit-mode re-entry); `cfg` supplies initial values (CFG when editing, DEFAULT_SETTINGS when fresh), `is_edit` governs what else differs (data-root's blank default, Save's label, the extra Close button). See solutions/concepts/19-relocated-code-rationale.md, "utils/wizard.R:66".
wiz_shared_form <- function(cfg, is_edit) {
  tagList(
    # Deliberately at the very top with muted styling and an explicit
    # "does not import" caption, so it can't be mistaken for a way to load
    # sample data -- it only ever writes the current in-wizard settings to a
    # file, and is never triggered by Save.
    div(class = "wiz-export-bar",
        style = "text-align:right; padding:6px 0; margin-bottom:12px; border-bottom:1px solid #ddd;",
        actionButton("wiz_export_samples", "⤓ Export settings to txt", class = "btn-outline-secondary btn-sm"),
        tags$div(class = "text-muted", style = "font-size:11px; margin-top:2px;",
                 "Writes the current settings below to a samples.txt file. Export only, does not import.")
    ),
    wellPanel(
      h4("Genome"),
      p(class = "text-muted", "GTF path and OrgDb package are required; TxDb and the two transcript-lookup caches are built automatically the first time this genome is used, and are never shown here -- they're pure derived files, not settings."),
      selectizeInput("wiz_genome_preset",
                  wiz_help("Quick fill from a bundled preset:",
                            "Pre-fills the fields below from a genome already set up on this server. You can still edit them afterward."),
                  choices = c("(choose one)" = "",
                              setNames(names(GENOME_PRESETS), vapply(GENOME_PRESETS, function(g) g$label, character(1))))),
      fluidRow(
        column(6, textInput("wiz_gtf_path",
                             wiz_help("GTF path:", "Path to a gene annotation file (.gtf or .gtf.gz), e.g. downloaded from GENCODE or Ensembl. Defines every gene/transcript this project can show."),
                             value = cfg$gtf_path %||% "", width = "100%")),
        column(3, textInput("wiz_orgdb",
                             wiz_help("OrgDb package:", "The Bioconductor annotation package for this species' gene symbols, e.g. org.Mm.eg.db (mouse) or org.Hs.eg.db (human). Must already be installed."),
                             value = cfg$orgdb %||% "")),
        column(3, textInput("wiz_genome_label",
                             wiz_help("Genome label:", "Must be the exact UCSC genome-build id for this genome, e.g. hg38, mm10, or mm39 -- used both as a display label and to look up real chromosome lengths from UCSC (the GTF itself doesn't contain them). Falls back to the OrgDb package name if left blank, which will not be a valid UCSC id."),
                             value = cfg$genome_label %||% ""))
      )
    ),
    wellPanel(
      h4("Data root"),
      p("Every sample file path below is entered relative to this directory."),
      textInput("wiz_data_root", "Data root directory (absolute path):",
                value = if (is_edit) cfg$data_dir else "", width = "100%"),
      uiOutput("wiz_data_root_status")
    ),
    wellPanel(
      h4("Samples"),
      actionButton("wiz_add_sample", "+ Add sample", class = "btn-primary"),
      hr(),
      uiOutput("wiz_samples_ui")
    ),
    wellPanel(
      tags$details(
        tags$summary(tags$b("Advanced settings")),
        p(class = "text-muted", style = "margin-top:8px;", "Every value below has a sensible default and can be changed later by re-opening this wizard."),

        h5("Layout & display"),
        fluidRow(
          column(3, selectizeInput("wiz_font", "Font:", choices = c("DejaVu Sans", "Helvetica", "Arial"), selected = cfg$font %||% "DejaVu Sans")),
          column(3, numericInput("wiz_fontsize", "Font size:", value = cfg$fontsize %||% 8, min = 4, max = 24)),
          column(3, numericInput("wiz_page_w",
                                  wiz_help("Page width (in):", "Width of the exported figure, in inches. Track heights below adapt to whatever's left after labels and margins."),
                                  value = cfg$page_w %||% 9, min = 3, max = 30)),
          column(3, selectizeInput("wiz_chr_fmt",
                                 wiz_help("Chromosome format:", "How chromosome names are displayed on the figure's axis."),
                                 choices = c("Dot (chr.19)" = "dot", "Space (chr 19)" = "space", "None (chr19)" = "plain"),
                                 selected = cfg$chr_fmt %||% "dot")),
          column(3, checkboxInput("wiz_show_track_axes",
                                   wiz_help("Show track value axes", "Draws a numeric value axis on coverage (bigWig) and junction/sashimi tracks. On by default."),
                                   value = cfg$show_track_axes %||% TRUE))
        ),

        hr(),
        h5("Coverage tracks (BigWig: ATAC / ChIP-seq / CUT&RUN / RNA-seq)"),
        fluidRow(
          column(3, numericInput("wiz_atac_height",
                                  wiz_help("Track height (in):", "Height in inches of each BigWig coverage track (ATAC, ChIP-seq, CUT&RUN, RNA-seq)."),
                                  value = cfg$atac_height %||% 0.85, min = 0.3, max = 8, step = 0.1))
        ),

        hr(),
        h5("Hi-C"),
        fluidRow(
          column(3, selectizeInput("wiz_hic_style",
                                 wiz_help("Display style:", "How the Hi-C contact matrix is drawn: a full rectangle, or a triangle (upper-left half)."),
                                 choices = c("Rectangle" = "rectangle", "Triangle" = "triangle"), selected = cfg$hic_style %||% "rectangle")),
          column(3, textInput("wiz_hic_norm",
                               wiz_help("Normalisation:", "Hi-C matrix normalisation method understood by the .hic/.cool/.mcool file, e.g. KR, VC, SCALE. Leave as KR unless your data needs a different one."),
                               value = cfg$hic_norm %||% "KR")),
          column(3, numericInput("wiz_zrange_cap",
                                  wiz_help("Colour scale cap:", "Caps the Hi-C colour scale at this contact-count value. Leave blank to scale automatically per plot."),
                                  value = wiz_clean_optional_numeric(cfg$zrange_cap) %||% NA, min = 1, step = 10)),
          column(3, numericInput("wiz_hic_height",
                                  wiz_help("Track height (in):", "Height in inches of each Hi-C track. Leave blank to size automatically based on the display style."),
                                  value = wiz_clean_optional_numeric(cfg$hic_height) %||% NA, min = 0.5, max = 15, step = 0.5))
        ),
        fluidRow(
          column(3, selectizeInput("wiz_hic_res_mode",
                                 wiz_help("Resolution mode:", "Auto picks the coarsest resolution that keeps the view under the bin cap; Finest always uses the highest-resolution matrix available (slower for wide views)."),
                                 choices = c("Auto" = "auto", "Finest" = "finest"), selected = cfg$hic_res_mode %||% "auto")),
          column(3, numericInput("wiz_hic_max_bins_per_view",
                                  wiz_help("Max bins per view:", "Upper limit on how many Hi-C bins are drawn at once -- used by Auto resolution mode to avoid overly dense plots."),
                                  value = cfg$hic_max_bins_per_view %||% 500, min = 50, max = 5000, step = 50))
        ),

        hr(),
        h5("Junction / arc tracks"),
        fluidRow(
          column(4, checkboxInput("wiz_arc_clip_partial",
                                   wiz_help("Clip partial arcs", "When an arc's anchors straddle the edge of the plotted region (one anchor visible, the other outside it), truncate the arc at the boundary instead of drawing it extending past the edge. Off by default. Long-range BEDPE-sourced pairs (e.g. Hi-C loop calls) are more likely to need this than short RNA-seq splice junctions -- see the BED12/BEDPE intake notes for why."),
                                   value = cfg$arc_clip_partial %||% FALSE))
        ),

        hr(),
        h5("Comparison (diff) tracks"),
        fluidRow(
          column(3, tags$label("Downregulated colour:", title = "Colour for the most downregulated end of a comparison track's scale."), br(),
                 tags$input(type = "color", id = "wiz_diff_colour_lo", value = cfg$diff_colour_lo %||% "#0033CC",
                            style = "width:60px; height:34px; padding:2px; border:1px solid #ccc; border-radius:4px; cursor:pointer; display:block;",
                            oninput = "Shiny.setInputValue('wiz_diff_colour_lo', this.value, {priority: 'event'});")),
          column(3, tags$label("Upregulated colour:", title = "Colour for the most upregulated end of a comparison track's scale."), br(),
                 tags$input(type = "color", id = "wiz_diff_colour_hi", value = cfg$diff_colour_hi %||% "#CC0000",
                            style = "width:60px; height:34px; padding:2px; border:1px solid #ccc; border-radius:4px; cursor:pointer; display:block;",
                            oninput = "Shiny.setInputValue('wiz_diff_colour_hi', this.value, {priority: 'event'});")),
          column(3, numericInput("wiz_junction_opacity_cap",
                                  wiz_help("Junctions opacity cap (log2FC):", "log2 fold-change value at which a junction arc reaches full opacity in a comparison track. Smaller values make modest changes more visible."),
                                  value = cfg$junction_opacity_cap %||% 2, min = 0.1, step = 0.5)),
          column(3, numericInput("wiz_min_junction_support",
                                  wiz_help("Min arc support:", "Minimum pooled nSupport for a paired-end arc (BED12 splice junction or BEDPE loop call) to be drawn. Default 0 (show everything) -- a BEDPE score column's scale is caller-specific (e.g. a HiChIP loop caller's continuous significance score, not a read count), so raise this only once you know your own data's scale."),
                                  value = cfg$min_junction_support %||% 0, min = 0, step = 1))
        ),
        fluidRow(
          column(3, tags$label("Reference colour:", title = "Fixed colour always used for the reference condition in a 2-sample comparison track, regardless of that condition's own colour."), br(),
                 tags$input(type = "color", id = "wiz_reference_colour", value = cfg$reference_colour %||% "#4D4D4D",
                            style = "width:60px; height:34px; padding:2px; border:1px solid #ccc; border-radius:4px; cursor:pointer; display:block;",
                            oninput = "Shiny.setInputValue('wiz_reference_colour', this.value, {priority: 'event'});")),
          column(3, tags$label("Comparison colour:", title = "Fixed colour always used for the comparison condition in a 2-sample comparison track, regardless of that condition's own colour."), br(),
                 tags$input(type = "color", id = "wiz_comparison_colour", value = cfg$comparison_colour %||% "#56B4E9",
                            style = "width:60px; height:34px; padding:2px; border:1px solid #ccc; border-radius:4px; cursor:pointer; display:block;",
                            oninput = "Shiny.setInputValue('wiz_comparison_colour', this.value, {priority: 'event'});"))
        ),

        hr(),
        h5("Gene view"),
        fluidRow(
          column(3, numericInput("wiz_gene_flank_pct",
                                  wiz_help("Flanking (%):", "Percent of the gene's length added to each side when jumping to a gene, so nearby regulatory regions stay in view."),
                                  value = cfg$gene_flank_pct %||% 15, min = 0, max = 200, step = 5)),
          column(3, numericInput("wiz_min_display_tpm",
                                  wiz_help("Min display TPM:", "Per-condition mean-TPM floor a transcript must clear to be individually named/displayed; below it, a transcript folds into 'Other'. Default 1."),
                                  value = cfg$min_display_tpm %||% 1, min = 0, step = 0.5)),
          column(4, checkboxInput("wiz_uniprot_reviewed_only",
                                   wiz_help("UniProt: reviewed only", "Restrict the Gene info panel's UniProt link-out to manually reviewed (Swiss-Prot) entries only."),
                                   value = cfg$uniprot_reviewed_only %||% FALSE))
        )
      )
    ),
    uiOutput("wiz_completion_error"),
    if (is_edit) tagList(
      actionButton("wiz_save", "Save", class = "btn-success btn-lg"),
      actionButton("wiz_close", "Close without saving", class = "btn-default btn-lg", style = "margin-left:10px;")
    ) else actionButton("wiz_save", "Save and launch", class = "btn-success btn-lg")
  )
}

wizard_ui <- function() {
  is_edit <- file.exists(SETTINGS_JSON)
  fluidPage(
    tags$head(tags$title(if (is_edit) "PlotTrellis: edit project settings" else "PlotTrellis: first-time setup")),
    tags$style(HTML("
      .wiz-row { padding: 4px 0; border-bottom: 1px solid #eee; }
      .wiz-row-error { background: #fdecea; }
      .wiz-error { color: #c62828; margin-bottom: 8px; }
      .wiz-ok    { color: #2e7d32; }
      /* Add-sample modal (modalDialog size='l') is Bootstrap's fixed 900px
         .modal-lg -- widen it so the instance-row selectize boxes have more
         room. #shiny-modal is the one live modal at a time, so this only
         affects a modal-lg while it's showing. */
      #shiny-modal .modal-lg { width: 1100px; max-width: 95vw; }
      /* selectize's single-select box doesn't truncate long selected paths by
         default -- they overflow past the visible input. Clip with ellipsis
         instead. */
      .selectize-input, .selectize-input > input, .selectize-input > .item {
        overflow: hidden;
        text-overflow: ellipsis;
        white-space: nowrap;
      }
      /* The Genome/Advanced-settings collapsible sections (tags$details) --
         summary is the only clickable toggle, so it needs to look like one. */
      details > summary { cursor: pointer; font-weight: 600; margin-bottom: 4px; }
      details[open] > summary { margin-bottom: 10px; }
      .wiz-help-icon { user-select: none; }
    ")),
    h2(if (is_edit) "PlotTrellis: edit project settings" else "PlotTrellis: first-time project setup"),
    p(if (is_edit)
        "Editing your existing project configuration. Save to restart back into the app with the new settings, or close without saving to discard any changes."
      else
        "No project_settings.json found at this deployment's expected path. Let's configure your data before launching the app."),
    uiOutput("wiz_body")
  )
}

wizard_server <- function(input, output, session) {

  is_edit <- file.exists(SETTINGS_JSON)
  cfg     <- if (is_edit) CFG else DEFAULT_SETTINGS

  wiz_samples     <- reactiveVal(if (is_edit) cfg$samples else list())  # committed samples: list(condition=,colour=,data=list(list(type=,name=,path=)))
  modal_instances <- reactiveVal(list())  # bw/hic/bed/paired-end rows in the currently-open add-sample modal
  modal_row_id    <- reactiveVal(0)
  # NULL = the modal is adding a brand-new sample; an integer = the 1-based
  # index (into wiz_samples()) of the sample currently being edited, so
  # wiz_confirm_sample can replace it in place instead of appending a
  # duplicate. Set by wiz_add_sample/wiz_edit_sample below, cleared on
  # confirm -- a re-click of either always resets it, so a stale value can't
  # survive a Cancel.
  modal_edit_index <- reactiveVal(NULL)
  # Per-row/per-field error messages from the last failed wiz_confirm_sample
  # attempt, so a validation error can highlight the specific offending
  # row/file instead of only listing errors in one combined banner. Keyed by
  # row id (character) for modal_row_errors, by "quantdir" for
  # modal_field_errors. Both reset to list() on a genuinely fresh modal open.
  modal_row_errors   <- reactiveVal(list())
  modal_field_errors <- reactiveVal(list())
  # NULL = landing screen (choose a path); "import" / "manual" = which entry
  # path was chosen; "edit" = re-entered from the running app to edit an
  # existing project (Phase C.3), which skips the landing screen entirely --
  # there's no "start fresh vs. import" choice to make when a project already
  # exists. Every non-NULL mode renders the same wiz_shared_form() below an
  # optional import box -- see wiz_shared_form()'s comment.
  wiz_mode <- reactiveVal(if (is_edit) "edit" else NULL)

  observeEvent(input$wiz_mode_import, wiz_mode("import"))
  observeEvent(input$wiz_mode_manual, wiz_mode("manual"))
  observeEvent(input$wiz_mode_back, wiz_mode(NULL), ignoreInit = TRUE)

  # Snapshot of every field this form can change, taken once at edit-mode
  # entry -- the baseline is_dirty() compares live values against, for the
  # Close-without-saving confirm guard (Phase C.3 steps 12/13). Not used
  # (and not needed) in fresh-setup mode: there's no prior configuration to
  # accidentally discard.
  original_snapshot <- if (is_edit) {
    list(
      data_root = cfg$data_dir,
      samples   = cfg$samples,
      flat      = cfg[c("font", "fontsize", "page_w", "chr_fmt", "show_track_axes", "hic_norm", "hic_style",
                         "zrange_cap", "atac_height", "hic_height", "hic_res_mode",
                         "hic_max_bins_per_view", "diff_colour_lo", "diff_colour_hi",
                         "junction_opacity_cap", "min_junction_support", "reference_colour", "comparison_colour",
                         "gene_flank_pct", "min_display_tpm", "uniprot_reviewed_only", "arc_clip_partial",
                         "gtf_path", "orgdb", "genome_label")]
    )
  } else NULL

  output$wiz_body <- renderUI({
    mode <- wiz_mode()
    if (is.null(mode)) {
      return(wellPanel(
        h4("Get started"),
        p("Import an existing samples.txt to prefill everything below, or set up a new project by hand."),
        fluidRow(
          column(6, actionButton("wiz_mode_import", "Import samples.txt", class = "btn-primary btn-lg", width = "100%")),
          column(6, actionButton("wiz_mode_manual", "Set up manually", class = "btn-default btn-lg", width = "100%"))
        )
      ))
    }
    tagList(
      if (!identical(mode, "edit"))
        actionLink("wiz_mode_back", "\u2190 Back", style = "display:inline-block; margin-bottom:10px;"),
      if (identical(mode, "import")) wellPanel(
        h4("Import samples.txt"),
        fluidRow(
          column(8, textInput("wiz_import_path", "Full path to samples.txt (on this server):", value = "", width = "100%")),
          column(4, style = "margin-top: 24px;", actionButton("wiz_import_load", "Import", class = "btn-default"))
        ),
        uiOutput("wiz_import_status")
      ),
      wiz_shared_form(cfg, is_edit)
    )
  })

  # Quick-fill the genome identity fields from a bundled GENOME_PRESETS
  # entry -- pure convenience, nothing is stored under the preset's id; Save
  # always writes the plain field values. GENOME_PRESETS carries no
  # gtf_path (2026-09-01, portability -- see utils/genomes.R's header), so
  # the GTF path field is cleared, not filled, on every preset pick: always
  # the user's own file, regardless of which preset (if any) was chosen.
  observeEvent(input$wiz_genome_preset, {
    preset <- GENOME_PRESETS[[input$wiz_genome_preset]]
    if (is.null(preset)) return()
    updateTextInput(session, "wiz_gtf_path", value = "")
    updateTextInput(session, "wiz_orgdb", value = preset$orgdb)
    updateTextInput(session, "wiz_genome_label", value = preset$genome_label)
  }, ignoreInit = TRUE)

  data_root <- reactive({
    dr <- trimws(input$wiz_data_root %||% "")
    if (nchar(dr) == 0) NULL else dr
  })

  output$wiz_data_root_status <- renderUI({
    dr <- data_root()
    if (is.null(dr)) return(NULL)
    if (dir.exists(dr)) tags$span(class = "wiz-ok", "\u2713 directory exists")
    else tags$span(class = "wiz-error", "\u2717 directory does not exist")
  })

  # ── Add-sample modal ──────────────────────────────────────────────────────

  # `error` NULL means a fresh modal open (fields blank/default); non-NULL means a re-display after a failed wiz_confirm_sample validation (fields seeded from current inputs instead, since Shiny never clears an unmounted input, so `error`'s presence is the only way to tell the two apart).
  add_sample_modal <- function(error = NULL) {
    dr    <- data_root()
    fresh <- is.null(error)
    seed  <- function(input_id, default) if (fresh) default else (isolate(input[[input_id]]) %||% default)
    field_errors <- if (fresh) list() else modal_field_errors()
    # Editing an existing sample (modal_edit_index() set) prefills name/
    # colour/quantdir from that sample. Adding a brand-new sample prefills
    # only the quantdir NAME (not path) -- feature 32, 2026-09-02 revision
    # (user's call): scans every already-committed sample for a quantdir
    # entry, most recent match wins, rather than only the very last sample --
    # so a naming convention set early on isn't lost just because the most
    # recent sample happens to have no quantdir at all. Nothing
    # sample-specific (condition name, colour, actual paths) is prefilled
    # into a new sample.
    edit_idx       <- modal_edit_index()
    editing_sample <- if (!is.null(edit_idx)) wiz_samples()[[edit_idx]] else NULL
    carry_quantdir <- if (is.null(editing_sample)) {
      quantdirs <- Filter(Negate(is.null), lapply(wiz_samples(), function(s) Find(function(d) identical(d$type, "quantdir"), s$data)))
      if (length(quantdirs) > 0) quantdirs[[length(quantdirs)]] else NULL
    } else NULL
    quantdir_inst  <- if (!is.null(editing_sample)) Find(function(d) identical(d$type, "quantdir"), editing_sample$data) else carry_quantdir
    default_name    <- editing_sample$condition %||% ""
    default_colour  <- editing_sample$colour %||% WIZARD_DEFAULT_COLOUR
    default_qname   <- quantdir_inst$name %||% "Quant"
    default_qpath   <- if (!is.null(editing_sample)) (quantdir_inst$path %||% "") else ""
    # Explicit "" choice: with only real file paths in `choices`,
    # selected="" matches no <option>, and
    # shiny's underlying native <select> then falls back to auto-selecting
    # the first real option before selectize.js upgrades it -- so a field the
    # user never touched silently pointed at whatever file happened to be
    # first. A real blank option gives selected="" something to bind to.
    quantdir_choices <- c("", wizard_scoped_choices(dr, "quantdir"))
    modalDialog(
      title = if (is.null(edit_idx)) "Add sample" else "Edit sample", size = "l", easyClose = FALSE,
      if (!is.null(error)) div(class = "wiz-error", error),
      fluidRow(
        column(6, textInput("wiz_sample_name", "Condition (unique):", value = seed("wiz_sample_name", default_name))),
        column(6, tags$label("Colour:"), br(),
               # shiny.js has no input binding for type="color" (confirmed against the
               # installed shiny.js; app.R's own s_diff_colour_lo/hi settings inputs hit
               # this same gap and work around it identically) -- without the oninput
               # handler below, input$wiz_sample_colour never updates from its initial value.
               tags$input(type = "color", id = "wiz_sample_colour", value = seed("wiz_sample_colour", default_colour),
                          style = "width:60px; height:34px; padding:2px; border:1px solid #ccc; border-radius:4px; cursor:pointer; display:block;",
                          oninput = "Shiny.setInputValue('wiz_sample_colour', this.value, {priority: 'event'});"))
      ),
      hr(),
      fluidRow(
        column(3, actionButton("wiz_add_bw",         "+ BigWig instance",     class = "btn-sm btn-default")),
        column(3, actionButton("wiz_add_hic",        "+ Hi-C instance",       class = "btn-sm btn-default")),
        column(3, actionButton("wiz_add_bed",        "+ BED instance",        class = "btn-sm btn-default")),
        column(3, actionButton("wiz_add_paired_end", "+ Paired-end instance", class = "btn-sm btn-default"))
      ),
      uiOutput("wiz_instance_rows"),
      hr(),
      p(class = "text-muted", style = "font-size:11px;", "Optional, one per sample:"),
      fluidRow(class = if (is.null(field_errors$quantdir)) "wiz-row" else "wiz-row wiz-row-error",
        column(2, tags$b("QUANT")),
        column(4, textInput("wiz_quantdir_name", NULL, value = seed("wiz_quantdir_name", default_qname))),
        column(6, selectizeInput("wiz_quantdir_path", NULL, choices = quantdir_choices,
                                  selected = seed("wiz_quantdir_path", default_qpath),
                                  options = list(create = TRUE, placeholder = "path relative to data root"))),
        if (!is.null(field_errors$quantdir)) column(12, tags$div(class = "wiz-error", style = "font-size:11px; margin-top:4px;", field_errors$quantdir))
      ),
      footer = tagList(modalButton("Cancel"), actionButton("wiz_confirm_sample", if (is.null(edit_idx)) "Add sample" else "Save changes", class = "btn-primary"))
    )
  }

  # `samples` is a list to prefill rows from: keep_path=TRUE for one sample being edited, FALSE for the deduped union of every committed sample's type:name pairs (feature 32); quantdir isn't included here, see add_sample_modal(). See solutions/concepts/19-relocated-code-rationale.md, "utils/wizard.R:458".
  wiz_prefill_rows <- function(samples, keep_path) {
    rows <- list()
    seen <- character(0)
    for (s in samples) {
      for (d in s$data) {
        if (identical(d$type, "quantdir")) next
        key <- paste0(d$type, ":", d$name)  # ":" isn't valid in a WIZARD_DATA_TYPES key, so this can't collide across types
        if (key %in% seen) next
        seen <- c(seen, key)
        rid <- modal_row_id() + 1
        modal_row_id(rid)
        rows[[length(rows) + 1]] <- list(id = rid, type = d$type, name = d$name,
                                          path = if (keep_path) d$path else "")
      }
    }
    rows
  }

  observeEvent(input$wiz_add_sample, {
    modal_edit_index(NULL)
    modal_instances(wiz_prefill_rows(wiz_samples(), keep_path = FALSE))
    modal_row_errors(list())
    modal_field_errors(list())
    # modal_row_id is NOT reset here: its ids are
    # input-element ids (wiz_inst_name_<id>/wiz_inst_path_<id>) that live
    # for the whole Shiny session even after their row is removed from the
    # DOM -- Shiny never clears input$wiz_inst_name_1 just because the
    # element disappeared. Resetting to 0 on every "+ Add sample" click meant
    # a fresh row on sample 2 reused id 1 and silently inherited sample 1's
    # already-committed values. Keeping the counter monotonic for the whole
    # session means an id is retired for good once used, so this can't recur.
    showModal(add_sample_modal())
  })

  observeEvent(input$wiz_edit_sample, {
    idx <- input$wiz_edit_sample
    modal_edit_index(idx)
    modal_instances(wiz_prefill_rows(list(wiz_samples()[[idx]]), keep_path = TRUE))
    modal_row_errors(list())
    modal_field_errors(list())
    showModal(add_sample_modal())
  }, ignoreInit = TRUE)

  add_modal_row <- function(type) {
    rid <- modal_row_id() + 1
    modal_row_id(rid)
    modal_instances(c(modal_instances(), list(list(id = rid, type = type))))
  }
  observeEvent(input$wiz_add_bw,         add_modal_row("bw"))
  observeEvent(input$wiz_add_hic,        add_modal_row("hic"))
  observeEvent(input$wiz_add_bed,        add_modal_row("bed"))
  observeEvent(input$wiz_add_paired_end, add_modal_row("paired-end"))

  observeEvent(input$wiz_remove_row, {
    modal_instances(Filter(function(r) r$id != input$wiz_remove_row, modal_instances()))
  }, ignoreInit = TRUE)

  output$wiz_instance_rows <- renderUI({
    rows <- modal_instances()
    dr <- data_root()
    errs <- modal_row_errors()
    if (length(rows) == 0) return(p(class = "text-muted", "No BigWig / Hi-C / BED / Paired-end instances added yet."))
    tagList(lapply(rows, function(r) {
      name_id <- paste0("wiz_inst_name_", r$id)
      path_id <- paste0("wiz_inst_path_", r$id)
      row_err <- errs[[as.character(r$id)]]
      fluidRow(class = if (is.null(row_err)) "wiz-row" else "wiz-row wiz-row-error",
        column(2, tags$b(toupper(r$type))),
        # r$name/r$path (return-pass fix): a row added via wiz_prefill_rows()
        # (edit-sample prefill, or a new sample carrying over the previous
        # one's instance layout) carries its prefill value here, since its
        # input id is freshly minted and isolate(input[[...]]) is NULL until
        # the user actually touches it. A row added via "+ <type> instance"
        # has no r$name/r$path (NULL), so this falls through to "" exactly
        # as before.
        column(4, textInput(name_id, NULL, value = isolate(input[[name_id]]) %||% r$name %||% "", placeholder = "instance name")),
        # "" prepended (return-pass fix, see add_sample_modal()): otherwise
        # selected="" matches no real <option> and the browser silently
        # falls back to the first matching file.
        column(5, selectizeInput(path_id, NULL, choices = c("", wizard_scoped_choices(dr, r$type)),
                                  selected = isolate(input[[path_id]]) %||% r$path %||% "",
                                  options = list(create = TRUE, placeholder = "path relative to data root"))),
        column(1, actionButton(paste0("wiz_rm_btn_", r$id), "\u2715", class = "btn-xs",
                                onclick = sprintf("Shiny.setInputValue('wiz_remove_row', %d, {priority:'event'})", r$id))),
        if (!is.null(row_err)) column(12, tags$div(class = "wiz-error", style = "font-size:11px; margin-top:4px;", row_err))
      )
    }))
  })

  observeEvent(input$wiz_confirm_sample, {
    edit_idx <- modal_edit_index()
    name <- trimws(input$wiz_sample_name %||% "")
    existing_names <- vapply(wiz_samples(), function(s) s$condition, character(1))
    # Editing a sample must not collide with its OWN pre-existing name --
    # excluding it here (rather than just allowing name == old name as a
    # special case) also lets the uniqueness check below stay a plain %in%.
    if (!is.null(edit_idx)) existing_names <- existing_names[-edit_idx]
    errs <- character(0)
    row_errs   <- list()  # row id (character) -> message, for wiz_instance_rows highlighting
    field_errs <- list()  # "quantdir" -> message, for add_sample_modal() highlighting

    if (nchar(name) == 0) errs <- c(errs, "Condition is required.")
    else if (name %in% existing_names) errs <- c(errs, sprintf("A condition named '%s' already exists.", name))

    data <- list()
    for (r in modal_instances()) {
      nm  <- trimws(input[[paste0("wiz_inst_name_", r$id)]] %||% "")
      pth <- trimws(input[[paste0("wiz_inst_path_", r$id)]] %||% "")
      if (nchar(pth) == 0) next  # blank instance = not declared, silently dropped (per phase-C.1.md step 2)
      if (nchar(nm) == 0) {
        msg <- sprintf("A %s instance has a path but no name.", r$type)
        errs <- c(errs, msg); row_errs[[as.character(r$id)]] <- msg
        next
      }
      if (wizard_name_has_bad_dot(nm)) {
        msg <- sprintf("Instance name '%s' contains a bare '\u00b7' (surround it with spaces, or remove it). This breaks track-label formatting.", nm)
        errs <- c(errs, msg); row_errs[[as.character(r$id)]] <- msg
        next
      }
      if (any(vapply(data, function(d) d$type == r$type && d$name == nm, logical(1)))) {
        msg <- sprintf("Duplicate %s instance name '%s'.", r$type, nm)
        errs <- c(errs, msg); row_errs[[as.character(r$id)]] <- msg
        next
      }
      # is_real_file(), not file.exists(): file.exists() is TRUE for
      # directories too, which would otherwise be accepted here only to fail
      # (or, for paired-end, crash) at first render.
      if (!is.null(data_root()) && !is_real_file(file.path(data_root(), pth))) {
        msg <- sprintf("Path for %s:%s does not resolve to a file under the data root: %s", r$type, nm, pth)
        errs <- c(errs, msg); row_errs[[as.character(r$id)]] <- msg
      } else if (identical(r$type, "paired-end") && !is.null(data_root()) &&
                 !is_paired_end_structural(file.path(data_root(), pth))) {
        # Cheap fail-fast check (utils/samples_txt_io.R): must be BED12 (12
        # cols) or BEDPE (6-10 cols) -- read_junction_file()'s real dispatch,
        # solutions/concepts/18-bed12-bedpe-arc-format.md.
        msg <- sprintf("Paired-end file for %s is not BED12 (12 columns) or BEDPE (6-10 columns): %s", nm, pth)
        errs <- c(errs, msg); row_errs[[as.character(r$id)]] <- msg
      }
      data[[length(data) + 1]] <- list(type = r$type, name = nm, path = pth)
    }

    quantdir_path <- trimws(input$wiz_quantdir_path %||% "")
    if (nchar(quantdir_path) > 0) {
      quantdir_name <- trimws(input$wiz_quantdir_name %||% "")
      if (nchar(quantdir_name) == 0) quantdir_name <- "Quant"
      if (!is.null(data_root()) && !is_real_file(file.path(data_root(), quantdir_path, "quant.sf"))) {
        msg <- sprintf("No quant.sf found under: %s", quantdir_path)
        errs <- c(errs, msg); field_errs$quantdir <- msg
      }
      data[[length(data) + 1]] <- list(type = "quantdir", name = quantdir_name, path = quantdir_path)
    }

    if (length(errs) > 0) {
      modal_row_errors(row_errs)
      modal_field_errors(field_errs)
      showModal(add_sample_modal(error = paste(errs, collapse = " | ")))
      return()
    }

    colour <- input$wiz_sample_colour %||% WIZARD_DEFAULT_COLOUR
    new_sample <- list(condition = name, colour = colour, data = data)
    if (is.null(edit_idx)) {
      wiz_samples(c(wiz_samples(), list(new_sample)))
    } else {
      samples <- wiz_samples()
      samples[[edit_idx]] <- new_sample
      wiz_samples(samples)
    }
    modal_edit_index(NULL)
    removeModal()
  })

  # ── Committed-sample list + edit/remove ───────────────────────────────────
  output$wiz_samples_ui <- renderUI({
    samples <- wiz_samples()
    if (length(samples) == 0) return(p(class = "text-muted", "No samples added yet."))
    tagList(lapply(seq_along(samples), function(i) {
      s <- samples[[i]]
      inst_labels <- if (length(s$data) > 0) vapply(s$data, function(d) paste0(d$type, ":", d$name), character(1)) else character(0)
      fluidRow(class = "wiz-row",
        column(1, tags$span(style = sprintf("display:inline-block;width:14px;height:14px;background:%s;border-radius:2px;", s$colour %||% "#ccc"))),
        column(3, tags$b(s$condition)),
        column(5, paste(inst_labels, collapse = ", ")),
        column(3,
               actionButton(paste0("wiz_edit_sample_btn_", i), "Edit", class = "btn-xs btn-default",
                             onclick = sprintf("Shiny.setInputValue('wiz_edit_sample', %d, {priority:'event'})", i)),
               actionButton(paste0("wiz_rm_sample_btn_", i), "Remove", class = "btn-xs btn-default",
                             onclick = sprintf("Shiny.setInputValue('wiz_remove_sample', %d, {priority:'event'})", i)))
      )
    }))
  })

  observeEvent(input$wiz_remove_sample, {
    wiz_samples(wiz_samples()[-input$wiz_remove_sample])
  }, ignoreInit = TRUE)

  # ── Import from samples.txt ───────────────────────────────────────────────
  output$wiz_import_status <- renderUI(NULL)

  observeEvent(input$wiz_import_load, {
    path <- trimws(input$wiz_import_path %||% "")
    if (nchar(path) == 0) {
      output$wiz_import_status <- renderUI(div(class = "wiz-error", "Enter a path to samples.txt first."))
      return()
    }
    if (!is_real_file(path)) {
      output$wiz_import_status <- renderUI(div(class = "wiz-error", sprintf("No such file: %s", path)))
      return()
    }
    res <- tryCatch(
      # defaults = DEFAULT_SETTINGS (Phase C.3 decision 6): lets any
      # `# <key>: value` line in the file be parsed back into
      # `res$overrides` and applied to the Other-settings fields below --
      # since the standalone genome-config task, that includes the 3 genome
      # identity fields (gtf_path/orgdb/genome_label), which are plain
      # DEFAULT_SETTINGS keys with no separate `# genome:` line or id
      # validation needed. Cache paths (txdb_path/tx_name_lookup_path/
      # tx_meta_lookup_path) are no longer settings at all (2026-08-25) --
      # a stray `# txdb_path: ...` line from an old export is simply ignored.
      import_samples_txt(path, defaults = DEFAULT_SETTINGS),
      error = function(e) list(errors = conditionMessage(e))
    )
    if (length(res$errors) > 0) {
      output$wiz_import_status <- renderUI(
        div(class = "wiz-error", tags$b("Import failed:"), tags$ul(lapply(res$errors, tags$li)))
      )
      return()
    }
    if (!is.null(res$data_root)) updateTextInput(session, "wiz_data_root", value = res$data_root)
    wiz_samples(res$samples)

    # Apply any `# <key>: value` overrides onto the Other-settings inputs
    # (Phase C.3 decision 6). The four colour fields are plain
    # <input type="color"> elements with no shiny.js input binding (see
    # add_sample_modal()'s comment on wiz_sample_colour) -- there's no
    # update*Input() for those, so a colour override from an imported file is
    # a known, minor limitation: it's parsed but not reflected in the swatch.
    ov <- res$overrides
    if (!is.null(ov$gtf_path))            updateTextInput(session, "wiz_gtf_path", value = ov$gtf_path)
    if (!is.null(ov$orgdb))               updateTextInput(session, "wiz_orgdb", value = ov$orgdb)
    if (!is.null(ov$genome_label))        updateTextInput(session, "wiz_genome_label", value = ov$genome_label)
    if (!is.null(ov$font))         updateSelectizeInput(session, "wiz_font", selected = ov$font)
    if (!is.null(ov$fontsize))     updateNumericInput(session, "wiz_fontsize", value = ov$fontsize)
    if (!is.null(ov$page_w))       updateNumericInput(session, "wiz_page_w", value = ov$page_w)
    if (!is.null(ov$chr_fmt))      updateSelectizeInput(session, "wiz_chr_fmt", selected = ov$chr_fmt)
    if (!is.null(ov$show_track_axes)) updateCheckboxInput(session, "wiz_show_track_axes", value = ov$show_track_axes)
    if (!is.null(ov$hic_norm))     updateTextInput(session, "wiz_hic_norm", value = ov$hic_norm)
    if (!is.null(ov$hic_style))    updateSelectizeInput(session, "wiz_hic_style", selected = ov$hic_style)
    if (!is.null(ov$zrange_cap))   updateNumericInput(session, "wiz_zrange_cap", value = ov$zrange_cap)
    if (!is.null(ov$atac_height))  updateNumericInput(session, "wiz_atac_height", value = ov$atac_height)
    if (!is.null(ov$hic_height))   updateNumericInput(session, "wiz_hic_height", value = ov$hic_height)
    if (!is.null(ov$hic_res_mode)) updateSelectizeInput(session, "wiz_hic_res_mode", selected = ov$hic_res_mode)
    if (!is.null(ov$hic_max_bins_per_view)) updateNumericInput(session, "wiz_hic_max_bins_per_view", value = ov$hic_max_bins_per_view)
    if (!is.null(ov$junction_opacity_cap))  updateNumericInput(session, "wiz_junction_opacity_cap", value = ov$junction_opacity_cap)
    if (!is.null(ov$min_junction_support))  updateNumericInput(session, "wiz_min_junction_support", value = ov$min_junction_support)
    if (!is.null(ov$gene_flank_pct))        updateNumericInput(session, "wiz_gene_flank_pct", value = ov$gene_flank_pct)
    if (!is.null(ov$min_display_tpm))       updateNumericInput(session, "wiz_min_display_tpm", value = ov$min_display_tpm)
    if (!is.null(ov$uniprot_reviewed_only)) updateCheckboxInput(session, "wiz_uniprot_reviewed_only", value = ov$uniprot_reviewed_only)
    if (!is.null(ov$arc_clip_partial))      updateCheckboxInput(session, "wiz_arc_clip_partial", value = ov$arc_clip_partial)

    status_lines <- list(div(class = "wiz-ok", sprintf("Imported %d sample(s).", length(res$samples))))
    if (length(res$warnings) > 0)
      status_lines <- c(status_lines, list(div(class = "wiz-error", tags$b("Warnings:"), tags$ul(lapply(res$warnings, tags$li)))))
    output$wiz_import_status <- renderUI(tagList(status_lines))
  })

  # ── Export samples.txt ────────────────────────────────────────────────────

  # A separate action from Save, never triggered automatically (Phase C.3 decision 2).
  wiz_collect_flat_settings <- function() {
    list(
      gtf_path            = wiz_clean_optional_text(input$wiz_gtf_path),
      orgdb               = wiz_clean_optional_text(input$wiz_orgdb),
      genome_label        = wiz_clean_optional_text(input$wiz_genome_label),
      font        = trimws(input$wiz_font %||% "DejaVu Sans"),
      fontsize    = input$wiz_fontsize %||% 8,
      page_w      = input$wiz_page_w %||% 9,
      chr_fmt     = input$wiz_chr_fmt %||% "dot",
      show_track_axes = input$wiz_show_track_axes %||% TRUE,
      hic_norm    = trimws(input$wiz_hic_norm %||% "KR"),
      hic_style   = input$wiz_hic_style %||% "rectangle",
      zrange_cap  = wiz_clean_optional_numeric(input$wiz_zrange_cap),
      atac_height = input$wiz_atac_height %||% 0.85,
      hic_height  = wiz_clean_optional_numeric(input$wiz_hic_height),
      hic_res_mode = input$wiz_hic_res_mode %||% "auto",
      hic_max_bins_per_view = input$wiz_hic_max_bins_per_view %||% 500,
      # The 4 colour fields below are plain <input type="color"> elements
      # with no shiny.js input binding (see add_sample_modal()'s comment) --
      # input$wiz_* only updates once the user actually touches that specific
      # swatch, unlike every numericInput/selectInput above (which DO get
      # their `value=`/`selected=` reported back on initial render). Without
      # the cfg$... middle fallback, an untouched swatch would silently
      # revert to the hardcoded default instead of the field's real current
      # value -- both on Save (discarding a customized colour) and in
      # is_dirty()'s comparison (spuriously flagging a change that never
      # happened).
      diff_colour_lo    = input$wiz_diff_colour_lo    %||% cfg$diff_colour_lo    %||% "#0033CC",
      diff_colour_hi    = input$wiz_diff_colour_hi    %||% cfg$diff_colour_hi    %||% "#CC0000",
      junction_opacity_cap = input$wiz_junction_opacity_cap %||% 2,
      min_junction_support = input$wiz_min_junction_support %||% 0,
      reference_colour  = input$wiz_reference_colour  %||% cfg$reference_colour  %||% "#4D4D4D",
      comparison_colour = input$wiz_comparison_colour %||% cfg$comparison_colour %||% "#56B4E9",
      gene_flank_pct    = input$wiz_gene_flank_pct %||% 15,
      min_display_tpm   = input$wiz_min_display_tpm %||% 1,
      uniprot_reviewed_only = input$wiz_uniprot_reviewed_only %||% FALSE,
      arc_clip_partial  = input$wiz_arc_clip_partial %||% FALSE
    )
  }

  observeEvent(input$wiz_export_samples, {
    dr <- data_root()
    showModal(modalDialog(
      title = "Export settings to txt", easyClose = TRUE,
      textInput("wiz_export_path", "Write to (full path):",
                value = if (is.null(dr)) "" else file.path(dr, "samples.txt"), width = "100%"),
      uiOutput("wiz_export_status"),
      footer = tagList(modalButton("Cancel"), actionButton("wiz_confirm_export", "Export", class = "btn-primary"))
    ))
  })

  output$wiz_export_status <- renderUI(NULL)

  observeEvent(input$wiz_confirm_export, {
    path <- trimws(input$wiz_export_path %||% "")
    if (nchar(path) == 0) {
      output$wiz_export_status <- renderUI(div(class = "wiz-error", "Enter a path first."))
      return()
    }
    settings <- c(
      list(data_root = data_root(), samples = wiz_samples()),
      wiz_collect_flat_settings()
    )
    ok <- tryCatch({
      writeLines(export_samples_txt(settings, defaults = DEFAULT_SETTINGS), path)
      TRUE
    }, error = function(e) {
      output$wiz_export_status <- renderUI(div(class = "wiz-error", sprintf("Failed to write: %s", conditionMessage(e))))
      FALSE
    })
    if (ok) output$wiz_export_status <- renderUI(div(class = "wiz-ok", sprintf("Wrote %s.", path)))
  })

  # ── Dirty guard (edit mode only, Phase C.3 steps 12/13) ───────────────────
  is_dirty <- reactive({
    if (!is_edit) return(FALSE)
    if (!identical(data_root(), original_snapshot$data_root)) return(TRUE)
    if (!identical(wiz_samples(), original_snapshot$samples)) return(TRUE)
    flat <- wiz_collect_flat_settings()
    any(vapply(names(flat), function(k) settings_value_differs(flat[[k]], original_snapshot$flat[[k]]), logical(1)))
  })

  wiz_do_close <- function() {
    # Clear so a later re-click of "Edit project settings" from the main app
    # re-enters the wizard on its own terms rather than being permanently
    # stuck routing back here.
    Sys.unsetenv("PLOTTRELLIS_EDIT_SETTINGS")
    showModal(modalDialog(
      title = "Closing", easyClose = FALSE, footer = NULL,
      p("Returning to the app without saving, a new browser tab will open in a few seconds.")
    ))
    message("PlotTrellis: settings unchanged, returning to the app.")
    later::later(function() stopApp(returnValue = "restart"), 0.2)
  }

  observeEvent(input$wiz_close, {
    if (is_dirty()) {
      showModal(modalDialog(
        title = "Discard changes?", easyClose = FALSE,
        p("You have unsaved changes. Closing now will discard them and return to the app with its previous settings."),
        footer = tagList(modalButton("Keep editing"), actionButton("wiz_confirm_close", "Discard and close", class = "btn-danger"))
      ))
    } else {
      wiz_do_close()
    }
  })

  observeEvent(input$wiz_confirm_close, {
    removeModal()
    wiz_do_close()
  }, ignoreInit = TRUE)

  # ── Save ───────────────────────────────────────────────────────────────────
  output$wiz_completion_error <- renderUI(NULL)

  observeEvent(input$wiz_save, {
    dr <- data_root()
    errs <- character(0)
    if (is.null(dr)) errs <- c(errs, "Data root is required.")
    else if (!dir.exists(dr)) errs <- c(errs, "Data root directory does not exist.")
    if (length(wiz_samples()) == 0) errs <- c(errs, "Add at least one sample before saving.")
    if (is.null(wiz_clean_optional_text(input$wiz_gtf_path))) errs <- c(errs, "GTF path is required.")
    if (is.null(wiz_clean_optional_text(input$wiz_orgdb))) errs <- c(errs, "OrgDb package is required.")

    if (length(errs) > 0) {
      output$wiz_completion_error <- renderUI(div(class = "wiz-error", paste(errs, collapse = " | ")))
      return()
    }

    settings <- c(
      list(
        data_dir        = dr,
        samples         = wiz_samples()
      ),
      wiz_collect_flat_settings()
    )
    # null = "null": jsonlite's toJSON() default serializes a NULL list element
    # (e.g. zrange_cap/hic_height left blank for "auto") as {} rather than the
    # JSON literal null. On reload, fromJSON() then reads that back as an empty
    # list, not NULL, and clean_optional_numeric()'s is.finite() crashes on it.
    writeLines(jsonlite::toJSON(settings, auto_unbox = TRUE, pretty = TRUE, null = "null"), SETTINGS_JSON)

    # See wiz_do_close()'s comment -- same reasoning applies to Save.
    Sys.unsetenv("PLOTTRELLIS_EDIT_SETTINGS")

    output$wiz_completion_error <- renderUI(NULL)
    showModal(modalDialog(
      title = "Settings saved", easyClose = FALSE, footer = NULL,
      p("project_settings.json has been written."),
      p("Restarting with your configured project, a new browser tab will open in a few seconds.")
    ))
    message("PlotTrellis: settings saved, launching the app.")
    later::later(function() stopApp(returnValue = "restart"), 0.2)
  })
}
