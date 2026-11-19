# PlotTrellis composer: one page with a Compose section and an Output section.
# Figures render with render_all() in a background R process
# (R/viz-render-worker.R), so a long render can be cancelled; the R script the page
# shows and downloads re-runs the same functions on the same entries and
# settings (assemble_download_script()).
#
# The exported entry points (run_plottrellis() and friends) are in
# R/api-entry-points.R; this file holds the project-settings defaults/reader and the
# composer app builder (.pt_composer_app()). .pt_composer_app() builds the
# startup context (R/app-context.R) and the page (R/app-ui.R), and its
# server calls one top-level function per block (the R/app-*.R block files); the
# closure keeps only the names tests read through it.

# The display settings and their defaults: how any track is drawn or filtered.
# These are the settings of the composer's Figure settings panel, and the ones
# render_all() fills in when a caller leaves them out (.pt_fill_settings()).
# NULL means "decided when the figure is drawn".
.pt_display_defaults <- function() list(
  font        = "DejaVu Sans",
  fontsize    = 8L,     # points
  page_w      = 9,      # inches
  # the chromosome name on the genome axis: "dot" (chr.19), "space" (chr 19),
  # "plain" (chr19)
  chr_fmt     = "dot",
  # The value axis on coverage, coverage-difference and junction tracks, which
  # show no other value scale (a Hi-C track has its colour legend, Transcript
  # expression its own TPM axis).
  show_track_axes = TRUE,
  # the figure's title line, when the figure has a title
  show_title      = TRUE,
  # the genomic scale bar on a figure's first track, for the types that draw one
  show_scale_bar  = TRUE,
  hic_norm    = "KR",
  hic_style   = "rectangle",
  hic_colour_scale = "linear", # "linear" | "log"
  # top of the Hi-C colour scale (Colour max): "max", "p95" or "p90" of the
  # contacts drawn, or a number; NULL = Auto, "p95" on a linear colour scale
  # and "max" on a log one
  zrange_cap  = NULL,
  atac_height     = 0.85,   # inches, a coverage track
  # inches; NULL = half the track width (rectangle) or all of it (square)
  hic_height      = NULL,
  hic_res_mode    = "auto", # "auto" | "finest" | "fixed"
  # fewest bins across the view that the automatic Hi-C resolution accepts
  hic_max_bins_per_view = 500,
  # bin size in bp for hic_res_mode "fixed"; NULL = not set
  hic_resolution  = NULL,
  gene_flank_pct  = 15,     # % of the gene's length shown on each side of it
  # Least pooled support an arc (a BED12 splice junction or a BEDPE loop) needs
  # to be drawn. 0 shows every arc: a BEDPE score's scale depends on the caller
  # that wrote it, so no other value is safe for every file.
  min_junction_support = 0,
  # Mean TPM a transcript must reach in a sample to be drawn and named on its
  # own, instead of being folded into "Other".
  min_display_tpm = 1,
  # the Gene info panel's UniProt link: Swiss-Prot entries only if TRUE
  uniprot_reviewed_only = FALSE,
  # Junction tracks: TRUE cuts an arc with one end outside the region at the
  # region's edge, FALSE draws it running past the edge. Long-range BEDPE pairs
  # (Hi-C loops) are the case that needs TRUE; splice junctions rarely reach
  # that far.
  arc_clip_partial = FALSE,
  # inches, a junction track's arcs per sample, under its heading (a mirrored
  # comparison: 2x)
  junction_height  = 0.85,

  # ---- Two-sample comparison tracks only ----
  diff_colour_lo  = "#0033CC",  # vivid blue: lower in Comparison
  diff_colour_hi  = "#CC0000",  # vivid red: higher in Comparison
  # Mirrored junction difference: the log2 fold change drawn fully opaque. A
  # ceiling, so one sparse arc's large ratio cannot stretch the scale and wash
  # out every other arc.
  junction_opacity_cap = 2,
  # Coverage difference: tint of the shared (overlap) layer, 0 = white to 1 =
  # the full reference colour.
  diff_shared_alpha = 0.15,
  # The two role colours of every comparison track: Reference and Comparison are
  # always drawn in these, whatever the two samples' own colours, so which one
  # is the baseline reads at a glance.
  reference_colour  = "#4D4D4D", # neutral grey
  comparison_colour = "#56B4E9"  # colour-blind-safe blue
)

# Every project setting with its default, as one flat list (reading the settings
# file, the wizard and the composer all rely on it being flat): the core import
# settings, which the wizard sets (.PT_CORE_SETTINGS), then the display
# settings. A function, not a top-level constant: app-main.R collates before
# genome-main.R, so a top-level list(orgdb = GENOME_PRESETS$...) would fail at
# install.
.pt_default_settings <- function() c(
  list(
    # Blank until the wizard is given real paths, which it requires before Save.
    gtf_path     = "",
    orgdb        = GENOME_PRESETS$mm39$orgdb,
    genome_label = GENOME_PRESETS$mm39$genome_label,
    data_dir     = "",
    # The samples live in the settings file itself (samples.txt is an
    # import/export format, never read at startup). Each is list(condition =,
    # colour =, data = list(list(type =, name =, path =))), the shape
    # import_samples_txt() and export_samples_txt() use; `condition` alone
    # identifies a sample, so a replicate is a separate entry with its own
    # condition.
    samples      = list()
  ),
  .pt_display_defaults()
)

.pt_read_settings <- function(SETTINGS_JSON, DEFAULT_SETTINGS) {
  if (!file.exists(SETTINGS_JSON)) return(DEFAULT_SETTINGS)
  # simplifyDataFrame = FALSE: jsonlite's default simplification turns
  # `samples[]` (an array of same-shaped objects) into a data.frame instead
  # of a list of lists, which crashes load_samples_from_settings()'s
  # `s$condition`/`s$data` access (`$` on an atomic column vector). Scalars
  # are unaffected -- only array-of-objects simplification is disabled.
  parsed <- tryCatch(jsonlite::fromJSON(SETTINGS_JSON,
                                        simplifyDataFrame = FALSE),
                     error = function(e) {
    # Point at the project_samples.txt the wizard's Save keeps beside it.
    txt <- file.path(dirname(SETTINGS_JSON), "project_samples.txt")
    stop("Can't read ", SETTINGS_JSON, " (", conditionMessage(e), "). ",
         if (file.exists(txt))
           paste0("Rebuild it from the project_samples.txt saved beside it: parse_samples_txt(",
                  encodeString(txt, quote = "\""), ", project_dir = ",
                  encodeString(dirname(SETTINGS_JSON), quote = "\""),
                  ", overwrite = TRUE)")
         else "Delete it and run the setup wizard again, or rebuild it from a samples.txt with parse_samples_txt(..., overwrite = TRUE).",
         call. = FALSE)
  })
  cfg <- modifyList(DEFAULT_SETTINGS, parsed)
  # modifyList()'s default (keep.null = FALSE) recurses into any key that's a
  # list on BOTH sides instead of overwriting it outright -- `samples` is
  # exactly such a key (DEFAULT_SETTINGS$samples = list()), and a JSON array
  # always parses to an UNNAMED list, so that recursive merge finds no names
  # to copy and silently discards every sample. Restore it directly instead
  # of trusting modifyList's recursive merge with it.
  cfg$samples <- if (is.null(parsed$samples)) list() else parsed$samples
  .pt_clean_display_settings(cfg, SETTINGS_JSON)
}

# A hand-edited display setting of the wrong type or length (page_w = "wide",
# fontsize = []), or a page_w of 0 or below, would stop every figure later with
# an R error: each one falls back to its default instead, with one warning
# naming them all. A box the app saved empty (null, NA, or toJSON()'s "NA")
# takes the default without a warning. Unset optional numbers (NULL default) may
# also be {}, which the render settings treat as automatic.
.pt_clean_display_settings <- function(cfg, SETTINGS_JSON) {
  defaults <- .pt_display_defaults()
  bad <- character(0)
  for (k in names(defaults)) {
    v <- cfg[[k]]; d <- defaults[[k]]
    if (is.null(v) ||
        (length(v) == 1 && !is.list(v) && (is.na(v) || identical(v, "NA")))) {
      cfg[k] <- list(d)
      next
    }
    if (is.null(d)) {
      if (is.list(v)) next
      # Colour max: also "max", "p95" or "p90" (.pt_hic_zrange_cap()).
      if (k == "zrange_cap") {
        cap <- .pt_hic_zrange_cap(v, bad = "bad")
        if (identical(cap, "bad")) { bad <- c(bad, k); cap <- NULL }
        cfg[k] <- list(cap)
        next
      }
      if (!(is.numeric(v) && length(v) == 1 && is.finite(v))) {
        bad <- c(bad, k); cfg[k] <- list(NULL)
      }
      next
    }
    ok <- length(v) == 1 && !is.list(v) &&
      (if (is.numeric(d)) is.numeric(v) && is.finite(v)
       else if (is.logical(d)) is.logical(v) else is.character(v))
    # A page width of 0 or below stops the whole figure (a narrow one only its
    # tracks, each with its own message).
    if (ok && k == "page_w") ok <- v > 0
    if (!ok) { bad <- c(bad, k); cfg[[k]] <- d }
  }
  if (length(bad) > 0)
    # immediate: the app's runApp() would hold it until the app closes
    warning("Ignoring unusable value", if (length(bad) > 1) "s", " in ",
            SETTINGS_JSON, " for ", paste(bad, collapse = ", "),
            "; using the default", if (length(bad) > 1) "s", ".",
            call. = FALSE, immediate. = TRUE)
  cfg
}

# Shared setup for every entry point (R/api-entry-points.R): creates project_dir and
# reads its settings merged onto the defaults.
.pt_project <- function(project_dir) {
  # The project folder may not exist yet -- on a
  # genuinely first run (nothing has ever written here before) the wizard's
  # Save handler's writeLines(SETTINGS_JSON) fails with "No such file or
  # directory" otherwise.
  dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)
  # Absolute, so the cache/TxDb paths that end up in the render settings (and so
  # in a downloaded script, pt_save() after setwd(), ...) don't depend on the
  # working directory.
  project_dir <- normalizePath(project_dir, winslash = "/")
  SETTINGS_JSON <- file.path(project_dir, "project_settings.json")
  DEFAULT_SETTINGS <- .pt_default_settings()
  # A project's caches, genome files included, go to its own cache/ folder.
  list(settings_json = SETTINGS_JSON, defaults = DEFAULT_SETTINGS,
       cfg = .pt_read_settings(SETTINGS_JSON, DEFAULT_SETTINGS),
       cache_dir = file.path(project_dir, "cache"))
}

# A chromosome's length from the chrom sizes (CHR_LENGTHS), or `default` for
# one they lack (a GTF scaffold UCSC does not list): `[[` on the named vector
# stops there.
.pt_chr_length <- function(lengths, chrom, default) {
  if (length(chrom) == 1L && !is.na(chrom) && chrom %in% names(lengths))
    lengths[[chrom]] else default
}

# Builds (does not run) the composer app for one project. Its startup below
# (genome migration, TxDb load, ...) runs on each call.
.pt_composer_app <- function(SETTINGS_JSON, CFG, RENDER_POOL = NULL) {

# CFG/SETTINGS_JSON come from .pt_project(), see R/api-entry-points.R. RENDER_POOL:
# the launcher's render processes (.pt_render_pool(), R/viz-render-worker.R), or
# NULL: one per session.

# Startup: everything loaded once per launch (R/app-context.R).
ctx <- .pt_composer_context(SETTINGS_JSON, CFG)
# The main screen (R/app-screen-server.R) unless the classic one is asked for
# (.pt_screen()); the rest of this function is the classic composer.
if (identical(.pt_screen(), "new"))
  return(.pt_main_app(ctx, SETTINGS_JSON, RENDER_POOL))
# Tests read the app's startup state through the server's environment
# (environment(app$serverFuncSource())): these names stay as locals here.
TRACK_ITEMS <- ctx$TRACK_ITEMS
RENDER_CFG  <- ctx$RENDER_CFG
START_VIEW  <- ctx$START_VIEW
PROJECT_DIR <- ctx$PROJECT_DIR
CACHE_DIR   <- ctx$CACHE_DIR
track_split <- function(conditions)
  .pt_track_split(conditions, TRACK_ITEMS, START_VIEW)

# ── UI ────────────────────────────────────────────────────────────────────────
# The page (R/app-ui.R); its stylesheet and script are
# inst/www/composer.css and composer.js.
ui <- .pt_composer_ui(ctx)

# ── Server ────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {
  # repeat-click gate for the "pt-once" buttons (R/app-once-gate.R)
  .pt_once_server(input, session)

  # This session's background render process (R/viz-render-worker.R), from
  # RENDER_POOL (the launcher's, .pt_render_pool()) when there is one.
  # shiny::testServer()'s mock session renders in-process instead, so its tests
  # see the figure as soon as they click.
  render_worker <- if (!inherits(session, "MockShinySession"))
    if (is.null(RENDER_POOL)) .pt_render_worker() else RENDER_POOL$take()
  session$onSessionEnded(function() {
    if (is.null(render_worker)) return()
    if (is.null(RENDER_POOL)) render_worker$close()
    else RENDER_POOL$give_back(render_worker)
  })

  # selected= is always passed explicitly, as "" rather than character(0) --
  # server-mode selectize defaults to ALL_GENES's first alphabetical entry if
  # given a zero-length selected, so "" (a value absent from ALL_GENES) is
  # the only way to start with nothing shown. No gene is preselected unless the
  # preloaded start view has one.
  # A start gene named only in the GTF is not among the OrgDb's symbols.
  updateSelectizeInput(session, "gene_name",
                       choices = union(ctx$ALL_GENES, ctx$START_VIEW$gene),
                       server = TRUE, selected = ctx$START_VIEW$gene %||% "")

  # ── Mode switch: single/two/compose is one radio (input$view_mode), read
  # directly -- no derived reactiveVal/output needed, same as input$region_mode.
  # Coordinates is selectable in every sample-count mode, so there's no
  # selection to reset on a mode switch.

  state <- reactiveValues(
    gene_region          = NULL,
    gene_error           = NULL,
    gene_note            = NULL,
    svg_path             = NULL,
    rcode                = NULL,
    plot_error           = NULL,
    start_parse_err      = NULL,
    end_parse_err        = NULL,
    last_region          = NULL,
    # track labels and title of the last figure, so zoom and Apply rebuild it
    # with the live settings
    last_labels          = NULL,
    last_title           = "",
    rendering            = FALSE
  )

  # Figure settings: saved and live settings, the panel, Save and Revert
  # (R/app-settings.R).
  fs <- .pt_figure_settings_server(input, output, RENDER_CFG, SETTINGS_JSON)
  saved_cfg       <- fs$saved_cfg
  live_cfg        <- fs$live_cfg
  settings_dirty  <- fs$settings_dirty
  settings_ui_gen <- fs$settings_ui_gen
  display_keys    <- fs$display_keys

  session$onSessionEnded(function() {
    isolate({
      if (!is.null(state$svg_path) && file.exists(state$svg_path))
        unlink(state$svg_path)
    })
  })

  # Edit project settings, Reset, Clear all caches, Close app:
  # R/app-project-actions.R.
  unsaved_settings_note <- function()
    .pt_unsaved_settings_note(settings_dirty())
  .pt_project_actions_server(input, session, unsaved_settings_note, PROJECT_DIR, SETTINGS_JSON, CACHE_DIR)

  # Renders read live_cfg() (Figure settings, above); core lookups read
  # RENDER_CFG directly.

  # Coordinates: chromosome, slider and text boxes, region validity
  # (R/app-region.R).
  rg <- .pt_region_server(input, output, session, state, ctx$CHR_LENGTHS)
  chr_len_r    <- rg$chr_len_r
  region_valid <- rg$region_valid

  # Gene lookup: the gene box to a region in state (R/app-gene.R).
  .pt_gene_lookup_server(input, output, session, state, ctx, saved_cfg)

  # Gene / protein info: Output page, about the rendered figure's gene
  # (R/app-gene-info-panel.R).
  .pt_gene_info_server(output, session, state, ctx$txdb, ctx$CDS_BY_TX,
    genome = RENDER_CFG[c("orgdb_pkg", "tx_name_lookup_path",
                          "tx_meta_lookup_path", "genome_label", "cache_dir")],
    reviewed_only = reactive(input$s_uniprot_reviewed_only %||%
                               saved_cfg()$uniprot_reviewed_only))

  # Conditions: Reference/Comparison exclusion, the picked conditions, the track
  # pool (R/app-conditions.R).
  cs <- .pt_conditions_server(input, output, session, ctx)
  sel_prev             <- cs$sel_prev
  selected_conditions  <- cs$selected_conditions
  has_conditions_ready <- cs$has_conditions_ready
  figure_pool_key      <- cs$figure_pool_key

  # Stage gating: resolved region, locks, invalid tracks, Generate
  # (R/app-gate.R).
  gt <- .pt_generate_gate_server(input, output, session, state, ctx,
                                 region_valid, has_conditions_ready)
  resolved_region       <- gt$resolved_region
  has_mode              <- gt$has_mode
  invalid_figure_tracks <- gt$invalid_figure_tracks

  # Generate, Cancel, zoom and Apply to figure (R/app-render.R).
  .pt_render_server(input, output, session, state, ctx, render_worker, live_cfg,
                    resolved_region, has_mode)

  # The Output section's outputs and downloads (R/app-output.R).
  .pt_output_server(output, state)
}

shinyApp(ui = ui, server = server)

}
