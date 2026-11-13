# PlotTrellis composer: one page with a Compose section and an Output section. Figures
# render in-process with render_all(); the R script the page shows and downloads re-runs
# the same functions on the same entries and settings (assemble_download_script()).
#
# The exported entry points (run_plottrellis() and friends) are in R/entry-points.R; this file
# holds the project-settings defaults/reader and the composer app builder (.pt_composer_app()).
# The ui, the server and every helper in .pt_composer_app() share its scope as closures.

# The display settings and their defaults: how any track is drawn or filtered. These are the
# settings of the composer's Figure settings panel, and the ones render_all() fills in when
# a caller leaves them out (.pt_fill_settings()). NULL means "decided when the figure is drawn".
.pt_display_defaults <- function() list(
  font        = "DejaVu Sans",
  fontsize    = 8L,     # points
  page_w      = 9,      # inches
  chr_fmt     = "dot",  # the chromosome name on the genome axis: "dot" (chr.19), "space" (chr 19), "plain" (chr19)
  # The value axis on coverage, coverage-difference and junction tracks, which show no other
  # value scale (a Hi-C track has its colour legend, Transcript expression its own TPM axis).
  show_track_axes = TRUE,
  show_title      = TRUE,   # the figure's title line, when the figure has a title
  show_scale_bar  = TRUE,   # the genomic scale bar on a figure's first track, for the types that draw one
  hic_norm    = "KR",
  hic_style   = "rectangle",
  zrange_cap  = NULL,   # top of the Hi-C colour scale; NULL = the largest value in view
  atac_height     = 0.85,   # inches, a coverage track
  hic_height      = NULL,   # inches; NULL = half the track width (rectangle) or all of it (square)
  hic_res_mode    = "auto", # "auto" | "finest" | "fixed"
  hic_max_bins_per_view = 500, # fewest bins across the view that the automatic Hi-C resolution accepts
  hic_resolution  = NULL,   # bin size in bp for hic_res_mode "fixed"; NULL = not set
  gene_flank_pct  = 15,     # % of the gene's length shown on each side of it
  # Least pooled support an arc (a BED12 splice junction or a BEDPE loop) needs to be drawn.
  # 0 shows every arc: a BEDPE score's scale depends on the caller that wrote it, so no
  # other value is safe for every file.
  min_junction_support = 0,
  # Mean TPM a transcript must reach in a sample to be drawn and named on its own, instead
  # of being folded into "Other".
  min_display_tpm = 1,
  uniprot_reviewed_only = FALSE, # the Gene info panel's UniProt link: Swiss-Prot entries only if TRUE
  # Junction tracks: TRUE cuts an arc with one end outside the region at the region's edge,
  # FALSE draws it running past the edge. Long-range BEDPE pairs (Hi-C loops) are the case
  # that needs TRUE; splice junctions rarely reach that far.
  arc_clip_partial = FALSE,
  junction_height  = 0.85,  # inches, a junction track's arcs per sample, under its heading (a mirrored comparison: 2x)

  # ---- Two-sample comparison tracks only ----
  diff_colour_lo  = "#0033CC",  # vivid blue: lower in Comparison
  diff_colour_hi  = "#CC0000",  # vivid red: higher in Comparison
  # Mirrored junction difference: the log2 fold change drawn fully opaque. A ceiling, so one
  # sparse arc's large ratio cannot stretch the scale and wash out every other arc.
  junction_opacity_cap = 2,
  # Coverage difference: tint of the shared (overlap) layer, 0 = white to 1 = the full
  # reference colour.
  diff_shared_alpha = 0.15,
  # The two role colours of every comparison track: Reference and Comparison are always
  # drawn in these, whatever the two samples' own colours, so which one is the baseline
  # reads at a glance.
  reference_colour  = "#4D4D4D", # neutral grey
  comparison_colour = "#56B4E9"  # colour-blind-safe blue
)

# Every project setting with its default, as one flat list (reading the settings file,
# the wizard and the composer all rely on it being flat): the core import settings, which
# the wizard sets (.PT_CORE_SETTINGS), then the display settings. A function, not a
# top-level constant: app.R collates before genomes.R, so a top-level
# list(orgdb = GENOME_PRESETS$...) would fail at install.
.pt_default_settings <- function() c(
  list(
    # Blank until the wizard is given real paths, which it requires before Save.
    gtf_path     = "",
    orgdb        = GENOME_PRESETS$mm39$orgdb,
    genome_label = GENOME_PRESETS$mm39$genome_label,
    data_dir     = "",
    # The samples live in the settings file itself (samples.txt is an import/export
    # format, never read at startup). Each is list(condition =, colour =,
    # data = list(list(type =, name =, path =))), the shape import_samples_txt() and
    # export_samples_txt() use; `condition` alone identifies a sample, so a replicate is
    # a separate entry with its own condition.
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
  parsed <- tryCatch(jsonlite::fromJSON(SETTINGS_JSON, simplifyDataFrame = FALSE), error = function(e) {
    # Point at the project_samples.txt the wizard's Save keeps beside it.
    txt <- file.path(dirname(SETTINGS_JSON), "project_samples.txt")
    stop("Can't read ", SETTINGS_JSON, " (", conditionMessage(e), "). ",
         if (file.exists(txt)) paste0("Rebuild it from the project_samples.txt saved beside it: parse_samples_txt(",
                                      encodeString(txt, quote = "\""), ", project_dir = ",
                                      encodeString(dirname(SETTINGS_JSON), quote = "\""), ", overwrite = TRUE)")
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
  cfg
}

# Shared setup for every entry point (R/entry-points.R): creates project_dir and
# reads its settings merged onto the defaults.
.pt_project <- function(project_dir) {
  # tools::R_user_dir() only returns the path, it doesn't create it -- on a
  # genuinely first run (nothing has ever written here before) the wizard's
  # Save handler's writeLines(SETTINGS_JSON) fails with "No such file or
  # directory" otherwise.
  dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)
  # Absolute, so the cache/TxDb paths that end up in the render settings (and so in
  # a downloaded script, pt_save() after setwd(), ...) don't depend on the working directory.
  project_dir <- normalizePath(project_dir, winslash = "/")
  SETTINGS_JSON <- file.path(project_dir, "project_settings.json")
  DEFAULT_SETTINGS <- .pt_default_settings()
  # A project's caches go to its own cache/ folder (genome files already beside the GTF are
  # used where they are, see .genome_cache_path()).
  list(settings_json = SETTINGS_JSON, defaults = DEFAULT_SETTINGS,
       cfg = .pt_read_settings(SETTINGS_JSON, DEFAULT_SETTINGS),
       cache_dir = file.path(project_dir, "cache"))
}

# Builds (does not run) the composer app for one project. Its startup below
# (genome migration, TxDb load, ...) runs on each call.
.pt_composer_app <- function(SETTINGS_JSON, CFG) {

# CFG/SETTINGS_JSON come from .pt_project(), see R/entry-points.R.

# ── Startup: load heavy objects once ──────────────────────────────────────────
# The genome is the flat fields gtf_path, orgdb and genome_label on CFG. An older
# project_settings.json with a genome_id and no orgdb is migrated: orgdb and genome_label
# come from GENOME_PRESETS by that id and are written back, so this runs once. An id that
# isn't a preset is an error (never a silent default): fill in the Genome section by hand.
if (is.null(CFG$orgdb) && !is.null(CFG$genome_id)) {
  preset <- GENOME_PRESETS[[CFG$genome_id]]
  if (is.null(preset))
    stop("project_settings.json has an old-format genome_id ('", CFG$genome_id,
         "') that isn't a known built-in preset, so it can't be auto-migrated. ",
         "Re-open Settings and fill in the Genome section (GTF path, OrgDb, genome label) by hand.", call. = FALSE)
  CFG$genome_id    <- NULL
  CFG$orgdb        <- preset$orgdb
  CFG$genome_label <- preset$genome_label
  writeLines(jsonlite::toJSON(CFG, auto_unbox = TRUE, pretty = TRUE, null = "null"), SETTINGS_JSON)
  message("Migrated project_settings.json's old genome_id ('", preset$label, "') to flat genome fields.")
}

# Cache paths are never read from or written to CFG: every launch re-derives the
# genome's cache paths from gtf_path (ensure_genome_ready()). ORGDB_PKG below is a
# local of this app; the figure and lookup functions get the genome's packages and
# paths as arguments (RENDER_CFG, the render settings). Caches go to the project's
# cache/ folder, beside its project_settings.json; genome files already beside the GTF
# are used where they are (.genome_cache_path()).
CACHE_DIR <- file.path(dirname(SETTINGS_JSON), "cache")
GENOME <- ensure_genome_ready(gtf_path = CFG$gtf_path, orgdb = CFG$orgdb, genome_label = CFG$genome_label,
                              cache_dir = CACHE_DIR)
ORGDB_PKG           <- GENOME$orgdb

message("Loading PlotTrellis...")
# The TxDb and its CDS index come from .pt_annotation()'s session memo, so every
# render and the Gene info panel share this one copy instead of opening their own.
ANNOT    <- .pt_annotation(GENOME)
txdb     <- ANNOT$txdb

# Built once here (not per gene render) since GenomicFeatures::cdsBy() walks
# the whole TxDb regardless of which genes are actually asked for. The Gene
# info table's UniProt canonical-isoform check (uniprot_canonical.R) and the
# protein-domains track (through the same memo) look up a handful of
# transcripts' translated CDS lengths per gene out of this.
CDS_BY_TX <- ANNOT$cds_by_tx()

# Chromosome lengths: a GTF carries none, so they are never read from the TxDb's
# seqinfo. They come from UCSC through GENOME$chrom_sizes_path (build_chrom_sizes(),
# built and cached by ensure_genome_ready() above), never from a hardcoded table: one
# species' lengths would be silently wrong under any other genome.
CHR_LENGTHS <- readRDS(GENOME$chrom_sizes_path)

SAMPLES_DF  <- load_samples_from_settings(CFG$samples, CFG$data_dir)

# Every sample is selectable in every view_mode (One/Two/Multi), whatever data it has.
# Which tracks a sample, or a pair of samples, offers is decided in track_detection.R
# (detect_generic(), detect_pair_generic() and the modules' own detect_<key>()), which
# give no item rather than an error when a data type is missing.
AVAILABLE_CONDITIONS <- vapply(SAMPLES_DF, function(r) r$display, character(1))

# Every sample, with no pre-filter: detect_pair_generic() checks the data on both sides
# of each pair, so e.g. two ATAC-only samples still get an ATAC coverage comparison.
TRACK_ITEMS <- c(
  build_track_items(SAMPLES_DF),
  build_comparison_track_items(SAMPLES_DF)
)

ALL_GENES <- sort(AnnotationDbi::keys(get(ORGDB_PKG, envir = asNamespace(ORGDB_PKG)), keytype = "SYMBOL"))

# ── Helpers ────────────────────────────────────────────────────────────────────
# The one cfg every render and lookup reads, computed once from CFG: the flat CFG keys
# plus the genome-derived values (orgdb_pkg, genome_label, tx_*_lookup_path). The wizard
# edits the core settings (with a restart), and the composer's Figure settings panel
# layers live values for the other keys on top of it (saved_cfg()/live_cfg() in server()).
RENDER_CFG <- .pt_render_cfg(CFG, GENOME, CACHE_DIR)

# Bolds whatever precedes the first " · " in a label (or the whole label if it has none); format-driven, not an enumerated word list, since a track's label can start with any user-given instance name.
format_track_label <- function(label) {
  sub("^([^\u00B7]+?)( \u00B7 |$)", "<strong>\\1</strong>\\2", label)
}

# Gene query (symbol in any case, or Ensembl gene ID; see resolve_gene_query()) -> coords,
# the canonical symbol every track/panel keys on, and a user-facing error.
resolve_gene_region <- function(query, txdb, flank_pct)
  .pt_resolve_gene_region(query, txdb, flank_pct, ORGDB_PKG, ALL_GENES, GENOME$genome_label)

# The Generate button; its "generate-btn" class puts it under the page script's spinner/double-click gate, which matches that class, not an id.
generate_button_ui <- function(id, enabled) {
  if (enabled)
    actionButton(id, "Generate figure \u2192", class = "btn-success btn-lg generate-btn")
  else
    tags$button("Generate figure \u2192",
                class = "btn btn-success btn-lg generate-btn", disabled = "disabled")
}

# Shared collapsible ".plottrellis-section" panel shell; the app-wide "hide.bs.collapse"/"show.bs.collapse" JS handlers key off any ".plottrellis-section .collapse" element (not a specific id), so this same shell works for every section.
collapsible_panel <- function(section_id, title, ..., wrapper_id = NULL, wrapper_style = NULL, collapsed = FALSE) {
  div(id = wrapper_id, class = paste("panel panel-default plottrellis-section", if (collapsed) "collapsed-panel"), style = wrapper_style,
    div(class = "panel-heading",
      tags$h4(class = "panel-title",
        tags$a(href = paste0("#", section_id, "-body"), `data-toggle` = "collapse",
          tags$span(class = "collapse-arrow", if (collapsed) "\u25B8" else "\u25BE"), paste0(" ", title)
        )
      )
    ),
    div(id = paste0(section_id, "-body"), class = paste("collapse", if (!collapsed) "in"),
      div(class = "panel-body", ...)
    )
  )
}

# Each sample's own colour, else PALETTE in sample order.
DEFAULT_SAMPLE_COLS <- .pt_sample_colours(SAMPLES_DF)

# Optional preloaded starting view (pt_preload()): the UI starts with this mode,
# condition(s), gene and figure tracks already picked, so Generate is the only
# click. NULL (no start_view in the settings) keeps the usual blank start.
START_VIEW <- if (is.list(CFG$start_view)) lapply(CFG$start_view, unlist) else NULL
# Keys that no longer exist (a sample renamed since) are dropped, not shown, after
# renaming the junction track's old keys. Order is kept, and so are repeats of
# Gene annotation, the one track a figure can hold several times; any other
# track is kept once (a pt_figure() can repeat one, the picker can't).
if (!is.null(START_VIEW)) {
  sv_keys <- as.character(.pt_renamed_junctions_key(START_VIEW$tracks))
  START_VIEW$tracks <- sv_keys[sv_keys %in% names(TRACK_ITEMS) & (sv_keys == "genes" | !duplicated(sv_keys))]
}
start_condition <- function(i) {
  v <- START_VIEW$conditions
  if (length(v) >= i) v[[i]] else ""
}

# Main chromosomes first (chr1..chrN, chrX, chrY, chrM), then every
# alt/random/fix/chrUn contig alphabetically below them -- kept, not filtered.
sort_chrs <- function(chrs) {
  is_std <- base::grepl("^chr[0-9]+$", chrs)
  nums   <- suppressWarnings(as.integer(base::sub("^chr", "", chrs)))
  sex_mt <- intersect(c("chrX", "chrY", "chrM"), chrs)
  c(chrs[is_std][order(nums[is_std])], sex_mt, sort(setdiff(chrs[!is_std], sex_mt)))
}
fmt_num <- function(x, digits = 2) {
  s <- formatC(round(x, digits), format = "f", digits = digits)
  s <- sub("0+$", "", s); sub("\\.$", "", s)
}
format_bp <- function(bp) {
  if (bp >= 1e6) paste0(fmt_num(bp / 1e6), "Mb")
  else if (bp >= 1e3) paste0(fmt_num(bp / 1e3), "Kb")
  else paste0(bp, "bp")
}
format_coord <- function(n) {
  if (n >= 1e6) paste0(fmt_num(n / 1e6, 3), "Mb")
  else if (n >= 1e3) paste0(fmt_num(n / 1e3, 3), "Kb")
  else paste0(n, " bp")
}

# Parse human-readable bp strings: "17000000", "17Mb", "17.56Mb", "15kb", "150KB"
parse_bp <- function(s) {
  s <- trimws(s)
  if (grepl("^[0-9]+$", s)) return(as.integer(s))
  m <- regmatches(s, regexec("^([0-9]*\\.?[0-9]+)\\s*(kb|mb)$", s, ignore.case = TRUE))[[1]]
  if (length(m) < 3) return(NA_integer_)
  num  <- as.numeric(m[2])
  unit <- toupper(m[3])
  if (unit == "KB") as.integer(round(num * 1e3))
  else if (unit == "MB") as.integer(round(num * 1e6))
  else NA_integer_
}

# Picker keys -> render_all() entries (.pt_build_entries()).
build_entries <- function(figure_labels, col_map, settings)
  .pt_build_entries(figure_labels, TRACK_ITEMS, col_map, settings)

# ── Conditions and track pool ─────────────────────────────────────────────────

track_labels <- function(keys) {
  setNames(lapply(keys, function(k) HTML(format_track_label(TRACK_ITEMS[[k]]$label))), keys)
}

# Alphabetical (case-insensitive) by TRACK_ITEMS' label, split on " · ", so tracks group by assay name, then by sample; never applied to figure_keys, whose order is deliberately curated.
sort_keys_by_label <- function(keys) {
  label_parts <- function(label) {
    label <- tolower(label)
    idx <- regexpr(" \u00B7 ", label, fixed = TRUE)
    if (idx == -1) return(c(label, ""))
    c(substr(label, 1, idx - 1),
      substr(label, idx + attr(idx, "match.length"), nchar(label)))
  }
  parts <- vapply(keys, function(k) label_parts(TRACK_ITEMS[[k]]$label), character(2))
  keys[order(parts[1, ], parts[2, ])]
}

# Every sample-count mode starts with just the gene-annotation track (or the preloaded start view's tracks), a plain, predictable starting point; anything else goes in from Available tracks explicitly.
default_figure_keys <- function(conditions)
  if (!is.null(START_VIEW) && identical(conditions, START_VIEW$conditions)) START_VIEW$tracks else "genes"

# Item pool + default figure/available split for 1-2 conditions, filtered by each item's own `sample` field; fixed_item()s are found by scanning TRACK_ITEMS itself, never a hardcoded type list, so a new sample-free track type shows up here without an edit.
track_split <- function(conditions) {
  cond_items <- Filter(function(sp) isTRUE(sp$sample %in% conditions), TRACK_ITEMS)
  comparison_keys <- character(0)
  if (length(conditions) > 1) {
    comparison_keys <- names(TRACK_ITEMS)[vapply(TRACK_ITEMS, function(sp) {
      !is.null(sp$samples) && length(sp$samples) == 2 && identical(sp$samples, conditions)
    }, logical(1))]
  }
  # is.null(sp$samples) is required alongside is.na(sp$sample) here -- every
  # comparison item (detect_pair_generic()) also sets sample = NA (it uses
  # `samples`, plural, instead), so without this every comparison track of
  # EVERY condition pair (not just the one selected above) would land in
  # "fixed" and appear whatever the conditions.
  fixed_keys  <- names(Filter(function(sp) isTRUE(is.na(sp$sample)) && is.null(sp$samples), TRACK_ITEMS))
  pool_keys   <- c(names(cond_items), comparison_keys, fixed_keys)
  figure_keys <- default_figure_keys(conditions)
  figure_keys <- figure_keys[figure_keys %in% pool_keys]
  # Gene annotation stays in Available whatever Figure tracks holds: a figure can
  # hold it many times (the composer JS keeps exactly one in Available).
  avail_keys  <- sort_keys_by_label(union(setdiff(pool_keys, figure_keys), intersect("genes", pool_keys)))
  list(available = track_labels(avail_keys), figure = track_labels(figure_keys))
}


# Multi's "Available tracks" pool: TRACK_ITEMS minus EVERY two-sample comparison
# item: a comparison is a Reference-vs-Comparison concept that Two-sample mode
# offers for its own pair, and offering every pair's again here would only be
# clutter. Gene annotation is in it too, as in every mode (Figure tracks also
# starts with one).
# Sorted alphabetically by track label (sort_keys_by_label()) so the "All" filter
# shows a stable, name-grouped order (ATAC vs RNA-seq vs Junctions, etc.)
# rather than TRACK_ITEMS' own sample-then-type insertion order.
COMPOSE_AVAILABLE_KEYS <- sort_keys_by_label(Filter(function(k) {
  sp <- TRACK_ITEMS[[k]]
  is.null(sp$samples) || length(sp$samples) != 2
}, names(TRACK_ITEMS)))

# Choices for the Condition/Reference/Comparison selects: every AVAILABLE_CONDITIONS entry is always offered plus a blank placeholder; collisions (picking the value already sat in the other box) are resolved by clearing the other box server-side, not by removing options here.
condition_choices <- function() {
  c(setNames("", "(select)"), setNames(AVAILABLE_CONDITIONS, AVAILABLE_CONDITIONS))
}

inline_svg <- function(svg_path) {
  lines    <- readLines(svg_path, warn = FALSE)
  svg_text <- paste(lines, collapse = "\n")
  svg_text <- base::sub("<svg ", "<svg style='width:100%;height:auto;' ", svg_text)
  # Capped + centred, not full-bleed: navbarPage is fluid, so an uncapped
  # width:100% stretches the SVG to however wide the monitor is -- since
  # height:auto preserves aspect ratio, that also inflates the height,
  # pushing tracks near the bottom (genes/transcripts) well past one
  # viewport. Capping keeps the figure at a legible, comparable size
  # regardless of screen width; narrower screens still shrink it via the
  # width:100% below the cap. No max-height/overflow-y here, deliberately: a
  # capped-height scroll box makes tall figures feel clipped/hidden.
  div(style = "width:100%; max-width:1400px; margin:0 auto; overflow-x:auto;", HTML(svg_text))
}

# JS map from track key → type, embedded in page for client-side filter
TRACK_TYPE_JS <- jsonlite::toJSON(
  setNames(lapply(TRACK_ITEMS, function(x) x$type), names(TRACK_ITEMS)),
  auto_unbox = TRUE
)
# JS map from track key -> broken: a declared file that is missing (the `broken`
# field track_detection.R sets on per-sample items, checked once at startup) gets
# red-boxed in the track picker, the same look as the gene-required flag.
# isTRUE() turns the common case (an item with no `broken` field at all, e.g.
# genes, transcripts, comparison tracks) from NULL into FALSE, rather than
# needing every module to declare it.
TRACK_BROKEN_JS <- jsonlite::toJSON(
  setNames(lapply(TRACK_ITEMS, function(x) isTRUE(x$broken)), names(TRACK_ITEMS)),
  auto_unbox = TRUE
)
# JS map from track type -> the id of the filter button that lists it, from each
# module's info_<key>() `filter` (.pt_type_filters()): e.g. the comparison types share
# "comparison" ("Comparisons", Two-sample mode only), and the gene annotation and
# transcript-structure tracks share "annotation" ("Annotations").
TRACK_FILTER_JS <- jsonlite::toJSON(.pt_type_filters(), auto_unbox = TRUE)
# Gene-required track types (needs_gene in their info_<key>()): from coordinates
# they can only draw an empty-state message. Kept as a plain R vector too (not just the JSON string below) so
# server() can run the same red-box check the client already does, to gate Generate
# on it.
GENE_REQUIRED_TYPES <- .pt_gene_required_types()
GENE_REQUIRED_TYPES_JS <- jsonlite::toJSON(GENE_REQUIRED_TYPES)

# ── UI ─────────────────────────────────────────────────────────────────────────
ui <- navbarPage(
  title = "PlotTrellis",
  id    = "main_tabs",

  tags$head(tags$style(HTML("
    :root {
      --bio-font-base:    15px;
      --bio-font-small:   13px;
      --bio-font-label:   16px;
      --bio-font-heading: 19px;
      --bio-font-code:    12px;
    }
    body  { font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif;
            font-size: var(--bio-font-base); }
    label { font-size: var(--bio-font-label); }
    h4    { font-size: var(--bio-font-heading); font-weight: 600;
            margin-top: 8px; margin-bottom: 6px; }
    .shiny-input-container .radio label,
    .shiny-input-container .checkbox label { font-size: var(--bio-font-label); }
    .shiny-input-container .radio input[type='radio'] {
      transform: scale(1.35); transform-origin: center; margin-right: 4px; cursor: pointer;
    }
    .zoom-group-label { font-size: var(--bio-font-label); font-weight:600; color:#555;
                        margin-bottom:5px; letter-spacing:0.03em; }
    .rank-list-container { min-height: 60px; max-height: 42vh; overflow-y: auto; }
    .rank-list-item { cursor: grab; }
    .rank-list-item:active { cursor: grabbing; }
    .track-filters { margin-bottom:10px; display:flex; flex-wrap:wrap; gap:4px; }
    .track-filters .btn { flex:1 1 0; margin-right:0; }
    .track-filters .btn.active { background-color:#337ab7; color:#fff; border-color:#2e6da4; }
    /* Gene-dependent tracks (Transcripts (filtered)/Transcripts/Transcripts
       annotated (filtered)/Canonical protein annotation/Transcript expression/
       Junctions mirrored diff (TPM-scaled), see GENE_REQUIRED_TYPES_JS)
       while region mode is Coordinates, not Gene name -- flagged, not hidden,
       since a track can already be sitting in
       Figure tracks when the mode switches away from Gene name; Generate stays
       disabled while one is there. Applies in every sample-count mode:
       applyGeneRequiredFlag() scopes to the single #tracks-block shared by
       all three modes. */
    /* A declared file that is missing (moved/deleted) gets the same red box
       as the gene-required flag above, in every view mode's picker. */
    .rank-list-item.track-gene-required,
    .rank-list-item.track-broken {
      background-color: #fdecea !important;
      border: 1px dashed #d9534f !important;
      opacity: 0.7;
    }
    #loading-overlay {
      position: fixed; top: 0; left: 0; width: 100%; height: 100%;
      background: rgba(255,255,255,0.78); z-index: 9999;
      text-align: center; padding-top: 38vh; box-sizing: border-box;
    }
    .plottrellis-spinner {
      display: inline-block; width: 54px; height: 54px;
      border: 5px solid #ddd; border-top-color: #337ab7;
      border-radius: 50%; animation: plottrellis-spin 0.75s linear infinite;
    }
    @keyframes plottrellis-spin { to { transform: rotate(360deg); } }
    /* Single-page layout: hide the nav tab links, keep the PlotTrellis title bar */
    .navbar-nav { display: none !important; }
    /* Collapsible section panels */
    .plottrellis-section { margin-bottom: 12px; border-radius: 6px; }
    .plottrellis-section > .panel-heading {
      background-color: #f0f4f8; border-bottom: 1px solid #dce3ea;
      border-radius: 6px 6px 0 0; padding: 8px 18px; cursor: pointer;
    }
    .plottrellis-section.collapsed-panel > .panel-heading {
      border-radius: 6px;
    }
    .plottrellis-section .panel-heading a {
      color: #2c3e50; text-decoration: none; font-size: 16px;
      font-weight: 600; display: block;
    }
    .plottrellis-section .panel-heading a:hover { color: #337ab7; }
    .plottrellis-section > .panel-body { padding: 12px 18px; }
    hr { margin-top: 10px; margin-bottom: 10px; }
    #topbar-actions {
      position: fixed; top: 10px; right: 15px; z-index: 9998;
      display: flex; align-items: center; gap: 8px;
    }
    /* region_mode radios (Gene / Coordinates), dimmed by
       setConditionGateLocked() while mode+conditions aren't both resolved
       yet. */
    .shiny-input-container .radio-inline.disabled-radio-label {
      opacity: 0.45; cursor: not-allowed;
    }
    /* Track picker (Available/Figure rank-lists), dimmed by setPickerLocked()
       while mode+conditions+region aren't all resolved yet. */
    .rank-list-locked {
      opacity: 0.45; cursor: not-allowed;
    }
    /* Add all / Clear all, at the bottom of each list's box (addListButtons()),
       kept in view when a long list scrolls the box. */
    .track-list-action { align-self: center; margin-top: 8px; position: sticky; bottom: 0; }
  "))),
  tags$head(tags$script(HTML(paste0("
    var TRACK_TYPE_MAP = ", TRACK_TYPE_JS, ";
    var TRACK_BROKEN_MAP = ", TRACK_BROKEN_JS, ";
    var TRACK_FILTER_MAP = ", TRACK_FILTER_JS, ";
    var GENE_REQUIRED_TYPES = ", GENE_REQUIRED_TYPES_JS, ";

    // Track-filter buttons + bucket are scoped per [data-available-id] block.
    function applyTrackFilter(blockEl) {
      var $block  = $(blockEl);
      var availId = $block.attr('data-available-id');
      var f = $block.find('.filter-btn.active').data('filter') || 'all';
      $block.find('#' + availId + ' .rank-list-item').each(function() {
        var key  = $(this).attr('data-rank-id') || '';
        var type = TRACK_TYPE_MAP[key] || '';
        var matches = f === 'all' || TRACK_FILTER_MAP[type] === f;
        if (matches) {
          $(this).show();
        } else {
          $(this).hide();
        }
      });
    }

    function applyAllTrackFilters() {
      $('[data-available-id]').each(function() { applyTrackFilter(this); });
    }

    // Flags, rather than hides, gene-required tracks while the region mode is
    // Coordinates: one may already be in Figure tracks. A track that is also
    // broken keeps applyBrokenFlag()'s tooltip, whichever of the two runs first.
    function applyGeneRequiredFlag() {
      var isCoords = $('input[name=\"region_mode\"]:checked').val() === 'coords';
      $('#tracks-block .rank-list-item').each(function() {
        var key      = $(this).attr('data-rank-id') || '';
        var type     = TRACK_TYPE_MAP[key] || '';
        var flagIt   = isCoords && GENE_REQUIRED_TYPES.indexOf(type) !== -1;
        $(this).toggleClass('track-gene-required', flagIt);
        if (TRACK_BROKEN_MAP[key]) return;
        $(this).attr('title', flagIt ? \"Needs 'Gene name' region mode -- won't render from coordinates\" : '');
      });
    }

    // Every .rank-list-item in the block, Figure tracks included (a broken track
    // can already be placed there), not just the Available list. Never cleared:
    // TRACK_BROKEN_MAP is fixed at page load.
    function applyBrokenFlag() {
      $('[data-available-id]').each(function() {
        $(this).find('.rank-list-item').each(function() {
          var key    = $(this).attr('data-rank-id') || '';
          var broken = TRACK_BROKEN_MAP[key] || false;
          if (broken) {
            $(this).addClass('track-broken');
            $(this).attr('title', \"This track's file path is missing or unreadable\");
          }
        });
      });
    }

    // Locks the region_mode radios and the gene box until a sample-count mode and
    // its condition(s) are chosen. Driven by the server's condition_lock message,
    // since whether the conditions are ready depends on inputs this script has
    // no simple DOM check for. The Coordinates panel needs no lock of its own: it
    // can't be reached until region_mode can be clicked.
    function setConditionGateLocked(locked) {
      var $regionModeInputs = $('input[name=\"region_mode\"]');
      $regionModeInputs.prop('disabled', locked);
      $regionModeInputs.closest('.radio-inline').toggleClass('disabled-radio-label', locked);
      var geneSel = document.getElementById('gene_name');
      if (geneSel && geneSel.selectize) {
        if (locked) geneSel.selectize.disable(); else geneSel.selectize.enable();
      }
    }

    // The track picker locks (both Available and Figure tracks, not just
    // visually dimmed) whenever the server-side region_lock message says the
    // mode+conditions+region gate isn't satisfied yet -- see the
    // `region_lock` custom message handler below for how this actually
    // disables dragging via each rank-list's own Sortable instance.
    function setPickerLocked(locked) {
      setupRepeatableTrack();
      addListButtons();
      $('#add_all_tracks, #clear_all_tracks').prop('disabled', locked);
      ['available_tracks', 'figure_tracks'].forEach(function(id) {
        var el = document.getElementById(id);
        if (!el) return;
        var inst = window.Sortable && window.Sortable.get(el);
        if (inst) inst.option('disabled', locked);
        $(el).toggleClass('rank-list-locked', locked);
      });
    }

    // Gene annotation ('genes') can sit in Figure tracks any number of times, every
    // one the same track, while Available always holds exactly one: dragging it
    // out of Available leaves one behind (Sortable's clone pull), and one dropped
    // on Available that already has one is removed (a drop's 'add' event runs
    // before the list reports itself in 'sort'). bucket_list() gives each list
    // its group as a plain name, so the clone pull is set here on the live
    // Sortable instance, once (update_rank_list() keeps the instances).
    function isRepeatableTrack(el) { return el.getAttribute('data-rank-id') === 'genes'; }
    function holdsRepeatableTrack(list, except) {
      return Array.prototype.some.call(list.children, function(c) { return c !== except && isRepeatableTrack(c); });
    }
    function setupRepeatableTrack() {
      var avail = document.getElementById('available_tracks');
      var inst = avail && window.Sortable && Sortable.get(avail);
      if (!inst || inst.ptRepeatable) return;
      inst.ptRepeatable = true;
      inst.option('group', { name: inst.options.group.name,
        pull: function(to, from, dragEl) { return isRepeatableTrack(dragEl) ? 'clone' : true; }, put: true });
      var onAdd = inst.options.onAdd;
      inst.option('onAdd', function(evt) {
        if (isRepeatableTrack(evt.item) && holdsRepeatableTrack(avail, evt.item)) evt.item.parentNode.removeChild(evt.item);
        if (onAdd) onAdd.call(this, evt);
      });
    }

    // Moves tracks of one list, in the order given, to the bottom of the other
    // (Available -> Figure, Figure -> Available), then lets both lists report
    // themselves through sortable's own onSort handlers, as after a drag (and
    // set sortable's empty-list class, which a drag sets on move). A locked
    // picker ignores it, like drag. Gene annotation as a drag does: a copy goes
    // to Figure tracks, and one sent back to an Available that has one is removed.
    // Returns the last track placed, or null.
    function moveTracks(fromId, items) {
      var from = document.getElementById(fromId);
      var to = document.getElementById(fromId === 'available_tracks' ? 'figure_tracks' : 'available_tracks');
      var src = from && window.Sortable && Sortable.get(from), dst = to && window.Sortable && Sortable.get(to);
      if (!src || !dst || src.option('disabled') || dst.option('disabled') || items.length === 0) return null;
      var last = null;
      items.forEach(function(el) {
        var item = el;
        if (isRepeatableTrack(el) && fromId === 'available_tracks') item = el.cloneNode(true);
        else if (isRepeatableTrack(el) && holdsRepeatableTrack(to, el)) item = null;
        if (item) { to.appendChild(item); last = item; } else from.removeChild(el);
      });
      [from, to].forEach(function(el) { Sortable.utils.toggleClass(el, 'rank-list-empty', el.children.length === 0); });
      [src, dst].forEach(function(s) { if (s.options.onSort) s.options.onSort.call(s, {}); });
      return last;
    }

    // Double-click moves one track.
    $(document).on('dblclick', '#available_tracks .rank-list-item, #figure_tracks .rank-list-item', function() {
      if (window.getSelection) window.getSelection().removeAllRanges();
      var placed = moveTracks(this.parentNode.id, [this]);
      if (placed) placed.scrollIntoView({ block: 'nearest' });
    });

    // 'Add all' at the bottom of Available tracks (every track its filter shows,
    // in its order) and 'Clear all' at the bottom of Figure tracks (every track),
    // each centred inside its list's box, below the list: update_rank_list()
    // replaces only the list, so the buttons stay. Same moves as double-click.
    // Disabled until the server's first region_lock message (setPickerLocked()).
    function addListButtons() {
      [['available_tracks', 'add_all_tracks', 'Add all'], ['figure_tracks', 'clear_all_tracks', 'Clear all']].forEach(function(b) {
        var box = document.getElementById('rank-list-' + b[0]);
        if (!box || document.getElementById(b[1])) return;
        var btn = document.createElement('button');
        btn.type = 'button';
        btn.id = b[1];
        btn.className = 'btn btn-default track-list-action';
        btn.textContent = b[2];
        btn.disabled = true;
        box.appendChild(btn);
      });
    }
    $(document).on('click', '#add_all_tracks', function() {
      var shown = Array.prototype.filter.call(document.getElementById('available_tracks').children,
                                              function(el) { return el.style.display !== 'none'; });
      moveTracks('available_tracks', shown);
    });
    $(document).on('click', '#clear_all_tracks', function() {
      moveTracks('figure_tracks', Array.prototype.slice.call(document.getElementById('figure_tracks').children));
    });

    $(document).ready(function() {
      addListButtons();
      setTimeout(function() { setupRepeatableTrack(); applyAllTrackFilters(); applyGeneRequiredFlag(); applyBrokenFlag(); }, 400);
    });
    $(document).on('shiny:inputchanged', function(e) {
      if (e.name === 'available_tracks' || e.name === 'figure_tracks') {
        setTimeout(applyAllTrackFilters, 80);
      }
      if (e.name === 'available_tracks' || e.name === 'figure_tracks' || e.name === 'region_mode') {
        setTimeout(applyGeneRequiredFlag, 80);
        setTimeout(applyBrokenFlag, 80);
      }
    });
    $(document).on('click', '.filter-btn', function() {
      var $block = $(this).closest('[data-available-id]');
      $block.find('.filter-btn').removeClass('active');
      $(this).addClass('active');
      applyTrackFilter($block[0]);
    });

    // A button of comparisons only (Comparisons) exists in two-sample mode only
    // (see the conditionalPanel around it, .pt_filter_buttons()) -- if it was active
    // when the user switches to single/compose, it's now hidden with no way to click
    // back to 'All', so reset the filter for them instead of leaving it stuck.
    $(document).on('change', 'input[name=\"view_mode\"]', function() {
      if ($(this).val() === 'two') return;
      var $block = $('#tracks-block');
      var $cmp   = $block.find('.shiny-panel-conditional .filter-btn');
      if ($cmp.hasClass('active')) {
        $cmp.removeClass('active');
        $block.find('.filter-btn[data-filter=\"all\"]').addClass('active');
        applyTrackFilter($block[0]);
      }
    });

    // Tooltips: init on load and after Shiny re-renders
    $(function() { $('[data-toggle=\"tooltip\"]').tooltip(); });
    $(document).on('shiny:value', function() {
      setTimeout(function() { $('[data-toggle=\"tooltip\"]').tooltip(); }, 250);
    });

    // -- Loading overlay --
    // Injected into <body> directly so position:fixed works as a true viewport overlay
    // (navbarPage wraps non-tabPanel children as <li> items which breaks fixed positioning)
    $(document).ready(function() {
      $('body').append(
        '<div id=\"loading-overlay\" style=\"display:none;\">' +
        '<div class=\"plottrellis-spinner\"></div>' +
        '<p style=\"margin-top:12px;font-size:14px;color:#555;\">Rendering\\u2026</p>' +
        '</div>'
      );
    });

    var plottrellisRendering = false;

    var ZOOM_BTN_IDS = ['zoom_in_15x','zoom_in_3x','zoom_in_10x',
                        'zoom_out_15x','zoom_out_3x','zoom_out_10x'];

    function showSpinner() {
      plottrellisRendering = true;
      document.getElementById('loading-overlay').style.display = 'block';
      ZOOM_BTN_IDS.forEach(function(id) {
        var el = document.getElementById(id);
        if (el) { el.disabled = true; el.classList.add('disabled'); }
      });
    }

    function hideSpinner() {
      if (!plottrellisRendering) return;
      plottrellisRendering = false;
      document.getElementById('loading-overlay').style.display = 'none';
      ZOOM_BTN_IDS.forEach(function(id) {
        var el = document.getElementById(id);
        if (el) { el.disabled = false; el.classList.remove('disabled'); }
      });
    }

    // Capture-phase gate: runs before Shiny's own click binding (which is
    // bound directly on the button), so a click that arrives while a render
    // is already in flight never reaches Shiny at all, it can't queue up
    // and fire a second render after the first one finishes.
    document.addEventListener('click', function(e) {
      // .generate-btn, not a hardcoded #generate id -- Generate and Apply to
      // figure both carry this class, so neither render skips the
      // spinner/double-click gate.
      var gateSelector = '.generate-btn, #' + ZOOM_BTN_IDS.join(', #');
      var el = e.target.closest(gateSelector);
      if (!el) return;
      if (plottrellisRendering) {
        e.stopImmediatePropagation();
        e.preventDefault();
        return;
      }
      showSpinner();
    }, true);

    // Primary signal: Shiny has delivered new SVG or error content to the DOM.
    // Double rAF gives the browser one paint cycle so SVG is visible before we uncover it.
    $(document).on('shiny:value', function(e) {
      if (!plottrellisRendering) return;
      if (e.name === 'svg_display_ui' || e.name === 'render_error_ui') {
        requestAnimationFrame(function() { requestAnimationFrame(hideSpinner); });
      }
    });

    // Fallback: explicit signal from R (fires after on.exit in do_render)
    Shiny.addCustomMessageHandler('hide_spinner', function(x) { hideSpinner(); });

    // Server-side region_lock drives the track
    // picker's lock state (see setPickerLocked() above and the
    // has_mode()/has_conditions_ready()/resolved_region()-driven observe()
    // in server()).
    Shiny.addCustomMessageHandler('region_lock', function(x) { setPickerLocked(x.locked); });

    // Server-side condition_lock drives the
    // region_mode/gene_name lock state (see setConditionGateLocked() above
    // and the has_mode()/has_conditions_ready()-driven observe() in
    // server()).
    Shiny.addCustomMessageHandler('condition_lock', function(x) { setConditionGateLocked(x.locked); });

    // Safety net: any uncaught Shiny error also dismisses the spinner
    $(document).on('shiny:error', function() { hideSpinner(); });

    // Show output section, expand it, scroll to it after first render.
    Shiny.addCustomMessageHandler('show_output', function(x) {
      var sec = document.getElementById('output-section');
      if (sec) sec.style.display = 'block';
      $('#output-body').collapse('show');
      setTimeout(function() {
        if (sec) sec.scrollIntoView({behavior: 'smooth'});
      }, 350);
    });

    // Keep collapse arrow in sync with open/close state
    $(document).on('hide.bs.collapse', '.plottrellis-section .collapse', function() {
      $(this).closest('.plottrellis-section').find('.collapse-arrow').text('\\u25B8');
      $(this).closest('.plottrellis-section').addClass('collapsed-panel');
    });
    $(document).on('show.bs.collapse', '.plottrellis-section .collapse', function() {
      $(this).closest('.plottrellis-section').find('.collapse-arrow').text('\\u25BE');
      $(this).closest('.plottrellis-section').removeClass('collapsed-panel');
    });

    // Coord text inputs: only fire Shiny update on blur or Enter (not every keystroke)
    $(document).on('shiny:bound', function(e) {
      var id = e.target.id;
      if (id !== 'region_start_text' && id !== 'region_end_text') return;
      var el = e.target;
      $(el).off('.textInputBinding');  // remove default keyup/input reactive listeners
      $(el).on('blur.coordInput', function() {
        Shiny.setInputValue(id, el.value, {priority: 'event'});
      });
      $(el).on('keydown.coordInput', function(evt) {
        if (evt.which === 13) {
          evt.preventDefault();
          Shiny.setInputValue(id, el.value, {priority: 'event'});
        }
      });
    });
  ")))),

  tabPanel("",

    # Edit settings / clear caches / reset / close buttons, fixed to the top-right
    # navbar area via CSS. Reset is its own button: it is a separate, more
    # destructive action than editing settings.
    div(id = "topbar-actions",
      actionButton("open_settings", "Edit project settings \u25B8",
                   class = "btn-default btn-sm"),
      actionButton("clear_caches_btn", "Clear all caches...",
                   class = "btn-warning btn-sm"),
      actionButton("reset_project_btn", "Reset project configuration...",
                   class = "btn-danger btn-sm"),
      actionButton("close_app_btn", "Close app",
                   class = "btn-default btn-sm")
    ),

    # ── Compose: one 3-way mode radio ─────────────────────────────────────────────
    # "single"/"two" (One/Two) pick 1-2 conditions from dropdowns and narrow the
    # track pool to them; "compose" (Multi) offers every sample's tracks, with the
    # track picker as its only sample selector.
    collapsible_panel("compose", "Compose",

      h4("Number of samples"),
      radioButtons("view_mode", NULL,
                   choices  = c("One" = "single",
                                "Two" = "two",
                                "Multi" = "compose"),
                   selected = START_VIEW$view_mode %||% character(0), inline = TRUE),
      hr(),

      # Condition selection -- single/two-sample modes only; track-compose has
      # no dropdown, the track picker below is the only sample selector there.
      # Placed above Genomic region so Reference/Comparison (or Condition) is
      # chosen before searching for a gene. All start unselected:
      # region_mode/gene_name and the track picker stay locked until these
      # resolve too, not just until a mode is chosen (see
      # has_conditions_ready() in server()).
      conditionalPanel("input.view_mode === 'single'",
        fluidRow(column(4,
          selectizeInput("sel_condition", "Condition:",
                      choices = condition_choices(),
                      selected = if (identical(START_VIEW$view_mode, "single")) start_condition(1) else "")
        ))
      ),
      conditionalPanel("input.view_mode === 'two'",
        fluidRow(
          column(4,
            selectizeInput("sel_reference", "Reference:",
                        choices = condition_choices(),
                        selected = if (identical(START_VIEW$view_mode, "two")) start_condition(1) else "")
          ),
          column(4,
            selectizeInput("sel_comparison", "Comparison:",
                        choices = condition_choices(),
                        selected = if (identical(START_VIEW$view_mode, "two")) start_condition(2) else "")
          )
        )
      ),

      # Region (full width)
      fluidRow(
        column(12,
          h4("Genomic region"),
          radioButtons("region_mode", NULL,
                       choices  = c("Gene name" = "gene", "Coordinates" = "coords"),
                       selected = if (is.null(START_VIEW$gene)) character(0) else "gene", inline = TRUE),
          hr(style = "margin:8px 0;"),

          conditionalPanel("input.region_mode === 'gene'",
            fluidRow(
              column(4,
                selectizeInput("gene_name", "Gene symbol or Ensembl ID:",
                               choices = NULL,
                               # maxOptions caps how many server-side search
                               # matches selectize will ever render, not a
                               # "load more on scroll" limit -- a small cap
                               # silently truncates most typed queries well
                               # before the user could scroll to a match.
                               # create = TRUE lets an Ensembl ID (not in the symbol list) be
                               # entered; resolve_gene_region() maps it to its symbol.
                               options = list(placeholder = "e.g. Pax5", maxOptions = 1000, create = TRUE))
              ),
              column(6, br(), uiOutput("gene_result_ui"))
            )
          ),

          # Coordinates panel, shown whenever region_mode is 'coords' -- selectable
          # in every sample-count mode; tracks that need a gene symbol can't draw
          # from coordinates and are red-boxed in the picker.
          conditionalPanel("input.region_mode === 'coords'",
            fluidRow(
              column(3, uiOutput("chr_select_ui")),
              column(9, uiOutput("range_slider_ui"))
            ),
            fluidRow(
              column(3, uiOutput("start_num_ui")),
              column(3, uiOutput("end_num_ui")),
              column(6, br(), uiOutput("region_error_ui"))
            )
          )
        )
      ),

      hr(),

      # Track assignment (full width)
      fluidRow(
        column(12,
          h4("Track assignment"),

          div(id = "tracks-block", `data-available-id` = "available_tracks",
            div(class = "track-filters",
              tags$button("All",         class = "btn btn-lg btn-default filter-btn active", `data-filter` = "all"),
              # One button per filter group the modules declare (info_<key>()'s `filter`).
              # "Comparisons" shows in Two-sample mode only: one-sample has none, and
              # track-compose's pool (COMPOSE_AVAILABLE_KEYS) excludes every comparison
              # item outright, so there is nothing for it to select there either.
              .pt_filter_buttons()
            ),

            bucket_list(
              header      = NULL,
              group_name  = "tracks",
              orientation = "horizontal",
              add_rank_list(
                text     = "Available tracks (drag \u2192 to figure)",
                labels   = track_labels(COMPOSE_AVAILABLE_KEYS),
                input_id = "available_tracks"
              ),
              add_rank_list(
                text     = "Figure tracks (drag to reorder)",
                labels   = list(genes = HTML(format_track_label(TRACK_ITEMS[["genes"]]$label))),
                input_id = "figure_tracks"
              )
            )
          )
        )
      ),

      # Every non-core setting, live for this session. Changes apply
      # on Generate, zoom or Apply to figure (no automatic re-render).
      collapsible_panel("figure_settings", "Figure settings", collapsed = TRUE,
        .pt_plot_title_input(),
        p(class = "text-muted", "Changes apply on the next Generate, zoom or Apply to figure, for this session only. Save as project default keeps them for next time."),
        uiOutput("display_settings_ui"),
        tags$div(style = "margin-top:10px; display:flex; gap:8px; align-items:center; flex-wrap:wrap;",
          uiOutput("apply_settings_ui", inline = TRUE),
          actionButton("save_settings", "Save as project default", class = "btn-default"),
          actionButton("revert_settings", "Revert to saved", class = "btn-default"),
          uiOutput("settings_dirty_ui", inline = TRUE)
        )
      ),
      tags$div(style = "height:8px;"),
      fluidRow(
        column(4, uiOutput("generate_ui")),
        column(8, uiOutput("generate_warning_ui"))
      ),
      tags$div(style = "height:16px;")
    ), # /collapsible_panel: compose

    # ── Output section (hidden until first render) ───────────────────────────────
    # Full-width, no inner scrollbox -- the figure renders at its natural size
    # and the whole page scrolls to reveal more of it (a max-height/overflow-y
    # box makes tall figures feel clipped/hidden).
    collapsible_panel("output", "Output",
      uiOutput("render_error_ui"),
      uiOutput("zoom_ui"),
      uiOutput("svg_display_ui"),
      tags$div(style = "margin-top:10px;",
        downloadButton("dl_svg", "Download SVG", class = "btn-lg"),
        tags$span(style = "margin-left:6px;",
          `data-toggle` = "tooltip", `data-placement` = "top",
          title = "The exact R script used to produce this figure. Run it in a clean R session to reproduce identically.",
          downloadButton("dl_rcode", "Download R script", class = "btn-lg")
        )
      ),
      .pt_gene_info_ui(collapsible_panel),
      tags$div(style = "margin-top:12px;",
        collapsible_panel("rcode", "R script", collapsed = TRUE,
          div(style = "max-height: 60vh; overflow-y: auto; font-size: var(--bio-font-code);",
            verbatimTextOutput("rcode_output"))
        )
      ),
      wrapper_id = "output-section", wrapper_style = "display:none;"
    )

  ) # /tabPanel
)

# ── Server ─────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  # selected= is always passed explicitly, as "" rather than character(0) --
  # server-mode selectize defaults to ALL_GENES's first alphabetical entry if
  # given a zero-length selected, so "" (a value absent from ALL_GENES) is
  # the only way to start with nothing shown. No gene is preselected unless the
  # preloaded start view has one.
  updateSelectizeInput(session, "gene_name", choices = ALL_GENES, server = TRUE, selected = START_VIEW$gene %||% "")

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
    last_labels          = NULL,  # track labels and title of the last figure, so zoom
    last_title           = "",    # and Apply rebuild it with the live settings
    rendering            = FALSE
  )

  # ── Figure settings ─────────────────────────────────────────────────────────
  # saved_cfg: the project's saved settings (RENDER_CFG, then whatever Save as
  # project default writes). live_cfg: saved_cfg with the composer panel's
  # s_<key> inputs on top. Every render reads live_cfg(); core keys (genome,
  # data root, samples) and the derived paths always come from saved_cfg.
  # Assigned with x[names(v)] <- v, not modifyList(), which would drop the
  # NULL ("auto") values of zrange_cap/hic_height.
  saved_cfg <- reactiveVal(RENDER_CFG)
  settings_ui_gen <- reactiveVal(0)  # bumped to re-render the panel from saved_cfg()
  display_keys <- setdiff(names(.pt_default_settings()), .PT_CORE_SETTINGS)
  live_cfg <- reactive({
    x <- saved_cfg()
    v <- .pt_collect_display_settings(input, "s_", x)
    x[names(v)] <- v
    x
  })
  settings_dirty <- reactive({
    l <- live_cfg(); sv <- saved_cfg()
    any(vapply(display_keys, function(k) settings_value_differs(l[[k]], sv[[k]]), logical(1)))
  })

  output$display_settings_ui <- renderUI({
    settings_ui_gen()
    ui <- .pt_display_settings_ui(isolate(saved_cfg()), "s_")
    tagList(ui$display, .pt_settings_group("Comparison (diff) tracks", ui$diff))
  })
  # Rendered while the panel is still collapsed, so every s_ input reports
  # its value from the start.
  outputOptions(output, "display_settings_ui", suspendWhenHidden = FALSE)

  output$settings_dirty_ui <- renderUI({
    if (settings_dirty()) tags$span(class = "text-warning", "Unsaved changes (this session only).")
  })

  # Save as project default: writes only the non-core keys into the saved
  # project_settings.json (re-read, so core keys and samples stay exactly as
  # saved), and refreshes project_samples.txt the way the wizard's Save does.
  observeEvent(input$save_settings, {
    live <- live_cfg()
    ok <- tryCatch({
      # Never write a project from defaults (e.g. reset from another tab meanwhile).
      if (!file.exists(SETTINGS_JSON)) stop(SETTINGS_JSON, " no longer exists; restart the app.", call. = FALSE)
      x <- .pt_read_settings(SETTINGS_JSON, .pt_default_settings())
      x[display_keys] <- live[display_keys]
      .pt_write_project_files(x, SETTINGS_JSON)
      TRUE
    }, error = function(e) {
      showNotification(paste("Couldn't save the settings:", conditionMessage(e)), type = "error", duration = NULL)
      FALSE
    })
    if (ok) {
      saved_cfg(live)
      showNotification("Saved as the project default.", type = "message")
    }
  })

  # Revert to saved: re-render the panel from saved_cfg(); the browser then
  # reports the re-rendered inputs' values.
  observeEvent(input$revert_settings, settings_ui_gen(settings_ui_gen() + 1))

  # Shown in the Edit project settings / Clear caches / Close app confirms.
  unsaved_settings_note <- function() {
    if (isTRUE(settings_dirty()))
      p(class = "text-warning", strong("Your Figure settings changes are not saved and will be lost."),
        " Use Save as project default first to keep them.")
  }

  session$onSessionEnded(function() {
    isolate({
      if (!is.null(state$svg_path) && file.exists(state$svg_path))
        unlink(state$svg_path)
    })
  })

  # ── Settings (re-enter wizard) ────────────────────────────────────────────────
  # Editing settings stops this app and restarts into the wizard, a separate app
  # (run_plottrellis()'s loop), since everything here is built from the settings
  # the wizard edits -- confirm first so this isn't a silent surprise
  # mid-session, rather than restarting immediately on click.
  observeEvent(input$open_settings, {
    showModal(modalDialog(
      title = "Edit project settings?", easyClose = FALSE,
      p("This restarts the app into the setup wizard, pre-filled with your current settings. Any in-progress figure/selection in this session will be lost."),
      unsaved_settings_note(),
      footer = tagList(modalButton("Cancel"), actionButton("confirm_edit_settings", "Continue", class = "btn-primary"))
    ))
  })

  observeEvent(input$confirm_edit_settings, {
    # Same "show feedback, then delay stopApp()" pattern as
    # confirm_reset_project below -- an immediate stopApp() can kill the
    # process before the modal's removal/replacement actually
    # reaches the browser.
    showModal(modalDialog(
      title = "Restarting into the setup wizard", easyClose = FALSE, footer = NULL,
      p(strong("You can close this tab now."), " The setup wizard opens in a new tab from RStudio in a few seconds. ",
        "If nothing opens, go back to RStudio: the app's address is in the Console.")
    ))
    later::later(function() stopApp(returnValue = "edit_settings"), 0.2)
  }, ignoreInit = TRUE)

  observeEvent(input$reset_project_btn, {
    showModal(modalDialog(
      title = "Reset project configuration?", easyClose = FALSE,
      p("This deletes project_settings.json and restarts the app into first-time setup, which offers to restore it from project_samples.txt (kept). Your own samples.txt is not touched."),
      footer = tagList(modalButton("Cancel"), actionButton("confirm_reset_project", "Delete and restart", class = "btn-danger"))
    ))
  })

  observeEvent(input$confirm_reset_project, {
    if (file.exists(SETTINGS_JSON)) file.remove(SETTINGS_JSON)
    showModal(modalDialog(
      title = "Project configuration reset", easyClose = FALSE, footer = NULL,
      p("project_settings.json has been deleted."),
      p("Restarting into first-time setup."),
      p(strong("You can close this tab now."), " The app reopens in a new tab from RStudio in a few seconds. ",
        "If nothing opens, go back to RStudio: the app's address is in the Console.")
    ))
    later::later(function() stopApp(returnValue = "restart"), 0.2)
  }, ignoreInit = TRUE)

  # ── Close app ────────────────────────────────────────────────────────────────
  # Same confirm-then-delayed-stopApp() pattern as reset/edit-settings, but
  # with a returnValue that ISN'T "restart"/"edit_settings", so
  # run_plottrellis()'s loop returns instead of relaunching: this actually
  # ends the app, not just the tab. Manual only -- no auto-stop on tab close/
  # disconnect, since a flaky connection dropping and reconnecting would
  # otherwise kill the app out from under the user.
  observeEvent(input$close_app_btn, {
    showModal(modalDialog(
      title = "Close app?", easyClose = FALSE,
      p("This stops the R process running this app. You'll need to relaunch it to use PlotTrellis again."),
      unsaved_settings_note(),
      footer = tagList(modalButton("Cancel"), actionButton("confirm_close_app", "Close app", class = "btn-primary"))
    ))
  })

  observeEvent(input$confirm_close_app, {
    showModal(modalDialog(
      title = "App closed", easyClose = FALSE, footer = NULL,
      p("You can close this browser tab now.")
    ))
    later::later(function() stopApp(returnValue = "closed"), 0.2)
  }, ignoreInit = TRUE)

  # ── Clear all caches ─────────────────────────────────────────────────────────
  # The disk caches the dialog lists (the chromosome-sizes file is kept);
  # project_settings.json is left alone, so the restart goes straight back into the
  # composer, and each cache rebuilds on next use.
  observeEvent(input$clear_caches_btn, {
    showModal(modalDialog(
      title = "Clear all caches?", easyClose = FALSE,
      p("This deletes every cached derived file this app writes:"),
      tags$ul(
        tags$li("Genome TxDb / transcript-name / transcript-meta lookups"),
        tags$li("Named-transcript display cache"),
        tags$li("UniProt canonical-isoform cache"),
        tags$li("Protein-domain annotation cache (all genes)")
      ),
      p("project_settings.json and samples.txt are not touched. Everything above rebuilds automatically (one-time cost) the next time it's needed."),
      unsaved_settings_note(),
      footer = tagList(modalButton("Cancel"), actionButton("confirm_clear_caches", "Delete and restart", class = "btn-warning"))
    ))
  })

  observeEvent(input$confirm_clear_caches, {
    cache_root <- .pt_cache_root(RENDER_CFG$cache_dir)
    cache_files <- c(GENOME$txdb_path, GENOME$tx_name_lookup_path, GENOME$tx_meta_lookup_path,
                      file.path(cache_root, "gene_display_transcripts_cache.rds"),
                      file.path(cache_root, "uniprot_canonical_cache.rds"))
    for (f in cache_files) if (!is.null(f) && file.exists(f)) file.remove(f)
    pd_dir <- file.path(cache_root, "protein_domains")
    if (dir.exists(pd_dir)) unlink(pd_dir, recursive = TRUE)

    showModal(modalDialog(
      title = "Caches cleared", easyClose = FALSE, footer = NULL,
      p("All caches have been deleted."),
      p("Restarting with your current project settings."),
      p(strong("You can close this tab now."), " The app reopens in a new tab from RStudio in a few seconds. ",
        "If nothing opens, go back to RStudio: the app's address is in the Console.")
    ))
    later::later(function() stopApp(returnValue = "restart"), 0.2)
  }, ignoreInit = TRUE)

  # Renders read live_cfg() (Figure settings, above); core lookups read
  # RENDER_CFG directly.

  # ── Region: chromosome selector ─────────────────────────────────────────────
  output$chr_select_ui <- renderUI({
    chrs <- sort_chrs(names(CHR_LENGTHS))
    selectizeInput("region_chr", "Chromosome:", choices = chrs, selected = chrs[1])
  })

  chr_len_r <- reactive({
    req(input$region_chr)
    CHR_LENGTHS[[input$region_chr]] %||% 200000000L
  })

  output$range_slider_ui <- renderUI({
    len <- chr_len_r()
    sliderInput("region_range", "Region:",
                min = 1L, max = len, value = c(1L, len),
                step = 1000L, width = "100%", sep = ",")
  })

  output$start_num_ui <- renderUI({
    textInput("region_start_text", "Start:",
              value = format_coord(1L), width = "100%",
              placeholder = "e.g. 17Mb, 15.5kb, 17000000")
  })

  output$end_num_ui <- renderUI({
    textInput("region_end_text", "End:",
              value = format_coord(chr_len_r()), width = "100%",
              placeholder = "e.g. 17.49Mb, 17490kb")
  })

  # Slider → text: update display when slider moves
  observeEvent(input$region_range, {
    req(input$region_range)
    updateTextInput(session, "region_start_text", value = format_coord(input$region_range[1]))
    updateTextInput(session, "region_end_text",   value = format_coord(input$region_range[2]))
  }, ignoreInit = TRUE)

  # Text → slider: parse on blur/Enter (JS suppresses auto-fire; only explicit events reach here)
  observeEvent(input$region_start_text, {
    txt <- input$region_start_text
    if (is.null(txt) || nchar(trimws(txt)) == 0) return()
    req(input$region_range)
    v <- parse_bp(txt)
    if (is.na(v)) {
      state$start_parse_err <- "Cannot parse start. Accepts: 17000000, 17Mb, 15.5kb (case-insensitive)"
      return()
    }
    state$start_parse_err <- NULL
    if (!isTRUE(input$region_range[1] == v))
      updateSliderInput(session, "region_range", value = c(v, input$region_range[2]))
  }, ignoreInit = TRUE, ignoreNULL = TRUE)

  observeEvent(input$region_end_text, {
    txt <- input$region_end_text
    if (is.null(txt) || nchar(trimws(txt)) == 0) return()
    req(input$region_range)
    v <- parse_bp(txt)
    if (is.na(v)) {
      state$end_parse_err <- "Cannot parse end. Accepts: 17000000, 17Mb, 15.5kb (case-insensitive)"
      return()
    }
    state$end_parse_err <- NULL
    if (!isTRUE(input$region_range[2] == v))
      updateSliderInput(session, "region_range", value = c(input$region_range[1], v))
  }, ignoreInit = TRUE, ignoreNULL = TRUE)

  region_valid <- reactive({
    req(input$region_chr, input$region_range)
    len   <- CHR_LENGTHS[[input$region_chr]] %||% Inf
    start <- as.integer(input$region_range[1])
    end   <- as.integer(input$region_range[2])
    if (is.na(start) || start < 1)  return(list(valid = FALSE, msg = "Start must be \u2265 1"))
    if (is.na(end) || end <= start) return(list(valid = FALSE, msg = "End must be > start"))
    if (end > len)                  return(list(valid = FALSE, msg = "End exceeds chromosome length"))
    list(valid = TRUE, chrom = input$region_chr, start = start, end = end)
  })

  output$region_error_ui <- renderUI({
    if (is.null(input$region_mode) || input$region_mode != "coords") return(NULL)
    parse_errs <- c(state$start_parse_err, state$end_parse_err)
    if (length(parse_errs) > 0)
      return(div(class = "alert alert-warning", paste(parse_errs, collapse = " | ")))
    v <- tryCatch(region_valid(), error = function(e) NULL)
    if (!is.null(v) && !isTRUE(v$valid))
      div(class = "alert alert-danger", v$msg)
  })

  # ── Region: gene lookup (shared by all 3 modes) ─────────────────────────────
  # The symbol this observer itself just put in the box, so the re-run that
  # update triggers keeps a typed Ensembl ID's own region rather than
  # re-resolving the symbol to whichever of its IDs comes first.
  gene_echo <- NULL
  # The text the current gene was resolved from, so a flank-only re-run
  # (Figure settings) re-resolves a typed Ensembl ID, not its echoed symbol.
  gene_query <- NULL
  observe({
    sym <- trimws(input$gene_name %||% "")
    # Read first, so the echo run below still depends on it.
    pct <- (input$s_gene_flank_pct %||% saved_cfg()$gene_flank_pct) / 100
    if (identical(sym, gene_echo)) { gene_echo <<- NULL; return() }
    gene_echo <<- NULL
    if (nchar(sym) == 0) { state$gene_symbol <- NULL; gene_query <<- NULL; return() }
    q <- if (!is.null(gene_query) && identical(sym, isolate(state$gene_symbol))) gene_query else sym
    gene_query <<- q
    res <- resolve_gene_region(q, txdb, pct)
    state$gene_region <- res$region
    state$gene_error  <- res$error
    state$gene_symbol <- res$symbol
    state$gene_note   <- res$note
    # An Ensembl ID or wrong-case symbol resolved: show the canonical symbol,
    # which is what every track and panel below looks genes up by.
    if (!is.null(res$region) && !identical(res$symbol, sym)) {
      gene_echo <<- res$symbol
      updateSelectizeInput(session, "gene_name", choices = ALL_GENES, server = TRUE, selected = res$symbol)
    }
  })

  observeEvent(input$region_mode, { state$gene_error <- NULL })

  output$gene_result_ui <- renderUI({
    if (!is.null(state$gene_error))
      div(class = "alert alert-danger", state$gene_error)
    else if (!is.null(state$gene_region)) {
      r <- state$gene_region
      div(class = "alert alert-success",
          strong("Found: "),
          paste0(r$chrom, ":", format_coord(r$start), " \u2013 ", format_coord(r$end),
                 " (", format_bp(r$end - r$start), ")"),
          if (!is.null(state$gene_note)) tagList(br(), state$gene_note))
    }
  })

  # Gene / protein info: Output page, about the rendered figure's gene (R/gene_info_panel.R).
  .pt_gene_info_server(output, session, state, txdb, CDS_BY_TX,
    genome = RENDER_CFG[c("orgdb_pkg", "tx_name_lookup_path", "tx_meta_lookup_path", "genome_label", "cache_dir")],
    reviewed_only = reactive(input$s_uniprot_reviewed_only %||% saved_cfg()$uniprot_reviewed_only))

  # ── Two-sample reference/comparison mutual exclusion ──
  # Both boxes always offer every condition; the only guard is against picking the same one in
  # both -- sel_prev tracks each box's last-settled value so we can tell which box just changed
  # (vs. one only just relabeled by this same observer) and clear the OTHER box on a collision.
  # ignoreInit = TRUE since the initial choices/markers are already baked into the UI; both
  # start "" (neither box preselects a condition), matching sel_reference/sel_comparison's
  # own defaults.
  sel_prev <- reactiveValues(control = "", treatment = "")

  observeEvent(list(input$sel_reference, input$sel_comparison), {
    ctrl <- input$sel_reference    %||% ""
    trt  <- input$sel_comparison  %||% ""

    if (!identical(ctrl, sel_prev$control) && nchar(ctrl) > 0 && identical(ctrl, trt)) {
      trt <- ""  # control just changed to collide with treatment -- clear treatment
    } else if (!identical(trt, sel_prev$treatment) && nchar(trt) > 0 && identical(trt, ctrl)) {
      ctrl <- ""  # treatment just changed to collide with control -- clear control
    }

    sel_prev$control   <- ctrl
    sel_prev$treatment <- trt

    updateSelectizeInput(session, "sel_reference",
      choices = condition_choices(), selected = ctrl)
    updateSelectizeInput(session, "sel_comparison",
      choices = condition_choices(), selected = trt)
  }, ignoreInit = TRUE)

  # ── Which condition(s) are currently selected, mode-dependent ──
  selected_conditions <- reactive({
    if (identical(input$view_mode, "two")) c(input$sel_reference, input$sel_comparison)
    else input$sel_condition
  })

  # Whether the condition-selection stage is fully resolved for the current
  # mode: Compose has no condition dropdown at all (the track picker is its
  # own sample selector), so it's trivially ready; single/two need their
  # dropdown(s) actually filled in, since none of sel_condition/
  # sel_reference/sel_comparison preselect a condition. Shared by the region/gene lock, the track-pool
  # populator, and Generate below so all three agree on the same stage
  # boundary.
  has_conditions_ready <- reactive({
    if (identical(input$view_mode, "compose")) TRUE
    else if (identical(input$view_mode, "two"))
      nchar(input$sel_reference %||% "") > 0 && nchar(input$sel_comparison %||% "") > 0
    else if (identical(input$view_mode, "single")) nchar(input$sel_condition %||% "") > 0
    else FALSE
  })

  # Re-scopes and re-defaults the track picker when the mode or the selected
  # condition(s) change, never on a region/gene change (that would wipe what the
  # user has placed; the picker is locked by region_lock instead). figure_pool_key
  # is the (mode, conditions) the pool was last built from: a trigger that leaves
  # it unchanged is skipped. Runs at startup too (no ignoreInit): the picker's
  # static markup holds Multi's pool, which the other modes must replace.
  figure_pool_key <- reactiveVal(NULL)

  observeEvent(list(input$view_mode, input$sel_condition, input$sel_reference, input$sel_comparison), {
    if (identical(input$view_mode, "compose")) {
      key   <- "compose"
      split <- list(
        available = track_labels(COMPOSE_AVAILABLE_KEYS),
        figure    = list(genes = HTML(format_track_label(TRACK_ITEMS[["genes"]]$label)))
      )
      if (identical(START_VIEW$view_mode, "compose"))
        split <- list(available = track_labels(setdiff(COMPOSE_AVAILABLE_KEYS, setdiff(START_VIEW$tracks, "genes"))),
                      figure    = track_labels(START_VIEW$tracks))
    } else {
      conds <- selected_conditions()
      if (length(conds) == 0 || anyNA(conds) || any(nchar(conds) == 0)) {
        key   <- "not_ready"
        split <- list(available = list(), figure = list())
      } else {
        key   <- paste(input$view_mode, paste(conds, collapse = "|"), sep = "::")
        split <- track_split(conds)
      }
    }
    if (identical(key, figure_pool_key())) return()
    figure_pool_key(key)
    update_rank_list("available_tracks", labels = split$available, session = session)
    update_rank_list("figure_tracks",    labels = split$figure,    session = session)
  })

  # ── Region resolution: identical across all 3 sample-count modes ───────────
  # Coordinates is selectable in every mode, so which region resolves depends
  # only on region_mode -- never on view_mode. gene_name is only meaningful in
  # 'gene' region_mode; tracks that need a gene symbol (GENE_REQUIRED_TYPES)
  # can't draw from coordinates in any mode, which the picker's red-boxing
  # handles, not this resolution logic.
  resolved_region <- reactive({
    if (!is.null(input$region_mode) && input$region_mode == "coords") {
      v <- tryCatch(region_valid(), error = function(e) NULL)
      if (!is.null(v) && isTRUE(v$valid))
        list(rgn = list(chrom = v$chrom, start = v$start, end = v$end), gene_name = NA_character_)
      else NULL
    } else if (!is.null(state$gene_region)) {
      list(rgn = state$gene_region, gene_name = state$gene_symbol %||% input$gene_name)
    } else NULL
  })

  # Whether a sample-count mode has been chosen at all (view_mode starts
  # unselected) -- Generate and the locks below all stay gated on this in
  # addition to conditions/region, so the intended stage order (mode ->
  # conditions -> region/gene -> tracks -> generate) holds whatever else is
  # already set.
  has_mode <- reactive({ !is.null(input$view_mode) })

  # region_mode/gene_name lock until mode AND conditions are both resolved,
  # so a gene can't be searched before Reference/Comparison (or Condition)
  # are chosen. Mirrors the 'region_lock' message
  # below (same server-driven pattern), consumed by setConditionGateLocked().
  observe({
    session$sendCustomMessage("condition_lock", list(locked = !(has_mode() && has_conditions_ready())))
  })

  # The track picker (both Available and Figure tracks) locks whenever
  # mode+conditions+region aren't all ready -- covers the "nothing chosen
  # yet" case and the mid-composition one (switching to Coordinates with no
  # coordinates entered locks the existing Figure tracks too, since it's
  # driven by the same single flag rather than two separate lock layers).
  # The track-pool populator above also waits for the conditions.
  observe({
    session$sendCustomMessage("region_lock", list(locked = !(has_mode() && has_conditions_ready() && !is.null(resolved_region()))))
  })

  # Figure-tracks keys that can't draw as things stand -- same two conditions
  # the client already red-boxes (the track-gene-required/track-broken CSS
  # classes, applyGeneRequiredFlag()/applyBrokenFlag() above), reimplemented
  # server-side so Generate is disabled while one is there.
  invalid_figure_tracks <- reactive({
    keys <- input$figure_tracks %||% character(0)
    if (length(keys) == 0) return(character(0))
    is_coords <- identical(input$region_mode, "coords")
    Filter(function(k) {
      item <- TRACK_ITEMS[[k]]
      !is.null(item) && (isTRUE(item$broken) || (is_coords && item$type %in% GENE_REQUIRED_TYPES))
    }, keys)
  })

  # ── Generate button (disabled until mode + conditions + region + tracks,
  # and no invalid track sitting in Figure tracks) ───────────────────────────
  output$generate_ui <- renderUI({
    labels     <- input$figure_tracks %||% character(0)
    has_tracks <- length(labels) > 0
    has_region <- !is.null(resolved_region())
    has_valid_tracks <- length(invalid_figure_tracks()) == 0
    generate_button_ui("generate", has_mode() && has_conditions_ready() && has_region && has_tracks && has_valid_tracks)
  })

  # Message beside Generate naming the red-boxed track(s) blocking it --
  # mirrors the .track-gene-required/.track-broken title tooltip text
  # (applyGeneRequiredFlag()/applyBrokenFlag() above) but as visible text,
  # since a disabled button with no explanation is a cryptic dead end.
  output$generate_warning_ui <- renderUI({
    bad <- invalid_figure_tracks()
    if (length(bad) == 0) return(NULL)
    names_str <- paste(vapply(bad, function(k) TRACK_ITEMS[[k]]$label %||% k, character(1)), collapse = ", ")
    div(class = "alert alert-warning", style = "margin-bottom:0;",
        strong("Remove or fix before generating: "), names_str)
  })

  # Render helper: render_all() straight to an svglite device, in-process (Generate, zoom and Apply to figure); the displayed/downloadable R script is assemble_download_script() on the same entries/region_params/settings.
  do_render <- function(region_params, entries, settings_base, title_str = "", labels = NULL) {
    state$rendering <- TRUE
    on.exit({
      state$rendering <- FALSE
      session$sendCustomMessage("hide_spinner", list())
    }, add = TRUE)
    settings <- c(settings_base, list(
      title       = title_str,
      output_path = tempfile(fileext = ".svg")
    ))
    state$plot_error <- NULL
    tryCatch({
      dev_before <- dev.cur()
      # .pt_quiet(): plotgardener prints a "text[text2]"-style line per grob, hundreds per figure.
      .pt_quiet(render_all(entries, region_params, settings,
                 device_fn = function(w, h) svglite(settings$output_path, width = w, height = h)))
      if (dev.cur() != dev_before) dev.off()
      .pt_tidy_svg(settings$output_path, unlock_text = FALSE)
      if (!is.null(state$svg_path) && file.exists(state$svg_path))
        unlink(state$svg_path)
      state$svg_path         <- settings$output_path
      state$rcode             <- assemble_download_script(entries, region_params, settings)
      state$last_region      <- region_params
      state$last_labels      <- labels
      state$last_title       <- title_str
      state$rendered_gene    <- region_params$gene  # Gene / protein info follows the figure
      session$sendCustomMessage("show_output", list())
    }, error = function(e) {
      if (dev.cur() != dev_before) try(dev.off(), silent = TRUE)
      state$plot_error <- conditionMessage(e)
      session$sendCustomMessage("show_output", list())
    })
  }

  # Shared by Generate, zoom and Apply to figure: build_entries() + settings_base + do_render(), from a region/labels/gene-name/title already resolved by the caller.
  render_from_selection <- function(rgn, labels, gene_name, title_str) {
    cfg     <- live_cfg()
    entries <- build_entries(labels, DEFAULT_SAMPLE_COLS, cfg)

    region_params <- list(
      chrom = rgn$chrom,
      start = as.integer(rgn$start),
      end   = as.integer(rgn$end),
      gene  = gene_name
    )
    settings_base <- .pt_settings_base(cfg)

    do_render(region_params, entries, settings_base, title_str, labels)
  }

  # ── Generate: build specs from UI then render ────────────────────────────────
  observeEvent(input$generate, {
    req(!isTRUE(state$rendering))
    req(has_mode())
    labels <- input$figure_tracks %||% character(0)
    req(length(labels) > 0)

    rr <- resolved_region()
    req(!is.null(rr))
    rgn       <- rr$rgn
    gene_name <- rr$gene_name

    if (identical(input$view_mode, "two")) {
      req(nchar(input$sel_reference %||% "") > 0, nchar(input$sel_comparison %||% "") > 0)
    }

    render_from_selection(rgn, labels, gene_name, trimws(input$plot_title %||% ""))
  })

  # ── Zoom: recentre and rescale the last rendered region ─────────────────────
  # Rebuilds entries and settings from live_cfg(), keeping the
  # last figure's tracks, gene and title.
  zoom_to <- function(factor) {
    req(!isTRUE(state$rendering))
    req(!is.null(state$last_region), !is.null(state$last_labels))
    r         <- state$last_region
    mid       <- (r$start + r$end) / 2
    hw        <- max((r$end - r$start) / 2 / factor, 500L)
    clen      <- CHR_LENGTHS[[r$chrom]] %||% 2e8
    new_start <- as.integer(max(1L,   round(mid - hw)))
    new_end   <- as.integer(min(clen, round(mid + hw)))
    new_rgn   <- list(chrom = r$chrom, start = new_start, end = new_end)
    updateSliderInput(session, "region_range",      value = c(new_start, new_end))
    updateTextInput(session,   "region_start_text", value = format_coord(new_start))
    updateTextInput(session,   "region_end_text",   value = format_coord(new_end))
    render_from_selection(new_rgn, state$last_labels, r$gene, state$last_title)
  }
  observeEvent(input$zoom_in_15x,  zoom_to(1.5))
  observeEvent(input$zoom_in_3x,   zoom_to(3))
  observeEvent(input$zoom_in_10x,  zoom_to(10))
  observeEvent(input$zoom_out_15x, zoom_to(1 / 1.5))
  observeEvent(input$zoom_out_3x,  zoom_to(1 / 3))
  observeEvent(input$zoom_out_10x, zoom_to(0.1))

  # Apply to figure: re-render the last figure, same region, with the live settings.
  # A generate-btn, so the spinner/double-click gate covers it; shown only once a
  # figure exists.
  output$apply_settings_ui <- renderUI({
    req(!is.null(state$last_region))
    actionButton("apply_settings", "Apply to figure", class = "btn-primary generate-btn")
  })
  observeEvent(input$apply_settings, zoom_to(1))

  # ── Output page ──────────────────────────────────────────────────────────────
  output$zoom_ui <- renderUI({
    req(!is.null(state$last_region))
    tags$div(
      style = "display:flex; justify-content:center; align-items:flex-start; gap:28px; padding:8px 0 12px;",
      tags$div(
        style = "text-align:center;",
        tags$div(class = "zoom-group-label",
                 "Zoom in"),
        actionButton("zoom_in_15x", "1.5\u00D7", class = "btn-default"),
        actionButton("zoom_in_3x",  "3\u00D7",   class = "btn-default"),
        actionButton("zoom_in_10x", "10\u00D7",  class = "btn-default")
      ),
      tags$div(
        style = "width:1px; background:#ddd; align-self:stretch; margin-top:22px;"
      ),
      tags$div(
        style = "text-align:center;",
        tags$div(class = "zoom-group-label",
                 "Zoom out"),
        actionButton("zoom_out_15x", "1.5\u00D7", class = "btn-default"),
        actionButton("zoom_out_3x",  "3\u00D7",   class = "btn-default"),
        actionButton("zoom_out_10x", "10\u00D7",  class = "btn-default")
      )
    )
  })

  output$render_error_ui <- renderUI({
    req(state$plot_error)
    div(class = "alert alert-danger",
        strong("Error: "), state$plot_error,
        br(),
        tags$pre(style="font-size:11px;margin-top:8px;", state$plot_error))
  })

  output$svg_display_ui <- renderUI({
    req(state$svg_path, file.exists(state$svg_path))
    inline_svg(state$svg_path)
  })

  output$rcode_output <- renderText({
    state$rcode %||% ""
  })

  output$dl_svg <- downloadHandler(
    filename = function() paste0("plottrellis_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".svg"),
    content  = function(file) {
      req(state$svg_path, file.exists(state$svg_path))
      file.copy(state$svg_path, file)
      .pt_tidy_svg(file)
    }
  )

  output$dl_rcode <- downloadHandler(
    filename = function() paste0("plottrellis_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".R"),
    content  = function(file) {
      req(state$rcode)
      writeLines(state$rcode, file)
    }
  )

}

shinyApp(ui = ui, server = server)

}
