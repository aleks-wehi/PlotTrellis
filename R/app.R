# PlotTrellis -- Multi-modal genomic visualisation
# 2-page app: Compose | Output  (project settings in session-only modal)
# Rendering pipeline: 4-script assembly → source() → SVG
# Single-script principle: downloaded R script is exactly what was sourced.
#
# Package entry point. Everything below used to be a Shiny app.R sourced
# directly by Shiny (which lands ui/server-supporting code in separate
# envs, hence the old file's APP_DIR/SCRIPTS_DIR/SETTINGS_JSON/
# DEFAULT_SETTINGS/CFG assign(envir = globalenv()) calls and its many
# source(file.path(APP_DIR, ...)) lines) -- now that it's one function body,
# ui/server and every helper below share this function's own lexical scope
# by ordinary closure, and every former source()d file is just another
# function in this package's namespace, loaded once at package-load time.
# `project_dir` is where this project's project_settings.json (and its
# downloaded-script/cache state) lives -- defaults to a per-user config dir
# so multiple projects can each pass their own.
#' Launch the PlotTrellis Shiny app
#'
#' Starts the PlotTrellis genomic track composer: a 2-page Shiny app
#' (Compose | Output) that lets you build plotgardener-based visualisations
#' without writing code, and download a standalone R script that reproduces
#' the exact figure offline. Blocks until the app is closed; on first run (or
#' when re-entering project settings), automatically launches the setup
#' wizard first and continues straight into the composer once it's saved --
#' no separate re-launch needed.
#'
#' @param project_dir Directory holding this project's `project_settings.json`
#'   and its cache/downloaded-script state. Defaults to a per-user config
#'   directory (via [tools::R_user_dir()]) so multiple projects can each pass
#'   their own.
#' @param ... Passed through to [shiny::runApp()] (e.g. `launch.browser`,
#'   `port`, `host`).
#' @return Invisibly, the return value of the final [shiny::runApp()] call
#'   (e.g. once the browser tab is closed) -- not a `shiny.appobj`, since this
#'   function itself runs the app rather than returning one to run.
#' @examples
#' if (interactive()) {
#'   PlotTrellisApp(project_dir = tempdir())
#' }
#' @export
PlotTrellisApp <- function(project_dir = tools::R_user_dir("PlotTrellis", which = "config"), ...) {

# tools::R_user_dir() only returns the path, it doesn't create it -- on a
# genuinely first run (nothing has ever written here before) the wizard's
# Save handler's writeLines(SETTINGS_JSON) fails with "No such file or
# directory" otherwise.
dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)
SETTINGS_JSON <- file.path(project_dir, "project_settings.json")

DEFAULT_SETTINGS <- list(
  # rationale: see solutions/concepts/19-relocated-code-rationale.md, "app.R:54"
  gtf_path            = "",
  orgdb               = GENOME_PRESETS$mm39$orgdb,
  genome_label        = GENOME_PRESETS$mm39$genome_label,
  # Same "portable blank default, wizard requires a real path before Save"
  # pattern as gtf_path above -- the wizard's "Data root directory" field
  # (utils/wizard.R) already defaults to "" for a new project regardless of
  # this value, so a hardcoded deployment-specific path here was dead for
  # the normal flow and only a latent portability trap via CFG <-
  # modifyList(DEFAULT_SETTINGS, parsed) if a saved settings file were ever
  # missing this field.
  data_dir            = "",
  # Phase C.3: samples live embedded here, not behind a samples_txt file
  # pointer -- samples.txt is now a pure import/export interchange format
  # (utils/samples_txt_io.R), never read at startup. Each entry:
  # list(condition=, colour=, data=list(list(type=,name=,path=))) -- same
  # shape import_samples_txt()/export_samples_txt() already use. `condition`
  # is the sole identifying field (no sample_id/replicate); a replicate is
  # just a separate entry with its own distinct condition string.
  samples     = list(),
  font        = "DejaVu Sans",
  fontsize    = 8L,
  page_w      = 9,
  chr_fmt     = "dot",
  # Per-track value axis (feature: toggleable axes), drawn on bw/bw_diff_stacked/
  # sashimi tracks -- the ones with no other value scale already shown (unlike
  # hic's annoHeatmapLegend or tpm_usage's own always-on hand-drawn axis, which
  # this setting leaves untouched). On by default: reading real values off a
  # track without asking a bioinformatician is the whole point of this app.
  show_track_axes = TRUE,
  hic_norm    = "KR",
  hic_style   = "rectangle",
  zrange_cap  = NULL,   # NULL = auto (colour scale runs to the observed max, see viz/hic.R)
  atac_height     = 0.85,
  hic_height      = NULL,   # NULL = auto (half track width for rectangle, full for square)
  hic_res_mode    = "auto", # "auto" | "finest"
  hic_max_bins_per_view = 500, # bin-count ceiling for Hi-C auto-resolution selection
  gene_flank_pct  = 15,     # % of gene length added to each side of gene view
  # Used by bw_diff_stacked's gain/loss colours (VIZ_CONVENTIONS.md §6)
  diff_colour_lo  = "#0033CC",  # vivid blue, highly downregulated
  diff_colour_hi  = "#CC0000",  # vivid red,  highly upregulated
  # Junctions mirrored diff track (BUGS.md #4): hard ceiling on the opacity
  # scale (log2FC that maps to full opacity), so one small/sparse arc's
  # outsized ratio can't stretch the scale and wash out every other arc's
  # opacity.
  junction_opacity_cap = 2,
  # Minimum pooled nSupport for a paired-end arc (BED12 splice junction or
  # BEDPE loop call) to be drawn. Default 0 (show everything) -- unlike a
  # real read count, a BEDPE score column's scale/meaning is caller-specific
  # (e.g. FitHiChIP writes a continuous significance score, not a read
  # count), so there's no safe nonzero default. See utils/transcript_filter.R.
  min_junction_support = 0,
  # Per-condition mean-TPM floor a transcript must clear to be individually
  # named/displayed (vs. folded into "Other"). Default 1 -- see
  # utils/transcript_filter.R's comment for why this is a plain absolute
  # expression floor rather than a statistical test.
  min_display_tpm = 1,
  # Fixed ROLE colour for every 2-sample comparison track (VIZ_CONVENTIONS.md §3.2)
  # -- Reference/Comparison always render in these two colours regardless of
  # either condition's own identity/PALETTE colour, so "which one is the
  # baseline" reads at a glance instead of depending on an arbitrary sample hue.
  # THESE TWO DEFAULTS ARE DUPLICATED (by necessity -- same pattern as
  # diff_colour_lo/diff_colour_hi below) in utils/wizard.R's colour <input>
  # `value=` and in RENDER_CFG's fallback chain -- if you change either
  # hex here, change it in those two places too. Full rationale + the exact
  # checklist of every spot these values live: VIZ_CONVENTIONS.md §3.2.
  reference_colour  = "#4D4D4D", # neutral grey, always the Reference condition
  comparison_colour = "#56B4E9", # colour-blind-safe blue, always the Comparison condition
  uniprot_reviewed_only = FALSE, # UniProt link-out: Swiss-Prot only if TRUE
  # Junction/arc tracks (feature 25, solutions/concepts/18-bed12-bedpe-arc-format.md):
  # whether an arc with only ONE anchor inside the plotted region is
  # truncated at the plot boundary (TRUE, plotPairsArches(clip=)) or drawn
  # extending past it (FALSE, matches the shipped sashimi look unchanged).
  # Off by default -- BEDPE-sourced long-range pairs (e.g. Hi-C loops) are
  # the case expected to actually need this; splice junctions rarely span
  # far enough past a gene's view for it to matter.
  arc_clip_partial = FALSE
)

read_settings <- function() {
  if (!file.exists(SETTINGS_JSON)) return(DEFAULT_SETTINGS)
  # simplifyDataFrame = FALSE: jsonlite's default simplification turns
  # `samples[]` (an array of same-shaped objects) into a data.frame instead
  # of a list of lists, which crashes load_samples_from_settings()'s
  # `s$condition`/`s$data` access (`$` on an atomic column vector). Scalars
  # are unaffected -- only array-of-objects simplification is disabled.
  parsed <- jsonlite::fromJSON(SETTINGS_JSON, simplifyDataFrame = FALSE)
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

# Loops internally so the wizard-to-app transition needs no manual
# re-launch: CFG/WIZARD_EDIT_MODE are re-read each pass since a completed
# wizard Save (or the main app's "Edit project settings" restart) changes
# both on disk/env between iterations. Mirrors the restart contract
# wizard.R already relies on (stopApp(returnValue = "restart"/"edit_settings")),
# now consumed here directly instead of only by the pre-package deployment's
# own external loop.
repeat {

CFG <- read_settings()

# Phase C.1 step 1: no project_settings.json yet -> guided first-run wizard
# instead of the normal app. Phase C.3 decision 1 extends this branch to also
# cover re-entering the wizard to EDIT an existing settings file: the main
# app's "Edit project settings" button sets this env var and restarts --
# same process-restart requirement as first-run, see utils/wizard.R's header
# comment for why. wizard_server() itself decides prefill-vs-blank-landing
# purely from file.exists(SETTINGS_JSON), not this flag -- this flag only
# controls which of the two cases routes here at all.
WIZARD_EDIT_MODE <- identical(Sys.getenv("PLOTTRELLIS_EDIT_SETTINGS", "0"), "1")
app_obj <- if (!file.exists(SETTINGS_JSON) || WIZARD_EDIT_MODE) {

.set_wizard_globals(SETTINGS_JSON, CFG, DEFAULT_SETTINGS)
shinyApp(ui = wizard_ui(), server = wizard_server)

} else {


# DEFAULT_SETTINGS/read_settings()/CFG now defined unconditionally, above the
# wizard-vs-main-app branch (Phase C.3 step D) -- see that block's own
# comment for why.

# ── Startup: load heavy objects once ──────────────────────────────────────────
# Active genome config (standalone genome-config task, 2026-08-18): flat
# fields on CFG, no more genome_id/GENOME_TABLE lookup. One-time migration
# for a pre-existing project_settings.json still carrying the old genome_id
# field with no orgdb (this deployment's own project, before this change) --
# fill in the missing fields from GENOME_PRESETS by that id and persist the
# migrated flat fields back to disk, so it only ever runs once. A genome_id
# that isn't a known preset can't be auto-migrated (loud error, matching the
# "not silently defaulted" convention) -- re-open Settings and fill in the
# Genome section by hand.
if (is.null(CFG$orgdb) && !is.null(CFG$genome_id)) {
  preset <- GENOME_PRESETS[[CFG$genome_id]]
  if (is.null(preset))
    stop("project_settings.json has an old-format genome_id ('", CFG$genome_id,
         "') that isn't a known built-in preset, so it can't be auto-migrated. ",
         "Re-open Settings and fill in the Genome section (GTF path, OrgDb, genome label) by hand.")
  CFG$genome_id    <- NULL
  CFG$orgdb        <- preset$orgdb
  CFG$genome_label <- preset$genome_label
  writeLines(jsonlite::toJSON(CFG, auto_unbox = TRUE, pretty = TRUE, null = "null"), SETTINGS_JSON)
  message("Migrated project_settings.json's old genome_id ('", preset$label, "') to flat genome fields.")
}

# Cache paths are NEVER read from/written to CFG (2026-08-25, see utils/genomes.R's
# header) -- every launch re-derives ORGDB_PKG/TX_NAME_LOOKUP_PATH/TX_META_LOOKUP_PATH
# fresh from gtf_path; this file's own ALL_GENES lookup, utils/gene_lookup.R, and
# utils/transcript_annotation.R read these globals directly, while
# utils/transcript_meta.R/utils/transcript_filter.R (never sourced) get the
# same-named globals from scripts/01_imports.R instead, defined once per assembled
# script exactly like ORGDB_PKG below.
GENOME <- ensure_genome_ready(gtf_path = CFG$gtf_path, orgdb = CFG$orgdb, genome_label = CFG$genome_label)
if (GENOME$rebuilt) message("Built and cached this genome's TxDb/lookup files (next to its GTF).")
ORGDB_PKG           <- GENOME$orgdb
TX_NAME_LOOKUP_PATH <- GENOME$tx_name_lookup_path
TX_META_LOOKUP_PATH <- GENOME$tx_meta_lookup_path
library(ORGDB_PKG, character.only = TRUE)

message("Loading TxDb...")
txdb     <- AnnotationDbi::loadDb(GENOME$txdb_path)
ASSEMBLY <- assembly(
  Genome = GENOME$genome_label, TxDb = txdb, OrgDb = ORGDB_PKG,
  gene.id.column = "ENSEMBL", display.column = "SYMBOL"
)
message("TxDb loaded.")

# Built once here (not per gene render) since GenomicFeatures::cdsBy() walks
# the whole TxDb regardless of which genes are actually asked for -- only
# needed for the Gene info table's UniProt canonical-isoform check
# (utils/uniprot_canonical.R), which looks up a handful of transcripts'
# translated CDS lengths per gene out of this.
message("Building CDS-by-transcript index...")
CDS_BY_TX <- GenomicFeatures::cdsBy(txdb, by = "tx", use.names = TRUE)
message("CDS index built.")

# Chr lengths -- a GTF carries no seqlength info of its own (confirmed empty
# on this app's own GENCODE fixtures; see utils/build_chrom_sizes.R), so this
# is never read from the TxDb's own seqinfo. Always sourced from UCSC via
# GENOME$chrom_sizes_path (built/cached by ensure_genome_ready() above),
# never a hardcoded constant (Phase F.1.3 step 7 / Phase F.1 Q6: a silent
# fallback to one species' lengths would give wrong lengths under any other
# genome with no error).
CHR_LENGTHS <- readRDS(GENOME$chrom_sizes_path)

SAMPLES_DF  <- load_samples_from_settings(CFG$samples, CFG$data_dir)

# Every sample is selectable in every view_mode ("One"/"Two"/"Multi") -- no
# RNA-seq-eligibility pre-filter (Phase C.2 step 7, dropping the old
# GR_ELIGIBLE gate per solutions/phases/phase-C.2.md). Per-condition track
# availability is each viz/*.R module's own detect()/detect_all(), and
# per-pair availability is each viz_comparisons/*.R module's own
# detect_pair()/detect_pair_all() -- both already degrade gracefully to
# "no item" rather than erroring, so nothing here needs to pre-check RNA-seq
# (or any other data type) presence before offering a condition.
AVAILABLE_CONDITIONS <- vapply(SAMPLES_DF, function(r) r$display, character(1))

# Every sample, not just RNA-seq-eligible ones (Phase C.2 step 7 fix): each
# pair-type's detect_pair_generic() already checks per-instance availability
# on both sides of the pair, so a pre-filter here only meant e.g. two
# ATAC-only samples could never be offered a bw_diff_stacked ATAC comparison
# even though nothing else prevented it.
TRACK_ITEMS <- c(
  build_track_items(SAMPLES_DF),
  build_comparison_track_items(SAMPLES_DF)
)

message("Loading gene symbols...")
ALL_GENES <- sort(AnnotationDbi::keys(get(ORGDB_PKG, envir = asNamespace(ORGDB_PKG)), keytype = "SYMBOL"))
message("Ready.")

# ── Helpers ────────────────────────────────────────────────────────────────────
PALETTE <- c("#1f78b4","#33a02c","#e31a1c","#ff7f00","#6a3d9a",
             "#a6761d","#e6ab02","#66a61e","#e7298a","#7570b3")

# %||% is defined once, in boilerplate.R (NA-aware) -- see phase-G.1.md's
# log for why the two other duplicate definitions that used to live here and
# in samples_txt_io.R were removed rather than left to collide on load order.

# Treat NULL/NA/non-finite/list alike as "unset" (list: jsonlite::toJSON() used to round-trip a NULL setting as {}, fixed at the write site in wizard.R's wiz_save, but a settings.json written before that fix still reads {} back as an empty list via fromJSON()); shared by RENDER_CFG's zrange_cap and hic_height, both falling back to auto.
clean_optional_numeric <- function(v) {
  if (is.null(v) || is.list(v) || is.na(v) || !is.finite(v)) NULL else v
}

# Every render/lookup call site's single "cfg" object (Phase C.3 decision 7,
# replaces the old cur_settings() reactive + settings_modal()'s session-only
# input$s_... overrides). CFG is the sole source of truth again -- all
# editing now goes through the wizard's restart round trip (utils/wizard.R),
# so this only needs computing once, not per-render. Bundles a few
# genome-derived globals (orgdb_pkg/genome_label/tx_*_lookup_path) alongside
# the flat CFG keys purely so callers keep getting one object with both.
RENDER_CFG <- list(
  txdb_path   = GENOME$txdb_path,
  gtf_path    = CFG$gtf_path,
  samples     = CFG$samples,
  orgdb_pkg           = ORGDB_PKG,
  genome_label        = GENOME$genome_label,
  tx_name_lookup_path = TX_NAME_LOOKUP_PATH,
  tx_meta_lookup_path = TX_META_LOOKUP_PATH,
  data_dir    = CFG$data_dir,
  font        = trimws(CFG$font %||% "DejaVu Sans"),
  fontsize    = CFG$fontsize %||% 8L,
  page_w      = CFG$page_w %||% 9,
  chr_fmt     = CFG$chr_fmt %||% "dot",
  show_track_axes = CFG$show_track_axes %||% TRUE,
  hic_norm    = trimws(CFG$hic_norm %||% "KR"),
  hic_style   = CFG$hic_style %||% "rectangle",
  zrange_cap  = clean_optional_numeric(CFG$zrange_cap),
  atac_height  = CFG$atac_height %||% 0.85,
  hic_height   = clean_optional_numeric(CFG$hic_height),
  hic_res_mode = CFG$hic_res_mode %||% "auto",
  hic_max_bins_per_view = CFG$hic_max_bins_per_view %||% 500,
  diff_colour_lo  = trimws(CFG$diff_colour_lo %||% "#0033CC"),
  diff_colour_hi  = trimws(CFG$diff_colour_hi %||% "#CC0000"),
  junction_opacity_cap = CFG$junction_opacity_cap %||% 2,
  min_junction_support = CFG$min_junction_support %||% 0,
  min_display_tpm = CFG$min_display_tpm %||% 1,
  reference_colour  = trimws(CFG$reference_colour %||% "#4D4D4D"),
  comparison_colour = trimws(CFG$comparison_colour %||% "#56B4E9"),
  gene_flank_pct  = CFG$gene_flank_pct %||% 15,
  uniprot_reviewed_only = CFG$uniprot_reviewed_only %||% FALSE,
  arc_clip_partial = CFG$arc_clip_partial %||% FALSE
)

# Bolds whatever precedes the first " · " in a label (or the whole label if it has none); format-driven, not an enumerated word list (Phase C.0 decision 14), since bw.R's tracks can start with an arbitrary user-given instance name an enumerated list couldn't anticipate.
format_track_label <- function(label) {
  sub("^([^\u00B7]+?)( \u00B7 |$)", "<strong>\\1</strong>\\2", label)
}

# Gene symbol -> coords + user-facing error, shared by Compose's and Gene Report's gene-name lookup observers so the two can never disagree on error wording.
resolve_gene_region <- function(sym, txdb, flank_pct) {
  result <- gene_to_coords(sym, txdb = txdb, flank_pct = flank_pct, orgdb_pkg = ORGDB_PKG)
  if (is.null(result))
    list(region = NULL, error = paste0("'", sym, "' not found in ", GENOME$genome_label, " (e.g. Pax5, not PAX5)"))
  else
    list(region = result, error = NULL)
}

# Shared by Compose's and Gene Report's Generate buttons; both carry "generate-btn" so the spinner/double-click JS gate (which matches on that class, not a hardcoded id) covers both.
generate_button_ui <- function(id, enabled) {
  if (enabled)
    actionButton(id, "Generate figure \u2192", class = "btn-success btn-lg generate-btn")
  else
    tags$button("Generate figure \u2192",
                class = "btn btn-success btn-lg generate-btn", disabled = "disabled")
}

# Shared collapsible ".plottrellis-section" panel shell; the app-wide "hide.bs.collapse"/"show.bs.collapse" JS handlers key off any ".plottrellis-section .collapse" element (not a specific id), so this same shell works for Compose, Gene Report, and Output alike.
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

UNIQUE_SAMPLES <- unique(vapply(SAMPLES_DF, function(r) r$display, character(1)))
# Use colour column from Samples_super.txt when present; fall back to PALETTE
DEFAULT_SAMPLE_COLS <- setNames(
  vapply(seq_along(UNIQUE_SAMPLES), function(i) {
    s   <- UNIQUE_SAMPLES[i]
    row <- Find(function(r) r$display == s, SAMPLES_DF)
    if (!is.null(row) && !is.na(row$colour %||% NA_character_)) return(row$colour)
    PALETTE[((i - 1L) %% length(PALETTE)) + 1L]
  }, character(1)),
  UNIQUE_SAMPLES
)

# Main chromosomes first (chr1..chrN, chrX, chrY, chrM), then every
# alt/random/fix/chrUn contig alphabetically below them -- kept, not filtered
# (Phase I.1 step 6 decision).
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

# Converts TRACK_ITEMS entries (picker selections, utils/track_detection.R's shape) into pipeline/render_all()'s entry shape; height/data are always NULL here since render_all()'s settings_<viz>() functions compute height themselves.
build_entries <- function(figure_labels, col_map, settings) {
  Filter(Negate(is.null), lapply(figure_labels, function(lbl) {
    sp <- TRACK_ITEMS[[lbl]]
    if (is.null(sp)) return(NULL)

    # Comparison items carry `samples` = c(reference_display, comparison_display)
    # instead of a single colour_key. control_colour/treatment_colour are
    # resolved HERE, from two fixed per-project ROLE colours (settings$
    # reference_colour/comparison_colour), never from either condition's own
    # identity/PALETTE colour -- "which one is the baseline" reads the same
    # way in every comparison track regardless of the samples involved.
    if (!is.null(sp$samples) && length(sp$samples) == 2) {
      known <- c("key", "label", "type", "sample", "file", "colour_key", "samples",
                 "control_label", "treatment_label", "control_path", "treatment_path")
      extra <- sp[setdiff(names(sp), known)]  # e.g. junctions_mirrored_diff_tpm_scaled's control_quant/treatment_quant
      return(c(list(
        viz = sp$type, name = sp$name,
        control_label = sp$control_label, treatment_label = sp$treatment_label,
        control_path  = sp$control_path,  treatment_path  = sp$treatment_path,
        control_colour   = settings$reference_colour  %||% "#4D4D4D",
        treatment_colour = settings$comparison_colour %||% "#56B4E9",
        height = NULL, data = NULL
      ), extra))
    }

    col <- if (!is.na(sp$colour_key %||% NA)) col_map[[sp$colour_key]] else NA_character_
    entry <- list(viz = sp$type, height = NULL, data = NULL)
    if (sp$type == "tpm_usage") {
      entry <- c(entry, list(sample = sp$sample, colour = col, quant_path = sp$quant_path, junction_path = sp$junction_path))
    } else if (sp$type == "sashimi") {
      entry <- c(entry, list(sample = sp$sample, colour = col, name = sp$name, path = sp$junction_path))
    } else if (sp$type %in% c("bw", "hic", "bed")) {
      entry <- c(entry, list(sample = sp$sample, colour = col, name = sp$name, path = sp$file))
      if (sp$type == "hic") entry$hic_style <- settings$hic_style
    }
    entry
  }))
}

# ── Gene Report: condition list + per-condition track pool ────────────────────
# AVAILABLE_CONDITIONS now computed at startup, above (build_comparison_
# track_items() needs it too).

track_labels <- function(keys) {
  setNames(lapply(keys, function(k) HTML(format_track_label(TRACK_ITEMS[[k]]$label))), keys)
}

# Alphabetical (case-insensitive) by TRACK_ITEMS' label, split on " · "; never applied to figure_keys, whose order is deliberately curated. Full rationale: solutions/concepts/19-relocated-code-rationale.md, "app.R:511".
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

# Just the gene-annotation track for every sample-count mode (bug fix, 2026-09-15: this used to also default in tpm_usage/bw/sashimi/transcripts_annotated/protein_domains/diff tracks, swapping in a whole different stack on every mode/condition change instead of a plain, predictable starting point); anything else goes in from Available tracks explicitly.
default_figure_keys <- function(conditions) "genes"

# Item pool + default figure/available split for 1-2 conditions, filtered by each item's own `sample` field; fixed_item()s are found by scanning TRACK_ITEMS itself, never a hardcoded type list. Full rationale: solutions/concepts/19-relocated-code-rationale.md, "app.R:545".
track_split <- function(conditions) {
  cond_items <- Filter(function(sp) isTRUE(sp$sample %in% conditions), TRACK_ITEMS)
  comparison_keys <- character(0)
  if (length(conditions) > 1) {
    comparison_keys <- names(TRACK_ITEMS)[vapply(TRACK_ITEMS, function(sp) {
      !is.null(sp$samples) && length(sp$samples) == 2 && identical(sp$samples, conditions)
    }, logical(1))]
  }
  # is.null(sp$samples) is required alongside is.na(sp$sample) here -- every
  # viz_comparisons/*.R item ALSO sets sample = NA (it uses `samples`,
  # plural, instead), so without this every pairwise comparison track for
  # EVERY condition pair (not just the one selected above) would leak into
  # "fixed" and appear regardless of conditions -- confirmed empirically
  # (comparison_keys' pair-matching filter was being bypassed entirely: the
  # fixed bucket alone pulled in all ~470 comparison items across all pairs).
  fixed_keys  <- names(Filter(function(sp) isTRUE(is.na(sp$sample)) && is.null(sp$samples), TRACK_ITEMS))
  pool_keys   <- c(names(cond_items), comparison_keys, fixed_keys)
  figure_keys <- intersect(default_figure_keys(conditions), pool_keys)
  avail_keys  <- sort_keys_by_label(setdiff(pool_keys, figure_keys))
  list(available = track_labels(avail_keys), figure = track_labels(figure_keys))
}


# Multi/track-compose's "Available tracks" pool: TRACK_ITEMS minus "genes"
# (pinned straight to Figure) minus EVERY comparison-type ("diff") item.
# Supersedes design decision F (solutions/mvp/04-merge-app-modes.md), which
# kept all DE/DTU-backed comparison items in Compose's pool -- "Comparisons"
# is a 2-sample control-vs-treatment concept, pairwise mode already surfaces
# plenty of diff tracks for it, and Compose offering the same items again
# (for every pair, not just one) was redundant clutter, not a feature.
# Sorted alphabetically by track label (sort_keys_by_label()) so the "All" filter
# shows a stable, name-grouped order (ATAC vs RNA-seq vs Junctions, etc.)
# rather than TRACK_ITEMS' own sample-then-type insertion order.
COMPOSE_AVAILABLE_KEYS <- sort_keys_by_label(Filter(function(k) {
  sp <- TRACK_ITEMS[[k]]
  if (identical(k, "genes")) return(FALSE)
  is.null(sp$samples) || length(sp$samples) != 2
}, names(TRACK_ITEMS)))

# Choices for the Gene Report's Control/Treatment selects: every AVAILABLE_CONDITIONS entry is always offered plus a blank placeholder; collisions (picking the value already sat in the other box) are resolved by clearing the other box server-side, not by removing options here.
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
  # width:100% below the cap. No max-height/overflow-y here (deliberately --
  # that was tried before and made tall figures feel clipped/hidden).
  div(style = "width:100%; max-width:1400px; margin:0 auto; overflow-x:auto;", HTML(svg_text))
}

# JS map from track key → type, embedded in page for client-side filter
TRACK_TYPE_JS <- jsonlite::toJSON(
  setNames(lapply(TRACK_ITEMS, function(x) x$type), names(TRACK_ITEMS)),
  auto_unbox = TRUE
)
# JS map from track key -> broken (Phase C.2 step 8): a declared-but-
# unresolvable path (viz/bw.R, viz/hic.R, viz/bed.R, viz/sashimi.R,
# viz/tpm_usage.R's own `broken` field, computed live from SAMPLES_DF at
# startup) gets red-boxed in the track picker, same treatment as the
# existing gene-required styling. isTRUE() coerces the common case (a track
# type with no `broken` field at all, e.g. genes/transcripts/hic comparison
# tracks) from NULL to FALSE, rather than needing every module to declare it.
TRACK_BROKEN_JS <- jsonlite::toJSON(
  setNames(lapply(TRACK_ITEMS, function(x) isTRUE(x$broken)), names(TRACK_ITEMS)),
  auto_unbox = TRUE
)
# The "Comparisons" filter button (Gene Report only) matches any of these
# distinct per-modality comparison types at once, not a single literal type
# string -- each viz_comparisons/*.R module has its own VIZ_TYPE (so
# assemble.R's per-type render dispatch works), but the button groups them
# for the biologist under one label.
COMPARISON_VIZ_TYPES_JS <- jsonlite::toJSON(
  c("bw_diff_stacked", "junctions_mirrored_diff", "junctions_mirrored_diff_tpm_scaled")
)
# "Annotations" (Gene Report + Compose) groups gene annotation with all
# transcript-structure tracks (filtered, unfiltered, and filtered+annotated)
# under one label -- there is no standalone "Transcripts" filter button.
ANNOTATION_VIZ_TYPES_JS <- jsonlite::toJSON(c("genes", "transcripts", "transcripts_all", "transcripts_annotated", "protein_domains"))
# Gene-required track types -- the pipeline/ vizzes that stop() on a missing
# gene (settings_transcripts_all/transcripts/transcripts_annotated/
# protein_domains/tpm_usage load via load_transcripts(), settings_junctions_
# mirrored_diff_tpm_scaled needs filtered_transcripts() for its TPM scale).
# Kept as a plain R vector too (not just the JSON string below) so server()
# can run the same red-box check the client already does, to gate Generate
# on it (E.2 follow-up, 2026-09-08).
GENE_REQUIRED_TYPES <- c("transcripts_all", "transcripts", "transcripts_annotated",
                          "protein_domains", "tpm_usage", "junctions_mirrored_diff_tpm_scaled")
GENE_REQUIRED_TYPES_JS <- jsonlite::toJSON(GENE_REQUIRED_TYPES)

# settings_modal() removed (Phase C.3 decision 7) -- every field it exposed
# now lives in utils/wizard.R's "Other settings" section instead, since
# editing goes through the wizard's restart+Save round trip, not a
# session-only, never-written-to-disk override. Per-sample colours are set
# per-sample in the wizard (remove + re-add to change one) rather than via a
# separate live swatch section here.

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
    .section-head  { font-size: var(--bio-font-heading); font-weight:600; color:#555; margin-bottom:6px; }
    .bucket-label  { font-size: var(--bio-font-label); font-weight:600; margin-bottom:4px; }
    .zoom-group-label { font-size: var(--bio-font-label); font-weight:600; color:#555;
                        margin-bottom:5px; letter-spacing:0.03em; }
    .rank-list-container { min-height: 60px; max-height: 42vh; overflow-y: auto; }
    .rank-list-item { cursor: grab; }
    .rank-list-item:active { cursor: grabbing; }
    .track-filters { margin-bottom:10px; display:flex; flex-wrap:wrap; gap:4px; }
    .track-filters .btn { flex:1 1 0; margin-right:0; }
    .track-filters .btn.active { background-color:#337ab7; color:#fff; border-color:#2e6da4; }
    /* Gene-dependent tracks (Transcripts (filtered)/Transcripts/Transcripts
       annotated (filtered)/Canonical protein annotation/Transcript counts/
       Junctions mirrored diff (TPM-scaled), see app.R's GENE_REQUIRED_TYPES_JS)
       while region mode is Coordinates, not Gene name -- flagged, not hidden,
       since a track can already be sitting in
       Figure tracks when the mode switches away from Gene name; Generate will
       still hard-fail on these if left in place. Applies in every sample-count
       mode, not just Compose -- Phase E.2 step 6 (Q6) lifted the old
       Coordinates-is-compose-only restriction, and applyGeneRequiredFlag()
       scopes to the single #tracks-block shared by all three modes. */
    /* Phase C.2 step 8: a declared-but-unresolvable path (moved/deleted file)
       gets the same red-box treatment as the gene-required flag above --
       deliberately reusing the same look rather than inventing a new one,
       since Phase E hasn't decided error-state visuals yet (sign-off
       decision 3, solutions/phase-logs/phase-C.2.md). Applies in every view
       mode's picker, same as track-gene-required above. */
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
       yet (Phase E.2 step 5; conditions added to the gate 2026-09-08). */
    .shiny-input-container .radio-inline.disabled-radio-label {
      opacity: 0.45; cursor: not-allowed;
    }
    /* Track picker (Available/Figure rank-lists), dimmed by setPickerLocked()
       while mode+conditions+region aren't all resolved yet (Phase E.2 step 6
       + Q7; conditions added to the gate 2026-09-08). */
    .rank-list-locked {
      opacity: 0.45; cursor: not-allowed;
    }
  "))),
  tags$head(tags$script(HTML(paste0("
    var TRACK_TYPE_MAP = ", TRACK_TYPE_JS, ";
    var TRACK_BROKEN_MAP = ", TRACK_BROKEN_JS, ";
    var COMPARISON_VIZ_TYPES = ", COMPARISON_VIZ_TYPES_JS, ";
    var ANNOTATION_VIZ_TYPES = ", ANNOTATION_VIZ_TYPES_JS, ";
    var GENE_REQUIRED_TYPES = ", GENE_REQUIRED_TYPES_JS, ";

    // Track-filter buttons + bucket are scoped per [data-available-id] block,
    // so Compose's and Gene Report's independent pickers never interfere.
    function applyTrackFilter(blockEl) {
      var $block  = $(blockEl);
      var availId = $block.attr('data-available-id');
      var f = $block.find('.filter-btn.active').data('filter') || 'all';
      $block.find('#' + availId + ' .rank-list-item').each(function() {
        var key  = $(this).attr('data-rank-id') || '';
        var type = TRACK_TYPE_MAP[key] || '';
        var matches = f === 'all' || type === f ||
          (f === 'comparison'  && COMPARISON_VIZ_TYPES.indexOf(type) !== -1) ||
          (f === 'annotation'  && ANNOTATION_VIZ_TYPES.indexOf(type) !== -1);
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

    // rationale: see solutions/concepts/19-relocated-code-rationale.md, \"app.R:812\"
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

    // rationale: see solutions/concepts/19-relocated-code-rationale.md, \"app.R:841\"
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

    // rationale: see solutions/concepts/19-relocated-code-rationale.md, \"app.R:871\"
    function setConditionGateLocked(locked) {
      var $regionModeInputs = $('input[name=\"region_mode\"]');
      $regionModeInputs.prop('disabled', locked);
      $regionModeInputs.closest('.radio-inline').toggleClass('disabled-radio-label', locked);
      var geneSel = document.getElementById('gene_name');
      if (geneSel && geneSel.selectize) {
        if (locked) geneSel.selectize.disable(); else geneSel.selectize.enable();
      }
    }

    // Phase E.2 step 6 + Q7 (conditions added to the gate 2026-09-08): the
    // track picker locks (both Available and Figure tracks, not just
    // visually dimmed) whenever the server-side region_lock message says the
    // mode+conditions+region gate isn't satisfied yet -- see the
    // `region_lock` custom message handler below for how this actually
    // disables dragging via each rank-list's own Sortable instance.
    function setPickerLocked(locked) {
      ['available_tracks', 'figure_tracks'].forEach(function(id) {
        var el = document.getElementById(id);
        if (!el) return;
        var inst = window.Sortable && window.Sortable.get(el);
        if (inst) inst.option('disabled', locked);
        $(el).toggleClass('rank-list-locked', locked);
      });
    }

    // toggleModeGate()/condition-gate state no longer needs a client-side
    // trigger here: setConditionGateLocked() (like setPickerLocked()) is
    // driven entirely by the server's `condition_lock` message, which fires
    // on its own whenever view_mode/sel_condition/sel_reference/
    // sel_comparison change (see has_conditions_ready() in server()).
    $(document).ready(function() {
      setTimeout(function() { applyAllTrackFilters(); applyGeneRequiredFlag(); applyBrokenFlag(); }, 400);
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

    // The Comparisons button only exists in two-sample mode (see the
    // conditionalPanel around it) -- if it was active when the user switches
    // to single/compose, it's now hidden with no way to click back to 'All',
    // so reset the filter for them instead of leaving it stuck.
    $(document).on('change', 'input[name=\"view_mode\"]', function() {
      if ($(this).val() === 'two') return;
      var $block = $('#tracks-block');
      var $cmp   = $block.find('.filter-btn[data-filter=\"comparison\"]');
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
      // .generate-btn, not a hardcoded #generate id -- Compose's and Gene
      // Report's Generate buttons both carry this class so neither mode's
      // render silently skips the spinner/double-click gate.
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

    // Phase E.2 step 6 + Q7: server-side region_lock drives the track
    // picker's lock state (see setPickerLocked() above and the
    // has_mode()/has_conditions_ready()/resolved_region()-driven observe()
    // in server()).
    Shiny.addCustomMessageHandler('region_lock', function(x) { setPickerLocked(x.locked); });

    // E.2 follow-up (2026-09-08): server-side condition_lock drives the
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

    # Settings/reset buttons, fixed to top-right navbar area via CSS. Reset
    # used to live inside settings_modal()'s "Danger zone" section; moved up
    # here as its own button now that settings_modal() is gone (Phase C.3
    # decision 7) -- it's a separate, more destructive action from editing
    # settings, not a sub-item of it.
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

    # ── Explore: Compose / Gene Report, unified behind one 3-way mode radio ─────
    # See solutions/mvp/04-merge-app-modes.md for the design rationale --
    # "single"/"two" are the old Gene Report modes (data-layer-locked pool,
    # narrowed to 1-2 conditions); "compose" is the old Track Composer
    # (unrestricted multi-sample pool, no auto-defaulted tracks).
    collapsible_panel("compose", "Compose",

      h4("Number of samples"),
      radioButtons("view_mode", NULL,
                   choices  = c("One" = "single",
                                "Two" = "two",
                                "Multi" = "compose"),
                   selected = character(0), inline = TRUE),
      hr(),

      # Condition selection -- single/two-sample modes only; track-compose has
      # no dropdown, the track picker below is the only sample selector there.
      # Placed above Genomic region so control/treatment (or condition) is
      # chosen before searching for a gene. All start unselected (E.2
      # follow-up, 2026-09-08) -- region_mode/gene_name and the track picker
      # now stay locked until these resolve too, not just until a mode is
      # chosen (see has_conditions_ready() in server()).
      conditionalPanel("input.view_mode === 'single'",
        fluidRow(column(4,
          selectizeInput("sel_condition", "Condition:",
                      choices = condition_choices(), selected = "")
        ))
      ),
      conditionalPanel("input.view_mode === 'two'",
        fluidRow(
          column(4,
            selectizeInput("sel_reference", "Reference:",
                        choices = condition_choices(), selected = "")
          ),
          column(4,
            selectizeInput("sel_comparison", "Comparison:",
                        choices = condition_choices(),
                        selected = "")
          )
        )
      ),

      # Region (full width)
      fluidRow(
        column(12,
          h4("Genomic region"),
          radioButtons("region_mode", NULL,
                       choices  = c("Gene name" = "gene", "Coordinates" = "coords"),
                       selected = character(0), inline = TRUE),
          hr(style = "margin:8px 0;"),

          conditionalPanel("input.region_mode === 'gene'",
            fluidRow(
              column(4,
                selectizeInput("gene_name", "Gene symbol:",
                               choices = NULL,
                               # maxOptions caps how many server-side search
                               # matches selectize will ever render, not a
                               # "load more on scroll" limit -- 10 meant most
                               # typed queries silently truncated well before
                               # the user could scroll to a match.
                               options = list(placeholder = "e.g. Pax5", maxOptions = 1000))
              ),
              column(6, br(), uiOutput("gene_result_ui"))
            )
          ),

          # Coordinates panel, shown whenever region_mode is 'coords' -- selectable
          # in every sample-count mode as of Phase E.2 step 6 (Q6); tracks that need
          # a gene symbol simply won't render from coordinates (Phase E.3 red-boxing).
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
              # Two-sample mode only -- "Comparisons" is a control/treatment
              # concept; one-sample has none, and track-compose's pool
              # (COMPOSE_AVAILABLE_KEYS) now excludes every comparison-type
              # item outright, so there is nothing for this filter to select
              # there either (supersedes design decision F,
              # solutions/mvp/04-merge-app-modes.md).
              conditionalPanel("input.view_mode === 'two'",
                tags$button("Comparisons", class = "btn btn-lg btn-default filter-btn", `data-filter` = "comparison")
              ),
              tags$button("Annotations", class = "btn btn-lg btn-default filter-btn", `data-filter` = "annotation"),
              tags$button("BED",         class = "btn btn-lg btn-default filter-btn", `data-filter` = "bed"),
              # Phase C.0: ATAC/ChIP-seq/CUT&RUN/RNA-seq coverage collapsed
              # into one generic bw VIZ_TYPE -- a track's modality is now in
              # its NAME (bold heading text), not a distinguishable type, so
              # one "Coverage" button (data-filter = "bw") replaces the old
              # separate "ATAC"/"RNA-seq" buttons and also newly covers
              # ChIP-seq/CUT&RUN, which never had a dedicated button before
              # (12-data-intake-spec.md's known gap).
              tags$button("Coverage",    class = "btn btn-lg btn-default filter-btn", `data-filter` = "bw"),
              tags$button("HiC",         class = "btn btn-lg btn-default filter-btn", `data-filter` = "hic"),
              tags$button("Sashimi",     class = "btn btn-lg btn-default filter-btn", `data-filter` = "sashimi"),
              tags$button("TPMs",        class = "btn btn-lg btn-default filter-btn", `data-filter` = "tpm_usage")
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

      fluidRow(column(4,
        textInput("plot_title", "Title (optional):",
                  placeholder = "Leave blank to auto-generate", width = "100%")
      )),
      tags$div(style = "height:8px;"),
      fluidRow(
        column(4, uiOutput("generate_ui")),
        column(8, uiOutput("generate_warning_ui"))
      ),
      tags$div(style = "height:16px;"),

      # rationale: see solutions/concepts/19-relocated-code-rationale.md, "app.R:1262"
      conditionalPanel("(input.view_mode === 'single' || input.view_mode === 'two') && input.region_mode === 'gene'",
        collapsible_panel("gene_protein_info", "Gene / protein info", collapsed = TRUE,
          tabsetPanel(type = "pills",
            tabPanel("Gene info", uiOutput("gr_extras_ui")),
            tabPanel("Canonical protein annotation", uiOutput("gr_protein_ui"))
          )
        )
      )
    ), # /collapsible_panel: compose

    # ── Output section (hidden until first render) ───────────────────────────────
    # Full-width, no inner scrollbox -- the figure renders at its natural size
    # and the whole page scrolls to reveal more of it (biologist ask: the old
    # max-height/overflow-y box made tall figures feel clipped/hidden).
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
  # the only way to start with nothing shown. No gene preselected as of the
  # E.2 follow-up (2026-09-08) -- previously defaulted to "Pax5".
  updateSelectizeInput(session, "gene_name", choices = ALL_GENES, server = TRUE, selected = "")

  # ── Mode switch: single/two/compose is now one radio (input$view_mode),
  # read directly -- no derived reactiveVal/output needed, same pattern
  # already used for input$region_mode elsewhere in this file. Coordinates
  # is selectable in every sample-count mode as of Phase E.2 step 6 (Q6), so
  # unlike before there's no stale-selection reset to make here on mode switch.

  # Per-sample colour state (session-only; defaults from PALETTE)
  session_colours <- do.call(reactiveValues, as.list(DEFAULT_SAMPLE_COLS))

  observe({
    for (s in UNIQUE_SAMPLES) {
      safe_id <- paste0("scolour_", gsub("[^A-Za-z0-9]", "_", s))
      val     <- input[[safe_id]]
      if (!is.null(val) && nchar(val) > 0) session_colours[[s]] <- val
    }
  })

  colour_map_r <- function(figure_labels) {
    # Comparison items carry `samples` (both condition display names) but,
    # unlike single-sample items' `colour_key`, these are NOT fed into the key
    # universe here -- control_colour/treatment_colour are fixed per-project
    # ROLE colours resolved directly from settings in build_entries(), not
    # looked up per-sample. Only single-sample `colour_key`s need a live
    # session colour.
    keys <- unlist(lapply(figure_labels, function(lbl) {
      sp <- TRACK_ITEMS[[lbl]]
      if (is.null(sp)) return(character(0))
      sp$colour_key
    }))
    unique_keys <- unique(keys[!is.na(keys)])
    cols <- vapply(seq_along(unique_keys), function(i) {
      k <- unique_keys[i]
      v <- isolate(session_colours[[k]])
      if (is.null(v)) PALETTE[((i - 1L) %% length(PALETTE)) + 1L] else v
    }, character(1))
    setNames(cols, unique_keys)
  }

  state <- reactiveValues(
    gene_region          = NULL,
    gene_error           = NULL,
    region               = NULL,
    svg_path             = NULL,
    rcode                = NULL,
    plot_error           = NULL,
    start_parse_err      = NULL,
    end_parse_err        = NULL,
    last_region          = NULL,
    last_entries         = NULL,
    last_settings_base   = NULL,
    rendering            = FALSE
  )

  session$onSessionEnded(function() {
    isolate({
      if (!is.null(state$svg_path) && file.exists(state$svg_path))
        unlink(state$svg_path)
    })
  })

  # ── Settings (re-enter wizard) ────────────────────────────────────────────────
  # Phase C.3 decision 1/Extra check 3: editing settings now restarts the
  # whole process (the wizard needs SAMPLES_DF/txdb/etc. to not exist yet, see
  # utils/wizard.R's header comment) -- confirm first so this isn't a silent
  # surprise mid-session, rather than restarting immediately on click.
  observeEvent(input$open_settings, {
    showModal(modalDialog(
      title = "Edit project settings?", easyClose = FALSE,
      p("This restarts the app into the setup wizard, pre-filled with your current settings. Any in-progress figure/selection in this session will be lost."),
      footer = tagList(modalButton("Cancel"), actionButton("confirm_edit_settings", "Continue", class = "btn-primary"))
    ))
  })

  observeEvent(input$confirm_edit_settings, {
    Sys.setenv(PLOTTRELLIS_EDIT_SETTINGS = "1")
    # Same "show feedback, then delay stopApp()" pattern as
    # confirm_reset_project below -- an immediate stopApp() can kill the
    # process before the modal's removal/replacement actually
    # reaches the browser.
    showModal(modalDialog(
      title = "Restarting into the setup wizard", easyClose = FALSE, footer = NULL,
      p("A new browser tab will open in a few seconds.")
    ))
    later::later(function() stopApp(returnValue = "edit_settings"), 0.2)
  }, ignoreInit = TRUE)

  observeEvent(input$reset_project_btn, {
    showModal(modalDialog(
      title = "Reset project configuration?", easyClose = FALSE,
      p("This deletes project_settings.json and restarts the app into first-time setup. samples.txt itself is not touched."),
      footer = tagList(modalButton("Cancel"), actionButton("confirm_reset_project", "Delete and restart", class = "btn-danger"))
    ))
  })

  observeEvent(input$confirm_reset_project, {
    if (file.exists(SETTINGS_JSON)) file.remove(SETTINGS_JSON)
    showModal(modalDialog(
      title = "Project configuration reset", easyClose = FALSE, footer = NULL,
      p("project_settings.json has been deleted."),
      p("Restarting into first-time setup, a new browser tab will open in a few seconds.")
    ))
    later::later(function() stopApp(returnValue = "restart"), 0.2)
  }, ignoreInit = TRUE)

  # ── Close app ────────────────────────────────────────────────────────────────
  # Same confirm-then-delayed-stopApp() pattern as reset/edit-settings, but
  # with a returnValue that ISN'T "restart"/"edit_settings", so launch.R's
  # repeat loop breaks instead of relaunching: this actually ends the R
  # process, not just the tab. Manual only -- no auto-stop on tab close/
  # disconnect, since a flaky connection dropping and reconnecting would
  # otherwise kill the app out from under the user.
  observeEvent(input$close_app_btn, {
    showModal(modalDialog(
      title = "Close app?", easyClose = FALSE,
      p("This stops the R process running this app. You'll need to relaunch it to use PlotTrellis again."),
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
  # rationale: see solutions/concepts/19-relocated-code-rationale.md, "app.R:1454"
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
      footer = tagList(modalButton("Cancel"), actionButton("confirm_clear_caches", "Delete and restart", class = "btn-warning"))
    ))
  })

  observeEvent(input$confirm_clear_caches, {
    cache_files <- c(GENOME$txdb_path, GENOME$tx_name_lookup_path, GENOME$tx_meta_lookup_path,
                      GENE_DISPLAY_TRANSCRIPTS_CACHE_PATH, UNIPROT_CANONICAL_CACHE_PATH)
    for (f in cache_files) if (!is.null(f) && file.exists(f)) file.remove(f)
    if (dir.exists(PROTEIN_DOMAINS_CACHE_DIR)) unlink(PROTEIN_DOMAINS_CACHE_DIR, recursive = TRUE)

    showModal(modalDialog(
      title = "Caches cleared", easyClose = FALSE, footer = NULL,
      p("All caches have been deleted."),
      p("Restarting with your current project settings, a new browser tab will open in a few seconds.")
    ))
    later::later(function() stopApp(returnValue = "restart"), 0.2)
  }, ignoreInit = TRUE)

  # cur_settings() removed (Phase C.3 decision 7) -- every call site below
  # now reads the top-level RENDER_CFG list directly instead.

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
  observe({
    sym <- trimws(input$gene_name %||% "")
    if (nchar(sym) == 0) return()
    pct <- RENDER_CFG$gene_flank_pct / 100
    res <- resolve_gene_region(sym, txdb, pct)
    state$gene_region <- res$region
    state$gene_error  <- res$error
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
                 " (", format_bp(r$end - r$start), ")"))
    }
  })

  # rationale: see solutions/concepts/19-relocated-code-rationale.md, "app.R:1434"
  output$gr_extras_ui <- renderUI({
    if (is.null(state$gene_region)) return(NULL)
    if (!is.null(input$region_mode) && input$region_mode != "gene") return(NULL)
    gene <- trimws(input$gene_name %||% "")
    if (nchar(gene) == 0) return(NULL)

    ann <- tryCatch(transcript_annotations(gene, txdb), error = function(e) NULL)

    if (!is.null(ann) && nrow(ann) > 0) {
      # UniProt check is a live network call (see utils/uniprot_canonical.R) --
      # failure degrades to NA/"unavailable" per row, never blanks out the
      # GTF-derived columns above, which don't depend on it.
      uc <- tryCatch(uniprot_canonical_match(gene, ann$tx_id, CDS_BY_TX, ann$ccds_id, orgdb = RENDER_CFG$orgdb_pkg), error = function(e) NULL)
      if (is.null(uc) || nrow(uc) != nrow(ann)) {
        uc <- data.frame(uniprot_accession = rep(NA_character_, nrow(ann)),
                          uniprot_match = rep(NA, nrow(ann)),
                          uniprot_unique = rep(NA, nrow(ann)),
                          uniprot_match_count = rep(NA_integer_, nrow(ann)),
                          uniprot_source = rep(NA_character_, nrow(ann)),
                          uniprot_canonical_isoform_id = rep(NA_character_, nrow(ann)),
                          uniprot_canonical_isoform_name = rep(NA_character_, nrow(ann)),
                          uniprot_row_isoform_id = rep(NA_character_, nrow(ann)),
                          uniprot_row_isoform_name = rep(NA_character_, nrow(ann)))
      }
      ann[names(uc)] <- uc
    }

    tag_pill <- function(text) tags$span(
      style = paste("display:inline-block; background:#eef2f7; color:#345; border-radius:3px;",
                     "padding:1px 6px; margin-right:4px; font-size:0.9em;"),
      text
    )

    uniprot_source_label <- c(
      ensembl_xref = "direct UniProt-Ensembl cross-reference",
      ccds_xref    = "UniProt-CCDS cross-reference",
      length_match = "protein length match only, not a direct cross-reference"
    )

    tagList(
      tags$a(href = uniprot_search_url(gene, RENDER_CFG$orgdb_pkg, RENDER_CFG$uniprot_reviewed_only), target = "_blank",
             class = "btn btn-default btn-sm",
             paste0("UniProt: ", gene, " (opens in new tab)")),
      if (!is.null(ann) && nrow(ann) > 0) {
        tags$div(style = "margin-top:8px; font-size: var(--bio-font-code, 12px);",
          tags$div(style = "display:flex; gap:10px; font-weight:600; color:#666; border-bottom:1px solid #ccc; padding-bottom:3px; margin-bottom:2px;",
            tags$span(style = "flex:0 0 90px;", "Transcript"),
            tags$span(style = "flex:0 0 60px;", "Type"),
            tags$span(style = "flex:1;", "Tags")
          ),
          lapply(seq_len(nrow(ann)), function(i) {
            src <- ann$uniprot_source[i]

            # The specific UniProt sequence this whole check compares
            # against -- an entry-level fact (isoform id + name, or just the
            # bare accession for a single-isoform entry), true regardless of
            # whether THIS transcript carries a direct tag for it.
            canon_id   <- ann$uniprot_canonical_isoform_id[i]
            canon_name <- ann$uniprot_canonical_isoform_name[i]
            canon_desc <- if (is.na(canon_id)) NA_character_ else
              if (is.na(canon_name)) canon_id else sprintf("isoform %s (%s)", canon_name, canon_id)

            uni_text <- if (is.na(ann$uniprot_match[i])) {
              "unknown (UniProt lookup unavailable)"
            } else if (ann$uniprot_match[i] && ann$uniprot_unique[i]) {
              sprintf("yes, this is UniProt's canonical %s (via %s)", canon_desc, uniprot_source_label[[src]])
            } else if (ann$uniprot_match[i]) {
              sprintf("matches UniProt's canonical %s (via %s), but %d transcripts of this gene share that exact sequence, so it doesn't uniquely identify this one",
                      canon_desc, uniprot_source_label[[src]], ann$uniprot_match_count[i])
            } else if (!is.na(src) && src != "length_match") {
              row_id   <- ann$uniprot_row_isoform_id[i]
              row_name <- ann$uniprot_row_isoform_name[i]
              row_desc <- if (is.na(row_id)) ann$uniprot_accession[i] else
                if (is.na(row_name)) row_id else sprintf("isoform %s (%s)", row_name, row_id)
              sprintf("no, this is %s, confirmed by %s; UniProt's canonical is %s",
                      row_desc, uniprot_source_label[[src]], canon_desc)
            } else if (!is.na(ann$uniprot_accession[i])) {
              sprintf("no; UniProt's canonical is %s", canon_desc)
            } else {
              "no reviewed UniProt entry found"
            }

            row_tags <- character(0)
            if (isTRUE(ann$ensembl_canonical[i])) row_tags <- c(row_tags, "Ensembl canonical")
            if (!is.na(ann$appris_tier[i])) row_tags <- c(row_tags, paste("APPRIS", ann$appris_tier[i]))
            if (isTRUE(ann$uniprot_match[i])) {
              row_tags <- c(row_tags, if (isTRUE(ann$uniprot_unique[i])) "UniProt canonical" else "UniProt (shared match)")
            }

            tags$details(
              style = "border-bottom: 1px solid #eee; padding: 3px 0;",
              tags$summary(
                style = "display:flex; gap:10px;",
                tags$span(style = "flex:0 0 90px; font-weight:600;", ann$transcript_name[i]),
                tags$span(style = "flex:0 0 60px;", ann$coding_status[i]),
                tags$span(style = "flex:1;", lapply(row_tags, tag_pill))
              ),
              tags$div(style = "margin: 4px 0 4px 14px; color:#444;",
                tags$div(tags$strong("Ensembl canonical: "), if (ann$ensembl_canonical[i]) "yes" else "no"),
                tags$div(tags$strong("APPRIS: "), if (is.na(ann$appris_tier[i])) "none" else ann$appris_tier[i]),
                tags$div(tags$strong("CCDS: "), if (is.na(ann$ccds_id[i])) "none" else ann$ccds_id[i]),
                tags$div(tags$strong("UniProt canonical: "), uni_text),
                tags$div(style = "margin-top:3px; font-style:italic; color:#888;",
                  "Source: GENCODE annotation tags, UniProt cross-reference (see UniProt canonical line above for which kind).")
              )
            )
          })
        )
      }
    )
  })

  # ── Canonical protein annotation: UniProt domain/region/motif deep dive ─────
  # Same progressive-disclosure style as gr_extras_ui above (one <details> per
  # feature), sits to its right. Uses gene_protein_domains()
  # (utils/protein_domains.R) directly -- the same resolver the
  # "Canonical protein annotation" TRACK calls, but here read straight from
  # its per-gene RDS cache rather than re-inlined into a generated script, so
  # this is a single cheap read after the first view of a gene. One row per
  # domain_id (a multi-exon feature collapsed to one row, not one row per
  # exon segment) -- matches how the track itself groups them.
  output$gr_protein_ui <- renderUI({
    if (is.null(state$gene_region)) return(NULL)
    if (!is.null(input$region_mode) && input$region_mode != "gene") return(NULL)
    gene <- trimws(input$gene_name %||% "")
    if (nchar(gene) == 0) return(NULL)

    pd <- tryCatch(
      gene_protein_domains(gene, txdb, RENDER_CFG$gtf_path, RENDER_CFG$orgdb_pkg, RENDER_CFG$genome_label),
      error = function(e) list(status = "error", message = conditionMessage(e))
    )

    pd_msg <- switch(pd$status,
      no_uniprot_entry     = "No reviewed UniProt entry found for this gene",
      no_domains           = "No annotated domain/region/motif features for this gene's canonical isoform",
      no_canonical_transcript = "UniProt's canonical isoform could not be matched to any annotated transcript for this gene",
      ok = NULL,
      "Protein domain lookup failed (network issue?) -- will retry on next render"
    )

    method_label <- if (identical(pd$status, "ok")) switch(pd$method,
      ccds_id_join                    = "via CCDS ID",
      length_match_fallback           = "via protein-length match (CCDS not tagged in GTF yet)",
      length_match_fallback_ambiguous = "via protein-length match (ambiguous, CDS shared by multiple transcripts)",
      pd$method
    ) else NULL

    tagList(
      if (!is.null(pd_msg))
        div(style = "color:#888; font-style:italic;", pd_msg)
      else {
        domains <- pd$domains
        dids <- unique(domains$domain_id)
        dids <- dids[order(vapply(dids, function(d) min(domains$start[domains$domain_id == d]), numeric(1)))]

        tagList(
          div(style = "font-size:0.9em; color:#666; margin-bottom:6px;",
              sprintf("canonical: %s (%s)", pd$canonical_tx, method_label)),
          tags$div(style = "font-size: var(--bio-font-code, 12px);",
            tags$div(style = "display:flex; gap:10px; font-weight:600; color:#666; border-bottom:1px solid #ccc; padding-bottom:3px; margin-bottom:2px;",
              tags$span(style = "flex:0 0 14px;", ""),
              tags$span(style = "flex:1;", "Feature"),
              tags$span(style = "flex:0 0 60px;", "Exons")
            ),
            lapply(dids, function(did) {
              seg  <- domains[domains$domain_id == did, ]
              seg  <- seg[order(seg$start), ]
              ty   <- seg$type[1]
              col  <- PD_TYPE_COLOURS[[ty]]
              nseg <- nrow(seg)

              tags$details(
                style = "border-bottom: 1px solid #eee; padding: 3px 0;",
                tags$summary(
                  style = "display:flex; gap:10px; align-items:center;",
                  tags$span(style = sprintf(
                    "flex:0 0 12px; width:12px; height:12px; border-radius:2px; background:%s; display:inline-block;",
                    col
                  )),
                  tags$span(style = "flex:1;", seg$note[1]),
                  tags$span(style = "flex:0 0 60px; color:#888;", if (nseg > 1) sprintf("%d exons", nseg) else "")
                ),
                tags$div(style = "margin: 4px 0 4px 22px; color:#444;",
                  tags$div(tags$strong("Type: "), ty),
                  tags$div(tags$strong("Position: "), sprintf("aa %d-%d", seg$aa_start[1], seg$aa_end[1])),
                  tags$div(tags$strong("Exon segments: "), nseg)
                )
              )
            })
          )
        )
      }
    )
  })

  # ── Two-sample reference/comparison mutual exclusion ──
  # Both boxes always offer every condition; the only guard is against picking the same one in
  # both -- sel_prev tracks each box's last-settled value so we can tell which box just changed
  # (vs. one only just relabeled by this same observer) and clear the OTHER box on a collision.
  # ignoreInit = TRUE since the initial choices/markers are already baked into the UI; both
  # start "" (E.2 follow-up, 2026-09-08: Reference no longer preselects a default condition),
  # matching sel_reference/sel_comparison's own defaults.
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

  # ── Gene Report: which condition(s) are currently selected, mode-dependent ──
  selected_conditions <- reactive({
    if (identical(input$view_mode, "two")) c(input$sel_reference, input$sel_comparison)
    else input$sel_condition
  })

  # Whether the condition-selection stage is fully resolved for the current
  # mode (E.2 follow-up, 2026-09-08): Compose has no condition dropdown at
  # all (the track picker is its own sample selector), so it's trivially
  # ready; single/two need their dropdown(s) actually filled in, since none
  # of sel_condition/sel_reference/sel_comparison default to a preselected
  # condition any more. Shared by the region/gene lock, the track-pool
  # populator, and Generate below so all three agree on the same stage
  # boundary.
  has_conditions_ready <- reactive({
    if (identical(input$view_mode, "compose")) TRUE
    else if (identical(input$view_mode, "two"))
      nchar(input$sel_reference %||% "") > 0 && nchar(input$sel_comparison %||% "") > 0
    else if (identical(input$view_mode, "single")) nchar(input$sel_condition %||% "") > 0
    else FALSE
  })

  # rationale: see solutions/concepts/19-relocated-code-rationale.md, "app.R:1879"
  figure_pool_key <- reactiveVal(NULL)

  observeEvent(list(input$view_mode, input$sel_condition, input$sel_reference, input$sel_comparison), {
    if (identical(input$view_mode, "compose")) {
      key   <- "compose"
      split <- list(
        available = track_labels(COMPOSE_AVAILABLE_KEYS),
        figure    = list(genes = HTML(format_track_label(TRACK_ITEMS[["genes"]]$label)))
      )
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
  # Phase E.2 Q6: Coordinates is selectable in every mode now, not just
  # Compose, so which region resolves depends only on region_mode -- never
  # on view_mode. gene_name is only meaningful in 'gene' region_mode; tracks
  # that need a gene symbol (GENE_REQUIRED_TYPES) simply won't render from
  # coordinates in any mode (Q6's own note -- handled by Phase E.3's
  # red-boxing, not by this resolution logic).
  resolved_region <- reactive({
    if (!is.null(input$region_mode) && input$region_mode == "coords") {
      v <- tryCatch(region_valid(), error = function(e) NULL)
      if (!is.null(v) && isTRUE(v$valid))
        list(rgn = list(chrom = v$chrom, start = v$start, end = v$end), gene_name = NA_character_)
      else NULL
    } else if (!is.null(state$gene_region)) {
      list(rgn = state$gene_region, gene_name = input$gene_name)
    } else NULL
  })

  # Whether a sample-count mode has been chosen at all (Phase E.2 step 4:
  # view_mode starts unselected) -- Generate and the locks below all stay
  # gated on this in addition to conditions/region, so the intended stage
  # order (mode -> conditions -> region/gene -> tracks -> generate) actually
  # holds even though state$gene_region can already be non-NULL before any
  # mode is chosen (gene_name has no startup default any more, but a user
  # could still type a gene while unselected before this gate existed).
  has_mode <- reactive({ !is.null(input$view_mode) })

  # E.2 follow-up (2026-09-08): region_mode/gene_name lock until mode AND
  # conditions are both resolved -- previously this only waited on mode
  # (Phase E.2 step 5), which let a user search a gene before Reference/
  # Comparison (or Condition) were chosen. Mirrors the 'region_lock' message
  # below (same server-driven pattern), consumed by setConditionGateLocked().
  observe({
    session$sendCustomMessage("condition_lock", list(locked = !(has_mode() && has_conditions_ready())))
  })

  # Phase E.2 step 6 + Q7: the track picker (both Available and Figure
  # tracks) locks whenever mode+conditions+region aren't all ready -- covers
  # the baseline "nothing chosen yet" case and Q7's mid-composition case
  # (switching to Coordinates with an already-placed gene-required track
  # locks the existing Figure-tracks pool too, since it's driven by the same
  # single flag rather than two separate lock layers). Conditions added to
  # the gate in the E.2 follow-up (2026-09-08), alongside the track-pool
  # populator above also now requiring them before showing any tracks.
  observe({
    session$sendCustomMessage("region_lock", list(locked = !(has_mode() && has_conditions_ready() && !is.null(resolved_region()))))
  })

  # Figure-tracks keys that would hard-fail on Generate as things stand --
  # same two conditions the client already red-boxes (Phase E.3's
  # track-gene-required/track-broken CSS classes, applyGeneRequiredFlag()/
  # applyBrokenFlag() above), reimplemented server-side so Generate can gate
  # on it too (E.2 follow-up, 2026-09-08: previously red-boxing was purely a
  # visual warning that didn't stop a doomed render).
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
  # since a disabled button with no explanation is exactly the kind of
  # cryptic dead-end Phase E.4's "no code, no crash" bar is meant to avoid.
  output$generate_warning_ui <- renderUI({
    bad <- invalid_figure_tracks()
    if (length(bad) == 0) return(NULL)
    names_str <- paste(vapply(bad, function(k) TRACK_ITEMS[[k]]$label %||% k, character(1)), collapse = ", ")
    div(class = "alert alert-warning", style = "margin-bottom:0;",
        strong("Remove or fix before generating: "), names_str)
  })

  # Render helper: render_all() straight to an svglite device (called by generate and zoom), in-process, no assemble+source subprocess; the downloadable/displayed R script is a separate, literal concatenation of pipeline/ (assemble_download_script()) built from the same entries/region_params/settings.
  do_render <- function(region_params, entries, settings_base, title_str = "") {
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
      render_all(entries, region_params, settings,
                 device_fn = function(w, h) svglite(settings$output_path, width = w, height = h))
      if (dev.cur() != dev_before) dev.off()
      if (!is.null(state$svg_path) && file.exists(state$svg_path))
        unlink(state$svg_path)
      state$svg_path         <- settings$output_path
      state$rcode             <- assemble_download_script(entries, region_params, settings)
      state$last_region      <- region_params
      state$last_entries     <- entries
      state$last_settings_base <- settings_base
      session$sendCustomMessage("show_output", list())
    }, error = function(e) {
      if (dev.cur() != dev_before) try(dev.off(), silent = TRUE)
      state$plot_error <- conditionMessage(e)
      session$sendCustomMessage("show_output", list())
    })
  }

  # Shared by Compose's and Gene Report's Generate handlers: build_entries() + settings_base + do_render(), from a region/labels/gene-name/title already resolved by the caller (the two modes only differ in how they resolve those).
  render_from_selection <- function(rgn, labels, gene_name, title_str) {
    cfg     <- RENDER_CFG
    col_map <- colour_map_r(labels)
    entries <- build_entries(labels, col_map, cfg)

    region_params <- list(
      chrom = rgn$chrom,
      start = as.integer(rgn$start),
      end   = as.integer(rgn$end),
      gene  = gene_name
    )
    settings_base <- list(
      page_w    = cfg$page_w,
      font      = cfg$font,
      fontsize  = as.integer(cfg$fontsize),
      hic_norm  = cfg$hic_norm,
      zrange    = if (is.null(cfg$zrange_cap) || is.na(cfg$zrange_cap))
                    NULL else c(0, cfg$zrange_cap),
      hic_style = cfg$hic_style,
      hic_res_mode = cfg$hic_res_mode,
      hic_height   = cfg$hic_height,
      hic_max_bins_per_view = cfg$hic_max_bins_per_view,
      atac_height = cfg$atac_height,
      chr_fmt   = cfg$chr_fmt,
      show_track_axes = cfg$show_track_axes,
      txdb_path = cfg$txdb_path,
      gtf_path  = cfg$gtf_path,
      orgdb_pkg           = cfg$orgdb_pkg,
      genome_label        = cfg$genome_label,
      tx_name_lookup_path = cfg$tx_name_lookup_path,
      tx_meta_lookup_path = cfg$tx_meta_lookup_path,
      data_dir = cfg$data_dir,
      samples  = cfg$samples,
      diff_colour_lo  = cfg$diff_colour_lo,
      diff_colour_hi  = cfg$diff_colour_hi,
      junction_opacity_cap = cfg$junction_opacity_cap,
      min_junction_support = cfg$min_junction_support,
      min_display_tpm = cfg$min_display_tpm,
      reference_colour  = cfg$reference_colour,
      comparison_colour = cfg$comparison_colour,
      arc_clip_partial  = cfg$arc_clip_partial
    )

    do_render(region_params, entries, settings_base, title_str)
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
  zoom_to <- function(factor) {
    req(!isTRUE(state$rendering))
    req(!is.null(state$last_region), !is.null(state$last_entries))
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
    do_render(new_rgn, state$last_entries, state$last_settings_base)
  }
  observeEvent(input$zoom_in_15x,  zoom_to(1.5))
  observeEvent(input$zoom_in_3x,   zoom_to(3))
  observeEvent(input$zoom_in_10x,  zoom_to(10))
  observeEvent(input$zoom_out_15x, zoom_to(1 / 1.5))
  observeEvent(input$zoom_out_3x,  zoom_to(1 / 3))
  observeEvent(input$zoom_out_10x, zoom_to(0.1))

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

result <- shiny::runApp(app_obj, ...)
if (!identical(result, "restart") && !identical(result, "edit_settings")) return(invisible(result))

}
}
