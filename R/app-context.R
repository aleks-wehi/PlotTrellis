# The composer's startup (R/app-main.R): everything .pt_composer_app() loads once per
# launch (genome migration, TxDb, chromosome lengths, samples, track items, gene
# list, render settings, start view), returned as one list that the UI builder
# and the server blocks read. The body is unindented, as it was inside the app
# closure.
.pt_composer_context <- function(SETTINGS_JSON, CFG) {

# ── Startup: load heavy objects once ──────────────────────────────────────────
# The genome is the flat fields gtf_path, orgdb and genome_label on CFG. An
# older project_settings.json with a genome_id and no orgdb is migrated: orgdb
# and genome_label come from GENOME_PRESETS by that id and are written back, so
# this runs once. An id that isn't a preset is an error (never a silent
# default): fill in the Genome section by hand.
if (is.null(CFG$orgdb) && !is.null(CFG$genome_id)) {
  preset <- GENOME_PRESETS[[CFG$genome_id]]
  if (is.null(preset))
    stop("project_settings.json has an old-format genome_id ('", CFG$genome_id,
         "') that isn't a known built-in preset, so it can't be auto-migrated. ",
         "Re-open Settings and fill in the Genome section (GTF path, OrgDb, genome label) by hand.",
         call. = FALSE)
  CFG$genome_id    <- NULL
  CFG$orgdb        <- preset$orgdb
  CFG$genome_label <- preset$genome_label
  writeLines(jsonlite::toJSON(CFG, auto_unbox = TRUE, pretty = TRUE,
                              null = "null", digits = NA), SETTINGS_JSON)
  message("Migrated project_settings.json's old genome_id ('", preset$label,
          "') to flat genome fields.")
}

# Cache paths are never read from or written to CFG: every launch re-derives the
# genome's cache paths from gtf_path (ensure_genome_ready()). ORGDB_PKG in the
# list below is a local of this app; the figure and lookup functions get the
# genome's packages and paths as arguments (RENDER_CFG, the render settings).
# Caches go to the project's cache/ folder, beside its project_settings.json
# (.genome_cache_path()).
PROJECT_DIR <- dirname(SETTINGS_JSON)
CACHE_DIR <- file.path(PROJECT_DIR, "cache")
.pt_require_genome_label(CFG$genome_label, SETTINGS_JSON)
GENOME <- ensure_genome_ready(gtf_path = CFG$gtf_path, orgdb = CFG$orgdb,
                              genome_label = CFG$genome_label,
                              cache_dir = CACHE_DIR)
ORGDB_PKG           <- GENOME$orgdb

message("Loading PlotTrellis...")
# The TxDb and its CDS index come from .pt_annotation()'s session memo, so every
# render and the Gene info panel share this one copy instead of opening their
# own.
ANNOT    <- .pt_annotation(GENOME)
txdb     <- ANNOT$txdb

# Built once here (not per gene render) since GenomicFeatures::cdsBy() walks
# the whole TxDb regardless of which genes are actually asked for. The Gene
# info table's UniProt canonical-isoform check (genome-uniprot-canonical.R) and the
# protein-domains track (through the same memo) look up a handful of
# transcripts' translated CDS lengths per gene out of this.
CDS_BY_TX <- ANNOT$cds_by_tx()

# Chromosome lengths: a GTF carries none, so they are never read from the TxDb's
# seqinfo. They come from UCSC through GENOME$chrom_sizes_path
# (build_chrom_sizes(), built and cached by ensure_genome_ready() above), never
# from a hardcoded table: one species' lengths would be silently wrong under any
# other genome.
CHR_LENGTHS <- readRDS(GENOME$chrom_sizes_path)

SAMPLES_DF  <- load_samples_from_settings(CFG$samples, CFG$data_dir)

# Every sample is selectable in every view_mode (One/Two/Multi), whatever data
# it has. Which tracks a sample, or a pair of samples, offers is decided in
# viz-track-detection.R (detect_generic(), detect_pair_generic() and the modules'
# own detect_<key>()), which give no item rather than an error when a data type
# is missing.
AVAILABLE_CONDITIONS <- vapply(SAMPLES_DF, function(r) r$display, character(1))

# Every sample, with no pre-filter: detect_pair_generic() checks the data on
# both sides of each pair, so e.g. two ATAC-only samples still get an ATAC
# coverage comparison.
TRACK_ITEMS <- c(
  build_track_items(SAMPLES_DF),
  build_comparison_track_items(SAMPLES_DF)
)

ALL_GENES <- sort(AnnotationDbi::keys(get(ORGDB_PKG,
                                          envir = asNamespace(ORGDB_PKG)),
                                      keytype = "SYMBOL"))

# The one cfg every render and lookup reads, computed once from CFG: the flat
# CFG keys plus the genome-derived values (orgdb_pkg, genome_label,
# tx_*_lookup_path). The wizard edits the core settings (with a restart), and
# the composer's Figure settings panel layers live values for the other keys on
# top of it (saved_cfg()/live_cfg() in server()).
RENDER_CFG <- .pt_render_cfg(CFG, GENOME, CACHE_DIR)

# Each sample's own colour, else PALETTE in sample order.
DEFAULT_SAMPLE_COLS <- .pt_sample_colours(SAMPLES_DF)

# Optional preloaded starting view (pt_preload()): the UI starts with this mode,
# condition(s), gene and figure tracks already picked, so Generate is the only
# click. NULL (no start_view in the settings) keeps the usual blank start.
START_VIEW <- if (is.list(CFG$start_view))
  lapply(CFG$start_view, unlist) else NULL
# Keys that no longer exist (a sample renamed since) are dropped, not shown,
# after renaming the junction track's old keys. Order is kept, and so are
# repeats of Gene annotation, the one track a figure can hold several times; any
# other track is kept once (a pt_figure() can repeat one, the picker can't).
if (!is.null(START_VIEW)) {
  sv_keys <- as.character(.pt_renamed_junctions_key(START_VIEW$tracks))
  START_VIEW$tracks <- sv_keys[sv_keys %in% names(TRACK_ITEMS) &
                                 (sv_keys == "genes" | !duplicated(sv_keys))]
}

# Multi's Available pool, sorted by label (.pt_compose_available_keys()).
COMPOSE_AVAILABLE_KEYS <- .pt_compose_available_keys(TRACK_ITEMS)

# Gene-required track types (needs_gene in their info_<key>()): from coordinates
# they can only draw an empty-state message. Kept as a plain R vector too (not
# just the JSON string the page carries) so server() can run the same red-box
# check the client already does, to gate Generate on it.
GENE_REQUIRED_TYPES <- .pt_gene_required_types()

list(SETTINGS_JSON = SETTINGS_JSON, CFG = CFG, PROJECT_DIR = PROJECT_DIR,
     CACHE_DIR = CACHE_DIR,
     GENOME = GENOME, ORGDB_PKG = ORGDB_PKG, txdb = txdb, CDS_BY_TX = CDS_BY_TX,
     CHR_LENGTHS = CHR_LENGTHS, SAMPLES_DF = SAMPLES_DF,
     AVAILABLE_CONDITIONS = AVAILABLE_CONDITIONS,
     TRACK_ITEMS = TRACK_ITEMS, ALL_GENES = ALL_GENES, RENDER_CFG = RENDER_CFG,
     DEFAULT_SAMPLE_COLS = DEFAULT_SAMPLE_COLS, START_VIEW = START_VIEW,
     COMPOSE_AVAILABLE_KEYS = COMPOSE_AVAILABLE_KEYS,
     GENE_REQUIRED_TYPES = GENE_REQUIRED_TYPES)
}
