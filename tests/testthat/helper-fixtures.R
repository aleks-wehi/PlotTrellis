# Shared bootstrap + fixture table for the golden-image regression suite.
# testthat sources every tests/testthat/helper-*.R file once before running
# any test-*.R file (its own convention) -- nothing in here should be
# re-sourced by an individual test file.
#
# Fixture-path design (solutions/phase-logs/phase-D.3.md's
# "regeneration-friendly" decision, 2026-09-01): every dataset's file list
# lives in ONE entry of the FIXTURES list below, with paths relative to
# that dataset's own data_dir, exactly like project_settings.json's own
# `samples` shape. Adding a new dataset, or a new file to an existing one
# (e.g. the user's incoming additional junction/HiC data), is a one-line
# addition here, never a change to a test-*.R file.
#
# Package version (Phase G.1): render_all()/gene_to_coords()/etc. are already
# in the PlotTrellis namespace via test_check("PlotTrellis") loading the
# package under test -- no source()/APP_DIR bootstrap needed any more.

suppressMessages({
  library(plotgardener)
  library(svglite)
  library(AnnotationDbi)
  library(GenomicFeatures)
  library(GenomicRanges)
  library(rtracklayer)
  library(ggplot2)
  library(org.Mm.eg.db)
  library(org.Hs.eg.db)
})

# Real fixture data (genome caches, ATAC/HiC/junction/quant files) lives in a
# gitignored data/ folder (the package root by default; in the monorepo, its
# root data/ via PT_TEST_DATA_DIR), not shipped with the package
# (far too large, and licensed test data) -- every test_that() below starts
# with skip_if_not(PT_TEST_DATA_AVAILABLE, ...) so `R CMD check` on a machine
# without this data (e.g. R-universe's build farm) reports these as skipped,
# not failed. PT_TEST_DATA_DIR overrides the location (the monorepo keeps the
# fixtures at its root data/, see scripts/sync-plottrellis-to-public.sh).
DATA_ROOT <- Sys.getenv("PT_TEST_DATA_DIR")
if (!nzchar(DATA_ROOT)) DATA_ROOT <- testthat::test_path("..", "..", "data")
PT_TEST_DATA_AVAILABLE <- dir.exists(DATA_ROOT)

# ---- Genome fixtures --------------------------------------------------------
# Real D.2 annotation caches (data/annotations/gencode.vM25*/gencode.v46*).
# ensure_genome_ready() is idempotent (utils/genomes.R) -- the TxDb .sqlite
# already exists from D.2, so only the tx_name/tx_meta lookups get built
# here, once, the first time this helper runs.
if (PT_TEST_DATA_AVAILABLE) {
  GENOMES <- list(
    mouse = ensure_genome_ready(
      gtf_path     = file.path(DATA_ROOT, "annotations", "gencode.vM25.annotation.gtf.gz"),
      orgdb        = "org.Mm.eg.db",
      genome_label = "mm10"
    ),
    human = ensure_genome_ready(
      gtf_path     = file.path(DATA_ROOT, "annotations", "gencode.v46.annotation.gtf.gz"),
      orgdb        = "org.Hs.eg.db",
      genome_label = "hg38"
    )
  )
  TXDBS <- list(
    mouse = AnnotationDbi::loadDb(GENOMES$mouse$txdb_path),
    human = AnnotationDbi::loadDb(GENOMES$human$txdb_path)
  )
} else {
  GENOMES <- NULL
  TXDBS <- NULL
}

# ---- Per-dataset fixture table ----------------------------------------------
# `samples`: RAW shape (condition/colour/data, paths relative to data_dir) --
# identical to project_settings.json's own embedded `samples[]`, and what
# load_samples_from_settings()/gene_display_transcripts() expect. Not every
# condition carries every data type -- only what actually landed (see
# data/mouse_mammary_test/README.md, data/human_encode_test/README.md):
#   - Mouse `loops/LP_*.bedpe`/`ML_*.bedpe` are genuine 2-byte empty source
#     files, so LP/ML have no junction entry here; Ba+P6 (both real BEDPE)
#     cover the mouse-side junction/BEDPE fixtures instead.
#   - Mouse has no per-condition RNA-seq bigWig (D.2 only produced Salmon
#     quant.sf, no coverage track) -- no "rna-seq"-named bw instance for
#     mouse, only "atac".
#   - Human junctions/*.bed are real BED12 splice-junction files (STAR/
#     TopHat convention) -- the human side exercises read_junction_file()'s
#     BED12 path, mouse's BEDPE loops exercise the other.
FIXTURES <- if (!PT_TEST_DATA_AVAILABLE) NULL else list(
  mouse = list(
    data_dir = file.path(DATA_ROOT, "mouse_mammary_test"),
    gene     = "Elf5",
    genome   = "mouse",
    samples = list(
      list(condition = "Ba", colour = "#E41A1C", data = list(
        list(type = "bw",       name = "atac",      path = "atac/Ba_ATAC.bw"),
        list(type = "hic",      name = "hic",       path = "hic/Ba.hic"),
        list(type = "junction", name = "Junctions", path = "loops/Ba_OmniC_FitHiChIP_10kb_FDR01.bedpe"),
        list(type = "quantdir", name = "Quant",     path = "quant/Basal_rep1")
      )),
      list(condition = "LP", colour = "#377EB8", data = list(
        list(type = "bw",       name = "atac",  path = "atac/LP_ATAC.bw"),
        list(type = "hic",      name = "hic",   path = "hic/LP.hic"),
        list(type = "quantdir", name = "Quant", path = "quant/LP_rep1")
      )),
      list(condition = "P6", colour = "#4DAF4A", data = list(
        list(type = "bw",       name = "atac",      path = "atac/P6_ATAC.bw"),
        list(type = "hic",      name = "hic",       path = "hic/P6.hic"),
        list(type = "junction", name = "Junctions", path = "loops/P6_OmniC_FitHiChIP_10kb_FDR01.bedpe")
      )),
      list(condition = "P7", colour = "#984EA3", data = list(
        list(type = "bw",       name = "atac",      path = "atac/P7_ATAC.bw"),
        list(type = "hic",      name = "hic",       path = "hic/P7.hic"),
        list(type = "junction", name = "Junctions", path = "loops/P7_OmniC_FitHiChIP_10kb_FDR01.bedpe")
      ))
    )
  ),
  human = list(
    data_dir = file.path(DATA_ROOT, "human_encode_test"),
    gene     = "PTPRC",
    genome   = "human",
    samples = list(
      list(condition = "GM12878", colour = "#E41A1C", data = list(
        list(type = "bw",       name = "atac",       path = "atac/GM12878_ATAC_foldchange.bigWig"),
        list(type = "bw",       name = "h3k27ac",    path = "atac/GM12878_H3K27ac_foldchange.bigWig"),
        list(type = "bed",      name = "peaks",      path = "atac/GM12878_ATAC_IDRpeaks.bed.gz"),
        list(type = "hic",      name = "hic",        path = "hic/GM12878_CTCF_ChIAPET_rep1.hic"),
        list(type = "junction", name = "Junctions",  path = "junctions/GM12878_rep1.bed"),
        list(type = "quantdir", name = "Quant",      path = "quant/GM12878_rep1")
      )),
      list(condition = "K562", colour = "#377EB8", data = list(
        list(type = "bw",       name = "atac",       path = "atac/K562_ATAC_foldchange.bigWig"),
        list(type = "bw",       name = "h3k27ac",    path = "atac/K562_H3K27ac_foldchange.bigWig"),
        list(type = "bed",      name = "peaks",      path = "atac/K562_ATAC_IDRpeaks.bed.gz"),
        list(type = "hic",      name = "hic",        path = "hic/K562_CTCF_ChIAPET_rep1.hic"),
        list(type = "junction", name = "Junctions",  path = "junctions/K562_rep1.bed"),
        list(type = "quantdir", name = "Quant",      path = "quant/K562_rep1")
      ))
    )
  )
)

# sample_row list (display-based, absolute paths) for a dataset -- what the
# entry builders below actually consume. Built via the app's own
# load_samples_from_settings(), so this fixture harness resolves paths
# exactly the way a real project_settings.json-backed render does (including
# the type="junction" -> "paired-end" rename load_samples_from_settings()
# applies).
test_rows <- function(dataset) {
  fx <- FIXTURES[[dataset]]
  rows <- load_samples_from_settings(fx$samples, fx$data_dir)
  names(rows) <- vapply(rows, `[[`, character(1), "display")
  rows
}

# region_params for a dataset's default fixture gene, or an explicit
# gene/flank override. Uses the real, shared gene_to_coords()
# (utils/gene_lookup.R) directly -- BUGS.md #17's gene_id-versioning bug
# that used to break this against the vM25/v46 fixture TxDbs is fixed at
# the source (utils/genomes.R's ensure_genome_ready() now strips GENCODE
# version suffixes from a freshly-built TxDb's gene table), so no local
# workaround is needed here any more.
test_region <- function(dataset, gene = NULL, flank_pct = 0.15) {
  fx   <- FIXTURES[[dataset]]
  gene <- if (is.null(gene)) fx$gene else gene
  coords <- gene_to_coords(gene, TXDBS[[fx$genome]], flank_pct = flank_pct, orgdb_pkg = GENOMES[[fx$genome]]$orgdb)
  if (is.null(coords)) stop("test_region(): no coords for gene '", gene, "' in dataset '", dataset, "'")
  list(chrom = coords$chrom, start = coords$start, end = coords$end, gene = gene)
}

# Minimal settings list for a dataset, with any field overridable (e.g. a
# non-default atac_height/hic_style/hic_res_mode fixture). settings$samples/
# data_dir are the RAW shape (see FIXTURES above) -- required by
# tpm_usage.R/transcripts.R/transcripts_annotated.R's gene_display_transcripts()
# calls, which do their own load_samples_from_settings() internally.
test_settings <- function(dataset, overrides = list()) {
  fx <- FIXTURES[[dataset]]
  gm <- GENOMES[[fx$genome]]
  base <- list(
    txdb_path = gm$txdb_path, orgdb_pkg = gm$orgdb, genome_label = gm$genome_label,
    tx_name_lookup_path = gm$tx_name_lookup_path, tx_meta_lookup_path = gm$tx_meta_lookup_path,
    gtf_path = gm$gtf_path, data_dir = fx$data_dir, samples = fx$samples,
    title = "", output_path = tempfile(fileext = ".svg"),
    page_w = 9, font = "DejaVu Sans", fontsize = 8L, chr_fmt = "dot",
    hic_norm = "KR", hic_max_bins_per_view = 500, zrange = NULL,
    atac_height = 0.85, hic_style = "rectangle", hic_res_mode = "auto",
    diff_colour_lo = "#0033CC", diff_colour_hi = "#CC0000",
    junction_opacity_cap = 2, arc_clip_partial = FALSE, show_track_axes = TRUE
  )
  modifyList(base, overrides)
}

# ---- render_all() entry constructors ----------------------------------------
# app.R's build_entries() job, reimplemented here since this harness
# deliberately doesn't source app.R (Shiny-reactive code): converts a
# sample_row (test_rows()) into the entry shape render_all() consumes.
#
# entry$viz vs data-instance `type`: a single underlying FILE FORMAT ("bw")
# can back more than one VISUALISATION (a plain bw track, and a bw-vs-bw
# comparison track). `viz` is the only field render_all()'s dispatch tables
# read; for a plain single-file track it's also the data-instance type
# (bw/hic/bed), except sashimi/junctions_mirrored_diff*, which read
# "paired-end" instances under a viz name that's never literally that.
build_entry <- function(type, sample_row) {
  if (type %in% c("genes", "transcripts_all", "transcripts", "transcripts_annotated", "protein_domains")) {
    return(list(viz = type, height = NULL, data = NULL))
  }
  if (type == "tpm_usage") {
    return(list(
      viz = "tpm_usage", sample = sample_row$display, colour = sample_row$colour %||% "#1B9E77",
      quant_path = sample_data_path(sample_row, "quantdir"),
      junction_path = sample_data_path(sample_row, "paired-end", name = "Junctions"),
      height = NULL, data = NULL
    ))
  }
  data_type <- if (type == "sashimi") "paired-end" else type
  inst <- Find(function(d) d$type == data_type, sample_row$data)
  if (is.null(inst)) stop("build_entry(): no '", data_type, "' instance for sample '", sample_row$display, "'")
  list(
    viz = type, sample = sample_row$display, colour = sample_row$colour %||% "#1B9E77",
    name = inst$name, path = inst$path, height = NULL, data = NULL
  )
}

build_comparison_entry <- function(type, control_row, treatment_row) {
  data_type <- if (type == "bw_diff_stacked") "bw" else "paired-end"
  ctrl_inst  <- Find(function(d) d$type == data_type, control_row$data)
  treat_inst <- Find(function(d) d$type == data_type, treatment_row$data)
  if (is.null(ctrl_inst) || is.null(treat_inst)) {
    stop("build_comparison_entry(): '", data_type, "' instance missing on control or treatment sample")
  }
  entry <- list(
    viz = type, name = ctrl_inst$name,
    control_label = control_row$display, treatment_label = treatment_row$display,
    control_path = ctrl_inst$path, treatment_path = treat_inst$path,
    control_colour = control_row$colour %||% "#4D4D4D", treatment_colour = treatment_row$colour %||% "#56B4E9",
    height = NULL, data = NULL
  )
  if (type == "junctions_mirrored_diff_tpm_scaled") {
    entry$control_quant   <- sample_data_path(control_row, "quantdir")
    entry$treatment_quant <- sample_data_path(treatment_row, "quantdir")
  }
  entry
}

# One entry for a single-sample viz. `settings` only matters for "hic"
# (production's build_entries() attaches settings$hic_style onto the entry
# itself -- settings_hic()/render_hic() never read settings$hic_style
# directly, see pipeline/viz/hic.R). An explicit `height` sets entry$height
# up front, same override precedence render_all()'s settings_<viz>()
# functions already give a non-NULL entry$height.
single_track <- function(type, sample_row, extra = list(), height = NULL, settings = NULL) {
  entry <- build_entry(type, sample_row)
  if (type == "hic") entry$hic_style <- settings$hic_style %||% "rectangle"
  if (!is.null(height)) entry$height <- height
  modifyList(entry, extra)
}

pair_track <- function(type, control_row, treatment_row, extra = list(), height = NULL) {
  entry <- build_comparison_entry(type, control_row, treatment_row)
  if (!is.null(height)) entry$height <- height
  modifyList(entry, extra)
}

# vdiffr-ready zero-arg closure for one figure: render_all() draws to
# whichever device is already open when called, same convention vdiffr's
# writer/comparer already expects.
#
# Phase D.4 step 8: optional perf_id/module (both default NULL, so every
# pre-existing call site is unaffected). When perf_id is supplied, the
# returned closure times its own draw (system.time() around the real
# render_all() call) and records it via record_perf() -- which also
# asserts against _perf_baseline.csv's recorded value for this exact
# perf_id, if one exists yet, right there in whichever test_that() block
# calls fig(). Scoped 1:1 with named vdiffr baselines only (perf_id == that
# fixture's vdiffr title) -- a render used only for a numeric assertion (no
# doppelganger) intentionally passes no perf_id.
test_fig <- function(track_specs, region_params, settings, perf_id = NULL, module = NULL) {
  raw_fig <- function() render_all(track_specs, region_params, settings)
  if (is.null(perf_id)) return(raw_fig)
  function() {
    elapsed <- system.time(val <- raw_fig())[["elapsed"]]
    record_perf(perf_id, module, elapsed)
    invisible(val)
  }
}

# ---- Performance-baseline instrumentation (Phase D.4 step 8) ---------------
# _perf_baseline.csv (columns: fixture, module, seconds) lives alongside
# _snaps/ and is committed the same way -- written once, the first time this
# suite runs clean with every correct-render fixture instrumented (see
# test-zzz-perf-baseline.R, which runs last alphabetically among test-*.R
# files and does the actual writing). PERF_BASELINE is loaded once, here, at
# helper-source time -- NULL on the very first run (nothing to assert
# against yet), the real data.frame on every run after that. test_path()
# resolves relative to tests/testthat/ under both devtools::test() and
# R CMD check, so no package-root path is needed here.
PERF_BASELINE_PATH <- testthat::test_path("_perf_baseline.csv")
PERF_BASELINE <- if (file.exists(PERF_BASELINE_PATH)) {
  utils::read.csv(PERF_BASELINE_PATH, stringsAsFactors = FALSE)
} else {
  NULL
}
PERF_TIMES <- new.env()
PERF_TIMES$records <- list()

# Per Phase D Q3 (solutions/concepts/11-phase-dry-run-decisions.md): the
# regression ceiling is 2x the recorded baseline seconds, not an absolute
# ceiling or any other multiplier. Asserted here (inside the helper) rather
# than at each call site -- testthat attributes an expectation raised from
# inside a helper function to whichever test_that() block is currently
# executing, the same convention every other assertion helper in this file
# already relies on.
record_perf <- function(fixture, module, seconds) {
  PERF_TIMES$records[[length(PERF_TIMES$records) + 1]] <- list(
    fixture = fixture, module = module, seconds = seconds
  )
  if (!is.null(PERF_BASELINE)) {
    row <- PERF_BASELINE[PERF_BASELINE$fixture == fixture, ]
    if (nrow(row) == 1) {
      testthat::expect_true(
        seconds <= 2 * row$seconds[1],
        info = sprintf(
          "fixture '%s' took %.2fs, more than 2x its %.2fs baseline",
          fixture, seconds, row$seconds[1]
        )
      )
    }
  }
}
