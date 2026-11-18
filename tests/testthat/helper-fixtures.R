# Shared bootstrap + fixture table for the golden-image regression suite.
# testthat sources every tests/testthat/helper-*.R file once before running
# any test-*.R file (its own convention) -- nothing in here should be
# re-sourced by an individual test file.
#
# Fixture paths: every dataset's file list
# lives in ONE entry of the FIXTURES list below, with paths relative to
# that dataset's own data_dir, exactly like project_settings.json's own
# `samples` shape. Adding a new dataset, or a new file to an existing one
# (e.g. the user's incoming additional junction/HiC data), is a one-line
# addition here, never a change to a test-*.R file.
#
# render_all()/gene_to_coords()/etc. are already in the PlotTrellis namespace
# via test_check("PlotTrellis") loading the package under test.

suppressMessages({
  library(plotgardener)
  library(svglite)
  library(AnnotationDbi)
  library(GenomicFeatures)
  library(GenomicRanges)
  library(rtracklayer)
})

# Real fixture data (genome caches, ATAC/HiC/junction/quant files) lives in a
# gitignored data/ folder (the package root by default; in the monorepo, its
# root data/ via PT_TEST_DATA_DIR), not shipped with the package
# (far too large, and licensed test data) -- every test_that() that uses it
# starts with skip_if_not(PT_TEST_DATA_AVAILABLE, ...) or
# skip_if_not(PT_TEST_MOUSE_AVAILABLE, ...) (see below), so `R CMD check` on a
# machine without this data (e.g. R-universe's build farm) reports these as
# skipped, not failed. PT_TEST_DATA_DIR overrides the location.
DATA_ROOT <- Sys.getenv("PT_TEST_DATA_DIR")
if (!nzchar(DATA_ROOT)) DATA_ROOT <- testthat::test_path("..", "..", "data")
PT_TEST_DATA_AVAILABLE <- dir.exists(DATA_ROOT)

# Gene-level caches (transcripts shown per gene, UniProt, protein domains) of the fixture
# renders: a fresh folder per test session, so a code change never passes on stale entries.
TEST_CACHE_DIR <- file.path(tempdir(), "pt-test-cache")

# Without the full data, the mouse tests run on fixtures/mini/: the same mouse
# files cut down to the Elf5 region (no Hi-C), committed with the package and
# rebuilt by data-raw/build-test-fixtures.R. Copied to a temporary folder
# first (its genome caches then go in a fresh test-genome-cache there too). Tests that need
# only the mouse data skip on PT_TEST_MOUSE_AVAILABLE; the rest (Hi-C, human,
# genome-wide gene search) keep skipping on PT_TEST_DATA_AVAILABLE.
# Only when NOT_CRAN is "true" (devtools::test(), the node script): a plain
# `R CMD check`, as on R-universe, keeps skipping all of them as before. Source
# tree only: fixtures/ is in .Rbuildignore, like the golden images it is
# checked against, so a check of the built tarball never turns this on.
PT_TEST_MOUSE_AVAILABLE <- PT_TEST_DATA_AVAILABLE ||
  (identical(Sys.getenv("NOT_CRAN"), "true") &&
     dir.exists(testthat::test_path("fixtures", "mini")) &&
     requireNamespace("org.Mm.eg.db", quietly = TRUE) &&
     requireNamespace("vdiffr", quietly = TRUE))
if (!PT_TEST_DATA_AVAILABLE && PT_TEST_MOUSE_AVAILABLE) {
  DATA_ROOT <- tempfile("pt-mini-data-")
  dir.create(DATA_ROOT)
  file.copy(list.files(testthat::test_path("fixtures", "mini"), full.names = TRUE), DATA_ROOT, recursive = TRUE)
}

# ---- Genome fixtures --------------------------------------------------------
# The fixture genomes' caches (TxDb, lookups, chromosome sizes for data/annotations/
# gencode.vM25 and gencode.v46) live in one test project cache, data/test-genome-cache by
# default (PT_TEST_GENOME_CACHE overrides), kept between runs: a full TxDb build takes minutes.
# ensure_genome_ready() is idempotent (genome-main.R) -- whatever caches are missing get built
# there, once. Caches an older version left beside the fixture GTFs seed it (hard-linked
# where possible, so a write to a test genome file would change the original too), so the first run neither rebuilds nor fetches chromosome sizes. With no seeds,
# the first run builds them; run that once before a sharded run.
TEST_GENOME_CACHE <- Sys.getenv("PT_TEST_GENOME_CACHE", file.path(DATA_ROOT, "test-genome-cache"))
test_genome <- function(gtf_name, orgdb, genome_label) {
  gtf_path <- file.path(DATA_ROOT, "annotations", gtf_name)
  for (suffix in c("_txdb.sqlite", "_tx_names.rds", "_tx_meta.rds", "_chrom_sizes.rds")) {
    seed <- file.path(dirname(gtf_path), paste0(sub("\\.gtf(\\.gz)?$", "", gtf_name), suffix))
    dest <- .genome_cache_path(gtf_path, suffix, TEST_GENOME_CACHE)
    if (file.exists(seed) && !file.exists(dest)) {
      # Hard-linked (no extra disk; copied where links can't be made), then renamed into
      # place, so concurrent test shards never see half a file.
      tmp <- tempfile(tmpdir = dirname(dest))
      if (suppressWarnings(file.link(seed, tmp)) || file.copy(seed, tmp)) file.rename(tmp, dest) else unlink(tmp)
    }
  }
  ensure_genome_ready(gtf_path = gtf_path, orgdb = orgdb, genome_label = genome_label,
                      cache_dir = TEST_GENOME_CACHE)
}
# Gives a new project on a fixture genome its genome caches without a rebuild: links the
# test genome cache's files into the project's own cache folder, where the app looks.
seed_project_genome <- function(project_dir, genome = GENOMES$mouse) {
  suffixes <- c(txdb_path = "_txdb.sqlite", tx_name_lookup_path = "_tx_names.rds",
                tx_meta_lookup_path = "_tx_meta.rds", chrom_sizes_path = "_chrom_sizes.rds")
  for (k in names(suffixes)) {
    dest <- .genome_cache_path(genome$gtf_path, suffixes[[k]], file.path(project_dir, "cache"))
    src <- normalizePath(genome[[k]])  # a relative link target would dangle
    if (!file.exists(dest) && !file.symlink(src, dest)) file.copy(src, dest)
  }
  invisible(project_dir)
}
if (PT_TEST_DATA_AVAILABLE) {
  GENOMES <- list(
    mouse = test_genome("gencode.vM25.annotation.gtf.gz", "org.Mm.eg.db", "mm10"),
    human = test_genome("gencode.v46.annotation.gtf.gz", "org.Hs.eg.db", "hg38")
  )
  TXDBS <- list(
    mouse = AnnotationDbi::loadDb(GENOMES$mouse$txdb_path),
    human = AnnotationDbi::loadDb(GENOMES$human$txdb_path)
  )
} else if (PT_TEST_MOUSE_AVAILABLE) {
  # The mini fixtures: txdbmaker's two notices for any GENCODE GTF, as in test-example-data.R.
  GENOMES <- list(mouse = withCallingHandlers(
    test_genome("gencode.vM25.annotation.gtf.gz", "org.Mm.eg.db", "mm10"),
    warning = function(w) {
      if (grepl("\"phase\" metadata column|genome version information", conditionMessage(w)))
        invokeRestart("muffleWarning")
    }))
  TXDBS <- list(mouse = AnnotationDbi::loadDb(GENOMES$mouse$txdb_path))
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
#   - Mouse has no per-condition RNA-seq bigWig (only Salmon quant.sf, no
#     coverage track) -- no "rna-seq"-named bw instance for
#     mouse, only "atac".
#   - Human junctions/*.bed are real BED12 splice-junction files (STAR/
#     TopHat convention) -- the human side exercises read_junction_file()'s
#     BED12 path, mouse's BEDPE loops exercise the other.
FIXTURES <- if (!PT_TEST_MOUSE_AVAILABLE) NULL else list(
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
# The mini fixtures have no human data and no Hi-C: the mouse rows keep their
# hic entries (paths that don't exist there), so only tests that skip on
# PT_TEST_DATA_AVAILABLE may use them.
if (!PT_TEST_DATA_AVAILABLE) FIXTURES$human <- NULL

# Whether a test drawing `type` from `dataset` can run here: the mini fixtures
# cover the mouse dataset except Hi-C. For loops over datasets and track types.
fixture_available <- function(dataset, type) {
  PT_TEST_DATA_AVAILABLE || (PT_TEST_MOUSE_AVAILABLE && dataset == "mouse" && type != "hic")
}

# The protein_domains track fetches UniProt on a cold cache. The node has
# network; a machine running on the mini fixtures may not, so there it skips
# instead of failing its golden image. Checked once per run, 5 s timeout.
.uniprot_reachable <- local({
  ok <- NULL
  function() {
    if (is.null(ok)) {
      old <- options(timeout = 5)
      on.exit(options(old))
      ok <<- tryCatch({
        con <- url("https://rest.uniprot.org/uniprotkb/search?query=Elf5&size=1&fields=accession")
        open(con, "rb")
        close(con)
        TRUE
      }, error = function(e) FALSE, warning = function(w) FALSE)
    }
    ok
  }
})
skip_if_mini_offline <- function() {
  if (!PT_TEST_DATA_AVAILABLE && !.uniprot_reachable())
    testthat::skip("UniProt not reachable (the protein_domains track fetches it)")
}

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
# (genome-gene-lookup.R) directly; it matches the fixture TxDbs' gene ids because
# ensure_genome_ready() strips GENCODE version suffixes from a freshly built
# TxDb's gene table.
test_region <- function(dataset, gene = NULL, flank_pct = 0.15) {
  fx   <- FIXTURES[[dataset]]
  gene <- if (is.null(gene)) fx$gene else gene
  coords <- gene_to_coords(gene, TXDBS[[fx$genome]], flank_pct = flank_pct, orgdb_pkg = GENOMES[[fx$genome]]$orgdb)
  if (is.null(coords)) stop("test_region(): no coords for gene '", gene, "' in dataset '", dataset, "'")
  list(chrom = coords$chrom, start = coords$start, end = coords$end, gene = gene)
}

# Minimal settings list for a dataset, with any field overridable (e.g. a
# non-default atac_height/hic_style/hic_res_mode fixture). settings$samples/
# data_dir are the RAW shape (see FIXTURES above) -- required by the
# transcript tracks' gene_display_transcripts() calls, which do their own
# load_samples_from_settings() internally.
test_settings <- function(dataset, overrides = list()) {
  fx <- FIXTURES[[dataset]]
  gm <- GENOMES[[fx$genome]]
  base <- list(
    txdb_path = gm$txdb_path, orgdb_pkg = gm$orgdb, genome_label = gm$genome_label,
    tx_name_lookup_path = gm$tx_name_lookup_path, tx_meta_lookup_path = gm$tx_meta_lookup_path,
    gtf_path = gm$gtf_path, data_dir = fx$data_dir, samples = fx$samples,
    title = "", output_path = tempfile(fileext = ".svg"),
    page_w = 9, font = "DejaVu Sans", fontsize = 8L, chr_fmt = "dot",
    hic_norm = "KR", hic_max_bins_per_view = 500, zrange_cap = NULL,
    atac_height = 0.85, hic_style = "rectangle", hic_res_mode = "auto",
    diff_colour_lo = "#0033CC", diff_colour_hi = "#CC0000",
    junction_opacity_cap = 2, arc_clip_partial = FALSE, show_track_axes = TRUE,
    cache_dir = TEST_CACHE_DIR
  )
  modifyList(base, overrides)
}

# ---- render_all() entry constructors ----------------------------------------
# What .pt_build_entries() (api-figure-helpers.R) does for the app, here from a
# sample_row (test_rows()) instead of the app's track items: converts it
# into the entry shape render_all() consumes.
#
# entry$viz vs data-instance `type`: a single underlying FILE FORMAT ("bw")
# can back more than one VISUALISATION (a plain bw track, and a bw-vs-bw
# comparison track). `viz` is the only field render_all() dispatches on; for a plain single-file track it's also the data-instance type
# (bw/hic/bed), except junctions/junctions_mirrored_diff*, which read
# "paired-end" instances under a viz name that's never literally that. The data
# type and the kind come from the module's own info_<key>(); a type that is not a
# declared track of the kind asked for is an error.
build_entry <- function(type, sample_row) {
  info <- .pt_viz_info(type)
  if (identical(info$kind, "fixed")) {
    return(list(viz = type, height = NULL, data = NULL))
  }
  if (type == "tpm_usage") {
    return(list(
      viz = "tpm_usage", sample = sample_row$display, colour = sample_row$colour %||% "#1B9E77",
      quant_path = sample_data_path(sample_row, "quantdir"),
      height = NULL, data = NULL
    ))
  }
  if (!identical(info$kind, "sample")) stop("build_entry(): '", type, "' is not a per-sample or sample-free track type")
  data_type <- info$data_type
  inst <- Find(function(d) d$type == data_type, sample_row$data)
  if (is.null(inst)) stop("build_entry(): no '", data_type, "' instance for sample '", sample_row$display, "'")
  list(
    viz = type, sample = sample_row$display, colour = sample_row$colour %||% "#1B9E77",
    name = inst$name, path = inst$path, height = NULL, data = NULL
  )
}

build_comparison_entry <- function(type, control_row, treatment_row) {
  info <- .pt_viz_info(type)
  if (!identical(info$kind, "comparison")) stop("build_comparison_entry(): '", type, "' is not a comparison track type")
  data_type <- info$data_type
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
# (.pt_build_entries() copies settings$hic_style onto the entry itself --
# settings_hic()/render_hic() never read settings$hic_style directly,
# see viz-hic.R). An explicit `height` sets entry$height
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
# Optional perf_id/module (both default NULL). When perf_id is supplied, the
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
    # Registered before rendering, so test-zzz-perf-baseline.R can tell a
    # render that crashed (started, never recorded) from one that finished.
    PERF_TIMES$started <- union(PERF_TIMES$started, perf_id)
    elapsed <- system.time(val <- raw_fig())[["elapsed"]]
    record_perf(perf_id, module, elapsed)
    invisible(val)
  }
}

# ---- Performance-baseline instrumentation ---------------------------------
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
PERF_TIMES$started <- character(0)

# The regression ceiling is a multiple of the recorded baseline seconds, not an
# absolute ceiling: 2x (a loaded login node can exceed it; a compute node
# doesn't). Asserted here (inside the helper) rather
# than at each call site -- testthat attributes an expectation raised from
# inside a helper function to whichever test_that() block is currently
# executing, the same convention every other assertion helper in this file
# already relies on.
record_perf <- function(fixture, module, seconds) {
  PERF_TIMES$records[[length(PERF_TIMES$records) + 1]] <- list(
    fixture = fixture, module = module, seconds = seconds
  )
  # Timed against the full data on the compute node only: the mini fixtures
  # draw less and run on any machine, where a 2x ceiling means nothing.
  if (!is.null(PERF_BASELINE) && PT_TEST_DATA_AVAILABLE) {
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
