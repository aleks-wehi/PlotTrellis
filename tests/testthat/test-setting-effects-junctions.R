# Phase I.3 step 28 (stage B): per-setting effect tests for the junction/arc
# settings (min_junction_support, arc_clip_partial, and junction_opacity_cap
# and show_track_axes on the mirrored-diff tracks). Each test renders the
# same figure at two values of one setting and asserts a specific difference.
#
# Every test marked "no fixture data" uses a synthetic TxDb + junction files
# in tempdir(), so it also runs on R-universe. The TPM-scaled test needs the
# human fixture (quant.sf + junctions) and skips without it.

synth_txdb_path <- local({
  path <- NULL
  function() {
    if (is.null(path)) {
      txdb <- suppressWarnings(txdbmaker::makeTxDb(
        transcripts = data.frame(tx_id = 1L, tx_name = "ENSMUST00000000001", tx_chrom = "chr1",
                                 tx_strand = "+", tx_start = 2000L, tx_end = 8000L),
        splicings = data.frame(tx_id = 1L, exon_rank = 1:2, exon_start = c(2000L, 7000L), exon_end = c(3000L, 8000L)),
        genes = data.frame(tx_id = 1L, gene_id = "ENSMUSG00000000001"),
        chrominfo = data.frame(chrom = "chr1", length = 100000L, is_circular = FALSE)
      ))
      path <<- tempfile("synth-txdb-", fileext = ".sqlite")
      AnnotationDbi::saveDb(txdb, path)
    }
    path
  }
})

SYNTH_REGION <- list(chrom = "chr1", start = 1000L, end = 10000L, gene = NA_character_)

synth_settings <- function(overrides = list()) {
  modifyList(list(
    txdb_path = synth_txdb_path(), orgdb_pkg = "org.Mm.eg.db", genome_label = "synthetic",
    page_w = 9, font = "DejaVu Sans", fontsize = 8L, chr_fmt = "dot", show_track_axes = TRUE
  ), overrides)
}

# BED12 junctions: one line per (intron_start, intron_end, reads); each has two
# 10bp anchor blocks either side of the intron.
bed12_file <- function(junctions) {
  f <- tempfile("synth-", fileext = ".bed")
  writeLines(vapply(junctions, function(j) {
    s <- j[1] - 10L; e <- j[2] + 10L
    sprintf("chr1\t%d\t%d\tj\t%d\t+\t%d\t%d\t0\t2\t10,10\t0,%d", s, e, j[3], s, e, e - s - 10L)
  }, character(1)), f)
  f
}

render_svg <- function(items, settings, region) {
  f <- tempfile(fileext = ".svg")
  on.exit(unlink(f), add = TRUE)
  dev <- NULL
  on.exit(if (!is.null(dev) && dev %in% grDevices::dev.list()) grDevices::dev.off(dev), add = TRUE, after = FALSE)
  utils::capture.output(render_all(items, region, settings, device_fn = function(w, h) {
    svglite::svglite(f, width = w, height = h)
    dev <<- grDevices::dev.cur()
  }))
  grDevices::dev.off(dev)
  dev <- NULL
  paste(readLines(f), collapse = "\n")
}

# plotPairsArches() draws each arc as one <polygon>, and nothing else in these figures is one.
n_arcs <- function(svg) lengths(regmatches(svg, gregexpr("<polygon", svg, fixed = TRUE)))
has_text <- function(svg, label) grepl(paste0(">", label, "</text>"), svg, fixed = TRUE)

skip_no_synth <- function() {
  # R-universe Windows checks render these tracks without the SVG text these
  # tests read (2026-10-27); they run on Linux and macOS.
  skip_on_os("windows")
  skip_if_not_installed("org.Mm.eg.db")
  skip_if_not_installed("txdbmaker")
}

sashimi_track <- function(path) list(viz = "sashimi", sample = "S1", colour = "#1B9E77", name = "Junctions",
                                     path = path, height = NULL, data = NULL)

test_that("min_junction_support filters sashimi arcs (no fixture data)", {
  skip_no_synth()
  track <- sashimi_track(bed12_file(list(c(3000L, 7000L, 12L), c(4000L, 6000L, 2L), c(4500L, 5500L, 30L))))
  svg_at <- function(ms) render_svg(list(track), synth_settings(list(min_junction_support = ms)), SYNTH_REGION)
  expect_identical(n_arcs(svg_at(0)), 3L)
  expect_identical(n_arcs(svg_at(10)), 2L)
  above_all <- svg_at(31)
  expect_identical(n_arcs(above_all), 0L)
  expect_true(has_text(above_all, "No data"))
})

test_that("min_junction_support filters junctions_mirrored_diff arcs (no fixture data)", {
  skip_no_synth()
  track <- list(viz = "junctions_mirrored_diff", name = "Junctions", control_label = "A", treatment_label = "B",
                control_path   = bed12_file(list(c(3000L, 7000L, 12L), c(4000L, 6000L, 2L))),
                treatment_path = bed12_file(list(c(3000L, 7000L, 30L), c(4000L, 6000L, 5L))),
                control_colour = "#4D4D4D", treatment_colour = "#56B4E9", height = NULL, data = NULL)
  svg_at <- function(ms) render_svg(list(track), synth_settings(list(min_junction_support = ms)), SYNTH_REGION)
  expect_identical(n_arcs(svg_at(0)), 4L)
  expect_identical(n_arcs(svg_at(10)), 2L)
})

test_that("arc_clip_partial drops sashimi arcs that overhang the window (no fixture data)", {
  skip_no_synth()
  # Intron 510-5000 starts left of the 1000-10000 window: the partial arc.
  track <- sashimi_track(bed12_file(list(c(510L, 5000L, 30L), c(3000L, 7000L, 12L), c(4000L, 6000L, 2L))))
  svg_at <- function(clip) render_svg(list(track), synth_settings(list(arc_clip_partial = clip)), SYNTH_REGION)
  expect_identical(n_arcs(svg_at(FALSE)), 3L)
  expect_identical(n_arcs(svg_at(TRUE)), 2L)
})

test_that("arc_clip_partial = FALSE keeps a partial arc from a BEDPE file too (no fixture data)", {
  skip("BUGS.md: read_junction_file()'s BEDPE branch drops any pair whose intron overhangs the window, so arc_clip_partial has no effect on BEDPE input, see step 28 stage B report")
  skip_no_synth()
  bedpe <- tempfile("synth-", fileext = ".bedpe")
  writeLines(c("chr1\t500\t510\tchr1\t5000\t5010\tj0\t30",
               "chr1\t3000\t3010\tchr1\t7000\t7010\tj1\t12"), bedpe)
  svg_at <- function(clip) render_svg(list(sashimi_track(bedpe)), synth_settings(list(arc_clip_partial = clip)), SYNTH_REGION)
  expect_identical(n_arcs(svg_at(FALSE)), 2L)
  expect_identical(n_arcs(svg_at(TRUE)), 1L)
})

# Value-axis tick labels are the only text drawn left of the track column
# (track_x = 0.8in = 57.6pt).
left_of_track_texts <- function(svg) {
  t <- regmatches(svg, gregexpr("<text[^>]*>[^<]*</text>", svg))[[1]]
  t[as.numeric(sub("^<text x='([0-9.]+)'.*", "\\1", t)) < 0.8 * 72]
}

test_that("show_track_axes toggles the value axis on junctions_mirrored_diff (no fixture data)", {
  skip_no_synth()
  track <- list(viz = "junctions_mirrored_diff", name = "Junctions", control_label = "A", treatment_label = "B",
                control_path   = bed12_file(list(c(3000L, 7000L, 12L))),
                treatment_path = bed12_file(list(c(3000L, 7000L, 30L))),
                control_colour = "#4D4D4D", treatment_colour = "#56B4E9", height = NULL, data = NULL)
  on  <- render_svg(list(track), synth_settings(), SYNTH_REGION)
  off <- render_svg(list(track), synth_settings(list(show_track_axes = FALSE)), SYNTH_REGION)
  expect_gte(length(left_of_track_texts(on)), 3)
  expect_length(left_of_track_texts(off), 0)
})

# One test for both TPM-scaled settings: settings_<viz>() (the slow quant.sf
# pass, which reads neither setting) runs once, then render_all() reuses its
# entry via a mock, so each extra render costs only the drawing.
test_that("junction_opacity_cap and show_track_axes change the TPM-scaled track (human dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  rows          <- test_rows("human")
  region_params <- test_region("human")
  spec          <- pair_track("junctions_mirrored_diff_tpm_scaled", rows[["GM12878"]], rows[["K562"]])
  base          <- test_settings("human")
  prepared <- settings_junctions_mirrored_diff_tpm_scaled(list(), spec, base, base$page_w - 1.8, region_params)
  local_mocked_bindings(settings_junctions_mirrored_diff_tpm_scaled = function(dict, entry, ...) prepared)
  svg_at <- function(...) render_svg(list(spec), modifyList(base, list(...)), region_params)
  opacities <- function(svg) {
    polys <- regmatches(svg, gregexpr("<polygon[^>]*>", svg))[[1]]
    op <- sub(".*fill-opacity: ([0-9.]+);.*", "\\1", polys)
    as.numeric(ifelse(grepl("fill-opacity", polys, fixed = TRUE), op, "1"))
  }
  def <- svg_at(junction_opacity_cap = 2)
  low <- svg_at(junction_opacity_cap = 0.5)
  expect_true(has_text(def, "2.00 log2FC"))
  expect_true(has_text(low, "0.50 log2FC"))
  op_def <- opacities(def)
  op_low <- opacities(low)
  expect_gt(length(op_def), 0)
  expect_identical(length(op_low), length(op_def))
  # A lower cap saturates sooner: every arc at least as opaque, some strictly more.
  expect_true(all(op_low >= op_def))
  expect_true(any(op_low > op_def))

  expect_gte(length(left_of_track_texts(def)), 3)
  expect_length(left_of_track_texts(svg_at(show_track_axes = FALSE)), 0)
})
