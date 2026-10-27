# Phase I.3 step 28 (stage B): per-setting effect tests for the page/display
# settings (font, fontsize, page_w, chr_fmt, show_track_axes, atac_height on
# the bw_diff_stacked path). Each test renders the same figure at two values
# of one setting and asserts a specific difference, so it fails if the
# setting stops changing the output.
#
# NO FIXTURE DATA NEEDED: every test here builds a one-transcript TxDb and
# tiny bigWig/BEDPE files in tempdir(), so these also run on R-universe.
# Only org.Mm.eg.db (Suggests) is required, for render_all()'s assembly().

synth <- local({
  cache <- new.env()
  function() {
    if (is.null(cache$txdb_path)) {
      txdb <- suppressWarnings(txdbmaker::makeTxDb(
        transcripts = data.frame(tx_id = 1L, tx_name = "ENSMUST00000000001", tx_chrom = "chr1",
                                 tx_strand = "+", tx_start = 2000L, tx_end = 8000L),
        splicings = data.frame(tx_id = 1L, exon_rank = 1:2, exon_start = c(2000L, 7000L), exon_end = c(3000L, 8000L)),
        genes = data.frame(tx_id = 1L, gene_id = "ENSMUSG00000000001"),
        chrominfo = data.frame(chrom = "chr1", length = 100000L, is_circular = FALSE)
      ))
      cache$txdb_path <- tempfile("synth-txdb-", fileext = ".sqlite")
      AnnotationDbi::saveDb(txdb, cache$txdb_path)
      bw <- function(phase) {
        pos <- seq(1000L, 9990L, 10L)
        gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(pos, width = 10L),
                                     score = 1 + 9 * abs(sin(pos / 500 + phase)),
                                     seqinfo = GenomeInfoDb::Seqinfo("chr1", 100000L))
        f <- tempfile("synth-", fileext = ".bw")
        rtracklayer::export.bw(gr, f)
        f
      }
      cache$bw_a <- bw(0)
      cache$bw_b <- bw(1)
      cache$bedpe <- tempfile("synth-", fileext = ".bedpe")
      writeLines(c("chr1\t3000\t3010\tchr1\t7000\t7010\tj1\t12",
                   "chr1\t4000\t4010\tchr1\t6000\t6010\tj2\t5"), cache$bedpe)
    }
    as.list(cache)
  }
})

synth_settings <- function(overrides = list()) {
  modifyList(list(
    txdb_path = synth()$txdb_path, orgdb_pkg = "org.Mm.eg.db", genome_label = "synthetic",
    page_w = 9, font = "DejaVu Sans", fontsize = 8L, chr_fmt = "dot", show_track_axes = TRUE,
    atac_height = 0.85, diff_colour_lo = "#0033CC", diff_colour_hi = "#CC0000"
  ), overrides)
}

SYNTH_REGION <- list(chrom = "chr1", start = 1000L, end = 10000L, gene = NA_character_)

synth_svg <- function(items, settings, region = SYNTH_REGION) {
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

svg_texts <- function(svg) regmatches(svg, gregexpr("<text[^>]*>[^<]*</text>", svg))[[1]]
text_of <- function(svg, label) {
  t <- svg_texts(svg)
  t[grepl(paste0(">", label, "</text>"), t, fixed = TRUE)]
}

synth_bw_track <- function() list(viz = "bw", sample = "S1", colour = "#1B9E77", name = "atac",
                                  path = synth()$bw_a, height = NULL, data = NULL)
synth_diff_track <- function() list(viz = "bw_diff_stacked", name = "atac", control_label = "A", treatment_label = "B",
                                    control_path = synth()$bw_a, treatment_path = synth()$bw_b,
                                    control_colour = "#4D4D4D", treatment_colour = "#56B4E9", height = NULL, data = NULL)
synth_sashimi_track <- function() list(viz = "sashimi", sample = "S1", colour = "#1B9E77", name = "Junctions",
                                       path = synth()$bedpe, height = NULL, data = NULL)

skip_no_synth <- function() {
  # R-universe Windows checks render these tracks without the SVG text these
  # tests read (2026-10-27); they run on Linux and macOS.
  skip_on_os("windows")
  skip_if_not_installed("org.Mm.eg.db")
  skip_if_not_installed("txdbmaker")
}

# svglite writes the font the system resolved the family to (e.g. "Courier"
# becomes "Nimbus Mono PS" on Milton), so this compares families rather than
# expecting a literal name. "mono" resolves to a monospace face on any
# fontconfig system, never to the default's sans face.
test_that("font changes the font-family of figure text (no fixture data)", {
  skip_no_synth()
  family <- function(svg) {
    t <- text_of(svg, "chr.1")
    expect_length(t, 1)
    sub(".*font-family: ([^;]*);.*", "\\1", t)
  }
  def <- family(synth_svg(list(), synth_settings()))
  alt <- family(synth_svg(list(), synth_settings(list(font = "mono"))))
  expect_false(identical(def, alt))
})

test_that("fontsize changes text size (no fixture data)", {
  skip_no_synth()
  def <- synth_svg(list(), synth_settings())
  big <- synth_svg(list(), synth_settings(list(fontsize = 14L)))
  expect_match(text_of(def, "chr.1"), "font-size: 8.00px", fixed = TRUE)
  expect_match(text_of(big, "chr.1"), "font-size: 14.00px", fixed = TRUE)
})

test_that("page_w changes the page width (no fixture data)", {
  skip_no_synth()
  width_pt <- function(svg) as.numeric(sub(".*<svg[^>]* width='([0-9.]+)pt'.*", "\\1", svg))
  expect_equal(width_pt(synth_svg(list(), synth_settings())), 9 * 72)
  expect_equal(width_pt(synth_svg(list(), synth_settings(list(page_w = 12)))), 12 * 72)
})

test_that("chr_fmt changes the genome-axis chromosome label (no fixture data)", {
  skip_no_synth()
  labels <- function(fmt) {
    svg <- synth_svg(list(), synth_settings(list(chr_fmt = fmt)))
    c(dot = length(text_of(svg, "chr.1")), space = length(text_of(svg, "chr 1")), plain = length(text_of(svg, "chr1")))
  }
  expect_identical(labels("dot"),   c(dot = 1L, space = 0L, plain = 0L))
  expect_identical(labels("space"), c(dot = 0L, space = 1L, plain = 0L))
  expect_identical(labels("plain"), c(dot = 0L, space = 0L, plain = 1L))
})

# Value-axis tick labels are the only text drawn left of the track column
# (track_x = 0.8in = 57.6pt); the genome axis starts AT track_x.
left_of_track_texts <- function(svg) {
  t <- svg_texts(svg)
  x <- as.numeric(sub("^<text x='([0-9.]+)'.*", "\\1", t))
  t[x < 0.8 * 72]
}

test_that("show_track_axes toggles the value axis on bw, bw_diff_stacked and sashimi (no fixture data)", {
  skip_no_synth()
  for (track in list(synth_bw_track(), synth_diff_track(), synth_sashimi_track())) {
    on  <- synth_svg(list(track), synth_settings())
    off <- synth_svg(list(track), synth_settings(list(show_track_axes = FALSE)))
    expect_gte(length(left_of_track_texts(on)), 2, label = paste(track$viz, "axis labels with show_track_axes = TRUE"))
    expect_length(left_of_track_texts(off), 0)
  }
})

test_that("atac_height sizes the bw_diff_stacked track (no fixture data)", {
  skip_no_synth()
  height_pt <- function(svg) as.numeric(sub(".*<svg[^>]* height='([0-9.]+)pt'.*", "\\1", svg))
  def  <- height_pt(synth_svg(list(synth_diff_track()), synth_settings()))
  tall <- height_pt(synth_svg(list(synth_diff_track()), synth_settings(list(atac_height = 1.6))))
  expect_equal(tall - def, (1.6 - 0.85) * 72, tolerance = 1e-6)
  # Default layout, no title (none given): top 0.1 + track 0.85 + gap 0.1 + genome axis 0.35.
  expect_equal(def, (0.1 + 0.85 + 0.1 + GENOME_AXIS_RESERVE) * 72, tolerance = 1e-6)
})
