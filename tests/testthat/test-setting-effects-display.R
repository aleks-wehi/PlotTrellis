# Per-setting effect tests for the page/display
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
synth_junctions_track <- function() list(viz = "junctions", sample = "S1", colour = "#1B9E77", name = "Junctions",
                                       path = synth()$bedpe, height = NULL, data = NULL)

skip_no_synth <- function() {
  # R-universe Windows checks render these tracks without the SVG text these
  # tests read; they run on Linux and macOS.
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

# A bed track reads the first three columns of any BED file and draws plain boxes: extra
# columns and BED12 blocks make no difference (the vignette's data-type table says so).
# Not covered: gappedPeak-style files with a decimal column, which have a known problem.
test_that("bed draws one plain box per interval from a BED3, BED6 or BED12 file (no fixture data)", {
  skip_no_synth()
  bed_file <- function(lines) { f <- tempfile(fileext = ".bed"); writeLines(lines, f); f }
  spans <- rbind(c(2000, 2500), c(5000, 5800), c(9000, 9500))
  bed3  <- bed_file(sprintf("chr1\t%d\t%d", spans[, 1], spans[, 2]))
  bed6  <- bed_file(sprintf("chr1\t%d\t%d\tp%d\t5\t+", spans[, 1], spans[, 2], 1:3))
  bed12 <- bed_file(sprintf("chr1\t%d\t%d\tp%d\t5\t+\t%d\t%d\t0\t2\t100,100\t0,%d",
                            spans[, 1], spans[, 2], 1:3, spans[, 1], spans[, 2], spans[, 2] - spans[, 1] - 100))
  outside <- bed_file("chr1\t100\t200")
  on.exit(unlink(c(bed3, bed6, bed12, outside)), add = TRUE)
  n_boxes <- function(path) {
    track <- list(viz = "bed", sample = "S1", name = "peaks", colour = "#1B9E77", path = path, height = NULL, data = NULL)
    svg <- suppressWarnings(synth_svg(list(track), synth_settings()))
    rects <- regmatches(svg, gregexpr("<rect[^>]*>", svg))[[1]]
    sum(grepl("fill: #1B9E77", rects, fixed = TRUE))  # the track's own colour: the interval boxes
  }
  expect_identical(n_boxes(bed3), 3L)
  expect_identical(n_boxes(bed6), 3L)
  expect_identical(n_boxes(bed12), 3L)  # blocks are not drawn: still one box per interval
  expect_identical(n_boxes(outside), 0L)
})

# ---- Fonts: each offered font renders, a missing one falls back ------------

# The Font choices offered in the composer's Figure settings, read from the UI
# itself so a font added there is tested without editing this file.
offered_fonts <- function() {
  html <- as.character(.pt_display_settings_ui(list(), "s_")$display)
  select <- regmatches(html, regexpr('(?s)<select [^>]*id="s_font"[^>]*>.*?</select>', html, perl = TRUE))
  sub('value="([^"]*)"', "\\1", regmatches(select, gregexpr('value="[^"]*"', select))[[1]])
}

# The face svglite writes for `font` when the device resolves it on its own,
# outside render_all(): what every figure label should say.
resolved_family <- function(font) {
  f <- tempfile(fileext = ".svg")
  on.exit(unlink(f))
  svglite::svglite(f, width = 1, height = 1)
  grid::grid.text("x", gp = grid::gpar(fontfamily = font))
  grDevices::dev.off()
  sub(".*font-family: ([^;]*);.*", "\\1", grep("<text", readLines(f), value = TRUE)[1])
}

svg_families <- function(svg) sub(".*font-family: ([^;]*);.*", "\\1", svg_texts(svg))
svg_labels <- function(svg) sub(".*>([^<]*)</text>", "\\1", svg_texts(svg))

# Genes, coverage, junction and stacked-diff tracks plus the genome axis,
# scale bar, value axes and title: every kind of text the synthetic data draws.
synth_all_tracks <- function() {
  list(list(viz = "genes", height = NULL, data = NULL), synth_bw_track(), synth_junctions_track(), synth_diff_track())
}

# render_all() draws a failed track as an error line and only says so in a
# message, so a clean render means no warning and no "track failed" message.
expect_clean_render <- function(expr) {
  expect_no_warning(expect_no_message(expr, message = "track failed"))
}

# macOS's png device (Quartz) warns 'no font could be found for family "<font>"'
# for each label it draws in a family the system lacks, then draws it in its
# own fallback face. That is the fallback a missing font should get, and macOS
# lacks DejaVu Sans too, so that one warning, for the font being rendered, is
# muffled; any other warning still fails the render.
muffle_font_fallback <- function(expr, font) {
  withCallingHandlers(expr, warning = function(w) {
    if (identical(conditionMessage(w), sprintf('no font could be found for family "%s"', font)))
      invokeRestart("muffleWarning")
  })
}

# Renders through the same devices as pt_save(): .pdf and .png must open,
# draw and close with no error, warning or failed track, and write a
# non-empty file.
expect_renders_on_pdf_and_png <- function(settings) {
  for (ext in c("pdf", "png")) {
    f <- tempfile(fileext = paste0(".", ext))
    dev <- switch(ext,
      pdf = function(w, h) grDevices::cairo_pdf(f, width = w, height = h),
      png = function(w, h) grDevices::png(f, width = w, height = h, units = "in", res = 72))
    before <- grDevices::dev.cur()
    expect_clean_render(muffle_font_fallback(
      utils::capture.output(render_all(synth_all_tracks(), SYNTH_REGION, settings, device_fn = dev)),
      settings$font
    ))
    if (grDevices::dev.cur() != before) grDevices::dev.off()
    expect_gt(file.size(f), 0)
    unlink(f)
  }
}

test_that("every font offered in Figure settings renders to SVG, PDF and PNG, all labels in it (no fixture data)", {
  skip_no_synth()
  fonts <- offered_fonts()
  expect_true(all(c("DejaVu Sans", "Helvetica", "Arial") %in% fonts))
  for (font in fonts) {
    settings <- synth_settings(list(font = font, title = "Title"))
    svg <- expect_clean_render(synth_svg(synth_all_tracks(), settings))
    fams <- svg_families(svg)
    expect_gt(length(fams), 10)
    expect_identical(unique(fams), resolved_family(font), label = paste("SVG text families for", font))
    expect_renders_on_pdf_and_png(settings)
  }
})

# A font the system doesn't have (a hand-edited project_settings.json or a
# samples.txt "# font:" line can name any) must not break the figure: every
# device falls back to the face the system uses for an unknown family.
test_that("a font missing from the system falls back to the system's face and still draws every label (no fixture data)", {
  skip_no_synth()
  missing <- "PlotTrellis Missing Font"
  def <- synth_svg(synth_all_tracks(), synth_settings(list(title = "Title")))
  settings <- synth_settings(list(font = missing, title = "Title"))
  svg <- expect_clean_render(synth_svg(synth_all_tracks(), settings))
  fallback <- resolved_family(missing)
  expect_false(identical(fallback, paste0('"', missing, '"')))
  expect_identical(unique(svg_families(svg)), fallback)
  expect_identical(svg_labels(svg), svg_labels(def))
  expect_renders_on_pdf_and_png(settings)
})

# Every track type with fixture data, at a face distinct from the default:
# catches a track that draws some text without ctx$font.
font_cases <- list(
  list(type = "genes", dataset = "mouse", sample = "Ba"),
  list(type = "bw", dataset = "mouse", sample = "Ba"),
  list(type = "hic", dataset = "mouse", sample = "Ba"),
  list(type = "protein_domains", dataset = "mouse", sample = "Ba"),
  list(type = "junctions", dataset = "human", sample = "GM12878"),
  list(type = "tpm_usage", dataset = "mouse", sample = "Ba"),
  list(type = "transcripts", dataset = "mouse", sample = "Ba"),
  list(type = "transcripts_all", dataset = "mouse", sample = "Ba"),
  list(type = "transcripts_annotated", dataset = "mouse", sample = "Ba"),
  list(type = "bed", dataset = "human", sample = "GM12878"),
  list(type = "bw_diff_stacked", dataset = "mouse", control = "Ba", treatment = "LP"),
  list(type = "junctions_mirrored_diff", dataset = "mouse", control = "Ba", treatment = "P6"),
  list(type = "junctions_mirrored_diff_tpm_scaled", dataset = "human", control = "GM12878", treatment = "K562")
)

for (cs in font_cases) {
  local({
    cs <- cs
    test_that(paste0("every ", cs$type, " label follows the font setting (", cs$dataset, " dataset)"), {
      skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
      rows <- test_rows(cs$dataset)
      spec <- if (is.null(cs$sample)) pair_track(cs$type, rows[[cs$control]], rows[[cs$treatment]])
              else single_track(cs$type, rows[[cs$sample]])
      settings <- test_settings(cs$dataset, list(font = "mono", title = "Title"))
      region <- test_region(cs$dataset)
      fams <- svg_families(expect_clean_render(synth_svg(list(spec), settings, region = region)))
      # More text than the title and genome axis alone: the track drew labels of its own.
      expect_gt(length(fams), length(svg_texts(synth_svg(list(), settings, region = region))))
      expect_identical(unique(fams), resolved_family("mono"))
    })
  })
}

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

test_that("show_track_axes toggles the value axis on bw, bw_diff_stacked and junctions (no fixture data)", {
  skip_no_synth()
  for (track in list(synth_bw_track(), synth_diff_track(), synth_junctions_track())) {
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
