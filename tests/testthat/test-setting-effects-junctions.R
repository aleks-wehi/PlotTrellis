# Per-setting effect tests for the junction/arc
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
  # tests read; they run on Linux and macOS.
  skip_on_os("windows")
  skip_if_not_installed("org.Mm.eg.db")
  skip_if_not_installed("txdbmaker")
}

junctions_track <- function(path) list(viz = "junctions", sample = "S1", colour = "#1B9E77", name = "Junctions",
                                     path = path, height = NULL, data = NULL)

test_that("min_junction_support filters junction arcs (no fixture data)", {
  skip_no_synth()
  track <- junctions_track(bed12_file(list(c(3000L, 7000L, 12L), c(4000L, 6000L, 2L), c(4500L, 5500L, 30L))))
  svg_at <- function(ms) render_svg(list(track), synth_settings(list(min_junction_support = ms)), SYNTH_REGION)
  expect_identical(n_arcs(svg_at(0)), 3L)
  expect_identical(n_arcs(svg_at(10)), 2L)
  above_all <- svg_at(31)
  expect_identical(n_arcs(above_all), 0L)
  # All filtered out: the track as usual with no arcs, not a message.
  expect_true(has_text(above_all, " \u00B7 S1 (max_h=0 reads)"))
  expect_false(grepl("No data", above_all, fixed = TRUE))
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

test_that("arc_clip_partial drops junction arcs that overhang the window (no fixture data)", {
  skip_no_synth()
  # Intron 510-5000 starts left of the 1000-10000 window: the partial arc.
  track <- junctions_track(bed12_file(list(c(510L, 5000L, 30L), c(3000L, 7000L, 12L), c(4000L, 6000L, 2L))))
  svg_at <- function(clip) render_svg(list(track), synth_settings(list(arc_clip_partial = clip)), SYNTH_REGION)
  expect_identical(n_arcs(svg_at(FALSE)), 3L)
  expect_identical(n_arcs(svg_at(TRUE)), 2L)
})

test_that("arc_clip_partial = FALSE keeps a partial arc from a BEDPE file too (no fixture data)", {
  skip_no_synth()
  bedpe <- tempfile("synth-", fileext = ".bedpe")
  # j0 starts left of the 1000-10000 window (the partial arc); j2 lies wholly left of it.
  writeLines(c("chr1\t500\t510\tchr1\t5000\t5010\tj0\t30",
               "chr1\t3000\t3010\tchr1\t7000\t7010\tj1\t12",
               "chr1\t100\t110\tchr1\t300\t310\tj2\t8"), bedpe)
  svg_at <- function(clip) render_svg(list(junctions_track(bedpe)), synth_settings(list(arc_clip_partial = clip)), SYNTH_REGION)
  expect_identical(n_arcs(svg_at(FALSE)), 2L)
  expect_identical(n_arcs(svg_at(TRUE)), 1L)
})

test_that("read_junction_file() keeps a BEDPE arc with an end in the window, drops one with none", {
  bedpe <- tempfile("synth-", fileext = ".bedpe")
  on.exit(unlink(bedpe))
  writeLines(c("chr1\t500\t510\tchr1\t5000\t5010\tpartial_left\t30",
               "chr1\t3000\t3010\tchr1\t12000\t12010\tpartial_right\t12",
               "chr1\t4000\t4010\tchr1\t6000\t6010\tinside\t5",
               "chr1\t400\t410\tchr1\t20000\t20010\tspans_window\t500",
               "chr1\t100\t110\tchr1\t300\t310\toutside\t8",
               "chr1\t5000\t5010\tchr1\t500\t510\tswapped_anchors\t7",
               "chr1\t100\t110\tchr1\t1000\t1010\ton_left_edge\t9",     # intron ends exactly at the window start
               "chr1\t9990\t10000\tchr1\t15000\t15010\ton_right_edge\t11",  # intron starts exactly at its end
               "chr1\t990\t1000\tchr1\t20000\t20010\tstarts_on_left_edge\t13"), bedpe)  # drawn: runs from the window start
  got <- read_junction_file(bedpe, "chr1", 1000L, 10000L)
  # Not 500 (spans), 8 (outside), 9 or 11 (touch the edge from outside, so nothing is drawn).
  expect_identical(sort(got$nSupport), c(5L, 7L, 12L, 13L, 30L))
})

# A synthetic HiCCUPS loop file (merged_loops.bedpe / postprocessed_pixels_<res>.bedpe layout:
# 24 tab-separated columns, a column-name header, a juicer_tools line, "." in score) and the
# 10-column BEDPE the HiCCUPS-to-BEDPE conversion makes of it: anchors as they are, score =
# `observed`, comment lines and rows without a numeric start dropped, sorted, named, strands ".".
HICCUPS_HEADER <- paste(c("#chr1", "x1", "x2", "chr2", "y1", "y2", "name", "score", "strand1", "strand2",
                          "color", "observed", "expectedBL", "expectedDonut", "expectedH", "expectedV",
                          "fdrBL", "fdrDonut", "fdrH", "fdrV", "numCollapsed", "centroid1", "centroid2",
                          "radius"), collapse = "\t")
hiccups_row <- function(c1, x1, x2, c2, y1, y2, observed) {
  paste(c(c1, x1, x2, c2, y1, y2, ".", ".", ".", ".", "0,255,255", observed, "27.05", "27.64", "21.74",
          "30.42", "0.0045", "0.0065", "5.3E-5", "0.0071", "4", x2, y1, "7500"), collapse = "\t")
}
hiccups_file <- function(chr = "chr1") {
  f <- tempfile("synth-hiccups-", fileext = ".bedpe")
  writeLines(c(HICCUPS_HEADER,
               "# juicer_tools version 2.20.00",
               hiccups_row(chr, 4000, 4010, chr, 6000, 6010, "5.0"),   # inside the window
               hiccups_row(chr, 500, 510, chr, 5000, 5010, "30.0"),    # partial, left end outside
               hiccups_row(chr, 3000, 3010, "chr2", 7000, 7010, "40.0"), # trans: dropped
               "# a stray comment line",
               hiccups_row(chr, 100, 110, chr, 300, 310, "8.0"),       # wholly left of the window
               hiccups_row(chr, "NA", 3010, chr, 7000, 7010, "99.0"),  # no numeric start: skipped
               hiccups_row(chr, 3000, 3010, chr, 7000, 7010, "12.5")), # inside, not a whole count
             f)
  f
}
hiccups_converted_file <- function() {
  f <- tempfile("synth-converted-", fileext = ".bedpe")
  writeLines(c("chr1\t100\t110\tchr1\t300\t310\tS_merged_loop1\t8.0\t.\t.",
               "chr1\t500\t510\tchr1\t5000\t5010\tS_merged_loop2\t30.0\t.\t.",
               "chr1\t3000\t3010\tchr1\t7000\t7010\tS_merged_loop3\t12.5\t.\t.",
               "chr1\t3000\t3010\tchr2\t7000\t7010\tS_merged_loop4\t40.0\t.\t.",
               "chr1\t4000\t4010\tchr1\t6000\t6010\tS_merged_loop5\t5.0\t.\t."), f)
  f
}
arcs_sorted <- function(d) {
  d <- d[order(d$intron_start, d$intron_end), ]
  rownames(d) <- NULL
  d
}

test_that("read_junction_file() reads a raw HiCCUPS loop file as its converted BEDPE (no fixture data)", {
  raw <- hiccups_file()
  bare <- hiccups_file(chr = "1")  # HiCCUPS run on a genome without "chr" names
  conv <- hiccups_converted_file()
  on.exit(unlink(c(raw, bare, conv)))
  expect_true(is_paired_end_structural(raw))
  want <- arcs_sorted(read_junction_file(conv, "chr1", 1000L, 10000L))
  expect_identical(want$nSupport, c(30, 12.5, 5))  # the observed column: not 8 (outside), 40 (trans), 99 (bad start)
  expect_identical(arcs_sorted(read_junction_file(raw, "chr1", 1000L, 10000L)), want)
  expect_identical(arcs_sorted(read_junction_file(bare, "chr1", 1000L, 10000L)), want)
  expect_identical(nrow(read_junction_file(raw, "chr3", 1000L, 10000L)), 0L)
})

test_that("a raw HiCCUPS loop file draws the same arcs as its converted BEDPE (no fixture data)", {
  skip_no_synth()
  raw <- hiccups_file()
  conv <- hiccups_converted_file()
  on.exit(unlink(c(raw, conv)))
  svg_of <- function(path, ms) render_svg(list(junctions_track(path)), synth_settings(list(min_junction_support = ms)), SYNTH_REGION)
  for (ms in c(0, 10)) expect_identical(svg_of(raw, ms), svg_of(conv, ms))
  expect_identical(n_arcs(svg_of(raw, 0)), 3L)
  expect_identical(n_arcs(svg_of(raw, 10)), 2L)  # the min-support filter acts on `observed`
})

test_that("a BEDPE with a HiCCUPS-style header but no observed column still reads as a BEDPE (no fixture data)", {
  f <- tempfile("synth-bedpe-", fileext = ".bedpe")
  on.exit(unlink(f))
  writeLines(c("#chr1\tx1\tx2\tchr2\ty1\ty2\tname\tscore\tstrand1\tstrand2",
               "chr1\t4000\t4010\tchr1\t6000\t6010\tl1\t5\t.\t."), f)
  expect_true(is_paired_end_structural(f))
  expect_identical(read_junction_file(f, "chr1", 1000L, 10000L)$nSupport, 5L)
})

test_that("HiCCUPS numbers take the types the converted BEDPE reads with (no fixture data)", {
  f <- tempfile("synth-hiccups-", fileext = ".bedpe")
  on.exit(unlink(f))
  writeLines(c(HICCUPS_HEADER, hiccups_row("chr1", 4000, 4010, "chr1", 6000, 6010, "7"),
               hiccups_row("chr1", 4000, 4010, "chr1", 3000000000, 3000000010, "9")), f)
  got <- read_hiccups_loops(f)
  expect_identical(got$start1, c(4000L, 4000L))
  expect_identical(got$start2, c(6000, 3e9))  # past the integer range: double, as read.delim reads it
  expect_identical(got$nSupport, c(7L, 9L))
})

test_that("a HiCCUPS line cut before its second chromosome is dropped, not an all-NA arc (no fixture data)", {
  f <- tempfile("synth-hiccups-", fileext = ".bedpe")
  on.exit(unlink(f))
  writeLines(c(HICCUPS_HEADER, hiccups_row("chr2", 4000, 4010, "chr2", 6000, 6010, "7"), "chr2\t9000\t90"), f)
  expect_identical(nrow(read_junction_file(f, "chr1", 1000L, 10000L)), 0L)
  expect_identical(nrow(pool_junctions(list(list(junction_path = f)), "chr1", 1000L, 10000L, 0)), 0L)
})

# A junction BED6 (chrom, start, end, name or ".", count, strand), the interval being the intron,
# and the BEDPE with the same introns between 10 bp anchors.
bed6_rows <- function(chr) c(
  paste(chr, 3000, 7000, ".", 12, "+", sep = "\t"),  # inside the window
  paste(chr, 4000, 6000, ".", 2, "-", sep = "\t"),   # inside, few reads
  paste(chr, 500, 5000, ".", 30, "*", sep = "\t"),   # partial, left end outside
  paste(chr, 100, 300, ".", 8, "+", sep = "\t"),     # wholly left of the window
  paste("2", 3000, 7000, ".", 40, "+", sep = "\t"))  # another chromosome
bed6_file <- function(chr = "1") {
  f <- tempfile("synth-junctions-", fileext = ".bed")
  writeLines(bed6_rows(chr), f)
  f
}
bed6_as_bedpe_file <- function() {
  f <- tempfile("synth-junctions-", fileext = ".bedpe")
  writeLines(c("chr1\t2990\t3000\tchr1\t7000\t7010\tj1\t12",
               "chr1\t3990\t4000\tchr1\t6000\t6010\tj2\t2",
               "chr1\t490\t500\tchr1\t5000\t5010\tj3\t30",
               "chr1\t90\t100\tchr1\t300\t310\tj4\t8",
               "chr2\t2990\t3000\tchr2\t7000\t7010\tj5\t40"), f)
  f
}

test_that("read_junction_file() reads a junction BED6: the interval is the intron, the score its reads (no fixture data)", {
  bare <- bed6_file("1")
  prefixed <- bed6_file("chr1")
  named_file <- function(name) {
    f <- tempfile("synth-junctions-", fileext = ".bed")
    writeLines(sub("\t\\.\t", paste0("\t", name, "\t"), bed6_rows("1")), f)
    f
  }
  # Names in column 4 instead of ".", including ones that look like chromosomes.
  named <- vapply(c("JUNC1", "1", "chr1:3001-7000:+", "CHRNA1"), named_file, character(1))
  on.exit(unlink(c(bare, prefixed, named)))
  expect_true(is_paired_end_structural(bare))
  want <- data.frame(intron_start = c(500L, 3000L, 4000L), intron_end = c(5000L, 7000L, 6000L),
                     nSupport = c(30L, 12L, 2L))
  for (f in c(bare, prefixed, named)) expect_identical(arcs_sorted(read_junction_file(f, "chr1", 1000L, 10000L)), want, info = f)
})

test_that("a junction BED6 with a strand that is not + - * or . is an error naming it, not an empty track (no fixture data)", {
  dots <- tempfile("synth-junctions-", fileext = ".bed")
  padded <- tempfile("synth-junctions-", fileext = ".bed")
  blank <- tempfile("synth-junctions-", fileext = ".bed")
  on.exit(unlink(c(dots, padded, blank)))
  writeLines(sub("\t[-+*]$", "\t.", bed6_rows("1")), dots)    # every strand "."
  writeLines(sub("$", " ", bed6_rows("1")), padded)             # trailing space after the strand
  rows <- bed6_rows("1")
  rows[2] <- sub("\t[-+*]$", "\t", rows[2])                   # one strand left empty
  writeLines(rows, blank)
  expect_identical(nrow(read_junction_file(dots, "chr1", 1000L, 10000L)), 3L)
  expect_identical(nrow(read_junction_file(padded, "chr1", 1000L, 10000L)), 3L)
  expect_error(read_junction_file(blank, "chr1", 1000L, 10000L), "column 6 must be a strand.*found ''")
})

test_that("a 6-column BEDPE still reads as BEDPE: a number in column 6 (no fixture data)", {
  f <- tempfile("synth-bedpe-", fileext = ".bedpe")
  bare <- tempfile("synth-bedpe-", fileext = ".bedpe")
  trans <- tempfile("synth-bedpe-", fileext = ".bedpe")
  on.exit(unlink(c(f, bare, trans)))
  writeLines(c("chr1\t490\t500\tchr1\t5000\t5010", "chr1\t3000\t3010\tchr2\t7000\t7010"), f)
  writeLines(c("1\t490\t500\t1\t5000\t5010", "1\t3990\t4000\t1\t6000\t6010"), bare)
  writeLines(c("1\t2990\t3000\t2\t7000\t7010", "1\t3990\t4000\tX\t6000\t6010"), trans)  # trans pairs only
  expect_identical(nrow(read_junction_file(trans, "chr1", 1000L, 10000L)), 0L)
  expect_identical(arcs_sorted(read_junction_file(f, "chr1", 1000L, 10000L)),
                   data.frame(intron_start = 500L, intron_end = 5000L, nSupport = 1L))
  expect_identical(arcs_sorted(read_junction_file(bare, "chr1", 1000L, 10000L)),
                   data.frame(intron_start = c(500L, 4000L), intron_end = c(5000L, 6000L), nSupport = c(1L, 1L)))
})

test_that("a junction BED6 draws the same arcs as the BEDPE with the same introns (no fixture data)", {
  skip_no_synth()
  bed6 <- bed6_file()
  bedpe <- bed6_as_bedpe_file()
  on.exit(unlink(c(bed6, bedpe)))
  svg_of <- function(path, ms) render_svg(list(junctions_track(path)), synth_settings(list(min_junction_support = ms)), SYNTH_REGION)
  for (ms in c(0, 10)) expect_identical(svg_of(bed6, ms), svg_of(bedpe, ms))
  expect_identical(n_arcs(svg_of(bed6, 0)), 3L)
  expect_identical(n_arcs(svg_of(bed6, 10)), 2L)  # the min-support filter acts on the score column
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

# The value axis is the vertical line left of the track (track_x - 0.08 in = 51.84pt): it spans
# the arcs, so its length is the height the arcs were given.
axis_spans <- function(svg) {
  lines <- regmatches(svg, gregexpr("<line x1='51.84' [^>]*x2='51.84'[^>]*>", svg))[[1]]
  y <- function(a) as.numeric(sub(sprintf(".*\\b%s='([0-9.]+)'.*", a), "\\1", lines))
  span <- abs(y("y2") - y("y1"))
  list(top = pmin(y("y1"), y("y2"))[span > 5], length = span[span > 5])
}

test_that("junction_height: a single-sample track's arcs are that tall, each half of a mirrored comparison too (no fixture data)", {
  skip_no_synth()
  single <- junctions_track(bed12_file(list(c(3000L, 7000L, 12L))))
  mirrored <- list(viz = "junctions_mirrored_diff", name = "Junctions", control_label = "A", treatment_label = "B",
                   control_path   = bed12_file(list(c(3000L, 7000L, 12L))),
                   treatment_path = bed12_file(list(c(3000L, 7000L, 30L))),
                   control_colour = "#4D4D4D", treatment_colour = "#56B4E9", height = NULL, data = NULL)
  for (jh in list(NULL, 1.4)) {
    svg <- render_svg(list(single, mirrored), synth_settings(list(junction_height = jh)), SYNTH_REGION)
    jh <- jh %||% 0.85
    axes <- axis_spans(svg)
    expect_equal(axes$length, c(jh, 2 * jh) * 72, tolerance = 0.02 / 60)
    # Each track is its heading band (0.3 in) plus its arcs; tracks are 0.1 in apart.
    expect_equal(diff(axes$top), (jh + 0.1 + 0.3) * 72, tolerance = 0.02 / 60)
    # The mirrored track's 0 tick is exactly halfway: one junction_height per sample.
    t <- left_of_track_texts(svg)
    zero_y <- as.numeric(sub(".* y='([0-9.]+)'.*", "\\1", t[grepl(">0</text>$", t)]))
    max_y  <- as.numeric(sub(".* y='([0-9.]+)'.*", "\\1", t[grepl(">30</text>$", t)]))
    expect_length(max_y, 2)
    expect_equal(diff(sort(c(max_y, zero_y[length(zero_y)]))), rep(jh * 72, 2), tolerance = 0.02 / 60)
  }
})

test_that("mirrored_junction_arches(): opacity is the capped |log2 fold change| of all reads; min_support and treat_scale", {
  treat <- data.frame(intron_start = c(100L, 500L), intron_end = c(1000L, 2000L), nSupport = c(3L, 30L))
  ctrl  <- data.frame(intron_start = 100L, intron_end = 1000L, nSupport = 12L)
  a <- mirrored_junction_arches(treat, ctrl, "chr1", min_support = 0, cap = 2)
  expect_identical(a$treat$nSupport, c(3, 30))
  expect_identical(a$treat[c("chrom1", "start1", "end1", "start2", "end2")],
                   data.frame(chrom1 = "chr1", start1 = c(100L, 500L), end1 = c(101L, 501L),
                              start2 = c(1000L, 2000L), end2 = c(1001L, 2001L)))
  expect_equal(a$treat$opacity, c(abs(log2(4 / 13)) / 2, 1))  # 500-2000 is absent from ctrl: capped
  expect_equal(a$ctrl$opacity, abs(log2(4 / 13)) / 2)
  # Junctions under min_support aren't drawn, but still count in the other sample's fold change.
  a <- mirrored_junction_arches(treat, ctrl, "chr1", min_support = 10, cap = 2)
  expect_identical(a$treat$start1, 500L)
  expect_equal(a$ctrl$opacity, abs(log2(4 / 13)) / 2)
  # treat_scale multiplies Comparison's counts in the heights and the fold change.
  a <- mirrored_junction_arches(treat, ctrl, "chr1", min_support = 0, cap = 2, treat_scale = 4)
  expect_identical(a$treat$nSupport, c(12, 120))
  expect_identical(a$ctrl$nSupport, 12)
  expect_equal(a$ctrl$opacity, 0)
  # A sample with no junctions gives zero arcs, not an error.
  a <- mirrored_junction_arches(treat[0, ], ctrl, "chr1", min_support = 0, cap = 2)
  expect_identical(nrow(a$treat), 0L)
  expect_identical(a$treat$opacity, numeric(0))
  expect_equal(a$ctrl$opacity, 1)
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
  prepared <- settings_junctions_mirrored_diff_tpm_scaled(list(), spec, .pt_fill_settings(base), base$page_w - 1.8,
                                                          region_params)
  expect_identical(prepared$entry$height, 0.3 + 2 * 0.85)  # heading band + a junction_height per sample
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

test_that("TPM-scaled junctions: a quant file whose TPM column isn't numbers is an error, not a scale from 0 (no fixture data)", {
  # filtered_transcripts() flags such a file (and counts it as 0 TPM) instead of stopping.
  local_mocked_bindings(
    load_junctions = function(...) NULL,
    .pt_annotation = function(settings) NULL,
    filtered_transcripts = function(gene, samples, anno) structure(data.frame(), total_tpm = 0,
      tpm_not_numeric = if (samples[[1]]$sample_id == "LP") "LP_rep1/quant.sf" else character(0))
  )
  run <- function(ctrl, treat) settings_junctions_mirrored_diff_tpm_scaled(list(),
    list(control_label = ctrl, treatment_label = treat), list(gtf_path = NULL), 7.2, list(gene = "Elf5"))
  msg <- "^the TPM column of a quant file is not all numbers \\(corrupt file\\?\\): LP_rep1/quant.sf$"
  expect_error(run("LP", "Ba"), msg)
  expect_error(run("Ba", "LP"), msg)
  expect_no_error(run("Ba", "Ba"))
})

