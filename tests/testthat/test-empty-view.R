# What a track draws when it has nothing to show (W59), for every track type with a data file:
#   - a view with nothing in it (the file reads fine): the track as usual, heading and axes
#     included, nothing plotted, no message;
#   - a file that is missing or can't be read: the red "<key> (<sample>) track failed: <reason>";
#   - a chromosome the file lacks: "Chromosome not found".
# "No data" is drawn by no track. Synthetic files in tempdir() (no fixture data), read off the
# SVG text. Hi-C's empty view needs a real .hic file: test-soft-fail.R covers it.

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

# The view: chr1:1,000-10,000. Every "far" file has its data at chr1:50,000-51,000 only.
VIEW <- list(chrom = "chr1", start = 1000L, end = 10000L, gene = NA_character_)

view_settings <- function() list(
  txdb_path = synth_txdb_path(), orgdb_pkg = "org.Mm.eg.db", genome_label = "synthetic",
  page_w = 9, font = "DejaVu Sans", fontsize = 8L, chr_fmt = "dot", show_track_axes = TRUE
)

synth_file <- function(ext, lines = character(0)) {
  f <- tempfile("w59-", fileext = ext)
  writeLines(lines, f)
  if (!length(lines)) file.create(f)  # zero bytes
  f
}
bw_file <- function(chrom = "chr1") {
  f <- tempfile("w59-", fileext = ".bw")
  gr <- GenomicRanges::GRanges(chrom, IRanges::IRanges(50001, 51000), score = 5)
  GenomeInfoDb::seqlengths(gr) <- setNames(100000, chrom)
  rtracklayer::export(gr, f)
  f
}
far_bed       <- function() synth_file(".bed", "chr1\t50000\t51000\tpeak\t0\t+")
far_junctions <- function() synth_file(".bed", "chr1\t49990\t51010\tj\t5\t+\t49990\t51010\t0\t2\t10,10\t0,1010")
missing_file  <- function(ext) file.path(tempdir(), paste0("w59-not-there", ext))

sample_track <- function(viz, name, path) list(viz = viz, sample = "S1", colour = "#1B9E77", name = name,
                                               path = path, height = NULL, data = NULL)
pair_of <- function(viz, name, ctrl, treat) list(viz = viz, name = name, control_label = "A", treatment_label = "B",
                                                 control_path = ctrl, treatment_path = treat,
                                                 control_colour = "#4D4D4D", treatment_colour = "#56B4E9",
                                                 height = NULL, data = NULL)

render_view <- function(items, region = VIEW) {
  f <- tempfile(fileext = ".svg")
  on.exit(unlink(f), add = TRUE)
  suppressWarnings(suppressMessages(utils::capture.output(
    render_all(items, region, view_settings(), device_fn = function(w, h) svglite::svglite(f, width = w, height = h)))))
  grDevices::dev.off()
  paste(readLines(f), collapse = "\n")
}
has_text <- function(svg, label) grepl(paste0(">", label, "</text>"), svg, fixed = TRUE)
# The grey line plotPairsArches() draws under its arcs (.pt_arc_baseline() when there are none).
n_baselines <- function(svg) lengths(regmatches(svg, gregexpr("<line[^>]*stroke: #BEBEBE", svg)))

skip_no_synth <- function() {
  # R-universe Windows checks render these tracks without the SVG text these tests read.
  skip_on_os("windows")
  skip_if_not_installed("org.Mm.eg.db")
  skip_if_not_installed("txdbmaker")
}

test_that("an empty view draws each track as usual: heading and axes, no message (no fixture data)", {
  skip_no_synth()
  svg <- render_view(list(
    sample_track("bw", "bw_track", bw_file()),
    sample_track("bed", "bed_track", far_bed()),
    sample_track("bed", "empty_bed_file", synth_file(".bed")),
    sample_track("junctions", "junction_track", far_junctions()),
    sample_track("junctions", "empty_junction_file", synth_file(".bed")),
    pair_of("bw_diff_stacked", "bw_pair", bw_file(), bw_file()),
    pair_of("junctions_mirrored_diff", "junction_pair", far_junctions(), far_junctions())
  ))
  expect_false(grepl("No data", svg, fixed = TRUE))
  expect_false(grepl("track failed", svg, fixed = TRUE))
  expect_false(grepl("Chromosome not found", svg, fixed = TRUE))
  for (heading in c("bw_track", "bed_track", "empty_bed_file", "junction_track", "empty_junction_file",
                    "bw_pair diff (stacked)", "junction_pair mirrored diff"))
    expect_true(has_text(svg, heading), label = heading)
  # The junction tracks keep their usual heading detail, at 0 reads.
  expect_identical(lengths(regmatches(svg, gregexpr(">\\s*· S1 \\(max_h=0 reads\\)<", svg))), 2L)
  # Value axes on the 0-1 scale (bigWig, stacked diff), and no arcs, only the baseline of each
  # junction track (two single, one mirrored).
  expect_true(has_text(svg, "0.8"))
  expect_false(grepl("<polygon", svg, fixed = TRUE))
  expect_identical(n_baselines(svg), 3L)
})

test_that("an empty sample keeps the shared range of the others (no fixture data)", {
  skip_no_synth()
  near <- bw_file()
  # Signal in view for one sample: the other, empty, draws the same axis (0 to 5 * 10/9).
  gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(4001, 5000), score = 5)
  GenomeInfoDb::seqlengths(gr) <- c(chr1 = 100000)
  rtracklayer::export(gr, near)
  n_text <- function(svg, label) lengths(regmatches(svg, gregexpr(paste0(">", label, "</text>"), svg, fixed = TRUE)))
  one  <- render_view(list(sample_track("bw", "atac", near)))
  both <- render_view(list(sample_track("bw", "atac", near), sample_track("bw", "atac", bw_file())))
  expect_identical(n_text(one, "5"), 1L)
  expect_identical(n_text(both, "5"), 2L)  # the empty track's axis reaches 5 too
  expect_identical(n_text(both, "0.8"), 0L)  # not the 0-1 scale of a view empty everywhere
})

test_that("a missing or unreadable data file fails its track in red; the rest still draws (no fixture data)", {
  skip_no_synth()
  garbage <- synth_file(".bw", "this is not a bigWig file")
  bad_junctions <- synth_file(".bed", "chr1\t100\t200")  # 3 columns: neither BED12 nor BEDPE
  svg <- render_view(list(
    sample_track("bw", "ok_track", bw_file()),
    sample_track("bw", "atac", missing_file(".bw")),
    sample_track("bw", "atac", garbage),
    sample_track("bed", "peaks", missing_file(".bed")),
    sample_track("hic", "hic", missing_file(".hic")),
    sample_track("junctions", "junctions", missing_file(".junc.bed")),
    sample_track("junctions", "junctions", bad_junctions),
    pair_of("bw_diff_stacked", "atac", bw_file(), missing_file(".bw")),
    pair_of("junctions_mirrored_diff", "junctions", far_junctions(), missing_file(".junc.bed"))
  ))
  failed <- regmatches(svg, gregexpr(">[^<]*track failed: [^<]*<", svg))[[1]]
  expect_length(failed, 8)
  expect_true(any(grepl("^>bw \\(S1\\) track failed: file not found: w59-not-there\\.bw \\(in ", failed)))
  expect_true(any(grepl(paste0("^>bw \\(S1\\) track failed: could not read ", basename(garbage), ": ", basename(garbage),
                               " is not a big wig file<"), failed)))
  expect_true(any(grepl("^>bed \\(S1\\) track failed: file not found: w59-not-there\\.bed", failed)))
  expect_true(any(grepl("^>hic \\(S1\\) track failed: file not found: w59-not-there\\.hic", failed)))
  expect_true(any(grepl("^>junctions \\(S1\\) track failed: file not found: w59-not-there\\.junc\\.bed", failed)))
  expect_true(any(grepl(paste0("^>junctions \\(S1\\) track failed: could not read ", basename(bad_junctions),
                               ": Unrecognised junction/arc file"), failed)))
  expect_identical(sum(grepl("^>bw_diff_stacked track failed: file not found", failed)), 1L)
  expect_identical(sum(grepl("^>junctions_mirrored_diff track failed: file not found", failed)), 1L)
  # In red, and the healthy track still drew.
  expect_identical(lengths(regmatches(svg, gregexpr("fill: #CC0000;[^>]*>[^<]*track failed", svg))), 8L)
  expect_true(has_text(svg, "ok_track"))
  expect_false(grepl("No data", svg, fixed = TRUE))
})

test_that("an unreadable BED file fails its track (no fixture data)", {
  skip_no_synth()
  skip_on_cran()
  bed <- far_bed()
  Sys.chmod(bed, "000")
  on.exit(Sys.chmod(bed, "644"), add = TRUE)
  skip_if(file.access(bed, 4) == 0, "this user can read a file with no read permission (root?)")
  svg <- render_view(list(sample_track("bed", "peaks", bed)))
  expect_match(svg, paste0("bed (S1) track failed: could not read ", basename(bed), ": "), fixed = TRUE)
})

test_that("a chromosome the file lacks stays as text (no fixture data)", {
  skip_no_synth()
  svg <- render_view(list(
    sample_track("bw", "atac", bw_file("chr2")),
    sample_track("bed", "peaks", synth_file(".bed", "chr2\t50000\t51000\tpeak\t0\t+")),
    pair_of("bw_diff_stacked", "atac", bw_file("chr2"), bw_file())
  ))
  expect_identical(lengths(regmatches(svg, gregexpr(">Chromosome not found</text>", svg, fixed = TRUE))), 3L)
  expect_false(grepl("track failed", svg, fixed = TRUE))
})

test_that("gene-only tracks without a gene say so, not \"No data\" (no fixture data)", {
  skip_no_synth()
  tpm <- list(viz = "tpm_usage", sample = "S1", colour = "#1B9E77", name = "quant", quant_path = tempdir(),
              height = NULL, data = NULL)
  scaled <- utils::modifyList(pair_of("junctions_mirrored_diff_tpm_scaled", "junctions", far_junctions(), far_junctions()),
                              list(control_quant = tempdir(), treatment_quant = tempdir()))
  svg <- render_view(list(tpm, scaled))
  expect_identical(lengths(regmatches(svg, gregexpr(">No gene selected for this region</text>", svg, fixed = TRUE))), 2L)
  expect_false(grepl("No data", svg, fixed = TRUE))
  expect_false(grepl("track failed", svg, fixed = TRUE))
})

test_that("a failed first track keeps the figure's scale bar (no fixture data)", {
  skip_no_synth()
  # ">1Kb<" is the scale bar's label and the genome axis's start: twice with the scale bar.
  n_1kb <- function(svg) lengths(regmatches(svg, gregexpr(">1Kb</text>", svg, fixed = TRUE)))
  expect_identical(n_1kb(render_view(list(sample_track("bw", "atac", bw_file())))), 2L)
  failed_in_settings <- render_view(list(sample_track("bw", "atac", missing_file(".bw"))))
  expect_match(failed_in_settings, "track failed", fixed = TRUE)
  expect_identical(n_1kb(failed_in_settings), 2L)
  failed_in_render <- render_view(list(pair_of("bw_diff_stacked", "atac", bw_file(), missing_file(".bw"))))
  expect_match(failed_in_render, "track failed", fixed = TRUE)
  expect_identical(n_1kb(failed_in_render), 2L)
})

test_that("a BED file with only header and comment lines is an empty track (no fixture data)", {
  skip_no_synth()
  svg <- render_view(list(sample_track("bed", "header_only", synth_file(".bed", c("track name=peaks", "# no peaks")))))
  expect_true(has_text(svg, "header_only"))
  expect_false(grepl("Chromosome not found|track failed|No data", svg))
})

test_that("odd data paths fail with a plain reason (no fixture data)", {
  expect_error(.pt_check_file(NA_character_), "^no data file given$")
  expect_error(.pt_check_file(NULL), "^no data file given$")
  expect_error(.pt_check_file(tempdir()), "is a folder, not a data file")
  expect_error(.pt_check_file(missing_file(".bw")), "^file not found: w59-not-there\\.bw \\(in ")
  expect_silent(.pt_check_file(synth_file(".bed")))
})

test_that("a missing quant file is named in the error (no fixture data)", {
  local_mocked_bindings(gtf_transcript_structure = function(...) data.frame(transcript_name = character(0)))
  expect_error(filtered_transcripts("Gene1", list(list(sample_id = "S1", quant_path = missing_file(".sf"))), NULL),
               "^file not found: w59-not-there\\.sf")
})

test_that("a sample's missing quant.sf fails Transcript expression with the sample and file named (bundled example)", {
  skip_if_not_installed("org.Mm.eg.db")
  skip_on_cran()  # builds the example's TxDb
  txt <- plottrellis_example()
  file.rename(file.path(dirname(txt), "quant", "LP_rep1", "quant.sf"), file.path(dirname(txt), "quant", "LP_rep1", "moved.sf"))
  proj <- suppressWarnings(suppressMessages(pt_project(tempfile("w59-project-"), txt)))
  svg <- tempfile(fileext = ".svg")
  suppressWarnings(suppressMessages(pt_save(pt_figure(proj, "Elf5", c("transcript_counts@Basal", "genes")), svg)))
  text <- paste(readLines(svg), collapse = "\n")
  # Basal's track fails on LP's file (the transcripts shown are picked from every sample), so the
  # reason names the sample; the folder's end survives a long temp dir.
  expect_match(text, "tpm_usage (Basal) track failed: sample LP: file not found: quant.sf (in ", fixed = TRUE)
  expect_match(text, "LP_rep1)", fixed = TRUE)
  expect_false(grepl("cannot open the connection", text, fixed = TRUE))
})

# y (pt from the page top) of every point of every <polygon> in an SVG, one vector per polygon.
polygon_ys <- function(svg) {
  pts <- regmatches(svg, gregexpr("<polygon points='[^']*'", svg))[[1]]
  lapply(pts, function(p) as.numeric(sub("^[^,]*,", "", strsplit(sub("^<polygon points='", "", sub("'$", "", p)), " ")[[1]])))
}
# The first track's slot, in pt: render_all() starts it 0.1 in from the top (no title) and it is
# atac_height tall.
first_slot_pt <- function() 72 * c(0.1, 0.1 + .pt_display_defaults()$atac_height)

test_that("signal that is never above 0 is drawn inside its track, not mistaken for an empty view (no fixture data)", {
  skip_no_synth()
  neg <- tempfile("w59-", fileext = ".bw")
  gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(c(2001, 5001), c(3000, 6000)), score = c(-1, -3))
  GenomeInfoDb::seqlengths(gr) <- c(chr1 = 100000)
  rtracklayer::export(gr, neg)
  slot <- first_slot_pt()
  for (items in list(list(sample_track("bw", "atac", neg)), list(pair_of("bw_diff_stacked", "atac", neg, neg)))) {
    svg <- render_view(items)
    ys <- polygon_ys(svg)
    expect_gt(length(ys), 0)
    # On a 0-1 scale the negative signal would be drawn below the track (clipped away); on its
    # own range every point is inside the slot, and the signal has height.
    for (y in ys) expect_true(all(y >= slot[1] - 0.5 & y <= slot[2] + 0.5))
    expect_gt(max(vapply(ys, function(y) diff(range(y)), 0)), 5)
    expect_false(has_text(svg, "0.8"))       # not on the empty view's 0-1 scale
    expect_false(grepl("track failed", svg, fixed = TRUE))
  }
  # Compare: a view with no signal at all draws no polygon and the 0-1 scale.
  empty <- render_view(list(sample_track("bw", "atac", bw_file())))
  expect_length(polygon_ys(empty), 0)
  expect_true(has_text(empty, "0.8"))
})

test_that("a view of explicit zeros is an empty view: the 0-1 axis is drawn (no fixture data)", {
  skip_no_synth()
  # deepTools bamCoverage writes zero-score bins rather than leaving gaps.
  zeros <- function() {
    f <- tempfile("w59-", fileext = ".bw")
    gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(seq(1, 99901, by = 100), width = 100), score = 0)
    GenomeInfoDb::seqlengths(gr) <- c(chr1 = 100000)
    rtracklayer::export(gr, f)
    f
  }
  n_text <- function(svg, label) lengths(regmatches(svg, gregexpr(paste0(">", label, "</text>"), svg, fixed = TRUE)))
  expect_identical(n_text(render_view(list(sample_track("bw", "atac", zeros()))), "0.8"), 1L)
  # Shared range: a zero sample next to one with no signal in view, both on the 0-1 axis.
  expect_identical(n_text(render_view(list(sample_track("bw", "atac", zeros()), sample_track("bw", "atac", bw_file()))), "0.8"), 2L)
  expect_identical(n_text(render_view(list(pair_of("bw_diff_stacked", "atac", zeros(), zeros()))), "0.8"), 1L)
})
