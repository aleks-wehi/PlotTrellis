# Per-setting effect tests for the transcript /
# TPM tracks (min_display_tpm through a real transcripts-track render,
# fontsize on the tpm_usage track's text). test-caching.R already
# covers min_display_tpm at the gene_display_transcripts() level;
# the render test here also catches the setting not reaching the track.
# Needs the mouse fixture (Elf5, quant.sf).

render_svg <- function(items, region, settings) {
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
text_label <- function(t) sub(".*>([^<]*)</text>$", "\\1", t)
text_size  <- function(t) as.numeric(sub(".*font-size: ([0-9.]+)px.*", "\\1", t))

test_that("min_display_tpm changes which transcripts the transcripts track draws (mouse dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  rows   <- test_rows("mouse")
  region <- test_region("mouse")
  spec   <- single_track("transcripts", rows[["Ba"]])
  base   <- test_settings("mouse")
  # Fresh cache folder per value, so a stale cache entry can't mask or fake the effect.
  def  <- text_label(svg_texts(render_svg(list(spec), region,
                                          modifyList(base, list(min_display_tpm = 1, cache_dir = tempfile("pt-cache-"))))))
  none <- text_label(svg_texts(render_svg(list(spec), region,
                                          modifyList(base, list(min_display_tpm = 1e9, cache_dir = tempfile("pt-cache-"))))))
  expect_true(any(grepl("^Elf5-", def)))
  expect_false("No retained transcripts for this gene in view" %in% def)
  expect_false(any(grepl("^Elf5-", none)))
  expect_true("No retained transcripts for this gene in view" %in% none)
})

test_that("fontsize scales the tpm_usage track's bar labels and legend (mouse dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  rows   <- test_rows("mouse")
  region <- test_region("mouse")
  spec   <- single_track("tpm_usage", rows[["Ba"]])
  sizes <- function(fontsize) {
    t <- svg_texts(render_svg(list(spec), region, test_settings("mouse", overrides = list(fontsize = fontsize))))
    lab <- text_label(t)
    list(bar = text_size(t[grepl("^[0-9]+%$", lab)]), legend = text_size(t[grepl("^Elf5-", lab)]))
  }
  def <- sizes(8L)
  big <- sizes(16L)
  expect_gt(length(def$bar), 0)
  expect_gt(length(def$legend), 0)
  # Both are one point under the fontsize setting.
  expect_equal(def$bar, rep(7, length(def$bar)), tolerance = 0.01)
  expect_equal(def$legend, rep(7, length(def$legend)), tolerance = 0.01)
  # unique(): bigger text drops a % label that no longer fits its segment.
  expect_equal(unique(big$bar), 15, tolerance = 0.01)
  expect_equal(big$legend, rep(15, length(big$legend)), tolerance = 0.01)
})

# The Transcript expression panel's parts in a rendered SVG (svglite units: 72 per inch).
svg_attr <- function(el, a) as.numeric(sub(sprintf(".*\\b%s='([-0-9.]+)(px)?'.*", a), "\\1", el))
svg_num  <- function(x) as.numeric(gsub(",", "", x))
# Bar segments are the only white-outlined rects; legend keys the only filled rects without an outline.
svg_bars <- function(svg) regmatches(svg, gregexpr("<rect[^>]*stroke: #FFFFFF;[^>]*>", svg))[[1]]
svg_keys <- function(svg) regmatches(svg, gregexpr("<rect x=[^>]*stroke: none; fill: #[^>]*>", svg))[[1]]
svg_lines <- function(svg) regmatches(svg, gregexpr("<line [^>]*>", svg))[[1]]
# The TPM axis as px = intercept + slope * TPM, from its tick labels (centred numbers).
svg_tpm_axis <- function(t) {
  ticks <- t[grepl("^[0-9.,]+$", text_label(t)) & grepl("text-anchor='middle'", t)]
  ticks <- ticks[order(svg_attr(ticks, "x"))]
  list(ticks = ticks, fit = stats::lm(svg_attr(ticks, "x") ~ svg_num(text_label(ticks))))
}
tpm_px <- function(axis, tpm) sum(stats::coef(axis$fit) * c(1, tpm))

# The Transcript expression bars against their own TPM axis.
test_that("tpm_usage bars end at their total on a readable TPM axis, and % labels fit their segment (mouse dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  rows <- test_rows("mouse")
  # Default settings, and a narrow page with bigger text, where the axis is widened past the
  # usual 1/9 to fit the largest total's label and fewer ticks fit.
  for (ov in list(list(), list(page_w = 6, fontsize = 12L))) {
    what <- if (length(ov)) paste(names(ov), unlist(ov), collapse = ", ") else "default"
    svg  <- render_svg(list(single_track("tpm_usage", rows[["Ba"]]), single_track("tpm_usage", rows[["LP"]])),
                       test_region("mouse"), test_settings("mouse", overrides = ov))
    t   <- svg_texts(svg)
    lab <- text_label(t)

    axis <- svg_tpm_axis(t)
    ticks <- axis$ticks
    expect_gte(length(ticks), 2)
    # Neighbouring tick labels don't overlap.
    half <- svg_attr(ticks, "textLength") / 2
    expect_true(all(diff(svg_attr(ticks, "x")) > half[-1] + half[-length(half)]), label = what)

    bars <- svg_bars(svg)
    totals <- lab[grepl("^[0-9.,]+ TPM$", lab)]
    expect_length(totals, 2)
    rows_y <- sort(unique(svg_attr(bars, "y")))
    expect_length(rows_y, 2)
    for (k in seq_along(rows_y)) {
      seg <- bars[svg_attr(bars, "y") == rows_y[k]]
      left  <- min(svg_attr(seg, "x"))
      right <- max(svg_attr(seg, "x") + svg_attr(seg, "width"))
      total <- svg_num(sub(" TPM$", "", totals[k]))
      expect_lt(abs(left - tpm_px(axis, 0)), 0.5, label = what)
      # Totals are printed rounded (integers from 10 up, 2 significant digits below): allow
      # half a TPM (or 1 px).
      expect_lt(abs(right - tpm_px(axis, total)), max(1, 0.5 * stats::coef(axis$fit)[[2]]), label = what)
    }
    if (!length(ov)) default_svg <- svg
  }

  t   <- svg_texts(default_svg)
  lab <- text_label(t)
  bars <- svg_bars(default_svg)

  # Every % label lies inside one segment of its own row.
  pct <- t[grepl("^[0-9]+%$", lab)]
  expect_gt(length(pct), 0)
  for (p in pct) {
    x <- svg_attr(p, "x"); y <- svg_attr(p, "y"); w <- svg_attr(p, "textLength")
    inside <- vapply(bars, function(b) {
      bx <- svg_attr(b, "x"); by <- svg_attr(b, "y")
      bx <= x - w / 2 && x + w / 2 <= bx + svg_attr(b, "width") && by <= y && y <= by + svg_attr(b, "height")
    }, logical(1))
    expect_true(any(inside), info = text_label(p))
  }
})

test_that("tpm_usage: total TPM 0 draws '0 TPM'; a gene not in the quant file or a non-number total is a track error (mouse dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  rows <- test_rows("mouse")
  settings <- test_settings("mouse")
  region <- test_region("mouse")
  elf5_tx <- gtf_transcript_structure("Elf5", .pt_annotation(settings))$tx_id_stripped
  quant <- read.delim(single_track("tpm_usage", rows[["LP"]])$quant_path)
  is_elf5 <- strip_version(quant$Name) %in% elf5_tx
  expect_gt(sum(is_elf5), 0)
  # LP's real quant.sf, edited, as the track's own quant file.
  labels_with_quant <- function(q) {
    path <- tempfile(fileext = ".sf")
    on.exit(unlink(path))
    write.table(q, path, sep = "\t", quote = FALSE, row.names = FALSE)
    text_label(svg_texts(render_svg(list(single_track("tpm_usage", rows[["LP"]], extra = list(quant_path = path))),
                                    region, settings)))
  }
  # An error is drawn in the sample's row, on two lines when long: joined here.
  joined <- function(lab) paste(lab, collapse = " ")

  zero_q <- quant
  zero_q$TPM[is_elf5] <- 0
  lab <- labels_with_quant(zero_q)
  expect_true("0 TPM" %in% lab)
  expect_false(grepl("failed:", joined(lab)))
  expect_true(all(c("Transcript expression", "LP", "Elf5-201", "TPM", "0") %in% lab))  # drawn as usual

  lab <- labels_with_quant(quant[!is_elf5, ])
  expect_match(joined(lab), paste("failed: none of Elf5's transcripts are in this sample's quant file",
                                  "(is the GTF the one the quant files were made with?)"), fixed = TRUE)
  expect_false("0 TPM" %in% lab)
  expect_true(all(c("Transcript expression", "LP") %in% lab))  # the row keeps its heading and label

  # One stray word anywhere reads the TPM column as text: an error, not every transcript at 0 TPM.
  text_q <- quant
  text_q$TPM[which(!is_elf5)[1]] <- "abc"
  lab <- labels_with_quant(text_q)
  expect_match(joined(lab), "failed: the TPM column of this sample's quant file is not all numbers (corrupt file?): ", fixed = TRUE)
  expect_false("0 TPM" %in% lab)
  # ...and only that sample's track: LP's corrupt file must not fail Ba's tracks, which also
  # read every condition's quant file to pick the transcripts shown.
  lp_bad <- tempfile("LP_bad")
  dir.create(lp_bad)
  on.exit(unlink(lp_bad, recursive = TRUE), add = TRUE)
  write.table(text_q, file.path(lp_bad, "quant.sf"), sep = "\t", quote = FALSE, row.names = FALSE)
  bad_settings <- settings
  bad_settings$cache_dir <- tempfile("pt-cache-")
  lp <- which(vapply(bad_settings$samples, function(x) identical(x$condition, "LP"), logical(1)))
  q <- which(vapply(bad_settings$samples[[lp]]$data, function(d) d$type == "quantdir", logical(1)))
  # Quant folders are relative to data_dir.
  up <- strrep("../", length(strsplit(normalizePath(settings$data_dir), "/")[[1]]) - 1)
  bad_settings$samples[[lp]]$data[[q]]$path <- paste0(up, sub("^/", "", normalizePath(lp_bad)))
  expect_true(file.exists(file.path(settings$data_dir, bad_settings$samples[[lp]]$data[[q]]$path, "quant.sf")))
  lab <- text_label(svg_texts(render_svg(list(single_track("tpm_usage", rows[["Ba"]]), single_track("transcripts", rows[["Ba"]])),
                                         region, bad_settings)))
  expect_false(grepl("failed:", joined(lab)))
  expect_true(all(c("Ba", "0", "TPM") %in% lab) && any(grepl("^Elf5-", lab)))  # Ba's tracks really drew
  # ...and that result isn't cached: once LP's file is fixed, LP's own transcripts count again.
  expect_length(.gene_display_transcripts_cache_load(file.path(bad_settings$cache_dir, "gene_display_transcripts_cache.rds")), 0)
})

test_that("tpm_usage_sample(): a non-number total is an error; named transcripts keep their order and the rest is 'Other'", {
  tx <- data.frame(transcript_name = c("G-201", "G-202", "G-203"), mean_tpm = c(5, 1, 4))
  quant <- function(total) structure(tx, total_tpm = total, n_quantified = 3L, tpm_not_numeric = character(0))
  local_mocked_bindings(filtered_transcripts = function(...) quant(10))
  expect_identical(tpm_usage_sample("G", "S", "quant.sf", c("G-203", "G-201"), "Other (n=1)", NULL),
                   c("G-203" = 4, "G-201" = 5, "Other (n=1)" = 1))
  expect_identical(tpm_usage_sample("G", "S", "quant.sf", c("G-201", "G-202", "G-203"), NULL, NULL),
                   c("G-201" = 5, "G-202" = 1, "G-203" = 4))
  # "inf" in a corrupt quant.sf reads as Inf and sums to Inf.
  for (bad in list(Inf, NA_real_, NULL)) {
    local_mocked_bindings(filtered_transcripts = function(...) quant(bad))
    expect_error(tpm_usage_sample("G", "S", "quant.sf", "G-201", "Other (n=2)", NULL),
                 "the TPM values for G in the quant file are not numbers (corrupt file?)", fixed = TRUE)
  }
})

test_that("a sample that can't be drawn doesn't change the shared TPM axis of the other tracks", {
  local_mocked_bindings(
    load_transcripts = function(entry, ...) { entry$data <- list(named = "G-201", n_annotated = 1L); entry },
    .pt_annotation = function(settings) NULL,
    filtered_transcripts = function(gene, samples, ...) {
      total <- if (samples[[1]]$sample_id == "bad") Inf else 40
      structure(data.frame(transcript_name = "G-201", mean_tpm = total), total_tpm = total,
                n_quantified = 1L, tpm_not_numeric = character(0))
    })
  settings <- .pt_fill_settings(list())  # the settings render_all() would pass on
  run <- function(samples, range) {
    rows <- lapply(samples, function(s) list(viz = "tpm_usage", sample = s))
    suppressMessages(settings_tpm_usage(list(tpm_range = range), list(viz = "tpm_usage", rows = rows),
                                        settings, 7.2, list(gene = "G")))
  }
  out <- run(c("ok", "bad"), 500)
  expect_identical(out$dict$tpm_range, 500)
  expect_null(out$entry$rows[[1]]$error)
  expect_match(out$entry$rows[[2]]$error, "^failed: the TPM values for G")  # the row is labelled "bad" already
  expect_identical(run(c("ok", "bad"), 10)$dict$tpm_range, 40)
  expect_message(settings_tpm_usage(list(), list(viz = "tpm_usage", rows = list(list(viz = "tpm_usage", sample = "bad"))),
                                    settings, 7.2, list(gene = "G")),
                 "PlotTrellis: tpm_usage (bad) track failed", fixed = TRUE)
})

test_that("tpm_usage total labels stay clear of the legend at large text and narrow pages (mouse dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  rows <- test_rows("mouse")
  specs <- list(single_track("tpm_usage", rows[["Ba"]]), single_track("tpm_usage", rows[["LP"]]))
  for (ov in list(list(fontsize = 16L), list(page_w = 6, fontsize = 12L), list(page_w = 5))) {
    svg <- render_svg(specs, test_region("mouse"), test_settings("mouse", overrides = ov))
    t <- svg_texts(svg)
    totals <- t[grepl("^[0-9.,]+ TPM$", text_label(t))]
    keys <- svg_keys(svg)
    expect_length(totals, 2)
    expect_gt(length(keys), 0)
    expect_lt(max(svg_attr(totals, "x") + svg_attr(totals, "textLength")), min(svg_attr(keys, "x")),
              label = paste(names(ov), unlist(ov), collapse = ", "))
  }
  # Too narrow for bars next to the legend: an error, not a backwards axis.
  lab <- text_label(svg_texts(render_svg(specs, test_region("mouse"), test_settings("mouse", overrides = list(page_w = 3)))))
  expect_true(any(grepl("tpm_usage (Ba, LP) track failed: not enough width", lab, fixed = TRUE)))
})

# Adjacent Transcript expression tracks are one panel.
test_that(".pt_group_row_panels() folds each run of adjacent tpm_usage tracks into one entry, in order", {
  tpm <- function(s) list(viz = "tpm_usage", sample = s)
  bw  <- list(viz = "bw", sample = "A")
  out <- .pt_group_row_panels(list(tpm("A"), tpm("B"), bw, bw, tpm("C")))
  expect_identical(vapply(out, `[[`, "", "viz"), c("tpm_usage", "bw", "bw", "tpm_usage"))
  expect_identical(out[[1]]$rows, list(tpm("A"), tpm("B")))
  expect_identical(out[[2]], bw)
  expect_identical(out[[4]]$rows, list(tpm("C")))
  # Every panel also knows all of the figure's rows of its type.
  expect_identical(out[[1]]$figure_rows, list(tpm("A"), tpm("B"), tpm("C")))
  expect_identical(out[[4]]$figure_rows, out[[1]]$figure_rows)
  expect_identical(.pt_group_row_panels(list()), list())
})

test_that("adjacent tpm_usage tracks draw as one panel: one heading, legend and axis, a row per sample in figure order (mouse dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  rows <- test_rows("mouse")
  settings <- test_settings("mouse")
  region <- test_region("mouse")
  ba <- single_track("tpm_usage", rows[["Ba"]]); lp <- single_track("tpm_usage", rows[["LP"]])
  count <- function(lab, what) sum(lab == what)
  label_y <- function(t, what) svg_attr(t[text_label(t) == what], "y")

  svg <- render_svg(list(ba, lp), region, settings)
  t <- svg_texts(svg); lab <- text_label(t)
  expect_identical(count(lab, "Transcript expression"), 1L)
  expect_identical(count(lab, "TPM"), 1L)
  expect_identical(count(lab, "Elf5-201"), 1L)
  expect_identical(count(lab, "0"), 1L)  # one axis
  expect_lt(label_y(t, "Ba"), label_y(t, "LP"))
  # Each sample's label, bar and total share a row.
  bars <- svg_bars(svg)
  bar_mid <- sort(unique(svg_attr(bars, "y"))) + unique(svg_attr(bars, "height")) / 2
  expect_length(bar_mid, 2)
  totals <- t[grepl("^[0-9.,]+ TPM$", lab)]
  expect_equal(svg_attr(totals, "y"), c(label_y(t, "Ba"), label_y(t, "LP")))
  expect_equal(bar_mid, c(label_y(t, "Ba"), label_y(t, "LP")) - 2.55, tolerance = 0.5 / 30)  # text baseline vs centre
  # Figure order, not sample order.
  t <- svg_texts(render_svg(list(lp, ba), region, settings))
  expect_gt(label_y(t, "Ba"), label_y(t, "LP"))

  # Two tracks merged are much shorter than two separate panels, which share one TPM range.
  page_h <- function(svg) as.numeric(sub(".*<svg[^>]* height='([0-9.]+)pt'.*", "\\1", svg))
  apart <- render_svg(list(ba, single_track("bw", rows[["Ba"]]), lp), region, settings)
  t <- svg_texts(apart); lab <- text_label(t)
  expect_identical(count(lab, "Transcript expression"), 2L)
  expect_identical(count(lab, "TPM"), 2L)
  zero <- t[lab == "0" & grepl("text-anchor='middle'", t)]
  five <- t[lab == "500" & grepl("text-anchor='middle'", t)]
  expect_length(zero, 2); expect_length(five, 2)
  expect_equal(diff(svg_attr(zero, "x")), 0)
  expect_equal(diff(svg_attr(five, "x")), 0)
  # Panels with different numbers of rows and a long sample name still share the scale
  # (human PTPRC: 9 legend entries, so the legend's columns matter).
  h <- test_rows("human")
  gm <- single_track("tpm_usage", h[["GM12878"]])
  k562 <- modifyList(single_track("tpm_usage", h[["K562"]]), list(sample = "K562 with a much longer sample name"))
  t <- svg_texts(render_svg(list(gm, single_track("bw", h[["GM12878"]]), gm, k562), test_region("human"), test_settings("human")))
  ticks <- t[grepl("^[0-9.,]+$", text_label(t)) & grepl("text-anchor='middle'", t)]
  by_panel <- split(ticks, svg_attr(ticks, "y"))
  expect_length(by_panel, 2)
  expect_identical(text_label(by_panel[[1]]), text_label(by_panel[[2]]))
  expect_equal(svg_attr(by_panel[[1]], "x"), svg_attr(by_panel[[2]], "x"))
  expect_gt(length(by_panel[[1]]), 2)
  # A long name on a narrow page is cut short; it doesn't fail the figure's panels.
  long <- modifyList(k562, list(sample = "K562_RNAseq_rep1_batch2_x"))
  lab <- text_label(svg_texts(render_svg(list(gm, single_track("bw", h[["GM12878"]]), long), test_region("human"),
                                         test_settings("human", overrides = list(page_w = 5)))))
  expect_false(any(grepl("track failed", lab)))
  expect_identical(sum(grepl("^[0-9.,]+ TPM$", lab)), 2L)
  expect_identical(sum(grepl("^K562_RNA.*\\.\\.\\.$", lab)), 1L)
  # Each panel's slot is as tall as what is drawn in it, also where another panel's long name
  # changes the legend's columns (fontsize 12 on an 8 in page: two columns, one with the name).
  st <- test_settings("human", overrides = list(fontsize = 12L, page_w = 8))
  heading_gap <- function(last) {
    t <- svg_texts(render_svg(list(gm, single_track("bw", h[["GM12878"]]), last), test_region("human"), st))
    list(gap = diff(svg_attr(t[text_label(t) == "Transcript expression"], "y")),
         keys = unique(grep("^(PTPRC-|Other )", text_label(t), value = TRUE)))
  }
  short <- heading_gap(modifyList(k562, list(sample = "K562")))
  long <- heading_gap(k562)
  height <- function(samples) transcript_counts_layout(1, short$keys, 8 - 1.8, 12, samples)$height
  expect_gt(height(c("GM12878", k562$sample)), height(c("GM12878", "K562")) + 0.3)
  expect_equal(long$gap - short$gap, 72 * (height(c("GM12878", k562$sample)) - height(c("GM12878", "K562"))),
               tolerance = 1e-6)

  # The second row adds one bar and one gap to the page: 0.26 in.
  expect_equal(page_h(svg) - page_h(render_svg(list(ba), region, settings)), 72 * 0.26, tolerance = 1e-6)

  # A sample that can't be drawn says so in its own row; the other row and the rest still draw.
  svg <- render_svg(list(ba, modifyList(lp, list(quant_path = single_track("tpm_usage", test_rows("human")[["K562"]])$quant_path))),
                    region, settings)
  t <- svg_texts(svg); lab <- text_label(t)
  failed <- t[grepl("fill: #CC0000", t, fixed = TRUE)]
  expect_identical(paste(text_label(failed), collapse = " "),
                   "failed: none of Elf5's transcripts are in this sample's quant file (is the GTF the one the quant files were made with?)")
  expect_identical(count(lab, "Transcript expression"), 1L)
  expect_true(all(c("Ba", "LP", "TPM", "Elf5-201") %in% lab))
  expect_length(unique(svg_attr(svg_bars(svg), "y")), 1)
  expect_equal(svg_attr(failed, "y")[1] < label_y(t, "LP") + 6 && svg_attr(failed, "y")[1] > label_y(t, "Ba"), TRUE)
  # On a narrow page the two lines still start with the reason, not with the track's name.
  t <- svg_texts(render_svg(list(ba, modifyList(lp, list(quant_path = single_track("tpm_usage", test_rows("human")[["K562"]])$quant_path))),
                            region, test_settings("mouse", overrides = list(page_w = 5))))
  failed <- text_label(t[grepl("fill: #CC0000", t, fixed = TRUE)])
  expect_length(failed, 2)
  expect_match(paste(failed, collapse = " "), "^failed: none of Elf5's transcripts.*\\.\\.\\.$")
})

test_that(".pt_track_error_label() names a row panel by its samples, the first three of many", {
  panel <- function(samples) list(viz = "tpm_usage", rows = lapply(samples, function(s) list(viz = "tpm_usage", sample = s)))
  label <- function(entry) suppressMessages(.pt_track_error_label(entry, simpleError("boom")))
  expect_identical(label(panel("A")), "tpm_usage (A) track failed: boom")
  expect_identical(label(panel(c("A", "B", "C"))), "tpm_usage (A, B, C) track failed: boom")
  expect_identical(label(panel(LETTERS[1:12])), "tpm_usage (A, B, C, +9 more) track failed: boom")
  expect_identical(label(list(viz = "bw", sample = "A")), "bw (A) track failed: boom")
  expect_identical(label(list(viz = "genes")), "genes track failed: boom")
})

test_that(".pt_track_error_label() cuts a long reason in the middle, keeping its start and the path's end", {
  deep <- file.path(paste(rep("a-deep-project-folder", 12), collapse = "/"), "quant", "LP_rep1")
  e <- simpleError(paste0("file not found: quant.sf (in ", deep, ")"))
  lab <- suppressMessages(.pt_track_error_label(list(viz = "tpm_usage", sample = "Basal"), e))
  reason <- sub("^tpm_usage \\(Basal\\) track failed: ", "", lab)
  expect_identical(nchar(reason), 140L)
  expect_match(reason, "^file not found: quant\\.sf \\(in a-deep-project-folder/")
  expect_match(reason, "/a-deep-project-folder/quant/LP_rep1\\)$")
  expect_match(reason, "...", fixed = TRUE)
  # A reason of up to 140 characters is kept whole.
  short <- strrep("x", 140)
  expect_identical(suppressMessages(.pt_track_error_label(list(viz = "genes"), simpleError(short))),
                   paste("genes track failed:", short))
})

# draw_transcript_counts() on its own page, no fixture data: `tpm` rows are samples.
counts_svg <- function(tpm, ..., width = 7.2, fontsize = 8, label_samples = rownames(tpm)) {
  f <- tempfile(fileext = ".svg")
  on.exit(unlink(f), add = TRUE)
  h <- transcript_counts_layout(nrow(tpm), colnames(tpm), width, fontsize, label_samples)$height
  svglite::svglite(f, width = width + 1, height = h + 1)
  dev <- grDevices::dev.cur()
  on.exit(if (dev %in% grDevices::dev.list()) grDevices::dev.off(dev), add = TRUE, after = FALSE)
  utils::capture.output({
    plotgardener::pageCreate(width = width + 1, height = h + 1, default.units = "inches", showGuides = FALSE)
    draw_transcript_counts(tpm, x = 0.5, y = 0.5, width = width, fontsize = fontsize, fontfamily = "DejaVu Sans",
                           label_samples = label_samples, ...)
  })
  grDevices::dev.off(dev)
  paste(readLines(f), collapse = "\n")
}
counts_matrix <- function(n_samples, n_tx) {
  matrix(seq_len(n_samples * n_tx) * 10, n_samples, n_tx,
         dimnames = list(paste0("S", seq_len(n_samples)), sprintf("GENE-2%02d", seq_len(n_tx))))
}

test_that("draw_transcript_counts(): bars are 0.2 in thick whatever the number of rows, in a panel as tall as its rows need", {
  skip_on_os("windows")
  for (n in 1:4) {
    tpm <- counts_matrix(n, 3)
    svg <- counts_svg(tpm, transcript_colours = TX_PALETTE[1:3])
    bars <- svg_bars(svg)
    expect_length(bars, n * 3)
    expect_equal(unique(svg_attr(bars, "height")), 0.2 * 72)
    # Rows are 0.06 in apart and nothing of the panel is drawn outside its layout height.
    expect_equal(diff(sort(unique(svg_attr(bars, "y")))), rep(0.26 * 72, n - 1))
    lay <- transcript_counts_layout(n, colnames(tpm), 7.2)
    expect_equal(lay$height, lay$heading_h + n * 0.2 + (n - 1) * 0.06 + lay$axis_h)
    top <- 0.5 * 72; bottom <- (0.5 + lay$height) * 72
    t <- svg_texts(svg)
    expect_true(all(svg_attr(t, "y") > top & svg_attr(t, "y") <= bottom + 0.01))
    expect_true(all(svg_attr(bars, "y") >= top & svg_attr(bars, "y") + svg_attr(bars, "height") <= bottom))
    # Bars end at their totals on the axis.
    axis <- svg_tpm_axis(t)
    for (i in seq_len(n)) {
      seg <- bars[svg_attr(bars, "y") == sort(unique(svg_attr(bars, "y")))[i]]
      expect_equal(max(svg_attr(seg, "x") + svg_attr(seg, "width")), tpm_px(axis, sum(tpm[i, ])), tolerance = 0.5 / 300)
    }
  }
  expect_lt(transcript_counts_layout(1, colnames(counts_matrix(1, 3)), 7.2)$height, 0.7)
})

test_that("draw_transcript_counts(): no horizontal gridline through the rows", {
  skip_on_os("windows")
  svg <- counts_svg(rbind(counts_matrix(1, 3), Low = c(0.5, 0.1, 0.2)), transcript_colours = TX_PALETTE[1:3])
  lines <- svg_lines(svg)
  expect_gt(length(lines), 0)
  expect_true(all(svg_attr(lines, "x1") == svg_attr(lines, "x2")))
  expect_true("0.8 TPM" %in% text_label(svg_texts(svg)))
})

test_that("draw_transcript_counts(): the legend wraps into columns inside the panel and the bars stay clear of it", {
  skip_on_os("windows")
  for (n in 1:3) {
    tpm <- counts_matrix(n, 9)
    svg <- counts_svg(tpm, transcript_colours = c(TX_PALETTE, TX_OTHER_COLOUR))
    keys <- svg_keys(svg)
    t <- svg_texts(svg); lab <- text_label(t)
    expect_length(keys, 9)
    lay <- transcript_counts_layout(n, colnames(tpm), 7.2)
    expect_true(all(svg_attr(keys, "y") >= 0.5 * 72 & svg_attr(keys, "y") + svg_attr(keys, "height") <= (0.5 + lay$height) * 72))
    expect_identical(length(unique(svg_attr(keys, "x"))), 2L)  # whatever the number of rows
    # Legend labels end at the panel's right edge; total labels end left of the legend.
    legend <- t[grepl("^GENE-", lab)]
    expect_equal(max(svg_attr(legend, "x") + svg_attr(legend, "textLength")), (0.5 + 7.2) * 72, tolerance = 0.5 / 500)
    totals <- t[grepl(" TPM$", lab)]
    expect_length(totals, n)
    expect_lt(max(svg_attr(totals, "x") + svg_attr(totals, "textLength")), min(svg_attr(keys, "x")))
  }
})

test_that("transcript_counts_layout(): the legend takes fewer columns and the panel grows when the bars would be squeezed", {
  keys <- c(sprintf("PTPRC-2%02d", 1:8), "Other (n=15)")
  rows_h <- function(lay, n) lay$heading_h + n * lay$bar_h + (n - 1) * lay$gap + lay$axis_h
  # Room to spare: two columns beside one row, and the panel is only as tall as its row.
  wide <- transcript_counts_layout(1, keys, 7.2)
  expect_identical(wide$legend_cols, 2)
  expect_equal(wide$height, rows_h(wide, 1))
  # A 5 in page: one column, and the panel is as tall as its 9 legend entries.
  narrow <- transcript_counts_layout(1, keys, 3.2)
  expect_identical(narrow$legend_cols, 1)
  expect_equal(narrow$height, 9 * narrow$key_pitch)
  expect_gt(narrow$height, rows_h(narrow, 1))
  # Fontsize 16 on the default page: three columns would fit the height but squeeze the bars.
  big <- transcript_counts_layout(1, keys, 7.2, fontsize = 16)
  expect_identical(big$legend_cols, 2)
  expect_identical(transcript_counts_layout(2, character(0), 7.2)$legend_cols, 1)  # no legend (no gene)
  # The column count doesn't depend on the number of rows (panels in one figure line up)...
  expect_identical(transcript_counts_layout(4, keys, 7.2)$legend_cols, 2)
  # ...and a long sample name counts as part of the squeeze: one column where two fit before.
  expect_identical(transcript_counts_layout(1, keys, 6.2, fontsize = 12)$legend_cols, 2)
  expect_identical(transcript_counts_layout(1, keys, 6.2, fontsize = 12, label_samples = strrep("x", 34))$legend_cols, 1)
  # The label column stops growing at 40% of the width, so a longer name changes nothing more.
  expect_identical(transcript_counts_layout(1, keys, 7.2, label_samples = strrep("x", 60)),
                   transcript_counts_layout(1, keys, 7.2, label_samples = strrep("x", 200)))
})

test_that("draw_transcript_counts(): 9 transcripts draw at large text and on narrow pages, everything clear of the legend", {
  skip_on_os("windows")
  # (fontsize, page width): the track is 1.8 in narrower than the page.
  # Also with a long sample name, which narrows the bars further.
  # Every (fontsize, page width) here draws with a short name, so it has to with a long one.
  for (case in list(c(8, 5), c(12, 6), c(16, 7), c(16, 8), c(16, 9), c(24, 9), c(8, 9),
                    c(12, 6, 20), c(12, 8, 34), c(12, 9, 34), c(12, 7, 34), c(8, 5, 25), c(10, 6, 30),
                    c(16, 7, 16), c(16, 8, 30), c(24, 9, 34))) for (n in 1:2) {
    what <- sprintf("fontsize %d, page_w %d, %d row(s), sample name %s", case[1], case[2], n,
                    if (length(case) > 2) paste(case[3], "characters") else "short")
    tpm <- counts_matrix(n, 9)
    colnames(tpm) <- c(sprintf("PTPRC-2%02d", 1:8), "Other (n=15)")
    if (length(case) > 2) {
      rownames(tpm)[1] <- substr("A sample name far longer than usual", 1, case[3])
      tpm <- tpm / 10  # two-digit totals, as in the real case (PTPRC) these settings failed on
    }
    svg <- counts_svg(tpm, transcript_colours = c(TX_PALETTE, TX_OTHER_COLOUR), width = case[2] - 1.8, fontsize = case[1])
    keys <- svg_keys(svg); bars <- svg_bars(svg)
    t <- svg_texts(svg); lab <- text_label(t)
    expect_length(keys, 9)
    expect_length(bars, n * 9)
    lay <- transcript_counts_layout(n, colnames(tpm), case[2] - 1.8, case[1], rownames(tpm))
    expect_true(all(svg_attr(keys, "y") >= 0.5 * 72 - 0.01 &
                    svg_attr(keys, "y") + svg_attr(keys, "height") <= (0.5 + lay$height) * 72 + 0.01), label = what)
    # Total labels and tick labels end left of the legend.
    left_of_legend <- t[grepl(" TPM$", lab) | (grepl("^[0-9.,]+$", lab) & grepl("text-anchor='middle'", t))]
    right <- svg_attr(left_of_legend, "x") + svg_attr(left_of_legend, "textLength") /
      ifelse(grepl("text-anchor='middle'", left_of_legend), 2, 1)
    expect_lt(max(right), min(svg_attr(keys, "x")), label = what)
    # Row labels (cut short where the name is too long for its column) end left of the bars.
    names_el <- t[abs(svg_attr(t, "x") - 0.55 * 72) < 0.01 & lab != "Transcript expression"]
    expect_length(names_el, n)
    expect_true(all(startsWith(rownames(tpm), sub("\\.\\.\\.$", "", text_label(names_el)))), label = what)
    expect_lt(max(svg_attr(names_el, "x") + svg_attr(names_el, "textLength")), min(svg_attr(bars, "x")), label = what)
  }
})

test_that("draw_transcript_counts(): panels given the same label_samples and tpm_max share one TPM scale", {
  skip_on_os("windows")
  keys <- c(sprintf("PTPRC-2%02d", 1:8), "Other (n=15)")
  one <- counts_matrix(1, 9); three <- counts_matrix(3, 9)
  rownames(three)[2] <- "A much longer sample name than usual"
  colnames(one) <- colnames(three) <- keys
  samples <- c(rownames(one), rownames(three))
  scale_of <- function(tpm) {
    svg <- counts_svg(tpm, transcript_colours = c(TX_PALETTE, TX_OTHER_COLOUR), tpm_max = 1000, label_samples = samples)
    ticks <- svg_tpm_axis(svg_texts(svg))$ticks
    list(labels = text_label(ticks), x = svg_attr(ticks, "x"), bars_left = min(svg_attr(svg_bars(svg), "x")))
  }
  expect_identical(scale_of(one), scale_of(three))
  expect_gt(length(scale_of(one)$x), 2)
})

test_that("draw_transcript_counts(): the sample label column widens for a long sample name", {
  skip_on_os("windows")
  bars_left <- function(svg) min(svg_attr(svg_bars(svg), "x"))
  tpm <- counts_matrix(2, 3)
  expect_equal(bars_left(counts_svg(tpm, transcript_colours = TX_PALETTE[1:3])), (0.5 + 1.4) * 72)
  rownames(tpm)[2] <- "A much longer sample name than usual"
  svg <- counts_svg(tpm, transcript_colours = TX_PALETTE[1:3])
  t <- svg_texts(svg)
  name <- t[text_label(t) == rownames(tpm)[2]]
  expect_gt(bars_left(svg), svg_attr(name, "x") + svg_attr(name, "textLength") + 0.1 * 72)

  # ...up to 40% of the panel's width: a name longer than that is cut short.
  long <- paste(rep("A much longer sample name than usual", 3), collapse = " ")
  rownames(tpm)[2] <- long
  svg <- counts_svg(tpm, transcript_colours = TX_PALETTE[1:3])
  expect_equal(bars_left(svg), (0.5 + 0.4 * 7.2) * 72)
  t <- svg_texts(svg); lab <- text_label(t)
  name <- t[grepl("^A much longer", lab)]
  expect_length(name, 1)
  expect_match(text_label(name), "\\.\\.\\.$")
  expect_true(startsWith(long, sub("\\.\\.\\.$", "", text_label(name))))
  expect_gt(nchar(text_label(name)), 30)
  expect_lt(svg_attr(name, "x") + svg_attr(name, "textLength"), bars_left(svg) - 0.1 * 72)
  expect_true("S1" %in% lab)  # a name that fits is left alone
  # On a narrow panel the column gives way to the bars, down to its usual 1.4 in.
  keys9 <- c(sprintf("PTPRC-2%02d", 1:8), "Other (n=15)")
  tpm9 <- counts_matrix(1, 9) / 10
  dimnames(tpm9) <- list("K562_RNAseq_rep1_batch2_x", keys9)
  svg <- counts_svg(tpm9, transcript_colours = c(TX_PALETTE, TX_OTHER_COLOUR), width = 3.2)
  expect_equal(bars_left(svg), (0.5 + 1.4) * 72)
  expect_length(svg_bars(svg), 9)
})

test_that("draw_transcript_counts(): a row's error replaces its bar; with no row to draw there is no axis or legend", {
  skip_on_os("windows")
  tpm <- counts_matrix(2, 3)
  msg <- paste("tpm_usage (S2) track failed:", paste(rep("a long reason", 8), collapse = " "))
  svg <- counts_svg(tpm, transcript_colours = TX_PALETTE[1:3], errors = c(NA, msg))
  t <- svg_texts(svg); lab <- text_label(t)
  bars <- svg_bars(svg)
  expect_length(unique(svg_attr(bars, "y")), 1)
  expect_identical(sum(grepl(" TPM$", lab)), 1L)
  failed <- t[grepl("a long reason", lab)]
  expect_length(failed, 2)  # wrapped onto two lines
  expect_identical(paste(text_label(failed), collapse = " "), msg)
  expect_true(all(grepl("fill: #CC0000", failed)))
  # Both lines stay left of the legend and inside the failed row's band.
  expect_lt(max(svg_attr(failed, "x") + svg_attr(failed, "textLength")), min(svg_attr(svg_keys(svg), "x")))
  row2_top <- min(svg_attr(bars, "y")) + 0.26 * 72
  expect_true(all(svg_attr(failed, "y") > row2_top - 2 & svg_attr(failed, "y") < row2_top + 0.2 * 72 + 4))

  # Too long for two lines on a narrow panel: cut, still left of the legend.
  svg <- counts_svg(tpm, transcript_colours = TX_PALETTE[1:3], errors = c(NA, msg), width = 4.2, fontsize = 12)
  t <- svg_texts(svg)
  failed <- t[grepl("fill: #CC0000", t, fixed = TRUE)]
  expect_length(failed, 2)
  expect_true(startsWith(msg, sub("\\.\\.\\.$", "", paste(text_label(failed), collapse = " "))))
  # The cut is at the end of a word, unless over 12 characters of the next word fit.
  svg <- counts_svg(tpm, transcript_colours = TX_PALETTE[1:3], width = 4.2, fontsize = 12,
                    errors = c(NA, paste(strrep("y", 300), "ab", strrep("x", 300))))
  t <- svg_texts(svg)
  expect_match(text_label(t[grepl("fill: #CC0000", t, fixed = TRUE)])[2], "^ab x{13,}\\.\\.\\.$")
  for (width in c(4.2, 4.4, 4.6)) {
    svg <- counts_svg(tpm, transcript_colours = TX_PALETTE[1:3], width = width, fontsize = 12,
                      errors = c(NA, paste("failed:", paste(rep("transcripts", 12), collapse = " "))))
    t <- svg_texts(svg)
    expect_match(text_label(t[grepl("fill: #CC0000", t, fixed = TRUE)])[2], "(^| )transcripts\\.\\.\\.$")
  }
  expect_lt(max(svg_attr(failed, "x") + svg_attr(failed, "textLength")), min(svg_attr(svg_keys(svg), "x")))
  # A single word too long for a line (a path) is cut too.
  svg <- counts_svg(tpm, transcript_colours = TX_PALETTE[1:3], errors = c(NA, paste0("failed: /", strrep("folder/", 40))))
  t <- svg_texts(svg)
  failed <- t[grepl("fill: #CC0000", t, fixed = TRUE)]
  expect_length(failed, 2)
  expect_lt(max(svg_attr(failed, "x") + svg_attr(failed, "textLength")), min(svg_attr(svg_keys(svg), "x")))

  svg <- counts_svg(tpm, transcript_colours = TX_PALETTE[1:3], errors = c("first failed", "second failed"))
  lab <- text_label(svg_texts(svg))
  expect_setequal(lab, c("Transcript expression", "S1", "S2", "first failed", "second failed"))
  expect_error(counts_svg(unname(tpm), transcript_colours = TX_PALETTE[1:3]), "needs row names")
  expect_length(svg_bars(svg), 0)
  expect_length(svg_keys(svg), 0)
})
