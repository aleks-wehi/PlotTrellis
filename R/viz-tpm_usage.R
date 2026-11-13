info_tpm_usage <- function() list(
  kind = "sample", data_type = "quantdir", order = 50, label = "Transcript expression", api = "transcript_expression",
  filter = c(tpm_usage = "TPMs"), row_panel = TRUE, needs_gene = TRUE
)

# One Transcript expression track per sample with a quantification, whatever its name.
detect_tpm_usage <- function(sample_row, info) {
  quant <- sample_data_path(sample_row, "quantdir")
  if (is.na(quant)) return(list())
  safe_s <- gsub("[^A-Za-z0-9._-]", "_", sample_row$display)
  list(list(
    key = paste0("tpm_usage_", safe_s),
    label = paste0(info$label, " \u00B7 ", sample_row$display),
    type = "tpm_usage", sample = sample_row$display, file = NA_character_, colour_key = sample_row$display,
    quant_path = quant,
    broken = !file.exists(quant)
  ))
}

# Transcript expression: one horizontal stacked bar per sample, each segment one transcript's TPM.
# Adjacent Transcript expression tracks reach this module as ONE entry with a row per track
# (row_panel in info_tpm_usage() above) and draw as one panel: the heading once, a bar per sample,
# one legend, one TPM axis. Separate panels in one figure draw on the same TPM scale:
# dict[["tpm_range"]] is the largest total among all of the figure's Transcript expression rows,
# and every panel sizes its label column from all of their sample names (entry$figure_rows).
# A sample name too long for that column is cut short ("...").
# The panel's height follows from its rows and legend; a height preset on a track is not used.
settings_tpm_usage <- function(dict, entry, settings, track_w, region_params) {
  gene <- region_params$gene
  if (!.pt_no_gene(gene)) {
    entry <- load_transcripts(entry, region_params, settings)
    named <- entry$data$named
    n_other <- entry$data$n_annotated - length(named)
    # The gene's transcripts not shown by name are one "Other" segment.
    other <- if (n_other > 0) sprintf("Other (n=%d)", n_other)
    entry$data$keys <- c(named, other)
  }
  entry$height <- transcript_counts_layout(length(entry$rows), entry$data$keys, track_w, settings$fontsize,
                                           label_samples = vapply(entry$figure_rows, `[[`, "", "sample"))$height
  if (.pt_no_gene(gene)) return(list(dict = dict, entry = entry))

  anno <- .pt_annotation(settings)
  entry$rows <- lapply(entry$rows, function(row) {
    # A sample that can't be drawn keeps its row and says why there; the others still draw.
    # The console gets the usual "tpm_usage (<sample>) track failed: <why>"; the row, which is
    # labelled with the sample already, only "failed: <why>".
    tryCatch({
      row$tpm <- tpm_usage_sample(gene, row$sample, row$quant_path, named, other, anno)
      row
    }, error = function(e) {
      row$error <- sub("^.*? track failed: ", "failed: ", .pt_track_error_label(row, e), perl = TRUE)
      row
    })
  })
  totals <- vapply(entry$rows, function(row) sum(row$tpm), numeric(1))
  dict[["tpm_range"]] <- max(dict[["tpm_range"]], totals)
  list(dict = dict, entry = entry)
}

# One sample's TPM for each named transcript, in `named` order, then (when `other` is its
# label, not NULL) the summed TPM of the gene's transcripts not shown by name. Stops when there
# is nothing to draw a bar from; a gene found in the quant file at 0 TPM is a real result and
# gives all zeros.
tpm_usage_sample <- function(gene, sample, quant_path, named, other, anno) {
  tx <- filtered_transcripts(gene = gene, samples = list(list(sample_id = sample, quant_path = quant_path)),
                             anno = anno)
  if (length(attr(tx, "tpm_not_numeric")))
    stop("the TPM column of this sample's quant file is not all numbers (corrupt file?): ",
         attr(tx, "tpm_not_numeric")[1], call. = FALSE)
  if (identical(attr(tx, "n_quantified"), 0L))
    stop("none of ", gene, "'s transcripts are in this sample's quant file (is the GTF the one the quant files were made with?)", call. = FALSE)
  total <- attr(tx, "total_tpm")
  if (is.null(total) || !is.finite(total))
    stop("the TPM values for ", gene, " in the quant file are not numbers (corrupt file?)", call. = FALSE)

  tpm <- vapply(named, function(name) sum(tx$mean_tpm[tx$transcript_name == name]), numeric(1))
  if (!is.null(other)) tpm[[other]] <- max(0, total - sum(tpm))
  tpm
}

render_tpm_usage <- function(entry, dict, ctx) {
  if (.pt_no_gene(ctx$region_params$gene)) {
    .pt_status_text(ctx, .pd_status_text("no_gene"))
    return(invisible(NULL))
  }
  rows <- entry$rows
  named <- entry$data$named
  keys <- entry$data$keys
  colours <- unname(transcript_colours(named)[named])
  if (length(keys) > length(named)) colours <- c(colours, TX_OTHER_COLOUR)

  # A row with an error has no TPMs: zeros keep its place, its message is drawn instead.
  tpm <- do.call(rbind, lapply(rows, function(row) row$tpm %||% rep(0, length(keys))))
  dimnames(tpm) <- list(vapply(rows, `[[`, "", "sample"), keys)
  draw_transcript_counts(
    tpm, x = ctx$track_x, y = ctx$y, width = ctx$track_w,
    transcript_colours = colours, sample_colours = vapply(rows, `[[`, "", "colour"),
    tpm_max = dict[["tpm_range"]], errors = vapply(rows, function(row) row$error %||% NA_character_, ""),
    label_samples = vapply(entry$figure_rows, `[[`, "", "sample"),
    fontsize = ctx$fontsize, fontfamily = ctx$font, fontcolor = ctx$fontcol
  )
  invisible(NULL)
}

# Sizes (inches) of a Transcript expression panel `width` wide, with `n` sample rows and a legend
# entry per name in `transcripts`: a heading line, `n` bars `bar_h` thick and `gap` apart, the
# TPM axis below them, and the legend in `legend_cols` columns. Text-sized parts grow with
# `fontsize`.
#
# The legend takes as few columns as fit beside a single row, so a short panel doesn't grow to
# hold it. Where that would leave the bars under 40% of the width right of the sample labels
# (`label_samples`), it takes fewer columns and the panel grows instead. Neither choice depends
# on `n`, so panels with the same transcripts, width and `label_samples` line their bars up.
# Text widths are estimated from the label lengths here: the height is needed before there is a
# page to measure text on. Pass the `label_samples` that draw_transcript_counts() will get.
transcript_counts_layout <- function(n, transcripts, width, fontsize = 8, label_samples = character(0)) {
  line <- (fontsize - 1) / 72  # height of the panel's label text
  est_w <- function(labels) 0.55 * line * max(nchar(labels), 0)
  heading_h <- line + 0.08
  bar_h <- max(0.2, 2 * line)
  gap <- 0.06
  axis_h <- 0.07 + (fontsize - 2) / 72 + (fontsize - 1.5) / 72  # tick marks, tick labels, "TPM"
  rows_h <- function(n) heading_h + n * bar_h + (n - 1) * gap + axis_h

  # Sample label column, left of the bars: at most 40% of the width (longer labels are cut).
  label_w <- max(1.4, min(est_w(label_samples) + 0.25, 0.4 * width))
  key_pitch <- 1.15 * line
  legend_col_w <- 0.75 * line + 0.05 + est_w(transcripts) + 0.12
  legend_cols <- max(1, ceiling(length(transcripts) / floor(rows_h(1) / key_pitch)))
  while (legend_cols > 1 && width - label_w - legend_cols * legend_col_w < 0.4 * (width - label_w))
    legend_cols <- legend_cols - 1
  legend_h <- ceiling(length(transcripts) / legend_cols) * key_pitch
  list(heading_h = heading_h, bar_h = bar_h, gap = gap, axis_h = axis_h,
       key_pitch = key_pitch, legend_cols = legend_cols, height = max(rows_h(n), legend_h))
}

# Draws a Transcript expression panel on the current plotgardener page, with its top left corner
# at (x, y) inches, `width` wide and
# transcript_counts_layout(nrow(tpm), colnames(tpm), width, fontsize, label_samples)$height tall.
#
#   tpm                 numeric matrix of TPMs: one row per sample (rownames are the row labels),
#                       one column per transcript (colnames are the legend labels). Bars are
#                       drawn top to bottom in row order, segments left to right in column order.
#   transcript_colours  one fill colour per column of `tpm`
#   sample_colours      one label colour per row of `tpm`
#   tpm_max             the largest total the axis has to fit (pass the same value to several
#                       panels to give them the same axis)
#   errors              optional, one per row: NA, or a message to draw in place of that bar
#   label_samples       the sample labels the label column is sized for. Give several panels
#                       the same set (all of their row labels) and the same `tpm_max`, and
#                       their bars line up on one TPM scale.
#
# The label column is as wide as its longest label needs, up to 40% of `width` and to what
# leaves the bars room; a label too long for it is cut short ("...").
#
# Needs only grid and plotgardener. Each bar ends at its sample's total on the TPM axis; a
# segment is labelled with its share of the bar when that is at least 5% and the label fits.
draw_transcript_counts <- function(tpm, x, y, width, transcript_colours, sample_colours = "black",
                                   tpm_max = max(rowSums(tpm)), errors = NULL, label_samples = rownames(tpm),
                                   heading = "Transcript expression",
                                   fontsize = 8, fontfamily = "", fontcolor = "black") {
  stopifnot("`tpm` needs row names (samples) and column names (transcripts)" =
              !is.null(rownames(tpm)) && !is.null(colnames(tpm)))
  n <- nrow(tpm)
  lay <- transcript_counts_layout(n, colnames(tpm), width, fontsize, label_samples)
  if (is.null(errors)) errors <- rep(NA_character_, n)
  drawn <- is.na(errors)
  totals <- rowSums(tpm)

  text_w <- function(labels, size, face = "plain") vapply(labels, function(label)
    grid::convertWidth(grid::grobWidth(grid::textGrob(label, gp = grid::gpar(
      fontsize = size, fontface = face, fontfamily = fontfamily))), "inches", valueOnly = TRUE), numeric(1))
  text <- function(labels, tx, ty, just, size = fontsize - 1, colour = fontcolor, face = "plain")
    plotText(label = labels, x = tx, y = ty, just = just, default.units = "inches",
             fontsize = size, fontface = face, fontfamily = fontfamily, fontcolor = colour)
  # No scientific notation; whole numbers from 10 TPM up, 2 significant digits below.
  tpm_label <- function(v) vapply(v, function(t) format(if (t >= 10) round(t) else signif(t, 2),
                                                        big.mark = ",", scientific = FALSE), "")
  # `label` cut short ("...") where it is wider than `max_w`; with `words`, back at the end of
  # the last whole word, unless over 12 characters of the next word fit.
  cut_text <- function(label, max_w, words = FALSE) {
    if (text_w(label, fontsize - 1) <= max_w + 1e-6) return(label)
    full <- label
    while (nchar(label) > 1 && text_w(paste0(label, "..."), fontsize - 1) > max_w) label <- substr(label, 1, nchar(label) - 1)
    if (words && substr(full, nchar(label) + 1, nchar(label) + 1) != " ") label <- sub(" [^ ]{1,12}$", "", label)
    paste0(trimws(label), "...")
  }

  # Legend sizes: the layout's columns, filled top to bottom.
  keys <- colnames(tpm)
  key <- 0.75 * (fontsize - 1) / 72
  pitch <- lay$key_pitch
  per_col <- ceiling(length(keys) / lay$legend_cols)
  key_col <- (seq_along(keys) - 1) %/% per_col + 1
  key_row <- (seq_along(keys) - 1) %% per_col
  col_w <- key + 0.05 + tapply(text_w(keys, fontsize - 1), key_col, max) + 0.12
  legend_w <- sum(col_w) - 0.12
  legend_h <- per_col * pitch

  # Rows: the heading line on top, then one bar per sample. Left of the bars, the sample label
  # column: as wide as its longest label needs, up to 40% of the width and to what leaves the
  # bars their least width (room for 1.5 total labels).
  bar_y <- y + lay$heading_h + (seq_len(n) - 1) * (lay$bar_h + lay$gap)
  mid_y <- bar_y + lay$bar_h / 2
  bars_bottom <- bar_y[n] + lay$bar_h
  label_room <- text_w(paste(tpm_label(tpm_max), "TPM"), fontsize - 1) + 0.05
  min_bars_w <- max(0.5, 1.5 * label_room)
  label_w <- max(1.4, min(max(text_w(label_samples, fontsize - 1), 0) + 0.25, 0.4 * width,
                          width - legend_w - 0.15 - min_bars_w))
  bars_x <- x + label_w
  text(heading, x + 0.05, y + 0.02, c("left", "top"), face = "bold")
  text(vapply(rownames(tpm), cut_text, "", max_w = label_w - 0.25), x + 0.05, mid_y, c("left", "center"),
       colour = sample_colours)

  # An error message takes its bar's place: on two lines if one is too long, each cut ("...")
  # where it would run past `max_w`.
  draw_errors <- function(max_w) for (i in which(!drawn)) {
    words <- strsplit(errors[i], " ", fixed = TRUE)[[1]]
    first <- max(1, sum(cumsum(text_w(paste0(words, " "), fontsize - 1)) <= max_w))
    lines <- c(paste(words[seq_len(first)], collapse = " "), paste(words[-seq_len(first)], collapse = " "))
    lines <- vapply(lines[nzchar(lines)], cut_text, "", max_w = max_w, words = TRUE)
    text(lines, bars_x, mid_y[i] + (seq_along(lines) - (length(lines) + 1) / 2) * fontsize / 72,
         c("left", "center"), colour = "#CC0000")
  }
  if (!any(drawn)) {
    draw_errors(width - label_w)
    return(invisible(NULL))
  }

  # Legend, against the panel's right edge.
  key_x <- (x + width - legend_w + cumsum(c(0, col_w)))[key_col]
  # Centred on the bars; kept inside the panel when taller than them.
  key_y <- min(max((bar_y[1] + bars_bottom - legend_h) / 2, y), y + lay$height - legend_h) + key_row * pitch
  plotRect(x = key_x, y = key_y + (pitch - key) / 2, width = key, height = key, just = c("left", "top"),
           default.units = "inches", fill = transcript_colours, linecolor = NA)
  text(keys, key_x + key + 0.05, key_y + pitch / 2, c("left", "center"))

  # TPM axis: runs 1/9 past the largest total, further if that bar's total label (drawn
  # 0.05 in after the bar) needs the room to stay clear of the legend.
  bars_w <- width - label_w - legend_w - 0.15
  if (bars_w < min_bars_w - 1e-9)
    stop("not enough width for the bars between the sample names and the legend; make the page wider or the font smaller", call. = FALSE)
  tpm_range <- max(1, tpm_max) * max(10 / 9, 1 / (1 - label_room / bars_w))
  tpm_x <- function(v) bars_x + v / tpm_range * bars_w
  # Fewer ticks on a narrow panel, so their labels never overlap; none whose label would
  # reach the legend.
  for (n_ticks in 5:1) {
    ticks <- pretty(c(0, tpm_range), n = n_ticks)
    ticks <- ticks[ticks >= 0 & ticks <= tpm_range]
    ticks <- ticks[tpm_x(ticks) + text_w(tpm_label(ticks), fontsize - 2) / 2 <= bars_x + bars_w + 0.1]
    if (length(ticks) < 2 ||
        (ticks[2] - ticks[1]) / tpm_range * bars_w >= max(text_w(tpm_label(ticks), fontsize - 2)) + 0.05) break
  }
  plotSegments(x0 = tpm_x(ticks), y0 = bar_y[1], x1 = tpm_x(ticks), y1 = bars_bottom,
               default.units = "inches", linecolor = "grey92", lwd = 0.75)
  plotSegments(x0 = tpm_x(ticks), y0 = bars_bottom, x1 = tpm_x(ticks), y1 = bars_bottom + 0.04,
               default.units = "inches", linecolor = "grey50", lwd = 0.5)
  text(tpm_label(ticks), tpm_x(ticks), bars_bottom + 0.05, c("center", "top"), size = fontsize - 2)
  text("TPM", bars_x + bars_w / 2, bars_bottom + lay$axis_h, c("center", "bottom"), size = fontsize - 1.5)

  # Bars: each segment starts where the ones before it end.
  left <- do.call(rbind, lapply(seq_len(n), function(i) cumsum(tpm[i, ]) - tpm[i, ]))
  seg_w <- tpm / tpm_range * bars_w
  seg <- which(tpm > 0 & drawn, arr.ind = TRUE)
  if (nrow(seg))
    plotRect(x = tpm_x(left[seg]), y = bar_y[seg[, 1]], width = seg_w[seg], height = lay$bar_h,
             just = c("left", "top"), default.units = "inches",
             fill = transcript_colours[seg[, 2]], linecolor = "white", lwd = 1.3)
  share <- array(sprintf("%.0f%%", 100 * tpm / totals), dim(tpm))
  seg <- which(drawn & tpm >= 0.05 * totals & tpm > 0 &
               seg_w >= text_w(share, fontsize - 1, "bold") + 0.04, arr.ind = TRUE)
  if (nrow(seg))
    text(share[seg], tpm_x(left[seg]) + seg_w[seg] / 2, mid_y[seg[, 1]],
         c("center", "center"), colour = "white", face = "bold")
  text(paste(tpm_label(totals[drawn]), "TPM"), tpm_x(totals[drawn]) + 0.05, mid_y[drawn], c("left", "center"))
  draw_errors(bars_w)
  invisible(NULL)
}
