settings_transcripts_annotated <- function(dict, entry, settings, track_w, region_params) {
  entry <- load_transcripts(entry, region_params, settings)
  tx_all <- entry$data$tx_all
  kept <- tx_all[tx_all$transcript_name %in% entry$data$named, ]
  kept <- kept[order(kept$transcript_name), ]
  entry$data$kept <- kept
  entry$data$n_excluded <- nrow(tx_all) - nrow(kept)
  entry$data$colours <- transcript_colours(kept$transcript_name)
  if (is.null(entry$height)) entry$height <- 2.0
  list(dict = dict, entry = entry)
}

render_transcripts_annotated <- function(entry, dict, ctx) {
  plotText(label = "Transcripts annotated (filtered)", x = ctx$track_x + 0.05, y = ctx$y + 0.02,
           just = c("left", "top"), fontsize = ctx$fontsize - 1, fontface = "bold",
           fontfamily = ctx$font, fontcolor = ctx$fontcol)

  if (entry$data$n_excluded > 0) {
    header_text_h <- grid::convertHeight(grid::grobHeight(grid::textGrob("Xg", gp = grid::gpar(
      fontsize = ctx$fontsize - 1, fontface = "bold", fontfamily = ctx$font))), "inches", valueOnly = TRUE)
    plotText(label = sprintf("+%d not shown", entry$data$n_excluded), x = ctx$track_x + ctx$track_w,
             y = ctx$y + 0.02 + header_text_h / 2, just = c("right", "center"),
             fontsize = ctx$fontsize - 2, fontfamily = ctx$font, fontcolor = "grey50")
  }

  header_h <- 0.2
  kept <- entry$data$kept
  n_rows <- nrow(kept)

  if (n_rows == 0) {
    plotText(label = "No retained transcripts for this gene in view", x = ctx$track_x + ctx$track_w / 2,
             y = ctx$y + header_h + (ctx$height - header_h) / 2, just = "center",
             fontsize = ctx$fontsize - 1, fontfamily = ctx$font, fontcolor = ctx$fontcol)
    return(invisible(NULL))
  }

  row_h <- (ctx$height - header_h) / n_rows
  row_y <- ctx$y + header_h + (seq_len(n_rows) - 1) * row_h
  label_dy <- row_h - grid::convertHeight(grid::unit(1, "mm"), "inches", valueOnly = TRUE)

  # coding_status -> (text colour, tag background); "coding" draws no tag.
  status_colour <- c(NMD = "grey40", "no CDS" = "grey40", "non-coding" = "grey40", unknown = "#c76e00")
  status_bg     <- c(NMD = "#e8e8e8", "no CDS" = "#e8e8e8", "non-coding" = "#e8e8e8", unknown = "#fbe4c8")
  swatch_size <- 0.09
  tag_pad_x   <- 0.045
  tag_pad_y   <- 0.03

  txn_text_w <- function(label, fontsize, fontfamily, fontface = "plain") {
    grid::convertWidth(grid::grobWidth(grid::textGrob(label, gp = grid::gpar(fontsize = fontsize, fontfamily = fontfamily, fontface = fontface))), "inches", valueOnly = TRUE)
  }
  txn_text_h <- function(label, fontsize, fontfamily, fontface = "plain") {
    grid::convertHeight(grid::grobHeight(grid::textGrob(label, gp = grid::gpar(fontsize = fontsize, fontfamily = fontfamily, fontface = fontface))), "inches", valueOnly = TRUE)
  }

  colours <- entry$data$colours
  for (j in seq_len(n_rows)) {
    plotTranscripts(params = ctx$region, transcriptFilter = kept$transcript_id[j], labels = NULL,
                     x = ctx$track_x, y = row_y[j], width = ctx$track_w, height = row_h,
                     boxHeight = grid::unit(2, "mm"), fill = colours[[kept$transcript_name[j]]], colorbyStrand = FALSE)

    row_yc <- row_y[j] + label_dy
    status <- kept$coding_status[j]
    name   <- kept$transcript_name[j]
    face   <- if (kept$is_canonical[j]) "bold" else "plain"

    if (status != "coding") {
      tag_w <- txn_text_w(status, ctx$fontsize - 1.5, ctx$font) + 2 * tag_pad_x
      tag_h <- txn_text_h("Xg", ctx$fontsize - 1.5, ctx$font) + 2 * tag_pad_y
      tag_x <- ctx$track_x - 0.05 - tag_w
      plotRect(x = tag_x, y = row_yc, width = tag_w, height = tag_h, just = c("left", "center"),
               default.units = "inches", fill = unname(status_bg[status]), linecolor = NA)
      plotText(label = status, x = tag_x + tag_pad_x, y = row_yc, just = c("left", "center"),
               fontsize = ctx$fontsize - 1.5, fontfamily = ctx$font,
               fontface = if (status == "unknown") "bold" else "plain", fontcolor = unname(status_colour[status]))
    }

    cur_x <- ctx$track_x + ctx$track_w + 0.05
    plotRect(x = cur_x, y = row_yc, width = swatch_size, height = swatch_size, just = c("left", "center"),
             default.units = "inches", fill = colours[[name]], linecolor = NA)
    cur_x <- cur_x + swatch_size + 0.06

    plotText(label = name, x = cur_x, y = row_yc, just = c("left", "center"),
             fontsize = ctx$fontsize - 1.5, fontfamily = ctx$font, fontface = face, fontcolor = ctx$fontcol)
  }
  invisible(NULL)
}
