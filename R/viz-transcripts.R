settings_transcripts <- function(dict, entry, settings, track_w, region_params) {
  entry <- load_transcripts(entry, region_params, settings)
  tx_all <- entry$data$tx_all
  kept <- tx_all[tx_all$transcript_name %in% entry$data$named, ]
  kept <- kept[order(kept$transcript_name), ]
  entry$data$kept <- kept
  entry$data$n_excluded <- nrow(tx_all) - nrow(kept)
  entry$data$colours <- transcript_colours(kept$transcript_name)
  if (is.null(entry$height)) entry$height <- 1.6
  list(dict = dict, entry = entry)
}

render_transcripts <- function(entry, dict, ctx) {
  plotText(label = "Transcripts (filtered)", x = ctx$track_x + 0.05, y = ctx$y + 0.02, just = c("left", "top"),
           fontsize = ctx$fontsize - 1, fontface = "bold", fontfamily = ctx$font, fontcolor = ctx$fontcol)

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
  colours <- entry$data$colours

  for (j in seq_len(n_rows)) {
    plotTranscripts(params = ctx$region, transcriptFilter = kept$transcript_id[j], labels = NULL,
                     x = ctx$track_x, y = row_y[j], width = ctx$track_w, height = row_h,
                     boxHeight = grid::unit(2, "mm"), fill = colours[[kept$transcript_name[j]]], colorbyStrand = FALSE)
    plotText(label = kept$transcript_name[j], x = ctx$track_x + ctx$track_w + 0.05, y = row_y[j] + label_dy,
             just = c("left", "center"), fontsize = ctx$fontsize - 1.5, fontfamily = ctx$font, fontcolor = ctx$fontcol)
  }
  invisible(NULL)
}
