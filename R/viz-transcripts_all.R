settings_transcripts_all <- function(dict, entry, settings, track_w, region_params) {
  entry <- load_transcripts(entry, region_params, settings)
  tx <- entry$data$tx_all[order(entry$data$tx_all$transcript_name), ]
  entry$data$tx <- tx
  entry$data$colours <- transcript_colours(tx$transcript_name)
  if (is.null(entry$height)) {
    n <- nrow(tx)
    entry$height <- if (n > 0) 0.2 + n * 0.16 else 2.4  # row_h below must match this 0.16
  }
  list(dict = dict, entry = entry)
}

render_transcripts_all <- function(entry, dict, ctx) {
  plotText(label = "Transcripts", x = ctx$track_x + 0.05, y = ctx$y + 0.02, just = c("left", "top"),
           fontsize = ctx$fontsize - 1, fontface = "bold", fontfamily = ctx$font, fontcolor = ctx$fontcol)
  header_h <- 0.2
  tx <- entry$data$tx
  n_rows <- nrow(tx)

  if (n_rows == 0) {
    plotText(label = "No annotated transcripts for this gene in view", x = ctx$track_x + ctx$track_w / 2,
             y = ctx$y + header_h + (ctx$height - header_h) / 2, just = "center",
             fontsize = ctx$fontsize - 1, fontfamily = ctx$font, fontcolor = ctx$fontcol)
    return(invisible(NULL))
  }

  row_h <- 0.16  # must match settings_transcripts_all's height formula
  row_y <- ctx$y + header_h + (seq_len(n_rows) - 1) * row_h
  label_dy <- row_h - grid::convertHeight(grid::unit(1, "mm"), "inches", valueOnly = TRUE)
  colours <- entry$data$colours

  for (j in seq_len(n_rows)) {
    plotTranscripts(params = ctx$region, transcriptFilter = tx$transcript_id[j], labels = NULL,
                     x = ctx$track_x, y = row_y[j], width = ctx$track_w, height = row_h,
                     boxHeight = grid::unit(2, "mm"), fill = colours[[tx$transcript_name[j]]], colorbyStrand = FALSE)
    plotText(label = tx$transcript_name[j], x = ctx$track_x + ctx$track_w + 0.05, y = row_y[j] + label_dy,
             just = c("left", "center"), fontsize = ctx$fontsize - 1.5, fontfamily = ctx$font, fontcolor = ctx$fontcol)
  }
  invisible(NULL)
}
