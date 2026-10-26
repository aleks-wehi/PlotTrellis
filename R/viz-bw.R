settings_bw <- function(dict, entry, settings, track_w, region_params) {
  entry$data <- load_bw_file(entry$path, region_params)
  if (is.null(entry$height)) entry$height <- settings$atac_height %||% 0.85

  key <- paste0("bw_range_", gsub("[^A-Za-z0-9._]", "_", entry$name))
  this_max <- if (length(entry$data$gr) > 0) max(entry$data$gr$score, na.rm = TRUE) else 0
  cur <- dict[[key]]
  dict[[key]] <- if (is.null(cur)) this_max else max(cur, this_max)
  list(dict = dict, entry = entry)
}

render_bw <- function(entry, dict, ctx) {
  key   <- paste0("bw_range_", gsub("[^A-Za-z0-9._]", "_", entry$name))  # must match settings_bw's key
  raw   <- dict[[key]]
  range <- if (!is.null(raw) && is.finite(raw) && raw > 0) c(0, round(raw * (10 / 9), 4)) else NULL

  if (isTRUE(entry$data$failed)) {
    label <- if (isTRUE(entry$data$no_chrom)) "Chromosome not found" else "No data"
    no_data_text(ctx$track_x, ctx$y, ctx$track_w, ctx$height, ctx$fontsize - 1, ctx$font, ctx$fontcol, label = label)
  } else {
    plotSignal(
      data = entry$path, params = ctx$region,
      chrom = entry$data$chrom, chromstart = ctx$region_params$start, chromend = ctx$region_params$end,
      x = ctx$track_x, y = ctx$y, width = ctx$track_w, height = ctx$height,
      fill = entry$colour, linecolor = entry$colour, range = range
    )
    mixed_text(
      list(list(text = entry$name, colour = ctx$fontcol, fontface = "bold"),
           list(text = paste0(" \u00B7 ", entry$sample), colour = entry$colour)),
      x = ctx$track_x + 0.05, y = ctx$y + 0.08, fontsize = ctx$fontsize - 1, fontfamily = ctx$font
    )
    if (isTRUE(ctx$show_track_axes) && !is.null(range)) {
      draw_value_axis(range, ctx$track_x - 0.08, ctx$y, ctx$height, ctx$fontsize - 2, ctx$font, ctx$fontcol)
    }
  }
  if (ctx$is_first) draw_scale_bar(ctx)
}
