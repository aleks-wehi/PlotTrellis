info_bw <- function() list(
  kind = "sample", data_type = "bw", order = 10, filter = c(bw = "Coverage"), scale_bar = TRUE
)

settings_bw <- function(dict, entry, settings, track_w, region_params) {
  entry$data <- load_bw_file(entry$path, region_params)
  if (is.null(entry$height)) entry$height <- settings$atac_height

  key <- .pt_shared_range_key("bw_range_", entry)
  # -Inf: no signal in view at all. `bw_negative_`: some track sharing the range has a score
  # below 0 in view (render_bw() then keeps plotgardener's own range rather than 0-1).
  scores <- entry$data$gr$score
  this_max <- if (length(scores) > 0) max(scores, na.rm = TRUE) else -Inf
  cur <- dict[[key]]
  dict[[key]] <- if (is.null(cur)) this_max else max(cur, this_max)
  neg_key <- .pt_shared_range_key("bw_negative_", entry)
  dict[[neg_key]] <- isTRUE(dict[[neg_key]]) || any(scores < 0, na.rm = TRUE)
  list(dict = dict, entry = entry)
}

render_bw <- function(entry, dict, ctx) {
  key   <- .pt_shared_range_key("bw_range_", entry)
  raw   <- dict[[key]]
  # Signal above 0: 0 to its shared maximum. None above 0 and none below (no signal in view, or
  # only explicit zeros): an empty 0-1 scale, axis included. Some below 0 and none above (all
  # negative, say): plotgardener's own range (NULL), which draws it.
  negative <- isTRUE(dict[[.pt_shared_range_key("bw_negative_", entry)]])
  range <- if (!is.null(raw) && is.finite(raw) && raw > 0) c(0, round(raw * (10 / 9), 4)) else if (!negative) c(0, 1)

  if (isTRUE(entry$data$no_chrom)) {
    .pt_status_text(ctx, "Chromosome not found")
  } else {
    plotSignal(
      data = entry$path, params = ctx$region,
      chrom = entry$data$chrom, chromstart = ctx$region_params$start, chromend = ctx$region_params$end,
      x = ctx$track_x, y = ctx$y, width = ctx$track_w, height = ctx$height,
      fill = entry$colour, linecolor = entry$colour, range = range
    )
    .pt_heading(ctx, entry)
    if (isTRUE(ctx$show_track_axes) && !is.null(range)) {
      draw_value_axis(range, ctx$track_x - 0.08, ctx$y, ctx$height, ctx$fontsize - 2, ctx$font, ctx$fontcol)
    }
  }
  if (ctx$is_first) draw_scale_bar(ctx)
}
