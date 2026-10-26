settings_hic <- function(dict, entry, settings, track_w, region_params) {
  entry <- load_hic(entry, region_params, settings)
  if (is.null(entry$height)) {
    sty <- entry$hic_style %||% "rectangle"
    entry$height <- settings$hic_height %||% (if (sty == "square") track_w else round(track_w / 2, 2))
  }

  key <- paste0("hic_zrange_", gsub("[^A-Za-z0-9._]", "_", entry$name))
  ov <- settings$zrange
  if (!is.null(ov) && length(ov) == 2 && !anyNA(ov)) {
    dict[[key]] <- list(override = ov)
    return(list(dict = dict, entry = entry))
  }
  this_max <- if (isTRUE(entry$data$failed) || nrow(entry$data$matrix) == 0) 0 else max(entry$data$matrix$counts, na.rm = TRUE)
  cur <- dict[[key]]
  cur_max <- if (is.list(cur)) NULL else cur
  dict[[key]] <- if (is.null(cur_max)) this_max else max(cur_max, this_max)
  list(dict = dict, entry = entry)
}

render_hic <- function(entry, dict, ctx) {
  if (isTRUE(entry$data$no_chrom)) {
    no_data_text(ctx$track_x, ctx$y, ctx$track_w, ctx$height, ctx$fontsize - 1, ctx$font, ctx$fontcol, label = "Chromosome not found")
    if (ctx$is_first) draw_scale_bar(ctx)
    return(invisible(NULL))
  }

  key    <- paste0("hic_zrange_", gsub("[^A-Za-z0-9._]", "_", entry$name))  # must match settings_hic's key
  raw    <- dict[[key]]
  # NULL (plotgardener's own default) when there's no positive finite max,
  # not c(0, 0) -- a degenerate zrange makes plotHicRectangle() error out.
  zrange <- if (is.list(raw)) raw$override else if (is.finite(raw %||% 0) && (raw %||% 0) > 0) c(0, round(raw, 4)) else NULL
  sty    <- entry$hic_style %||% "rectangle"

  # plotHicRectangle/Square/Triangle do eval(match.call()[[1]]) internally, so
  # calling through an aliased variable breaks dispatch -- call by literal name.
  hic_call <- function(...) switch(sty,
    rectangle = plotHicRectangle(...),
    square    = plotHicSquare(...),
    plotHicTriangle(...)
  )

  tryCatch({
    hic_plot <- hic_call(
      data = entry$path, params = ctx$region, chrom = entry$data$chrom,
      resolution = entry$data$resolution, zrange = zrange, norm = ctx$hic_norm,
      x = ctx$track_x, y = ctx$y, width = ctx$track_w, height = ctx$height
    )
    annoHeatmapLegend(
      plot = hic_plot, x = ctx$track_x + ctx$track_w + 0.15, y = ctx$y,
      width = 0.2, height = min(ctx$height * 0.6, 1.5),
      just = c("left", "top"), default.units = "inches",
      fontsize = ctx$fontsize, fontcolor = ctx$fontcol, fontfamily = ctx$font
    )
    mixed_text(
      list(list(text = entry$name, colour = ctx$fontcol, fontface = "bold"),
           list(text = paste0(" \u00B7 ", entry$sample), colour = entry$colour)),
      x = ctx$track_x + 0.05, y = ctx$y + 0.08, fontsize = ctx$fontsize - 1, fontfamily = ctx$font
    )
  }, error = function(e) {
    no_data_text(ctx$track_x, ctx$y, ctx$track_w, ctx$height, ctx$fontsize - 1, ctx$font, ctx$fontcol)
  })
  if (ctx$is_first) draw_scale_bar(ctx)
}
