info_hic <- function() list(
  kind = "sample", data_type = "hic", order = 20, filter = c(hic = "HiC"), scale_bar = TRUE,
  entry_settings = "hic_style"
)

settings_hic <- function(dict, entry, settings, track_w, region_params) {
  entry <- load_hic(entry, region_params, settings)
  if (is.null(entry$height)) {
    sty <- entry$hic_style %||% "rectangle"
    entry$height <- settings$hic_height %||% (if (sty == "square") track_w else round(track_w / 2, 2))
  }

  key <- .pt_shared_range_key("hic_zrange_", entry)
  ov <- settings$zrange
  if (!is.null(ov) && length(ov) == 2 && !anyNA(ov)) {
    dict[[key]] <- list(override = ov)
    return(list(dict = dict, entry = entry))
  }
  this_max <- if (isTRUE(entry$data$failed) || nrow(entry$data$matrix) == 0 || all(is.na(entry$data$matrix$counts))) 0 else max(entry$data$matrix$counts, na.rm = TRUE)
  cur <- dict[[key]]
  cur_max <- if (is.list(cur)) NULL else cur
  dict[[key]] <- if (is.null(cur_max)) this_max else max(cur_max, this_max)
  list(dict = dict, entry = entry)
}

render_hic <- function(entry, dict, ctx) {
  if (isTRUE(entry$data$no_chrom)) {
    .pt_status_text(ctx, "Chromosome not found")
    if (ctx$is_first) draw_scale_bar(ctx)
    return(invisible(NULL))
  }

  key    <- .pt_shared_range_key("hic_zrange_", entry)
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

  draw <- function(data) hic_call(
    data = data, params = ctx$region, chrom = entry$data$chrom,
    resolution = entry$data$resolution, zrange = zrange, norm = ctx$hic_norm,
    x = ctx$track_x, y = ctx$y, width = ctx$track_w, height = ctx$height
  )
  # No contacts on the chromosome (load_hic()): drawn from the empty matrix, as the file would
  # fail to read; plotgardener's warnings about the empty input are its own, not the user's.
  hic_plot <- if (isTRUE(entry$data$no_contacts)) suppressWarnings(draw(entry$data$matrix)) else draw(entry$path)
  # Colour legend laid out like plotgardener's annoHeatmapLegend() (max label, bar,
  # min label; same 100 colours), but drawn as vector rects: that one is a raster.
  # There is no colour range to show only when no Hi-C track of this name has counts in
  # view (and no colour range is set): the track is then its heading alone.
  if (!is.null(hic_plot$zrange) && !is.null(hic_plot$color_palette)) {
    leg_x <- ctx$track_x + ctx$track_w + 0.15
    leg_h <- min(ctx$height * 0.6, 1.5)
    lab_h <- grid::convertHeight(grid::grobHeight(grid::textGrob("0", gp = grid::gpar(
      fontsize = ctx$fontsize, fontfamily = ctx$font))), "inches", valueOnly = TRUE)
    plotText(format(hic_plot$zrange[2], scientific = FALSE, digits = 1), x = leg_x + 0.1, y = ctx$y,
             just = c("center", "top"), fontsize = ctx$fontsize, fontcolor = ctx$fontcol, fontfamily = ctx$font)
    plotText(format(hic_plot$zrange[1], scientific = FALSE, digits = 1), x = leg_x + 0.1, y = ctx$y + leg_h,
             just = c("center", "bottom"), fontsize = ctx$fontsize, fontcolor = ctx$fontcol, fontfamily = ctx$font)
    .pt_gradient_bar(rev(mapColors(seq(hic_plot$zrange[1], hic_plot$zrange[2], length.out = 100), hic_plot$color_palette)),
                     x = leg_x, y = ctx$y + 1.5 * lab_h, width = 0.2, height = leg_h - 3 * lab_h)
  }
  .pt_heading(ctx, entry)
  if (ctx$is_first) draw_scale_bar(ctx)
}
