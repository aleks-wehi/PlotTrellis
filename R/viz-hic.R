info_hic <- function() list(
  kind = "sample", data_type = "hic", order = 20, filter = c(hic = "HiC"),
  scale_bar = TRUE,
  entry_settings = "hic_style"
)

settings_hic <- function(dict, entry, settings, track_w, region_params) {
  if (is.null(entry$height)) {
    sty <- entry$hic_style %||% "rectangle"
    entry$height <- settings$hic_height %||%
      (if (sty == "square") track_w else round(track_w / 2, 2))
  }
  # A colour range with no width (a colour cap of 0) can't be drawn.
  ov <- settings$zrange
  if (!is.null(ov) && length(ov) == 2 && !anyNA(ov) && !(ov[2] > ov[1]))
    .pt_track_stop(if (ov[1] == 0)
                     sprintf("the Hi-C colour cap must be above 0 (it is %s).",
                             format(ov[2]))
                   else sprintf(
                     "the Hi-C colour range must go up (it is %s to %s).",
                     format(ov[1]), format(ov[2])
                   ),
                   entry$height)
  # A track that fails to load keeps its height (render_all() reads `pt_height`
  # off the error), so the figure's layout doesn't change when it fails.
  entry <- tryCatch(load_hic(entry, region_params, settings),
                    error = function(e) {
    e$pt_height <- entry$height
    stop(e)
  })

  # Colour scale (Colour scale setting): Linear, as Juicebox and plotgardener
  # draw, or Log, as HiGlass and cooler do. render_hic() reads it off the loaded
  # entry.
  log_scale <- identical(settings$hic_colour_scale, "log")
  entry$data$log_scale <- log_scale
  counts <- if (isTRUE(entry$data$failed) || is.null(entry$data$matrix))
    numeric(0) else entry$data$matrix$counts
  counts <- counts[is.finite(counts)]
  # A log scale starts at the smallest value above 0 in view (HiGlass's rule),
  # shared like the top.
  if (log_scale && any(counts > 0)) {
    low_key <- .pt_shared_range_key("hic_zlow_", entry)
    dict[[low_key]] <- min(dict[[low_key]], counts[counts > 0])
  }

  key <- .pt_shared_range_key("hic_zrange_", entry)
  ov <- settings$zrange
  if (!is.null(ov) && length(ov) == 2 && !anyNA(ov)) {
    dict[[key]] <- list(override = ov)
    return(list(dict = dict, entry = entry))
  }
  # The colour scale's top when no Colour max is set. Linear: the 95th
  # percentile of the contacts in view, as Juicebox does (the strongest 5% share
  # the darkest colour); the maximum, plotgardener's default, is a few pixels on
  # the diagonal, and the rest of the map came out almost white. Log: the
  # maximum, as HiGlass does (a log scale already spreads the faint contacts
  # out).
  this_max <- if (!length(counts)) 0 else if (log_scale) max(counts) else
    stats::quantile(counts, HIC_AUTO_TOP_QUANTILE, names = FALSE)
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
  # NULL (plotgardener's own default) when there's no positive finite top,
  # not c(0, 0) -- a degenerate zrange makes plotHicRectangle() error out.
  zrange <- if (is.list(raw)) raw$override else
    if (is.finite(raw %||% 0) && (raw %||% 0) > 0)
      c(0, .pt_hic_round_top(raw)) else NULL
  log_scale <- isTRUE(entry$data$log_scale)
  # A log scale can't start at 0: it starts at the smallest value above 0 in
  # view, rounded down (a thousandth of the top when that is not below the top,
  # e.g. a Colour max under every value).
  if (log_scale && !is.null(zrange) && zrange[1] <= 0) {
    low <- dict[[.pt_shared_range_key("hic_zlow_", entry)]]
    zrange[1] <- if (!is.null(low) && low < zrange[2])
      .pt_hic_round_top(low, up = FALSE) else zrange[2] / 1000
  }
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
    colorTrans = if (log_scale) "log10" else "linear",
    x = ctx$track_x, y = ctx$y, width = ctx$track_w, height = ctx$height
  )
  # No contacts on the chromosome (load_hic()): drawn from the empty matrix, as
  # the file would fail to read; plotgardener's warnings about the empty input
  # are its own, not the user's.
  hic_plot <- if (isTRUE(entry$data$no_contacts))
    suppressWarnings(draw(entry$data$matrix)) else draw(entry$path)
  # Colour legend laid out like plotgardener's annoHeatmapLegend() (max label,
  # bar, min label; same 100 colours), but drawn as vector rects: that one is a
  # raster. There is no colour range to show only when no Hi-C track of this
  # name has counts in view (and no colour range is set): the track is then its
  # heading alone.
  if (!is.null(hic_plot$zrange) && !is.null(hic_plot$color_palette)) {
    leg_x <- ctx$track_x + ctx$track_w + 0.15
    leg_h <- min(ctx$height * 0.6, 1.5)
    lab_h <- grid::convertHeight(
      grid::grobHeight(grid::textGrob("0", gp = grid::gpar(
        fontsize = ctx$fontsize, fontfamily = ctx$font
      ))),
      "inches", valueOnly = TRUE
    )
    plotText(format(signif(hic_plot$zrange[2], 3), big.mark = ",",
                    scientific = FALSE),
             x = leg_x + 0.1, y = ctx$y,
             just = c("center", "top"), fontsize = ctx$fontsize,
             fontcolor = ctx$fontcol, fontfamily = ctx$font)
    plotText(format(signif(hic_plot$zrange[1], 3), big.mark = ",",
                    scientific = FALSE),
             x = leg_x + 0.1, y = ctx$y + leg_h,
             just = c("center", "bottom"), fontsize = ctx$fontsize,
             fontcolor = ctx$fontcol, fontfamily = ctx$font)
    .pt_gradient_bar(rev(mapColors(seq(hic_plot$zrange[1], hic_plot$zrange[2],
                                       length.out = 100),
                                   hic_plot$color_palette)),
                     x = leg_x, y = ctx$y + 1.5 * lab_h, width = 0.2,
                     height = leg_h - 3 * lab_h)
  }
  # The bin size drawn, under the colour key: it changes with the view
  # (Resolution setting) and is what one pixel of the map stands for. A log
  # colour scale says so below it.
  drawn_res <- hic_plot$resolution %||% entry$data$resolution
  if (!is.null(drawn_res) && is.finite(drawn_res)) {
    notes <- c(paste(.pt_hic_bin_label(drawn_res), "bins"),
               if (log_scale) "log scale")
    for (k in seq_along(notes))
      plotText(notes[k], x = ctx$track_x + ctx$track_w + 0.25,
               y = ctx$y + min(ctx$height * 0.6, 1.5) + 0.05 + (k - 1) * 0.16,
               just = c("center", "top"), fontsize = ctx$fontsize - 1,
               fontcolor = ctx$fontcol, fontfamily = ctx$font)
  }
  .pt_heading(ctx, entry)
  if (ctx$is_first) draw_scale_bar(ctx)
}

# Quantile of the contacts in view that tops the colour scale when no Colour max
# is set.
HIC_AUTO_TOP_QUANTILE <- 0.95

# An end of the colour scale as drawn and written into the script: 4 significant
# digits, rounded up for the top so it is never below the value it stands for (4
# decimals put 0.00012345 at 0.0001), down for a log scale's bottom. The colour
# key labels it to 3 significant digits.
.pt_hic_round_top <- function(x, up = TRUE) {
  step <- 10^(floor(log10(x)) - 3)
  r <- (if (up) ceiling else floor)(signif(x / step, 12)) * step
  signif(r, 4)
}

# A bin size as people write it: 500 bp, 5 kb, 2.5 Mb.
.pt_hic_bin_label <- function(bp) {
  f <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
  if (bp >= 1e6) paste(f(bp / 1e6), "Mb") else
    if (bp >= 1e3) paste(f(bp / 1e3), "kb") else paste(f(bp), "bp")
}
