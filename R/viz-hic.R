info_hic <- function() list(
  kind = "sample", data_type = "hic", order = 20, filter = c(hic = "HiC"),
  scale_bar = TRUE,
  entry_settings = "hic_style"
)

settings_hic <- function(dict, entry, settings, track_w, region_params) {
  sty <- entry$hic_style %||% "rectangle"
  if (is.null(entry$height)) {
    entry$height <- settings$hic_height %||%
      (if (sty == "square") track_w else round(track_w / 2, 2))
  }
  # A Colour max that is not one (a script's typo), or a number of 0 or below,
  # can't be drawn.
  cap <- .pt_hic_zrange_cap(settings$zrange_cap)
  if (length(cap) == 1 && is.na(cap))
    .pt_track_stop(paste0("the Hi-C colour max must be \"max\", \"p95\", ",
                          "\"p90\" or a number (it is ",
                          paste(vapply(settings$zrange_cap, format, ""),
                                collapse = ", "),
                          ")."), entry$height)
  if (is.numeric(cap) && !(cap > 0))
    .pt_track_stop(sprintf("the Hi-C colour cap must be above 0 (it is %s).",
                           format(cap)), entry$height)
  # Rectangle draws a strip as high as the track around the view, so it reads
  # the contacts out to that far beyond each end of the view (as plotgardener
  # does to draw it): its height over its width, in view widths.
  ratio <- entry$height / track_w
  pad <- if (sty == "rectangle")
    ratio * (region_params$end - region_params$start) else 0
  # A track that fails to load keeps its height (render_all() reads `pt_height`
  # off the error), so the figure's layout doesn't change when it fails.
  entry <- tryCatch(load_hic(entry, region_params, settings, pad),
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
    numeric(0) else .pt_hic_drawn_counts(entry$data$matrix, sty,
                                         region_params$start,
                                         region_params$end,
                                         entry$data$resolution, ratio)
  counts <- counts[is.finite(counts)]
  # A log scale starts at the smallest value above 0 drawn (HiGlass's rule),
  # shared like the top.
  if (log_scale && any(counts > 0)) {
    low_key <- .pt_shared_range_key("hic_zlow_", entry)
    dict[[low_key]] <- min(dict[[low_key]], counts[counts > 0])
  }

  key <- .pt_shared_range_key("hic_zrange_", entry)
  if (is.numeric(cap)) {
    dict[[key]] <- list(override = c(0, cap))
    return(list(dict = dict, entry = entry))
  }
  # The colour scale's top, from the contacts drawn (Colour max setting). Auto
  # (NULL): on Linear the 95th percentile, as Juicebox does (the strongest 5%
  # share the darkest colour; the maximum, plotgardener's default, is a few
  # pixels on the diagonal, and the rest of the map came out almost white); on
  # Log the maximum, as HiGlass does (a log scale already spreads the faint
  # contacts out).
  stat <- if (is.character(cap)) cap else if (log_scale) "max" else "p95"
  this_max <- if (!length(counts)) 0 else if (stat == "max") max(counts) else
    stats::quantile(counts, HIC_TOP_QUANTILES[[stat]], names = FALSE)
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

# The Colour max choices that are a quantile of the contacts drawn.
HIC_TOP_QUANTILES <- c(p95 = 0.95, p90 = 0.90)

# The Hi-C Colour max setting (zrange_cap) as a track reads it: NULL for Auto
# (also blank, NA, "auto", or {} from an older settings file), "max", "p95" or
# "p90", or a number (numeric text, as samples.txt and the Colour max box give
# it, is read as one). `bad` for anything else.
.pt_hic_zrange_cap <- function(v, bad = NA) {
  if (is.null(v) || is.list(v)) return(NULL)
  if (length(v) != 1) return(bad)
  if (is.na(v)) return(NULL)
  if (is.numeric(v)) return(if (is.finite(v)) v else bad)
  v <- tolower(trimws(as.character(v)))
  if (v %in% c("", "auto", "na")) return(NULL)
  if (v %in% c("max", names(HIC_TOP_QUANTILES))) return(v)
  n <- suppressWarnings(as.numeric(v))
  if (is.finite(n)) n else bad
}

# The values of the contacts a Hi-C track draws, from the matrix load_hic()
# read (x <= y: the upper half): Rectangle shows those within `ratio` (its
# height over its width) view widths of the diagonal whose midpoint is in the
# view, from a matrix read that far beyond the view; Square every contact of
# the view's square, those off the diagonal twice (it draws both halves);
# Triangle the square's upper half. A bin counts as drawn when its centre is.
.pt_hic_drawn_counts <- function(m, style, start, end, res, ratio) {
  if (!nrow(m)) return(numeric(0))
  if (identical(style, "square")) return(c(m$counts, m$counts[m$x != m$y]))
  if (identical(style, "rectangle")) {
    mid <- (m$x + m$y + res) / 2
    d <- m$y - m$x
    return(m$counts[mid >= start & mid <= end & d >= 0 &
                      d <= 2 * ratio * (end - start)])
  }
  m$counts
}

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

# The top of each Hi-C colour scale a figure drew, from render_all()'s `dict`,
# named by its shared key (.pt_shared_range_key()), rounded as drawn. A Colour
# max typed as a number is left out: it is its own top.
.pt_hic_colour_tops <- function(dict) {
  v <- Filter(function(x) is.numeric(x) && length(x) == 1 && is.finite(x) &&
                x > 0, dict[grepl("^hic_zrange_", names(dict))])
  vapply(v, .pt_hic_round_top, numeric(1))
}
