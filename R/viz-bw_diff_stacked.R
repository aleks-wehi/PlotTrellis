settings_bw_diff_stacked <- function(dict, entry, settings, track_w, region_params) {
  if (is.null(entry$height)) entry$height <- settings$bw_diff_stacked_height %||% settings$atac_height %||% 0.85
  list(dict = dict, entry = entry)
}

render_bw_diff_stacked <- function(entry, dict, ctx) {
  ctrl  <- load_bw_file(entry$control_path,  ctx$region_params)
  treat <- load_bw_file(entry$treatment_path, ctx$region_params)

  if (isTRUE(ctrl$no_chrom) || isTRUE(treat$no_chrom) || ctrl$chrom != treat$chrom) {
    entry$data <- list(ctrl = GRanges(), treat = GRanges(), failed = TRUE, no_chrom = TRUE, chrom = NA_character_)
  } else {
    entry$data <- list(
      ctrl = ctrl$gr, treat = treat$gr, failed = ctrl$failed || treat$failed,
      no_chrom = FALSE, chrom = ctrl$chrom
    )
  }

  gmax <- max(
    if (length(entry$data$treat) > 0) max(entry$data$treat$score, na.rm = TRUE) else 0,
    if (length(entry$data$ctrl)  > 0) max(entry$data$ctrl$score,  na.rm = TRUE) else 0
  )
  entry$range <- if (is.finite(gmax) && gmax > 0) c(0, gmax * (10 / 9)) else NULL

  if (isTRUE(entry$data$failed)) {
    label <- if (isTRUE(entry$data$no_chrom)) "Chromosome not found" else "No data"
    no_data_text(ctx$track_x, ctx$y, ctx$track_w, ctx$height, ctx$fontsize - 1, ctx$font, ctx$fontcol, label = label)
    if (ctx$is_first) draw_scale_bar(ctx)
    return(invisible(NULL))
  }

  # Pre-blends a colour with the page background at `alpha` into an OPAQUE
  # hex, so the top cap layer fully occludes the two layers under it.
  faint_solid <- function(colour, alpha, bg = "#FFFFFF") {
    cf  <- grDevices::col2rgb(colour) / 255
    bgf <- grDevices::col2rgb(bg) / 255
    m   <- alpha * cf + (1 - alpha) * bgf
    grDevices::rgb(m[1], m[2], m[3])
  }

  hi <- ctx$diff_colour_hi
  lo <- ctx$diff_colour_lo
  ch <- entry$data$chrom

  # Layer 1 (bottom): treatment, full height, GAIN colour.
  plotSignal(data = entry$treatment_path, params = ctx$region,
             chrom = ch, chromstart = ctx$region_params$start, chromend = ctx$region_params$end,
             x = ctx$track_x, y = ctx$y,
             width = ctx$track_w, height = ctx$height, fill = hi, linecolor = hi, range = entry$range)
  # Layer 2: control, full height, LOSS colour -- occludes layer 1 wherever control > treatment.
  plotSignal(data = entry$control_path, params = ctx$region,
             chrom = ch, chromstart = ctx$region_params$start, chromend = ctx$region_params$end,
             x = ctx$track_x, y = ctx$y,
             width = ctx$track_w, height = ctx$height, fill = lo, linecolor = lo, range = entry$range)

  # Layer 3 (top): min(control, treatment), opaque faint reference colour --
  # always exactly as tall as whichever of layers 1/2 must be fully erased.
  which_gr <- GRanges(ch, IRanges(ctx$region_params$start, ctx$region_params$end))
  min_gr <- {
    # Single-bracket `[ch]`, not `[[ch]]` -- keeps a length-1 RleList so
    # as(pmin(...), "GRanges") can recover seqnames (no method for a bare Rle).
    rle_ctrl  <- tryCatch(rtracklayer::import(entry$control_path,  which = which_gr, as = "RleList")[ch], error = function(e) NULL)
    rle_treat <- tryCatch(rtracklayer::import(entry$treatment_path, which = which_gr, as = "RleList")[ch], error = function(e) NULL)
    if (!is.null(rle_ctrl) && !is.null(rle_treat)) {
      g <- as(pmin(rle_ctrl, rle_treat), "GRanges")
      g[GenomicRanges::start(g) < ctx$region_params$end & GenomicRanges::end(g) > ctx$region_params$start]
    } else GRanges()
  }
  plotSignal(data = min_gr, params = ctx$region,
             chrom = ch, chromstart = ctx$region_params$start, chromend = ctx$region_params$end,
             x = ctx$track_x, y = ctx$y,
             width = ctx$track_w, height = ctx$height,
             fill = faint_solid(entry$control_colour, 0.15), linecolor = faint_solid(entry$control_colour, 0.15),
             range = entry$range)

  mixed_text(
    list(list(text = paste0(entry$name, " diff (stacked)"), colour = ctx$fontcol, fontface = "bold"),
         list(text = " · ", colour = ctx$fontcol),
         list(text = entry$treatment_label, colour = entry$treatment_colour),
         list(text = " (Cmp) vs ", colour = ctx$fontcol),
         list(text = entry$control_label, colour = entry$control_colour),
         list(text = " (Ref)", colour = ctx$fontcol)),
    x = ctx$track_x + 0.05, y = ctx$y + 0.08, fontsize = ctx$fontsize - 1, fontfamily = ctx$font
  )
  if (isTRUE(ctx$show_track_axes) && !is.null(entry$range)) {
    draw_value_axis(entry$range, ctx$track_x - 0.08, ctx$y, ctx$height, ctx$fontsize - 2, ctx$font, ctx$fontcol)
  }
  plotLegend(
    legend = c("Gain", "Loss"), fill = c(hi, lo), border = FALSE, bg = NA,
    x = ctx$track_x + ctx$track_w + 0.1, y = ctx$y, width = 0.8, height = 0.45,
    just = c("left", "top"), default.units = "inches", fontsize = ctx$fontsize, fontfamily = ctx$font
  )
  if (ctx$is_first) draw_scale_bar(ctx)
}
