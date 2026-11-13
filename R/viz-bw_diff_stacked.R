info_bw_diff_stacked <- function() list(
  kind = "comparison", data_type = "bw", order = 10, label = "{name} diff (stacked)", api = "{name}_diff",
  filter = c(comparison = "Comparisons"), scale_bar = TRUE
)

settings_bw_diff_stacked <- function(dict, entry, settings, track_w, region_params) {
  if (is.null(entry$height)) entry$height <- settings$atac_height
  # Project-wide setting, threaded through the entry since ctx has no field for it: how strongly
  # the shared (overlap) layer is tinted with the reference colour, 0 = white, 1 = full colour.
  entry$shared_alpha <- settings$diff_shared_alpha
  list(dict = dict, entry = entry)
}

render_bw_diff_stacked <- function(entry, dict, ctx) {
  ctrl  <- load_bw_file(entry$control_path,  ctx$region_params)
  treat <- load_bw_file(entry$treatment_path, ctx$region_params)

  if (isTRUE(ctrl$no_chrom) || isTRUE(treat$no_chrom) || ctrl$chrom != treat$chrom) {
    .pt_status_text(ctx, "Chromosome not found")
    if (ctx$is_first) draw_scale_bar(ctx)
    return(invisible(NULL))
  }
  entry$data <- list(ctrl = ctrl$gr, treat = treat$gr, chrom = ctrl$chrom)

  gmax <- max(
    if (length(entry$data$treat) > 0) max(entry$data$treat$score, na.rm = TRUE) else 0,
    if (length(entry$data$ctrl)  > 0) max(entry$data$ctrl$score,  na.rm = TRUE) else 0
  )
  # Signal above 0: 0 to its maximum. None above 0 and none below (no signal in view, or only
  # explicit zeros): an empty 0-1 scale, axis included. Some below 0 and none above:
  # plotgardener's own range (NULL), which draws it.
  negative <- any(entry$data$treat$score < 0, entry$data$ctrl$score < 0, na.rm = TRUE)
  entry$range <- if (is.finite(gmax) && gmax > 0) c(0, gmax * (10 / 9)) else if (!negative) c(0, 1)

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
    rle_ctrl  <- tryCatch(rtracklayer::import(entry$control_path,  which = which_gr, as = "RleList")[ch],
                          error = function(e) .pt_read_error(entry$control_path, e))
    rle_treat <- tryCatch(rtracklayer::import(entry$treatment_path, which = which_gr, as = "RleList")[ch],
                          error = function(e) .pt_read_error(entry$treatment_path, e))
    if (!is.null(rle_ctrl) && !is.null(rle_treat)) {
      g <- as(pmin(rle_ctrl, rle_treat), "GRanges")
      g[GenomicRanges::start(g) < ctx$region_params$end & GenomicRanges::end(g) > ctx$region_params$start]
    } else GRanges()
  }
  plotSignal(data = min_gr, params = ctx$region,
             chrom = ch, chromstart = ctx$region_params$start, chromend = ctx$region_params$end,
             x = ctx$track_x, y = ctx$y,
             width = ctx$track_w, height = ctx$height,
             fill = faint_solid(entry$control_colour, entry$shared_alpha),
             linecolor = faint_solid(entry$control_colour, entry$shared_alpha),
             range = entry$range)

  .pt_heading_pair(ctx, entry, paste0(entry$name, " diff (stacked)"))
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
