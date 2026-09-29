settings_genome_axis <- function(dict, entry, settings, track_w, region_params) {
  entry$height <- GENOME_AXIS_RESERVE
  list(dict = dict, entry = entry)
}

render_genome_axis <- function(entry, dict, ctx) {
  axis_y <- ctx$y
  tick_h <- 0.04
  plotSegments(ctx$track_x, axis_y, ctx$track_x + ctx$track_w, axis_y, default.units = "inches", lwd = 1, linecolor = ctx$fontcol)
  plotSegments(ctx$track_x, axis_y, ctx$track_x, axis_y + tick_h, default.units = "inches", lwd = 1, linecolor = ctx$fontcol)
  plotSegments(ctx$track_x + ctx$track_w, axis_y, ctx$track_x + ctx$track_w, axis_y + tick_h, default.units = "inches", lwd = 1, linecolor = ctx$fontcol)

  rp <- ctx$region_params
  chrom_label <- switch(ctx$chr_fmt,
    dot   = sub("^chr", "chr.", rp$chrom),
    space = sub("^chr", "chr ", rp$chrom),
    rp$chrom)
  plotText(format_bp_generic(rp$start, 3, space_before_bp = TRUE), x = ctx$track_x, y = axis_y + tick_h + 0.02,
           just = c("left", "top"), fontsize = ctx$fontsize, fontfamily = ctx$font, fontcolor = ctx$fontcol)
  plotText(chrom_label, x = ctx$track_x + ctx$track_w / 2, y = axis_y + tick_h + 0.02,
           just = c("center", "top"), fontsize = ctx$fontsize, fontfamily = ctx$font, fontcolor = ctx$fontcol)
  plotText(format_bp_generic(rp$end, 3, space_before_bp = TRUE), x = ctx$track_x + ctx$track_w, y = axis_y + tick_h + 0.02,
           just = c("right", "top"), fontsize = ctx$fontsize, fontfamily = ctx$font, fontcolor = ctx$fontcol)
  invisible(NULL)
}
