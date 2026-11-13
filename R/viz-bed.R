info_bed <- function() list(kind = "sample", data_type = "bed", order = 30, filter = c(bed = "BED"))

settings_bed <- function(dict, entry, settings, track_w, region_params) {
  if (is.null(entry$height)) entry$height <- 0.32
  list(dict = dict, entry = entry)
}

render_bed <- function(entry, dict, ctx) {
  entry <- load_bed(entry, ctx$region_params)
  if (isTRUE(entry$data$no_chrom)) {
    .pt_status_text(ctx, "Chromosome not found")
    return(invisible(NULL))
  }
  line_h <- 0.16
  line_y <- ctx$y + ctx$height - line_h
  if (!isTRUE(entry$data$empty)) plotRanges(
    data = entry$path, params = ctx$region,
    chrom = entry$data$chrom, chromstart = ctx$region_params$start, chromend = ctx$region_params$end,
    x = ctx$track_x, y = line_y, width = ctx$track_w, height = line_h,
    fill = entry$colour, linecolor = NA, boxHeight = grid::unit(1.2, "mm"), collapse = TRUE
  )
  mixed_text(
    list(list(text = entry$name, colour = ctx$fontcol, fontface = "bold"),
         list(text = paste0(" \u00B7 ", entry$sample), colour = entry$colour)),
    x = ctx$track_x + 0.05, y = ctx$y + 0.02, fontsize = ctx$fontsize - 1, fontfamily = ctx$font
  )
  invisible(NULL)
}
