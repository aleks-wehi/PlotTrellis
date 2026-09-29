settings_genes <- function(dict, entry, settings, track_w, region_params) {
  if (is.null(entry$height)) entry$height <- 0.65
  list(dict = dict, entry = entry)
}

render_genes <- function(entry, dict, ctx) {
  plotGenes(params = ctx$region, x = ctx$track_x, y = ctx$y, width = ctx$track_w, height = ctx$height, fontsize = ctx$fontsize)
  invisible(NULL)
}
