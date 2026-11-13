info_genes <- function() list(
  kind = "fixed", order = 10, label = "Gene annotation", filter = c(annotation = "Annotations")
)

settings_genes <- function(dict, entry, settings, track_w, region_params) {
  if (is.null(entry$height)) entry$height <- 0.65
  list(dict = dict, entry = entry)
}

render_genes <- function(entry, dict, ctx) {
  # plotGenes() has no font argument; its labels take the family of the viewport
  # they are drawn in, and it fits them using that same font.
  grid::pushViewport(grid::viewport(gp = grid::gpar(fontfamily = ctx$font)))
  on.exit(grid::upViewport())
  plotGenes(params = ctx$region, x = ctx$track_x, y = ctx$y, width = ctx$track_w, height = ctx$height, fontsize = ctx$fontsize)
  invisible(NULL)
}
