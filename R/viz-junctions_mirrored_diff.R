info_junctions_mirrored_diff <- function() list(
  kind = "comparison", data_type = "paired-end", order = 20, label = "{name} mirrored diff",
  api = "{name}_diff", filter = c(comparison = "Comparisons")
)

# Mirrored junction comparison: Comparison's arcs above a shared baseline, Reference's
# mirrored below it, each half junction_height tall under the heading band.
settings_junctions_mirrored_diff <- function(dict, entry, settings, track_w, region_params) {
  if (is.null(entry$height)) entry$height <- JUNCTION_HEADING_H + 2 * settings$junction_height
  list(dict = dict, entry = entry)
}

render_junctions_mirrored_diff <- function(entry, dict, ctx) {
  # Read before the call: a read error raised inside merge()'s lazily evaluated argument would
  # reach the track's failure message wrapped in S4 dispatch text.
  treat <- load_junctions(entry$treatment_path, ctx$region_params)
  ctrl  <- load_junctions(entry$control_path,   ctx$region_params)
  arches <- mirrored_junction_arches(
    treat = treat, ctrl = ctrl,
    chrom = ctx$region_params$chrom, min_support = ctx$min_junction_support, cap = ctx$junction_opacity_cap
  )
  render_mirrored_junction_track(arches, " mirrored diff", " (Cmp) vs ", entry, ctx)
}

# A mirrored junction track from its arcs: the arcs, axis and legend below the heading band,
# then the heading "<name><title> . <Comparison><cmp_note><Reference> (Ref)". Shared with the
# TPM-scaled track, which differs in its arcs and heading text only.
render_mirrored_junction_track <- function(arches, title, cmp_note, entry, ctx) {
  draw_mirrored_junctions(
    arches, region = ctx$region, x = ctx$track_x, y = ctx$y + JUNCTION_HEADING_H,
    width = ctx$track_w, height = ctx$height - JUNCTION_HEADING_H,
    treat_colour = entry$treatment_colour, ctrl_colour = entry$control_colour, cap = ctx$junction_opacity_cap,
    show_axis = isTRUE(ctx$show_track_axes), fontsize = ctx$fontsize, fontfamily = ctx$font, fontcolor = ctx$fontcol
  )
  .pt_heading_pair(ctx, entry, paste0(entry$name, title), cmp_note)
  invisible(NULL)
}

# The two samples' junctions (as read by load_junctions(): intron_start, intron_end, nSupport)
# as arcs for draw_mirrored_junctions(): list(treat, ctrl), each junction_arches() plus an
# `opacity` column from 0 to 1.
#
# Opacity is the junction's |log2 fold change| in read count between the samples, full at
# `cap`. The fold change uses every junction read (pseudocount 1; a junction missing from one
# sample counts 0 there); only junctions with at least `min_support` reads in a sample are
# drawn for it. `treat_scale` multiplies Comparison's counts first, in both the arc heights
# and the fold change (1 = raw counts).
mirrored_junction_arches <- function(treat, ctrl, chrom, min_support, cap, treat_scale = 1) {
  merged <- merge(treat, ctrl, by = c("intron_start", "intron_end"), all = TRUE, suffixes = c("_treat", "_ctrl"))
  merged$nSupport_treat[is.na(merged$nSupport_treat)] <- 0
  merged$nSupport_ctrl[is.na(merged$nSupport_ctrl)]   <- 0
  abs_lfc <- abs(log2((merged$nSupport_treat * treat_scale + 1) / (merged$nSupport_ctrl + 1)))

  arches <- function(junctions, scale) {
    junctions <- junctions[junctions$nSupport >= min_support, ]
    a <- junction_arches(junctions, chrom)
    a$nSupport <- a$nSupport * scale
    in_merged <- match(paste(junctions$intron_start, junctions$intron_end),
                       paste(merged$intron_start, merged$intron_end))
    a$opacity <- pmin(abs_lfc[in_merged] / cap, 1)
    a
  }
  list(treat = arches(treat, treat_scale), ctrl = arches(ctrl, 1))
}

# Draws mirrored_junction_arches() on the current plotgardener page: `arches$treat` above the
# baseline and `arches$ctrl` below it, on one shared count scale, each arc in its sample's
# colour at its own opacity; with the count axis on the left (show_axis) and the opacity
# legend to the right of the track. (x, y) is the top left corner of the arcs in inches and
# `height` covers both halves. `region` is the pgParams() genomic region.
draw_mirrored_junctions <- function(arches, region, x, y, width, height, treat_colour, ctrl_colour, cap,
                                    show_axis = TRUE, fontsize = 8, fontfamily = "", fontcolor = "black") {
  half_h <- height / 2
  mid_y  <- y + half_h
  bottom <- y + height
  max_count <- max(arches$treat$nSupport, arches$ctrl$nSupport, 0)
  if (max_count <= 0) max_count <- 1

  half <- function(a, colour, flip, top) {
    # Guarded on nrow() > 0: plotPairsArches errors on zero-row data when fill/linecolor are vectors.
    if (nrow(a) == 0) return(invisible(NULL))
    fill <- vapply(a$opacity, function(o) grDevices::adjustcolor(colour, alpha.f = o), character(1))
    plotPairsArches(
      data = a, params = region,
      archHeight = "nSupport", range = c(0, max_count), flip = flip,
      fill = fill, linecolor = fill, alpha = 1, baseline = TRUE,
      x = x, y = top, width = width, height = half_h
    )
  }
  half(arches$treat, treat_colour, flip = FALSE, top = y)
  half(arches$ctrl,  ctrl_colour,  flip = TRUE,  top = mid_y)
  # Either half's arcs draw the shared baseline; with none in view at all, draw it alone.
  if (nrow(arches$treat) == 0 && nrow(arches$ctrl) == 0) .pt_arc_baseline(x, mid_y, width)

  text <- function(label, tx, ty, just, size = fontsize - 2, face = "plain")
    plotText(label, x = tx, y = ty, just = just, fontsize = size, fontface = face,
             fontfamily = fontfamily, fontcolor = fontcolor)

  # Count axis: 0 on the baseline, the shared maximum at the top and bottom.
  if (show_axis) {
    axis_x <- x - 0.08
    ticks  <- c(y, mid_y, bottom)
    plotSegments(axis_x, y, axis_x, bottom, default.units = "inches", lwd = 1, linecolor = fontcolor)
    plotSegments(axis_x - 0.04, ticks, axis_x, ticks, default.units = "inches", lwd = 1, linecolor = fontcolor)
    max_label <- sprintf("%d", as.integer(round(max_count)))
    text(c(max_label, "0", max_label), axis_x - 0.06, ticks, c("right", "center"))
  }

  # Opacity legend: a black gradient from full opacity (|log2FC| >= cap) down to 0, as vector
  # rects (.pt_gradient_bar(), boilerplate.R).
  legend_h <- min(height * 0.6, 1.0)
  legend_x <- x + width + 0.2
  legend_y <- mid_y - legend_h / 2
  text("Opacity", legend_x + 0.09, legend_y - 0.08, c("center", "bottom"), size = fontsize - 1, face = "bold")
  .pt_gradient_bar(grDevices::rgb(0, 0, 0, alpha = seq(1, 0, length.out = 100)),
                   x = legend_x, y = legend_y, width = 0.18, height = legend_h)
  text(sprintf("%.2f log2FC", cap), legend_x + 0.22, legend_y, c("left", "top"))
  text("(full)", legend_x + 0.22, legend_y + 0.11, c("left", "top"))
  text("0", legend_x + 0.22, legend_y + legend_h, c("left", "bottom"))
  invisible(NULL)
}
