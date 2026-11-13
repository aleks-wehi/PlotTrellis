info_junctions <- function() list(
  kind = "sample", data_type = "paired-end", order = 40, filter = c(junctions = "Junctions"),
  path_field = "junction_path"
)

# Junction tracks keep their heading in a band of this height (inches) above the arcs. Below
# it a single-sample track's arcs take junction_height, a mirrored comparison's twice that
# (one junction_height for each sample).
JUNCTION_HEADING_H <- 0.3

# Junctions (one row each: intron_start, intron_end, nSupport) as plotPairsArches() input:
# an arc per junction between its two ends, with nSupport (read count) as its height.
junction_arches <- function(junctions, chrom) {
  chrom <- rep(chrom, nrow(junctions))
  data.frame(
    chrom1 = chrom, start1 = junctions$intron_start, end1 = junctions$intron_start + 1L,
    chrom2 = chrom, start2 = junctions$intron_end,   end2 = junctions$intron_end + 1L,
    nSupport = junctions$nSupport
  )
}

settings_junctions <- function(dict, entry, settings, track_w, region_params) {
  entry$data <- list(raw = load_junctions(entry$path, region_params))

  if (is.null(entry$height)) entry$height <- JUNCTION_HEADING_H + settings$junction_height

  # Project-wide setting, threaded through the entry since ctx has no field for it.
  entry$clip <- settings$arc_clip_partial

  min_support <- settings$min_junction_support
  raw <- entry$data$raw
  filtered <- raw[raw$nSupport >= min_support, ]

  key <- .pt_shared_range_key("junctions_range_", entry)
  this_max <- if (nrow(filtered) > 0) max(filtered$nSupport) else 0
  cur <- dict[[key]]
  dict[[key]] <- if (is.null(cur)) this_max else max(cur, this_max)
  list(dict = dict, entry = entry)
}

render_junctions <- function(entry, dict, ctx) {
  pooled <- entry$data$raw
  filtered <- pooled[pooled$nSupport >= ctx$min_junction_support, ]

  key   <- .pt_shared_range_key("junctions_range_", entry)
  raw   <- dict[[key]]
  range <- if (!is.null(raw) && is.finite(raw) && raw > 0) c(0, raw) else c(0, 1)

  arcs_y <- ctx$y + JUNCTION_HEADING_H
  arcs_h <- ctx$height - JUNCTION_HEADING_H
  # No junctions in view (or none with min_junction_support reads): the track as usual, no arcs.
  if (nrow(filtered) > 0) {
    plotPairsArches(
      data = junction_arches(filtered, ctx$region_params$chrom), params = ctx$region,
      archHeight = "nSupport", range = range, flip = FALSE, clip = entry$clip,
      fill = entry$colour, linecolor = entry$colour, alpha = 0.55, baseline = TRUE,
      x = ctx$track_x, y = arcs_y, width = ctx$track_w, height = arcs_h
    )
  } else {
    .pt_arc_baseline(ctx$track_x, arcs_y + arcs_h, ctx$track_w)
  }

  # %s, not %d: a BEDPE source's score column can be non-integer (e.g. FitHiChIP).
  max_h <- if (nrow(filtered) > 0) max(filtered$nSupport) else 0
  label <- sprintf(" \u00B7 %s (max_h=%s reads)", entry$sample, format(max_h, digits = 3, scientific = FALSE))
  .pt_heading(ctx, entry, detail = label)

  if (isTRUE(ctx$show_track_axes)) {
    draw_value_axis(range, ctx$track_x - 0.08, arcs_y, arcs_h, ctx$fontsize - 2, ctx$font, ctx$fontcol, integer_ticks = TRUE)
  }
  invisible(NULL)
}
