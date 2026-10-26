settings_sashimi <- function(dict, entry, settings, track_w, region_params) {
  pooled <- tryCatch(load_junctions(entry$path, region_params), error = function(e) NULL)
  entry$data <- list(raw = if (is.null(pooled)) data.frame() else pooled)

  if (is.null(entry$height)) entry$height <- 0.85

  # Project-wide setting, threaded through the entry since ctx has no field for it.
  entry$clip <- settings$arc_clip_partial %||% FALSE

  min_support <- settings$min_junction_support %||% 0
  raw <- entry$data$raw
  filtered <- raw[raw$nSupport >= min_support, ]

  key <- paste0("sashimi_range_", gsub("[^A-Za-z0-9._]", "_", entry$name))
  this_max <- if (nrow(filtered) > 0) max(filtered$nSupport) else 0
  cur <- dict[[key]]
  dict[[key]] <- if (is.null(cur)) this_max else max(cur, this_max)
  list(dict = dict, entry = entry)
}

render_sashimi <- function(entry, dict, ctx) {
  # nrow(pooled) == 0 needs its own branch: a length-1 scalar chrom column
  # against otherwise-empty vector columns errors on plain recycling.
  sashimi_to_arches <- function(chrom, pooled) {
    if (nrow(pooled) == 0) {
      return(data.frame(chrom1 = character(0), start1 = integer(0), end1 = integer(0),
                         chrom2 = character(0), start2 = integer(0), end2 = integer(0),
                         nSupport = integer(0)))
    }
    data.frame(
      chrom1 = chrom, start1 = pooled$intron_start, end1 = pooled$intron_start + 1L,
      chrom2 = chrom, start2 = pooled$intron_end,   end2 = pooled$intron_end + 1L,
      nSupport = pooled$nSupport
    )
  }

  pooled <- entry$data$raw
  filtered <- pooled[pooled$nSupport >= ctx$min_junction_support, ]
  arches <- sashimi_to_arches(ctx$region_params$chrom, filtered)
  if (nrow(arches) == 0) {
    no_data_text(ctx$track_x, ctx$y, ctx$track_w, ctx$height, ctx$fontsize - 1, ctx$font, ctx$fontcol)
    return(invisible(NULL))
  }

  key   <- paste0("sashimi_range_", gsub("[^A-Za-z0-9._]", "_", entry$name))  # must match settings_sashimi's key
  raw   <- dict[[key]]
  range <- if (!is.null(raw) && is.finite(raw) && raw > 0) c(0, raw) else c(0, 1)

  plotPairsArches(
    data = arches, params = ctx$region,
    archHeight = "nSupport", range = range, flip = FALSE, clip = entry$clip,
    fill = entry$colour, linecolor = entry$colour, alpha = 0.55, baseline = TRUE,
    x = ctx$track_x, y = ctx$y, width = ctx$track_w, height = ctx$height
  )

  label <- sprintf(" \u00B7 %s (max_h=%d reads)", entry$sample, max(arches$nSupport))
  mixed_text(
    list(list(text = entry$name, colour = ctx$fontcol, fontface = "bold"), list(text = label, colour = entry$colour)),
    x = ctx$track_x + 0.05, y = ctx$y + 0.08, fontsize = ctx$fontsize - 1, fontfamily = ctx$font
  )

  if (isTRUE(ctx$show_track_axes)) {
    draw_value_axis(range, ctx$track_x - 0.08, ctx$y, ctx$height, ctx$fontsize - 2, ctx$font, ctx$fontcol, integer_ticks = TRUE)
  }
  invisible(NULL)
}
