info_junctions <- function() list(
  kind = "sample", data_type = "paired-end", order = 40,
  filter = c(junctions = "Junctions"),
  path_field = "junction_path"
)

# Mirrored comparison junction tracks keep their heading in a band of this
# height (inches) above the arcs, which take twice junction_height (one for each
# sample). A single-sample junction track has no band: like a coverage track,
# its arcs fill the whole junction_height slot and the heading is drawn over
# their top left.
JUNCTION_HEADING_H <- 0.3

# Junctions (one row each: intron_start, intron_end, nSupport) as
# plotPairsArches() input: an arc per junction between its two ends, with
# nSupport (read count) as its height.
junction_arches <- function(junctions, chrom) {
  chrom <- rep(chrom, nrow(junctions))
  data.frame(
    chrom1 = chrom,
    start1 = junctions$intron_start,
    end1 = junctions$intron_start + 1L,
    chrom2 = chrom,
    start2 = junctions$intron_end,
    end2 = junctions$intron_end + 1L,
    nSupport = junctions$nSupport
  )
}

# TRUE when a junction file's scores (column 5 of a BED12 or junction BED6) are
# all "." on the data lines .pt_arc_file_kind() looks at: read_junction_file()
# then counts every junction once, so the heights say nothing and the heading
# shows no max. Such a file is still "junctions" (it keeps its TPM-scaled
# comparison). A STAR SJ.out.tab's column 5 is the intron motif, never ".".
.pt_junction_scoreless <- function(path) {
  tryCatch({
    first <- suppressWarnings(.pt_arc_first_data(path))
    data <- suppressWarnings(readLines(path, n = first$skip + 50, warn = FALSE))
    data <- data[seq_along(data) > first$skip]
    data <- data[nzchar(trimws(data)) & !startsWith(data, "#")]
    length(data) > 0 &&
      all(trimws(vapply(strsplit(data, "\t"), `[`, "", 5)) %in% ".")
  }, error = function(e) FALSE)
}

# The junctions whose arcs plotPairsArches(clip = TRUE) draws: both ends in the
# window (junction_arches() makes each end a 1 bp anchor). With arc_clip_partial
# on, the junction tracks keep only these, so their scale and heading count only
# arcs that are drawn.
.pt_clip_partial <- function(junctions, region_params) {
  junctions[junctions$intron_start >= region_params$start &
              junctions$intron_end + 1 <= region_params$end, , drop = FALSE]
}

settings_junctions <- function(dict, entry, settings, track_w, region_params) {
  entry$data <- list(raw = load_junctions(entry$path, region_params))

  if (is.null(entry$height)) entry$height <- settings$junction_height

  # Project-wide setting, threaded through the entry since ctx has no field for
  # it.
  entry$clip <- settings$arc_clip_partial

  min_support <- settings$min_junction_support
  raw <- entry$data$raw
  filtered <- raw[raw$nSupport >= min_support, ]
  if (isTRUE(entry$clip)) filtered <- .pt_clip_partial(filtered, region_params)

  key <- .pt_shared_range_key("junctions_range_", entry)
  this_max <- if (nrow(filtered) > 0) max(filtered$nSupport) else 0
  cur <- dict[[key]]
  dict[[key]] <- if (is.null(cur)) this_max else max(cur, this_max)
  list(dict = dict, entry = entry)
}

render_junctions <- function(entry, dict, ctx) {
  pooled <- entry$data$raw
  filtered <- pooled[pooled$nSupport >= ctx$min_junction_support, ]
  if (isTRUE(entry$clip))
    filtered <- .pt_clip_partial(filtered, ctx$region_params)

  key   <- .pt_shared_range_key("junctions_range_", entry)
  raw   <- dict[[key]]
  range <- if (!is.null(raw) && is.finite(raw) && raw > 0) c(0, raw) else
    c(0, 1)

  arcs_y <- ctx$y
  arcs_h <- ctx$height
  # No junctions in view (or none with min_junction_support reads): the track as
  # usual, no arcs.
  if (nrow(filtered) > 0) {
    plotPairsArches(
      data = junction_arches(filtered, ctx$region_params$chrom),
      params = ctx$region,
      archHeight = "nSupport", range = range, flip = FALSE, clip = entry$clip,
      fill = entry$colour, linecolor = entry$colour, alpha = 0.55,
      baseline = TRUE,
      x = ctx$track_x, y = arcs_y, width = ctx$track_w, height = arcs_h
    )
  } else {
    .pt_arc_baseline(ctx$track_x, arcs_y + arcs_h, ctx$track_w)
  }

  # %s, not %d: a BEDPE source's score column can be non-integer (e.g.
  # FitHiChIP). The tallest arc's read count (splice junctions) or score
  # (loops); a file without scores has no height to show.
  max_h <- if (nrow(filtered) > 0) max(filtered$nSupport) else 0
  kind  <- .pt_arc_file_kind(entry$path)
  unscored <- identical(kind, "unscored") ||
    (identical(kind, "junctions") && .pt_junction_scoreless(entry$path))
  label <- if (unscored) sprintf(" \u00B7 %s", entry$sample) else
    sprintf(" \u00B7 %s (max=%s)", entry$sample,
            format(max_h, digits = 3, scientific = FALSE))
  .pt_heading(ctx, entry, detail = label)

  if (isTRUE(ctx$show_track_axes)) {
    draw_value_axis(range, ctx$track_x - 0.08, arcs_y, arcs_h, ctx$fontsize - 2,
                    ctx$font, ctx$fontcol, integer_ticks = TRUE)
  }
  invisible(NULL)
}
