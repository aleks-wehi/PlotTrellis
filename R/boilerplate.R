`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a

# Fixed height (inches) of the bottom genome axis, not derived from the font size.
GENOME_AXIS_RESERVE <- 0.35

nice_bp <- function(r) {
  m    <- 10^floor(log10(r * 0.15))
  cand <- m * c(1, 2, 5, 10)
  cand[which.min(abs(cand / r - 0.15))]
}
.fmt_num <- function(x, digits = 2) {
  s <- formatC(x, format = "f", digits = digits)
  s <- sub("0+$", "", s); sub("\\.$", "", s)
}
format_bp_generic <- function(n, digits = 2, space_before_bp = FALSE) {
  if (n >= 1e6) paste0(.fmt_num(n / 1e6, digits), "Mb")
  else if (n >= 1e3) paste0(.fmt_num(n / 1e3, digits), "Kb")
  else paste0(n, if (space_before_bp) " bp" else "bp")
}

# TRUE when `gene` names no gene: NULL, not exactly one value, NA, or a blank string.
.pt_no_gene <- function(gene) length(gene) != 1 || is.na(gene) || !nzchar(trimws(gene))

# One entry per TxDb file for the session: the open TxDb and, once asked for, its CDS index.
# A pure cache keyed by the file's path, so results never depend on call order.
.pt_txdb_memo <- new.env(parent = emptyenv())

# The annotation a render needs, from settings: the TxDb (opened once per file, reopened
# when the file changes, e.g. after Clear caches), the lookup paths, and cds_by_tx(), which
# builds the whole-transcriptome CDS index on first call and returns the same one after.
.pt_annotation <- function(settings) {
  path <- settings$txdb_path
  m <- .pt_txdb_memo[[path]]
  if (is.null(m) || !identical(m$mtime, file.mtime(path))) {
    m <- new.env(parent = emptyenv())
    m$mtime <- file.mtime(path)
    # Named as its package would name it: loadDb() called without a package name attaches
    # GenomicFeatures and what it depends on to the user's search path, noisily.
    loadNamespace("GenomicFeatures")
    m$txdb <- AnnotationDbi::loadDb(path, packageName = "GenomicFeatures")
    assign(path, m, envir = .pt_txdb_memo)
  }
  list(txdb = m$txdb, gtf_path = settings$gtf_path, orgdb_pkg = settings$orgdb_pkg,
       tx_name_lookup_path = settings$tx_name_lookup_path, tx_meta_lookup_path = settings$tx_meta_lookup_path,
       cds_by_tx = function() {
         if (is.null(m$cds)) m$cds <- GenomicFeatures::cdsBy(m$txdb, by = "tx", use.names = TRUE)
         m$cds
       })
}

mixed_text <- function(segments, x, y, fontsize, fontfamily) {
  seg_w <- vapply(segments, function(seg) {
    fface <- if (is.null(seg$fontface)) "plain" else seg$fontface
    grid::convertWidth(grid::grobWidth(grid::textGrob(seg$text, gp = grid::gpar(fontsize = fontsize, fontfamily = fontfamily, fontface = fface))), "inches", valueOnly = TRUE)
  }, numeric(1))
  text_h <- grid::convertHeight(grid::grobHeight(grid::textGrob("Xg", gp = grid::gpar(fontsize = fontsize, fontfamily = fontfamily))), "inches", valueOnly = TRUE)
  pad <- 0.03
  plotRect(x = x - pad, y = y - pad, width = sum(seg_w) + 2 * pad, height = text_h + 2 * pad,
           just = c("left", "top"), default.units = "inches", fill = "white", linecolor = NA)
  cur_x <- x
  for (si in seq_along(segments)) {
    seg   <- segments[[si]]
    fface <- if (is.null(seg$fontface)) "plain" else seg$fontface
    plotText(label = seg$text, x = cur_x, y = y, just = c("left", "top"),
             fontsize = fontsize, fontfamily = fontfamily, fontface = fface,
             fontcolor = seg$colour)
    cur_x <- cur_x + seg_w[si]
  }
}

# A per-sample track's heading at the slot's top left: `title` in bold, then `detail`
# (by default " \u00B7 <sample>") in the sample's identity colour.
.pt_heading <- function(ctx, entry, title = entry$name, detail = paste0(" \u00B7 ", entry$sample)) {
  mixed_text(
    list(list(text = title, colour = ctx$fontcol, fontface = "bold"),
         list(text = detail, colour = entry$colour)),
    x = ctx$track_x + 0.05, y = ctx$y + 0.08, fontsize = ctx$fontsize - 1, fontfamily = ctx$font
  )
}

# A comparison track's heading at the slot's top left:
# "<title> \u00B7 <Comparison><cmp_note><Reference> (Ref)", the title in bold and each
# sample in its role colour.
.pt_heading_pair <- function(ctx, entry, title, cmp_note = " (Cmp) vs ") {
  mixed_text(
    list(
      list(text = title, colour = ctx$fontcol, fontface = "bold"),
      list(text = " \u00B7 ", colour = ctx$fontcol),
      list(text = entry$treatment_label, colour = entry$treatment_colour),
      list(text = cmp_note, colour = ctx$fontcol),
      list(text = entry$control_label, colour = entry$control_colour),
      list(text = " (Ref)", colour = ctx$fontcol)
    ),
    x = ctx$track_x + 0.05, y = ctx$y + 0.08, fontsize = ctx$fontsize - 1, fontfamily = ctx$font
  )
}

# The `dict` key under which the tracks of one type with the same data name share a value
# range: `prefix`, then the data's name made safe for a list name. A module's settings and
# render functions both build the key here, so the two always agree.
.pt_shared_range_key <- function(prefix, entry) paste0(prefix, gsub("[^A-Za-z0-9._]", "_", entry$name))

draw_value_axis <- function(range, x, y_top, h, fontsize, fontfamily, fontcolor, integer_ticks = FALSE) {
  if (diff(range) <= 0) return(invisible(NULL))
  ticks <- pretty(range)
  if (integer_ticks) ticks <- unique(round(ticks))
  ticks <- ticks[ticks >= range[1] & ticks <= range[2]]
  if (length(ticks) == 0) return(invisible(NULL))
  ys <- y_top + h * (1 - (ticks - range[1]) / diff(range))
  plotSegments(x, y_top, x, y_top + h, default.units = "inches", lwd = 1, linecolor = fontcolor)
  for (i in seq_along(ticks)) {
    plotSegments(x - 0.04, ys[i], x, ys[i], default.units = "inches", lwd = 1, linecolor = fontcolor)
    plotText(label = ticks[i], x = x - 0.06, y = ys[i], just = c("right", "center"), fontsize = fontsize, fontfamily = fontfamily, fontcolor = fontcolor)
  }
}

no_data_text <- function(x, y, w, h, fontsize, fontfamily, fontcolor, label) {
  plotText(label = label, x = x + w / 2, y = y + h / 2, just = "center",
           fontsize = fontsize, fontfamily = fontfamily, fontcolor = fontcolor)
}

# A status message ("Chromosome not found", "No gene selected for this region") across a
# track's whole slot, in the track's text style: what a module draws when the track does not
# apply to the view. An empty view is not one: the track draws as usual with nothing in it.
.pt_status_text <- function(ctx, label) {
  no_data_text(ctx$track_x, ctx$y, ctx$track_w, ctx$height, ctx$fontsize - 1, ctx$font, ctx$fontcol, label = label)
}

# A data file that is missing or can't be read fails its track (render_all() draws the red
# "track failed: <reason>"); it is never drawn as an empty track. Stops unless `path` is a file.
.pt_check_file <- function(path) {
  if (length(path) != 1 || is.na(path) || !nzchar(path)) stop("no data file given", call. = FALSE)
  if (dir.exists(path)) stop(basename(path), " is a folder, not a data file (in ", dirname(path), ")", call. = FALSE)
  if (!file.exists(path)) stop("file not found: ", basename(path), " (in ", dirname(path), ")", call. = FALSE)
}

# Stops with a reader's error `e` for `path`, the file named by its base name so the reason
# survives the cut of a long track error message.
.pt_read_error <- function(path, e) {
  stop("could not read ", basename(path), ": ", gsub(path, basename(path), conditionMessage(e), fixed = TRUE),
       call. = FALSE)
}

# The grey line plotPairsArches() draws along the base of its arcs (its baseline.color and
# baseline.lwd defaults), from `x` to `x + width` at `y` (inches): an arc track with no arcs in
# view draws only this, as plotPairsArches() draws nothing at all for zero arcs.
.pt_arc_baseline <- function(x, y, width) {
  plotSegments(x, y, x + width, y, default.units = "inches", linecolor = "grey", lwd = 1)
}

# A vertical colour bar as stacked vector rects (first colour at the top), so a
# legend stays editable in an SVG instead of being a raster image. Colours with alpha
# are flattened onto white, as they look on the page; each rect overlaps the next a
# little so no seams show.
.pt_gradient_bar <- function(colours, x, y, width, height) {
  rgba <- grDevices::col2rgb(as.vector(colours), alpha = TRUE) / 255
  a <- rgba[4, ]
  flat <- grDevices::rgb(t(rgba[1:3, , drop = FALSE]) * a + (1 - a))
  n <- length(flat)
  step <- height / n
  plotRect(x = x, y = y + (seq_len(n) - 1) * step, width = width, height = c(rep(step * 1.5, n - 1), step),
           just = c("left", "top"), default.units = "inches", fill = flat, linecolor = NA)
  invisible(NULL)
}

# Every entry is a GROUP of equivalent spellings for one chromosome; a spelling
# with a capture group (e.g. "chr(\\d+)" paired with "(\\d+)") means the two
# spellings agree on the same digit, so resolve_chrom() only needs one group
# entry to cover every digit/X/Y chromosome instead of enumerating chr1..chr22.
CHROM_ALIASES <- list(
  c("chr([0-9XY]+)", "([0-9XY]+)"),
  c("chrM", "MT", "M", "chrMT")
)

# Resolves region_params$chrom against `vocab` (a file's own chromosome
# names). Returns vocab's own spelling, or NULL (soft "no data", never a
# crash) if nothing matches. vocab = NULL trusts the caller's chrom as-is.
resolve_chrom <- function(chrom, vocab) {
  if (is.null(vocab)) return(chrom)
  match_in_vocab <- function(cand) {
    hit <- vocab[tolower(vocab) == tolower(cand)]
    if (length(hit) > 0) hit[1] else NULL
  }
  direct <- match_in_vocab(chrom)
  if (!is.null(direct)) return(direct)

  for (group in CHROM_ALIASES) {
    for (spelling in group) {
      m <- regmatches(chrom, regexec(paste0("^", spelling, "$"), chrom, ignore.case = TRUE))[[1]]
      if (length(m) == 0) next
      captured <- if (length(m) > 1) m[2] else NA_character_
      for (other in setdiff(group, spelling)) {
        has_capture <- grepl("(", other, fixed = TRUE)
        cand <- if (!has_capture) other else if (!is.na(captured)) sub("\\([^)]*\\)", captured, other) else NA_character_
        if (!is.na(cand)) {
          hit <- match_in_vocab(cand)
          if (!is.null(hit)) return(hit)
        }
      }
    }
  }
  NULL
}

draw_scale_bar <- function(ctx) {
  sc <- ctx$scale
  by <- ctx$y + 0.12
  plotSegments(sc$bx0, by,        sc$bx1, by,        default.units = "inches", lwd = 1.5, linecolor = ctx$fontcol)
  plotSegments(sc$bx0, by - 0.05, sc$bx0, by + 0.05, default.units = "inches", lwd = 1.5, linecolor = ctx$fontcol)
  plotSegments(sc$bx1, by - 0.05, sc$bx1, by + 0.05, default.units = "inches", lwd = 1.5, linecolor = ctx$fontcol)
  plotText(sc$label, x = sc$x_mid, y = by - 0.02, just = c("center", "bottom"),
           fontsize = ctx$fontsize, fontfamily = ctx$font, fontcolor = ctx$fontcol)
}
