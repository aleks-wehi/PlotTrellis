#' The `settings_<key>()`/`render_<key>()` module contract
#'
#' A viz module is a pair of top-level functions, `settings_<key>()` and
#' `render_<key>()`, defined anywhere under `R/` (by convention in a
#' `viz-<key>.R` file). Adding a new pair is enough for both the live render
#' path (`render_all()`) and the downloaded-script assembly path
#' (`assemble_download_script()`, `download-script.R`) to pick it up
#' automatically -- see [viz-module-discovery] for how that lookup works.
#' There is no list of modules to update by hand.
#'
#' `settings_<key>(dict, entry, settings, track_w, region_params)` runs once
#' per track before any drawing happens. It loads/validates this track's data,
#' fills in `entry$height` if not already set, and can read or write shared
#' state via `dict` (e.g. a cross-track value range every same-type track
#' needs to agree on). It returns `list(dict = dict, entry = entry)`.
#'
#' `render_<key>(entry, dict, ctx)` does the actual plotgardener drawing for
#' one track, using the `entry`/`dict` from its matching `settings_<key>()`
#' call plus `ctx` (page geometry, region, fonts, per-render-call options)
#' built once by `render_all()` and passed to every track's render function.
#'
#' A module with only one half defined is silently excluded from
#' [.pt_viz_keys()] rather than erroring, so an in-progress module doesn't
#' break discovery for every other one.
#'
#' @keywords internal
#' @name viz-module-contract
NULL

#' Viz-module discovery
#'
#' Internal helpers backing the [viz-module-contract]: resolving a
#' `settings_<key>()`/`render_<key>()` pair by name, and enumerating every
#' key that has both halves defined.
#'
#' @keywords internal
#' @name viz-module-discovery
NULL

#' @rdname viz-module-discovery
.pt_dispatch_fn <- function(prefix, key) {
  # environment(sys.function()) rather than a hardcoded asNamespace("PlotTrellis")
  # -- resolves to the package namespace when running inside the package, and
  # to globalenv() when this function is pasted standalone into a downloaded
  # script (assemble_download_script(), download-script.R), where no
  # PlotTrellis namespace exists at all.
  ns <- environment(sys.function())
  fn <- get0(paste0(prefix, "_", key), envir = ns, mode = "function", inherits = FALSE)
  if (is.null(fn)) stop("PlotTrellis: no ", prefix, "_", key, "() is defined for viz type '", key, "'")
  fn
}

#' @rdname viz-module-discovery
# All viz keys that have BOTH a settings_<key>() and a render_<key>() defined
# -- the two halves of a complete viz module. A key with only one half is
# excluded rather than erroring here, since a module mid-development (only
# one function written so far) shouldn't break discovery for every other one.
.pt_viz_keys <- function() {
  ns <- environment(sys.function())
  settings_keys <- sub("^settings_", "", ls(ns, pattern = "^settings_"))
  render_keys   <- sub("^render_",   "", ls(ns, pattern = "^render_"))
  intersect(settings_keys, render_keys)
}

render_all <- function(items, region_params, settings, device_fn = NULL) {
  page_w  <- settings$page_w %||% 9
  track_x <- 0.8
  track_w <- page_w - track_x - 1.0
  first_y <- 0.5
  gap     <- 0.10
  font     <- settings$font %||% "DejaVu Sans"
  fontsize <- settings$fontsize %||% 8L
  fontcol  <- "black"

  txdb <- AnnotationDbi::loadDb(settings$txdb_path)
  asm  <- assembly(Genome = settings$genome_label, TxDb = txdb, OrgDb = settings$orgdb_pkg,
                   gene.id.column = "ENSEMBL", display.column = "SYMBOL")
  region <- pgParams(chrom = region_params$chrom, chromstart = region_params$start,
                     chromend = region_params$end, assembly = asm)

  # Built from only the caller-supplied items -- genome_axis (appended below) isn't a data track.
  title_items <- items
  items <- c(items, list(list(viz = "genome_axis", height = NULL)))

  dict <- list()
  for (i in seq_along(items)) {
    settings_fn <- .pt_dispatch_fn("settings", items[[i]]$viz)
    res <- settings_fn(dict, items[[i]], settings, track_w, region_params)
    dict <- res$dict
    items[[i]] <- res$entry
  }
  for (it in items) if (is.null(it$height)) stop("render_all(): settings_", it$viz, "() must set an entry height")

  heights <- vapply(items, `[[`, numeric(1), "height")
  y       <- first_y + c(0, cumsum(head(heights + gap, -1)))
  page_h  <- tail(y, 1) + tail(heights, 1)

  region_width <- region_params$end - region_params$start
  scale_len <- nice_bp(region_width)
  scale_w   <- scale_len * track_w / region_width
  scale_bx1 <- track_x + track_w - 0.05
  scale_bx0 <- scale_bx1 - scale_w
  scale <- list(bx0 = round(scale_bx0, 3), bx1 = round(scale_bx1, 3),
                x_mid = round((scale_bx0 + scale_bx1) / 2, 3), label = format_bp_generic(scale_len))

  title <- if (nchar(trimws(settings$title %||% "")) > 0) {
    settings$title
  } else {
    types    <- vapply(title_items, `[[`, character(1), "viz")
    types    <- types[types != "genes"]  # annotation overlay, not a data type
    gene     <- region_params$gene
    gene_lbl <- if (is.null(gene) || is.na(gene) || nchar(trimws(gene)) == 0) "" else paste0(gene, " \u00B7 ")
    paste0(
      gene_lbl, paste(toupper(unique(types)), collapse = " + "), " \u00B7 ",
      region_params$chrom, ":", format_bp_generic(region_params$start, 3), "\u2013", format_bp_generic(region_params$end, 3)
    )
  }

  if (!is.null(device_fn)) device_fn(page_w, page_h)  # caller opens its own SVG device sized to the actual content, e.g. svglite()

  pageCreate(width = page_w, height = page_h, default.units = "inches", showGuides = FALSE)
  plotText(title, x = track_x + track_w / 2, y = 0.2, just = "center",
           fontfamily = font, fontsize = fontsize + 2, fontface = "bold", fontcolor = fontcol)

  for (i in seq_along(items)) {
    entry <- items[[i]]
    ctx <- list(
      index = i, is_first = (i == 1L) && (entry$viz %in% SCALE_BAR_TYPES),
      y = y[[i]], height = heights[[i]],
      track_x = track_x, track_w = track_w, region = region,
      font = font, fontsize = fontsize, fontcol = fontcol,
      hic_norm = settings$hic_norm %||% "KR", show_track_axes = settings$show_track_axes %||% TRUE,
      scale = scale, region_params = region_params, chr_fmt = settings$chr_fmt %||% "dot",
      diff_colour_hi = settings$diff_colour_hi %||% "#E41A1C", diff_colour_lo = settings$diff_colour_lo %||% "#377EB8",
      min_junction_support = settings$min_junction_support %||% 0,
      junction_opacity_cap = settings$junction_opacity_cap %||% 2,
      txdb = txdb, gtf_path = settings$gtf_path, orgdb_pkg = settings$orgdb_pkg,
      genome_label = settings$genome_label
    )
    render_fn <- .pt_dispatch_fn("render", entry$viz)
    render_fn(entry, dict, ctx)
  }

  invisible(NULL)
}
