#' The `info_<key>()`/`settings_<key>()`/`render_<key>()` module contract
#'
#' A track type ("viz module") is a set of top-level functions named after its
#' key, by convention all in one `R/viz-<key>.R` file: the pair
#' `settings_<key>()`/`render_<key>()`, which draw it, and `info_<key>()`, which
#' declares what it is. They are found by name ([viz-module-discovery]), both by
#' `render_all()` and by the downloaded-script assembly (`download-script.R`).
#' Every list of track types the app and the code API use is derived from the
#' declarations: adding a track type is one new file plus a rebuild of
#' `R/sysdata.rda` (`Rscript --vanilla data-raw/build-source-text.R` in the
#' package folder), which holds the source text a downloaded script is pasted
#' from; without the rebuild, `pt_script()` stops for any figure using the new
#' track type.
#'
#' `settings_<key>(dict, entry, settings, track_w, region_params)` runs once per
#' track (once per panel for a row-panel type, below), in figure order, before
#' any drawing happens. It loads this track's data into `entry$data`, sets
#' `entry$height` (inches) if it is not already set, and can read or write
#' `dict`, a named list every track's settings call passes on to the next (e.g.
#' a value range all tracks of the type share). It returns
#' `list(dict = dict, entry = entry)`. `settings` has every display setting the
#' caller left unset filled from its default (`.pt_fill_settings()`).
#'
#' `render_<key>(entry, dict, ctx)` draws the track with plotgardener, from its
#' `entry`, the final `dict` and `ctx`, built by `render_all()` for each track:
#' its slot (`y`, `height`, `track_x`, `track_w`, inches from the page's top
#' left), `region` (plotgardener params with the assembly), `region_params`
#' (`chrom`, `start`, `end`, `gene`, the gene being `NA` or absent for a
#' coordinates view), `font`, `fontsize`, `fontcol`, `scale` and `is_first`
#' (below), `index`, the project-wide display options (`show_track_axes`,
#' `hic_norm`, `chr_fmt`, `diff_colour_hi`/`diff_colour_lo`,
#' `min_junction_support`, `junction_opacity_cap`) and the genome's annotation;
#' the `ctx` list in `render_all()` is the complete set.
#'
#' A view with nothing in it (the file is read fine but has no data in the region)
#' is not a special case: the track draws as usual (heading, axes, shared range),
#' with nothing plotted and no message. Besides drawing, a track ends in one of
#' three ways:
#' \describe{
#'   \item{The track does not apply to the view}{(no gene, `.pt_no_gene()`; a
#'     chromosome the file lacks): the module draws a message across its slot with
#'     `.pt_status_text(ctx, label)` and returns.}
#'   \item{The track can't be drawn}{(a data file missing or unreadable,
#'     `.pt_check_file()`, `.pt_read_error()`; a wrong format): the module calls
#'     `stop(..., call. = FALSE)`; an error in either function becomes an in-figure
#'     "<key> track failed: ..." message in red in that track's slot, and the rest
#'     of the figure still draws.}
#'   \item{The figure can't be drawn}{An error outside the modules (the genome's
#'     annotation, an unknown track type, no height set) stops `render_all()`.}
#' }
#'
#' The `entry` a module receives, by kind: a per-sample track has `viz` (the
#' key), `sample`, `colour` (the sample's identity colour), `name` (its data's
#' name in `samples.txt`, e.g. `atac`) and `path` (the data file), or the fields
#' its own `detect_<key>()` put on the track item (below); a two-sample
#' comparison has `viz`, `name`, `control_label`/`treatment_label`,
#' `control_path`/`treatment_path`, `control_colour`/`treatment_colour` (the
#' project's two role colours) and any fields its `enrich_<key>()` added; a
#' sample-free track has only `viz`. Each also has `height` (a preset, else
#' `NULL`) and `data` (`NULL` until its settings call).
#'
#' `info_<key>()` returns a plain list; only `kind` (and `data_type` for a track
#' with data) is required:
#' \describe{
#'   \item{`kind`}{`"sample"` (one track per sample that has the data),
#'     `"comparison"` (one track per ordered pair of samples that both have data
#'     of the same name) or `"fixed"` (one sample-free track).}
#'   \item{`data_type`}{The `samples.txt` data type it reads: one of
#'     `DATA_TYPE_VOCAB` (`bw`, `bed`, `hic`, `quantdir`, `paired-end`).}
#'   \item{`label`}{Its label in the app; `{name}` stands for the data's name.
#'     The sample(s) follow it. Default `"{name}"`, the key for a fixed track.}
#'   \item{`api`}{Its track name in [pt_tracks()]/[pt_figure()], `{name}` as in
#'     `label`. Default `"{name}"`, the key for a fixed track. When two types give
#'     the same name, the type detected later gets `_<key>` appended.}
#'   \item{`order`}{Its position in detection among tracks of its kind (a
#'     number, default last). It orders [pt_tracks()] and decides which of two
#'     same-named types keeps the plain name.}
#'   \item{`filter`}{The app's filter button that lists it, as
#'     `c(<button id> = "<button label>")`; types sharing an id share the button.
#'     A button whose tracks are all comparisons shows in Two-sample mode only.
#'     Default none (listed under All only).}
#'   \item{`scale_bar`}{`TRUE`: when it is the figure's first track, `ctx$is_first`
#'     is `TRUE` (unless the scale bar is switched off) and the module draws the
#'     scale bar with `draw_scale_bar(ctx)`. Default `FALSE`.}
#'   \item{`row_panel`}{`TRUE`: adjacent tracks of the type draw as one panel
#'     (below). Default `FALSE`.}
#'   \item{`needs_gene`}{`TRUE`: it draws only for a gene, not for coordinates;
#'     the app flags it in a coordinates view. Default `FALSE`.}
#'   \item{`path_field`}{The track item field holding a per-sample track's data
#'     path. It reaches the entry as `path`; an item made by a module's own
#'     `detect_<key>()` with no path there (`NA`, empty) gives an entry without `path`.
#'     Default `"file"`.}
#'   \item{`entry_settings`}{Names of project settings copied into each of its
#'     entries, so a downloaded script carries them per track. Default none.}
#' }
#' A per-sample track is offered for each of a sample's data of `data_type`; a
#' module with its own rule defines `detect_<key>(sample_row, info)`, returning
#' a list of track items. A comparison can add fields to its items with
#' `enrich_<key>(items, control_row, treatment_row)`, returning the items (an
#' empty list offers none). A declaration that is malformed (an unknown kind or
#' field, a missing data type, a filter button given two labels) is an error
#' naming the module, raised where tracks are listed (the app's start,
#' [pt_project()]). A module without
#' `info_<key>()` is offered nowhere, but still draws in a figure that names it
#' (the genome axis `render_all()` adds is one).
#'
#' Adjacent tracks of a `row_panel` type draw as one panel with a row per track:
#' such a run (a lone track too) reaches its module as a single entry,
#' `list(viz = <key>, rows = <the run's entries, in figure order>)`, plus
#' `figure_rows`, the rows of every panel of that type in the figure. A row panel
#' sizes itself: a height preset on one of its tracks is not used.
#'
#' A key with only one of `settings_<key>()`/`render_<key>()` is silently
#' excluded from [.pt_viz_keys()], so an in-progress module doesn't break
#' discovery for every other one; an `info_<key>()` without both of them is an
#' error.
#'
#' @keywords internal
#' @name viz-module-contract
NULL

#' Viz-module discovery
#'
#' Internal helpers backing the [viz-module-contract]: resolving a module's
#' functions by name, enumerating every key that has both halves defined, and
#' reading the modules' `info_<key>()` declarations.
#'
#' @keywords internal
#' @name viz-module-discovery
NULL

#' @rdname viz-module-discovery
# A module's `<prefix>_<key>()` (settings, render, info, detect, enrich), or NULL.
.pt_module_fn <- function(prefix, key) {
  # environment(sys.function()) rather than a hardcoded asNamespace("PlotTrellis")
  # -- resolves to the package namespace when running inside the package, and
  # to globalenv() when this function is pasted standalone into a downloaded
  # script (assemble_download_script(), download-script.R), where no
  # PlotTrellis namespace exists at all.
  ns <- environment(sys.function())
  get0(paste0(prefix, "_", key), envir = ns, mode = "function", inherits = FALSE)
}

#' @rdname viz-module-discovery
.pt_dispatch_fn <- function(prefix, key) {
  fn <- .pt_module_fn(prefix, key)
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

#' @rdname viz-module-discovery
# One TRUE/FALSE field of a module's info_<key>() declaration (`scale_bar`, `row_panel`),
# FALSE when the module has no declaration or leaves the field out. This is all drawing
# needs of a declaration, so a downloaded script carries the declarations but no checking of them.
.pt_viz_flag <- function(key, flag) {
  fn <- .pt_module_fn("info", key)
  !is.null(fn) && isTRUE(fn()[[flag]])
}

#' @rdname viz-module-discovery
# A module's info_<key>() declaration, checked, with every field filled in, plus
# `key` and `offered` (FALSE for a module without one: it is offered nowhere,
# but still draws when a figure names it).
.pt_viz_info <- function(key) {
  fn <- .pt_module_fn("info", key)
  info <- if (is.null(fn)) list() else fn()
  bad <- function(...) stop("PlotTrellis: info_", key, "() ", ..., call. = FALSE)
  if (!is.list(info)) bad("must return a list")
  fields <- c("kind", "data_type", "label", "api", "order", "filter", "scale_bar", "row_panel",
              "needs_gene", "path_field", "entry_settings")
  unknown <- setdiff(names(info), fields)
  if (length(unknown)) bad("has an unknown field '", unknown[1], "' (fields: ", paste(fields, collapse = ", "), ")")
  is_string <- function(x) is.character(x) && length(x) == 1 && !is.na(x) && nzchar(x)
  is_flag   <- function(x) is.logical(x) && length(x) == 1 && !is.na(x)
  if (!is.null(fn)) {
    if (!(is_string(info$kind) && info$kind %in% c("sample", "comparison", "fixed")))
      bad("field 'kind' must be \"sample\", \"comparison\" or \"fixed\"")
    if (info$kind != "fixed" && !(is_string(info$data_type) && info$data_type %in% DATA_TYPE_VOCAB))
      bad("field 'data_type' must be one of ", paste(DATA_TYPE_VOCAB, collapse = ", "))
  }
  fixed <- identical(info$kind, "fixed")
  info <- modifyList(list(kind = NA_character_, data_type = NA_character_,
                          label = if (fixed) key else "{name}", api = if (fixed) key else "{name}",
                          order = Inf, scale_bar = FALSE, row_panel = FALSE, needs_gene = FALSE,
                          path_field = "file", entry_settings = character(0)), info)
  for (f in c("label", "api", "path_field")) if (!is_string(info[[f]])) bad("field '", f, "' must be one string")
  for (f in c("scale_bar", "row_panel", "needs_gene")) if (!is_flag(info[[f]])) bad("field '", f, "' must be TRUE or FALSE")
  if (!(is.numeric(info$order) && length(info$order) == 1 && !is.na(info$order))) bad("field 'order' must be one number")
  if (!is.null(info$filter) && !(is_string(info$filter) && is_string(names(info$filter))))
    bad("field 'filter' must be c(<button id> = \"<button label>\")")
  if (!is.character(info$entry_settings)) bad("field 'entry_settings' must be a character vector")
  c(list(key = key, offered = !is.null(fn)), info)
}

#' @rdname viz-module-discovery
# Every offered module's declaration (.pt_viz_info()), named by key, in detection
# order: by `order`, then key.
.pt_viz_infos <- function() {
  keys <- .pt_viz_keys()
  # Projects, figures and start views kept from before the junction track was renamed
  # still say "sashimi"; their keys are read as "junctions" (.pt_renamed_junctions_key()).
  reserved <- grep("^sashimi(_|$)", keys, value = TRUE)
  if (length(reserved))
    stop("PlotTrellis: a track type's key can't be 'sashimi' or start with 'sashimi_' (", reserved[1],
         "): that prefix is read as the old name of the junctions track", call. = FALSE)
  ns <- environment(sys.function())
  declared <- sub("^info_", "", ls(ns, pattern = "^info_"))
  declared <- Filter(function(key) !is.null(.pt_module_fn("info", key)), declared)
  orphans <- setdiff(declared, keys)
  if (length(orphans))
    stop("PlotTrellis: info_", orphans[1], "() is defined but settings_", orphans[1], "() and render_",
         orphans[1], "() are not both", call. = FALSE)
  infos <- setNames(lapply(keys, .pt_viz_info), keys)
  infos <- Filter(function(info) info$offered, infos)
  .pt_check_filter_labels(infos)
  infos[order(vapply(infos, `[[`, numeric(1), "order"), names(infos), method = "radix")]
}

# Types that share a filter button (the same id in `filter`) must give it the same label.
.pt_check_filter_labels <- function(infos) {
  infos <- Filter(function(info) !is.null(info$filter), infos)
  ids <- vapply(infos, function(info) names(info$filter), character(1))
  labels <- vapply(infos, function(info) unname(info$filter), character(1))
  for (id in unique(ids)) {
    if (length(unique(labels[ids == id])) > 1)
      stop("PlotTrellis: filter button '", id, "' has two labels: ",
           paste(unique(labels[ids == id]), collapse = ", "), call. = FALSE)
  }
}

# Per-track error isolation: an error in one track's settings_/render_ call
# becomes an in-figure message naming the track and the reason, so the rest
# of the figure still renders. Also echoed to the console via message().
.pt_track_error_label <- function(entry, e) {
  # A row panel (.pt_group_row_panels() below) is named by its rows' samples, the first three of many.
  samples <- unlist(lapply(entry$rows, `[[`, "sample"))
  if (length(samples) > 3) samples <- c(samples[1:3], sprintf("+%d more", length(samples) - 3))
  sample <- entry$sample %||% paste(samples, collapse = ", ")
  who <- if (length(sample) == 1 && nzchar(sample)) sprintf("%s (%s)", entry$viz, sample) else entry$viz
  # iconv first: a message with invalid bytes would otherwise make gsub() fail inside this handler.
  msg <- trimws(gsub("\\s+", " ", iconv(conditionMessage(e), to = "UTF-8", sub = "?")))
  # A long one is cut in the middle: its start says what went wrong, its end is often the useful
  # end of a path (the sample's folder in "file not found: quant.sf (in <folder>)").
  if (nchar(msg) > 140) msg <- paste0(substr(msg, 1, 90), "...", substr(msg, nchar(msg) - 46, nchar(msg)))
  label <- sprintf("%s track failed: %s", who, msg)
  message("PlotTrellis: ", label)
  label
}
.pt_draw_track_error <- function(label, ctx) {
  # Blank the slot first: a render_ call can fail after drawing part of its track.
  plotRect(x = ctx$track_x, y = ctx$y, width = ctx$track_w, height = ctx$height,
           just = c("left", "top"), default.units = "inches", fill = "white", linecolor = NA)
  no_data_text(ctx$track_x, ctx$y, ctx$track_w, ctx$height, ctx$fontsize - 1, ctx$font, "#CC0000", label = label)
  # The figure's scale bar belongs to its first track: it stays when that track fails.
  if (isTRUE(ctx$is_first)) draw_scale_bar(ctx)
}

# Folds each run of adjacent tracks of a row-panel type (`row_panel` in its info_<key>(),
# see the viz-module-contract above) into one entry holding the run's entries as `rows`.
# Each such entry also gets `figure_rows`, the rows of every panel of its type in the
# figure, so separate panels can agree on a layout.
.pt_group_row_panels <- function(items) {
  out <- list()
  for (it in items) {
    n <- length(out)
    if (!.pt_viz_flag(it$viz, "row_panel")) out[[n + 1]] <- it
    else if (n > 0 && identical(out[[n]]$viz, it$viz)) out[[n]]$rows <- c(out[[n]]$rows, list(it))
    else out[[n + 1]] <- list(viz = it$viz, rows = list(it))
  }
  lapply(out, function(it) {
    if (!is.null(it$rows)) it$figure_rows <- Filter(function(row) identical(row$viz, it$viz), items)
    it
  })
}

# render_all()'s settings with each display setting left unset (NULL or NA) taken from its
# default (.pt_display_defaults()), so modules read settings$<name> as it is.
.pt_fill_settings <- function(settings) {
  defaults <- Filter(Negate(is.null), .pt_display_defaults())
  for (k in names(defaults)) settings[[k]] <- settings[[k]] %||% defaults[[k]]
  settings
}

# The family the figure's text is drawn in for the Font setting `font`. On macOS a family the
# system lacks (the default DejaVu Sans is not shipped there) is drawn as the generic "sans":
# Quartz (png()) measures and draws a missing family with no font at all (zero-width,
# undrawn labels), and svglite matches its bold face to an unrelated font. Elsewhere the
# family is used as it is: fontconfig already substitutes a close face on every device.
# Installed means the system's own fonts (Quartz can't use systemfonts' registered ones); the
# generic families are those quartzFonts() defines.
.pt_device_font <- function(font, sysname = Sys.info()[["sysname"]],
                            installed = function() systemfonts::system_fonts()$family) {
  if (!identical(sysname, "Darwin") || font %in% c("", "sans", "serif", "mono")) return(font)
  if (font %in% installed()) font else "sans"
}

render_all <- function(items, region_params, settings, device_fn = NULL) {
  settings <- .pt_fill_settings(settings)
  page_w  <- settings$page_w
  track_x <- 0.8
  track_w <- page_w - track_x - 1.0
  # The title and the scale bar can be switched off (both on by default). There
  # is no auto-title either: no title text means no title. Without one the tracks start near
  # the top of the page instead of below it.
  title      <- trimws(settings$title %||% "")
  show_title <- isTRUE(settings$show_title) && nzchar(title)
  first_y <- if (show_title) 0.5 else 0.1
  gap     <- 0.10
  font     <- .pt_device_font(settings$font)
  fontsize <- settings$fontsize
  fontcol  <- "black"

  anno <- .pt_annotation(settings)
  txdb <- anno$txdb
  asm  <- assembly(Genome = settings$genome_label, TxDb = txdb, OrgDb = settings$orgdb_pkg,
                   gene.id.column = "ENSEMBL", display.column = "SYMBOL")
  region <- pgParams(chrom = region_params$chrom, chromstart = region_params$start,
                     chromend = region_params$end, assembly = asm)

  items <- c(.pt_group_row_panels(items), list(list(viz = "genome_axis", height = NULL)))

  dict <- list()
  for (i in seq_along(items)) {
    settings_fn <- .pt_dispatch_fn("settings", items[[i]]$viz)
    res <- tryCatch(settings_fn(dict, items[[i]], settings, track_w, region_params), error = function(e) {
      entry <- items[[i]]
      entry$track_error <- .pt_track_error_label(entry, e)
      entry$height <- entry$height %||% 0.4
      list(dict = dict, entry = entry)
    })
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

  if (!is.null(device_fn)) device_fn(page_w, page_h)  # caller opens its own SVG device sized to the actual content, e.g. svglite()

  pageCreate(width = page_w, height = page_h, default.units = "inches", showGuides = FALSE)
  if (show_title) plotText(title, x = track_x + track_w / 2, y = 0.2, just = "center",
                           fontfamily = font, fontsize = fontsize + 2, fontface = "bold", fontcolor = fontcol)

  for (i in seq_along(items)) {
    entry <- items[[i]]
    ctx <- list(
      index = i, is_first = (i == 1L) && .pt_viz_flag(entry$viz, "scale_bar") && isTRUE(settings$show_scale_bar),
      y = y[[i]], height = heights[[i]],
      track_x = track_x, track_w = track_w, region = region,
      font = font, fontsize = fontsize, fontcol = fontcol,
      hic_norm = settings$hic_norm, show_track_axes = settings$show_track_axes,
      scale = scale, region_params = region_params, chr_fmt = settings$chr_fmt,
      diff_colour_hi = settings$diff_colour_hi, diff_colour_lo = settings$diff_colour_lo,
      min_junction_support = settings$min_junction_support,
      junction_opacity_cap = settings$junction_opacity_cap,
      txdb = txdb, tx_meta_lookup_path = settings$tx_meta_lookup_path, orgdb_pkg = settings$orgdb_pkg,
      genome_label = settings$genome_label, cds_by_tx = anno$cds_by_tx, cache_dir = settings$cache_dir
    )
    if (!is.null(entry$track_error)) {
      .pt_draw_track_error(entry$track_error, ctx)
      next
    }
    render_fn <- .pt_dispatch_fn("render", entry$viz)
    tryCatch(render_fn(entry, dict, ctx),
             error = function(e) .pt_draw_track_error(.pt_track_error_label(entry, e), ctx))
  }

  invisible(NULL)
}
