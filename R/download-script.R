# Download script assembly: no template, no token-filling. A downloaded
# script is exactly the subset of this package's own R/ source that one figure's
# render_all() call reaches, pasted verbatim (recovered from PT_SOURCE_TEXT, captured at
# build time -- see data-raw/build-source-text.R) below a short header, the library()
# calls and an editable data block (entries/region_params/settings), with one
# render_all() call at the end -- R/ is the single source of truth for what a viz DOES;
# this only decides what to paste and in what order.
#
# What to paste is discovered, not hand-listed: starting from render_all() and the
# info_<key>()/settings_<key>()/render_<key>() functions of the figure's own track types,
# codetools::findGlobals() walks each function's global references transitively (a function calling
# another package function pulls that function in too, and so on) until no new
# names are found. This is what lets adding a new viz-<key>.R file "just work"
# for a downloaded script with nothing to update here.
#
# Real R function *definitions* don't depend on execution order among
# themselves, so the discovered set is pasted in an order chosen for the reader: render_all(), each
# track type's info_/settings_/render_ functions, then the other functions and the constants.

# One object's source text, function or constant. Functions: exact original
# text (comments/formatting preserved), read from PT_SOURCE_TEXT
# (R/sysdata.rda, built by data-raw/build-source-text.R) rather than the
# live object's own `srcref` attribute -- an installed package doesn't carry
# srcrefs by default (needs R_KEEP_PKG_SOURCE=yes/--with-keep.source, which
# no ordinary install sets), so reading srcref at runtime would fail on
# every normal install. Non-function top-level values (e.g. GENOME_AXIS_RESERVE):
# never had a srcref either way, so `deparse()` is the fallback -- loses
# comments but is exactly value-equivalent, which is all correctness
# requires here.
.pt_source_text <- function(name, ns) {
  obj <- get(name, envir = ns)
  # A name that isn't syntactic (the infix operator `%||%`) needs backticks to be assigned to.
  lhs <- if (identical(make.names(name), name)) name else paste0("`", name, "`")
  if (is.function(obj)) {
    src <- PT_SOURCE_TEXT[[name]]
    if (is.null(src)) {
      stop("PlotTrellis: '", name, "' has no captured source text -- was ",
           "data-raw/build-source-text.R run after this function was added to R/?", call. = FALSE)
    }
    paste0(lhs, " <- ", src)
  } else if (is.environment(obj)) {
    # An environment (a session cache such as .pt_txdb_memo) deparses to
    # `<environment>`, which does not parse: the script starts it empty.
    paste0(lhs, " <- new.env(parent = emptyenv())")
  } else {
    paste0(lhs, " <- ", paste(deparse(obj), collapse = "\n"))
  }
}

# The track types a figure's render_all() call dispatches to: its entries' own, in
# figure order, plus the genome axis render_all() always appends.
.pt_script_keys <- function(entries) unique(c(vapply(entries, `[[`, character(1), "viz"), "genome_axis"))

# Transitive closure of every global name render_all() and the figure's own
# info_<key>()/settings_<key>()/render_<key>() functions reference, restricted to names that
# actually exist in this package's namespace (base/other-package functions are left as `::`-free
# bare calls, resolved at runtime by the assembled script's own library() calls).
# render_all() reaches a module only through .pt_dispatch_fn() and .pt_viz_flag(), by a
# computed name, so its functions are seeds; a row panel keeps its tracks' key, so the
# same functions draw it.
.pt_download_deps <- function(entries) {
  ns <- asNamespace("PlotTrellis")
  # all.names = TRUE -- otherwise every dot-prefixed internal helper
  # (.pt_dispatch_fn, .pd_*, .pt_no_gene, ...) is invisible to this
  # discovery and silently missing from the assembled script.
  ns_names <- ls(ns, all.names = TRUE)

  viz_keys <- .pt_script_keys(entries)
  seeds <- c(
    "render_all", "%||%", ".pt_tidy_svg", ".pt_quiet",
    paste0("info_", viz_keys), paste0("settings_", viz_keys), paste0("render_", viz_keys)
  )
  seeds <- intersect(seeds, ns_names)

  found <- character(0)
  frontier <- seeds
  while (length(frontier) > 0) {
    frontier <- setdiff(frontier, found)
    if (length(frontier) == 0) break
    found <- union(found, frontier)
    next_frontier <- character(0)
    for (nm in frontier) {
      obj <- get(nm, envir = ns)
      if (!is.function(obj)) next
      globals <- tryCatch(codetools::findGlobals(obj, merge = TRUE), error = function(e) character(0))
      next_frontier <- union(next_frontier, intersect(globals, ns_names))
    }
    frontier <- next_frontier
  }
  sort(found)
}

# Packages the pasted functions need, read off their own references. `attach`: the
# packages their bare names are imported from in this package's NAMESPACE, in import
# order (so a name two of them export resolves as it does in the package), without R's
# default packages, which Rscript attaches anyway. `need`: those plus every package they
# call as `pkg::`, without R's own base packages; the script checks each is installed.
.pt_script_packages <- function(deps, ns) {
  imports <- getNamespaceImports(ns)
  imports <- imports[nzchar(names(imports)) & vapply(imports, is.character, logical(1))]
  # Imported name -> its package; a later import of the same name wins, as in R.
  from <- character(0)
  for (p in names(imports)) from[names(imports[[p]])] <- p
  attach <- colons <- character(0)
  for (nm in deps) {
    obj <- get(nm, envir = ns)
    if (!is.function(obj)) next
    pkg <- from[setdiff(codetools::findGlobals(obj, merge = TRUE), deps)]
    attach <- c(attach, pkg[!is.na(pkg)])
    pd <- utils::getParseData(parse(text = PT_SOURCE_TEXT[[nm]], keep.source = TRUE))
    colons <- c(colons, pd$text[pd$token == "SYMBOL_PACKAGE"])
  }
  attach <- setdiff(intersect(names(imports), attach),
                    c("base", "methods", "datasets", "utils", "grDevices", "graphics", "stats"))
  need <- unique(c(attach, colons))
  is_base <- vapply(need, function(p) identical(utils::packageDescription(p, fields = "Priority"), "base"), logical(1))
  list(attach = attach, need = need[!is_base])
}

# R code for a value, laid out for reading: a list that fits in `width` from column `col`
# stays on one line, a longer one packs its fields onto lines under its first field, and a
# field that itself needs several lines starts a line of its own. Anything but a plain
# list (a vector, or a list with a class such as a data.frame) is written by deparse().
.pt_format_value <- function(x, col = 0, width = 100) {
  if (!is.list(x) || !all(names(attributes(x)) %in% "names")) {
    txt <- paste(deparse(x, width.cutoff = 500L), collapse = " ")
    # deparse() rounds doubles to 15 significant digits; 17 always read back exactly.
    if (!identical(eval(str2lang(txt), baseenv()), x)) {
      txt <- paste(deparse(x, width.cutoff = 500L, control = c("keepNA", "keepInteger", "niceNames",
                                                               "showAttributes", "digits17")), collapse = " ")
    }
    return(txt)
  }
  lhs <- .pt_field_names(x)
  inner <- col + nchar("list(")
  items <- vapply(seq_along(x), function(i) paste0(lhs[i], .pt_format_value(x[[i]], inner + nchar(lhs[i]), width)), "")
  lines <- character(0)
  for (it in items) {
    n <- length(lines)
    if (n > 0 && !grepl("\n", it) && !grepl("\n", lines[n]) && inner + nchar(lines[n]) + 2 + nchar(it) <= width) {
      lines[n] <- paste0(lines[n], ", ", it)
    } else {
      lines <- c(lines, it)
    }
  }
  paste0("list(", paste(lines, collapse = paste0(",\n", strrep(" ", inner))), ")")
}

# "name = " for each element of a list (backticked when not syntactic), "" when unnamed.
.pt_field_names <- function(x) {
  nms <- names(x) %||% rep("", length(x))
  ifelse(nzchar(nms), paste0(ifelse(make.names(nms) == nms, nms, paste0("`", nms, "`")), " = "), "")
}

# A list with one element per line: the figure's tracks, or its settings.
.pt_format_block <- function(x) {
  if (length(x) == 0) return("list()")
  lhs <- .pt_field_names(x)
  items <- vapply(seq_along(x), function(i) paste0(lhs[i], .pt_format_value(x[[i]], 2 + nchar(lhs[i]))), "")
  paste0("list(\n", paste0("  ", items, collapse = ",\n"), "\n)")
}

# A track entry without its NULL fields (an absent field reads as NULL too), except one
# whose name begins another field's name, which `$` would then partially match.
.pt_drop_null_fields <- function(entry) {
  drop <- vapply(names(entry), function(f) {
    is.null(entry[[f]]) && !any(startsWith(setdiff(names(entry), f), f))
  }, logical(1))
  entry[!drop]
}

# The header's one-line descriptions of the region and the tracks.
.pt_script_region_line <- function(region_params, settings) {
  gene <- region_params$gene
  where <- sprintf("%s:%s-%s", region_params$chrom, format(region_params$start, big.mark = ","),
                   format(region_params$end, big.mark = ","))
  if (is.character(gene) && length(gene) == 1 && !is.na(gene) && nzchar(gene)) where <- paste0(gene, ", ", where)
  if (!is.null(settings$genome_label)) where <- paste0(where, " (", settings$genome_label, ")")
  where
}
# Each track by its label in the app, e.g. "atac (Basal)", "Transcript expression (K562)".
.pt_script_tracks_line <- function(entries) {
  paste(vapply(entries, function(e) {
    what <- .pt_fill_name(.pt_viz_info(e$viz)$label, e$name %||% e$viz)
    who <- if (!is.null(e$treatment_label)) paste(e$treatment_label, "vs", e$control_label) else e$sample
    if (length(who) == 1 && !is.na(who) && nzchar(who)) paste0(what, " (", who, ")") else what
  }, character(1)), collapse = ", ")
}

# Main assembly function -- returns the downloadable script as a single string.
assemble_download_script <- function(entries, region_params, settings) {
  ns <- asNamespace("PlotTrellis")
  deps <- .pt_download_deps(entries)
  pkgs <- .pt_script_packages(deps, ns)
  # svglite writes the figure (the last lines); the OrgDb is looked up by name, never attached.
  pkgs$need <- union(pkgs$need, c("svglite", settings$orgdb_pkg,
                                  if (any(vapply(entries, function(e) identical(e$viz, "bed") && isTRUE(grepl("\\.gz$", e$path, ignore.case = TRUE)),
                                                 logical(1)))) "R.utils"))

  dl_entries <- lapply(entries, .pt_drop_null_fields)
  dl_settings <- settings
  dl_settings$output_path <- NULL  # tempfile path from the live app, not meaningful standalone
  # A setting a track type copies into each of its entries (entry_settings) is drawn from
  # the entry: the block shows it there only, where editing it has an effect.
  dl_settings[unique(unlist(lapply(.pt_viz_infos(), `[[`, "entry_settings")))] <- NULL
  # The order of the app's settings, then the render-only fields (paths, title).
  dl_settings <- dl_settings[c(intersect(names(.pt_default_settings()), names(dl_settings)),
                               setdiff(names(dl_settings), names(.pt_default_settings())))]
  gene <- region_params$gene
  output <- paste0("plottrellis_",
                   if (is.character(gene) && length(gene) == 1 && !is.na(gene) && nzchar(gene)) gene
                   else paste(region_params$chrom, region_params$start, region_params$end, sep = "_"),
                   ".svg")

  header <- c(
    "# PlotTrellis figure script: redraws one figure with plotgardener; PlotTrellis itself",
    "# need not be installed. Run it with `Rscript <this file>`.",
    # One line each, whatever a gene or sample is called.
    paste0("# Region: ", gsub("[\r\n]+", " ", .pt_script_region_line(region_params, settings))),
    paste0("# Tracks: ", gsub("[\r\n]+", " ", .pt_script_tracks_line(entries))),
    sprintf("# Written by PlotTrellis %s on %s.", utils::packageVersion("PlotTrellis"), format(Sys.Date())),
    paste0("# Needs the R packages: ", paste(pkgs$need, collapse = ", ")),
    "# It reads the data and annotation files named in `entries` and `settings`.",
    "#",
    "# To change the figure, edit the block after the library() calls: OUTPUT (the file",
    "# written), entries (the tracks, top to bottom), region_params (the region drawn) and",
    "# settings. Everything after that block is PlotTrellis's own code, unchanged, and the",
    "# last lines draw the figure."
  )
  libs <- c(
    paste0("for (pkg in ", .pt_format_value(pkgs$need), ")"),
    "  if (!requireNamespace(pkg, quietly = TRUE))",
    "    stop(\"This script needs the R package '\", pkg, \"'. Install it with BiocManager::install(\\\"\", pkg, \"\\\").\",",
    "         call. = FALSE)",
    "suppressPackageStartupMessages({",
    paste0("  library(", pkgs$attach, ")"),
    "})"
  )
  data_block <- c(
    "# ---- Edit here ---------------------------------------------------------------------",
    paste0("OUTPUT <- ", .pt_format_value(output), "  # written as SVG, whatever the name ends in"),
    "",
    "# The tracks, top to bottom. colour = \"#rrggbb\" sets the colour of a sample's track (a",
    "# comparison: control_colour, treatment_colour; a track without a sample has none);",
    "# height = <inches> sets a track's height (otherwise the settings below decide; a type",
    "# whose info_<key>() says row_panel sizes itself).",
    "# Tracks can be removed, reordered or repeated; in each list here, every item but the",
    "# last ends with a comma. This script holds the code of its own track types only",
    paste0("# (", paste(c(setdiff(.pt_script_keys(entries), "genome_axis"), "the genome axis"), collapse = ", "),
           "): a track of another type needs a script made from a figure that has one."),
    paste0("entries <- ", .pt_format_block(dl_entries)),
    "",
    "# The region drawn. `gene` is the gene a gene-level track (transcripts, protein domains)",
    "# draws; change the two together.",
    paste0("region_params <- ", .pt_format_value(region_params, nchar("region_params <- "))),
    "",
    "# The project's and the figure's settings (per-type heights are in inches).",
    paste0("settings <- ", .pt_format_block(dl_settings))
  )
  data_text <- paste(data_block, collapse = "\n")
  # The block is what the script draws from: it must give back exactly these values.
  chk <- new.env(parent = ns)
  eval(parse(text = data_text), chk)
  if (!identical(mget(c("entries", "region_params", "settings"), envir = chk),
                 list(entries = dl_entries, region_params = region_params, settings = dl_settings))) {
    stop("PlotTrellis: the script's data block does not reproduce this figure's values")
  }

  keys <- .pt_script_keys(entries)
  pairs <- intersect(as.vector(rbind(paste0("info_", keys), paste0("settings_", keys), paste0("render_", keys))),
                     deps)
  rest <- setdiff(deps, c("render_all", pairs))
  is_fn <- vapply(rest, function(nm) is.function(get(nm, envir = ns)), logical(1))
  pasted <- function(nms) {
    paste(vapply(nms, function(nm) paste(c(attr(PT_SOURCE_TEXT[[nm]], "comment"), .pt_source_text(nm, ns)),
                                         collapse = "\n"), character(1)), collapse = "\n\n")
  }
  pair_src <- vapply(keys, function(k) {
    fns <- intersect(paste0(c("info_", "settings_", "render_"), k), pairs)
    if (length(fns) == 0) "" else paste0("# ---- Track type ", k, " ----\n\n", pasted(fns))
  }, character(1))
  code <- c(
    "# ---- PlotTrellis code, unchanged: render_all() draws the figure ----",
    pasted("render_all"),
    pair_src[nzchar(pair_src)],
    "# ---- Helpers ----",
    pasted(rest[is_fn]),
    "# ---- Constants ----",
    pasted(rest[!is_fn])
  )

  footer <- c(
    "# ---- Draw the figure ----",
    "suppressPackageStartupMessages(.pt_quiet(render_all(",
    "  entries, region_params, settings,",
    "  device_fn = function(w, h) svglite::svglite(OUTPUT, width = w, height = h))))",
    "invisible(dev.off())",
    ".pt_tidy_svg(OUTPUT)",
    "cat(\"Figure written to\", normalizePath(OUTPUT), \"\\n\")"
  )

  paste(c(paste(header, collapse = "\n"), paste(libs, collapse = "\n"), data_text,
          code, paste(footer, collapse = "\n")), collapse = "\n\n")
}
