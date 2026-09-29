# Download script assembly: no template, no token-filling. A downloaded
# script is exactly the subset of this package's own R/ source that render_all()
# actually needs, pasted verbatim (recovered from PT_SOURCE_TEXT, captured at
# build time -- see data-raw/build-source-text.R) plus a small generated data
# block (deparsed entries/region_params/settings) and one render_all() call at
# the end -- R/ is the single source of truth for what a viz DOES; this only
# decides what to paste and in what order.
#
# What to paste is discovered, not hand-listed: starting from render_all() and
# every settings_<key>()/render_<key>()/load_<key>() function, codetools::findGlobals()
# walks each function's global references transitively (a function calling
# another package function pulls that function in too, and so on) until no new
# names are found. This is what lets adding a new viz-<key>.R file "just work"
# for a downloaded script with nothing to update here -- the old
# PIPELINE_FILE_ORDER/BOILERPLATE_UTIL_DEPS lists this replaces were a second,
# independently-maintained copy of "which files matter" that only stayed
# correct by luck.
#
# Real R function *definitions* don't depend on execution order among
# themselves (unlike render_all()'s old SETTINGS_FNS/RENDER_FNS, which had to
# be built last, as a side effect, from already-defined functions) -- so the
# discovered set is pasted in one fixed (alphabetical, for reproducibility)
# order with no dependency-ordering logic needed.

# Packages loaded by the assembled script -- kept as an explicit fixed list
# (not discovered) because it lives outside package function bodies with no
# global-reference footprint for findGlobals() to walk. This must stay in sync
# with what render_all() and the viz/loader functions actually call via `::`
# or bare names; DESCRIPTION's Imports: is the same list plus UI-only packages
# (shiny, sortable, jsonlite) this standalone script never needs.
.pt_download_script_libs <- c(
  "plotgardener", "rtracklayer", "GenomicRanges", "GenomeInfoDb",
  "AnnotationDbi", "strawr", "grid", "ggplot2", "svglite"
)

# One object's source text, function or constant. Functions: exact original
# text (comments/formatting preserved), read from PT_SOURCE_TEXT
# (R/sysdata.rda, built by data-raw/build-source-text.R) rather than the
# live object's own `srcref` attribute -- an installed package doesn't carry
# srcrefs by default (needs R_KEEP_PKG_SOURCE=yes/--with-keep.source, which
# no ordinary install sets), so reading srcref at runtime broke this for
# every normal install. Non-function top-level values (e.g. SCALE_BAR_TYPES):
# never had a srcref either way, so `deparse()` is the fallback -- loses
# comments but is exactly value-equivalent, which is all correctness
# requires here.
.pt_source_text <- function(name, ns) {
  obj <- get(name, envir = ns)
  # Backticks make the LHS valid for both ordinary names and non-syntactic
  # ones (infix operators like `%||%`, which a bare `%||% <- ...` fails to
  # parse) -- harmless for ordinary names, where backticked and bare are
  # identical.
  lhs <- paste0("`", name, "`")
  if (is.function(obj)) {
    src <- PT_SOURCE_TEXT[[name]]
    if (is.null(src)) {
      stop("PlotTrellis: '", name, "' has no captured source text -- was ",
           "data-raw/build-source-text.R run after this function was added to R/?")
    }
    paste0(lhs, " <- ", src)
  } else {
    paste0(lhs, " <- ", paste(deparse(obj), collapse = "\n"))
  }
}

# Transitive closure of every global name render_all() and the viz/loader
# functions reference, restricted to names that actually exist in this
# package's namespace (base/other-package functions are left as `::`-free bare
# calls, resolved at runtime by the assembled script's own library() calls).
.pt_download_deps <- function() {
  ns <- asNamespace("PlotTrellis")
  # all.names = TRUE -- otherwise every dot-prefixed internal helper
  # (.pt_dispatch_fn, .pd_*, .gene_missing, ...) is invisible to this
  # discovery and silently missing from the assembled script.
  ns_names <- ls(ns, all.names = TRUE)

  viz_keys <- .pt_viz_keys()
  seeds <- c(
    "render_all", "%||%",
    paste0("settings_", viz_keys), paste0("render_", viz_keys),
    ls(ns, pattern = "^load_")
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

# Main assembly function -- returns the downloadable script as a single string.
assemble_download_script <- function(entries, region_params, settings) {
  ns <- asNamespace("PlotTrellis")
  deps <- .pt_download_deps()
  body_src <- paste(vapply(deps, .pt_source_text, character(1), ns = ns), collapse = "\n\n")

  header <- paste(paste0("library(", .pt_download_script_libs, ")"), collapse = "\n")

  dl_settings <- settings
  dl_settings$output_path <- NULL  # tempfile path from the live app, not meaningful standalone

  data_block <- paste0(
    "entries <- ", paste(deparse(entries), collapse = "\n"), "\n\n",
    "region_params <- ", paste(deparse(region_params), collapse = "\n"), "\n\n",
    "settings <- ", paste(deparse(dl_settings), collapse = "\n")
  )

  footer <- paste(
    "library(settings$orgdb_pkg, character.only = TRUE)",
    'OUTPUT <- paste0("plottrellis_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".svg")',
    "render_all(entries, region_params, settings,",
    "           device_fn = function(w, h) svglite::svglite(OUTPUT, width = w, height = h))",
    "dev.off()",
    'cat("Written to", OUTPUT, "\\n")',
    sep = "\n"
  )

  paste(header, body_src, data_block, footer, sep = "\n\n")
}
