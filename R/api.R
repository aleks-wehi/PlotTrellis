# Code API (replaces phase-I.3 step 26's Tier 2): open a project once, then
# describe figures by gene + track names and save them. Builds the same
# entries/settings as the composer's Generate, through the shared helpers in
# R/figure_helpers.R. Design notes: solutions/spikes/bioinf-api/README.md.

# The bare-call project folder: visible, in the home directory, one per user.
.pt_default_dir <- function() file.path(path.expand("~"), "PlotTrellis")

#' Open a PlotTrellis project for scripting
#'
#' Loads a project's genome, samples and tracks once, so [pt_figure()] calls are quick. With
#' a `samples_txt`, it first creates the project in `project_dir` from it (the folder is
#' created if needed); without one, it opens the project already in `project_dir`.
#'
#' A project is re-imported automatically when its `samples.txt` is newer than the project's
#' `project_settings.json` (you edited it), so re-running a script after editing
#' `samples.txt` picks the edit up. Otherwise the existing project is opened as it is,
#' keeping any settings changed in the app. One folder holds one project: a different
#' `samples.txt` pointed at a folder that already has a project is an error.
#'
#' @param project_dir The project's folder: its settings and caches live here. Defaults to
#'   the `PlotTrellis` folder in your home directory, as for [run_plottrellis()].
#' @param samples_txt Path to a `samples.txt` describing the samples and genome, e.g.
#'   [pt_example()]. Leave out to open an existing project.
#' @param overwrite If `TRUE`, re-import `samples_txt` even if the project is newer
#'   (discarding settings changed in the app since), or replace a project made from a
#'   different `samples.txt`.
#' @return A `pt_project`. Print it to see its samples and the track names [pt_figure()]
#'   accepts.
#' @seealso [pt_figure()], [pt_tracks()], [pt_example()]
#' @examples
#' proj <- pt_project(file.path(tempdir(), "pt-example-human"), pt_example("human"))
#' proj
#' @export
pt_project <- function(project_dir = .pt_default_dir(), samples_txt = NULL, overwrite = FALSE) {
  .pt_import(project_dir, samples_txt, overwrite)
  cfg <- .pt_project(project_dir)$cfg
  if (!requireNamespace(cfg$orgdb, quietly = TRUE))
    stop("This project needs the Bioconductor package '", cfg$orgdb, "'. Install it with ",
         "BiocManager::install(\"", cfg$orgdb, "\").", call. = FALSE)
  genome <- ensure_genome_ready(gtf_path = cfg$gtf_path, orgdb = cfg$orgdb, genome_label = cfg$genome_label)
  # Attached as in .pt_composer_app(): the pipeline looks the OrgDb up by name.
  suppressPackageStartupMessages(library(genome$orgdb, character.only = TRUE))
  samples <- load_samples_from_settings(cfg$samples, cfg$data_dir)
  items   <- c(build_track_items(samples), build_comparison_track_items(samples))
  structure(list(
    dir        = normalizePath(project_dir),
    cfg        = cfg,
    genome     = genome,
    txdb       = AnnotationDbi::loadDb(genome$txdb_path),
    genes      = sort(AnnotationDbi::keys(get(genome$orgdb, envir = asNamespace(genome$orgdb)), keytype = "SYMBOL")),
    samples    = samples,
    colours    = .pt_sample_colours(samples),
    items      = items,
    tracks     = .pt_track_table(items),
    render_cfg = .pt_render_cfg(cfg, genome)
  ), class = "pt_project")
}

# Creates or refreshes project_dir from samples_txt (NULL: the project must
# already exist), shared by pt_project() and run_plottrellis(). A samples.txt
# may carry `# start_gene:` / `# start_tracks:` lines (track names as in
# pt_figure(), e.g. atac@Basal): the app then opens with that view preloaded.
.pt_import <- function(project_dir, samples_txt = NULL, overwrite = FALSE) {
  json <- file.path(project_dir, "project_settings.json")
  if (is.null(samples_txt)) {
    if (!file.exists(json))
      stop("No PlotTrellis project in '", project_dir, "' yet. Create one with ",
           "pt_project(project_dir, samples_txt = \"path/to/samples.txt\"), or run_plottrellis() to set ",
           "it up in the app.", call. = FALSE)
    return(invisible(project_dir))
  }
  if (!is_real_file(samples_txt)) stop("No such samples.txt: '", samples_txt, "'.", call. = FALSE)
  src <- normalizePath(samples_txt)
  # One folder holds one project: a different samples.txt (or a project set up
  # in the app) already there is an error, not silently opened or replaced.
  made_from <- if (file.exists(json)) jsonlite::read_json(json)$samples_txt_source %||% "the app's setup wizard"
  if (!is.null(made_from) && !identical(made_from, src) && !isTRUE(overwrite))
    stop("'", project_dir, "' already holds a project, made from ", made_from, ". Use another ",
         "project_dir for a separate project, or overwrite = TRUE to replace it.", call. = FALSE)
  if (is.null(made_from) || isTRUE(overwrite) || file.mtime(samples_txt) > file.mtime(json)) {
    if (!is.null(made_from)) message("Re-importing '", basename(samples_txt), "' into '", project_dir, "'.")
    parse_samples_txt(samples_txt, project_dir = project_dir, overwrite = TRUE)
    settings <- jsonlite::read_json(json)
    settings$samples_txt_source <- src
    writeLines(jsonlite::toJSON(settings, auto_unbox = TRUE, pretty = TRUE, null = "null"), json)
    hdr <- function(key) {
      ln <- grep(paste0("^#\\s*", key, ":"), readLines(samples_txt, warn = FALSE), value = TRUE)
      if (length(ln)) trimws(sub(paste0("^#\\s*", key, ":"), "", ln[1])) else NULL
    }
    gene <- hdr("start_gene")
    tracks <- hdr("start_tracks")
    if (!is.null(gene) && !is.null(tracks)) {
      cfg <- .pt_project(project_dir)$cfg
      samples <- load_samples_from_settings(cfg$samples, cfg$data_dir)
      items <- c(build_track_items(samples), build_comparison_track_items(samples))
      keys <- .pt_resolve_tracks(.pt_track_table(items), trimws(strsplit(tracks, ",")[[1]]), NULL)
      .pt_start_view(project_dir, items, keys, gene, names(.pt_sample_colours(samples)))
    }
  }
  invisible(project_dir)
}

#' @export
print.pt_project <- function(x, ...) {
  cat("PlotTrellis project:", x$dir, "\n")
  cat("Genome: ", x$genome$genome_label, " (", basename(x$cfg$gtf_path), ", ", x$genome$orgdb, ")\n", sep = "")
  cat("Samples:", paste(names(x$colours), collapse = ", "), "\n")
  cat("Tracks (for pt_figure(tracks = ...)):\n")
  tr <- x$tracks
  for (t in unique(tr$track)) {
    s <- tr$sample[tr$track == t]
    cat(sprintf("  %-22s %s\n", t, if (all(is.na(s))) "" else paste0("@ ", paste(s, collapse = ", "))))
  }
  invisible(x)
}

# One row per track item: the name pt_figure() accepts, its sample (NA for
# sample-free tracks, "B vs A" for comparisons), kind and TRACK_ITEMS key. Data
# tracks are named after samples.txt's own `type:name` column names.
.pt_track_table <- function(items) {
  cmp_suffix <- c(bw_diff_stacked = "_diff", junctions_mirrored_diff = "_diff",
                  junctions_mirrored_diff_tpm_scaled = "_diff_tpm")
  rows <- lapply(items, function(it) {
    kind <- if (!is.null(it$samples)) "comparison" else if (is.na(it$sample)) "fixed" else "sample"
    track <- switch(kind,
      fixed      = it$type,
      comparison = paste0(it$name, cmp_suffix[[it$type]]),
      sample     = if (it$type == "tpm_usage") "transcript_counts" else it$name)
    data.frame(track = track, type = it$type, kind = kind, key = it$key, label = it$label,
               sample = if (kind == "comparison") paste(it$samples[2], "vs", it$samples[1]) else it$sample,
               stringsAsFactors = FALSE)
  })
  tab <- do.call(rbind, rows)
  tab <- tab[order(match(tab$kind, c("fixed", "sample", "comparison"))), ]
  # Two data types under one column name (e.g. bw:atac and bed:atac): the
  # later types get a _<type> suffix so every name stays unambiguous.
  first_type <- tapply(tab$type, tab$track, `[`, 1)
  dup <- tab$kind != "fixed" & tab$type != first_type[tab$track]
  tab$track[dup] <- paste0(tab$track[dup], "_", tab$type[dup])
  rownames(tab) <- NULL
  tab[, c("track", "sample", "kind", "label", "key")]
}

#' List a project's tracks
#'
#' @param proj A [pt_project()].
#' @return A data frame with one row per available track: `track` (the name to use in
#'   [pt_figure()]), `sample` (`NA` for tracks that need no sample; `"B vs A"` for
#'   comparisons, B compared against reference A), `kind` (`"fixed"`, `"sample"` or
#'   `"comparison"`) and `label` (the app's label for it).
#' @examples
#' pt_tracks(pt_project(file.path(tempdir(), "pt-example-human"), pt_example("human")))
#' @export
pt_tracks <- function(proj) proj$tracks[, c("track", "sample", "kind", "label")]

# tracks -> TRACK_ITEMS keys, in order. A bare per-sample name expands over
# `par` (track-major, so e.g. adjacent transcript_counts tracks merge); a bare
# comparison needs exactly 2 there (second vs first). "name@sample" always
# picks exactly that one, with or without parallel.
.pt_resolve_tracks <- function(tab, tracks, par) {
  keys <- character(0)
  for (spec in tracks) {
    has_at <- grepl("@", spec, fixed = TRUE)
    name <- trimws(sub("@.*$", "", spec))
    smp  <- if (has_at) gsub("\\s+", " ", trimws(sub("^[^@]*@", "", spec))) else NA_character_
    hit <- tab[tab$track == name, ]
    if (nrow(hit) == 0)
      stop("Unknown track '", name, "'. This project has: ", paste(unique(tab$track), collapse = ", "),
           ". See pt_tracks(proj).", call. = FALSE)
    kind <- hit$kind[1]
    if (kind == "fixed") {
      if (has_at) stop("'", name, "' isn't per-sample; write just \"", name, "\".", call. = FALSE)
      keys <- c(keys, hit$key[1])
      next
    }
    want <- if (has_at) smp else if (kind == "sample") par else if (length(par) == 2) paste(par[2], "vs", par[1])
    if (is.null(want)) {
      ex <- hit$sample[1]
      stop(if (kind == "sample")
             paste0("'", name, "' needs a sample: write \"", name, "@", ex, "\", or use parallel = TRUE ",
                    "for every sample.")
           else
             paste0("'", name, "' compares 2 samples: write \"", name, "@", ex, "\" (B vs reference A), or ",
                    "pass exactly 2 samples to parallel, e.g. parallel = c(\"",
                    sub(".* vs ", "", ex), "\", \"", sub(" vs .*", "", ex), "\") (reference first)."),
           call. = FALSE)
    }
    for (w in want) {
      k <- hit$key[hit$sample == w]
      if (length(k) == 0) {
        if (!has_at && kind == "sample") {
          message("Note: no '", name, "' for ", w, ", skipped.")
          next
        }
        near <- hit$sample[tolower(hit$sample) == tolower(w)]
        stop("No '", name, "' track for '", w, "'.",
             if (length(near)) paste0(" Did you mean '", near[1], "'? (Sample names are case-sensitive.)"),
             " It has: ", paste(hit$sample, collapse = ", "), ".", call. = FALSE)
      }
      keys <- c(keys, k[1])
    }
  }
  keys
}

#' Describe a PlotTrellis figure
#'
#' Picks a gene (or region) and the tracks to show, in order, top to bottom. Nothing is
#' drawn until the figure is printed, saved with [pt_save()], or exported with
#' [pt_script()].
#'
#' Track names are the ones [pt_tracks()] lists: `genes`, `transcripts`,
#' `transcripts_all`, `transcripts_annotated` and `protein_domains` need no sample. The
#' others are named after your `samples.txt`'s columns (`bw:atac` gives `atac`;
#' `quantdir:` gives `transcript_counts`; comparisons get a `_diff` suffix, such as
#' `junctions_diff`) and need a sample, written `"atac@K562"` (or
#' `"junctions_diff@K562 vs GM12878"`, K562 against reference GM12878).
#'
#' With `parallel`, a per-sample name written without `@` stands for that track in
#' every sample (one after another), and a comparison without `@` compares the 2 samples
#' (the first is the reference). A name written with `@` always means just that one.
#'
#' @param proj A [pt_project()].
#' @param gene A gene symbol (any case) or Ensembl gene ID. The view is the gene plus the
#'   project's flank setting on each side.
#' @param tracks Character vector of track names, top to bottom.
#' @param parallel `FALSE` (every per-sample track says its sample with `@`), `TRUE`
#'   (bare names cover every sample, in `samples.txt` order), or a character vector of
#'   sample names to cover, in that order.
#' @param region Instead of `gene`: a region such as `"chr1:198,600,000-198,760,000"`.
#' @param title Figure title.
#' @return A `pt_figure`.
#' @seealso [pt_save()], [pt_script()]
#' @examples
#' proj <- pt_project(file.path(tempdir(), "pt-example-human"), pt_example("human"))
#' fig <- pt_figure(proj, "PSME2", c("genes", "transcripts", "atac", "junctions"), parallel = TRUE)
#' pt_save(fig, file.path(tempdir(), "PSME2.pdf"))
#'
#' # The same figure, one sample per track, mixed freely:
#' pt_figure(proj, "PSME2", c("genes", "atac@GM12878", "atac@K562", "junctions_diff@K562 vs GM12878"))
#' @export
pt_figure <- function(proj, gene = NULL, tracks = "genes", parallel = FALSE, region = NULL, title = "") {
  stopifnot(inherits(proj, "pt_project"))
  all_samples <- names(proj$colours)
  par <- if (isTRUE(parallel)) all_samples else if (is.character(parallel)) parallel else NULL
  bad <- setdiff(par, all_samples)
  if (length(bad))
    stop("Unknown sample(s) in parallel: ", paste(bad, collapse = ", "), ". This project has: ",
         paste(all_samples, collapse = ", "), ".", call. = FALSE)

  if (!is.null(gene) && !(is.character(gene) && length(gene) == 1))
    stop("gene must be one gene name or Ensembl ID. With region =, name the tracks: tracks = c(...).",
         call. = FALSE)
  if (!is.null(gene)) {
    res <- .pt_quiet(.pt_resolve_gene_region(gene, proj$txdb, proj$render_cfg$gene_flank_pct / 100,
                                             proj$genome$orgdb, proj$genes, proj$genome$genome_label))
    if (!is.null(res$error)) stop(res$error, call. = FALSE)
    rgn <- res$region
    rgn$gene <- res$symbol
  } else if (!is.null(region)) {
    m <- regmatches(region, regexec("^\\s*([^:\\s]+):([0-9,]+)-([0-9,]+)\\s*$", region, perl = TRUE))[[1]]
    if (length(m) != 4) stop("region must look like \"chr1:198,600,000-198,760,000\".", call. = FALSE)
    rgn <- list(chrom = m[2], start = gsub(",", "", m[3]), end = gsub(",", "", m[4]), gene = NA_character_)
  } else {
    stop("Pass gene = (a symbol or Ensembl ID) or region = (\"chr:start-end\").", call. = FALSE)
  }

  keys <- .pt_resolve_tracks(proj$tracks, tracks, par)
  structure(list(
    gene          = rgn$gene,
    keys          = keys,
    tracks        = proj$tracks$label[match(keys, proj$tracks$key)],
    entries       = .pt_build_entries(keys, proj$items, proj$colours, proj$render_cfg),
    region_params = list(chrom = rgn$chrom, start = as.integer(rgn$start), end = as.integer(rgn$end), gene = rgn$gene),
    settings      = c(.pt_settings_base(proj$render_cfg), list(title = title))
  ), class = "pt_figure")
}

#' @export
print.pt_figure <- function(x, ...) {
  r <- x$region_params
  cat("PlotTrellis figure: ", if (is.na(x$gene)) "" else paste0(x$gene, ", "),
      r$chrom, ":", format(r$start, big.mark = ","), "-", format(r$end, big.mark = ","), "\n", sep = "")
  cat(paste0("  ", x$tracks, "\n"), sep = "")
  viewer <- getOption("viewer")
  if (interactive() && is.function(viewer)) {
    f <- pt_save(x, tempfile(fileext = ".png"))
    viewer(f)
  } else {
    cat("Save it with pt_save(fig, \"figure.pdf\") (or .svg, .png).\n")
  }
  invisible(x)
}

#' Render a figure to a file
#'
#' @param fig A [pt_figure()].
#' @param file Output path ending in `.pdf`, `.svg` or `.png`. Its folder must exist.
#' @param dpi Resolution for `.png`.
#' @return `file`, invisibly.
#' @seealso [pt_figure()], [pt_script()]
#' @export
pt_save <- function(fig, file, dpi = 300) {
  stopifnot(inherits(fig, "pt_figure"))
  dev <- switch(tolower(tools::file_ext(file)),
    svg = function(w, h) svglite::svglite(file, width = w, height = h),
    pdf = function(w, h) grDevices::cairo_pdf(file, width = w, height = h),
    png = function(w, h) grDevices::png(file, width = w, height = h, units = "in", res = dpi),
    stop("pt_save() writes .pdf, .svg or .png files.", call. = FALSE))
  before <- grDevices::dev.cur()
  on.exit(if (grDevices::dev.cur() != before) grDevices::dev.off())
  .pt_quiet(render_all(fig$entries, fig$region_params, fig$settings, device_fn = dev))
  invisible(file)
}

# Hides plotgardener's per-grob "segments[segments2]" lines and AnnotationDbi's
# "select() returned 1:many" notes, hundreds per figure; other messages pass.
.pt_quiet <- function(expr) withCallingHandlers(expr, message = function(m) {
  if (grepl("^\\w+\\[\\w+\\]\\s*$|^'select\\(\\)' returned", conditionMessage(m))) invokeRestart("muffleMessage")
})

#' Export a figure as a standalone R script
#'
#' Writes the same self-contained script the app's Download R script button gives: it
#' redraws the figure with plotgardener alone (PlotTrellis need not be installed), and is
#' the place to start for edits the API doesn't offer.
#'
#' @param fig A [pt_figure()].
#' @param file Path for the `.R` script.
#' @return `file`, invisibly.
#' @export
pt_script <- function(fig, file) {
  stopifnot(inherits(fig, "pt_figure"))
  writeLines(assemble_download_script(fig$entries, fig$region_params, fig$settings), file)
  invisible(file)
}

#' @exportS3Method knitr::knit_print
knit_print.pt_figure <- function(x, ...) {
  f <- knitr::fig_path(".png")
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  knitr::knit_print(knitr::include_graphics(pt_save(x, f, dpi = 150)), ...)
}

#' Example samples.txt
#'
#' The path of a `samples.txt` shipped with PlotTrellis, for [run_plottrellis()] or
#' [pt_project()]. Its data is read where it is installed (a project keeps its caches in its
#' own folder), and the app opens on it with a gene and tracks already picked.
#'
#' - `"mammary"` (the default): mouse mammary basal cells (Basal) vs luminal progenitors (LP)
#'   (Milevskiy et al. 2023, GEO GSE227750) at Ehf/Elf5/Cat, mm10: ATAC and transcript counts.
#'   Elf5 is almost silent in Basal and highly expressed in LP. Needs `org.Mm.eg.db`.
#' - `"human"`: ENCODE GM12878 vs K562 at PSME1/PSME2, hg38: ATAC, splice junctions and
#'   transcript counts. Needs `org.Hs.eg.db`.
#'
#' @param dataset `"mammary"` or `"human"`.
#' @return The path of the example's `samples.txt`.
#' @examples
#' pt_example()
#' \dontrun{
#' run_plottrellis(project_dir = "PlotTrellis-example", samples_txt = pt_example())
#' }
#' @export
pt_example <- function(dataset = c("mammary", "human")) {
  dataset <- match.arg(dataset)
  system.file("extdata", if (dataset == "mammary") "example" else "example-human", "samples.txt",
              package = "PlotTrellis", mustWork = TRUE)
}

#' Open the app with a figure preloaded
#'
#' Makes `fig` the starting view of `proj`'s app: the next [run_plottrellis()] on this
#' project opens with the same gene, sample(s) and tracks already picked, so Generate
#' figure is the only click. Use it to hand a ready-to-go project to someone who doesn't
#' code. Saving new project settings in the app's wizard clears it.
#'
#' @param proj A [pt_project()].
#' @param fig A [pt_figure()] of `proj`, made with `gene =` (a `region =` view isn't
#'   preloaded).
#' @return `proj`'s folder, invisibly, for [run_plottrellis()].
#' @examples
#' proj <- pt_project(file.path(tempdir(), "pt-example-human"), pt_example("human"))
#' pt_preload(proj, pt_figure(proj, "PSME2", c("genes", "atac", "junctions"), parallel = TRUE))
#' \dontrun{
#' run_plottrellis(proj)
#' }
#' @export
pt_preload <- function(proj, fig) {
  stopifnot(inherits(proj, "pt_project"), inherits(fig, "pt_figure"))
  if (is.na(fig$gene)) stop("pt_preload() needs a figure made with gene =.", call. = FALSE)
  .pt_start_view(proj$dir, proj$items, fig$keys, fig$gene, names(proj$colours))
  invisible(proj$dir)
}

# The app's mode follows the sample count: 1 = One, 2 = Two (reference first,
# as a comparison track has it), more = Multi, which has no comparison tracks.
.pt_start_view <- function(dir, items, keys, gene, sample_order) {
  sel <- items[keys]
  cmp <- Filter(function(it) !is.null(it$samples), sel)
  conds <- if (length(cmp)) cmp[[1]]$samples else character(0)
  conds <- unique(c(conds, intersect(sample_order, unlist(lapply(sel, `[[`, "sample")))))
  mode <- if (length(conds) == 1) "single" else if (length(conds) == 2) "two" else "compose"
  if (mode == "compose" && length(cmp)) {
    warning("The app's Multi mode (3+ samples) has no comparison tracks; preloading without them.", call. = FALSE)
    keys <- setdiff(keys, names(cmp))
  }
  .pt_write_start_view(dir, mode, if (mode == "compose") character(0) else conds, gene, keys)
}

# Adds start_view (read by .pt_composer_app()'s START_VIEW) to a project's settings.
.pt_write_start_view <- function(dir, view_mode, conditions, gene, keys) {
  json <- file.path(dir, "project_settings.json")
  settings <- jsonlite::read_json(json)
  settings$start_view <- list(view_mode = view_mode, conditions = as.list(conditions), gene = gene,
                              tracks = as.list(keys))
  writeLines(jsonlite::toJSON(settings, auto_unbox = TRUE, pretty = TRUE, null = "null"), json)
}
