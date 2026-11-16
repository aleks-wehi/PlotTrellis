# Code API: open a project once, then describe figures by gene + track names
# and save them. Builds the same entries/settings as the composer's Generate,
# through the shared helpers in R/figure_helpers.R.

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
  p <- .pt_project(project_dir)
  cfg <- p$cfg
  if (!requireNamespace(cfg$orgdb, quietly = TRUE))
    stop("This project needs the Bioconductor package '", cfg$orgdb, "'. Install it with ",
         "BiocManager::install(\"", cfg$orgdb, "\").", call. = FALSE)
  .pt_require_genome_label(cfg$genome_label, p$settings_json)
  genome <- ensure_genome_ready(gtf_path = cfg$gtf_path, orgdb = cfg$orgdb, genome_label = cfg$genome_label,
                                cache_dir = p$cache_dir)
  samples <- load_samples_from_settings(cfg$samples, cfg$data_dir)
  items   <- c(build_track_items(samples), build_comparison_track_items(samples))
  structure(list(
    dir        = normalizePath(project_dir),
    cfg        = cfg,
    genome     = genome,
    txdb       = .pt_annotation(genome)$txdb,
    genes      = sort(AnnotationDbi::keys(get(genome$orgdb, envir = asNamespace(genome$orgdb)), keytype = "SYMBOL")),
    samples    = samples,
    colours    = .pt_sample_colours(samples),
    items      = items,
    tracks     = .pt_track_table(items),
    render_cfg = .pt_render_cfg(cfg, genome, p$cache_dir)
  ), class = "pt_project")
}

# Creates or refreshes project_dir from samples_txt (NULL: the project must
# already exist), shared by pt_project() and run_plottrellis(). A samples.txt
# may carry `# start_gene:` / `# start_tracks:` lines, and optionally `# start_samples:`,
# `# start_comparison:` and `# start_order:` (comma lists, as pt_figure()'s tracks, samples,
# comparison and order; "NA" for an empty samples slot): the app then opens with that
# view preloaded.
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
    # A start line that stops the import (below) must leave the project as it was: the two
    # project files are put back, with their dates, so the next pt_project() tries again
    # and says so, instead of opening the new settings without their start view.
    kept <- file.path(project_dir, c("project_settings.json", "project_samples.txt"))
    backup <- tempfile("pt-import-backup-", fileext = c(".json", ".txt"))
    had <- file.exists(kept)
    file.copy(kept[had], backup[had], copy.date = TRUE)
    on.exit(unlink(backup), add = TRUE)
    # Read before parse_samples_txt(): when samples_txt is the project's own copy (a rebuild
    # from it), that call rewrites the copy, and its start lines must come from the original.
    lines <- readLines(samples_txt, warn = FALSE)
    tryCatch({
      parse_samples_txt(samples_txt, project_dir = project_dir, overwrite = TRUE)
      settings <- jsonlite::read_json(json)
      settings$samples_txt_source <- src
      writeLines(jsonlite::toJSON(settings, auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA), json)
      sv <- .pt_txt_start_view(lines, .pt_project(project_dir)$cfg)
      if (!is.null(sv)) .pt_write_start_view(project_dir, sv$view_mode, unlist(sv$conditions), sv$gene, unlist(sv$tracks))
    }, error = function(e) {
      unlink(kept[!had])
      file.copy(backup[had], kept[had], overwrite = TRUE, copy.date = TRUE)
      stop(e)
    })
  }
  invisible(project_dir)
}

# The start view a samples.txt's start lines give (the lines are described above .pt_import()),
# for project settings `cfg` (its samples and data_dir); NULL without both start_gene and
# start_tracks. Shared by .pt_import() and the wizard's Import.
.pt_txt_start_view <- function(lines, cfg) {
  hdr <- function(key) {
    ln <- grep(paste0("^#\\s*", key, ":"), lines, value = TRUE)
    if (length(ln)) trimws(sub(paste0("^#\\s*", key, ":"), "", ln[1])) else NULL
  }
  gene <- hdr("start_gene")
  tracks <- hdr("start_tracks")
  if (is.null(gene) || is.null(tracks)) return(NULL)
  samples <- load_samples_from_settings(cfg$samples, cfg$data_dir)
  items <- c(build_track_items(samples), build_comparison_track_items(samples))
  split <- function(x) {
    if (is.null(x)) return(NULL)
    v <- trimws(strsplit(x, ",")[[1]])
    v[v == "NA"] <- NA_character_
    v
  }
  tab <- .pt_track_table(items)
  tracks <- split(tracks)
  start_samples <- split(hdr("start_samples"))
  # A track this project doesn't have (e.g. one a newer PlotTrellis removed) is left out of
  # the start view with one warning; the rest still opens. A quoted number or a name@sample
  # keeps its own error (each with its hint).
  known <- vapply(tracks, function(t) grepl("@", t, fixed = TRUE) || grepl("^[0-9]+$", t) ||
                    !inherits(tryCatch(.pt_track_names(tab, t), error = identity), "error"), logical(1))
  if (!all(known)) {
    warning("Leaving out of the start view the track", if (sum(!known) > 1) "s", " this project doesn't have: ",
            paste(tracks[!known], collapse = ", "), " (# start_tracks:).",
            call. = FALSE, immediate. = TRUE)  # immediate: run_plottrellis() would hold it until the app closes
    if (length(start_samples) == length(tracks)) start_samples <- start_samples[known]
    tracks <- tracks[known]
    if (!length(tracks)) return(NULL)
  }
  # "NA" for every track, all of them fixed (a copy written for such a view), means no samples
  # at all; with a per-sample track among them the NA stays, and its "needs a sample" error.
  fixed <- vapply(tracks, function(t) { k <- tab$kind[tab$track == t]; length(k) > 0 && all(k == "fixed") }, logical(1))
  if (!is.null(start_samples) && all(is.na(start_samples)) && all(fixed)) start_samples <- NULL
  keys <- .pt_resolve_tracks(tab, tracks, names(.pt_sample_colours(samples)),
                             start_samples, split(hdr("start_comparison")),
                             hdr("start_order") %||% "track")
  .pt_start_view_of(items, keys, gene, names(.pt_sample_colours(samples)))
}

# Start lines that give back the start view `sv` of project settings `cfg`, checked by reading
# them back with .pt_txt_start_view(); none when no set of them gives exactly that view (a
# Multi view of 2 samples, say). For the project_samples.txt copy, so rebuilding the project
# from it restores the view too.
.pt_start_view_lines <- function(sv, cfg) {
  samples <- load_samples_from_settings(cfg$samples, cfg$data_dir)
  tab <- .pt_track_table(c(build_track_items(samples), build_comparison_track_items(samples)))
  keys <- .pt_renamed_junctions_key(as.character(unlist(sv$tracks)))
  i <- match(keys, tab$key)
  if (!length(keys) || anyNA(i) || !is.character(sv$gene) || length(sv$gene) != 1) return(character(0))
  csv <- function(key, x) paste0("# ", key, ": ", paste(ifelse(is.na(x), "NA", x), collapse = ", "))
  conds <- as.character(unlist(sv$conditions))
  by_track <- csv("start_tracks", unique(tab$track[i]))
  # No start_samples line when no track of the view takes a sample.
  by_sample <- if (all(is.na(tab$sample[i]))) character(0) else csv("start_samples", tab$sample[i])
  forms <- list(c(csv("start_tracks", tab$track[i]), by_sample),
                c(by_track, csv("start_comparison", conds)),
                c(by_track, csv("start_comparison", conds), "# start_order: sample"))
  want <- .pt_start_view_entry(sv$view_mode, conds, sv$gene, keys)
  for (f in forms) {
    lines <- c(paste0("# start_gene: ", sv$gene), f)
    got <- tryCatch(suppressMessages(suppressWarnings(.pt_txt_start_view(lines, cfg))), error = function(e) NULL)
    if (identical(got, want)) return(lines)
  }
  character(0)
}

#' @export
print.pt_project <- function(x, ...) {
  cat("PlotTrellis project:", x$dir, "\n")
  cat("Genome: ", x$genome$genome_label, " (", basename(x$cfg$gtf_path), ", ", x$genome$orgdb, ")\n", sep = "")
  smp <- pt_samples(x)
  cat("Samples (pt_figure(samples = , comparison = ), by name or number):\n")
  cat(sprintf("  %2d %s\n", smp$id, smp$name), sep = "")
  tr <- pt_tracks(x)
  cat("Tracks (pt_figure(tracks = ), by name or number):\n")
  cat(sprintf("  %2d %-22s %s\n", tr$id, tr$track, tr$samples), sep = "")
  invisible(x)
}

# One row per track item: the name pt_figure() accepts, its sample (NA for
# sample-free tracks, "B vs A" for comparisons), kind and TRACK_ITEMS key. Data
# tracks are named after samples.txt's own `type:name` column names, through the
# `api` template of their module's info_<key>().
.pt_track_table <- function(items) {
  rows <- lapply(items, function(it) {
    kind <- if (!is.null(it$samples)) "comparison" else if (is.na(it$sample)) "fixed" else "sample"
    track <- .pt_fill_name(.pt_viz_info(it$type)$api, it$name %||% it$type)
    data.frame(track = track, type = it$type, kind = kind, key = it$key, label = it$label,
               sample = if (kind == "comparison") paste(it$samples[2], "vs", it$samples[1]) else it$sample,
               stringsAsFactors = FALSE)
  })
  tab <- do.call(rbind, rows)
  tab <- tab[order(match(tab$kind, c("fixed", "sample", "comparison"))), ]
  # Two types under one name (e.g. bw:atac and bed:atac): the later types
  # get a _<type> suffix so every name stays unambiguous.
  first_type <- tapply(tab$type, tab$track, `[`, 1)
  dup <- tab$type != first_type[tab$track]
  tab$track[dup] <- paste0(tab$track[dup], "_", tab$type[dup])
  rownames(tab) <- NULL
  tab[, c("track", "sample", "kind", "label", "key")]
}

#' List a project's tracks
#'
#' @param proj A [pt_project()].
#' @return A data frame with one row per track: `id` (its number, usable in
#'   [pt_figure()]'s `tracks` instead of the name), `track` (the name), `kind` (`"fixed"`:
#'   needs no sample; `"sample"`: drawn per sample; `"comparison"`: compares 2 samples),
#'   `label` (the app's label for it) and `samples` (the samples that have it; `"B vs A"`
#'   pairs for comparisons, B compared against reference A).
#' @examples
#' pt_tracks(pt_project(file.path(tempdir(), "pt-example-human"), pt_example("human")))
#' @export
pt_tracks <- function(proj) {
  tab <- proj$tracks
  nm  <- unique(tab$track)
  first <- match(nm, tab$track)
  data.frame(
    id      = seq_along(nm),
    track   = nm,
    kind    = tab$kind[first],
    label   = sub(" \u00B7 .*$", "", tab$label[first]),
    samples = vapply(nm, function(n) {
      s <- tab$sample[tab$track == n]
      if (all(is.na(s))) "" else paste(s, collapse = ", ")
    }, character(1), USE.NAMES = FALSE),
    stringsAsFactors = FALSE
  )
}

#' List a project's samples
#'
#' @param proj A [pt_project()].
#' @return A data frame with one row per sample, in `samples.txt` order: `id` (its number,
#'   usable in [pt_figure()]'s `samples` and `comparison` instead of the name), `name` and
#'   `colour`. Numbers follow the row order, so they change when `samples.txt` does: use
#'   names in scripts you keep.
#' @examples
#' pt_samples(pt_project(file.path(tempdir(), "pt-example-human"), pt_example("human")))
#' @export
pt_samples <- function(proj) {
  nm <- names(proj$colours)
  data.frame(id = seq_along(nm), name = nm, colour = unname(unlist(proj$colours)),
             stringsAsFactors = FALSE)
}

# tracks (names, or pt_tracks() ids) -> one track name each. Old names still work:
# transcript_counts (the quantification track) and <name>_sashimi (<name>_junctions).
.pt_track_names <- function(tab, tracks) {
  all_names <- unique(tab$track)
  listing <- paste(seq_along(all_names), all_names, collapse = ", ")
  if (!length(tracks))
    stop("tracks is empty: name at least one, e.g. tracks = \"genes\". See pt_tracks(proj).", call. = FALSE)
  if (is.numeric(tracks)) {
    bad <- tracks[is.na(tracks) | tracks != round(tracks) | tracks < 1 | tracks > length(all_names)]
    if (length(bad))
      stop("There is no track number ", bad[1], ". This project's tracks: ", listing,
           ". See pt_tracks(proj).", call. = FALSE)
    return(all_names[tracks])
  }
  if (!is.character(tracks) || anyNA(tracks))
    stop("tracks must be track names or numbers, e.g. tracks = c(\"genes\", \"", all_names[length(all_names)],
         "\") or tracks = c(1, ", length(all_names), "). See pt_tracks(proj).", call. = FALSE)
  vapply(tracks, function(spec) {
    name <- trimws(spec)
    if (name %in% all_names) return(name)
    if (identical(name, "transcript_counts") && any(startsWith(tab$key, "tpm_usage_")))
      return(tab$track[startsWith(tab$key, "tpm_usage_")][1])
    if (grepl("_sashimi$", name) && sub("_sashimi$", "_junctions", name) %in% all_names)
      return(sub("_sashimi$", "_junctions", name))
    stop("Unknown track '", name, "'. This project's tracks: ", listing, ". See pt_tracks(proj).",
         if (grepl("^[0-9]+$", name)) paste0(" '", name, "' is in quotes, so it is read as a name: write ",
                                             name, " without quotes (in a samples.txt, use names)."),
         call. = FALSE)
  }, character(1), USE.NAMES = FALSE)
}

# samples / comparison (names exactly as in samples.txt, or their row numbers) -> names.
.pt_sample_names <- function(x, all_samples, arg) {
  if (is.null(x)) return(NULL)
  listing <- paste(seq_along(all_samples), all_samples, collapse = ", ")
  if (is.logical(x) && all(is.na(x))) x <- rep(NA_character_, length(x))
  if (!length(x) || !(is.numeric(x) || is.character(x)))
    stop(arg, " must be sample names or numbers, e.g. ", arg, " = \"", all_samples[1], "\" or ",
         arg, " = 1",
         if (is.logical(x)) " (parallel = TRUE is gone: leave samples and comparison unset to draw every sample)",
         ". This project's samples: ", listing, ". See pt_samples(proj).", call. = FALSE)
  if (is.numeric(x)) {
    ok <- is.na(x) | (x == round(x) & x >= 1 & x <= length(all_samples))
    if (!all(ok))
      stop("There is no sample number ", x[!ok][1], " (in ", arg, "). This project's samples: ",
           listing, ". See pt_samples(proj).", call. = FALSE)
    return(all_samples[x])
  }
  bad <- x[!is.na(x) & !x %in% all_samples]
  if (length(bad))
    stop("Unknown sample '", bad[1], "' (in ", arg, "). Write the name exactly as in samples.txt, ",
         "or its number. This project's samples: ", listing, ". See pt_samples(proj).",
         if (grepl("^[0-9]+$", bad[1])) paste0(" '", bad[1], "' is in quotes, so it is read as a name: write ",
                                               bad[1], " without quotes (in a samples.txt, use names)."),
         call. = FALSE)
  x
}

# A tracks list still in the old "name@sample" form: stops, suggesting the new call only
# when one draws exactly the same tracks (checked by resolving it), else saying how.
.pt_stop_at_form <- function(tab, tracks, all_samples) {
  at  <- grepl("@", tracks, fixed = TRUE)
  nm  <- trimws(sub("@.*$", "", tracks))
  smp <- ifelse(at, gsub("\\s+", " ", trimws(sub("^[^@]*@", "", tracks))), NA_character_)
  kind <- tab$kind[match(nm, tab$track)]
  old <- vapply(seq_along(nm), function(i) {
    if (!at[i] && identical(kind[i], "fixed")) return(tab$key[tab$track == nm[i]][1])
    if (!at[i]) return(NA_character_)
    tab$key[tab$track == nm[i] & !is.na(tab$sample) & tab$sample == smp[i]][1]
  }, character(1))
  q <- function(x) paste0("c(", paste(ifelse(is.na(x), "NA", vapply(x, deparse, "")), collapse = ", "), ")")
  cands <- list()
  if (!anyNA(old)) {
    if (!any(kind == "comparison")) cands <- list(list(tracks = nm, samples = smp, order = "track"))
    pair <- unique(smp[kind == "comparison"])
    if (length(pair) == 1) {
      cmp <- c(sub(".* vs ", "", pair), sub(" vs .*", "", pair))
      # Each per-sample track once per run of adjacent ones (comparison = repeats it).
      run <- cumsum(kind != "sample")
      keep <- kind != "sample" | !duplicated(paste(run, nm))
      for (ord in c("track", "sample"))
        cands <- c(cands, list(list(tracks = nm[keep], comparison = cmp, order = ord)))
    }
  }
  same <- function(a) identical(old, tryCatch(
    suppressMessages(.pt_resolve_tracks(tab, a$tracks, all_samples, a$samples, a$comparison, a$order)),
    error = function(e) NULL))
  ok <- lapply(Filter(same, cands), function(a)
    paste0("tracks = ", q(a$tracks),
           if (!is.null(a$samples)) paste0(", samples = ", q(a$samples)),
           if (!is.null(a$comparison)) paste0(", comparison = ", q(a$comparison)),
           if (a$order == "sample") ", order = \"sample\""))
  stop("'", tracks[at][1], "': the name@sample form is gone. Write track names alone and give the ",
       "samples separately",
       if (length(ok)) paste0(". The same figure: ", ok[[1]]) else
         paste0(": samples = one sample for every track, or one per track (NA for tracks without ",
                "one); comparison = c(reference, other) for comparison tracks, which then draws ",
                "per-sample tracks for both. See ?pt_figure"),
       ". In a samples.txt, use # start_tracks: with # start_samples: or # start_comparison: lines.",
       call. = FALSE)
}

# tracks + samples/comparison/order -> TRACK_ITEMS keys, top to bottom. See pt_figure()
# for the rules. A per-sample track repeated over samples (no `samples`) skips, with a
# message, a sample that has no such track; one picked by `samples` errors instead.
.pt_resolve_tracks <- function(tab, tracks, all_samples, samples = NULL, comparison = NULL, order = "track") {
  if (!(is.character(order) && length(order) == 1 && order %in% c("track", "sample")))
    stop("order must be \"track\" (all of one track, then the next) or \"sample\" ",
         "(all of one sample's tracks, then the next sample's).", call. = FALSE)
  if (!is.null(samples) && !is.null(comparison))
    stop("Set samples or comparison, not both: samples = picks each track's sample; comparison = ",
         "repeats every per-sample track over its samples, and gives comparison tracks their 2.", call. = FALSE)
  if (is.character(tracks) && any(grepl("@", tracks, fixed = TRUE))) .pt_stop_at_form(tab, tracks, all_samples)
  nms   <- .pt_track_names(tab, tracks)
  kinds <- tab$kind[match(nms, tab$track)]
  samples    <- .pt_sample_names(samples, all_samples, "samples")
  comparison <- .pt_sample_names(comparison, all_samples, "comparison")
  if (anyNA(comparison))
    stop("comparison can't contain NA: list the samples to draw, e.g. comparison = c(\"",
         all_samples[1], "\", \"", all_samples[min(2, length(all_samples))], "\").", call. = FALSE)
  if (anyDuplicated(comparison))
    stop("comparison lists '", comparison[duplicated(comparison)][1], "' twice; list each sample once.", call. = FALSE)
  has_of <- function(n) tab$sample[tab$track == n]
  key_of <- function(n, s) tab$key[tab$track == n & !is.na(tab$sample) & tab$sample == s][1]

  cmp <- unique(nms[kinds == "comparison"])
  if (length(cmp)) {
    ex  <- has_of(cmp[1])[1]
    arg <- paste0("comparison = c(\"", sub(".* vs ", "", ex), "\", \"", sub(" vs .*", "", ex), "\")")
    if (!is.null(samples))
      stop("'", cmp[1], "' compares 2 samples, which samples = can't give: drop samples and set ",
           arg, " (reference first).", call. = FALSE)
    if (is.null(comparison))
      stop("'", cmp[1], "' compares 2 samples: set ", arg, " (reference first).", call. = FALSE)
    if (length(comparison) != 2)
      stop("'", cmp[1], "' compares 2 samples, so comparison needs exactly 2 (reference first); it has ",
           length(comparison), ".", call. = FALSE)
  }

  if (!is.null(samples)) {
    n <- length(nms)
    if (length(samples) == 1) {
      slot <- ifelse(kinds == "fixed", NA_character_, samples)
    } else if (length(samples) == n) {
      slot <- samples
      bad <- which(kinds == "fixed" & !is.na(slot))
      if (length(bad))
        stop("samples[", bad[1], "] is '", slot[bad[1]], "', but track ", bad[1], " ('", nms[bad[1]],
             "') has no sample: put NA there.", call. = FALSE)
    } else {
      stop("samples has ", length(samples), " entries for ", n, if (n == 1) " track: give 1." else
             paste0(" tracks: give 1 (that sample for every track) or ", n, " (one per track, NA for ",
                    "tracks without a sample, such as genes)."), call. = FALSE)
    }
    keys <- character(n)
    for (i in seq_len(n)) {
      if (kinds[i] == "fixed") { keys[i] <- tab$key[tab$track == nms[i]][1]; next }
      if (is.na(slot[i]))
        stop("Track ", i, " ('", nms[i], "') needs a sample, but samples[", if (length(samples) == 1) 1 else i,
             "] is NA.", call. = FALSE)
      keys[i] <- key_of(nms[i], slot[i])
      if (is.na(keys[i]))
        stop("No '", nms[i], "' track for '", slot[i], "'. Samples that have it: ",
             paste(has_of(nms[i]), collapse = ", "), ".", call. = FALSE)
    }
    if (all(kinds == "fixed"))
      warning("samples has no effect: none of these tracks needs a sample (", paste(unique(nms), collapse = ", "),
              "). It picks the sample of per-sample tracks; see pt_tracks().", call. = FALSE)
    return(keys)
  }

  if (!is.null(comparison) && all(kinds == "fixed"))
    warning("comparison has no effect: none of these tracks needs a sample (", paste(unique(nms), collapse = ", "),
            "). It picks the samples of per-sample and comparison tracks; see pt_tracks().", call. = FALSE)
  over <- comparison %||% all_samples
  one <- function(i, s) {
    k <- key_of(nms[i], s)
    if (is.na(k)) { message("Note: no '", nms[i], "' for ", s, ", skipped."); return(character(0)) }
    k
  }
  keys <- character(0)
  i <- 1L
  while (i <= length(nms)) {
    if (kinds[i] == "fixed") {
      keys <- c(keys, tab$key[tab$track == nms[i]][1])
    } else if (kinds[i] == "comparison") {
      pair <- paste(comparison[2], "vs", comparison[1])
      k <- key_of(nms[i], pair)
      if (is.na(k))
        stop("No '", nms[i], "' track for ", pair, ". Pairs that have it: ",
             paste(has_of(nms[i]), collapse = ", "), ".", call. = FALSE)
      keys <- c(keys, k)
    } else {
      # A run of per-sample tracks is one block: by track (each track over every sample)
      # or by sample (each sample's tracks, then the next sample's).
      j <- i
      while (j < length(nms) && kinds[j + 1] == "sample") j <- j + 1L
      block <- i:j
      if (order == "track") for (b in block) for (s in over) keys <- c(keys, one(b, s))
      else for (s in over) for (b in block) keys <- c(keys, one(b, s))
      i <- j
    }
    i <- i + 1L
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
#' `quantdir:` gives `transcript_expression` (`transcript_counts`, its old name, still works);
#' comparisons get a `_diff` suffix, such as `junctions_diff`).
#'
#' Which samples a per-sample track shows:
#' * neither `samples` nor `comparison`: every sample, in `samples.txt` order;
#' * `samples = "K562"`: that one sample for every track;
#' * `samples` as long as `tracks`: track i shows sample i, with `NA` for tracks that need
#'   no sample, e.g. `tracks = c("genes", "atac", "atac"), samples = c(NA, "K562", "GM12878")`;
#' * `comparison = c("GM12878", "K562")`: every per-sample track for those samples, in that
#'   order. A comparison track (such as `junctions_diff`) needs exactly 2 here and compares
#'   the second against the first, the reference.
#'
#' Tracks and samples can also be given by number: tracks by their `id` in [pt_tracks()],
#' samples by their row in `samples.txt` ([pt_samples()]). Numbers change when
#' `samples.txt` does, so use names in scripts you keep. Sample names must be written
#' exactly as in `samples.txt`.
#'
#' @param proj A [pt_project()].
#' @param gene A gene symbol (any case) or Ensembl gene ID. The view is the gene plus the
#'   project's flank setting on each side.
#' @param tracks Track names (or numbers, see [pt_tracks()]), top to bottom.
#' @param samples Which sample each per-sample track shows: one name (or number) for all of
#'   them, or one per track (`NA` for tracks without a sample). `NULL` (the default): every
#'   sample.
#' @param comparison Instead of `samples`: the samples to repeat every per-sample track over,
#'   in order; with a comparison track, exactly 2, the reference first.
#' @param region Instead of `gene`: a region such as `"chr1:198,600,000-198,760,000"`.
#' @param title Figure title.
#' @param order When a track is repeated over samples: `"track"` (all samples of one track,
#'   then the next track) or `"sample"` (all of one sample's tracks, then the next sample's).
#'   Grouping applies to per-sample tracks written next to each other: tracks without a
#'   sample and comparison tracks stay where they are written and split the groups, so write
#'   them first or last.
#' @param ... Not used; an old `parallel =` stops with what to write instead.
#' @return A `pt_figure`.
#' @seealso [pt_save()], [pt_script()], [pt_tracks()], [pt_samples()]
#' @examples
#' proj <- pt_project(file.path(tempdir(), "pt-example-human"), pt_example("human"))
#' fig <- pt_figure(proj, "PSME2", c("genes", "transcripts", "atac", "junctions"))
#' pt_save(fig, file.path(tempdir(), "PSME2.pdf"))
#'
#' # One sample for every track:
#' pt_figure(proj, "PSME2", c("genes", "atac", "junctions"), samples = "K562")
#'
#' # One sample per track, in any order:
#' pt_figure(proj, "PSME2", c("genes", "atac", "atac"), samples = c(NA, "K562", "GM12878"))
#'
#' # Two samples and their comparison (K562 against the reference GM12878), grouped by sample:
#' pt_figure(proj, "PSME2", c("atac", "junctions", "junctions_diff", "genes"),
#'           comparison = c("GM12878", "K562"), order = "sample")
#' @export
pt_figure <- function(proj, gene = NULL, tracks = "genes", samples = NULL, comparison = NULL,
                      region = NULL, title = "", order = "track", ...) {
  stopifnot(inherits(proj, "pt_project"))
  dots <- ...names() %||% rep("", ...length())
  if ("parallel" %in% dots)
    stop("parallel is gone: leave samples and comparison unset to draw every per-sample track for ",
         "every sample, or set comparison = c(...) to pick the samples and their order.", call. = FALSE)
  if (length(dots))
    stop("pt_figure() ",
         if (any(nzchar(dots))) paste0("has no argument ", paste0("'", dots[nzchar(dots)], "'", collapse = ", ")),
         if (any(nzchar(dots)) && !all(nzchar(dots))) ", and ",
         if (!all(nzchar(dots))) "takes 8 arguments; this call gives more",
         ". See ?pt_figure.", call. = FALSE)
  proj <- .pt_renamed_junctions_project(proj)

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
    err <- .pt_region_error(rgn$start, rgn$end)
    if (!is.null(err)) stop(err, call. = FALSE)
  } else {
    stop("Pass gene = (a symbol or Ensembl ID) or region = (\"chr:start-end\").", call. = FALSE)
  }

  keys <- .pt_resolve_tracks(proj$tracks, tracks, names(proj$colours), samples, comparison, order)
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
  fig$entries <- lapply(fig$entries, function(e) replace(e, "viz", .pt_renamed_junctions_key(e$viz)))
  dev <- switch(tolower(tools::file_ext(file)),
    svg = function(w, h) svglite::svglite(file, width = w, height = h),
    pdf = function(w, h) grDevices::cairo_pdf(file, width = w, height = h),
    png = function(w, h) grDevices::png(file, width = w, height = h, units = "in", res = dpi),
    stop("pt_save() writes .pdf, .svg or .png files.", call. = FALSE))
  # Closes only the device it opened (its number), then makes the caller's current
  # device current again.
  before <- grDevices::dev.cur()
  opened <- NULL
  on.exit({
    if (!is.null(opened) && opened %in% grDevices::dev.list()) grDevices::dev.off(opened)
    if (before %in% grDevices::dev.list()) grDevices::dev.set(before)
  })
  .pt_quiet(render_all(fig$entries, fig$region_params, fig$settings,
                       device_fn = function(w, h) { dev(w, h); opened <<- grDevices::dev.cur() }))
  if (!is.null(opened)) grDevices::dev.off(opened)
  opened <- NULL
  if (tolower(tools::file_ext(file)) == "svg") .pt_tidy_svg(file)
  invisible(file)
}

# Hides plotgardener's per-grob "segments[segments2]" lines, hundreds per figure,
# for the API and the app; other messages pass. (AnnotationDbi's "select()
# returned 1:many" notes are silenced at each select() call; the pattern stays
# as a backstop.)
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
  fig$entries <- lapply(fig$entries, function(e) replace(e, "viz", .pt_renamed_junctions_key(e$viz)))
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
#'   (Milevskiy et al. 2023, GEO GSE227750) at Ehf/Elf5/Cat, mm10: ATAC and transcript expression.
#'   Elf5 is almost silent in Basal and highly expressed in LP. Needs `org.Mm.eg.db`.
#' - `"human"`: ENCODE GM12878 vs K562 at PSME1/PSME2, hg38: ATAC, splice junctions and
#'   transcript expression. Needs `org.Hs.eg.db`.
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
#' code. Saving the app's wizard keeps it while the GTF, the OrgDb and the samples and
#' tracks it shows (by name) are all still there, and clears it otherwise. An Import in
#' the wizard replaces it by the imported samples.txt's start lines (none: no view); the
#' project_samples.txt copy beside the settings carries it as such lines.
#'
#' @param proj A [pt_project()].
#' @param fig A [pt_figure()] of `proj`, made with `gene =` (a `region =` view isn't
#'   preloaded).
#' @return `proj`'s folder, invisibly, for [run_plottrellis()].
#' @examples
#' proj <- pt_project(file.path(tempdir(), "pt-example-human"), pt_example("human"))
#' pt_preload(proj, pt_figure(proj, "PSME2", c("genes", "atac", "junctions")))
#' \dontrun{
#' run_plottrellis(proj)
#' }
#' @export
pt_preload <- function(proj, fig) {
  stopifnot(inherits(proj, "pt_project"), inherits(fig, "pt_figure"))
  if (is.na(fig$gene)) stop("pt_preload() needs a figure made with gene =.", call. = FALSE)
  proj <- .pt_renamed_junctions_project(proj)
  .pt_start_view(proj$dir, proj$items, .pt_renamed_junctions_key(fig$keys), fig$gene, names(proj$colours))
  invisible(proj$dir)
}

# The app's mode follows the sample count: 1 = One, 2 = Two (reference first,
# as a comparison track has it), more = Multi, which has no comparison tracks.
.pt_start_view <- function(dir, items, keys, gene, sample_order) {
  sv <- .pt_start_view_of(items, keys, gene, sample_order)
  .pt_write_start_view(dir, sv$view_mode, unlist(sv$conditions), sv$gene, unlist(sv$tracks))
}

# The start view itself, as project_settings.json holds it (see .pt_start_view()).
.pt_start_view_of <- function(items, keys, gene, sample_order) {
  sel <- items[keys]
  cmp <- Filter(function(it) !is.null(it$samples), sel)
  conds <- if (length(cmp)) cmp[[1]]$samples else character(0)
  conds <- unique(c(conds, intersect(sample_order, unlist(lapply(sel, `[[`, "sample")))))
  mode <- if (length(conds) == 1) "single" else if (length(conds) == 2) "two" else "compose"
  if (mode == "compose" && length(cmp)) {
    warning("The app's Multi mode (3+ samples) has no comparison tracks; preloading without them.", call. = FALSE)
    keys <- keys[!keys %in% names(cmp)]
  }
  .pt_start_view_entry(mode, if (mode == "compose") character(0) else conds, gene, keys)
}

.pt_start_view_entry <- function(view_mode, conditions, gene, keys)
  list(view_mode = view_mode, conditions = as.list(conditions), gene = gene, tracks = as.list(keys))

# Adds start_view (read by .pt_composer_app()'s START_VIEW) to a project's settings.
.pt_write_start_view <- function(dir, view_mode, conditions, gene, keys) {
  json <- file.path(dir, "project_settings.json")
  settings <- jsonlite::read_json(json)
  settings$start_view <- .pt_start_view_entry(view_mode, conditions, gene, keys)
  writeLines(jsonlite::toJSON(settings, auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA), json)
  # And the project_samples.txt copy, which then carries the view too (a failure is ignored).
  .pt_write_samples_copy(settings, json)
  invisible()
}

# The junction track's viz key used to be "sashimi". A saved start view can still hold
# its item keys ("sashimi_<name>_<sample>") and a saved pt_figure its entries' viz
# ("sashimi"); this gives their "junctions" equivalents. NULL (no tracks) stays NULL.
.pt_renamed_junctions_key <- function(keys) {
  if (is.null(keys)) return(NULL)
  sub("^sashimi(_|$)", "junctions\\1", keys)
}

# A pt_project made before that rename and kept (saveRDS()), with its track items and
# track table keyed as they are made now. A project made since comes back as it is.
.pt_renamed_junctions_project <- function(proj) {
  if (!any(grepl("^sashimi(_|$)", names(proj$items)))) return(proj)
  proj$items <- lapply(proj$items, function(it) {
    it$key  <- .pt_renamed_junctions_key(it$key)
    it$type <- .pt_renamed_junctions_key(it$type)
    it
  })
  names(proj$items) <- .pt_renamed_junctions_key(names(proj$items))
  proj$tracks$key <- .pt_renamed_junctions_key(proj$tracks$key)
  proj
}
