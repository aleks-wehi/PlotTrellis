# Code API: open a project once, then describe figures by gene + track names
# and save them. Builds the same entries/settings as the composer's Generate,
# through the shared helpers in R/api-figure-helpers.R.

# The bare-call project folder: visible, in the home directory, one per user.
.pt_default_dir <- function() file.path(path.expand("~"), "PlotTrellis")

#' Open a PlotTrellis project for scripting
#'
#' Loads a project's genome, samples and tracks once, so [pt_figure()] calls are
#' quick. With a `samples_txt`, it first creates the project in `project_dir`
#' from it (the folder is created if needed); without one, it opens the project
#' already in `project_dir`.
#'
#' A project is re-imported automatically when its `samples.txt` is newer than
#' the project's `project_settings.json` (you edited it), so re-running a script
#' after editing `samples.txt` picks the edit up. Otherwise the existing project
#' is opened as it is, keeping any settings changed in the app. One folder holds
#' one project: a different `samples.txt` pointed at a folder that already has a
#' project is an error.
#'
#' @param project_dir The project's folder: its settings and caches live here.
#'   Defaults to the `PlotTrellis` folder in your home directory, as for
#'   [run_plottrellis()].
#' @param samples_txt Path to a `samples.txt` describing the samples and genome,
#'   e.g. [pt_example()]. Leave out to open an existing project.
#' @param overwrite If `TRUE`, re-import `samples_txt` even if the project is
#'   newer (discarding settings changed in the app since), or replace a project
#'   made from a different `samples.txt`.
#' @return A `pt_project`. Print it to see its samples and the track names
#'   [pt_figure()] accepts.
#' @seealso [pt_figure()], [pt_tracks()], [pt_example()]
#' @examples
#' proj <- pt_project(file.path(tempdir(), "pt-example-human"),
#'                    pt_example("human"))
#' proj
#' @export
pt_project <- function(project_dir = .pt_default_dir(), samples_txt = NULL,
                       overwrite = FALSE) {
  .pt_import(project_dir, samples_txt, overwrite)
  p <- .pt_project(project_dir)
  cfg <- p$cfg
  if (!requireNamespace(cfg$orgdb, quietly = TRUE))
    stop("This project needs the Bioconductor package '", cfg$orgdb,
         "'. Install it with ",
         "BiocManager::install(\"", cfg$orgdb, "\").", call. = FALSE)
  .pt_require_genome_label(cfg$genome_label, p$settings_json)
  genome <- ensure_genome_ready(gtf_path = cfg$gtf_path, orgdb = cfg$orgdb,
                                genome_label = cfg$genome_label,
                                cache_dir = p$cache_dir)
  samples <- load_samples_from_settings(cfg$samples, cfg$data_dir)
  items   <- c(build_track_items(samples),
               build_comparison_track_items(samples))
  structure(list(
    dir        = normalizePath(project_dir),
    cfg        = cfg,
    genome     = genome,
    txdb       = .pt_annotation(genome)$txdb,
    genes      = sort(AnnotationDbi::keys(
      get(genome$orgdb, envir = asNamespace(genome$orgdb)), keytype = "SYMBOL"
    )),
    samples    = samples,
    colours    = .pt_sample_colours(samples),
    items      = items,
    tracks     = .pt_track_table(items),
    render_cfg = .pt_render_cfg(cfg, genome, p$cache_dir)
  ), class = "pt_project")
}

# Creates or refreshes project_dir from samples_txt (NULL: the project must
# already exist), shared by pt_project() and run_plottrellis(). The
# `# start_gene:`, `# start_tracks:`, ... lines an older samples.txt may carry
# (the app once opened on that view) are ignored, with one message.
.pt_import <- function(project_dir, samples_txt = NULL, overwrite = FALSE) {
  json <- file.path(project_dir, "project_settings.json")
  if (is.null(samples_txt)) {
    if (!file.exists(json))
      stop("No PlotTrellis project in '", project_dir,
           "' yet. Create one with ",
           "pt_project(project_dir, samples_txt = \"path/to/samples.txt\"), or run_plottrellis() to set ",
           "it up in the app.", call. = FALSE)
    return(invisible(project_dir))
  }
  if (!is_real_file(samples_txt))
    stop("No such samples.txt: '", samples_txt, "'.", call. = FALSE)
  src <- normalizePath(samples_txt)
  # One folder holds one project: a different samples.txt (or a project set up
  # in the app) already there is an error, not silently opened or replaced.
  old <- if (file.exists(json)) jsonlite::read_json(json)
  made_from <- if (!is.null(old))
    old$samples_txt_source %||% "the app's setup wizard"
  # The same bundled example from another PlotTrellis install (an R upgrade or
  # reinstall moves pt_example()'s path): the same project, re-imported so its
  # data paths follow the example to this install.
  moved <- !is.null(made_from) && !identical(made_from, src) &&
    .pt_same_example(made_from, src)
  if (!is.null(made_from) && !identical(made_from, src) && !moved &&
      !isTRUE(overwrite))
    stop("'", project_dir, "' already holds a project, made from ", made_from,
         ". Use another ",
         "project_dir for a separate project, or overwrite = TRUE to replace it.",
         call. = FALSE)
  # Edited since its import: newer than the time the import recorded (both
  # rounded alike). The settings' own time can't tell, since the app's Save as
  # project default rewrites them; a project imported before that record (or
  # saved by the wizard, which drops it) still compares with it.
  stamp <- function() sprintf("%.6f", as.numeric(file.mtime(samples_txt)))
  edited <- function() {
    if (is.character(old$samples_txt_mtime))
      as.numeric(stamp()) > as.numeric(old$samples_txt_mtime)
    else file.mtime(samples_txt) > file.mtime(json)
  }
  if (is.null(made_from) || isTRUE(overwrite) || moved || edited()) {
    if (!is.null(made_from))
      message("Re-importing '", basename(samples_txt), "' into '", project_dir,
              "'.")
    # An import that stops part way must leave the project as it was: the two
    # project files are put back, with their dates, so the next pt_project()
    # tries again and says so.
    kept <- file.path(project_dir,
                      c("project_settings.json", "project_samples.txt"))
    backup <- tempfile("pt-import-backup-", fileext = c(".json", ".txt"))
    had <- file.exists(kept)
    file.copy(kept[had], backup[had], copy.date = TRUE)
    on.exit(unlink(backup), add = TRUE)
    if (any(grepl("^#\\s*start_(gene|tracks|samples|comparison|order):",
                  readLines(samples_txt, warn = FALSE))))
      message("Ignoring the '# start_' lines of '", basename(samples_txt),
              "': the app opens with nothing picked. ",
              "A '# suggested_genes:' line lists genes for its gene box to suggest.")
    tryCatch({
      parse_samples_txt(samples_txt, project_dir = project_dir,
                        overwrite = TRUE)
      settings <- jsonlite::read_json(json)
      settings$samples_txt_source <- src
      # Taken after parse_samples_txt(): when samples_txt is the project's own
      # copy, that call rewrites it.
      settings$samples_txt_mtime <- stamp()
      writeLines(jsonlite::toJSON(settings, auto_unbox = TRUE, pretty = TRUE,
                                  null = "null", digits = NA),
                 json)
    }, error = function(e) {
      unlink(kept[!had])
      file.copy(backup[had], kept[had], overwrite = TRUE, copy.date = TRUE)
      stop(e)
    })
  }
  invisible(project_dir)
}

# TRUE when `old` and `new` (normalized paths) are the samples.txt of the same
# bundled example, `new` in this PlotTrellis install and `old` in any install
# or source tree of the package (.../PlotTrellis/extdata/<example>/, or
# .../PlotTrellis-pkg/inst/extdata/<example>/ in a source tree).
.pt_same_example <- function(old, new) {
  old <- gsub("\\\\", "/", old)
  new <- gsub("\\\\", "/", new)
  ext <- normalizePath(system.file("extdata", package = "PlotTrellis"),
                       winslash = "/", mustWork = FALSE)
  if (!nzchar(ext) || !startsWith(new, paste0(ext, "/"))) return(FALSE)
  example <- sub("/.*$", "", substring(new, nchar(ext) + 2))
  identical(new, paste0(ext, "/", example, "/samples.txt")) &&
    grepl(paste0("/PlotTrellis[^/]*(/inst)?/extdata/\\Q", example,
                 "\\E/samples\\.txt$"), old, perl = TRUE)
}

#' @export
print.pt_project <- function(x, ...) {
  cat("PlotTrellis project:", x$dir, "\n")
  cat("Genome: ", x$genome$genome_label, " (", basename(x$cfg$gtf_path), ", ",
      x$genome$orgdb, ")\n", sep = "")
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
    kind <- if (!is.null(it$samples))
      "comparison" else if (is.na(it$sample)) "fixed" else "sample"
    track <- .pt_fill_name(.pt_viz_info(it$type)$api, it$name %||% it$type)
    data.frame(track = track, type = it$type, kind = kind, key = it$key,
               label = it$label,
               sample = if (kind == "comparison")
                 paste(it$samples[2], "vs", it$samples[1]) else it$sample,
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
#'   [pt_figure()]'s `tracks` instead of the name), `track` (the name), `kind`
#'   (`"fixed"`: needs no sample; `"sample"`: drawn per sample; `"comparison"`:
#'   compares 2 samples), `label` (the app's label for it) and `samples` (the
#'   samples that have it; `"B vs A"` pairs for comparisons, B compared against
#'   reference A).
#' @examples
#' pt_tracks(pt_project(file.path(tempdir(), "pt-example-human"),
#'                      pt_example("human")))
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
#' @return A data frame with one row per sample, in `samples.txt` order: `id`
#'   (its number, usable in [pt_figure()]'s `samples` and `comparison` instead
#'   of the name), `name` and `colour`. Numbers follow the row order, so they
#'   change when `samples.txt` does: use names in scripts you keep.
#' @examples
#' pt_samples(pt_project(file.path(tempdir(), "pt-example-human"),
#'                       pt_example("human")))
#' @export
pt_samples <- function(proj) {
  nm <- names(proj$colours)
  data.frame(id = seq_along(nm), name = nm,
             colour = unname(unlist(proj$colours)),
             stringsAsFactors = FALSE)
}

# tracks (names, or pt_tracks() ids) -> one track name each. Old names still
# work: transcript_counts (the quantification track) and <name>_sashimi
# (<name>_junctions).
.pt_track_names <- function(tab, tracks) {
  all_names <- unique(tab$track)
  listing <- paste(seq_along(all_names), all_names, collapse = ", ")
  if (!length(tracks))
    stop("tracks is empty: name at least one, e.g. tracks = \"genes\". See pt_tracks(proj).",
         call. = FALSE)
  if (is.numeric(tracks)) {
    bad <- tracks[is.na(tracks) | tracks != round(tracks) | tracks < 1 |
                    tracks > length(all_names)]
    if (length(bad))
      stop("There is no track number ", bad[1],
           ". This project's tracks: ", listing,
           ". See pt_tracks(proj).", call. = FALSE)
    return(all_names[tracks])
  }
  if (!is.character(tracks) || anyNA(tracks))
    stop("tracks must be track names or numbers, e.g. tracks = c(\"genes\", \"",
         all_names[length(all_names)],
         "\") or tracks = c(1, ", length(all_names), "). See pt_tracks(proj).",
         call. = FALSE)
  vapply(tracks, function(spec) {
    name <- trimws(spec)
    if (name %in% all_names) return(name)
    if (identical(name, "transcript_counts") &&
        any(startsWith(tab$key, "tpm_usage_")))
      return(tab$track[startsWith(tab$key, "tpm_usage_")][1])
    if (grepl("_sashimi$", name) &&
        sub("_sashimi$", "_junctions", name) %in% all_names)
      return(sub("_sashimi$", "_junctions", name))
    stop("Unknown track '", name, "'. This project's tracks: ", listing,
         ". See pt_tracks(proj).",
         if (grepl("^[0-9]+$", name))
           paste0(" '", name,
                  "' is in quotes, so it is read as a name: write ",
                  name, " without quotes (in a samples.txt, use names)."),
         call. = FALSE)
  }, character(1), USE.NAMES = FALSE)
}

# samples / comparison (names exactly as in samples.txt, or their row numbers)
# -> names.
.pt_sample_names <- function(x, all_samples, arg) {
  if (is.null(x)) return(NULL)
  listing <- paste(seq_along(all_samples), all_samples, collapse = ", ")
  if (is.logical(x) && all(is.na(x))) x <- rep(NA_character_, length(x))
  if (!length(x) || !(is.numeric(x) || is.character(x)))
    stop(arg, " must be sample names or numbers, e.g. ", arg, " = \"",
         all_samples[1], "\" or ",
         arg, " = 1",
         if (is.logical(x))
           " (parallel = TRUE is gone: leave samples and comparison unset to draw every sample)",
         ". This project's samples: ", listing, ". See pt_samples(proj).",
         call. = FALSE)
  if (is.numeric(x)) {
    ok <- is.na(x) | (x == round(x) & x >= 1 & x <= length(all_samples))
    if (!all(ok))
      stop("There is no sample number ", x[!ok][1], " (in ", arg,
           "). This project's samples: ",
           listing, ". See pt_samples(proj).", call. = FALSE)
    return(all_samples[x])
  }
  bad <- x[!is.na(x) & !x %in% all_samples]
  if (length(bad))
    stop("Unknown sample '", bad[1], "' (in ", arg,
         "). Write the name exactly as in samples.txt, ",
         "or its number. This project's samples: ", listing,
         ". See pt_samples(proj).",
         if (grepl("^[0-9]+$", bad[1]))
           paste0(" '", bad[1],
                  "' is in quotes, so it is read as a name: write ",
                  bad[1], " without quotes (in a samples.txt, use names)."),
         call. = FALSE)
  x
}

# A tracks list still in the old "name@sample" form: stops, suggesting the new
# call only when one draws exactly the same tracks (checked by resolving it),
# else saying how.
.pt_stop_at_form <- function(tab, tracks, all_samples) {
  at  <- grepl("@", tracks, fixed = TRUE)
  nm  <- trimws(sub("@.*$", "", tracks))
  smp <- ifelse(at, gsub("\\s+", " ", trimws(sub("^[^@]*@", "", tracks))),
                NA_character_)
  kind <- tab$kind[match(nm, tab$track)]
  old <- vapply(seq_along(nm), function(i) {
    if (!at[i] && identical(kind[i], "fixed"))
      return(tab$key[tab$track == nm[i]][1])
    if (!at[i]) return(NA_character_)
    tab$key[tab$track == nm[i] & !is.na(tab$sample) & tab$sample == smp[i]][1]
  }, character(1))
  q <- function(x) paste0("c(", paste(ifelse(is.na(x), "NA",
                                             vapply(x, deparse, "")),
                                      collapse = ", "), ")")
  cands <- list()
  if (!anyNA(old)) {
    if (!any(kind == "comparison"))
      cands <- list(list(tracks = nm, samples = smp, order = "track"))
    pair <- unique(smp[kind == "comparison"])
    if (length(pair) == 1) {
      cmp <- c(sub(" vs .*", "", pair), sub(".* vs ", "", pair))
      # Each per-sample track once per run of adjacent ones (comparison =
      # repeats it).
      run <- cumsum(kind != "sample")
      keep <- kind != "sample" | !duplicated(paste(run, nm))
      for (ord in c("track", "sample"))
        cands <- c(cands, list(list(tracks = nm[keep], comparison = cmp,
                                    order = ord)))
    }
  }
  same <- function(a) identical(old, tryCatch(
    suppressMessages(.pt_resolve_tracks(tab, a$tracks, all_samples, a$samples,
                                        a$comparison, a$order)),
    error = function(e) NULL))
  ok <- lapply(Filter(same, cands), function(a)
    paste0("tracks = ", q(a$tracks),
           if (!is.null(a$samples)) paste0(", samples = ", q(a$samples)),
           if (!is.null(a$comparison))
             paste0(", comparison = ", q(a$comparison)),
           if (a$order == "sample") ", order = \"sample\""))
  stop("'", tracks[at][1],
       "': the name@sample form is gone. Write track names alone and give the ",
       "samples separately",
       if (length(ok)) paste0(". The same figure: ", ok[[1]]) else
         paste0(": samples = one sample for every track, or one per track (NA for tracks without ",
                "one); comparison = c(other, reference) for comparison tracks, which then draws ",
                "per-sample tracks for both. See ?pt_figure"),
       ".",
       call. = FALSE)
}

# tracks + samples/comparison/order -> TRACK_ITEMS keys, top to bottom. See
# pt_figure() for the rules. A per-sample track repeated over samples (no
# `samples`) skips, with a message, a sample that has no such track; one picked
# by `samples` errors instead.
.pt_resolve_tracks <- function(tab, tracks, all_samples, samples = NULL,
                               comparison = NULL, order = "track") {
  if (!(is.character(order) && length(order) == 1 &&
        order %in% c("track", "sample")))
    stop("order must be \"track\" (all of one track, then the next) or \"sample\" ",
         "(all of one sample's tracks, then the next sample's).", call. = FALSE)
  if (!is.null(samples) && !is.null(comparison))
    stop("Set samples or comparison, not both: samples = picks each track's sample; comparison = ",
         "repeats every per-sample track over its samples, and gives comparison tracks their 2.",
         call. = FALSE)
  if (is.character(tracks) && any(grepl("@", tracks, fixed = TRUE)))
    .pt_stop_at_form(tab, tracks, all_samples)
  nms   <- .pt_track_names(tab, tracks)
  kinds <- tab$kind[match(nms, tab$track)]
  samples    <- .pt_sample_names(samples, all_samples, "samples")
  comparison <- .pt_sample_names(comparison, all_samples, "comparison")
  if (anyNA(comparison))
    stop("comparison can't contain NA: list the samples to draw, e.g. comparison = c(\"",
         all_samples[1], "\", \"", all_samples[min(2, length(all_samples))],
         "\").", call. = FALSE)
  if (anyDuplicated(comparison))
    stop("comparison lists '", comparison[duplicated(comparison)][1],
         "' twice; list each sample once.", call. = FALSE)
  has_of <- function(n) tab$sample[tab$track == n]
  key_of <- function(n, s) tab$key[tab$track == n & !is.na(tab$sample) &
                                     tab$sample == s][1]

  cmp <- unique(nms[kinds == "comparison"])
  if (length(cmp)) {
    ex  <- has_of(cmp[1])[1]
    arg <- paste0("comparison = c(\"", sub(" vs .*", "", ex), "\", \"",
                  sub(".* vs ", "", ex), "\")")
    if (!is.null(samples))
      stop("'", cmp[1],
           "' compares 2 samples, which samples = can't give: drop samples and set ",
           arg, " (reference second).", call. = FALSE)
    if (is.null(comparison))
      stop("'", cmp[1], "' compares 2 samples: set ", arg,
           " (reference second).", call. = FALSE)
    if (length(comparison) != 2)
      stop("'", cmp[1],
           "' compares 2 samples, so comparison needs exactly 2 (reference second); it has ",
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
        stop("samples[", bad[1], "] is '", slot[bad[1]], "', but track ",
             bad[1], " ('", nms[bad[1]],
             "') has no sample: put NA there.", call. = FALSE)
    } else {
      stop("samples has ", length(samples), " entries for ", n,
           if (n == 1) " track: give 1." else
             paste0(" tracks: give 1 (that sample for every track) or ", n,
                    " (one per track, NA for ",
                    "tracks without a sample, such as genes)."), call. = FALSE)
    }
    keys <- character(n)
    for (i in seq_len(n)) {
      if (kinds[i] == "fixed") {
        keys[i] <- tab$key[tab$track == nms[i]][1]; next
      }
      if (is.na(slot[i]))
        stop("Track ", i, " ('", nms[i], "') needs a sample, but samples[",
             if (length(samples) == 1) 1 else i,
             "] is NA.", call. = FALSE)
      keys[i] <- key_of(nms[i], slot[i])
      if (is.na(keys[i]))
        stop("No '", nms[i], "' track for '", slot[i],
             "'. Samples that have it: ",
             paste(has_of(nms[i]), collapse = ", "), ".", call. = FALSE)
    }
    if (all(kinds == "fixed"))
      warning("samples has no effect: none of these tracks needs a sample (",
              paste(unique(nms), collapse = ", "),
              "). It picks the sample of per-sample tracks; see pt_tracks().",
              call. = FALSE)
    return(keys)
  }

  if (!is.null(comparison) && all(kinds == "fixed"))
    warning("comparison has no effect: none of these tracks needs a sample (",
            paste(unique(nms), collapse = ", "),
            "). It picks the samples of per-sample and comparison tracks; see pt_tracks().",
            call. = FALSE)
  over <- comparison %||% all_samples
  one <- function(i, s) {
    k <- key_of(nms[i], s)
    if (is.na(k)) {
      message("Note: no '", nms[i], "' for ", s, ", skipped.");
      return(character(0))
    }
    k
  }
  keys <- character(0)
  i <- 1L
  while (i <= length(nms)) {
    if (kinds[i] == "fixed") {
      keys <- c(keys, tab$key[tab$track == nms[i]][1])
    } else if (kinds[i] == "comparison") {
      pair <- paste(comparison[1], "vs", comparison[2])
      k <- key_of(nms[i], pair)
      if (is.na(k))
        stop("No '", nms[i], "' track for ", pair, ". Pairs that have it: ",
             paste(has_of(nms[i]), collapse = ", "), ".", call. = FALSE)
      keys <- c(keys, k)
    } else {
      # A run of per-sample tracks is one block: by track (each track over every
      # sample) or by sample (each sample's tracks, then the next sample's).
      j <- i
      while (j < length(nms) && kinds[j + 1] == "sample") j <- j + 1L
      block <- i:j
      if (order == "track")
        for (b in block) for (s in over) keys <- c(keys, one(b, s))
      else for (s in over) for (b in block) keys <- c(keys, one(b, s))
      i <- j
    }
    i <- i + 1L
  }
  keys
}

#' Describe a PlotTrellis figure
#'
#' Picks a gene (or region) and the tracks to show, in order, top to bottom.
#' Nothing is drawn until the figure is printed, saved with [pt_save()], or
#' exported with [pt_script()].
#'
#' Track names are the ones [pt_tracks()] lists: `genes`, `transcripts`,
#' `transcripts_all`, `transcripts_annotated` and `protein_domains` need no
#' sample. The others are named after your `samples.txt`'s columns (`bw:atac`
#' gives `atac`; `quantdir:` gives `transcript_expression` (`transcript_counts`,
#' its old name, still works); comparisons get a `_diff` suffix, such as
#' `junctions_diff`).
#'
#' Which samples a per-sample track shows:
#' * neither `samples` nor `comparison`: every sample, in `samples.txt` order;
#' * `samples = "K562"`: that one sample for every track;
#' * `samples` as long as `tracks`: track i shows sample i, with `NA` for tracks
#'   that need no sample, e.g.
#'   `tracks = c("genes", "atac", "atac"), samples = c(NA, "K562", "GM12878")`;
#' * `comparison = c("K562", "GM12878")`: every per-sample track for those
#'   samples, in that order. A comparison track (such as `junctions_diff`) needs
#'   exactly 2 here and compares the first against the second, the reference
#'   ("K562 vs GM12878", as its heading writes it).
#'
#' Tracks and samples can also be given by number: tracks by their `id` in
#' [pt_tracks()], samples by their row in `samples.txt` ([pt_samples()]).
#' Numbers change when `samples.txt` does, so use names in scripts you keep.
#' Sample names must be written exactly as in `samples.txt`.
#'
#' @param proj A [pt_project()].
#' @param gene A gene symbol (any case) or Ensembl gene ID. The view is the gene
#'   plus the project's flank setting on each side.
#' @param tracks Track names (or numbers, see [pt_tracks()]), top to bottom.
#' @param samples Which sample each per-sample track shows: one name (or number)
#'   for all of them, or one per track (`NA` for tracks without a sample).
#'   `NULL` (the default): every sample.
#' @param comparison Instead of `samples`: the samples to repeat every
#'   per-sample track over, in order; with a comparison track, exactly 2, the
#'   reference second.
#' @param region Instead of `gene`: a region such as
#'   `"chr1:198,600,000-198,760,000"`.
#' @param title Figure title.
#' @param order When a track is repeated over samples: `"track"` (all samples of
#'   one track, then the next track) or `"sample"` (all of one sample's tracks,
#'   then the next sample's). Grouping applies to per-sample tracks written next
#'   to each other: tracks without a sample and comparison tracks stay where
#'   they are written and split the groups, so write them first or last.
#' @param ... Not used; an old `parallel =` stops with what to write instead.
#' @return A `pt_figure`.
#' @seealso [pt_save()], [pt_script()], [pt_tracks()], [pt_samples()]
#' @examples
#' proj <- pt_project(file.path(tempdir(), "pt-example-human"),
#'                    pt_example("human"))
#' fig <- pt_figure(proj, "PSME2",
#'                  c("genes", "transcripts", "atac", "junctions"))
#' pt_save(fig, file.path(tempdir(), "PSME2.pdf"))
#'
#' # One sample for every track:
#' pt_figure(proj, "PSME2", c("genes", "atac", "junctions"), samples = "K562")
#'
#' # One sample per track, in any order:
#' pt_figure(proj, "PSME2", c("genes", "atac", "atac"),
#'           samples = c(NA, "K562", "GM12878"))
#'
#' # Two samples and their comparison (K562 against the reference GM12878), grouped by sample:
#' pt_figure(proj, "PSME2", c("atac", "junctions", "junctions_diff", "genes"),
#'           comparison = c("K562", "GM12878"), order = "sample")
#' @export
pt_figure <- function(proj, gene = NULL, tracks = "genes", samples = NULL,
                      comparison = NULL, region = NULL, title = "",
                      order = "track", ...) {
  stopifnot(inherits(proj, "pt_project"))
  dots <- ...names() %||% rep("", ...length())
  if ("parallel" %in% dots)
    stop("parallel is gone: leave samples and comparison unset to draw every per-sample track for ",
         "every sample, or set comparison = c(...) to pick the samples and their order.",
         call. = FALSE)
  if (length(dots))
    stop("pt_figure() ",
         if (any(nzchar(dots)))
           paste0("has no argument ",
                  paste0("'", dots[nzchar(dots)], "'", collapse = ", ")),
         if (any(nzchar(dots)) && !all(nzchar(dots))) ", and ",
         if (!all(nzchar(dots))) "takes 8 arguments; this call gives more",
         ". See ?pt_figure.", call. = FALSE)
  proj <- .pt_renamed_junctions_project(proj)

  if (!is.null(gene) && !(is.character(gene) && length(gene) == 1))
    stop("gene must be one gene name or Ensembl ID. With region =, name the tracks: tracks = c(...).",
         call. = FALSE)
  if (!is.null(gene)) {
    # The TxDb through .pt_annotation(), not proj$txdb: a project kept with
    # saveRDS() (or a knitr cache) comes back with a dead database connection;
    # this reopens it by its path.
    res <- .pt_quiet(.pt_resolve_gene_region(
      gene, .pt_annotation(proj$genome)$txdb,
      proj$render_cfg$gene_flank_pct / 100, proj$genome$orgdb, proj$genes,
      proj$genome$genome_label))
    if (!is.null(res$error)) stop(res$error, call. = FALSE)
    rgn <- res$region
    rgn$gene <- res$symbol
  } else if (!is.null(region)) {
    if (!(is.character(region) && length(region) == 1))
      stop("region must be one string, such as \"chr1:198,600,000-198,760,000\".",
           call. = FALSE)
    m <- regmatches(region,
                    regexec("^\\s*([^:\\s]+):([0-9,]+)-([0-9,]+)\\s*$",
                            region, perl = TRUE))[[1]]
    if (length(m) != 4)
      stop("region must look like \"chr1:198,600,000-198,760,000\".",
           call. = FALSE)
    # The chromosome in the project's own spelling ("2" or "MT" typed, "chr2" or
    # "chrM" in a chr-named project); as typed when the project has no such
    # chromosome.
    chroms <- tryCatch(names(readRDS(proj$genome$chrom_sizes_path)),
                       error = function(e) NULL)
    rgn <- list(chrom = resolve_chrom(m[2], chroms) %||% m[2],
                start = gsub(",", "", m[3]), end = gsub(",", "", m[4]),
                gene = NA_character_)
    err <- .pt_region_error(rgn$start, rgn$end)
    if (!is.null(err)) stop(err, call. = FALSE)
  } else {
    stop(
      "Pass gene = (a symbol or Ensembl ID) or region = (\"chr:start-end\").",
      call. = FALSE
    )
  }

  keys <- .pt_resolve_tracks(proj$tracks, tracks, names(proj$colours),
                             samples, comparison, order)
  structure(list(
    gene          = rgn$gene,
    project_dir   = proj$dir,
    keys          = keys,
    tracks        = proj$tracks$label[match(keys, proj$tracks$key)],
    entries       = .pt_build_entries(keys, proj$items, proj$colours,
                                      proj$render_cfg),
    region_params = list(chrom = rgn$chrom, start = as.integer(rgn$start),
                         end = as.integer(rgn$end), gene = rgn$gene),
    settings      = c(.pt_settings_base(proj$render_cfg), list(title = title))
  ), class = "pt_figure")
}

#' @export
print.pt_figure <- function(x, ...) {
  r <- x$region_params
  cat("PlotTrellis figure: ", if (is.na(x$gene)) "" else paste0(x$gene, ", "),
      r$chrom, ":", format(r$start, big.mark = ","), "-",
      format(r$end, big.mark = ","), "\n", sep = "")
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
#' @param file Output path ending in `.pdf`, `.svg` or `.png`. Its folder must
#'   exist.
#' @param dpi Resolution for `.png`.
#' @return `file`, invisibly.
#' @seealso [pt_figure()], [pt_script()]
#' @export
pt_save <- function(fig, file, dpi = 300) {
  stopifnot(inherits(fig, "pt_figure"))
  fig$entries <- lapply(fig$entries, function(e)
    replace(e, "viz", .pt_renamed_junctions_key(e$viz)))
  dev <- .pt_device(file, dpi, "pt_save()")
  # Closes only the device it opened (its number), then makes the caller's
  # current device current again.
  before <- grDevices::dev.cur()
  opened <- NULL
  on.exit({
    if (!is.null(opened) && opened %in% grDevices::dev.list())
      grDevices::dev.off(opened)
    if (before %in% grDevices::dev.list()) grDevices::dev.set(before)
  })
  .pt_quiet(render_all(fig$entries, fig$region_params, fig$settings,
                       device_fn = function(w, h) {
                         dev(w, h); opened <<- grDevices::dev.cur()
                       }))
  if (!is.null(opened)) grDevices::dev.off(opened)
  opened <- NULL
  if (tolower(tools::file_ext(file)) == "svg") .pt_tidy_svg(file)
  invisible(file)
}

# A function(w, h) opening the device for `file`, picked by its extension.
.pt_device <- function(file, dpi, caller) {
  switch(tolower(tools::file_ext(file)),
    svg = function(w, h) svglite::svglite(file, width = w, height = h),
    pdf = function(w, h) grDevices::cairo_pdf(file, width = w, height = h),
    png = function(w, h) grDevices::png(file, width = w, height = h,
                                        units = "in", res = dpi),
    stop(caller, " writes .pdf, .svg or .png files.", call. = FALSE))
}

#' Stitch figures into one file
#'
#' Puts several figures in one file, in `ncol` columns filled row by row
#' (one column by default: stacked top to bottom), for many genes in one
#' figure. Each figure keeps the size [pt_save()] gives it, at the top left of
#' its cell; a column is as wide as its widest figure, a row as tall as its
#' tallest. The `.svg` and `.pdf` stay vector with live text. A `.png` can be
#' at most 32767 pixels a side (about 109 inches at 300 dpi); for more, write
#' a `.pdf`, change `ncol` (more columns for a page too tall, fewer for one too
#' wide) or lower `dpi`.
#'
#' @param figs A list of [pt_figure()]s, in reading order.
#' @param file Output path ending in `.pdf`, `.svg` or `.png`. Its folder must
#'   exist.
#' @param ncol Number of columns. More than the number of figures gives one
#'   row.
#' @param dpi Resolution for `.png`.
#' @return `file`, invisibly.
#' @seealso [pt_save()], [pt_grob()] to compose a figure with other plots
#' @examples
#' proj <- pt_project(file.path(tempdir(), "pt-example-human"),
#'                    pt_example("human"))
#' figs <- list()
#' for (g in c("PSME2", "IRF9"))
#'   figs[[g]] <- pt_figure(proj, g, c("genes", "atac"), title = g)
#' pt_stitch(figs, file.path(tempdir(), "all_genes.pdf"))
#' # Side by side:
#' pt_stitch(figs, file.path(tempdir(), "all_genes_2col.pdf"), ncol = 2)
#' @export
pt_stitch <- function(figs, file, ncol = 1, dpi = 300) {
  if (inherits(figs, "pt_figure")) figs <- list(figs)
  bad <- which(!vapply(figs, inherits, logical(1), "pt_figure"))
  if (!length(figs) || length(bad))
    stop("figs must be a list of pt_figure()s",
         if (length(bad)) paste0("; item ", bad[1], " is not one"), ".",
         call. = FALSE)
  if (!is.numeric(ncol) || length(ncol) != 1 || is.na(ncol) || ncol < 1 ||
      ncol != round(ncol))
    stop("ncol must be a whole number, 1 or more.", call. = FALSE)
  ext <- tolower(tools::file_ext(file))
  dev <- .pt_device(file, dpi, "pt_stitch()")
  if (!dir.exists(dirname(file)))
    stop("The folder ", dirname(file), " does not exist.", call. = FALSE)
  grobs <- lapply(seq_along(figs), function(i)
    withCallingHandlers(.pt_grob(figs[[i]], ext, dpi), error = function(e) {
      gene <- figs[[i]]$gene
      stop("Figure ", i, if (isTRUE(!is.na(gene))) paste0(" (", gene, ")"), ": ",
           conditionMessage(e), call. = FALSE)
    }))
  sizes <- vapply(grobs, function(g) c(g$width, g$height), numeric(2))
  # Cells filled row by row; a column as wide as its widest figure, a row as
  # tall as its tallest.
  ncol <- min(ncol, length(grobs))
  col <- (seq_along(grobs) - 1) %% ncol + 1
  row <- (seq_along(grobs) - 1) %/% ncol + 1
  col_w <- vapply(seq_len(ncol), function(j) max(sizes[1, col == j]), 0)
  row_h <- vapply(seq_len(max(row)), function(r) max(sizes[2, row == r]), 0)
  left <- c(0, cumsum(col_w))[col]
  top <- sum(row_h) - c(0, cumsum(row_h))[row]
  if (ext == "png" && any(ceiling(c(sum(col_w), sum(row_h)) * dpi) > 32767))
    stop("Too big for one PNG at ", dpi, " dpi (at most 32767 pixels a side): ",
         "write a .pdf or .svg, change ncol, lower dpi, or stitch fewer figures.",
         call. = FALSE)
  before <- grDevices::dev.cur()
  opened <- NULL
  on.exit({
    if (!is.null(opened) && opened %in% grDevices::dev.list())
      grDevices::dev.off(opened)
    if (before %in% grDevices::dev.list()) grDevices::dev.set(before)
  })
  dev(sum(col_w), sum(row_h))
  opened <- grDevices::dev.cur()
  grid::grid.newpage()
  for (i in seq_along(grobs)) {
    grid::pushViewport(grid::viewport(
      x = left[i], y = top[i], width = sizes[1, i], height = sizes[2, i],
      just = c("left", "top"), default.units = "in"))
    grid::grid.draw(grobs[[i]])
    grid::popViewport()
  }
  grDevices::dev.off(opened)
  opened <- NULL
  if (ext == "svg") .pt_tidy_svg(file)
  invisible(file)
}

#' A figure as a grid grob, to compose with other plots
#'
#' Draws a figure once and keeps the drawing as a grid grob, for patchwork,
#' cowplot, gridExtra or a grid layout, e.g. beside a ggplot. The grob draws at
#' the figure's own size in inches (its `width` and `height`), centred in the
#' space it is given, and is never scaled: give it that much room, or it spills
#' over its neighbours. In patchwork, give its column its size in inches,
#' `plot_layout(widths = grid::unit(c(g$width, 1), c("in", "null")))`, and
#' `ggsave()` a page somewhat bigger than the figures (margins, tags).
#' [grid::grid.draw()] also draws a [pt_figure()] directly.
#'
#' R's own `pdf()` and `postscript()` devices (a plain `ggsave("both.pdf")`
#' uses `pdf()`) know only a few fonts, not the default DejaVu Sans: on them
#' the figure's text is drawn in the device's own font (usually Helvetica),
#' with a message. For the figure's own font, save with a cairo device, e.g.
#' `ggsave("both.pdf", device = cairo_pdf)`.
#'
#' @param fig A [pt_figure()].
#' @return A grob of class `pt_grob`, with the figure's `width` and `height` in
#'   inches.
#' @seealso [pt_stitch()] to put several figures in one file
#' @examples
#' proj <- pt_project(file.path(tempdir(), "pt-example-human"),
#'                    pt_example("human"))
#' g <- pt_grob(pt_figure(proj, "PSME2", c("genes", "atac"), title = "PSME2"))
#' c(g$width, g$height)
#' if (requireNamespace("patchwork", quietly = TRUE)) {
#'   p <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) + ggplot2::geom_point()
#'   both <- patchwork::wrap_elements(full = g) + p +
#'     patchwork::plot_layout(widths = grid::unit(c(g$width, 1), c("in", "null")))
#'   ggplot2::ggsave(file.path(tempdir(), "both.pdf"), both,
#'                   width = g$width + 4.5, height = g$height + 0.5)
#' }
#' @export
pt_grob <- function(fig) {
  if (!inherits(fig, "pt_figure"))
    stop("pt_grob() takes a pt_figure().", call. = FALSE)
  withCallingHandlers(.pt_grob(fig, "pdf", 300), error = function(e) {
    stop("Could not draw the figure", if (isTRUE(!is.na(fig$gene)))
      paste0(" of ", fig$gene), ": ", conditionMessage(e), call. = FALSE)
  })
}

# The figure drawn on a throwaway device of `ext`'s kind, as a pt_grob.
.pt_grob <- function(fig, ext, dpi) {
  rec <- .pt_capture(fig, ext, dpi)
  size <- attr(rec, "pt_size")
  grid::grob(steps = rec$children, width = size[1], height = size[2],
             gp_page = attr(rec, "pt_gp"),
             font = .pt_device_font(.pt_fill_settings(fig$settings)$font),
             cl = "pt_grob")
}

#' @exportS3Method grid::drawDetails
drawDetails.pt_grob <- function(x, recording) {
  # Named after the grob (grid numbers every grob), so it differs from any
  # other pt_grob's while both are drawn.
  vp <- paste0("pt_grob_", x$name)
  # R's pdf() and postscript() devices draw only the fonts in their own lists
  # (not the default DejaVu Sans; they would warn and draw nothing): the
  # figure's other fonts are drawn in the device's own font instead (Helvetica
  # unless its family is set). postscript() draws only its own family (and its
  # `fonts`), whatever postscriptFonts() lists, so on it every font is swapped.
  # Text keeps the places measured for the figure's own font.
  dev <- names(grDevices::dev.cur())
  known <- switch(dev, pdf = names(grDevices::pdfFonts()),
                  postscript = names(grDevices::postscriptFonts()))
  missing <- !is.null(known) && length(x$font) && nzchar(x$font) &&
    !x$font %in% known
  if (missing)
    message("R's ", dev, "() device has no ", x$font, ": the figure is drawn in ",
            "the device's font (usually Helvetica). For ", x$font, ", save with ",
            "a cairo device, e.g. ggsave(..., device = cairo_pdf).")
  if (missing || identical(dev, "postscript"))
    x$steps <- .pt_known_fonts(x$steps, if (dev == "pdf") known else character(0))
  # The device's defaults the figure was drawn with: settings of the viewport
  # around it (a cex, a bold face) would otherwise change its text. cex and lex
  # multiply down the viewports, so they are divided by the current ones; alpha
  # also multiplies but cannot go above 1, so a parent's alpha still applies.
  gp <- as.list(x$gp_page)
  gp$alpha <- NULL
  for (k in intersect(c("cex", "lex"), names(gp)))
    gp[[k]] <- gp[[k]] / grid::get.gpar(k)[[k]]
  grid::pushViewport(grid::viewport(width = x$width, height = x$height,
                                    default.units = "in", name = vp,
                                    gp = do.call(grid::gpar, gp)))
  # plotgardener finds its page by name from the root (seekViewport("page")):
  # removing this figure's viewports keeps the next figure out of its page.
  # Also after an error, so nothing drawn later lands inside the figure.
  on.exit({
    grid::seekViewport(vp)
    grid::popViewport()
  })
  # Replays the recorded steps one after another, as drawn: through
  # grid.draw() each step would reset the graphical parameters after it, and
  # the genes track's labels would lose the font of the viewport around them.
  for (step in x$steps) eval(step$expr, step$list, getNamespace("grid"))
}

# `x` (recorded drawing steps, grobs, viewports, gpars) with every font family
# not in `known` replaced by "", the device's own.
.pt_known_fonts <- function(x, known) {
  if (inherits(x, "gpar")) {
    fam <- x$fontfamily
    if (!is.null(fam)) x$fontfamily[!fam %in% known] <- ""
  } else if (is.list(x)) {
    for (i in seq_along(x))
      if (!is.null(x[[i]])) x[[i]] <- .pt_known_fonts(x[[i]], known)
  }
  x
}

#' @exportS3Method grid::widthDetails
widthDetails.pt_grob <- function(x) grid::unit(x$width, "in")

#' @exportS3Method grid::heightDetails
heightDetails.pt_grob <- function(x) grid::unit(x$height, "in")

# `g + p` with a ggplot `p` would otherwise go to ggplot2's method, which
# returns NULL without a word (and ggsave() of it writes a blank page).
# `p + g` stays ggplot2's (patchwork adds a grob); before ggplot2 4.0 it comes
# here.
#' @export
`+.pt_grob` <- function(e1, e2) {
  if (inherits(e1, "gg") && !inherits(e1, "pt_grob"))
    return(ggplot2::ggplot_add(e2, e1, deparse(substitute(e2))))
  stop("A pt_grob can't be added with +: use ",
       "patchwork::wrap_elements(full = g) + p (see ?pt_grob).", call. = FALSE)
}

# Picks `+.pt_grob` over the other side's method (ggplot2's S7 one included).
#' @export
chooseOpsMethod.pt_grob <- function(x, y, mx, my, cl, reverse) TRUE

#' @exportS3Method grid::grid.draw
grid.draw.pt_figure <- function(x, recording = TRUE) {
  grid::grid.draw(pt_grob(x), recording = recording)
}

# Draws a figure once on a throwaway device of the output's kind (so text is
# measured as in the final file) and keeps the recorded drawing steps, with the
# figure's size in inches as attribute "pt_size".
.pt_capture <- function(fig, ext, dpi) {
  fig$entries <- lapply(fig$entries, function(e)
    replace(e, "viz", .pt_renamed_junctions_key(e$viz)))
  tmp <- tempfile(fileext = paste0(".", ext))
  dev <- .pt_device(tmp, dpi, "pt_stitch()")
  before <- grDevices::dev.cur()
  opened <- NULL
  size <- NULL
  gp <- NULL
  on.exit({
    if (!is.null(opened) && opened %in% grDevices::dev.list())
      grDevices::dev.off(opened)
    if (before %in% grDevices::dev.list()) grDevices::dev.set(before)
    unlink(tmp)
  })
  .pt_quiet(render_all(fig$entries, fig$region_params, fig$settings,
                       device_fn = function(w, h) {
                         dev(w, h); opened <<- grDevices::dev.cur()
                         size <<- c(w, h)
                         gp <<- unclass(grid::get.gpar())
                       }))
  structure(grid::grid.grab(wrap.grobs = TRUE), pt_size = size, pt_gp = gp)
}

# Hides plotgardener's per-grob "segments[segments2]" lines, hundreds per
# figure, for the API and the app; other messages pass. (AnnotationDbi's
# "select() returned 1:many" notes are silenced at each select() call; the
# pattern stays as a backstop.) Also hides plotSignal()'s bin warnings for a
# short view (about 2 kb or less), which it adjusts for itself.
.pt_quiet <- function(expr) withCallingHandlers(expr, message = function(m) {
  if (grepl("^\\w+\\[\\w+\\]\\s*$|^'select\\(\\)' returned",
            conditionMessage(m)))
    invokeRestart("muffleMessage")
}, warning = function(w) {
  if (grepl("^Too many bins: adjusting|^Number of bins larger than plot length",
            conditionMessage(w)))
    invokeRestart("muffleWarning")
})

#' Export a figure as a standalone R script
#'
#' Writes the same self-contained script the app's Download R script button
#' gives: it redraws the figure with plotgardener alone (PlotTrellis need not be
#' installed), and is the place to start for edits the API doesn't offer.
#'
#' @param fig A [pt_figure()].
#' @param file Path for the `.R` script.
#' @return `file`, invisibly.
#' @export
pt_script <- function(fig, file) {
  stopifnot(inherits(fig, "pt_figure"))
  fig$entries <- lapply(fig$entries, function(e)
    replace(e, "viz", .pt_renamed_junctions_key(e$viz)))
  writeLines(assemble_download_script(fig$entries, fig$region_params,
                                      fig$settings,
                                      hic_tops = .pt_hic_figure_tops(fig)),
             file)
  invisible(file)
}

# The tops of the Hi-C colour scales `fig` draws (.pt_hic_colour_tops()), for
# its script: its Hi-C tracks drawn alone to a temporary file. NULL when it has
# none or they could not be drawn.
.pt_hic_figure_tops <- function(fig) {
  hic <- Filter(function(e) identical(e$viz, "hic"), fig$entries)
  if (!length(hic)) return(NULL)
  settings <- fig$settings
  settings$output_path <- tempfile(fileext = ".svg")
  on.exit(unlink(settings$output_path))
  res <- .pt_render_svg(hic, fig$region_params, settings)
  if (is.list(res)) res$hic_tops
}

#' @exportS3Method knitr::knit_print
knit_print.pt_figure <- function(x, ...) {
  # Outside a knit, fig_path() would be figure/-1.png in the working directory:
  # a temp file then.
  f <- if (isTRUE(getOption("knitr.in.progress"))) knitr::fig_path(".png")
       else tempfile(fileext = ".png")
  dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
  knitr::knit_print(knitr::include_graphics(pt_save(x, f, dpi = 150)), ...)
}

#' Example samples.txt
#'
#' The path of a `samples.txt` shipped with PlotTrellis, for [run_plottrellis()]
#' or [pt_project()]. Its data is read where it is installed (a project keeps
#' its caches in its own folder), and the app's gene box suggests its gene.
#'
#' - `"mammary"` (the default): mouse mammary basal cells (Basal) vs luminal
#'   progenitors (LP) (Milevskiy et al. 2023, GEO GSE227750) at five genes
#'   (Eya2, Lrrfip1, Acot7, Dzip1, Cxcl14), mm39, with every data type:
#'   RNA-seq coverage, splice junctions and transcript expression, ATAC and
#'   H3K36me3 (with its IgG control) coverage and peaks, and Hi-C (Omni-C).
#'   Basal cells mostly use Eya2's isoform from an upstream promoter, LP a
#'   shorter one. Needs `org.Mm.eg.db`.
#' - `"human"`: ENCODE GM12878 vs K562 at PSME1/PSME2, hg38: ATAC, splice
#'   junctions and transcript expression. Needs `org.Hs.eg.db`.
#'
#' @param dataset `"mammary"` or `"human"`.
#' @return The path of the example's `samples.txt`.
#' @examples
#' pt_example()
#' \dontrun{
#' run_plottrellis(project_dir = "PlotTrellis-example",
#'                 samples_txt = pt_example())
#' }
#' @export
pt_example <- function(dataset = c("mammary", "human")) {
  dataset <- match.arg(dataset)
  system.file("extdata",
              if (dataset == "mammary") "example" else "example-human",
              "samples.txt", package = "PlotTrellis", mustWork = TRUE)
}

# The junction track's viz key used to be "sashimi". A saved pt_figure can
# still hold its item keys ("sashimi_<name>_<sample>") and its entries' viz
# ("sashimi"); this gives their "junctions" equivalents. NULL (no tracks) stays
# NULL.
.pt_renamed_junctions_key <- function(keys) {
  if (is.null(keys)) return(NULL)
  sub("^sashimi(_|$)", "junctions\\1", keys)
}

# A pt_project made before that rename and kept (saveRDS()), with its track
# items and track table keyed as they are made now. A project made since comes
# back as it is.
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
