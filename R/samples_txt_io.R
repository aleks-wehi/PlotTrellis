# samples.txt import/export: one wide tab-separated table, one row per
# condition, with `type:name` column headers, plus a `# data_root:` line and
# `# <key>: value` lines for other settings. Pure functions, no Shiny: used by
# the wizard and parse_samples_txt().

# A cell read from samples.txt: NA or length 0 -> "". (%||%, from boilerplate.R,
# would leave character(0) as it is.)
chr0 <- function(x) if (length(x) == 0 || is.na(x)) "" else x

# Must resolve to an actual FILE, not a directory (file.exists() is TRUE for
# both) -- a directory here would crash is_paired_end_structural() or silently
# pass only to fail at render.
is_real_file <- function(path) file.exists(path) && !dir.exists(path)

DATA_TYPE_VOCAB <- c("bw", "bed", "hic", "quantdir", "paired-end")

# Settings keys that never round-trip through samples.txt's `# <key>: value`
# delta-header lines: `data_dir` already has its own dedicated `# data_root:`
# line; `samples` is structural, not a render-style setting. The genome fields
# (gtf_path/orgdb/genome_label) DO round-trip through the generic delta
# mechanism below (see SAMPLES_TXT_ALWAYS_EXPORT below for how export treats
# them differently). Everything else in DEFAULT_SETTINGS is fair game.
SAMPLES_TXT_DELTA_EXCLUDE <- c("data_dir", "samples")

# Genome identity fields ALWAYS get a `# key: value` line on export,
# regardless of whether the live value happens to equal DEFAULT_SETTINGS':
# delta-only export treats "matches the default" as "nothing to say", which is
# right for a cosmetic setting like page_w but wrong for identifying,
# must-travel-with-the-data information (a project on the default genome
# would otherwise export no genome at all). Still parsed on import via the
# same generic `# <key>: value` mechanism as every other key -- only the
# export side special-cases these.
SAMPLES_TXT_ALWAYS_EXPORT <- c("gtf_path", "orgdb", "genome_label")

# Coerces `val` to match `default_val`'s type (e.g. a logical default like
# uniprot_reviewed_only must come back logical, not "TRUE"/"FALSE"); malformed
# input coerces to NA, which the caller drops rather than injecting into live
# settings.
settings_value_coerce <- function(val, default_val) {
  if (is.logical(default_val)) return(suppressWarnings(as.logical(val)))
  if (is.numeric(default_val)) return(suppressWarnings(as.numeric(val)))
  val
}

# Whether `a` differs from `b` for export's delta-only decision: numeric values
# compare by value (all.equal, both coerced to double), since identical(8, 8L)
# is FALSE and would spuriously emit a delta line on every numericInput/jsonlite
# round trip; everything else compares by identical().
settings_value_differs <- function(a, b) {
  if (is.numeric(a) && is.numeric(b))
    return(!isTRUE(all.equal(as.numeric(a), as.numeric(b))))
  !identical(a, b)
}

# Cheap structural check only (a HiCCUPS header line, or the column count of the
# first data line, past any header lines, matching BED12's 12 or BED6/BEDPE/STAR
# SJ.out.tab's 6-10 and not a BED peaks file; not full semantic validation, see
# read_junction_file() in transcript_filter.R); caller must check is_real_file()
# first, this doesn't guard against `path` being a directory.
is_paired_end_structural <- function(path) {
  if (is_hiccups_header(readLines(path, n = 1, warn = FALSE))) return(TRUE)
  fields <- .pt_arc_first_data(path)$fields
  n_cols <- length(fields)
  n_cols == 12 || (n_cols >= 6 && n_cols <= 10 && !.pt_is_bed_peaks(fields))
}

# Parses samples.txt into list(errors=, warnings=, data_root=, overrides=,
# samples=) (samples[[i]]$data has the shape of a project's
# settings$samples[[i]]$data); errors are fatal and callers must check them
# before trusting samples, warnings are non-fatal broken paths, checked against
# the file's own `# data_root:` line.
import_samples_txt <- function(path, defaults = list()) {
  lines <- readLines(path, warn = FALSE)
  errors <- character(0)
  warnings <- character(0)

  data_root <- NULL
  overrides <- list()
  body_lines <- character(0)
  for (ln in lines) {
    if (grepl("^#\\s*data_root:", ln)) {
      data_root <- trimws(sub("^#\\s*data_root:", "", ln))
    } else if (grepl("^#\\s*([a-zA-Z_][a-zA-Z0-9_]*):(.*)$", ln)) {
      key <- trimws(sub("^#\\s*([a-zA-Z_][a-zA-Z0-9_]*):(.*)$", "\\1", ln))
      val <- trimws(sub("^#\\s*([a-zA-Z_][a-zA-Z0-9_]*):(.*)$", "\\2", ln))
      if (key %in% names(defaults) && !(key %in% SAMPLES_TXT_DELTA_EXCLUDE)) {
        coerced <- settings_value_coerce(val, defaults[[key]])
        if (!(length(coerced) == 1 && is.na(coerced)))
          overrides[[key]] <- coerced
      }
    } else if (grepl("^#", ln)) {
      next
    } else if (nchar(trimws(ln)) > 0) {
      body_lines <- c(body_lines, ln)
    }
  }

  # A relative `# data_root:`/`# gtf_path:` is relative to this file's own
  # folder, not R's working directory, so a samples.txt can travel with its
  # data. Absolute paths (or ~) and blanks are left exactly as written.
  base_dir <- normalizePath(dirname(path), winslash = "/")
  is_absolute <- function(p) grepl("^(/|~|[A-Za-z]:[/\\\\]|\\\\\\\\)", p)
  resolve <- function(p) {
    if (is.null(p) || !nzchar(p) || is_absolute(p)) return(p)
    normalizePath(file.path(base_dir, p), winslash = "/", mustWork = FALSE)
  }
  data_root <- resolve(data_root)
  if (!is.null(overrides$gtf_path))
    overrides$gtf_path <- resolve(overrides$gtf_path)

  if (length(body_lines) == 0) {
    return(list(errors = c(errors, "No header/data rows found"),
                warnings = warnings,
                data_root = data_root, overrides = overrides, samples = list()))
  }

  header <- strsplit(body_lines[1], "\t", fixed = TRUE)[[1]]
  fixed_cols <- c("condition", "colour")
  col_spec <- vector("list", length(header))
  for (j in seq_along(header)) {
    h <- header[j]
    if (h %in% fixed_cols) {
      col_spec[[j]] <- list(kind = "fixed", name = h)
    } else if (grepl(":", h, fixed = TRUE)) {
      parts <- strsplit(h, ":", fixed = TRUE)[[1]]
      type <- parts[1]
      name <- paste(parts[-1], collapse = ":")
      # Alias: older samples.txt files use "junction:<name>" for paired-end
      # data -- normalize so they keep loading without a manual edit.
      if (identical(type, "junction")) type <- "paired-end"
      if (!(type %in% DATA_TYPE_VOCAB)) {
        errors <- c(errors, sprintf(
          "Unrecognized type prefix '%s' in column '%s' (known types: %s)",
          type, h, paste(DATA_TYPE_VOCAB, collapse = ", ")))
        col_spec[[j]] <- list(kind = "bad")
      } else {
        col_spec[[j]] <- list(kind = "instance", type = type, name = name)
      }
    } else {
      errors <- c(errors, sprintf(
        "Malformed header column '%s' (not a known metadata column, and no 'type:name' prefix)",
        h))
      col_spec[[j]] <- list(kind = "bad")
    }
  }

  seen_conditions <- character(0)
  samples <- list()
  # per imported row: its condition and the instance columns it fills
  filled <- list()
  if (length(body_lines) > 1) {
    for (i in seq(2, length(body_lines))) {
      fields <- strsplit(body_lines[i], "\t", fixed = TRUE)[[1]]
      length(fields) <- length(header)  # pad short rows with NA

      cond <- trimws(chr0(fields[match("condition", header)]))
      if (nchar(cond) == 0) {
        errors <- c(errors, sprintf("Row %d: missing condition", i))
        next
      }
      if (cond %in% seen_conditions) {
        errors <- c(errors, sprintf("Duplicate condition '%s' (row %d)",
                                    cond, i))
        next
      }
      seen_conditions <- c(seen_conditions, cond)

      data <- list()
      for (j in seq_along(col_spec)) {
        spec <- col_spec[[j]]
        if (is.null(spec) || spec$kind != "instance") next
        val <- trimws(fields[j])
        if (is.na(val) || nchar(val) == 0) next
        if (!is.null(data_root)) {
          resolved <- file.path(data_root, val)
          if (spec$type == "quantdir") {
            resolved <- file.path(resolved, "quant.sf")
            ok <- is_real_file(resolved)
          } else {
            ok <- is_real_file(resolved)
          }
          # An absolute path is joined onto the data root too (here and at
          # render), so the file it names may exist and still not be found: say
          # so, not just "does not resolve".
          if (!ok) {
            warnings <- c(warnings, if (is_absolute(val)) sprintf(
              "Row %d (%s): path for '%s:%s' is an absolute path; write it relative to the data root (%s): %s",
              i, cond, spec$type, spec$name, data_root, val) else sprintf(
              "Row %d (%s): path for '%s:%s' does not resolve to a file: %s",
              i, cond, spec$type, spec$name, val))
          }
        }
        data[[length(data) + 1]] <- list(
          type = spec$type, name = spec$name, path = val
        )
        filled[[cond]] <- c(filled[[cond]], j)
      }

      col_val <- trimws(chr0(fields[match("colour", header)]))

      samples[[length(samples) + 1]] <- list(
        condition = cond,
        colour    = if (nchar(col_val) > 0) col_val else NA_character_,
        data      = data
      )
    }
  }

  # Never return partial samples alongside a structural error -- enforced at
  # the function boundary itself, not left to callers to remember to check
  # `errors` first. `warnings` (unresolved
  # paths) do NOT trigger this: those rows are still genuinely usable, just
  # missing one file, so every OTHER sample still comes in.
  if (length(errors) > 0) samples <- list()

  list(errors = errors, warnings = warnings, data_root = data_root,
       overrides = overrides, samples = samples,
       layout_warning = if (length(errors) == 0)
         .pt_layout_warning(header, col_spec, filled) else character(0))
}

# One column per sample instead of one per data type (`bw:Grans ATAC`,
# `bw:pDC ATAC`, each filled on one row): tracks are named after their column,
# and only tracks of the same name share a scale and get comparison tracks, so
# such a file imports but no two samples can be compared. Flagged when two or
# more columns of one type are each filled on exactly one row, on different
# rows, and have the same name once that row's sample name is taken out
# ("Grans ATAC", "pDC ATAC" -> "ATAC"). Different data on different samples
# (`bw:CUT&RUN` only on A, `bw:H3K27ac` only on B) is a normal file and isn't
# flagged. Not caught: replicate rows (Grans_r1, Grans_r2 filling
# `bw:Grans ATAC`), since those columns fill more than one row. Returns one
# message (character(1)) or character(0).
.pt_layout_warning <- function(header, col_spec, filled) {
  if (length(filled) < 2) return(character(0))
  row_of <- list()  # column index -> conditions filling it
  for (cond in names(filled)) for (j in filled[[cond]])
    row_of[[as.character(j)]] <- c(row_of[[as.character(j)]], cond)
  single <- Filter(function(j) length(row_of[[as.character(j)]]) == 1,
                   as.integer(names(row_of)))
  # The column's name without its sample's name, or NA if the name doesn't
  # contain it.
  base <- vapply(single, function(j) {
    name <- col_spec[[j]]$name
    cond <- row_of[[as.character(j)]]
    pattern <- paste0("(^|[ _.:-])",
                      gsub("([][{}()+*^$|\\\\.?])", "\\\\\\1", cond),
                      "($|[ _.:-])")
    stripped <- trimws(gsub("^[ _.:-]+|[ _.:-]+$", "",
                            sub(pattern, " ", name, ignore.case = TRUE)))
    if (nzchar(stripped) && !identical(stripped, name))
      stripped else NA_character_
  }, "")
  keep <- !is.na(base)
  if (!any(keep)) return(character(0))
  single <- single[keep]
  type_base <- paste0(vapply(col_spec[single], `[[`, "", "type"), ":",
                      base[keep])
  # Grouped ignoring case ("grans_atac", "PDC-ATAC"); the suggestion keeps the
  # first column's case.
  groups <- Filter(
    function(js) length(unique(unlist(row_of[as.character(js)]))) >= 2,
    split(single, tolower(type_base))
  )
  if (length(groups) == 0) return(character(0))
  flagged <- sort(unlist(groups, use.names = FALSE))
  firsts <- sort(vapply(groups, min, 0L))
  suggest <- type_base[match(firsts, single)]
  shown <- header[flagged]
  more <- if (length(shown) > 4)
    sprintf(" (and %d more)", length(shown) - 4) else ""
  sprintf(paste0(
    "Columns %s%s are each used by one sample only, so every sample gets its own track and ",
    "its own scale, and there are no comparison tracks. Use one column per data type, with a ",
    "row per sample: %s. See vignette(\"PlotTrellis\"), section \"Your own data in code ",
    "(writing a samples.txt)\"."),
    paste(head(shown, 4), collapse = ", "), more,
    paste(suggest, collapse = ", "))
}

# Mirrors import_samples_txt(): flattens settings$samples[] into a wide
# type:name-header samples.txt plus a `# data_root:` line and one
# `# <key>: value` delta line per non-default setting (only those, to keep the
# file short), except the 3 genome fields which always get a line
# (SAMPLES_TXT_ALWAYS_EXPORT); pass DEFAULT_SETTINGS as `defaults`.
export_samples_txt <- function(settings, defaults = list()) {
  samples <- settings$samples %||% list()
  delta_keys <- setdiff(intersect(names(settings), names(defaults)),
                        SAMPLES_TXT_DELTA_EXCLUDE)
  delta_lines <- character(0)
  for (key in delta_keys) {
    val <- settings[[key]]
    # A NULL/NA live value (e.g. a genome cache path not yet built) must not
    # emit "# key: NA" -- unlike a numeric NULL (already handled by
    # settings_value_differs()'s all.equal path), a character NA/NULL would
    # otherwise print the literal text "NA", which re-imports as the STRING
    # "NA" (settings_value_coerce() returns a character value as-is), not as
    # "still unset". Guarded against here for the character case
    # specifically.
    if (is.null(val) || (length(val) == 1 && is.na(val))) next
    if (key %in% SAMPLES_TXT_ALWAYS_EXPORT ||
        settings_value_differs(val, defaults[[key]])) {
      delta_lines <- c(delta_lines, paste0("# ", key, ": ", val))
    }
  }
  all_cols <- unique(unlist(lapply(samples, function(s) {
    vapply(s$data, function(d) paste0(d$type, ":", d$name), character(1))
  })))
  header <- c("condition", "colour", all_cols)

  rows <- vapply(samples, function(s) {
    lookup <- character(0)
    if (length(s$data) > 0) {
      lookup <- setNames(vapply(s$data, function(d) d$path, character(1)),
                         vapply(s$data, function(d) paste0(d$type, ":", d$name),
                                character(1)))
    }
    vals <- vapply(header, function(h) {
      switch(h,
        condition = s$condition %||% "",
        colour    = chr0(s$colour %||% NA_character_),
        chr0(unname(lookup[h]))
      )
    }, character(1))
    paste(vals, collapse = "\t")
  }, character(1))

  c(
    paste0("# data_root: ", settings$data_root %||% ""),
    delta_lines,
    paste(header, collapse = "\t"),
    rows
  )
}
