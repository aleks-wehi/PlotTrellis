# samples.txt import/export: one wide tab-separated table, one row per condition, with
# `type:name` column headers, plus a `# data_root:` line and `# <key>: value` lines for
# other settings. Pure functions, no Shiny: used by the wizard and parse_samples_txt().

# A cell read from samples.txt: NA or length 0 -> "". (%||%, from boilerplate.R, would
# leave character(0) as it is.)
chr0 <- function(x) if (length(x) == 0 || is.na(x)) "" else x

# Must resolve to an actual FILE, not a directory (file.exists() is TRUE for both) -- a directory here would crash is_paired_end_structural() or silently pass only to fail at render.
is_real_file <- function(path) file.exists(path) && !dir.exists(path)

DATA_TYPE_VOCAB <- c("bw", "bed", "hic", "quantdir", "paired-end")

# Settings keys that never round-trip through samples.txt's `# <key>: value`
# delta-header lines: `data_dir` already has its own dedicated
# `# data_root:` line; `samples` is structural, not a render-style setting.
# The genome fields (gtf_path/orgdb/genome_label) DO round-trip through the
# generic delta mechanism below (see SAMPLES_TXT_ALWAYS_EXPORT below for how
# export treats them differently). Everything else in DEFAULT_SETTINGS is fair game.
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

# Coerces `val` to match `default_val`'s type (e.g. a logical default like uniprot_reviewed_only must come back logical, not "TRUE"/"FALSE"); malformed input coerces to NA, which the caller drops rather than injecting into live settings.
settings_value_coerce <- function(val, default_val) {
  if (is.logical(default_val)) return(suppressWarnings(as.logical(val)))
  if (is.numeric(default_val)) return(suppressWarnings(as.numeric(val)))
  val
}

# Whether `a` differs from `b` for export's delta-only decision: numeric values compare by value (all.equal, both coerced to double), since identical(8, 8L) is FALSE and would spuriously emit a delta line on every numericInput/jsonlite round trip; everything else compares by identical().
settings_value_differs <- function(a, b) {
  if (is.numeric(a) && is.numeric(b)) return(!isTRUE(all.equal(as.numeric(a), as.numeric(b))))
  !identical(a, b)
}

# Cheap structural check only (a HiCCUPS header line, or a column count matching BED12's 12 or BED6/BEDPE's 6-10, not full semantic validation, see read_junction_file() in transcript_filter.R); caller must check is_real_file() first, this doesn't guard against `path` being a directory.
is_paired_end_structural <- function(path) {
  if (is_hiccups_header(readLines(path, n = 1, warn = FALSE))) return(TRUE)
  con <- file(path, "r")
  on.exit(close(con))
  repeat {
    ln <- readLines(con, n = 1, warn = FALSE)
    if (length(ln) == 0) return(FALSE)
    if (grepl("^#", ln) || nchar(trimws(ln)) == 0) next
    n_cols <- length(strsplit(ln, "\t", fixed = TRUE)[[1]])
    return(n_cols == 12 || (n_cols >= 6 && n_cols <= 10))
  }
}

# Parses samples.txt into list(errors=, warnings=, data_root=, overrides=, samples=) (samples[[i]]$data has the shape of a project's settings$samples[[i]]$data); errors are fatal and callers must check them before trusting samples, warnings are non-fatal broken paths, checked against the file's own `# data_root:` line.
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
        if (!(length(coerced) == 1 && is.na(coerced))) overrides[[key]] <- coerced
      }
    } else if (grepl("^#", ln)) {
      next
    } else if (nchar(trimws(ln)) > 0) {
      body_lines <- c(body_lines, ln)
    }
  }

  # A relative `# data_root:`/`# gtf_path:` is relative to this file's own folder, not R's
  # working directory, so a samples.txt can travel with its data.
  # Absolute paths (or ~) and blanks are left exactly as written.
  base_dir <- normalizePath(dirname(path), winslash = "/")
  resolve <- function(p) {
    if (is.null(p) || !nzchar(p) || grepl("^(/|~|[A-Za-z]:[/\\\\]|\\\\\\\\)", p)) return(p)
    normalizePath(file.path(base_dir, p), winslash = "/", mustWork = FALSE)
  }
  data_root <- resolve(data_root)
  if (!is.null(overrides$gtf_path)) overrides$gtf_path <- resolve(overrides$gtf_path)

  if (length(body_lines) == 0) {
    return(list(errors = c(errors, "No header/data rows found"), warnings = warnings,
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
        "Malformed header column '%s' (not a known metadata column, and no 'type:name' prefix)", h))
      col_spec[[j]] <- list(kind = "bad")
    }
  }

  seen_conditions <- character(0)
  samples <- list()
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
        errors <- c(errors, sprintf("Duplicate condition '%s' (row %d)", cond, i))
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
          if (!ok) {
            warnings <- c(warnings, sprintf(
              "Row %d (%s): path for '%s:%s' does not resolve to a file: %s",
              i, cond, spec$type, spec$name, val))
          }
        }
        data[[length(data) + 1]] <- list(type = spec$type, name = spec$name, path = val)
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
       overrides = overrides, samples = samples)
}

# Mirrors import_samples_txt(): flattens settings$samples[] into a wide type:name-header samples.txt plus a `# data_root:` line and one `# <key>: value` delta line per non-default setting (only those, to keep the file short), except the 3 genome fields which always get a line (SAMPLES_TXT_ALWAYS_EXPORT); pass DEFAULT_SETTINGS as `defaults`.
export_samples_txt <- function(settings, defaults = list()) {
  samples <- settings$samples %||% list()
  delta_keys <- setdiff(intersect(names(settings), names(defaults)), SAMPLES_TXT_DELTA_EXCLUDE)
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
    if (key %in% SAMPLES_TXT_ALWAYS_EXPORT || settings_value_differs(val, defaults[[key]])) {
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
                          vapply(s$data, function(d) paste0(d$type, ":", d$name), character(1)))
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
