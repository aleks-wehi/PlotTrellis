# rationale: see solutions/concepts/19-relocated-code-rationale.md, "utils/samples_txt_io.R:1"

# %||% is NA-aware and defined once, package-wide, in boilerplate.R (see
# phase-G.1.md's log) -- not a second local definition here. NA-specific
# defaulting for a possibly-NA character value still goes through chr0()
# below instead, since %||% would treat "" and NA the same way (both
# defaulted), which isn't always the desired distinction for a field read
# straight out of samples.txt.
chr0 <- function(x) if (length(x) == 0 || is.na(x)) "" else x

# Must resolve to an actual FILE, not a directory (file.exists() is TRUE for both) -- a directory here would crash is_paired_end_structural() or silently pass only to fail at render.
is_real_file <- function(path) file.exists(path) && !dir.exists(path)

DATA_TYPE_VOCAB <- c("bw", "bed", "hic", "quantdir", "paired-end")

# Settings keys that never round-trip through samples.txt's `# <key>: value`
# delta-header lines (Phase C.3 decision 6/7, updated by the standalone
# genome-config task): `data_dir` already has its own dedicated
# `# data_root:` line; `samples` is structural, not a render-style setting.
# Genome fields (gtf_path/orgdb/genome_label) are independent per-project
# settings with no genome_id/GENOME_TABLE indirection, so -- unlike before --
# they DO round-trip through the generic delta mechanism below (on the
# IMPORT side unchanged; see SAMPLES_TXT_ALWAYS_EXPORT below for how export
# treats them differently). Everything else in DEFAULT_SETTINGS is fair game.
SAMPLES_TXT_DELTA_EXCLUDE <- c("data_dir", "samples")

# Genome identity fields ALWAYS get a `# key: value` line on export,
# regardless of whether the live value happens to equal DEFAULT_SETTINGS'
# (found 2026-08-25: a project using this deployment's own default mm39
# genome silently exported NO genome info at all, since delta-only export
# treats "matches the default" as "nothing to say" -- correct for a cosmetic
# setting like page_w, wrong for identifying, must-travel-with-the-data
# information a biologist needs to see regardless of which machine produced
# the file). Still parsed on import via the same generic `# <key>: value`
# mechanism as every other key -- only the export side special-cases these.
SAMPLES_TXT_ALWAYS_EXPORT <- c("gtf_path", "orgdb", "genome_label")

# Coerces `val` to match `default_val`'s type (Phase C.3 step A, e.g. a logical default like uniprot_reviewed_only must come back logical, not "TRUE"/"FALSE"); malformed input coerces to NA, which the caller drops rather than injecting into live settings.
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

# Cheap structural check only (column count matches BED12's 12 or BEDPE's 6-10, not full semantic validation, see read_junction_file() in utils/transcript_filter.R and solutions/concepts/18-bed12-bedpe-arc-format.md, per phase-C.1.md step 5); caller must check is_real_file() first, this doesn't guard against `path` being a directory.
is_paired_end_structural <- function(path) {
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

# Parses samples.txt into list(errors=, warnings=, data_root=, overrides=, samples=) (samples[[i]]'s data shape matches utils/samples.R's sample$data); errors are fatal and callers must check them before trusting samples, warnings are non-fatal broken paths; data_root_override validates against the running CFG$data_dir instead of the file's own `# data_root:` line. See solutions/concepts/19-relocated-code-rationale.md, "utils/samples_txt_io.R:107".
import_samples_txt <- function(path, data_root_override = NULL, defaults = list()) {
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
  # working directory, so a samples.txt can travel with its data (phase-I.3 step 25).
  # Absolute paths (or ~) and blanks are left exactly as written. Inline base R only: this
  # function's text is pasted into every downloaded script (see .pt_download_deps()).
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

  effective_root <- data_root_override %||% data_root

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
      # Pre-rename alias: samples.txt files written before the BEDPE-support
      # generalization used "junction:<name>" -- normalize so they keep
      # loading without a manual edit (solutions/concepts/18-bed12-bedpe-arc-format.md).
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
        if (!is.null(effective_root)) {
          resolved <- file.path(effective_root, val)
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
  # paths, Phase C.3 decision 5b) do NOT trigger this: those rows are still
  # genuinely usable, just missing one file, so every OTHER sample still
  # comes in per Extra check 2a.
  if (length(errors) > 0) samples <- list()

  list(errors = errors, warnings = warnings, data_root = data_root,
       overrides = overrides, samples = samples)
}

# Mirrors import_samples_txt(): flattens settings$samples[] into a wide type:name-header samples.txt plus a `# data_root:` line and one `# <key>: value` delta line per non-default setting (Phase C.3 decision 6), except the 3 genome fields which always get a line (SAMPLES_TXT_ALWAYS_EXPORT); pass DEFAULT_SETTINGS as `defaults`. See solutions/concepts/19-relocated-code-rationale.md, "utils/samples_txt_io.R:269".
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
    # "still unset". Same failure shape as BUGS.md #13, guarded against here
    # for the character case specifically.
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
