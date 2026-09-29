# rationale: see solutions/concepts/19-relocated-code-rationale.md, "utils/samples.R:1"

# First matching instance's path for `type` (+ optional `name`), or NA if none declared -- a drop-in replacement for old fixed-column reads (e.g. sample_row$rnaseq_junction_path); types with multiple named instances (e.g. paired-end) use sample_data_instances() below instead (solutions/concepts/18-bed12-bedpe-arc-format.md).
sample_data_path <- function(sample_row, type, name = NULL) {
  inst <- Filter(function(d) d$type == type && (is.null(name) || d$name == name), sample_row$data)
  if (length(inst) == 0) return(NA_character_)
  inst[[1]]$path
}

# Every instance of `type` declared on a sample (0, 1, or many); used by modules supporting multiple named instances (viz/bw.R, viz/hic.R, viz/bed.R's/viz/sashimi.R's detect_all()).
sample_data_instances <- function(sample_row, type) {
  Filter(function(d) d$type == type, sample_row$data)
}

# Phase C.1 step 6: a samples.txt with any ":"-containing header is NEW-format, produced by the wizard's add-sample/export flow -- load via import_samples_txt() (utils/samples_txt_io.R) instead of the fixed-column parser below; Samples_super_merged.txt's own on-disk format is untouched.
load_samples_new_format <- function(samples_super_txt, data_dir) {
  res <- import_samples_txt(samples_super_txt, data_root_override = data_dir)
  if (length(res$errors) > 0) {
    stop("Failed to load '", samples_super_txt, "':\n", paste(res$errors, collapse = "\n"))
  }
  lapply(res$samples, function(s) {
    list(
      display = s$condition,
      data = lapply(s$data, function(d) list(
        type = d$type, name = d$name,
        path = if (is.na(d$path) || nchar(d$path) == 0) {
          NA_character_
        } else if (identical(d$type, "quantdir")) {
          file.path(data_dir, d$path, "quant.sf")
        } else {
          file.path(data_dir, d$path)
        }
      )),
      colour                = s$colour,
      rnaseq_condition_code = NA_character_,
      cell_type             = NA_character_,
      memory_naive          = NA_character_,
      resting_activation    = NA_character_
    )
  })
}

# Deliberately does NO file-existence validation (decision 5a) -- unlike load_samples_new_format(), which still validates because Import is a one-shot, human-supervised action worth flagging problems on immediately. Full rationale: solutions/concepts/19-relocated-code-rationale.md, "utils/samples.R:80".
load_samples_from_settings <- function(samples, data_dir) {
  lapply(samples, function(s) {
    list(
      display = s$condition,
      data = lapply(s$data, function(d) list(
        # Pre-rename alias: a project_settings.json saved before the
        # BEDPE-support generalization stored type="junction" -- normalize on
        # load so an already-saved project keeps working unedited (see
        # samples_txt_io.R's matching alias for the samples.txt import path).
        type = if (identical(d$type, "junction")) "paired-end" else d$type,
        name = d$name,
        path = if (is.null(d$path) || is.na(d$path) || nchar(d$path) == 0) {
          NA_character_
        } else if (identical(d$type, "quantdir")) {
          file.path(data_dir, d$path, "quant.sf")
        } else {
          file.path(data_dir, d$path)
        }
      )),
      colour                = s$colour,
      rnaseq_condition_code = NA_character_,
      cell_type             = NA_character_,
      memory_naive          = NA_character_,
      resting_activation    = NA_character_
    )
  })
}

load_samples <- function(samples_super_txt, data_dir) {
  header_line <- readLines(samples_super_txt, n = 1, warn = FALSE)
  if (any(grepl(":", strsplit(header_line, "\t", fixed = TRUE)[[1]], fixed = TRUE))) {
    return(load_samples_new_format(samples_super_txt, data_dir))
  }

  raw <- read.delim(samples_super_txt, header = TRUE, stringsAsFactors = FALSE,
                     check.names = FALSE, sep = "\t")
  raw$condition <- trimws(raw$condition)

  join_path <- function(rel) {
    rel <- trimws(rel)
    ifelse(!is.na(rel) & nchar(rel) > 0, file.path(data_dir, rel), NA_character_)
  }
  raw$atac_path            <- join_path(raw$atac_path)
  raw$chip_path            <- join_path(raw$chip_path)
  raw$cutrun_path          <- join_path(raw$cutrun_path)
  raw$hic_path             <- join_path(raw$hic_path)
  # bed_path: spike field for viz/bed.R (13-bed-support-spike-prompt.md).
  raw$bed_path             <- if ("bed_path" %in% names(raw)) join_path(raw$bed_path) else NA_character_
  raw$rnaseq_coverage_path <- join_path(raw$rnaseq_coverage_path)
  raw$rnaseq_junction_path <- join_path(raw$rnaseq_junction_path)

  quant_dir <- trimws(raw$rnaseq_quant_dir)
  raw$rnaseq_quant_path <- ifelse(
    !is.na(quant_dir) & nchar(quant_dir) > 0,
    file.path(data_dir, quant_dir, "quant.sf"), NA_character_
  )

  raw$rnaseq_condition_code <- trimws(raw$rnaseq_condition_code)
  raw$rnaseq_condition_code[nchar(raw$rnaseq_condition_code) == 0] <- NA_character_

  # Optional colour column -- valid hex codes become sample defaults
  raw$colour <- NA_character_
  if ("colour" %in% names(raw)) {
    rawcol <- trimws(raw$colour)
    valid  <- grepl("^#[0-9A-Fa-f]{6}$", rawcol)
    raw$colour[valid] <- rawcol[valid]
  }

  # Migration adapter: fixed columns -> `data` instances; names are fixed, human-readable strings chosen to match today's hardcoded track-label text exactly (decision 13). Full rationale: solutions/concepts/19-relocated-code-rationale.md, "utils/samples.R:160".
  build_data <- function(i) {
    specs <- list(
      list(type = "bw",       name = "ATAC",      path = raw$atac_path[i]),
      list(type = "bw",       name = "ChIP",      path = raw$chip_path[i]),
      list(type = "bw",       name = "CUT&RUN",   path = raw$cutrun_path[i]),
      list(type = "bw",       name = "RNA-seq",   path = raw$rnaseq_coverage_path[i]),
      list(type = "hic",      name = "HiC",       path = raw$hic_path[i]),
      list(type = "bed",      name = "BED",       path = raw$bed_path[i]),
      list(type = "paired-end", name = "Junctions", path = raw$rnaseq_junction_path[i]),
      list(type = "quantdir", name = "Quant",     path = raw$rnaseq_quant_path[i])
    )
    Filter(function(s) !is.na(s$path), specs)
  }

  lapply(seq_len(nrow(raw)), function(i) {
    list(
      display               = raw$condition[i],
      data                  = build_data(i),
      colour                = raw$colour[i],
      rnaseq_condition_code = raw$rnaseq_condition_code[i],
      cell_type             = raw$cell_type[i],
      memory_naive          = raw$memory_naive[i],
      resting_activation    = raw$resting_activation[i]
    )
  })
}
