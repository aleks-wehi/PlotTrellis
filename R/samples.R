# Sample rows: one per condition, each with `data`, a list of list(type=, name=, path=)
# instances (any number of named instances per type).

# First matching instance's path for `type` (+ optional `name`), or NA if none declared; types with multiple named instances (e.g. paired-end) use sample_data_instances() below instead.
sample_data_path <- function(sample_row, type, name = NULL) {
  inst <- Filter(function(d) d$type == type && (is.null(name) || d$name == name), sample_row$data)
  if (length(inst) == 0) return(NA_character_)
  inst[[1]]$path
}

# Every instance of `type` declared on a sample (0, 1, or many); used by detect_generic() and detect_pair_generic() (track_detection.R).
sample_data_instances <- function(sample_row, type) {
  Filter(function(d) d$type == type, sample_row$data)
}

# A project's settings$samples -> sample rows (display name, data instances with full paths,
# colour). Deliberately does no file-existence check: a missing file shows as a broken track
# (track_detection.R), while the samples.txt import warns about it up front.
load_samples_from_settings <- function(samples, data_dir) {
  lapply(samples, function(s) {
    list(
      display = s$condition,
      data = lapply(s$data, function(d) list(
        # Alias: an older project_settings.json stores paired-end data as
        # type="junction" -- normalize on load so such a project keeps working
        # unedited (see samples_txt_io.R's matching alias for the samples.txt
        # import path).
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
      colour                = s$colour
    )
  })
}
