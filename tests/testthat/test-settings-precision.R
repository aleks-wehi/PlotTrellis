# Bug 174: every writer of project_settings.json keeps a number's full precision (jsonlite's
# toJSON() rounds to 4 decimals by default). No fixture data needed.

# The example's samples.txt with absolute paths and one setting with many decimals.
precise_samples_txt <- function() {
  ex <- dirname(pt_example("mammary"))
  lines <- readLines(pt_example("mammary"))
  lines <- sub("^# data_root: .*", paste0("# data_root: ", ex), lines)
  lines <- sub("^# gtf_path: (.*)", paste0("# gtf_path: ", file.path(ex, "\\1")), lines)
  lines <- c("# junction_height: 1.234567", lines)
  out <- tempfile(fileext = ".txt")
  writeLines(lines, out)
  out
}

test_that("parse_samples_txt() and pt_project()'s re-import keep every decimal", {
  skip_if_not_installed("org.Mm.eg.db")
  json <- file.path(project <- tempfile("pt-precision-"), "project_settings.json")
  suppressMessages(suppressWarnings(parse_samples_txt(precise_samples_txt(), project_dir = project)))
  expect_identical(jsonlite::read_json(json)$junction_height, 1.234567)
  # .pt_import() (pt_project(), run_plottrellis(samples_txt =)) writes the file again.
  project2 <- tempfile("pt-precision-")
  suppressMessages(suppressWarnings(.pt_import(project2, precise_samples_txt())))
  expect_identical(jsonlite::read_json(file.path(project2, "project_settings.json"))$junction_height, 1.234567)
})

test_that("the wizard's and the app's save keep every decimal", {
  project <- tempfile("pt-precision-")
  dir.create(project)
  json <- file.path(project, "project_settings.json")
  .pt_write_project_files(list(data_dir = project, samples = list(), junction_height = 1.234567,
                               zrange_cap = 123.456789), json)
  expect_identical(jsonlite::read_json(json)$junction_height, 1.234567)
  expect_identical(jsonlite::read_json(json)$zrange_cap, 123.456789)
})
