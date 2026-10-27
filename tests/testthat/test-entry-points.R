# Phase I.3 step 25: the Tier 1 entry points (parse_samples_txt(),
# plottrellis_wizard(), plottrellis_composer(), run_plottrellis()) and the
# pieces they share with the wizard. No fixture data needed, so nothing here
# skips: these run on R-universe too. Temp files only.

# A throwaway project folder holding a samples.txt plus the files it names.
make_samples_project <- function(header_lines, rows = "A\t#FF0000\ta.bw",
                                 files = c("data/a.bw", "genes.gtf")) {
  dir <- tempfile("pt-proj-")
  dir.create(dir)
  for (f in files) {
    dir.create(dirname(file.path(dir, f)), recursive = TRUE, showWarnings = FALSE)
    file.create(file.path(dir, f))
  }
  txt <- file.path(dir, "samples.txt")
  writeLines(c(header_lines, "condition\tcolour\tbw:atac", rows), txt)
  list(dir = dir, txt = txt)
}

test_that("import_samples_txt() resolves relative data_root/gtf_path against the file's folder", {
  p <- make_samples_project(c("# data_root: data", "# gtf_path: genes.gtf"))
  # R's cwd (tests/testthat) is not the file's folder, so a cwd-relative
  # resolution would both return the wrong path and warn about a.bw.
  expect_false(identical(normalizePath(getwd()), normalizePath(p$dir)))
  expect_no_warning(res <- import_samples_txt(p$txt, defaults = list(gtf_path = "")))
  expect_length(res$errors, 0)
  expect_length(res$warnings, 0)
  expect_identical(res$data_root, normalizePath(file.path(p$dir, "data"), winslash = "/"))
  expect_identical(res$overrides$gtf_path, normalizePath(file.path(p$dir, "genes.gtf"), winslash = "/"))

  abs_root <- normalizePath(file.path(p$dir, "data"), winslash = "/")
  q <- make_samples_project(paste("# data_root:", abs_root))
  expect_identical(import_samples_txt(q$txt)$data_root, abs_root)

  blank <- make_samples_project("# data_root: ", files = character(0))
  expect_identical(import_samples_txt(blank$txt)$data_root, "")
})

test_that("wiz_save_problems() returns the wizard Save's exact messages, in order", {
  expect_identical(
    wiz_save_problems(NULL, list(), "", NULL),
    c("Data root is required.", "Add at least one sample before saving.",
      "GTF path is required.", "OrgDb package is required."))
  expect_identical(
    wiz_save_problems(file.path(tempdir(), "no-such-dir"), list(list(condition = "A")), "g.gtf", "org.Mm.eg.db"),
    "Data root directory does not exist.")
  expect_identical(
    wiz_save_problems(tempdir(), list(list(condition = "A")), "g.gtf", "org.Mm.eg.db"),
    character(0))
  # Bug 62: a genome label chromosome sizes can't be looked up by; blank is allowed.
  ok <- function(label) wiz_save_problems(tempdir(), list(list(condition = "A")), "g.gtf", "org.Mm.eg.db", label)
  expect_identical(ok("mm10"), character(0))
  expect_identical(ok(""), character(0))
  expect_identical(ok("m10"), "Genome label 'm10' isn't a known genome build id (e.g. mm10, mm39, hg38).")
})

test_that("parse_samples_txt() rejects a genome label that isn't a genome build id (bug 62)", {
  p <- make_samples_project(c("# data_root: data", "# gtf_path: genes.gtf", "# orgdb: org.Mm.eg.db", "# genome_label: m10"))
  out_dir <- tempfile("pt-out-")
  expect_error(parse_samples_txt(p$txt, project_dir = out_dir), "Genome label 'm10'", fixed = TRUE)
  expect_false(file.exists(file.path(out_dir, "project_settings.json")))
})

test_that("plottrellis_composer() errors, creating nothing, when the project has no settings", {
  dir <- tempfile("pt-none-")
  expect_error(plottrellis_composer(dir), "No project_settings.json.*parse_samples_txt")
  expect_false(dir.exists(dir))
})

# Runs run_plottrellis() with the two launchers mocked; each mock is a list of
# per-call functions (called in turn) and every call is recorded in order.
run_loop_with_mocks <- function(dir, wizard_steps, composer_steps) {
  old <- Sys.getenv("PLOTTRELLIS_EDIT_SETTINGS", NA)
  on.exit(if (is.na(old)) Sys.unsetenv("PLOTTRELLIS_EDIT_SETTINGS")
          else Sys.setenv(PLOTTRELLIS_EDIT_SETTINGS = old), add = TRUE)
  Sys.unsetenv("PLOTTRELLIS_EDIT_SETTINGS")
  calls <- character(0)
  n <- c(wizard = 0, composer = 0)
  step <- function(who, steps) function(project_dir, ...) {
    calls <<- c(calls, who)
    n[[who]] <<- n[[who]] + 1
    steps[[n[[who]]]](project_dir)
  }
  testthat::local_mocked_bindings(
    plottrellis_wizard   = step("wizard", wizard_steps),
    plottrellis_composer = step("composer", composer_steps)
  )
  result <- run_plottrellis(dir)
  list(calls = calls, result = result)
}

test_that("run_plottrellis() goes wizard -> composer on a first run with no manual step", {
  dir <- tempfile("pt-loop-")
  dir.create(dir)
  json <- file.path(dir, "project_settings.json")
  out <- run_loop_with_mocks(
    dir,
    wizard_steps   = list(function(d) { writeLines("{}", json); "restart" }),
    composer_steps = list(function(d) "closed"))
  expect_identical(out$calls, c("wizard", "composer"))
  expect_identical(out$result, "closed")
})

test_that("run_plottrellis() reopens the wizard after Edit project settings, then the composer", {
  dir <- tempfile("pt-loop-")
  dir.create(dir)
  writeLines("{}", file.path(dir, "project_settings.json"))
  out <- run_loop_with_mocks(
    dir,
    wizard_steps = list(function(d) { Sys.unsetenv("PLOTTRELLIS_EDIT_SETTINGS"); "restart" }),
    composer_steps = list(
      function(d) { Sys.setenv(PLOTTRELLIS_EDIT_SETTINGS = "1"); "edit_settings" },
      function(d) "closed"))
  expect_identical(out$calls, c("composer", "wizard", "composer"))
  expect_identical(out$result, "closed")
})

test_that("run_plottrellis() opens the wizard after Reset deletes the settings", {
  dir <- tempfile("pt-loop-")
  dir.create(dir)
  json <- file.path(dir, "project_settings.json")
  writeLines("{}", json)
  out <- run_loop_with_mocks(
    dir,
    wizard_steps   = list(function(d) "closed"),
    composer_steps = list(function(d) { file.remove(json); "restart" }))
  expect_identical(out$calls, c("composer", "wizard"))
  expect_identical(out$result, "closed")
})

test_that("PlotTrellisApp() passes everything through to run_plottrellis()", {
  got <- NULL
  testthat::local_mocked_bindings(run_plottrellis = function(...) { got <<- list(...); "closed" })
  expect_identical(PlotTrellisApp("x", port = 1), "closed")
  expect_identical(got, list(project_dir = "x", port = 1))
})

# A valid relative samples.txt, as parse_samples_txt() expects to see one.
good_header <- c("# data_root: data", "# gtf_path: genes.gtf", "# orgdb: org.Hs.eg.db")

test_that("parse_samples_txt() writes settings the composer reads back, samples included", {
  p <- make_samples_project(good_header)
  out_dir <- tempfile("pt-out-")
  written <- parse_samples_txt(p$txt, project_dir = out_dir)
  json <- file.path(out_dir, "project_settings.json")
  expect_true(file.exists(json))
  cfg <- .pt_read_settings(json, .pt_default_settings())
  expect_identical(cfg$data_dir, normalizePath(file.path(p$dir, "data"), winslash = "/"))
  expect_identical(cfg$gtf_path, normalizePath(file.path(p$dir, "genes.gtf"), winslash = "/"))
  expect_identical(cfg$orgdb, "org.Hs.eg.db")
  # Guards the modifyList() trap: samples must survive into the json.
  expect_length(cfg$samples, 1)
  expect_identical(cfg$samples[[1]]$condition, "A")
  expect_identical(cfg$samples[[1]]$data[[1]]$path, "a.bw")
  expect_equal(cfg$fontsize, 8)
  expect_identical(written$samples, cfg$samples)
})

test_that("parse_samples_txt() stops on a broken samples.txt, writing and creating nothing", {
  p <- make_samples_project(good_header, rows = c("A\t#FF0000\ta.bw", "A\t#00FF00\ta.bw"))
  out_dir <- tempfile("pt-out-")
  expect_error(parse_samples_txt(p$txt, project_dir = out_dir), "Duplicate condition 'A'")
  expect_false(dir.exists(out_dir))
})

test_that("parse_samples_txt() applies the wizard's checks (missing gtf_path)", {
  p <- make_samples_project(c("# data_root: data", "# orgdb: org.Hs.eg.db"))
  out_dir <- tempfile("pt-out-")
  expect_error(parse_samples_txt(p$txt, project_dir = out_dir), "GTF path is required.", fixed = TRUE)
  expect_false(file.exists(file.path(out_dir, "project_settings.json")))
})

test_that("parse_samples_txt() never silently replaces existing settings", {
  p <- make_samples_project(good_header)
  out_dir <- tempfile("pt-out-")
  parse_samples_txt(p$txt, project_dir = out_dir)
  json <- file.path(out_dir, "project_settings.json")
  before <- unname(tools::md5sum(json))

  q <- make_samples_project(good_header, rows = "B\t#00FF00\ta.bw")
  expect_error(parse_samples_txt(q$txt, project_dir = out_dir), "already exists")
  expect_identical(unname(tools::md5sum(json)), before)

  parse_samples_txt(q$txt, project_dir = out_dir, overwrite = TRUE)
  cfg <- .pt_read_settings(json, .pt_default_settings())
  expect_identical(cfg$samples[[1]]$condition, "B")
})

test_that("parse_samples_txt() warns about unresolved sample files but still saves", {
  p <- make_samples_project(good_header, rows = c("A\t#FF0000\ta.bw", "B\t#00FF00\tmissing.bw"))
  out_dir <- tempfile("pt-out-")
  expect_warning(parse_samples_txt(p$txt, project_dir = out_dir), "not found")
  json <- file.path(out_dir, "project_settings.json")
  expect_true(file.exists(json))
  expect_length(.pt_read_settings(json, .pt_default_settings())$samples, 2)
})

test_that("parse_samples_txt() saves NULL-default numeric settings as numbers, like the wizard", {
  # zrange_cap/hic_height default to NULL ("auto"), so the importer can't
  # infer their type and returns text; the composer's
  # clean_optional_numeric() would silently treat "50" as auto.
  p <- make_samples_project(c(good_header, "# zrange_cap: 50", "# hic_height: 2.5"))
  out_dir <- tempfile("pt-out-")
  parse_samples_txt(p$txt, project_dir = out_dir)
  cfg <- .pt_read_settings(file.path(out_dir, "project_settings.json"), .pt_default_settings())
  expect_identical(cfg$zrange_cap, 50L)
  expect_identical(cfg$hic_height, 2.5)
})

test_that("parse_samples_txt() refreshes project_samples.txt, so recovery restores the new samples", {
  # Reset and the broken-settings error both point at this copy (feature 43);
  # found by step 25's review: an overwrite left the old samples in it.
  p <- make_samples_project(good_header)
  out_dir <- tempfile("pt-out-")
  parse_samples_txt(p$txt, project_dir = out_dir)
  q <- make_samples_project(good_header, rows = "B\t#00FF00\ta.bw")
  parse_samples_txt(q$txt, project_dir = out_dir, overwrite = TRUE)
  mirror <- import_samples_txt(file.path(out_dir, "project_samples.txt"), defaults = .pt_default_settings())
  expect_length(mirror$errors, 0)
  expect_identical(vapply(mirror$samples, `[[`, "", "condition"), "B")
  expect_identical(mirror$data_root, normalizePath(file.path(q$dir, "data"), winslash = "/"))
})

test_that("parse_samples_txt() warns when an auto setting isn't a finite number", {
  p <- make_samples_project(c(good_header, "# hic_height: Inf"))
  out_dir <- tempfile("pt-out-")
  expect_warning(parse_samples_txt(p$txt, project_dir = out_dir), "hic_height: Inf")
  cfg <- .pt_read_settings(file.path(out_dir, "project_settings.json"), .pt_default_settings())
  expect_null(cfg$hic_height)
})

test_that("parse_samples_txt() says which genome fields came from the defaults", {
  p <- make_samples_project(c("# data_root: data", "# gtf_path: genes.gtf"))
  expect_message(parse_samples_txt(p$txt, project_dir = tempfile("pt-out-")), "orgdb = org.Mm.eg.db")
})
