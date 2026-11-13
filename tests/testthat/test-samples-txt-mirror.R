# The wizard's Save always writes a
# project_samples.txt beside project_settings.json; a first-time wizard that finds
# it opens on Import, prefilled; an unreadable project_settings.json says how
# to rebuild from it. Drives the real wizard_server() headlessly, no fixture
# data needed.

new_wizard_project <- function() {
  proj <- tempfile("pt-proj-")
  dir.create(file.path(proj, "data", "atac"), recursive = TRUE)
  file.create(file.path(proj, "data", "atac", "a.bw"))
  proj
}
drop_wizard_project <- function(proj) unlink(proj, recursive = TRUE)
# The wizard's server for the project in `proj`, as plottrellis_wizard() builds it.
project_wizard <- function(proj) {
  json <- file.path(proj, "project_settings.json")
  defaults <- .pt_default_settings()
  wizard_server(json, .pt_read_settings(json, defaults), defaults)
}

# Runs a wizard Save of one bw sample with the given data root, its font size
# imported from a samples.txt (the wizard has no display inputs). testServer() evaluates its block in the server's environment, so the
# values are spliced in.
wizard_save <- function(proj, data_root, fontsize = 11) {
  txt <- tempfile(fileext = ".txt")
  on.exit(unlink(txt), add = TRUE)
  writeLines(c(paste0("# fontsize: ", fontsize), "condition\tcolour\tbw:ATAC", "Ba\t#1B9E77\tatac/a.bw"), txt)
  eval(bquote(shiny::testServer(project_wizard(proj), {
    session$setInputs(wiz_import_path = .(txt), wiz_import_load = 1)
    session$setInputs(wiz_data_root = .(data_root), wiz_gtf_path = "/x/genes.gtf.gz",
                      wiz_orgdb = "org.Mm.eg.db", wiz_genome_label = "mm10")
    session$setInputs(wiz_save = 1)
  })))
}

test_that("wizard Save writes a project_samples.txt that re-imports to the same project", {
  proj <- new_wizard_project()
  on.exit(drop_wizard_project(proj), add = TRUE)
  dr <- file.path(proj, "data")
  wizard_save(proj, dr)
  json <- jsonlite::fromJSON(file.path(proj, "project_settings.json"), simplifyDataFrame = FALSE)
  txt  <- file.path(proj, "project_samples.txt")
  expect_true(file.exists(txt))
  expect_false(file.exists(file.path(proj, "samples.txt")))  # never a user-named file
  res <- import_samples_txt(txt, defaults = .pt_default_settings())
  expect_length(res$errors, 0)
  expect_same_path(res$data_root, json$data_dir)
  expect_equal(res$samples[[1]]$condition, json$samples[[1]]$condition)
  expect_equal(res$samples[[1]]$data[[1]]$path, json$samples[[1]]$data[[1]]$path)
  expect_same_path(res$overrides$gtf_path, json$gtf_path)
  expect_equal(as.numeric(res$overrides$fontsize), json$fontsize)  # a non-default setting survives
})

test_that("a relative data root is written as absolute, so the copy rebuilds from anywhere", {
  proj <- new_wizard_project()
  on.exit(drop_wizard_project(proj), add = TRUE)
  old_wd <- setwd(proj)
  on.exit(setwd(old_wd), add = TRUE, after = FALSE)
  wizard_save(proj, "data")
  res <- import_samples_txt(file.path(proj, "project_samples.txt"), defaults = .pt_default_settings())
  expect_length(res$errors, 0)
  expect_same_path(res$data_root, file.path(proj, "data"))
})

test_that("a first-time wizard that finds project_samples.txt opens on Import, prefilled", {
  proj <- new_wizard_project()
  on.exit(drop_wizard_project(proj), add = TRUE)
  writeLines("condition\tcolour\tbw:ATAC", file.path(proj, "project_samples.txt"))
  shiny::testServer(project_wizard(proj), {
    expect_equal(wiz_mode(), "import")
    expect_same_path(mirror_txt, file.path(proj, "project_samples.txt"))
  })
  unlink(file.path(proj, "project_samples.txt"))
  shiny::testServer(project_wizard(proj), expect_null(wiz_mode()))
})

test_that("an unreadable project_settings.json says how to rebuild from project_samples.txt", {
  proj <- tempfile("pt-proj-"); dir.create(proj)
  on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  json <- file.path(proj, "project_settings.json")
  writeLines("{ not json", json)
  expect_error(.pt_read_settings(json, .pt_default_settings()), "run the setup wizard again")
  writeLines("condition\tcolour", file.path(proj, "project_samples.txt"))
  expect_error(.pt_read_settings(json, .pt_default_settings()), "parse_samples_txt(\"", fixed = TRUE)
})

test_that("an edit-mode wizard Save keeps the saved display settings; no change is not dirty", {
  proj <- new_wizard_project()
  on.exit(drop_wizard_project(proj), add = TRUE)
  dr <- file.path(proj, "data")
  json <- file.path(proj, "project_settings.json")
  cfg <- utils::modifyList(.pt_default_settings(),
    list(data_dir = dr, gtf_path = "/x/genes.gtf.gz", fontsize = 14, zrange_cap = 50,
         reference_colour = "#123456"))
  cfg$samples <- list(list(condition = "Ba", colour = "#1B9E77",
                           data = list(list(type = "bw", name = "ATAC", path = "atac/a.bw"))))
  writeLines(jsonlite::toJSON(cfg, auto_unbox = TRUE, pretty = TRUE, null = "null"), json)
  shiny::testServer(project_wizard(proj), {
    session$setInputs(wiz_data_root = dr, wiz_gtf_path = "/x/genes.gtf.gz",
                      wiz_orgdb = "org.Mm.eg.db", wiz_genome_label = "mm39")
    expect_false(is_dirty())
    session$setInputs(wiz_save = 1)
  })
  saved <- jsonlite::fromJSON(json, simplifyDataFrame = FALSE)
  expect_equal(saved$fontsize, 14)
  expect_equal(saved$zrange_cap, 50)
  expect_equal(saved$reference_colour, "#123456")
})

test_that("display settings imported in the wizard are saved, numbers and colours included", {
  proj <- new_wizard_project()
  on.exit(drop_wizard_project(proj), add = TRUE)
  txt <- file.path(proj, "in.txt")
  writeLines(c("# zrange_cap: 50", "# diff_colour_lo: #123456", "# reference_colour: #654321",
               "# junction_height: 1.4", "# hic_height: big",
               "condition\tcolour\tbw:ATAC", "Ba\t#1B9E77\tatac/a.bw"), txt)
  shiny::testServer(project_wizard(proj), {
    session$setInputs(wiz_import_path = txt, wiz_import_load = 1)
    session$setInputs(wiz_data_root = file.path(proj, "data"), wiz_gtf_path = "/x/genes.gtf.gz",
                      wiz_orgdb = "org.Mm.eg.db", wiz_genome_label = "mm10")
    session$setInputs(wiz_save = 1)
  })
  saved <- jsonlite::fromJSON(file.path(proj, "project_settings.json"), simplifyDataFrame = FALSE)
  expect_identical(saved$zrange_cap, 50L)  # a number, not the text "50"
  expect_equal(saved$diff_colour_lo, "#123456")
  expect_equal(saved$reference_colour, "#654321")
  expect_equal(saved$junction_height, 1.4)
  expect_null(saved$hic_height)            # not a number: stays auto
  expect_equal(saved$fontsize, 8)          # not in the file: the default
})

test_that("the wizard form does not show the display settings", {
  html <- as.character(wiz_shared_form(.pt_default_settings(), TRUE))
  for (k in setdiff(names(.pt_default_settings()), .PT_CORE_SETTINGS))
    expect_false(grepl(sprintf('id="wiz_%s"', k), html, fixed = TRUE), info = k)
  expect_true(grepl('id="wiz_gtf_path"', html, fixed = TRUE))
})

test_that("the vignette's samples.txt example is tab-separated and imports as shown", {
  rmd <- testthat::test_path("..", "..", "vignettes", "PlotTrellis.Rmd")
  skip_if_not(file.exists(rmd), "vignette source not available (installed package)")
  lines <- readLines(rmd)
  # pandoc turns tabs into spaces unless told otherwise, and the rendered page's Copy button copies what it shows.
  expect_true(any(grepl("--preserve-tabs", lines[1:12], fixed = TRUE)))
  start <- grep("^# data_root: ", lines)[1]
  block <- lines[start:(start + grep("^```", lines[start:length(lines)])[1] - 2)]
  expect_true(any(grepl("\t", block[!startsWith(block, "#")])))
  expect_false(any(grepl("  ", block[!startsWith(block, "#")])))   # no runs of spaces in the table
  root <- tempfile("pt-vig-")
  files <- sub("^[^\t]*\t[^\t]*\t", "", block[!startsWith(block, "#") & grepl("\t", block)][-1])
  for (f in unlist(strsplit(files, "\t"))) {
    dir.create(file.path(root, dirname(f)), recursive = TRUE, showWarnings = FALSE)
    file.create(file.path(root, f))
  }
  txt <- tempfile(fileext = ".txt")
  writeLines(sub("^# data_root: .*", paste("# data_root:", root), block), txt)
  res <- import_samples_txt(txt, defaults = .pt_default_settings())
  expect_length(res$errors, 0)
  expect_length(res$samples, 2)
})
