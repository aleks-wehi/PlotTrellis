# The wizard's Import a samples.txt (.pt_wizard_import_server(), R/wizard-import.R) on its own:
# plain reactiveVals for what an import fills, a samples.txt in a temp folder, no genome. The
# field updates it sends are caught, not sent.

import_app <- function(input, output, session) {
  wiz_samples   <- reactiveVal(list())
  imported      <- reactiveVal(NULL)
  import_visual <- reactiveVal(list())
  .pt_wizard_import_server(input, output, session, wiz_samples, imported, import_visual, .pt_default_settings())
}

# Catches the field updates the code under test sends (id = value), for the rest of the
# calling test.
catch_updates <- function(env = parent.frame()) {
  sent <- new.env()
  local_mocked_bindings(updateTextInput = function(session, inputId, label = NULL, value = NULL, ...)
    assign(inputId, value, envir = sent), .env = env)
  sent
}

status_html <- function(out) if (is.null(out)) "" else paste(as.character(out$html), collapse = "")

new_import_project <- function() {
  root <- tempfile("pt-import-")
  dir.create(file.path(root, "atac"), recursive = TRUE)
  file.create(file.path(root, "atac", "a.bw"))
  root
}

test_that("a blank path or a missing file is said, and nothing is imported", {
  shiny::testServer(import_app, {
    expect_identical(status_html(output$wiz_import_status), "")
    session$setInputs(wiz_import_path = "  ", wiz_import_load = 1)
    expect_match(status_html(output$wiz_import_status), "Enter a path to samples.txt first.", fixed = TRUE)
    session$setInputs(wiz_import_path = file.path(tempdir(), "no-such-samples.txt"), wiz_import_load = 2)
    expect_match(status_html(output$wiz_import_status), "No such file:", fixed = TRUE)
    expect_length(wiz_samples(), 0)
    expect_null(imported())
    expect_identical(import_visual(), list())
  })
})

test_that("a file with errors imports nothing and lists them", {
  root <- new_import_project(); on.exit(unlink(root, recursive = TRUE), add = TRUE)
  txt <- file.path(root, "samples.txt")
  writeLines(c("condition\tcolour\tbw:ATAC", "\t#1B9E77\tatac/a.bw"), txt)
  sent <- catch_updates()
  shiny::testServer(import_app, {
    session$setInputs(wiz_import_path = txt, wiz_import_load = 1)
    expect_match(status_html(output$wiz_import_status), "Nothing imported. The file has these problems:", fixed = TRUE)
    expect_length(wiz_samples(), 0)
    expect_null(imported())
    expect_identical(import_visual(), list())
  })
  expect_length(ls(sent), 0)
})

test_that("an import fills the samples, the data folder and genome fields, and keeps its display settings", {
  root <- new_import_project(); on.exit(unlink(root, recursive = TRUE), add = TRUE)
  txt <- file.path(root, "samples.txt")
  lines <- c(paste("# data_root:", root), "# gtf_path: /x/genes.gtf.gz", "# orgdb: org.Mm.eg.db",
             "# genome_label: mm10", "# fontsize: 11", "# zrange_cap: 50",
             "condition\tcolour\tbw:ATAC", "Ba\t#1B9E77\tatac/a.bw")
  writeLines(lines, txt)
  sent <- catch_updates()
  shiny::testServer(import_app, {
    session$setInputs(wiz_import_path = txt, wiz_import_load = 1)
    expect_match(status_html(output$wiz_import_status), "Imported 1 sample(s).", fixed = TRUE)
    expect_no_match(status_html(output$wiz_import_status), "with warnings", fixed = TRUE)
    expect_identical(vapply(wiz_samples(), `[[`, "", "condition"), "Ba")
    expect_identical(wiz_samples()[[1]]$data[[1]][c("type", "name", "path")],
                     list(type = "bw", name = "ATAC", path = "atac/a.bw"))
    expect_identical(imported(), lines)
    # Display settings only; the genome fields go to their inputs, numbers stay numbers.
    expect_setequal(names(import_visual()), c("fontsize", "zrange_cap"))
    expect_equal(import_visual()$fontsize, 11)
    expect_equal(import_visual()$zrange_cap, 50)
  })
  expect_identical(mget(c("wiz_data_root", "wiz_gtf_path", "wiz_orgdb", "wiz_genome_label"), envir = sent),
                   list(wiz_data_root = root, wiz_gtf_path = "/x/genes.gtf.gz", wiz_orgdb = "org.Mm.eg.db",
                        wiz_genome_label = "mm10"))
})

test_that("an import with problems it can skip says so in an amber note", {
  root <- new_import_project(); on.exit(unlink(root, recursive = TRUE), add = TRUE)
  txt <- file.path(root, "samples.txt")
  writeLines(c(paste("# data_root:", root), "# hic_height: big",
               "condition\tcolour\tbw:ATAC", "Ba\t#1B9E77\tatac/missing.bw"), txt)
  sent <- catch_updates()
  shiny::testServer(import_app, {
    session$setInputs(wiz_import_path = txt, wiz_import_load = 1)
    html <- status_html(output$wiz_import_status)
    expect_match(html, "Imported 1 sample(s).", fixed = TRUE)
    expect_match(html, "Imported, with warnings:", fixed = TRUE)
    expect_match(html, "wiz-warn-box", fixed = TRUE)
    expect_match(html, "hic_height", fixed = TRUE)
    expect_null(import_visual()$hic_height)
  })
  # No genome lines in the file: those fields are left as they are.
  expect_identical(ls(sent), "wiz_data_root")
})
