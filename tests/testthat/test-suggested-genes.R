# The project setting suggested_genes (a samples.txt's "# suggested_genes:" line): how it
# travels through import, the project files, the wizard and export. What the gene box does
# with it is in test-main-api.R.

# A copy of the mammary example's samples.txt with `line` as its suggested_genes line in
# place of the example's own (NULL: none).
suggesting_example <- function(line = NULL) {
  txt <- example_copy("mammary")
  lines <- grep("^# suggested_genes:", readLines(txt), value = TRUE, invert = TRUE)
  at <- grep("^# genome_label:", lines)
  writeLines(append(lines, line, after = at), txt)
  txt
}

copy_lines <- function(proj) readLines(file.path(proj, "project_samples.txt"))

test_that("a samples.txt's suggested_genes line becomes the project's setting, as written", {
  skip_if_not_installed("org.Mm.eg.db")
  txt <- suggesting_example("# suggested_genes: Cat, ENSMUSG00000027186 ,elf5")
  proj <- tempfile("pt-sugg-")
  on.exit(unlink(c(proj, dirname(txt)), recursive = TRUE), add = TRUE)
  settings <- suppressMessages(parse_samples_txt(txt, project_dir = proj))
  expect_identical(settings$suggested_genes, "Cat, ENSMUSG00000027186 ,elf5")
  json <- jsonlite::read_json(file.path(proj, "project_settings.json"))
  expect_identical(json$suggested_genes, "Cat, ENSMUSG00000027186 ,elf5")
  # The project's own copy carries the line, so rebuilding from it keeps the list.
  expect_true("# suggested_genes: Cat, ENSMUSG00000027186 ,elf5" %in% copy_lines(proj))
  again <- tempfile("pt-sugg-")
  on.exit(unlink(again, recursive = TRUE), add = TRUE)
  rebuilt <- suppressMessages(parse_samples_txt(file.path(proj, "project_samples.txt"), project_dir = again))
  expect_identical(rebuilt$suggested_genes, settings$suggested_genes)
})

test_that("without the line the setting is blank and no line is written", {
  skip_if_not_installed("org.Mm.eg.db")
  txt <- suggesting_example(NULL)
  proj <- tempfile("pt-sugg-")
  on.exit(unlink(c(proj, dirname(txt)), recursive = TRUE), add = TRUE)
  settings <- suppressMessages(parse_samples_txt(txt, project_dir = proj))
  expect_identical(settings$suggested_genes, "")
  expect_false(any(grepl("suggested_genes", copy_lines(proj))))
  expect_identical(.pt_default_settings()$suggested_genes, "")
})

test_that("it is a core setting: no Figure settings input, and Save as project default keeps it", {
  skip_if_not_installed("org.Mm.eg.db")
  expect_true("suggested_genes" %in% .PT_CORE_SETTINGS)
  expect_false("suggested_genes" %in% names(.pt_display_defaults()))
  txt <- suggesting_example("# suggested_genes: Cat, Elf5")
  proj <- tempfile("pt-sugg-")
  on.exit(unlink(c(proj, dirname(txt)), recursive = TRUE), add = TRUE)
  suppressMessages(parse_samples_txt(txt, project_dir = proj))
  json <- file.path(proj, "project_settings.json")
  cfg <- .pt_read_settings(json, .pt_default_settings())
  expect_identical(cfg$suggested_genes, "Cat, Elf5")
  # The Figure settings block's Save, on render settings that do not carry the list.
  render_cfg <- cfg[setdiff(names(cfg), "suggested_genes")]
  shiny::testServer(function(input, output, session)
    fs <- .pt_figure_settings_server(input, output, render_cfg, json), {
    session$setInputs(s_fontsize = 12, save_settings = 1)
  })
  saved <- jsonlite::read_json(json)
  expect_identical(saved$fontsize, 12L)
  expect_identical(saved$suggested_genes, "Cat, Elf5")
  expect_true("# suggested_genes: Cat, Elf5" %in% copy_lines(proj))
})

test_that("the wizard carries it: the project's list, or the one of a samples.txt imported there", {
  cfg <- utils::modifyList(.pt_default_settings(), list(suggested_genes = "Cat, Elf5"))
  root <- tempfile("pt-sugg-wiz-")
  dir.create(file.path(root, "atac"), recursive = TRUE)
  file.create(file.path(root, "atac", "a.bw"))
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  txt <- file.path(root, "samples.txt")
  body <- c(paste("# data_root:", root), "# gtf_path: /x/genes.gtf.gz", "# orgdb: org.Mm.eg.db",
            "# genome_label: mm10", "condition\tcolour\tbw:atac", "A\t#1B9E77\tatac/a.bw")
  local_mocked_bindings(updateTextInput = function(...) invisible())
  out <- file.path(root, "exported.txt")
  shiny::testServer(function(input, output, session) {
    wiz_samples <- reactiveVal(list())
    import_visual <- reactiveVal(list())
    .pt_wizard_import_server(input, output, session, wiz_samples, import_visual, .pt_default_settings())
    blk <- .pt_wizard_export_server(input, output, function() root, wiz_samples, import_visual, cfg,
                                    .pt_default_settings())
  }, {
    # as opened: the project's
    expect_identical(blk$wiz_collect_flat_settings()$suggested_genes, "Cat, Elf5")
    # an imported file without the line leaves it (as it leaves the display settings)
    writeLines(body, txt)
    session$setInputs(wiz_import_path = txt, wiz_import_load = 1)
    expect_length(wiz_samples(), 1)
    expect_identical(blk$wiz_collect_flat_settings()$suggested_genes, "Cat, Elf5")
    # one with the line replaces it
    writeLines(append(body, "# suggested_genes: Ehf", after = 4), txt)
    session$setInputs(wiz_import_path = txt, wiz_import_load = 2)
    expect_identical(blk$wiz_collect_flat_settings()$suggested_genes, "Ehf")
    # and Export writes it
    session$setInputs(wiz_export_path = out, wiz_confirm_export = 1)
    expect_true("# suggested_genes: Ehf" %in% readLines(out))
  })
  # a project that never had one: blank, and no line exported
  shiny::testServer(function(input, output, session) {
    import_visual <- reactiveVal(list())
    blk <- .pt_wizard_export_server(input, output, function() root, function() list(), import_visual,
                                    .pt_default_settings(), .pt_default_settings())
  }, {
    expect_identical(blk$wiz_collect_flat_settings()$suggested_genes, "")
    session$setInputs(wiz_export_path = out, wiz_confirm_export = 1)
    expect_false(any(grepl("suggested_genes", readLines(out))))
  })
})

test_that("the app says at launch which suggested genes the annotation does not have", {
  skip_if_not_installed("org.Mm.eg.db")
  txt <- suggesting_example("# suggested_genes: Eya2, NoSuchGene, AlsoNot")  # Eya2: in the example
  proj <- tempfile("pt-sugg-")
  on.exit(unlink(c(proj, dirname(txt)), recursive = TRUE), add = TRUE)
  local_offline_chrom_sizes()
  p <- suppressMessages(pt_project(proj, txt))
  withr::local_options(PlotTrellis.stop_on_close = FALSE)
  prj <- .pt_project(proj)
  expect_warning(
    app <- suppressMessages(.pt_composer_app(prj$settings_json, prj$cfg)),
    "Not suggested, as this project's annotation has no such gene: NoSuchGene, AlsoNot (# suggested_genes:).",
    fixed = TRUE
  )
  expect_s3_class(app, "shiny.appobj")
})
