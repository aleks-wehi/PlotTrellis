# The wizard's Export these settings as a samples.txt (.pt_wizard_export_server(),
# R/wizard-export.R) on its own, with the settings it collects for Save and the dirty guard:
# plain reactives, a temp folder, no genome. The export dialog is caught, not shown.

# A server made of nothing but the export block; `root` is the Data folder, `visual` the
# display settings an import kept.
export_app <- function(root, samples = list(), visual = list(), cfg = .pt_default_settings()) {
  function(input, output, session) {
    data_root     <- reactive(root)
    wiz_samples   <- reactiveVal(samples)
    import_visual <- reactiveVal(visual)
    res <- .pt_wizard_export_server(input, output, data_root, wiz_samples, import_visual, cfg, .pt_default_settings())
  }
}

catch_modals <- function(env = parent.frame()) {
  got <- new.env(); got$shown <- character(0)
  local_mocked_bindings(showModal = function(ui, ...) got$shown <- c(got$shown, as.character(ui)), .env = env)
  got
}

status_html <- function(out) if (is.null(out)) "" else paste(as.character(out$html), collapse = "")

test_that("the collected settings: the genome fields (blank is NULL), the suggested genes, then cfg's display settings with imported ones on top", {
  cfg <- utils::modifyList(.pt_default_settings(), list(fontsize = 14, reference_colour = "#123456"))
  shiny::testServer(export_app(NULL, visual = list(fontsize = 11), cfg = cfg), {
    session$setInputs(wiz_gtf_path = " /x/genes.gtf.gz ", wiz_orgdb = "", wiz_genome_label = "mm10")
    flat <- res$wiz_collect_flat_settings()
    expect_identical(names(flat)[1:4], c("gtf_path", "orgdb", "genome_label", "suggested_genes"))
    expect_identical(flat$gtf_path, "/x/genes.gtf.gz")
    expect_null(flat$orgdb)
    expect_identical(flat$genome_label, "mm10")
    expect_equal(flat$fontsize, 11)                 # imported, over cfg's 14
    expect_identical(flat$reference_colour, "#123456")  # cfg's
    expect_setequal(names(flat)[-(1:4)], setdiff(names(.pt_default_settings()), .PT_CORE_SETTINGS))
    expect_identical(res$wiz_display_settings(), flat[-(1:4)])
  })
})

test_that("the export dialog offers samples.txt in the data folder, or a blank path with none", {
  got <- catch_modals()
  shiny::testServer(export_app("/data/proj"), session$setInputs(wiz_export_samples = 1))
  shiny::testServer(export_app(NULL), session$setInputs(wiz_export_samples = 1))
  expect_match(got$shown[[1]], 'id="wiz_export_path"[^>]*value="/data/proj/samples.txt"')
  expect_match(got$shown[[2]], 'id="wiz_export_path"[^>]*value=""')
})

test_that("Export writes a samples.txt that imports back to the same samples and settings", {
  root <- tempfile("pt-export-"); dir.create(file.path(root, "atac"), recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  file.create(file.path(root, "atac", "a.bw"))
  samples <- list(list(condition = "Ba", colour = "#1B9E77", data = list(list(type = "bw", name = "ATAC", path = "atac/a.bw"))))
  out <- file.path(root, "out.txt")
  shiny::testServer(export_app(root, samples, visual = list(fontsize = 11)), {
    expect_identical(status_html(output$wiz_export_status), "")
    session$setInputs(wiz_export_path = "  ", wiz_confirm_export = 1)
    expect_match(status_html(output$wiz_export_status), "Enter a path first.", fixed = TRUE)
    session$setInputs(wiz_gtf_path = "/x/genes.gtf.gz", wiz_orgdb = "org.Mm.eg.db", wiz_genome_label = "mm10")
    session$setInputs(wiz_export_path = out, wiz_confirm_export = 2)
    expect_match(status_html(output$wiz_export_status), sprintf("Wrote %s.", out), fixed = TRUE)
  })
  back <- import_samples_txt(out, defaults = .pt_default_settings())
  expect_length(back$errors, 0)
  expect_identical(back$data_root, root)
  expect_identical(back$samples[[1]]$condition, "Ba")
  expect_identical(back$samples[[1]]$data[[1]]$path, "atac/a.bw")
  expect_identical(back$overrides$gtf_path, "/x/genes.gtf.gz")
  expect_equal(as.numeric(back$overrides$fontsize), 11)
})

test_that("an Export that can't be written says why", {
  shiny::testServer(export_app(NULL), {
    suppressWarnings(session$setInputs(wiz_export_path = file.path(tempdir(), "no-such-folder", "s.txt"),
                                       wiz_confirm_export = 1))
    expect_match(status_html(output$wiz_export_status), "wiz-error.*Failed to write:")
  })
})
