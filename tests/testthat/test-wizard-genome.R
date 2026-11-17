# The wizard's Genome section (.pt_wizard_genome_server(), R/wizard-genome.R) on its own: the
# species preset and the live status lines under the GTF, OrgDb and genome label fields. A
# few-line GTF in a temp file, no genome data.

genome_app <- function(input, output, session) res <- .pt_wizard_genome_server(input, output, session)

status_html <- function(out) if (is.null(out)) "" else paste(as.character(out$html), collapse = "")

# A GENCODE-style mouse GTF of one gene, build GRCm39.
mouse_gtf <- function() {
  f <- tempfile(fileext = ".gtf")
  writeLines(c("##description: evidence-based annotation of the mouse genome (GRCm39), version M32 (Ensembl 109)",
               "##provider: GENCODE",
               'chr1\tHAVANA\tgene\t1\t100\t.\t+\t.\tgene_id "ENSMUSG00000000001.5"; gene_name "Gnai3";'), f)
  f
}

test_that("a status line marks ok, warn and error, and is NULL for blank", {
  expect_null(.pt_wizard_status_line(list(state = "blank", text = NULL)))
  ok <- as.character(.pt_wizard_status_line(list(state = "ok", text = "Fine.")))
  expect_match(ok, 'class="wiz-status wiz-ok"', fixed = TRUE)
  expect_match(ok, "✓ Fine.", fixed = TRUE)
  warn <- as.character(.pt_wizard_status_line(list(state = "warn", text = "Hmm.")))
  expect_match(warn, 'class="wiz-status wiz-warn"', fixed = TRUE)
  expect_match(warn, "⚠ Hmm.", fixed = TRUE)
  missing <- as.character(.pt_wizard_status_line(list(state = "missing", text = "Gone.")))
  expect_match(missing, 'class="wiz-status wiz-error"', fixed = TRUE)
  expect_match(missing, "✗ Gone.", fixed = TRUE)
})

test_that("a preset fills the OrgDb and genome label fields, never the GTF path; a blank one nothing", {
  sent <- list()
  # The build-hint placeholder update (bug 218) carries no value and is not recorded.
  local_mocked_bindings(updateTextInput = function(session, inputId, label = NULL, value = NULL, ...)
    if (!is.null(value)) sent[[length(sent) + 1]] <<- list(id = inputId, value = value))
  shiny::testServer(genome_app, {
    session$flushReact()
    session$setInputs(wiz_genome_preset = "")
    expect_length(sent, 0)
    session$setInputs(wiz_genome_preset = "mm10")
  })
  expect_identical(sent, list(list(id = "wiz_orgdb", value = "org.Mm.eg.db"),
                              list(id = "wiz_genome_label", value = "mm10")))
})

test_that("the GTF status line says what the file is, once the typing pauses", {
  gtf <- mouse_gtf(); on.exit(unlink(gtf), add = TRUE)
  shiny::testServer(genome_app, {
    session$setInputs(wiz_gtf_path = gtf)
    session$elapse(500)
    html <- status_html(output$wiz_gtf_status)
    expect_match(html, "wiz-ok", fixed = TRUE)
    expect_match(html, "Looks like a GENCODE GTF", fixed = TRUE)
    expect_identical(res$gtf_check()$species$orgdb, "org.Mm.eg.db")
    expect_identical(res$gtf_check()$build, "GRCm39")
    session$setInputs(wiz_gtf_path = file.path(tempdir(), "no-such.gtf"))
    session$elapse(500)
    expect_match(status_html(output$wiz_gtf_status), "wiz-error.*No file at this path")
    session$setInputs(wiz_gtf_path = "")
    session$elapse(500)
    expect_identical(status_html(output$wiz_gtf_status), "")
  })
})

test_that("the genome label line compares the label with the GTF's build", {
  gtf <- mouse_gtf(); on.exit(unlink(gtf), add = TRUE)
  shiny::testServer(genome_app, {
    session$setInputs(wiz_gtf_path = gtf, wiz_genome_label = "mm10")
    session$elapse(500)
    expect_match(status_html(output$wiz_genome_label_status), "The GTF is for GRCm39, whose UCSC name is mm39, not mm10.",
                 fixed = TRUE)
    session$setInputs(wiz_genome_label = "mm39")
    expect_match(status_html(output$wiz_genome_label_status), "Matches the GTF's build (GRCm39).", fixed = TRUE)
    session$setInputs(wiz_genome_label = "")
    expect_match(status_html(output$wiz_genome_label_status), "wiz-warn.*Needed")
  })
})

test_that("the OrgDb line says when the package isn't installed, or isn't the GTF's species", {
  gtf <- mouse_gtf(); on.exit(unlink(gtf), add = TRUE)
  shiny::testServer(genome_app, {
    session$setInputs(wiz_gtf_path = gtf, wiz_orgdb = "org.Xx.eg.db")
    session$elapse(500)
    expect_match(status_html(output$wiz_orgdb_status), "BiocManager::install('org.Xx.eg.db')", fixed = TRUE)
    session$setInputs(wiz_orgdb = "")
    expect_identical(status_html(output$wiz_orgdb_status), "")
  })
  skip_if_not_installed("org.Hs.eg.db")
  shiny::testServer(genome_app, {
    session$setInputs(wiz_gtf_path = gtf, wiz_orgdb = "org.Hs.eg.db")
    session$elapse(500)
    expect_match(status_html(output$wiz_orgdb_status), "the GTF's gene IDs are mouse: that species' package is org.Mm.eg.db",
                 fixed = TRUE)
  })
})
