# The setup wizard's live GTF and OrgDb checks (W29) and its help buttons. No fixture data:
# the bundled example GTF and small temp files only, so these run on R-universe too.

example_gtf <- system.file("extdata", "example", "annotation.gtf.gz", package = "PlotTrellis")

# A tiny GTF written with the given gene_id attribute values, one gene line each.
write_gtf <- function(gene_ids, header = character(0), ext = ".gtf") {
  f <- tempfile("wiz-gtf-", fileext = ext)
  body <- sprintf('chr1\tsrc\tgene\t%d\t%d\t.\t+\t.\tgene_id "%s"; gene_name "G%d";',
                  seq_along(gene_ids) * 1000, seq_along(gene_ids) * 1000 + 500, gene_ids, seq_along(gene_ids))
  writeLines(c(header, body), f)
  f
}

test_that("wiz_gtf_check() recognises the bundled GENCODE GTF and names the species", {
  skip_if(!nzchar(example_gtf))
  chk <- wiz_gtf_check(example_gtf)
  expect_identical(chk$state, "ok")
  expect_match(chk$text, "GENCODE")
  expect_match(chk$text, "GRCm39")
  expect_match(chk$text, "mouse")
  expect_identical(chk$species$orgdb, "org.Mm.eg.db")
})

test_that("wiz_gtf_check() flags a RefSeq-style GTF (symbol gene IDs) and says why", {
  f <- write_gtf(c("Nfatc1", "Coro2a"))
  chk <- wiz_gtf_check(f)
  expect_identical(chk$state, "warn")
  expect_match(chk$text, "Nfatc1", fixed = TRUE)
  expect_match(chk$text, "RefSeq", fixed = TRUE)
  expect_match(chk$text, "Ensembl IDs", fixed = TRUE)
  expect_null(chk$species)
})

test_that("wiz_gtf_check() reads plain Ensembl GTFs (bare IDs, #!genome-build header) and human IDs", {
  f <- write_gtf(c("ENSG00000141510", "ENSG00000012048"),
                 header = c("#!genome-build GRCh38.p14", "#!genome-version GRCh38"))
  chk <- wiz_gtf_check(f)
  expect_identical(chk$state, "ok")
  expect_match(chk$text, "Ensembl GTF, genome build GRCh38.p14 (human)", fixed = TRUE)
  expect_identical(chk$species$orgdb, "org.Hs.eg.db")
  expect_identical(chk$build, "GRCh38")
  # Ensembl IDs of another species: a warning (mouse and human only), no species for the OrgDb check.
  g <- write_gtf("ENSGALG00000000001")
  chk <- wiz_gtf_check(g)
  expect_identical(chk$state, "warn")
  expect_match(chk$text, "not mouse or human", fixed = TRUE)
  expect_null(chk$species)
})

test_that("wiz_gtf_check() flags GENCODE lift37 IDs, skips a non-Ensembl first line, ignores ref_gene_id", {
  lift <- write_gtf(c("ENSG00000290825.1_1", "ENSG00000223972.6_2"))
  chk <- wiz_gtf_check(lift)
  expect_identical(chk$state, "warn")
  expect_match(chk$text, "lift37", fixed = TRUE)
  expect_match(chk$text, "GENCODE 19", fixed = TRUE)
  # The species comes from the first Ensembl ID, not the first line.
  mixed <- write_gtf(c("Nfatc1", "ENSMUSG00000033016"))
  expect_identical(wiz_gtf_check(mixed)$species$orgdb, "org.Mm.eg.db")
  # StringTie-style lines: gene_id is MSTRG..., ref_gene_id the Ensembl one; the real gene_id wins.
  st <- tempfile(fileext = ".gtf")
  writeLines('chr1\tStringTie\ttranscript\t1000\t2000\t.\t+\t.\tgene_id "MSTRG.1"; transcript_id "MSTRG.1.1"; ref_gene_id "ENSG00000141510";', st)
  chk <- wiz_gtf_check(st)
  expect_identical(chk$state, "warn")
  expect_match(chk$text, "MSTRG.1", fixed = TRUE)
})

test_that("wiz_genome_label_check(): blank warns, the GTF's build must match, unknown build says nothing", {
  expect_identical(wiz_genome_label_check("")$state, "warn")
  expect_match(wiz_genome_label_check(NULL)$text, "first launch of a new GTF stops", fixed = TRUE)
  wrong <- wiz_genome_label_check("mm39", "GRCm38")
  expect_identical(wrong$state, "warn")
  expect_match(wrong$text, "GRCm38, whose UCSC name is mm10, not mm39", fixed = TRUE)
  expect_identical(wiz_genome_label_check("mm10", "GRCm38")$state, "ok")
  expect_identical(wiz_genome_label_check("hg38", NULL)$state, "blank")
  expect_identical(wiz_genome_label_check("hg38", "GRCz11")$state, "blank")
  # No build in the header: a label of the other species is still caught.
  mouse <- list(species = "mouse", orgdb = "org.Mm.eg.db")
  expect_match(wiz_genome_label_check("hg38", NULL, mouse)$text, "hg38 is a human build, but the GTF's gene IDs are mouse", fixed = TRUE)
  expect_identical(wiz_genome_label_check("mm10", NULL, mouse)$state, "blank")
})

test_that("wiz_gtf_check() takes the build from whichever header line names one", {
  f <- write_gtf("ENSG00000141510", header = c("##description: a subset", "#!genome-build GRCh38"))
  expect_identical(wiz_gtf_check(f)$build, "GRCh38")
})

test_that("the wizard's live status lines render from the inputs (debounced path)", {
  skip_if(!nzchar(example_gtf))
  json <- file.path(tempfile("wiz-live-"), "project_settings.json")
  dir.create(dirname(json))
  defaults <- .pt_default_settings()
  shiny::testServer(wizard_server(json, defaults, defaults), {
    session$setInputs(wiz_gtf_path = example_gtf, wiz_orgdb = "org.Xx.eg.db", wiz_genome_label = "mm10")
    session$elapse(500)
    expect_match(as.character(output$wiz_gtf_status$html), "GENCODE GTF", fixed = TRUE)
    expect_match(as.character(output$wiz_orgdb_status$html), "Not installed", fixed = TRUE)
    # The example GTF is GRCm39: mm10 is the wrong label for it.
    expect_match(as.character(output$wiz_genome_label_status$html), "UCSC name is mm39", fixed = TRUE)
    session$setInputs(wiz_genome_label = "mm39")
    expect_match(as.character(output$wiz_genome_label_status$html), "Matches the GTF", fixed = TRUE)
  })
})

test_that("the genome label's placeholder follows the GTF's build (bug 218)", {
  skip_if(!nzchar(example_gtf))
  hints <- character(0)
  local_mocked_bindings(updateTextInput = function(session, inputId, label = NULL, value = NULL, placeholder = NULL)
    if (identical(inputId, "wiz_genome_label") && !is.null(placeholder)) hints <<- c(hints, placeholder))
  json <- file.path(tempfile("wiz-hint-"), "project_settings.json")
  dir.create(dirname(json))
  defaults <- .pt_default_settings()
  shiny::testServer(wizard_server(json, defaults, defaults), {
    session$setInputs(wiz_gtf_path = example_gtf)
    session$elapse(500)
    expect_identical(tail(hints, 1), "mm39")  # the example GTF is GRCm39
    session$setInputs(wiz_gtf_path = write_gtf("ENSG00000000003.15", header = "#!genome-build GRCh38.p14"))
    session$elapse(500)
    expect_identical(tail(hints, 1), "hg38")
    session$setInputs(wiz_gtf_path = write_gtf("ENSMUSG00000000001.4"))  # no build in the header
    session$elapse(500)
    expect_identical(tail(hints, 1), "mm39")
    session$setInputs(wiz_gtf_path = write_gtf("ENSG00000000003.15"))
    session$elapse(500)
    expect_identical(tail(hints, 1), "hg38")  # no build: the species' current build
    session$setInputs(wiz_gtf_path = "")
    session$elapse(500)
    expect_identical(tail(hints, 1), "mm39")  # no GTF: the old default
  })
})

test_that("the GTF file is read once per pause in typing, not per keystroke", {
  skip_if(!nzchar(example_gtf))
  calls <- 0L
  real <- wiz_gtf_check
  local_mocked_bindings(wiz_gtf_check = function(path, ...) { calls <<- calls + 1L; real(path, ...) })
  json <- file.path(tempfile("wiz-deb-"), "project_settings.json")
  dir.create(dirname(json))
  defaults <- .pt_default_settings()
  shiny::testServer(wizard_server(json, defaults, defaults), {
    output$wiz_gtf_status  # the output exists before any input
    for (n in seq_len(nchar(example_gtf))[-1]) {
      session$setInputs(wiz_gtf_path = substr(example_gtf, 1, n))
      session$elapse(50)
    }
    before <- calls
    session$elapse(500)
    output$wiz_gtf_status
    expect_lte(before, 1L)
    expect_equal(calls, before + 1L)
  })
})

test_that("wiz_gtf_check(): blank, missing, wrong extension, no gene_id", {
  expect_identical(wiz_gtf_check("")$state, "blank")
  expect_identical(wiz_gtf_check(NULL)$state, "blank")
  expect_identical(wiz_gtf_check(file.path(tempdir(), "nope.gtf"))$state, "missing")
  # A folder is not a file.
  expect_identical(wiz_gtf_check(tempdir())$state, "missing")
  gff <- write_gtf("ENSG00000141510", ext = ".gff3")
  expect_identical(wiz_gtf_check(gff)$state, "warn")
  expect_match(wiz_gtf_check(gff)$text, "GFF3", fixed = TRUE)
  empty <- tempfile(fileext = ".gtf"); writeLines(c("# only a comment", ""), empty)
  chk <- wiz_gtf_check(empty)
  expect_identical(chk$state, "warn")
  expect_match(chk$text, "No gene_id", fixed = TRUE)
})

test_that("wiz_gtf_check() reads only the first lines of a large file", {
  f <- write_gtf(c(rep("ENSMUSG00000000001", 10), rep("Nfatc1", 10)))
  # The RefSeq-looking IDs sit after the first 10 lines: not seen with n_lines = 10, so OK.
  expect_identical(wiz_gtf_check(f, n_lines = 10)$state, "ok")
})

test_that("wiz_orgdb_check() says installed, not installed, or wrong species for the GTF", {
  expect_identical(wiz_orgdb_check("")$state, "blank")
  nope <- wiz_orgdb_check("org.Xx.eg.db")
  expect_identical(nope$state, "warn")
  expect_match(nope$text, "BiocManager::install('org.Xx.eg.db')", fixed = TRUE)
  skip_if_not_installed("org.Mm.eg.db")
  mouse <- list(species = "mouse", orgdb = "org.Mm.eg.db")
  expect_identical(wiz_orgdb_check("org.Mm.eg.db", mouse)$state, "ok")
  expect_match(wiz_orgdb_check("org.Mm.eg.db", mouse)$text, "matches the GTF (mouse)", fixed = TRUE)
  expect_identical(wiz_orgdb_check("org.Mm.eg.db", NULL)$state, "ok")
  human <- list(species = "human", orgdb = "org.Hs.eg.db")
  wrong <- wiz_orgdb_check("org.Mm.eg.db", human)
  expect_identical(wrong$state, "warn")
  expect_match(wrong$text, "org.Hs.eg.db", fixed = TRUE)
})

test_that("wiz_help() is a popover button carrying its text, not a title tooltip", {
  html <- as.character(wiz_help("Font size:", "Points, 4 to 24."))
  expect_match(html, "Font size:", fixed = TRUE)
  expect_match(html, '<button type="button" class="pt-help"', fixed = TRUE)
  expect_match(html, 'data-pt-help="Points, 4 to 24."', fixed = TRUE)
  expect_false(grepl("title=", html, fixed = TRUE))
  # The wizard page carries the popover script and stylesheet (inlined by includeScript/includeCSS).
  page <- as.character(wizard_ui(file.path(tempdir(), "no-such-settings.json")))
  expect_match(page, "pt-help-pinned", fixed = TRUE)
  expect_match(page, ".popover.pt-help-popover", fixed = TRUE)
  # So does the app's page (.pt_main_ui(), which builds it from no project data).
  app_page <- unlist(htmltools::renderTags(.pt_main_ui(list()))[c("head", "html")])
  app_page <- paste(app_page, collapse = "\n")
  expect_match(app_page, "pt-help-pinned", fixed = TRUE)
  expect_match(app_page, ".popover.pt-help-popover", fixed = TRUE)
})

test_that("import warnings are a note, not an alert; Save errors are an alert", {
  warn <- as.character(wiz_error_list("Imported, with warnings:", "x", kind = "warn"))
  expect_match(warn, 'role="status"', fixed = TRUE)
  expect_match(warn, "wiz-warn-box", fixed = TRUE)
  err <- as.character(wiz_error_list("Not saved yet. Please fix:", "x"))
  expect_match(err, 'role="alert"', fixed = TRUE)
})

test_that("the wizard's Save and Add-sample errors are lists, one problem per item", {
  html <- as.character(wiz_error_list("Not saved yet. Please fix:", c("Data root is required.", "GTF path is required.")))
  expect_match(html, "<li>Data root is required.</li>", fixed = TRUE)
  expect_match(html, "<li>GTF path is required.</li>", fixed = TRUE)
  expect_false(grepl(" | ", html, fixed = TRUE))
})

# W92: nobody reads more than one sentence, so a help popover is one short sentence. The Hi-C
# popovers (W93's) are left out here; the sample modal's help is covered by the same
# wiz_help() path but needs a session, so the composer's settings panel and the wizard's form
# stand for it.
test_that("help popovers are one short sentence (Hi-C settings excluded)", {
  defaults <- .pt_default_settings()
  ui <- .pt_display_settings_ui(defaults, "s_")
  html <- paste(as.character(ui$display), as.character(ui$diff),
                as.character(wiz_shared_form(defaults, is_edit = FALSE)),
                as.character(.pt_plot_title_input()))
  chunks <- strsplit(html, "<label", fixed = TRUE)[[1]]
  chunks <- chunks[grepl("data-pt-help=", chunks, fixed = TRUE)]
  expect_gt(length(chunks), 15)
  chunks <- chunks[!grepl("hic_|zrange_cap", substr(chunks, 1, 300))]
  texts <- sub('.*data-pt-help="([^"]*)".*', "\\1", chunks)
  # The colour swatches' labels carry a native title tooltip, held to the same rule.
  titles <- regmatches(html, gregexpr('<label[^>]*title="[^"]*"', html))[[1]]
  expect_gt(length(titles), 3)
  texts <- c(texts, sub('.*title="([^"]*)"$', "\\1", titles))
  expect_true(all(nchar(texts) <= 110), info = paste(texts[nchar(texts) > 110], collapse = " | "))
  # One sentence: no full stop followed by a new capitalised sentence.
  expect_false(any(grepl("\\. [A-Z]", texts)), info = paste(texts[grepl("\\. [A-Z]", texts)], collapse = " | "))
})
