# Characterisation of the setup wizard's server (wizard_server()): what each block does today,
# driven headlessly through shiny::testServer() with session$setInputs(), and read back from
# its outputs, its reactive values, the dialogs it shows and the files it writes. These pin
# the current behaviour (oddities included, marked "current behaviour, possibly a bug") so a
# restructuring of the server can be checked to change nothing. No fixture data needed.
#
# Only inputs, outputs and the server's reactive values are used here, never its internal
# helper functions (those are unit-tested in test-wizard-helpers.R).
#
# testServer() evaluates its block with the server's own variables in front of the test's, so
# test variables here avoid the server's names (data_root, cfg, is_edit, ...).

# A project folder with a data folder holding one file of each kind, two salmon folders, and a
# bigWig outside the data folder.
ch_project <- function() {
  proj <- tempfile("pt-wizch-")
  for (d in c("atac", "hic", "peaks", "junc", "quant/S1", "quant/S2"))
    dir.create(file.path(proj, "data", d), recursive = TRUE)
  file.create(file.path(proj, "data", c("atac/a.bw", "atac/b.bw", "hic/m.hic", "quant/S1/quant.sf", "quant/S2/quant.sf")))
  writeLines("chr1\t100\t200", file.path(proj, "data", "peaks", "p.bed"))
  writeLines(paste(c("chr1", "100", "900", "j1", "5", "+", "100", "900", "0", "2", "10,10", "0,790"), collapse = "\t"),
             file.path(proj, "data", "junc", "j.bed"))
  file.create(file.path(proj, "outside.bw"))
  proj
}
ch_json <- function(proj) file.path(proj, "project_settings.json")
# The wizard's server for `proj`, built as plottrellis_wizard() builds it.
ch_server <- function(proj) {
  defaults <- .pt_default_settings()
  wizard_server(ch_json(proj), .pt_read_settings(ch_json(proj), defaults), defaults)
}
# An existing project (edit mode): a settings file with two samples.
ch_edit_project <- function(extra = list()) {
  proj <- ch_project()
  s <- utils::modifyList(list(
    data_dir = file.path(proj, "data"), gtf_path = "/x/genes.gtf.gz", orgdb = "org.Mm.eg.db", genome_label = "mm10",
    samples = list(
      list(condition = "Ba", colour = "#1B9E77",
           data = list(list(type = "bw", name = "ATAC", path = "atac/a.bw"),
                       list(type = "quantdir", name = "Salmon", path = "quant/S1"))),
      list(condition = "Lp", colour = "#D95F02",
           data = list(list(type = "bw", name = "ATAC", path = "atac/b.bw"))))), extra)
  jsonlite::write_json(s, ch_json(proj), auto_unbox = TRUE, null = "null", pretty = TRUE)
  proj
}
# An output's HTML as one string ("" for a NULL render).
ch_html <- function(out) if (is.null(out)) "" else paste(as.character(out$html), collapse = "")
ch_conds <- function(samples) vapply(samples, function(s) s$condition, "")
# Records the wizard's dialogs, notifications, input updates and its scheduled stop, instead of
# sending them to a browser. The stop runs only when the test calls rec$later[[i]]$func().
ch_mocks <- function(env = parent.frame()) {
  rec <- new.env()
  rec$modals <- character(0); rec$removed <- 0L; rec$notes <- character(0)
  rec$updates <- list(); rec$later <- list(); rec$stopped <- list()
  local_mocked_bindings(
    showModal = function(ui, ...) rec$modals <- c(rec$modals, paste(as.character(ui), collapse = "")),
    removeModal = function(...) rec$removed <- rec$removed + 1L,
    showNotification = function(ui, ..., type = "default") rec$notes <- c(rec$notes, paste0(type, ": ", ui)),
    # A placeholder-only update (the genome label's build hint, bug 218) is not recorded.
    updateTextInput = function(session, inputId, label = NULL, value = NULL, placeholder = NULL)
      if (is.null(placeholder) || !is.null(value))
        rec$updates[[length(rec$updates) + 1]] <- list(id = inputId, value = value),
    stopApp = function(returnValue = invisible()) rec$stopped[[length(rec$stopped) + 1]] <- returnValue,
    .env = env)
  local_mocked_bindings(later = function(func, delay = 0, ...) {
    rec$later[[length(rec$later) + 1]] <- list(func = func, delay = delay)
    invisible(NULL)
  }, .package = "later", .env = env)
  rec
}
ch_last_modal <- function(rec) if (length(rec$modals)) rec$modals[[length(rec$modals)]] else ""
ch_updates <- function(rec) setNames(lapply(rec$updates, `[[`, "value"), vapply(rec$updates, `[[`, "", "id"))
# A samples.txt written to a temp file (removed when the calling test ends).
ch_txt <- function(lines, env = parent.frame()) {
  f <- tempfile("pt-wizch-", fileext = ".txt")
  writeLines(lines, f)
  withr::defer(unlink(f), envir = env)
  f
}
# Sets the genome and data folder inputs to values that pass Save's checks.
ch_good_genome <- function(session, proj)
  session$setInputs(wiz_data_root = file.path(proj, "data"), wiz_gtf_path = "/x/genes.gtf.gz",
                    wiz_orgdb = "org.Mm.eg.db", wiz_genome_label = "mm10")
# A tiny GTF: one gene line per id, with optional header lines; gzipped for a .gz name.
ch_gtf <- function(ids, header = character(0), ext = ".gtf", attr = 'gene_id "%s"; transcript_id "%s.1";', env = parent.frame()) {
  f <- tempfile("pt-wizch-gtf-", fileext = ext)
  body <- sprintf(paste0("chr1\tsrc\tgene\t%d\t%d\t.\t+\t.\t", attr), seq_along(ids) * 100, seq_along(ids) * 100 + 50, ids, ids)
  con <- if (grepl("\\.gz$", ext)) gzfile(f, "w") else file(f, "w")
  writeLines(c(header, body), con)
  close(con)
  withr::defer(unlink(f), envir = env)
  f
}

# ---- Mode switch and the page body -------------------------------------------------------

test_that("a first run with no settings and no copy opens on the landing choice", {
  proj <- ch_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  shiny::testServer(ch_server(proj), {
    expect_null(wiz_mode())
    expect_identical(wiz_samples(), list())
    html <- ch_html(output$wiz_body)
    expect_match(html, "What a project needs", fixed = TRUE)
    expect_match(html, 'id="wiz_mode_manual"', fixed = TRUE)
    expect_match(html, 'id="wiz_mode_import"', fixed = TRUE)
    expect_no_match(html, 'id="wiz_gtf_path"', fixed = TRUE)
    expect_no_match(html, "wiz_mode_back", fixed = TRUE)
  })
})

test_that("Set up by hand shows the form with a Back link and no import box", {
  proj <- ch_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  shiny::testServer(ch_server(proj), {
    session$setInputs(wiz_mode_manual = 1)
    expect_identical(wiz_mode(), "manual")
    html <- ch_html(output$wiz_body)
    expect_match(html, 'id="wiz_mode_back"', fixed = TRUE)
    expect_match(html, "Back to the start", fixed = TRUE)
    expect_no_match(html, 'id="wiz_import_path"', fixed = TRUE)
    expect_match(html, 'id="wiz_gtf_path"', fixed = TRUE)
    expect_match(html, "Save and open the app", fixed = TRUE)
    expect_no_match(html, 'id="wiz_close"', fixed = TRUE)
    # A first run's data folder starts blank, whatever the defaults say.
    expect_match(html, '<input id="wiz_data_root" type="text" class="shiny-input-text form-control" value=""', fixed = TRUE)
  })
})

test_that("Import a samples.txt shows the import box, blank, above the same form", {
  proj <- ch_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  shiny::testServer(ch_server(proj), {
    session$setInputs(wiz_mode_import = 1)
    expect_identical(wiz_mode(), "import")
    html <- ch_html(output$wiz_body)
    expect_match(html, 'id="wiz_mode_back"', fixed = TRUE)
    expect_match(html, '<input id="wiz_import_path" type="text" class="shiny-input-text form-control" value=""', fixed = TRUE)
    expect_match(html, "The file fills in the three sections below", fixed = TRUE)
    expect_match(html, 'id="wiz_import_load"', fixed = TRUE)
    expect_match(html, 'id="wiz_gtf_path"', fixed = TRUE)
  })
})

test_that("Back returns to the landing choice and clears Save's errors, but keeps the samples", {
  proj <- ch_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  rec <- ch_mocks()
  txt <- ch_txt(c("condition\tcolour\tbw:ATAC", "Ba\t#1B9E77\tatac/a.bw", "# fontsize: 11"))
  shiny::testServer(ch_server(proj), {
    session$setInputs(wiz_mode_import = 1)
    session$setInputs(wiz_import_path = txt, wiz_import_load = 1)
    session$setInputs(wiz_save = 1)
    expect_true(save_failed())
    session$setInputs(wiz_mode_back = 1)
    expect_null(wiz_mode())
    expect_false(save_failed())
    expect_null(write_error())
    expect_identical(ch_html(output$wiz_completion_error), "")
    # current behaviour, possibly a bug: Back keeps the imported samples and display settings,
    # so a manual setup started after it saves them.
    expect_identical(ch_conds(wiz_samples()), "Ba")
    expect_identical(import_visual(), list(fontsize = 11))
    session$setInputs(wiz_mode_manual = 1)
    expect_identical(wiz_mode(), "manual")
    expect_match(ch_html(output$wiz_samples_ui), "<b>Ba</b>", fixed = TRUE)
  })
})

test_that("an existing project opens in edit mode on its own settings, with no landing or Back", {
  proj <- ch_edit_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  shiny::testServer(ch_server(proj), {
    expect_identical(wiz_mode(), "edit")
    expect_identical(ch_conds(wiz_samples()), c("Ba", "Lp"))
    html <- ch_html(output$wiz_body)
    expect_no_match(html, "What a project needs", fixed = TRUE)
    expect_no_match(html, "wiz_mode_back", fixed = TRUE)
    expect_no_match(html, 'id="wiz_import_path"', fixed = TRUE)
    expect_match(html, "Save and return to the app", fixed = TRUE)
    expect_match(html, 'id="wiz_close"', fixed = TRUE)
    expect_match(html, sprintf('id="wiz_data_root" type="text" class="shiny-input-text form-control" value="%s"',
                               file.path(proj, "data")), fixed = TRUE)
    expect_match(html, 'id="wiz_gtf_path" type="text" class="shiny-input-text form-control" value="/x/genes.gtf.gz"', fixed = TRUE)
  })
})

test_that("a first run beside an earlier project_samples.txt opens on Import with that file filled in", {
  proj <- ch_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  mirror <- file.path(proj, "project_samples.txt")
  writeLines(c("condition\tcolour\tbw:ATAC", "Ba\t#1B9E77\tatac/a.bw"), mirror)
  shiny::testServer(ch_server(proj), {
    expect_identical(wiz_mode(), "import")
    html <- ch_html(output$wiz_body)
    expect_match(html, "A project_samples.txt saved with your last project was found", fixed = TRUE)
    expect_match(html, sprintf('id="wiz_import_path" type="text" class="shiny-input-text form-control" value="%s"', mirror), fixed = TRUE)
    # Nothing is imported until Import is clicked.
    expect_identical(wiz_samples(), list())
    # Back, then Import again: still prefilled with the copy.
    session$setInputs(wiz_mode_back = 1)
    expect_null(wiz_mode())
    session$setInputs(wiz_mode_import = 1)
    expect_match(ch_html(output$wiz_body), sprintf('value="%s"', mirror), fixed = TRUE)
  })
})

test_that("an existing project ignores a project_samples.txt beside it", {
  proj <- ch_edit_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  writeLines(c("condition\tcolour\tbw:ATAC", "Zz\t#1B9E77\tatac/a.bw"), file.path(proj, "project_samples.txt"))
  shiny::testServer(ch_server(proj), {
    expect_identical(wiz_mode(), "edit")
    expect_identical(ch_conds(wiz_samples()), c("Ba", "Lp"))
  })
})

test_that("a folder named project_samples.txt counts as the saved copy", {
  proj <- ch_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  dir.create(file.path(proj, "project_samples.txt"))
  shiny::testServer(ch_server(proj), {
    # current behaviour, possibly a bug: file.exists() is TRUE for a folder, so the wizard
    # opens on Import prefilled with a path Import then rejects.
    expect_identical(wiz_mode(), "import")
    session$setInputs(wiz_import_load = 1, wiz_import_path = file.path(proj, "project_samples.txt"))
    expect_match(ch_html(output$wiz_import_status), "No such file:", fixed = TRUE)
  })
})

# ---- Genome preset -----------------------------------------------------------------------

test_that("a genome preset fills the OrgDb and genome label, never the GTF path", {
  proj <- ch_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  rec <- ch_mocks()
  shiny::testServer(ch_server(proj), {
    session$flushReact()
    session$setInputs(wiz_genome_preset = "hg38")
    expect_identical(ch_updates(rec), list(wiz_orgdb = "org.Hs.eg.db", wiz_genome_label = "hg38"))
    session$setInputs(wiz_genome_preset = "mm10")
    expect_identical(ch_updates(rec)[3:4], list(wiz_orgdb = "org.Mm.eg.db", wiz_genome_label = "mm10"))
    # The blank choice and an unknown id change nothing.
    session$setInputs(wiz_genome_preset = "")
    session$setInputs(wiz_genome_preset = "dm6")
    expect_length(rec$updates, 4)
  })
})

test_that("a genome preset already set when the wizard starts is not applied", {
  proj <- ch_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  rec <- ch_mocks()
  shiny::testServer(ch_server(proj), {
    session$setInputs(wiz_genome_preset = "hg38")  # the first flush: ignored as the initial value
    expect_length(rec$updates, 0)
  })
})

# ---- Data folder -------------------------------------------------------------------------

test_that("the data folder: blank or spaces is unset, with no status line", {
  proj <- ch_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  shiny::testServer(ch_server(proj), {
    expect_null(data_root())
    expect_identical(ch_html(output$wiz_data_root_status), "")
    session$setInputs(wiz_data_root = "   ")
    expect_null(data_root())
    expect_identical(ch_html(output$wiz_data_root_status), "")
  })
})

test_that("the data folder: an existing one is found, trimmed, a trailing slash kept", {
  proj <- ch_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  dr <- file.path(proj, "data")
  shiny::testServer(ch_server(proj), {
    session$setInputs(wiz_data_root = paste0("  ", dr, "  "))
    expect_identical(data_root(), dr)
    expect_identical(ch_html(output$wiz_data_root_status), '<div class="wiz-status wiz-ok">\u2713 Folder found.</div>')
    session$setInputs(wiz_data_root = paste0(dr, "/"))
    expect_identical(data_root(), paste0(dr, "/"))
    expect_match(ch_html(output$wiz_data_root_status), "Folder found.", fixed = TRUE)
  })
})

test_that("the data folder: a missing one says so; a relative one is relative to R's working folder", {
  proj <- ch_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  shiny::testServer(ch_server(proj), {
    session$setInputs(wiz_data_root = file.path(proj, "nope"))
    expect_identical(ch_html(output$wiz_data_root_status),
                     '<div class="wiz-status wiz-error">\u2717 No folder at this path (it must be a path on the computer running R).</div>')
    withr::with_dir(proj, {
      session$setInputs(wiz_data_root = "data")
      expect_identical(data_root(), "data")  # kept as typed, not made absolute
      expect_match(ch_html(output$wiz_data_root_status), "Folder found.", fixed = TRUE)
    })
  })
})

# ---- Live GTF, OrgDb and genome label checks ---------------------------------------------

test_that("before anything is typed: no GTF or OrgDb line, the genome label warns", {
  proj <- ch_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  shiny::testServer(ch_server(proj), {
    expect_identical(gtf_check()$state, "blank")
    expect_identical(ch_html(output$wiz_gtf_status), "")
    expect_identical(ch_html(output$wiz_orgdb_status), "")
    expect_match(ch_html(output$wiz_genome_label_status), '<div class="wiz-status wiz-warn">\u26a0 Needed:', fixed = TRUE)
  })
})

test_that("a GTF path that isn't a file says so in red, after the pause in typing", {
  proj <- ch_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  shiny::testServer(ch_server(proj), {
    session$flushReact()
    session$setInputs(wiz_gtf_path = file.path(proj, "nope.gtf"))
    expect_identical(gtf_check()$state, "blank")  # not read yet: the path is debounced
    session$elapse(300)
    expect_identical(gtf_check()$state, "blank")
    session$elapse(500)
    expect_identical(gtf_check()$state, "missing")
    expect_identical(ch_html(output$wiz_gtf_status),
                     '<div class="wiz-status wiz-error">\u2717 No file at this path (it must be a path on the computer running R).</div>')
    # A folder is not a file either.
    session$setInputs(wiz_gtf_path = proj)
    session$elapse(500)
    expect_identical(gtf_check()$state, "missing")
  })
})

test_that("a file that isn't .gtf or .gtf.gz warns", {
  proj <- ch_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  f <- ch_gtf("ENSMUSG00000000001", ext = ".txt")
  shiny::testServer(ch_server(proj), {
    session$setInputs(wiz_gtf_path = f)
    session$elapse(500)
    expect_identical(ch_html(output$wiz_gtf_status), paste0(
      '<div class="wiz-status wiz-warn">\u26a0 Expected a .gtf or .gtf.gz file. ',
      'A GFF3 (.gff, .gff3) or other annotation format is not read.</div>'))
    expect_null(gtf_check()$species)
  })
})

test_that("a gzipped GENCODE mouse GTF: species and build found, OrgDb and label checked against them", {
  skip_if_not_installed("org.Mm.eg.db")
  proj <- ch_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  f <- ch_gtf(c("ENSMUSG00000000001.4", "ENSMUSG00000000003.15"), ext = ".gtf.gz",
              header = c("##description: evidence-based annotation of the mouse genome (GRCm38), version M25 (Ensembl 100)",
                         "##provider: GENCODE"))
  shiny::testServer(ch_server(proj), {
    session$setInputs(wiz_gtf_path = f, wiz_orgdb = "org.Mm.eg.db", wiz_genome_label = "mm10")
    session$elapse(500)
    expect_identical(gtf_check()$state, "ok")
    expect_identical(gtf_check()$build, "GRCm38")
    expect_identical(gtf_check()$species, list(species = "mouse", orgdb = "org.Mm.eg.db"))
    expect_identical(ch_html(output$wiz_gtf_status), paste0(
      '<div class="wiz-status wiz-ok">\u2713 Looks like a GENCODE GTF: mouse genome (GRCm38), version M25 (Ensembl 100). ',
      'Make sure it is the release your reads were quantified with.</div>'))
    expect_identical(ch_html(output$wiz_orgdb_status),
                     '<div class="wiz-status wiz-ok">\u2713 Installed, matches the GTF (mouse).</div>')
    expect_identical(ch_html(output$wiz_genome_label_status),
                     "<div class=\"wiz-status wiz-ok\">\u2713 Matches the GTF's build (GRCm38).</div>")
    # The label of another build of the same species, and the human OrgDb.
    session$setInputs(wiz_genome_label = "mm39", wiz_orgdb = "org.Xx.eg.db")
    expect_match(ch_html(output$wiz_genome_label_status), "The GTF is for GRCm38, whose UCSC name is mm10, not mm39.", fixed = TRUE)
    expect_match(ch_html(output$wiz_orgdb_status), "wiz-warn\">\u26a0 Not installed in this R.", fixed = TRUE)
  })
})

test_that("a plain human Ensembl GTF with the mouse OrgDb and no build line in the label's favour", {
  skip_if_not_installed("org.Mm.eg.db")
  proj <- ch_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  f <- ch_gtf(c("ENSG00000141510", "ENSG00000012048"), header = "#!genome-build GRCh38.p14")
  shiny::testServer(ch_server(proj), {
    session$setInputs(wiz_gtf_path = f, wiz_orgdb = "org.Mm.eg.db", wiz_genome_label = "hg19")
    session$elapse(500)
    expect_match(ch_html(output$wiz_gtf_status), "Looks like an Ensembl GTF, genome build GRCh38.p14 (human).", fixed = TRUE)
    expect_identical(ch_html(output$wiz_orgdb_status), paste0(
      "<div class=\"wiz-status wiz-warn\">\u26a0 Installed, but the GTF's gene IDs are human: ",
      "that species' package is org.Hs.eg.db.</div>"))
    expect_match(ch_html(output$wiz_genome_label_status), "The GTF is for GRCh38, whose UCSC name is hg38, not hg19.", fixed = TRUE)
  })
})

test_that("a GTF with no build line: the label is checked against the gene IDs' species only", {
  proj <- ch_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  f <- ch_gtf("ENSMUSG00000000001")
  shiny::testServer(ch_server(proj), {
    session$setInputs(wiz_gtf_path = f, wiz_genome_label = "hg38")
    session$elapse(500)
    expect_null(gtf_check()$build)
    expect_match(ch_html(output$wiz_gtf_status), "Looks like a GTF with Ensembl gene IDs (mouse).", fixed = TRUE)
    expect_match(ch_html(output$wiz_genome_label_status), "hg38 is a human build, but the GTF's gene IDs are mouse.", fixed = TRUE)
    session$setInputs(wiz_genome_label = "mm39")
    expect_identical(ch_html(output$wiz_genome_label_status), "")
  })
})

test_that("a GTF with transcript_id but no gene_id warns, and the OrgDb line has no species", {
  skip_if_not_installed("org.Mm.eg.db")
  proj <- ch_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  f <- ch_gtf(c("ENSMUSG00000000001", "ENSMUSG00000000002"), attr = 'transcript_id "%s"; gene_name "%s";')
  shiny::testServer(ch_server(proj), {
    session$setInputs(wiz_gtf_path = f, wiz_orgdb = "org.Mm.eg.db")
    session$elapse(500)
    expect_match(ch_html(output$wiz_gtf_status), "\u26a0 No gene_id attribute in the first 2000 lines. Is this a GTF file?", fixed = TRUE)
    expect_identical(ch_html(output$wiz_orgdb_status), '<div class="wiz-status wiz-ok">\u2713 Installed.</div>')
  })
})

test_that("a GTF with gene_id but no transcript_id still passes the live check", {
  proj <- ch_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  f <- ch_gtf("ENSMUSG00000000001", attr = 'gene_id "%s"; gene_name "%s";')
  shiny::testServer(ch_server(proj), {
    session$setInputs(wiz_gtf_path = f)
    session$elapse(500)
    expect_identical(gtf_check()$state, "ok")
  })
})

# ---- Add-sample dialog -------------------------------------------------------------------

test_that("+ Add a sample opens a blank Add dialog with the default colour and track name", {
  proj <- ch_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  rec <- ch_mocks()
  shiny::testServer(ch_server(proj), {
    session$setInputs(wiz_data_root = file.path(proj, "data"))
    session$setInputs(wiz_add_sample = 1)
    expect_length(rec$modals, 1)
    m <- ch_last_modal(rec)
    expect_match(m, "Add a sample", fixed = TRUE)
    expect_match(m, '<span class="action-label">Add sample</span>', fixed = TRUE)
    expect_match(m, 'id="wiz_sample_name" type="text" class="shiny-input-text form-control" value=""', fixed = TRUE)
    expect_match(m, 'id="wiz_sample_colour" value="#1f78b4"', fixed = TRUE)
    expect_match(m, "Shiny.setInputValue('wiz_sample_colour', \"#1f78b4\");", fixed = TRUE)
    expect_match(m, 'id="wiz_quantdir_name" type="text" class="shiny-input-text form-control" value="Quant"', fixed = TRUE)
    # The salmon folders under the data folder are offered, after a blank first choice.
    expect_match(m, '<option value="" selected></option>\n<option value="quant/S1">quant/S1</option>\n<option value="quant/S2">quant/S2</option>', fixed = TRUE)
    expect_null(modal_edit_name())
    expect_identical(modal_instances(), list())
    expect_identical(modal_field_errors(), list())
    expect_identical(ch_html(output$wiz_instance_rows),
                     '<p class="wiz-empty">No files yet. Use the buttons above to add one.</p>')
  })
})

test_that("the file-kind buttons add rows with ids that only go up; a row's picker lists its kind's files", {
  proj <- ch_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  rec <- ch_mocks()
  shiny::testServer(ch_server(proj), {
    session$setInputs(wiz_data_root = file.path(proj, "data"))
    session$setInputs(wiz_add_sample = 1)
    session$setInputs(wiz_add_bw = 1)
    session$setInputs(wiz_add_hic = 1)
    session$setInputs(wiz_add_bed = 1)
    session$setInputs(wiz_add_paired_end = 1)
    expect_identical(modal_instances(), list(list(id = 1, type = "bw"), list(id = 2, type = "hic"),
                                             list(id = 3, type = "bed"), list(id = 4, type = "paired-end")))
    html <- ch_html(output$wiz_instance_rows)
    for (lab in c("Coverage (bigWig)", "Hi-C map", "Intervals (BED)", "Arcs (junctions, loops)"))
      expect_match(html, sprintf('<div class="wiz-kind">%s</div>', lab), fixed = TRUE)
    for (i in 1:4) expect_match(html, sprintf('id="wiz_inst_path_%d"', i), fixed = TRUE)
    expect_match(html, '<select class="shiny-input-select form-control" id="wiz_inst_path_1"><option value="" selected></option>\n<option value="atac/a.bw">atac/a.bw</option>\n<option value="atac/b.bw">atac/b.bw</option></select>', fixed = TRUE)
    expect_match(html, '<option value="hic/m.hic">hic/m.hic</option>', fixed = TRUE)
    expect_match(html, "Shiny.setInputValue(&#39;wiz_remove_row&#39;, 4, {priority:&#39;event&#39;})", fixed = TRUE)
    # Cancel and open again: no rows, and the next row gets a new id, never a used one.
    session$setInputs(wiz_add_sample = 2)
    expect_identical(modal_instances(), list())
    session$setInputs(wiz_add_bw = 2)
    expect_identical(modal_instances(), list(list(id = 5, type = "bw")))
  })
})

test_that("a row's remove button takes out that row only; an unknown id changes nothing", {
  proj <- ch_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  rec <- ch_mocks()
  shiny::testServer(ch_server(proj), {
    session$setInputs(wiz_add_sample = 1)
    session$setInputs(wiz_add_bw = 1)
    session$setInputs(wiz_add_bed = 1)
    session$setInputs(wiz_add_bw = 2)
    session$setInputs(wiz_remove_row = 2)
    expect_identical(modal_instances(), list(list(id = 1, type = "bw"), list(id = 3, type = "bw")))
    session$setInputs(wiz_remove_row = 99)
    expect_identical(vapply(modal_instances(), `[[`, 0, "id"), c(1, 3))
  })
})

test_that("Add sample with a blank condition shows the dialog again with the error", {
  proj <- ch_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  rec <- ch_mocks()
  shiny::testServer(ch_server(proj), {
    session$setInputs(wiz_data_root = file.path(proj, "data"))
    session$setInputs(wiz_add_sample = 1)
    session$setInputs(wiz_sample_name = "   ", wiz_confirm_sample = 1)
    expect_identical(wiz_samples(), list())
    expect_length(rec$modals, 2)
    expect_identical(rec$removed, 0L)
    m <- ch_last_modal(rec)
    expect_match(m, "<b>Please fix before adding:</b>", fixed = TRUE)
    expect_match(m, "<li>Condition is required.</li>", fixed = TRUE)
    expect_match(m, 'role="alert"', fixed = TRUE)
  })
})

test_that("Add sample with a condition already used is refused; the typed values stay in the dialog", {
  proj <- ch_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  rec <- ch_mocks()
  shiny::testServer(ch_server(proj), {
    session$setInputs(wiz_data_root = file.path(proj, "data"))
    wiz_samples(list(list(condition = "Ba", colour = "#000000", data = list())))
    session$setInputs(wiz_add_sample = 1)
    session$setInputs(wiz_sample_name = " Ba ", wiz_sample_colour = "#ABCDEF", wiz_quantdir_name = "Sal",
                      wiz_confirm_sample = 1)
    expect_identical(ch_conds(wiz_samples()), "Ba")
    m <- ch_last_modal(rec)
    expect_match(m, "<li>A condition named 'Ba' already exists.</li>", fixed = TRUE)
    # A re-shown dialog keeps what was typed (the name untrimmed), not the defaults.
    expect_match(m, 'id="wiz_sample_name" type="text" class="shiny-input-text form-control" value=" Ba "', fixed = TRUE)
    expect_match(m, 'id="wiz_sample_colour" value="#ABCDEF"', fixed = TRUE)
    expect_match(m, 'id="wiz_quantdir_name" type="text" class="shiny-input-text form-control" value="Sal"', fixed = TRUE)
    # A fresh open goes back to the defaults although the inputs still hold the typed values.
    session$setInputs(wiz_add_sample = 2)
    m <- ch_last_modal(rec)
    expect_match(m, 'id="wiz_sample_name" type="text" class="shiny-input-text form-control" value=""', fixed = TRUE)
    expect_match(m, 'id="wiz_sample_colour" value="#1f78b4"', fixed = TRUE)
  })
})

test_that("a file row: a path with no name, a bare middle dot, a duplicate name each mark their row", {
  proj <- ch_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  rec <- ch_mocks()
  shiny::testServer(ch_server(proj), {
    session$setInputs(wiz_data_root = file.path(proj, "data"))
    session$setInputs(wiz_add_sample = 1)
    for (i in 1:4) session$setInputs(wiz_add_bw = i)
    session$setInputs(wiz_add_bed = 1)
    session$setInputs(
      wiz_inst_name_1 = "", wiz_inst_path_1 = "atac/a.bw",
      wiz_inst_name_2 = "ATAC\u00b7x", wiz_inst_path_2 = "atac/a.bw",
      wiz_inst_name_3 = "ATAC", wiz_inst_path_3 = "atac/a.bw",
      wiz_inst_name_4 = "ATAC", wiz_inst_path_4 = "atac/b.bw",
      wiz_inst_name_5 = "ATAC", wiz_inst_path_5 = "peaks/p.bed",  # same name, another kind: fine
      wiz_sample_name = "A", wiz_confirm_sample = 1)
    expect_identical(wiz_samples(), list())
    m <- ch_last_modal(rec)
    expect_match(m, "<li>A bw instance has a path but no name.</li>", fixed = TRUE)
    expect_match(m, "<li>Instance name 'ATAC\u00b7x' contains a bare '\u00b7' (surround it with spaces, or remove it). This breaks track-label formatting.</li>", fixed = TRUE)
    expect_match(m, "<li>Duplicate bw instance name 'ATAC'.</li>", fixed = TRUE)
    expect_no_match(m, "bed instance", fixed = TRUE)
    rows <- ch_html(output$wiz_instance_rows)
    expect_identical(lengths(regmatches(rows, gregexpr("wiz-row wiz-row-error", rows, fixed = TRUE))), 3L)
    expect_match(rows, '<div class="wiz-error" style="font-size:13px; margin-top:4px;">A bw instance has a path but no name.</div>', fixed = TRUE)
    # The rows re-render from the typed values.
    expect_match(rows, 'id="wiz_inst_name_4" type="text" class="shiny-input-text form-control" value="ATAC"', fixed = TRUE)
    expect_match(rows, '<option value="atac/b.bw" selected>atac/b.bw</option>', fixed = TRUE)
  })
})

test_that("a middle dot with spaces around it is fine; a row with no path is dropped silently", {
  proj <- ch_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  rec <- ch_mocks()
  shiny::testServer(ch_server(proj), {
    session$setInputs(wiz_data_root = file.path(proj, "data"))
    session$setInputs(wiz_add_sample = 1)
    session$setInputs(wiz_add_bw = 1)
    session$setInputs(wiz_add_bw = 2)
    session$setInputs(wiz_inst_name_1 = "ATAC \u00b7 rep1", wiz_inst_path_1 = "atac/a.bw",
                      wiz_inst_name_2 = "Unused", wiz_inst_path_2 = "  ",
                      wiz_sample_name = "A", wiz_confirm_sample = 1)
    expect_identical(wiz_samples(), list(list(condition = "A", colour = "#1f78b4",
      data = list(list(type = "bw", name = "ATAC \u00b7 rep1", path = "atac/a.bw")))))
  })
})

test_that("a file path that isn't a file under the data folder is refused", {
  proj <- ch_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  rec <- ch_mocks()
  shiny::testServer(ch_server(proj), {
    session$setInputs(wiz_data_root = file.path(proj, "data"))
    session$setInputs(wiz_add_sample = 1)
    session$setInputs(wiz_add_bw = 1); session$setInputs(wiz_add_bw = 2); session$setInputs(wiz_add_bw = 3)
    session$setInputs(wiz_inst_name_1 = "A1", wiz_inst_path_1 = "atac/none.bw",
                      wiz_inst_name_2 = "A2", wiz_inst_path_2 = "atac",               # a folder
                      wiz_inst_name_3 = "A3", wiz_inst_path_3 = file.path(proj, "data", "atac", "a.bw"),  # absolute
                      wiz_sample_name = "A", wiz_confirm_sample = 1)
    expect_identical(wiz_samples(), list())
    m <- ch_last_modal(rec)
    expect_match(m, "<li>Path for bw:A1 does not resolve to a file under the data root: atac/none.bw</li>", fixed = TRUE)
    expect_match(m, "<li>Path for bw:A2 does not resolve to a file under the data root: atac</li>", fixed = TRUE)
    # An absolute path is joined onto the data folder too, so even one inside it is refused.
    expect_match(m, "Path for bw:A3 does not resolve to a file under the data root", fixed = TRUE)
  })
})

test_that("a relative path climbing out of the data folder to a real file is accepted", {
  proj <- ch_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  rec <- ch_mocks()
  shiny::testServer(ch_server(proj), {
    session$setInputs(wiz_data_root = file.path(proj, "data"))
    session$setInputs(wiz_add_sample = 1)
    session$setInputs(wiz_add_bw = 1)
    session$setInputs(wiz_inst_name_1 = "ATAC", wiz_inst_path_1 = "../outside.bw",
                      wiz_sample_name = "A", wiz_confirm_sample = 1)
    # current behaviour, possibly a bug: "../" leaves the data folder, yet the path is accepted
    # because file.path(data folder, path) is a real file.
    expect_identical(wiz_samples()[[1]]$data[[1]]$path, "../outside.bw")
  })
})

test_that("with no data folder set, file and salmon paths are not checked at all", {
  proj <- ch_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  rec <- ch_mocks()
  shiny::testServer(ch_server(proj), {
    session$setInputs(wiz_add_sample = 1)
    session$setInputs(wiz_add_paired_end = 1)
    session$setInputs(wiz_inst_name_1 = "J", wiz_inst_path_1 = "no/such.bed",
                      wiz_quantdir_path = "no/such", wiz_quantdir_name = "",
                      wiz_sample_name = "A", wiz_confirm_sample = 1)
    expect_identical(wiz_samples(), list(list(condition = "A", colour = "#1f78b4", data = list(
      list(type = "paired-end", name = "J", path = "no/such.bed"),
      list(type = "quantdir", name = "Quant", path = "no/such")))))
    expect_identical(rec$removed, 1L)
  })
})

test_that("an arcs row needs a BED12, BED6/BEDPE or HiCCUPS file", {
  proj <- ch_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  rec <- ch_mocks()
  shiny::testServer(ch_server(proj), {
    session$setInputs(wiz_data_root = file.path(proj, "data"))
    session$setInputs(wiz_add_sample = 1)
    session$setInputs(wiz_add_paired_end = 1)
    session$setInputs(wiz_inst_name_1 = "J", wiz_inst_path_1 = "peaks/p.bed", wiz_sample_name = "A", wiz_confirm_sample = 1)
    expect_match(ch_last_modal(rec),
                 "<li>Paired-end file for J (peaks/p.bed) has 3 columns (arcs need 6 or more); if it is BED intervals, add it as Intervals (BED) instead.</li>",
                 fixed = TRUE)
    session$setInputs(wiz_inst_path_1 = "junc/j.bed", wiz_confirm_sample = 2)
    expect_identical(wiz_samples()[[1]]$data, list(list(type = "paired-end", name = "J", path = "junc/j.bed")))
  })
})

test_that("a salmon folder with no quant.sf marks the salmon row; a blank track name becomes Quant", {
  proj <- ch_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  rec <- ch_mocks()
  shiny::testServer(ch_server(proj), {
    session$setInputs(wiz_data_root = file.path(proj, "data"))
    session$setInputs(wiz_add_sample = 1)
    session$setInputs(wiz_quantdir_path = "atac", wiz_quantdir_name = " ", wiz_sample_name = "A", wiz_confirm_sample = 1)
    expect_identical(modal_field_errors(), list(quantdir = "No quant.sf found under: atac"))
    m <- ch_last_modal(rec)
    expect_match(m, '<div class="row wiz-row wiz-row-error">', fixed = TRUE)
    expect_match(m, "No quant.sf found under: atac</div>", fixed = TRUE)
    # Fixed: accepted, and the next fresh open has no field error.
    session$setInputs(wiz_quantdir_path = " quant/S2 ", wiz_confirm_sample = 2)
    expect_identical(wiz_samples()[[1]]$data, list(list(type = "quantdir", name = "Quant", path = "quant/S2")))
    session$setInputs(wiz_add_sample = 2)
    expect_identical(modal_field_errors(), list())
    expect_no_match(ch_last_modal(rec), "wiz-row-error", fixed = TRUE)
  })
})

test_that("a valid sample is appended as typed and trimmed, files in row order then salmon, the dialog closed", {
  proj <- ch_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  rec <- ch_mocks()
  shiny::testServer(ch_server(proj), {
    session$setInputs(wiz_data_root = file.path(proj, "data"))
    wiz_samples(list(list(condition = "First", colour = "#000000", data = list())))
    session$setInputs(wiz_add_sample = 1)
    session$setInputs(wiz_add_hic = 1)
    session$setInputs(wiz_add_bw = 1)
    session$setInputs(wiz_inst_name_1 = " HiC ", wiz_inst_path_1 = " hic/m.hic ",
                      wiz_inst_name_2 = "ATAC", wiz_inst_path_2 = "atac/a.bw",
                      wiz_quantdir_name = "Salmon", wiz_quantdir_path = "quant/S1",
                      wiz_sample_name = "  KO ", wiz_sample_colour = "#123456", wiz_confirm_sample = 1)
    expect_identical(wiz_samples()[[2]], list(condition = "KO", colour = "#123456", data = list(
      list(type = "hic", name = "HiC", path = "hic/m.hic"),
      list(type = "bw", name = "ATAC", path = "atac/a.bw"),
      list(type = "quantdir", name = "Salmon", path = "quant/S1"))))
    expect_identical(rec$removed, 1L)
    expect_null(modal_edit_name())
    expect_length(rec$modals, 1)  # no dialog re-shown
  })
})

test_that("a second Add carries the earlier samples' file rows (names, no paths) and the last salmon name", {
  proj <- ch_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  rec <- ch_mocks()
  shiny::testServer(ch_server(proj), {
    session$setInputs(wiz_data_root = file.path(proj, "data"))
    wiz_samples(list(
      list(condition = "A", colour = "#000000", data = list(list(type = "bw", name = "ATAC", path = "atac/a.bw"),
                                                         list(type = "quantdir", name = "Sal1", path = "quant/S1"))),
      list(condition = "B", colour = "#000000", data = list(list(type = "bw", name = "ATAC", path = "atac/b.bw"),
                                                         list(type = "hic", name = "HiC", path = "hic/m.hic"))),
      list(condition = "C", colour = "#000000", data = list())))
    session$setInputs(wiz_add_sample = 1)
    expect_identical(modal_instances(), list(list(id = 1, type = "bw", name = "ATAC", path = ""),
                                             list(id = 2, type = "hic", name = "HiC", path = "")))
    m <- ch_last_modal(rec)
    expect_match(m, "Add a sample", fixed = TRUE)
    # The salmon track name from the latest sample that has one; its path is not carried.
    expect_match(m, 'id="wiz_quantdir_name" type="text" class="shiny-input-text form-control" value="Sal1"', fixed = TRUE)
    expect_match(m, '<option value="" selected></option>', fixed = TRUE)
    rows <- ch_html(output$wiz_instance_rows)
    expect_match(rows, 'id="wiz_inst_name_1" type="text" class="shiny-input-text form-control" value="ATAC"', fixed = TRUE)
  })
})

test_that("Edit opens the sample's own values; saving keeps its place in the list", {
  proj <- ch_edit_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  rec <- ch_mocks()
  shiny::testServer(ch_server(proj), {
    session$setInputs(wiz_data_root = file.path(proj, "data"))
    wiz_samples(c(wiz_samples(), list(list(condition = "Zz", colour = "#999999", data = list()))))
    session$setInputs(wiz_edit_sample = "Ba")
    expect_identical(modal_edit_name(), "Ba")
    expect_identical(modal_instances(), list(list(id = 1, type = "bw", name = "ATAC", path = "atac/a.bw")))
    m <- ch_last_modal(rec)
    expect_match(m, "Edit sample", fixed = TRUE)
    expect_match(m, '<span class="action-label">Save changes</span>', fixed = TRUE)
    expect_match(m, 'id="wiz_sample_name" type="text" class="shiny-input-text form-control" value="Ba"', fixed = TRUE)
    expect_match(m, 'id="wiz_sample_colour" value="#1B9E77"', fixed = TRUE)
    expect_match(m, 'id="wiz_quantdir_name" type="text" class="shiny-input-text form-control" value="Salmon"', fixed = TRUE)
    expect_match(m, '<option value="quant/S1" selected>quant/S1</option>', fixed = TRUE)
    expect_match(ch_html(output$wiz_instance_rows), '<option value="atac/a.bw" selected>atac/a.bw</option>', fixed = TRUE)
    # Renamed and recoloured, the file kept (the browser sends the prefilled row's values).
    session$setInputs(wiz_inst_name_1 = "ATAC", wiz_inst_path_1 = "atac/a.bw",
                      wiz_quantdir_name = "Salmon", wiz_quantdir_path = "quant/S1",
                      wiz_sample_name = "Ba2", wiz_sample_colour = "#FF0000", wiz_confirm_sample = 1)
    expect_identical(ch_conds(wiz_samples()), c("Ba2", "Lp", "Zz"))
    expect_identical(wiz_samples()[[1]], list(condition = "Ba2", colour = "#FF0000", data = list(
      list(type = "bw", name = "ATAC", path = "atac/a.bw"), list(type = "quantdir", name = "Salmon", path = "quant/S1"))))
    expect_null(modal_edit_name())
  })
})

test_that("Edit: another sample's name is refused, its own is not; the rows' unsent values drop its files", {
  proj <- ch_edit_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  rec <- ch_mocks()
  shiny::testServer(ch_server(proj), {
    session$setInputs(wiz_data_root = file.path(proj, "data"))
    session$setInputs(wiz_edit_sample = "Lp")
    session$setInputs(wiz_sample_name = "Ba", wiz_confirm_sample = 1)
    expect_match(ch_last_modal(rec), "A condition named 'Ba' already exists.", fixed = TRUE)
    expect_match(ch_last_modal(rec), "Edit sample", fixed = TRUE)
    session$setInputs(wiz_sample_name = "Lp", wiz_confirm_sample = 2)
    expect_identical(ch_conds(wiz_samples()), c("Ba", "Lp"))
    # Here no browser sent the prefilled row's inputs, so the edited sample lost its file and
    # its colour came from the last typed colour input (none here: the default).
    expect_identical(wiz_samples()[[2]], list(condition = "Lp", colour = "#1f78b4", data = list()))
  })
})

test_that("Edit of a name not in the list does nothing; a confirm after its sample was removed warns", {
  proj <- ch_edit_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  rec <- ch_mocks()
  shiny::testServer(ch_server(proj), {
    session$flushReact()
    session$setInputs(wiz_edit_sample = "Nobody")
    expect_length(rec$modals, 0)
    expect_null(modal_edit_name())
    session$setInputs(wiz_edit_sample = "Lp")
    session$setInputs(wiz_remove_sample = "Lp")
    session$setInputs(wiz_sample_name = "Lp", wiz_confirm_sample = 1)
    expect_identical(rec$notes, "warning: Sample 'Lp' was removed, so these changes weren't saved.")
    expect_identical(rec$removed, 1L)
    expect_null(modal_edit_name())
    expect_identical(ch_conds(wiz_samples()), "Ba")
  })
})

test_that("an Edit opened before wiz_edit_sample's first flush is skipped as the initial value", {
  proj <- ch_edit_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  rec <- ch_mocks()
  shiny::testServer(ch_server(proj), {
    session$setInputs(wiz_edit_sample = "Ba")
    expect_length(rec$modals, 0)
    expect_null(modal_edit_name())
  })
})

# ---- Committed sample list ---------------------------------------------------------------

test_that("the sample list: empty message, then one row per sample with its swatch, files and buttons", {
  proj <- ch_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  shiny::testServer(ch_server(proj), {
    expect_identical(ch_html(output$wiz_samples_ui),
                     '<p class="wiz-empty">No samples yet. Add one per condition; you need at least one to save.</p>')
    wiz_samples(list(
      list(condition = "Ba", colour = "#1B9E77", data = list(list(type = "bw", name = "ATAC", path = "a.bw"),
                                                          list(type = "quantdir", name = "Quant", path = "q"))),
      list(condition = "No colour", colour = NULL, data = list())))
    session$flushReact()
    html <- ch_html(output$wiz_samples_ui)
    expect_match(html, '<div class="wiz-samples">', fixed = TRUE)
    expect_match(html, "background:#1B9E77;", fixed = TRUE)
    expect_match(html, "background:#ccc;", fixed = TRUE)
    expect_match(html, "<b>Ba</b>", fixed = TRUE)
    expect_match(html, '<div class="col-sm-5">bw:ATAC, quantdir:Quant</div>', fixed = TRUE)
    expect_match(html, '<div class="col-sm-5"></div>', fixed = TRUE)
    expect_match(html, 'id="wiz_edit_sample_btn_2"', fixed = TRUE)
    expect_match(html, 'id="wiz_rm_sample_btn_1"', fixed = TRUE)
    expect_match(html, "Shiny.setInputValue(&#39;wiz_edit_sample&#39;, &quot;No colour&quot;, {priority:&#39;event&#39;})", fixed = TRUE)
  })
})

test_that("Remove takes out the named sample only; an unknown name changes nothing", {
  proj <- ch_edit_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  shiny::testServer(ch_server(proj), {
    session$flushReact()
    session$setInputs(wiz_remove_sample = "Nobody")
    expect_identical(ch_conds(wiz_samples()), c("Ba", "Lp"))
    session$setInputs(wiz_remove_sample = "Ba")
    expect_identical(ch_conds(wiz_samples()), "Lp")
    expect_no_match(ch_html(output$wiz_samples_ui), "<b>Ba</b>", fixed = TRUE)
  })
})

# ---- Import ------------------------------------------------------------------------------

test_that("Import with no path, a missing file or a folder says so and imports nothing", {
  proj <- ch_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  shiny::testServer(ch_server(proj), {
    expect_identical(ch_html(output$wiz_import_status), "")
    session$setInputs(wiz_import_path = "  ", wiz_import_load = 1)
    expect_identical(ch_html(output$wiz_import_status), '<div class="wiz-error">Enter a path to samples.txt first.</div>')
    session$setInputs(wiz_import_path = file.path(proj, "none.txt"), wiz_import_load = 2)
    expect_identical(ch_html(output$wiz_import_status), sprintf('<div class="wiz-error">No such file: %s</div>', file.path(proj, "none.txt")))
    session$setInputs(wiz_import_path = proj, wiz_import_load = 3)
    expect_match(ch_html(output$wiz_import_status), "No such file:", fixed = TRUE)
    expect_identical(wiz_samples(), list())
    expect_identical(import_visual(), list())
  })
})

test_that("Import of a file with a bad row lists its problems and keeps the current samples", {
  proj <- ch_edit_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  rec <- ch_mocks()
  txt <- ch_txt(c("# fontsize: 12", "condition\tcolour\tfoo:ATAC", "Ba\t#1B9E77\tatac/a.bw", "Ba\t#1B9E77\tatac/b.bw"))
  shiny::testServer(ch_server(proj), {
    session$setInputs(wiz_import_path = txt, wiz_import_load = 1)
    html <- ch_html(output$wiz_import_status)
    expect_match(html, "<b>Nothing imported. The file has these problems:</b>", fixed = TRUE)
    expect_match(html, "<li>Unrecognized type prefix 'foo' in column 'foo:ATAC' (known types: bw, bed, hic, quantdir, paired-end)</li>", fixed = TRUE)
    expect_match(html, "<li>Duplicate condition 'Ba' (row 3)</li>", fixed = TRUE)
    expect_match(html, 'role="alert"', fixed = TRUE)
    expect_identical(ch_conds(wiz_samples()), c("Ba", "Lp"))
    expect_identical(import_visual(), list())
    expect_length(rec$updates, 0)
  })
})

test_that("Import of a good file: samples replaced, genome and folder filled, display settings kept for Save", {
  proj <- ch_edit_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  rec <- ch_mocks()
  txt <- ch_txt(c(paste0("# data_root: ", file.path(proj, "data")), "# gtf_path: /y/other.gtf.gz", "# orgdb: org.Hs.eg.db",
                  "# genome_label: hg38", "# fontsize: 12", "# zrange_cap: 50", "# show_title: FALSE", "# not_a_setting: 1",
                  "condition\tcolour\tbw:ATAC\tquantdir:Sal", "Ko\t#111111\tatac/a.bw\tquant/S1", "Wt\t\tatac/b.bw\t"))
  shiny::testServer(ch_server(proj), {
    session$setInputs(wiz_import_path = txt, wiz_import_load = 1)
    expect_identical(wiz_samples(), list(
      list(condition = "Ko", colour = "#111111", data = list(list(type = "bw", name = "ATAC", path = "atac/a.bw"),
                                                          list(type = "quantdir", name = "Sal", path = "quant/S1"))),
      list(condition = "Wt", colour = NA_character_, data = list(list(type = "bw", name = "ATAC", path = "atac/b.bw")))))
    expect_identical(ch_updates(rec), list(wiz_data_root = file.path(proj, "data"), wiz_gtf_path = "/y/other.gtf.gz",
                                           wiz_orgdb = "org.Hs.eg.db", wiz_genome_label = "hg38"))
    expect_identical(import_visual(), list(fontsize = 12, zrange_cap = 50, show_title = FALSE))
    expect_identical(ch_html(output$wiz_import_status), '<div class="wiz-ok">Imported 2 sample(s).</div>')
    expect_identical(wiz_mode(), "edit")  # the mode doesn't change
  })
})

test_that("Import with unresolved paths or an unusable auto number imports, with warnings", {
  proj <- ch_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  rec <- ch_mocks()
  txt <- ch_txt(c(paste0("# data_root: ", file.path(proj, "data")), "# hic_height: big",
                  "condition\tcolour\tbw:ATAC", "Ko\t#111111\tatac/none.bw"))
  shiny::testServer(ch_server(proj), {
    session$setInputs(wiz_import_path = txt, wiz_import_load = 1)
    expect_identical(ch_conds(wiz_samples()), "Ko")
    html <- ch_html(output$wiz_import_status)
    expect_match(html, '<div class="wiz-ok">Imported 1 sample(s).</div>', fixed = TRUE)
    expect_match(html, '<div class="wiz-warn wiz-warn-box" role="status">', fixed = TRUE)
    expect_match(html, "<b>Imported, with warnings:</b>", fixed = TRUE)
    expect_match(html, "<li>Row 2 (Ko): path for 'bw:ATAC' does not resolve to a file: atac/none.bw</li>", fixed = TRUE)
    expect_match(html, "<li>Ignoring '# hic_height: big' (not a number); using auto.</li>", fixed = TRUE)
    expect_identical(import_visual(), structure(list(), names = character(0)))
    # Only the folder is filled: the file names no genome.
    expect_identical(names(ch_updates(rec)), "wiz_data_root")
  })
})

test_that("Import resolves a relative data_root against the file's folder; no data_root fills nothing", {
  proj <- ch_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  rec <- ch_mocks()
  rel <- file.path(proj, "in.txt")
  writeLines(c("# data_root: data", "condition\tcolour\tbw:ATAC", "Ko\t#111111\tatac/a.bw"), rel)
  txt <- ch_txt(c("condition\tcolour\tbw:ATAC", "Ko\t#111111\tatac/a.bw"))
  shiny::testServer(ch_server(proj), {
    session$setInputs(wiz_import_path = rel, wiz_import_load = 1)
    expect_identical(ch_updates(rec), list(wiz_data_root = normalizePath(file.path(proj, "data"), winslash = "/")))
    session$setInputs(wiz_import_path = txt, wiz_import_load = 2)
    expect_length(rec$updates, 1)
    expect_identical(ch_html(output$wiz_import_status), '<div class="wiz-ok">Imported 1 sample(s).</div>')
  })
})

test_that("Import of a file with a column per sample warns about its layout", {
  proj <- ch_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  rec <- ch_mocks()
  txt <- ch_txt(c("condition\tcolour\tbw:Grans ATAC\tbw:pDC ATAC", "Grans\t#111111\tatac/a.bw\t", "pDC\t#222222\t\tatac/b.bw"))
  shiny::testServer(ch_server(proj), {
    session$setInputs(wiz_import_path = txt, wiz_import_load = 1)
    expect_match(ch_html(output$wiz_import_status), "<li>Columns bw:Grans ATAC, bw:pDC ATAC are each used by one sample only", fixed = TRUE)
    expect_identical(ch_conds(wiz_samples()), c("Grans", "pDC"))
  })
})

test_that("Import of a file with start lines takes its samples and no start lines as settings", {
  proj <- ch_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  rec <- ch_mocks()
  txt <- ch_txt(c("# start_gene: Elf5", "# start_tracks: ATAC, genes", "condition\tcolour\tbw:ATAC",
                  "Ba\t#1B9E77\tatac/a.bw", "Lp\t#D95F02\tatac/b.bw"))
  shiny::testServer(ch_server(proj), {
    session$setInputs(wiz_import_path = txt, wiz_import_load = 1)
    expect_identical(ch_conds(wiz_samples()), c("Ba", "Lp"))
    expect_identical(import_visual(), list())
    expect_identical(ch_html(output$wiz_import_status), '<div class="wiz-ok">Imported 2 sample(s).</div>')
  })
})

# ---- Export ------------------------------------------------------------------------------

test_that("Export opens a dialog proposing samples.txt in the data folder, or a blank path", {
  proj <- ch_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  rec <- ch_mocks()
  shiny::testServer(ch_server(proj), {
    session$setInputs(wiz_export_samples = 1)
    m <- ch_last_modal(rec)
    expect_match(m, "Export settings to txt", fixed = TRUE)
    expect_match(m, 'id="wiz_export_path" type="text" class="shiny-input-text form-control" value=""', fixed = TRUE)
    expect_match(m, 'id="wiz_confirm_export"', fixed = TRUE)
    expect_no_match(m, 'data-backdrop="static"', fixed = TRUE)  # easyClose
    session$setInputs(wiz_data_root = file.path(proj, "data"), wiz_export_samples = 2)
    expect_match(ch_last_modal(rec), sprintf('value="%s"', file.path(proj, "data", "samples.txt")), fixed = TRUE)
    expect_false(file.exists(file.path(proj, "data", "samples.txt")))  # nothing written by opening it
  })
})

test_that("Export with a blank path asks for one; a good path writes the form's samples.txt", {
  proj <- ch_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  rec <- ch_mocks()
  out <- file.path(proj, "out.txt")
  txt <- ch_txt(c("# fontsize: 12", "condition\tcolour\tbw:ATAC\tquantdir:Sal", "Ko\t#111111\tatac/a.bw\tquant/S1", "Wt\t#222222\tatac/b.bw\t"))
  shiny::testServer(ch_server(proj), {
    expect_identical(ch_html(output$wiz_export_status), "")
    session$setInputs(wiz_export_samples = 1)
    session$setInputs(wiz_export_path = " ", wiz_confirm_export = 1)
    expect_identical(ch_html(output$wiz_export_status), '<div class="wiz-error">Enter a path first.</div>')
    session$setInputs(wiz_import_path = txt, wiz_import_load = 1)
    session$setInputs(wiz_data_root = "data", wiz_gtf_path = " /x/genes.gtf.gz ", wiz_orgdb = "org.Mm.eg.db",
                      wiz_genome_label = "mm39")
    session$setInputs(wiz_export_path = out, wiz_confirm_export = 2)
    expect_identical(ch_html(output$wiz_export_status), sprintf('<div class="wiz-ok">Wrote %s.</div>', out))
    # The data folder as typed (not made absolute), the genome always, display settings that
    # differ from the defaults, then the table.
    expect_identical(readLines(out), c(
      "# data_root: data", "# gtf_path: /x/genes.gtf.gz", "# orgdb: org.Mm.eg.db", "# genome_label: mm39", "# fontsize: 12",
      "condition\tcolour\tbw:ATAC\tquantdir:Sal", "Ko\t#111111\tatac/a.bw\tquant/S1", "Wt\t#222222\tatac/b.bw\t"))
    expect_false(file.exists(ch_json(proj)))
  })
})

test_that("Export with blank genome fields leaves their lines out, and an unset folder is blank", {
  proj <- ch_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  rec <- ch_mocks()
  out <- file.path(proj, "out.txt")
  shiny::testServer(ch_server(proj), {
    session$setInputs(wiz_gtf_path = "", wiz_export_path = out, wiz_confirm_export = 1)
    expect_identical(readLines(out), c("# data_root: ", "condition\tcolour", character(0)))
  })
})

test_that("Export into a folder that doesn't exist says it failed", {
  proj <- ch_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  rec <- ch_mocks()
  out <- file.path(proj, "no", "such", "out.txt")
  shiny::testServer(ch_server(proj), {
    # current behaviour, possibly a bug: writeLines()'s "cannot open file" warning escapes the
    # observer alongside the shown error.
    expect_warning(session$setInputs(wiz_export_path = out, wiz_confirm_export = 1), "cannot open file")
    expect_identical(ch_html(output$wiz_export_status), '<div class="wiz-error">Failed to write: cannot open the connection</div>')
    expect_false(file.exists(out))
  })
})

test_that("Export in edit mode carries the project's non-default display settings", {
  proj <- ch_edit_project(list(page_w = 12, hic_norm = "VC")); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  rec <- ch_mocks()
  out <- file.path(proj, "out.txt")
  shiny::testServer(ch_server(proj), {
    ch_good_genome(session, proj)
    session$setInputs(wiz_export_path = out, wiz_confirm_export = 1)
    expect_identical(readLines(out)[1:7], c(
      paste0("# data_root: ", file.path(proj, "data")), "# gtf_path: /x/genes.gtf.gz", "# orgdb: org.Mm.eg.db",
      "# genome_label: mm10", "# page_w: 12", "# hic_norm: VC", "condition\tcolour\tbw:ATAC\tquantdir:Salmon"))
  })
})

# ---- Unsaved-changes guard and Close -----------------------------------------------------

test_that("a first run is never dirty", {
  proj <- ch_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  shiny::testServer(ch_server(proj), {
    expect_false(is_dirty())
    session$setInputs(wiz_data_root = proj, wiz_gtf_path = "/x/g.gtf")
    wiz_samples(list(list(condition = "A", colour = "#000000", data = list())))
    expect_false(is_dirty())
  })
})

test_that("edit mode: clean on the saved values, dirty after any change to a field or the samples", {
  proj <- ch_edit_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  shiny::testServer(ch_server(proj), {
    # current behaviour: before the browser sends the fields, the empty data folder differs
    # from the saved one.
    expect_true(is_dirty())
    ch_good_genome(session, proj)
    expect_false(is_dirty())
    session$setInputs(wiz_orgdb = "  org.Mm.eg.db  ")  # trimmed: the same
    expect_false(is_dirty())
    session$setInputs(wiz_genome_label = "mm39")
    expect_true(is_dirty())
    session$setInputs(wiz_genome_label = "mm10", wiz_data_root = paste0(file.path(proj, "data"), "/"))
    expect_true(is_dirty())  # a trailing slash is a different folder here
    session$setInputs(wiz_data_root = file.path(proj, "data"))
    expect_false(is_dirty())
    keep <- wiz_samples()
    wiz_samples(rev(keep))
    expect_true(is_dirty())
    wiz_samples(keep)
    expect_false(is_dirty())
  })
})

test_that("edit mode: an imported display setting makes the form dirty", {
  proj <- ch_edit_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  rec <- ch_mocks()
  shiny::testServer(ch_server(proj), {
    ch_good_genome(session, proj)
    expect_false(is_dirty())
    import_visual(list(fontsize = 8))
    expect_false(is_dirty())  # the same value as saved
    import_visual(list(fontsize = 9))
    expect_true(is_dirty())
  })
})

test_that("Close with nothing changed shows the closing note and stops with restart", {
  proj <- ch_edit_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  rec <- ch_mocks()
  shiny::testServer(ch_server(proj), {
    ch_good_genome(session, proj)
    expect_message(session$setInputs(wiz_close = 1), "PlotTrellis: settings unchanged, returning to the app.", fixed = TRUE)
    m <- ch_last_modal(rec)
    expect_match(m, "Closing", fixed = TRUE)
    expect_match(m, "Returning to the app without saving.", fixed = TRUE)
    expect_no_match(m, "modal-footer", fixed = TRUE)
    expect_length(rec$later, 1)
    expect_identical(rec$later[[1]]$delay, 0.2)
    expect_length(rec$stopped, 0)
    rec$later[[1]]$func()
    expect_identical(rec$stopped, list("restart"))
  })
  expect_false(file.exists(file.path(proj, "project_samples.txt")))  # nothing written
})

test_that("Close with changes asks first; Keep editing stays, Discard closes", {
  proj <- ch_edit_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  rec <- ch_mocks()
  before <- readLines(ch_json(proj))
  shiny::testServer(ch_server(proj), {
    ch_good_genome(session, proj)
    session$setInputs(wiz_genome_label = "hg38")
    session$setInputs(wiz_close = 1)
    m <- ch_last_modal(rec)
    expect_match(m, "Discard changes?", fixed = TRUE)
    expect_match(m, 'id="wiz_confirm_close"', fixed = TRUE)
    expect_match(m, "Keep editing", fixed = TRUE)
    expect_length(rec$later, 0)
    expect_message(session$setInputs(wiz_confirm_close = 1), "settings unchanged")
    expect_identical(rec$removed, 1L)
    expect_match(ch_last_modal(rec), "Returning to the app without saving.", fixed = TRUE)
    rec$later[[1]]$func()
    expect_identical(rec$stopped, list("restart"))
  })
  expect_identical(readLines(ch_json(proj)), before)
})

# ---- Save --------------------------------------------------------------------------------

test_that("Save with nothing filled in lists every problem, in order, and writes nothing", {
  proj <- ch_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  rec <- ch_mocks()
  shiny::testServer(ch_server(proj), {
    session$setInputs(wiz_mode_manual = 1)
    expect_false(save_failed())
    expect_identical(ch_html(output$wiz_completion_error), "")
    session$setInputs(wiz_save = 1)
    expect_true(save_failed())
    expect_identical(save_problems(), c("Data root is required.", "Add at least one sample before saving.",
                                        "GTF path is required.", "OrgDb package is required.", PT_GENOME_LABEL_REQUIRED))
    html <- ch_html(output$wiz_completion_error)
    expect_match(html, '<div class="wiz-error wiz-error-box" role="alert">', fixed = TRUE)
    expect_match(html, "<b>Not saved yet. Please fix:</b>", fixed = TRUE)
    expect_identical(lengths(regmatches(html, gregexpr("<li>", html, fixed = TRUE))), 5L)
    expect_length(rec$modals, 0)
    expect_length(rec$later, 0)
  })
  expect_false(file.exists(ch_json(proj)))
  expect_false(file.exists(file.path(proj, "project_samples.txt")))
})

test_that("Save's problems: a missing folder, an OrgDb not installed, an unknown genome label", {
  proj <- ch_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  rec <- ch_mocks()
  shiny::testServer(ch_server(proj), {
    wiz_samples(list(list(condition = "A", colour = "#000000", data = list())))
    session$setInputs(wiz_data_root = file.path(proj, "nope"), wiz_gtf_path = "/x/g.gtf",
                      wiz_orgdb = "org.Xx.eg.db", wiz_genome_label = "m10", wiz_save = 1)
    expect_identical(save_problems(), c(
      "Data root directory does not exist.",
      "OrgDb package 'org.Xx.eg.db' is not installed. Install it first, e.g. BiocManager::install('org.Xx.eg.db').",
      "Genome label 'm10' isn't a known genome build id (e.g. mm10, mm39, hg38)."))
    # A blank genome label has its own message.
    session$setInputs(wiz_genome_label = " ")
    expect_identical(save_problems()[3], PT_GENOME_LABEL_REQUIRED)
  })
  expect_false(file.exists(ch_json(proj)))
})

test_that("before any Save the problems are computed but not shown", {
  proj <- ch_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  shiny::testServer(ch_server(proj), {
    expect_length(save_problems(), 5)
    expect_identical(ch_html(output$wiz_completion_error), "")
  })
})

test_that("a good first Save writes the settings and the samples copy, says so, and stops with restart", {
  skip_if_not_installed("org.Mm.eg.db")
  proj <- ch_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  rec <- ch_mocks()
  shiny::testServer(ch_server(proj), {
    session$setInputs(wiz_mode_manual = 1)
    wiz_samples(list(list(condition = "Ba", colour = "#1B9E77", data = list(list(type = "bw", name = "ATAC", path = "atac/a.bw")))))
    session$setInputs(wiz_data_root = paste0(" ", file.path(proj, "data")), wiz_gtf_path = " /x/genes.gtf.gz",
                      wiz_orgdb = "org.Mm.eg.db ", wiz_genome_label = "mm10")
    expect_message(session$setInputs(wiz_save = 1), "PlotTrellis: settings saved, launching the app.", fixed = TRUE)
    expect_false(save_failed())
    expect_null(write_error())
    m <- ch_last_modal(rec)
    expect_match(m, "Settings saved", fixed = TRUE)
    expect_match(m, "project_settings.json has been written.", fixed = TRUE)
    expect_match(m, sprintf("A copy of your samples is in %s; it can rebuild the project", file.path(proj, "project_samples.txt")), fixed = TRUE)
    expect_no_match(m, "wiz-warn", fixed = TRUE)
    expect_length(rec$later, 1)
    expect_identical(rec$later[[1]]$delay, 0.2)
    rec$later[[1]]$func()
    expect_identical(rec$stopped, list("restart"))
  })
  saved <- jsonlite::fromJSON(ch_json(proj), simplifyDataFrame = FALSE)
  expect_identical(names(saved), c(
    "data_dir", "samples", "gtf_path", "orgdb", "genome_label", "suggested_genes", "font", "fontsize", "page_w", "chr_fmt",
    "show_track_axes", "show_title", "show_scale_bar", "hic_norm", "hic_style", "hic_colour_scale", "zrange_cap", "atac_height",
    "hic_height", "hic_res_mode", "hic_max_bins_per_view", "hic_resolution", "diff_colour_lo", "diff_colour_hi",
    "junction_opacity_cap", "min_junction_support", "reference_colour", "comparison_colour", "gene_flank_pct",
    "min_display_tpm", "uniprot_reviewed_only", "arc_clip_partial", "junction_height", "diff_shared_alpha"))
  expect_identical(saved$data_dir, file.path(proj, "data"))  # trimmed
  expect_identical(saved$gtf_path, "/x/genes.gtf.gz")
  expect_identical(saved$orgdb, "org.Mm.eg.db")
  expect_identical(saved$samples, list(list(condition = "Ba", colour = "#1B9E77",
                                            data = list(list(type = "bw", name = "ATAC", path = "atac/a.bw")))))
  d <- .pt_display_defaults()
  for (k in setdiff(names(d), c("zrange_cap", "hic_height", "hic_resolution"))) expect_equal(saved[[k]], d[[k]], info = k)
  raw <- readLines(ch_json(proj))
  for (k in c("zrange_cap", "hic_height", "hic_resolution")) expect_true(any(grepl(sprintf('"%s": null', k), raw, fixed = TRUE)), info = k)
  expect_identical(readLines(file.path(proj, "project_samples.txt")), c(
    paste0("# data_root: ", normalizePath(file.path(proj, "data"))), "# gtf_path: /x/genes.gtf.gz", "# orgdb: org.Mm.eg.db",
    "# genome_label: mm10", "condition\tcolour\tbw:ATAC", "Ba\t#1B9E77\tatac/a.bw"))
  expect_identical(sort(list.files(proj)), c("data", "outside.bw", "project_samples.txt", "project_settings.json"))
})

test_that("a Save after an Import writes the imported display settings over the defaults", {
  skip_if_not_installed("org.Mm.eg.db")
  proj <- ch_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  rec <- ch_mocks()
  txt <- ch_txt(c("# page_w: 11", "# chr_fmt: plain", "# hic_resolution: 5000", "condition\tcolour\tbw:ATAC", "Ba\t#1B9E77\tatac/a.bw"))
  shiny::testServer(ch_server(proj), {
    session$setInputs(wiz_mode_import = 1)
    session$setInputs(wiz_import_path = txt, wiz_import_load = 1)
    ch_good_genome(session, proj)
    expect_message(session$setInputs(wiz_save = 1), "settings saved")
  })
  saved <- jsonlite::fromJSON(ch_json(proj), simplifyDataFrame = FALSE)
  expect_equal(saved$page_w, 11)
  expect_identical(saved$chr_fmt, "plain")
  expect_equal(saved$hic_resolution, 5000)
  expect_false("start_view" %in% names(saved))
  expect_true(any(grepl("^# page_w: 11$", readLines(file.path(proj, "project_samples.txt")))))
})

test_that("an edit-mode Save keeps the samples.txt source recorded in the old settings", {
  skip_if_not_installed("org.Mm.eg.db")
  proj <- ch_edit_project(list(samples_txt_source = "/somewhere/samples.txt", page_w = 12))
  on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  rec <- ch_mocks()
  shiny::testServer(ch_server(proj), {
    ch_good_genome(session, proj)
    expect_message(session$setInputs(wiz_save = 1), "settings saved")
  })
  saved <- jsonlite::fromJSON(ch_json(proj), simplifyDataFrame = FALSE)
  expect_identical(saved$samples_txt_source, "/somewhere/samples.txt")
  expect_identical(utils::tail(names(saved), 1), "samples_txt_source")
  expect_equal(saved$page_w, 12)
  expect_identical(ch_conds(saved$samples), c("Ba", "Lp"))
})

test_that("a Save whose samples copy can't be written still saves, and says so", {
  skip_if_not_installed("org.Mm.eg.db")
  proj <- ch_edit_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  dir.create(file.path(proj, "project_samples.txt"))  # a folder in the copy's place
  rec <- ch_mocks()
  shiny::testServer(ch_server(proj), {
    ch_good_genome(session, proj)
    session$setInputs(wiz_genome_label = "mm39")
    expect_message(session$setInputs(wiz_save = 1), "settings saved")
    expect_match(ch_last_modal(rec), "Couldn't write project_samples.txt beside it (use Export samples.txt to save a copy yourself).", fixed = TRUE)
    expect_length(rec$later, 1)
  })
  expect_identical(jsonlite::fromJSON(ch_json(proj))$genome_label, "mm39")
})

test_that("a Save into a folder that can't be written keeps the session, and Back clears the message", {
  proj <- ch_project()
  on.exit({ Sys.chmod(proj, "0755"); unlink(proj, recursive = TRUE) }, add = TRUE)
  rec <- ch_mocks()
  Sys.chmod(proj, "0555")
  skip_if(file.access(proj, 2) == 0, "can't make a read-only folder here (running as root?)")
  shiny::testServer(ch_server(proj), {
    session$setInputs(wiz_mode_manual = 1)
    wiz_samples(list(list(condition = "Ba", colour = "#1B9E77", data = list())))
    ch_good_genome(session, proj)
    session$setInputs(wiz_save = 1)
    expect_match(write_error(), paste0("^Couldn't write ", ch_json(proj), " \\(.*Permission denied.*\\)\\. ",
                                       "Check that you can write to that folder, then Save again\\.$"))
    expect_false(save_failed())
    expect_length(rec$modals, 0)
    expect_length(rec$later, 0)
    expect_match(ch_html(output$wiz_completion_error), "<b>Not saved:</b>", fixed = TRUE)
    session$setInputs(wiz_mode_back = 1)
    expect_null(write_error())
    expect_identical(ch_html(output$wiz_completion_error), "")
  })
})

test_that("a failed Save then Back clears the list; the samples and fields stay", {
  proj <- ch_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  shiny::testServer(ch_server(proj), {
    session$setInputs(wiz_mode_manual = 1)
    wiz_samples(list(list(condition = "Ba", colour = "#1B9E77", data = list())))
    session$setInputs(wiz_data_root = file.path(proj, "data"), wiz_save = 1)
    expect_true(save_failed())
    expect_match(ch_html(output$wiz_completion_error), "GTF path is required.", fixed = TRUE)
    session$setInputs(wiz_mode_back = 1)
    expect_false(save_failed())
    expect_identical(ch_html(output$wiz_completion_error), "")
    expect_identical(data_root(), file.path(proj, "data"))
    expect_length(wiz_samples(), 1)
    # Manual again: still no list until the next Save.
    session$setInputs(wiz_mode_manual = 2)
    expect_identical(ch_html(output$wiz_completion_error), "")
    session$setInputs(wiz_save = 2)
    expect_match(ch_html(output$wiz_completion_error), "GTF path is required.", fixed = TRUE)
  })
})

test_that("a Save after a failed one clears the list once it succeeds", {
  skip_if_not_installed("org.Mm.eg.db")
  proj <- ch_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  rec <- ch_mocks()
  shiny::testServer(ch_server(proj), {
    wiz_samples(list(list(condition = "Ba", colour = "#1B9E77", data = list())))
    session$setInputs(wiz_save = 1)
    expect_true(save_failed())
    ch_good_genome(session, proj)
    expect_identical(ch_html(output$wiz_completion_error), "")  # the live re-check
    expect_true(save_failed())
    expect_message(session$setInputs(wiz_save = 2), "settings saved")
    expect_false(save_failed())
  })
  expect_true(file.exists(ch_json(proj)))
})
