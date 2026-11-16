# Bug 70 (part): the wizard's Save keeps the project's start view (pt_preload(), or a
# samples.txt's start_gene/start_tracks) while the genome, samples and tracks it refers to are
# unchanged, and drops it otherwise. Drives the real wizard_server() headlessly.

sv_project <- function() {
  proj <- tempfile("pt-sv-")
  dir.create(file.path(proj, "data", "atac"), recursive = TRUE)
  file.create(file.path(proj, "data", "atac", "a.bw"))
  proj
}
sv_wizard <- function(proj) {
  json <- file.path(proj, "project_settings.json")
  defaults <- .pt_default_settings()
  wizard_server(json, .pt_read_settings(json, defaults), defaults)
}
# A wizard Save. With `txt` (samples.txt lines) the samples are imported first, as on a first
# run; without, the wizard opens on the saved settings and only the genome/folder inputs are
# set (testServer() has no UI to prefill them), as an Edit project settings, then Save.
sv_save <- function(proj, gtf = "/x/genes.gtf.gz", txt = NULL) {
  f <- NULL
  if (!is.null(txt)) { f <- tempfile(fileext = ".txt"); writeLines(txt, f); on.exit(unlink(f), add = TRUE) }
  eval(bquote(shiny::testServer(sv_wizard(.(proj)), {
    if (!is.null(.(f))) session$setInputs(wiz_import_path = .(f), wiz_import_load = 1)
    session$setInputs(wiz_data_root = .(file.path(proj, "data")), wiz_gtf_path = .(gtf),
                      wiz_orgdb = "org.Mm.eg.db", wiz_genome_label = "mm10")
    session$setInputs(wiz_save = 1)
  })))
}
sv_settings <- function(proj) jsonlite::read_json(file.path(proj, "project_settings.json"))
BA <- c("condition\tcolour\tbw:ATAC", "Ba\t#1B9E77\tatac/a.bw")

# A project saved by the wizard with sample Ba, then given a start view on its ATAC track.
sv_preloaded <- function() {
  proj <- sv_project()
  sv_save(proj, txt = BA)
  keys <- names(build_track_items(load_samples_from_settings(sv_settings(proj)$samples, file.path(proj, "data"))))
  expect_true("bw_ATAC_Ba" %in% keys)
  .pt_write_start_view(proj, "single", "Ba", "Elf5", c("bw_ATAC_Ba", "genes"))
  proj
}

test_that("a wizard Save with nothing changed keeps the start view", {
  proj <- sv_preloaded()
  on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  before <- sv_settings(proj)$start_view
  expect_false(is.null(before))
  mirror <- file.path(proj, "project_samples.txt")
  unlink(mirror)
  sv_save(proj)
  expect_true(file.exists(mirror))  # the Save went through (it rewrites the mirror)
  expect_identical(sv_settings(proj)$start_view, before)
})

test_that("a wizard Save drops the start view when what it refers to changed", {
  # Another GTF: its gene may not be there.
  proj <- sv_preloaded()
  on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  sv_save(proj, gtf = "/x/other.gtf.gz")
  expect_null(sv_settings(proj)$start_view)

  # Any Import, even of the same samples: the file restates the project (this one has no
  # start lines), so the old view isn't brought back.
  proj3 <- sv_preloaded()
  on.exit(unlink(proj3, recursive = TRUE), add = TRUE)
  sv_save(proj3, txt = BA)
  expect_true(file.exists(file.path(proj3, "project_settings.json")))
  expect_null(sv_settings(proj3)$start_view)

  # Its sample renamed (re-imported as Bb): neither the sample nor its track is there.
  proj2 <- sv_preloaded()
  on.exit(unlink(proj2, recursive = TRUE), add = TRUE)
  sv_save(proj2, txt = c("condition\tcolour\tbw:ATAC", "Bb\t#1B9E77\tatac/a.bw"))
  expect_null(sv_settings(proj2)$start_view)
})

# Two samples, with the start lines pt_project() reads.
BA_LP_START <- c("# start_gene: Elf5", "# start_tracks: ATAC, genes",
                 "condition\tcolour\tbw:ATAC", "Ba\t#1B9E77\tatac/a.bw", "Lp\t#D95F02\tatac/a.bw")

test_that("the wizard's Import carries a samples.txt's start lines into the start view, as pt_project() does", {
  # First run: no project yet, the file imported in the wizard.
  proj <- sv_project()
  on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  sv_save(proj, txt = BA_LP_START)
  got <- sv_settings(proj)$start_view
  expect_identical(got, list(view_mode = "two", conditions = list("Ba", "Lp"), gene = "Elf5",
                             tracks = list("bw_ATAC_Ba", "bw_ATAC_Lp", "genes")))

  # The same file through pt_project()'s importer gives the same view.
  ref <- tempfile("pt-sv-ref-")
  on.exit(unlink(ref, recursive = TRUE), add = TRUE)
  txt <- tempfile(fileext = ".txt"); on.exit(unlink(txt), add = TRUE)
  writeLines(c(paste0("# data_root: ", file.path(proj, "data")), "# gtf_path: /x/genes.gtf.gz",
               "# orgdb: org.Mm.eg.db", "# genome_label: mm10", BA_LP_START), txt)
  .pt_import(ref, txt)
  expect_identical(jsonlite::read_json(file.path(ref, "project_settings.json"))$start_view, got)

  # A later Import replaces it: by its own start lines, then by none.
  sv_save(proj, txt = c("# start_gene: Cat", "# start_tracks: ATAC", "# start_samples: Lp", BA_LP_START[-(1:2)]))
  expect_identical(sv_settings(proj)$start_view, list(view_mode = "single", conditions = list("Lp"),
                                                      gene = "Cat", tracks = list("bw_ATAC_Lp")))
  sv_save(proj, txt = BA_LP_START[-(1:2)])
  expect_null(sv_settings(proj)$start_view)
})

test_that("start lines that don't fit the imported samples don't block the wizard's Save", {
  proj <- sv_project()
  on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  expect_message(sv_save(proj, txt = c("# start_gene: Elf5", "# start_tracks: DNase", BA_LP_START[-(1:2)])),
                 "start lines: .*DNase")
  expect_true(file.exists(file.path(proj, "project_settings.json")))
  expect_null(sv_settings(proj)$start_view)
})

test_that("writing a start view needs only the settings file writable, as before the copy carried it", {
  proj <- sv_project()
  on.exit({ Sys.chmod(proj, "0755"); unlink(proj, recursive = TRUE) }, add = TRUE)
  sv_save(proj, txt = BA)
  copy <- file.path(proj, "project_samples.txt")
  # A folder that can't be written (the settings file can): no temporary file to rename.
  Sys.chmod(proj, "0555")
  skip_if(file.access(proj, 2) == 0, "can't make a read-only folder here (running as root?)")
  expect_no_warning(.pt_write_start_view(proj, "single", "Ba", "Elf5", c("bw_ATAC_Ba", "genes")))
  expect_identical(sv_settings(proj)$start_view$gene, "Elf5")
  # A copy that can't be written: no warning, the view still saved.
  Sys.chmod(proj, "0755"); Sys.chmod(copy, "0444")
  expect_no_warning(.pt_write_start_view(proj, "single", "Ba", "Cat", c("bw_ATAC_Ba", "genes")))
  expect_identical(sv_settings(proj)$start_view$gene, "Cat")
})

test_that("project_samples.txt carries the start view, so rebuilding the project from it restores it", {
  # The bundled example (start lines in its samples.txt), then a Two view with a comparison
  # track (pt_preload()): each time the copy's start lines give the same view back.
  proj <- tempfile("pt-sv-ex-")
  on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  p <- suppressMessages(suppressWarnings(pt_project(proj, pt_example())))
  json <- file.path(proj, "project_settings.json")
  cmp <- pt_tracks(p)$track[pt_tracks(p)$kind == "comparison"][1]
  for (fig in list(NULL, pt_figure(p, "Elf5", c("genes", "atac", cmp), comparison = c("Basal", "LP")))) {
    if (!is.null(fig)) pt_preload(p, fig)
    s <- jsonlite::read_json(json)
    expect_false(is.null(s$start_view))
    copy <- readLines(file.path(proj, "project_samples.txt"))
    expect_identical(.pt_txt_start_view(copy, s), s$start_view)
  }
  # Reset, then Import the copy in the wizard (the genome inputs set by hand: testServer()
  # has no UI for the Import to fill).
  unlink(json)
  shiny::testServer(sv_wizard(proj), {
    session$setInputs(wiz_import_path = file.path(proj, "project_samples.txt"), wiz_import_load = 1)
    session$setInputs(wiz_data_root = s$data_dir, wiz_gtf_path = s$gtf_path,
                      wiz_orgdb = s$orgdb, wiz_genome_label = s$genome_label)
    session$setInputs(wiz_save = 1)
  })
  expect_identical(jsonlite::read_json(json)$start_view, s$start_view)
})

test_that("rebuilding a project from its own project_samples.txt keeps the view and the copy's start lines", {
  proj <- tempfile("pt-sv-rb-")
  on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  p <- suppressMessages(suppressWarnings(pt_project(proj, pt_example())))
  json <- file.path(proj, "project_settings.json")
  copy <- file.path(proj, "project_samples.txt")
  sv <- jsonlite::read_json(json)$start_view
  expect_false(is.null(sv))
  lines <- grep("^# start_", readLines(copy), value = TRUE)
  expect_true(length(lines) >= 2)
  # pt_project() from the copy (also run_plottrellis(samples_txt = <the copy>)).
  suppressMessages(suppressWarnings(pt_project(proj, copy, overwrite = TRUE)))
  expect_identical(jsonlite::read_json(json)$start_view, sv)
  expect_identical(grep("^# start_", readLines(copy), value = TRUE), lines)
  # parse_samples_txt() from the copy, the rebuild the unreadable-settings error suggests:
  # the copy keeps its start lines.
  suppressMessages(suppressWarnings(parse_samples_txt(copy, project_dir = proj, overwrite = TRUE)))
  expect_identical(grep("^# start_", readLines(copy), value = TRUE), lines)
})

test_that("a start view of fixed tracks only has no start_samples line, and an all-NA one reads without a warning", {
  cfg <- list(data_dir = file.path(sv_project(), "data"),
              samples = list(list(condition = "Lp", colour = "#1B9E77",
                                  data = list(list(type = "bw", name = "ATAC", path = "atac/a.bw")))))
  on.exit(unlink(dirname(cfg$data_dir), recursive = TRUE), add = TRUE)
  fixed <- .pt_start_view_entry("compose", character(0), "Elf5", c("genes", "transcripts_annotated"))
  lines <- .pt_start_view_lines(fixed, cfg)
  expect_identical(lines, c("# start_gene: Elf5", "# start_tracks: genes, transcripts_annotated"))
  expect_identical(.pt_txt_start_view(lines, cfg), fixed)
  # A copy written before this read "NA, NA" as samples and warned that they have no effect.
  expect_no_warning(got <- .pt_txt_start_view(c(lines, "# start_samples: NA, NA"), cfg))
  expect_identical(got, fixed)
  # Only for fixed tracks: a per-sample track with NA for its sample keeps its error.
  expect_error(.pt_txt_start_view(c("# start_gene: Elf5", "# start_tracks: genes, ATAC", "# start_samples: NA, NA"), cfg),
               "needs a sample")
})

test_that("a start view no set of start lines gives exactly is left out of the copy", {
  # Multi (compose) with 2 samples: read back, the lines would give Two.
  cfg <- list(data_dir = file.path(sv_project(), "data"),
              samples = lapply(c("Ba", "Lp"), function(n) list(condition = n, colour = "#1B9E77",
                data = list(list(type = "bw", name = "ATAC", path = "atac/a.bw")))))
  on.exit(unlink(dirname(cfg$data_dir), recursive = TRUE), add = TRUE)
  multi <- .pt_start_view_entry("compose", character(0), "Elf5", c("bw_ATAC_Ba", "bw_ATAC_Lp"))
  expect_identical(.pt_start_view_lines(multi, cfg), character(0))
  one <- .pt_start_view_entry("single", "Lp", "Elf5", c("genes", "bw_ATAC_Lp"))
  expect_identical(.pt_start_view_lines(one, cfg),
                   c("# start_gene: Elf5", "# start_tracks: genes, ATAC", "# start_samples: NA, Lp"))
  expect_identical(.pt_start_view_lines(.pt_start_view_entry("single", "Lp", "Elf5", "bw_gone_Lp"), cfg), character(0))
})

test_that(".pt_kept_start_view() needs every sample and track of the view", {
  old <- list(gtf_path = "/x/g.gtf.gz", orgdb = "org.Mm.eg.db",
              start_view = list(view_mode = "single", conditions = list("Ba"), gene = "Elf5",
                                tracks = list("bw_ATAC_Ba", "genes")))
  smp <- function(name) list(list(condition = "Ba", colour = "#1B9E77",
                                  data = list(list(type = "bw", name = name, path = "atac/a.bw"))))
  new <- list(gtf_path = "/x/g.gtf.gz", orgdb = "org.Mm.eg.db", data_dir = tempdir(), samples = smp("ATAC"))
  expect_identical(.pt_kept_start_view(old, new), old$start_view)
  with <- function(key, value) { x <- new; x[[key]] <- value; x }
  expect_null(.pt_kept_start_view(old, with("samples", smp("DNase"))))  # its track is gone
  expect_null(.pt_kept_start_view(old, with("orgdb", "org.Hs.eg.db")))
  expect_null(.pt_kept_start_view(old[c("gtf_path", "orgdb")], new))  # no view to keep
})

test_that("a wizard Save into a folder that can't be written shows a message, not a closed session", {
  proj <- sv_project()
  on.exit({ Sys.chmod(proj, "0755"); unlink(proj, recursive = TRUE) }, add = TRUE)
  f <- tempfile(fileext = ".txt"); writeLines(BA, f); on.exit(unlink(f), add = TRUE)
  Sys.chmod(proj, "0555")
  skip_if(file.access(proj, 2) == 0, "can't make a read-only folder here (running as root?)")
  shiny::testServer(sv_wizard(proj), {
    session$setInputs(wiz_import_path = f, wiz_import_load = 1)
    session$setInputs(wiz_data_root = file.path(proj, "data"), wiz_gtf_path = "/x/genes.gtf.gz",
                      wiz_orgdb = "org.Mm.eg.db", wiz_genome_label = "mm10")
    session$setInputs(wiz_save = 1)
    expect_match(paste(as.character(output$wiz_completion_error$html)), "Not saved:.*Couldn't write .*project_settings.json.*Permission denied")
  })
  expect_false(file.exists(file.path(proj, "project_settings.json")))
})

test_that("a write failure's message goes at the next Save, which shows its own result", {
  proj <- sv_project()
  on.exit({ Sys.chmod(proj, "0755"); unlink(proj, recursive = TRUE) }, add = TRUE)
  f <- tempfile(fileext = ".txt"); writeLines(BA, f); on.exit(unlink(f), add = TRUE)
  Sys.chmod(proj, "0555")
  skip_if(file.access(proj, 2) == 0, "can't make a read-only folder here (running as root?)")
  shiny::testServer(sv_wizard(proj), {
    shown <- function() paste(as.character(output$wiz_completion_error$html))
    session$setInputs(wiz_import_path = f, wiz_import_load = 1)
    session$setInputs(wiz_data_root = file.path(proj, "data"), wiz_gtf_path = "/x/genes.gtf.gz",
                      wiz_orgdb = "org.Mm.eg.db", wiz_genome_label = "mm10")
    session$setInputs(wiz_save = 1)
    expect_match(shown(), "Permission denied")
    # The folder fixed, a field emptied: the next Save lists that field, not the old write error.
    Sys.chmod(proj, "0755")
    session$setInputs(wiz_gtf_path = "")
    session$setInputs(wiz_save = 2)
    expect_match(shown(), "Not saved yet. Please fix")
    expect_no_match(shown(), "Permission denied")
    # Fixing the field clears the list as it is fixed (the live re-check).
    session$setInputs(wiz_gtf_path = "/x/genes.gtf.gz")
    expect_identical(shown(), character(0))
  })
})

test_that("a failed settings write leaves the saved project_settings.json as it was", {
  proj <- tempfile("pt-full-"); dir.create(proj)
  on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  json <- file.path(proj, "project_settings.json")
  writeLines('{"page_w": 7}', json)
  # A write that fails after the file is opened (as on a full disk; writeLines() opens, so
  # empties, its file before writing): the old file must survive and no temporary file stay.
  local_mocked_bindings(toJSON = function(...) 1L, .package = "jsonlite")
  expect_error(.pt_write_project_files(list(page_w = 9), json), "character")
  expect_identical(readLines(json), '{"page_w": 7}')
  expect_identical(list.files(proj), "project_settings.json")
})
