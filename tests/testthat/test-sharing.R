# Sharing between people: the sender gives the receiver a data folder and its samples.txt, a
# project folder, the project's project_samples.txt or a downloaded R script. Each test checks
# the receiver gets the sender's figure (the saved SVG, byte for byte) or an error naming what
# is missing. The sender's own folders are moved away first wherever the receiver should not
# need them, so a path still pointing at them fails the test. Helpers: helper-sharing.R.

test_that("a samples.txt with relative paths, copied with its data to another folder, draws the same figure there", {
  skip_if_not_installed("org.Hs.eg.db")
  local_example_chrom_sizes()
  sender <- share_sender()
  # The receiver's copy, at a different path (spaces again); the sender's data and project gone.
  home <- file.path(tempfile("pt-share-"), "receiver home")
  dir.create(file.path(home, "their data"), recursive = TRUE)
  stopifnot(all(file.copy(list.files(sender$data, full.names = TRUE), file.path(home, "their data"),
                          recursive = TRUE)))
  st <- file.path(home, "their data", "samples.txt")
  expect_true(any(readLines(st) == "# data_root: ."))
  # The check below can see the sender's paths: the sender's own figure holds them.
  expect_true(any(grepl(dirname(sender$home), share_strings(sender$fig), fixed = TRUE)))
  local_moved_away(sender$data)
  local_moved_away(sender$dir)
  # Opened from an unrelated working directory: it doesn't matter where R is.
  withr::local_dir(withr::local_tempdir())
  proj <- share_project(file.path(home, "PlotTrellis project"), st)
  fig <- share_figure(proj)
  expect_identical(share_svg(fig), sender$md5)
  # Nothing the receiver's project or figure holds (data, GTF, genome cache, cache folder, the
  # OrgDb) points back at the sender's folders.
  strs <- share_strings(proj, fig)
  expect_false(any(grepl(dirname(sender$home), strs, fixed = TRUE)))
  expect_true(any(startsWith(strs, normalizePath(home))))
  expect_identical(proj$cfg$orgdb, "org.Hs.eg.db")
})

test_that("a samples.txt one folder above its data (data_root: <folder>) travels too, opened by relative paths from any working directory", {
  skip_if_not_installed("org.Hs.eg.db")
  local_example_chrom_sizes()
  sender <- share_sender()
  top <- file.path(tempfile("pt-share-"), "shared folder")
  st <- share_copy_example(file.path(top, "my data"))
  txt <- readLines(st)
  txt <- sub("^# data_root: .*", "# data_root: my data", txt)
  txt <- sub("^# gtf_path: (.*)", "# gtf_path: my data/\\1", txt)
  writeLines(txt, file.path(top, "samples.txt"))
  unlink(st)
  # The receiver copies the whole folder, the original is moved away.
  recv <- file.path(tempfile("pt-share-"), "receiver home")
  dir.create(recv, recursive = TRUE)
  stopifnot(file.copy(top, recv, recursive = TRUE))
  local_moved_away(top)
  # samples_txt and project_dir given relative to the receiver's working directory, the figure
  # saved after moving to another one.
  withr::local_dir(recv)
  proj <- share_project("PlotTrellis project", file.path("shared folder", "samples.txt"))
  expect_same_path(proj$cfg$data_dir, file.path(recv, "shared folder", "my data"))
  fig <- share_figure(proj)
  setwd(withr::local_tempdir())
  expect_identical(share_svg(fig), sender$md5)
})

test_that("the default example, copied to another folder with the original out of reach, draws the installed example's figure", {
  skip_if_not_installed("org.Mm.eg.db")
  local_example_chrom_sizes()
  # The GTF's first gene may lie beyond some of the cut-down data: plotgardener then warns of
  # an empty window, which draws as normal (still read from the files, which must be found).
  withCallingHandlers(warning = function(w) {
    if (grepl("(No range data|Not enough data within range) to plot", conditionMessage(w)))
      invokeRestart("muffleWarning")
  }, {
    # The sender is the installed example itself; the receiver a copy of it elsewhere. (The
    # installed files can't be moved away, but the copy's project holds no path into them.)
    installed <- share_project(tempfile("pt-share-proj-"), pt_example())
    gene <- share_first_gene(installed)
    expect_false(is.na(gene))
    ref <- share_svg(share_figure(installed, gene))
    st <- share_copy_example(file.path(tempfile("pt-share-"), "receiver home", "mouse data"), "default")
    withr::local_dir(withr::local_tempdir())
    proj <- share_project(tempfile("pt-share-proj-"), st)
    fig <- share_figure(proj, gene)
    expect_false(any(grepl(dirname(pt_example()), share_strings(proj, fig), fixed = TRUE)))
    expect_identical(share_svg(fig), ref)
  })
})

test_that("a project folder copied to another place on the same filesystem opens there and draws the same figure, with its own caches", {
  skip_if_not_installed("org.Hs.eg.db")
  local_example_chrom_sizes()
  sender <- share_sender()
  dest <- file.path(tempfile("pt-share-"), "receiver home")
  dir.create(dest, recursive = TRUE)
  stopifnot(file.copy(sender$dir, dest, recursive = TRUE))
  copy <- file.path(dest, basename(sender$dir))
  # The sender's project is gone; the data stays where project_settings.json says it is.
  local_moved_away(sender$dir)
  withr::local_dir(withr::local_tempdir())
  msgs <- character()
  proj <- withCallingHandlers(suppressWarnings(pt_project(copy)),
                              message = function(m) {
                                msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage")
                              })
  # The genome cache came with the folder: nothing rebuilt.
  expect_false(any(grepl("Building TxDb", msgs, fixed = TRUE)))
  expect_same_path(proj$render_cfg$cache_dir, file.path(copy, "cache"))
  expect_true(all(startsWith(normalizePath(unlist(proj$genome[c("txdb_path", "tx_name_lookup_path",
                                                                 "tx_meta_lookup_path", "chrom_sizes_path")])),
                             normalizePath(copy))))
  fig <- share_figure(proj)
  expect_false(any(grepl(sender$dir, share_strings(proj, fig), fixed = TRUE)))
  expect_identical(share_svg(fig), sender$md5)
})

test_that("a project folder someone else can only read opens and draws the same figure", {
  skip_on_os("windows")
  skip_if_not_installed("org.Hs.eg.db")
  local_example_chrom_sizes()
  sender <- share_sender()
  ro <- file.path(tempfile("pt-share-"), "read-only project")
  dir.create(dirname(ro))
  stopifnot(file.copy(sender$dir, dirname(ro), recursive = TRUE),
            file.rename(file.path(dirname(ro), basename(sender$dir)), ro))
  everything <- c(ro, list.files(ro, recursive = TRUE, full.names = TRUE, include.dirs = TRUE))
  Sys.chmod(everything, "0555")
  withr::defer(Sys.chmod(everything, "0755"))
  skip_if(file.access(ro, 2) == 0, "can't make a read-only folder here (running as root?)")
  .pt_close_txdbs()
  proj <- pt_project(ro)
  expect_identical(share_svg(share_figure(proj)), sender$md5)
})

test_that("a project folder copied elsewhere with its data moved stops, naming the missing GTF; re-imported from the samples.txt beside the data it draws the same figure", {
  skip_if_not_installed("org.Hs.eg.db")
  local_example_chrom_sizes()
  sender <- share_sender()
  home <- file.path(tempfile("pt-share-"), "receiver home")
  dir.create(file.path(home, "their data"), recursive = TRUE)
  stopifnot(file.copy(sender$dir, home, recursive = TRUE),
            all(file.copy(list.files(sender$data, full.names = TRUE), file.path(home, "their data"),
                          recursive = TRUE)))
  copy <- file.path(home, basename(sender$dir))
  st <- file.path(home, "their data", "samples.txt")
  local_moved_away(sender$data)
  local_moved_away(sender$dir)
  # The paths in the copied project are the sender's: the error names the missing file.
  expect_error(pt_project(copy), paste("Genome GTF not found:", file.path(sender$data, "annotation.gtf.gz")),
               fixed = TRUE)
  # The same folder from the receiver's samples.txt: it holds the sender's project, so it is
  # replaced only with overwrite = TRUE (which also says so).
  expect_error(pt_project(copy, st), "already holds a project, made from .*overwrite = TRUE")
  expect_message(proj <- suppressWarnings(pt_project(copy, st, overwrite = TRUE)), "Re-importing")
  expect_same_path(proj$cfg$data_dir, file.path(home, "their data"))
  expect_identical(share_svg(share_figure(proj)), sender$md5)
})

test_that("project_samples.txt carries the sender's app settings: it works from another folder on the same filesystem, and elsewhere once its two paths are edited", {
  skip_if_not_installed("org.Hs.eg.db")
  local_example_chrom_sizes()
  sender <- share_sender()
  # The sender changes a setting in the app (its Save writes both project files) and draws
  # the figure again: a different one.
  mine <- tempfile("pt-share-proj-")
  stopifnot(file.copy(sender$dir, dirname(mine), recursive = TRUE),
            file.rename(file.path(dirname(mine), basename(sender$dir)), mine))
  json <- file.path(mine, "project_settings.json")
  settings <- jsonlite::read_json(json)
  settings$fontsize <- 11L
  .pt_write_project_files(settings, json)
  custom <- share_svg(share_figure(share_project(mine)))
  expect_false(identical(custom, sender$md5))
  shared <- file.path(withr::local_tempdir(), "from a colleague.txt")
  stopifnot(file.copy(file.path(mine, "project_samples.txt"), shared))
  txt <- readLines(shared)
  expect_true(any(txt == paste("# data_root:", normalizePath(sender$data))))
  expect_true("# fontsize: 11" %in% txt)

  # Same filesystem: the absolute paths in it still lead to the data.
  withr::local_dir(withr::local_tempdir())
  proj <- share_project(tempfile("pt-share-proj-"), shared)
  expect_identical(share_svg(share_figure(proj)), custom)

  # Another computer: the sender's folders aren't there. The import stops on the data folder.
  home <- file.path(tempfile("pt-share-"), "receiver home")
  dir.create(file.path(home, "their data"), recursive = TRUE)
  stopifnot(all(file.copy(list.files(sender$data, full.names = TRUE), file.path(home, "their data"),
                          recursive = TRUE)))
  local_moved_away(sender$data)
  local_moved_away(mine)
  expect_error(pt_project(tempfile("pt-share-proj-"), shared), "Data root directory does not exist", fixed = TRUE)
  # The receiver puts it beside their copy of the data and makes its two paths relative.
  txt <- sub("^# data_root: .*", "# data_root: .", txt)
  txt <- sub("^# gtf_path: .*", "# gtf_path: annotation.gtf.gz", txt)
  writeLines(txt, file.path(home, "their data", "from a colleague.txt"))
  proj <- share_project(tempfile("pt-share-proj-"), file.path(home, "their data", "from a colleague.txt"))
  expect_identical(share_svg(share_figure(proj)), custom)
})

test_that("a downloaded R script run from another folder names each data file it can't find once the data moved", {
  skip_on_cran()
  skip_if_not_installed("org.Hs.eg.db")
  local_example_chrom_sizes()
  sender <- share_sender()
  run_dir <- withr::local_tempdir()
  pt_script(sender$fig, file.path(run_dir, "figure.R"))
  # Its paths are the sender's, absolute (the project's genome cache is still there; a missing
  # one is rebuilt from the GTF, see test-download-script.R).
  expect_true(any(grepl(file.path(sender$data, "atac"), readLines(file.path(run_dir, "figure.R")), fixed = TRUE)))
  local_moved_away(sender$data)
  withr::local_dir(run_dir)
  out <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", "figure.R"), stdout = TRUE, stderr = TRUE,
                 env = paste0("R_LIBS=", paste(.libPaths(), collapse = .Platform$path.sep)))
  failed <- grep("track failed", out, value = TRUE)
  expect_gt(length(failed), 0)
  expect_true(all(grepl("file not found", failed, fixed = TRUE)))
  expect_true(any(grepl("GM12878_ATAC_foldchange.bigWig", failed, fixed = TRUE)))
})
