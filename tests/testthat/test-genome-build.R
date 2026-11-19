# A project's first open builds its TxDb from the GTF (ensure_genome_ready(), R/genome-main.R).
# That build gives no routine warnings that options(warn = 2) turns into errors, and never
# leaves a half-written TxDb at the cache path. Uses the bundled example, so needs no fixture
# data (org.Mm.eg.db only).

# The marker the "other process" below writes into its TxDb.
txdb_marker <- function(path) {
  con <- RSQLite::dbConnect(RSQLite::SQLite(), path)
  on.exit(RSQLite::dbDisconnect(con))
  RSQLite::dbGetQuery(con, "SELECT value FROM metadata WHERE name = 'saved by'")$value
}

test_that("a cold project opens under options(warn = 2) (bug 129)", {
  skip_if_not_installed("org.Mm.eg.db")
  proj <- tempfile("pt-warn2-")
  on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  old <- options(warn = 2)
  on.exit(options(old), add = TRUE)
  expect_s3_class(suppressMessages(pt_project(proj, pt_example())), "pt_project")
})

# The test above only sees the installed txdbmaker's wording; R-universe's newer txdbmaker
# says "genome assembly" where older ones say "genome version information".
test_that("the TxDb build drops both wordings of the missing-genome warning, and only those", {
  build_warning <- NULL
  local_mocked_bindings(makeTxDbFromGFF = function(...) {
    warning(build_warning)
    "a TxDb"
  }, .package = "txdbmaker")
  for (msg in c("genome version information is not available for this TxDb object",
                "genome assembly is not available for this TxDb object",
                paste("The \"phase\" metadata column contains non-NA values for features",
                      "of type Selenocysteine, stop_codon. This information was ignored."))) {
    build_warning <- msg
    expect_no_warning(expect_identical(.pt_make_txdb("x.gtf"), "a TxDb"))
  }
  build_warning <- "some CDS phases are missing or not between 0 and 2"
  expect_warning(.pt_make_txdb("x.gtf"), "some CDS phases", fixed = TRUE)
})

test_that("a TxDb build stopped part-way leaves no file at the cache path (bug 212)", {
  skip_if_not_installed("org.Mm.eg.db")
  cache_dir <- tempfile("pt-cache-")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  gtf <- file.path(dirname(pt_example()), "annotation.gtf.gz")
  txdb_path <- .genome_cache_path(gtf, "_txdb.sqlite", cache_dir)
  # The save writes part of a file, then the process stops (killed, disk full).
  local_mocked_bindings(saveDb = function(x, file) {
    writeLines("half a TxDb", file)
    stop("save interrupted")
  }, .package = "AnnotationDbi")
  expect_error(suppressWarnings(suppressMessages(
    ensure_genome_ready(gtf, "org.Mm.eg.db", "mm39", cache_dir))), "save interrupted")
  expect_false(file.exists(txdb_path))
  expect_identical(list.files(dirname(txdb_path)), character(0))
})

test_that("a TxDb another process saved first is kept, and no temporary file is left (bug 212)", {
  skip_if_not_installed("org.Mm.eg.db")
  cache_dir <- tempfile("pt-cache-")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  gtf <- file.path(dirname(pt_example()), "annotation.gtf.gz")
  txdb_path <- .genome_cache_path(gtf, "_txdb.sqlite", cache_dir)
  real_save <- AnnotationDbi::saveDb
  # While this process saves, another one finishes its own build of the same TxDb.
  local_mocked_bindings(saveDb = function(x, file) {
    real_save(x, txdb_path)
    con <- RSQLite::dbConnect(RSQLite::SQLite(), txdb_path)
    on.exit(RSQLite::dbDisconnect(con))
    RSQLite::dbExecute(con, "INSERT INTO metadata VALUES ('saved by', 'the other process')")
    real_save(x, file)
  }, .package = "AnnotationDbi")
  suppressWarnings(suppressMessages(ensure_genome_ready(gtf, "org.Mm.eg.db", "mm39", cache_dir)))
  expect_identical(txdb_marker(txdb_path), "the other process")
  expect_false(any(grepl("_txdb\\.sqlite-", list.files(dirname(txdb_path)))))
})

# A whole-genome GTF (over 10 MB) is built in a separate R process (bug 272); the bundled
# example is too small for that, so these tests ask for it.
test_that("a whole-genome build runs in a separate R process, with the same caches and messages (bug 272)", {
  skip_on_cran()  # starts R processes, as test-render-worker.R's tests
  skip_if_not_installed("org.Mm.eg.db")
  gtf <- file.path(dirname(pt_example()), "annotation.gtf.gz")
  here <- tempfile("pt-cache-")
  child <- tempfile("pt-cache-")
  on.exit(unlink(c(here, child), recursive = TRUE), add = TRUE)
  expect_false(.pt_genome_build_in_child(gtf, child))  # 100 kB: built here
  msgs_here <- character(0)
  g_here <- withCallingHandlers(
    suppressWarnings(ensure_genome_ready(gtf, "org.Mm.eg.db", "mm10", here)),
    message = function(m) {
      msgs_here <<- c(msgs_here, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  local_mocked_bindings(.pt_genome_build_in_child = function(gtf_path, cache_dir) TRUE)
  # A build in this process would stop here.
  local_mocked_bindings(makeTxDbFromGFF = function(...) stop("built in the test's own process"),
                        .package = "txdbmaker")
  msgs_child <- character(0)
  g_child <- withCallingHandlers(
    suppressWarnings(ensure_genome_ready(gtf, "org.Mm.eg.db", "mm10", child)),
    message = function(m) {
      msgs_child <<- c(msgs_child, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  expect_true(g_child$rebuilt)
  # the same messages; the separate process reads the gene names before it builds the TxDb
  expect_identical(sort(msgs_child), sort(msgs_here))
  expect_length(msgs_child, 2)
  expect_match(msgs_child[1], "^Reading the gene names in annotation\\.gtf\\.gz")
  expect_match(msgs_child[2], "^Building TxDb from .*example.annotation\\.gtf\\.gz")
  # the same cache files, the lookups byte for byte, the TxDb table for table
  rel <- function(dir) sort(list.files(dir, recursive = TRUE))
  expect_identical(rel(child), rel(here))
  for (k in c("tx_name_lookup_path", "tx_meta_lookup_path"))
    expect_identical(unname(tools::md5sum(g_child[[k]])), unname(tools::md5sum(g_here[[k]])))
  expect_identical(unname(tools::md5sum(.pt_gene_names_file(g_child$txdb_path))),
                   unname(tools::md5sum(.pt_gene_names_file(g_here$txdb_path))))
  tables <- function(path) {
    con <- RSQLite::dbConnect(RSQLite::SQLite(), path)
    on.exit(RSQLite::dbDisconnect(con))
    lapply(c("chrominfo", "transcript", "exon", "cds", "gene", "splicing"),
           function(t) RSQLite::dbReadTable(con, t))
  }
  expect_identical(tables(g_child$txdb_path), tables(g_here$txdb_path))
  # A second open finds the caches (built by either process) and builds nothing.
  expect_false(suppressMessages(ensure_genome_ready(gtf, "org.Mm.eg.db", "mm10", child))$rebuilt)
})

test_that("an error in the separate build process stops with the same message (bug 272)", {
  skip_on_cran()  # starts R processes, as test-render-worker.R's tests
  skip_if_not_installed("org.Mm.eg.db")
  row <- function(feature, attrs, start = 1, end = 100)
    paste("chr1", "src", feature, start, end, ".", "+", ".", attrs, sep = "\t")
  gtf <- tempfile(fileext = ".gtf")
  writeLines(c(row("gene", 'gene_id "Elf5"; gene_name "Elf5";'),
               row("transcript", 'gene_id "Elf5"; transcript_id "NM_1";'),
               row("exon", 'gene_id "Elf5"; transcript_id "NM_1"; exon_number "1";')), gtf)
  dirs <- c(tempfile("pt-cache-"), tempfile("pt-cache-"))
  on.exit(unlink(c(gtf, dirs), recursive = TRUE), add = TRUE)
  err <- function(dir) tryCatch(suppressWarnings(suppressMessages(
    ensure_genome_ready(gtf, "org.Mm.eg.db", "mm10", dir))), error = conditionMessage)
  here <- err(dirs[1])
  expect_match(here, "don't look like Ensembl IDs", fixed = TRUE)
  local_mocked_bindings(.pt_genome_build_in_child = function(gtf_path, cache_dir) TRUE)
  expect_identical(err(dirs[2]), here)
})

test_that("a build process that dies is an error saying so, and leaves no cache file (bug 272)", {
  skip_on_cran()  # starts R processes, as test-render-worker.R's tests
  skip_if_not_installed("org.Mm.eg.db")
  skip_if_not(file.exists(file.path("/proc", Sys.getpid(), "task")), "needs Linux /proc")
  gtf <- file.path(dirname(pt_example()), "annotation.gtf.gz")
  cache_dir <- tempfile("pt-cache-")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  kids <- function() {
    f <- sprintf("/proc/%d/task/%d/children", Sys.getpid(), Sys.getpid())
    scan(f, integer(), quiet = TRUE)
  }
  # Another process kills the build process (SIGKILL, as the out-of-memory killer does) a
  # second after it starts.
  killer <- callr::r_bg(function(parent, before) {
    f <- sprintf("/proc/%d/task/%d/children", parent, parent)
    for (i in 1:600) {
      now <- scan(f, integer(), quiet = TRUE)
      new <- setdiff(now, c(before, Sys.getpid()))
      if (length(new)) { Sys.sleep(1); tools::pskill(new, 9L); return(new) }
      Sys.sleep(0.1)
    }
  }, list(Sys.getpid(), kids()))
  on.exit(killer$kill(), add = TRUE)
  local_mocked_bindings(.pt_genome_build_in_child = function(gtf_path, cache_dir) TRUE)
  expect_error(suppressMessages(ensure_genome_ready(gtf, "org.Mm.eg.db", "mm10", cache_dir)),
               "the R process building it stopped before it finished (out of memory, or killed). The first open of a whole-genome GTF needs an R session with at least 10 GB of memory: start one that size and try again. It is built once per project; 8 GB is enough after that.",
               fixed = TRUE)
  killer$wait(5000)
  expect_length(killer$get_result(), 1)
  paths <- unlist(.pt_genome_paths(gtf, cache_dir))
  expect_false(any(file.exists(paths)))
})

test_that("a genome cache already built (by this or an earlier version) is read with no rebuild and no build process (bug 272)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "fixture data not available")
  # The fixture genomes' caches (helper-fixtures.R) were built by earlier versions; whole-genome GTFs.
  for (g in GENOMES) {
    expect_gt(file.size(g$gtf_path), 1e7)
    expect_false(.pt_genome_build_in_child(g$gtf_path, TEST_GENOME_CACHE))
    files <- unlist(.pt_genome_paths(g$gtf_path, TEST_GENOME_CACHE))
    before <- file.mtime(files)
    local_mocked_bindings(.pt_genome_files_child = function(...) stop("started a build process"))
    again <- ensure_genome_ready(g$gtf_path, g$orgdb, g$genome_label, TEST_GENOME_CACHE)
    expect_false(again$rebuilt)
    expect_identical(file.mtime(files), before)
  }
})
