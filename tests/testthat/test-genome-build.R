# A project's first open builds its TxDb from the GTF (ensure_genome_ready(), R/genomes.R).
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
  gtf <- file.path(dirname(pt_example()), "gencode.vM25.example.gtf.gz")
  txdb_path <- .genome_cache_path(gtf, "_txdb.sqlite", cache_dir)
  # The save writes part of a file, then the process stops (killed, disk full).
  local_mocked_bindings(saveDb = function(x, file) {
    writeLines("half a TxDb", file)
    stop("save interrupted")
  }, .package = "AnnotationDbi")
  expect_error(suppressWarnings(suppressMessages(
    ensure_genome_ready(gtf, "org.Mm.eg.db", "mm10", cache_dir))), "save interrupted")
  expect_false(file.exists(txdb_path))
  expect_identical(list.files(dirname(txdb_path)), character(0))
})

test_that("a TxDb another process saved first is kept, and no temporary file is left (bug 212)", {
  skip_if_not_installed("org.Mm.eg.db")
  cache_dir <- tempfile("pt-cache-")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  gtf <- file.path(dirname(pt_example()), "gencode.vM25.example.gtf.gz")
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
  suppressWarnings(suppressMessages(ensure_genome_ready(gtf, "org.Mm.eg.db", "mm10", cache_dir)))
  expect_identical(txdb_marker(txdb_path), "the other process")
  expect_false(any(grepl("_txdb\\.sqlite-", list.files(dirname(txdb_path)))))
})
