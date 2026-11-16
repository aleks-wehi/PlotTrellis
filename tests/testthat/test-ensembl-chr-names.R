# Bug 185: a TxDb built from an Ensembl GTF ("1", "MT") must answer in UCSC names ("chr1",
# "chrM") after it is saved and loaded again, including through select() and
# genes(filter =), which plotgardener's gene tracks query by chromosome. A GENCODE one
# ("chr1") must be left exactly as built.

# A GTF with one gene on each of `chroms`.
one_gene_gtf <- function(chroms) {
  gtf <- tempfile(fileext = ".gtf")
  writeLines(unlist(lapply(seq_along(chroms), function(i) {
    attrs <- sprintf('gene_id "ENSMUSG%011d"; transcript_id "ENSMUST%011d"; gene_name "G%d"; transcript_name "G%d-201";',
                     i, i, i, i)
    paste(chroms[i], "ensembl", c("gene", "transcript", "exon"), 1000, 2000, ".", "+", ".", attrs, sep = "\t")
  })), gtf)
  gtf
}

# The TxDb ensure_genome_ready() caches for that GTF, as a render loads it.
cached_txdb <- function(chroms) {
  skip_if_not_installed("org.Mm.eg.db")
  gtf <- one_gene_gtf(chroms)
  cache_dir <- tempfile()
  # Chromosome lengths already cached, so nothing is fetched from UCSC.
  saveRDS(c(chr1 = 100000L), .genome_cache_path(gtf, "_chrom_sizes.rds", cache_dir))
  genome <- suppressMessages(suppressWarnings(ensure_genome_ready(gtf, "org.Mm.eg.db", "mm10", cache_dir)))
  .pt_annotation(genome)$txdb
}

# A TxDb file saved with that GTF's own names (as a cache built before the fix is).
saved_txdb <- function(chroms) {
  path <- tempfile(fileext = ".sqlite")
  AnnotationDbi::saveDb(suppressWarnings(suppressMessages(txdbmaker::makeTxDbFromGFF(one_gene_gtf(chroms)))), path)
  path
}
chroms_of <- function(path) {
  con <- RSQLite::dbConnect(RSQLite::SQLite(), path)
  on.exit(RSQLite::dbDisconnect(con))
  sort(RSQLite::dbGetQuery(con, "SELECT chrom FROM chrominfo")$chrom)
}

test_that("an Ensembl-built TxDb answers in UCSC chromosome names once cached", {
  txdb <- cached_txdb(c("1", "MT", "JH584299.1"))
  expect_setequal(GenomeInfoDb::seqlevels(txdb), c("chr1", "chrM", "JH584299.1"))
  # What plotgardener's plotGenes()/plotTranscripts() ask for a region on chr1.
  expect_length(suppressMessages(GenomicFeatures::genes(txdb, filter = list(tx_chrom = "chr1"))), 1)
  sel <- suppressMessages(AnnotationDbi::select(txdb, keys = "ENSMUSG00000000001", keytype = "GENEID",
                                                columns = c("TXCHROM", "TXNAME")))
  expect_equal(sel$TXCHROM, "chr1")
  expect_equal(sel$TXNAME, "ENSMUST00000000001")
})

test_that("a GENCODE-built TxDb (chr names) keeps its names", {
  txdb <- cached_txdb(c("chr1", "chrM", "JH584299.1"))
  expect_setequal(GenomeInfoDb::seqlevels(txdb), c("chr1", "chrM", "JH584299.1"))
})

test_that("names with no UCSC equivalent are kept, without an error", {
  path <- saved_txdb(c("scaffold_1", "NC_000067.6"))
  expect_silent(.pt_txdb_ucsc_chroms(path))
  expect_equal(chroms_of(path), c("NC_000067.6", "scaffold_1"))
})

test_that("a TxDb file that can't be written is a warning, and opens as before", {
  skip_on_os("windows")
  path <- saved_txdb(c("1", "MT"))
  Sys.chmod(path, "0444")
  on.exit(Sys.chmod(path, "0644"))
  skip_if(file.access(path, 2) == 0, "running as a user who can write a read-only file")
  expect_warning(.pt_txdb_ucsc_chroms(path), "could not give the chromosomes")
  expect_equal(chroms_of(path), c("1", "MT"))
})

# Several processes opening one project at the same time: the first renames, the rest wait
# for its lock and then find the names done; one kept waiting too long gets a warning that
# says so, and the file is never renamed twice.
test_that("a renamer kept out by another process's lock warns that it is one, and is fine once the lock is gone", {
  skip_if_not_installed("txdbmaker")
  path <- saved_txdb(c("1", "X"))
  holder <- RSQLite::dbConnect(RSQLite::SQLite(), path)
  on.exit(try(RSQLite::dbDisconnect(holder), silent = TRUE), add = TRUE)
  RSQLite::dbExecute(holder, "BEGIN IMMEDIATE")
  expect_warning(res <- .pt_txdb_ucsc_chroms(path, busy_ms = 300), "another R process is renaming")
  expect_identical(res, "failed")
  expect_identical(chroms_of(path), c("1", "X"))
  RSQLite::dbExecute(holder, "ROLLBACK")
  expect_identical(.pt_txdb_ucsc_chroms(path), "renamed")
  expect_identical(chroms_of(path), c("chr1", "chrX"))
  expect_identical(.pt_txdb_ucsc_chroms(path), "unchanged")
})

test_that("a renamer that waits for another process finds the names done and leaves them", {
  skip_if_not_installed("txdbmaker")
  path <- saved_txdb(c("1", "X"))
  # Another R process takes the lock (and says so in `held`), renames the way the package
  # does, and commits 3 s later.
  held <- tempfile("held-")
  other <- callr::r_bg(function(path, held) {
    con <- RSQLite::dbConnect(RSQLite::SQLite(), path, synchronous = NULL)
    RSQLite::sqliteSetBusyHandler(con, 10000)
    RSQLite::dbExecute(con, "BEGIN IMMEDIATE")
    file.create(held)
    for (sql in c("UPDATE chrominfo SET chrom = ? WHERE chrom = ?",
                  "UPDATE transcript SET tx_chrom = ? WHERE tx_chrom = ?",
                  "UPDATE exon SET exon_chrom = ? WHERE exon_chrom = ?",
                  "UPDATE cds SET cds_chrom = ? WHERE cds_chrom = ?"))
      RSQLite::dbExecute(con, sql, params = list(c("chr1", "chrX"), c("1", "X")))
    Sys.sleep(3)
    RSQLite::dbExecute(con, "COMMIT")
    RSQLite::dbDisconnect(con)
    "committed"
  }, args = list(path, held))
  on.exit(try(other$kill(), silent = TRUE), add = TRUE)
  for (i in 1:200) { if (file.exists(held)) break; Sys.sleep(0.05) }   # until it holds the lock
  expect_true(file.exists(held))
  expect_no_warning(res <- .pt_txdb_ucsc_chroms(path, busy_ms = 10000))
  expect_identical(res, "done elsewhere")
  other$wait(10000)
  expect_identical(other$get_result(), "committed")
  expect_identical(chroms_of(path), c("chr1", "chrX"))
  con <- RSQLite::dbConnect(RSQLite::SQLite(), path); on.exit(RSQLite::dbDisconnect(con), add = TRUE)
  expect_false(any(grepl("chrchr", RSQLite::dbGetQuery(con, "SELECT tx_chrom FROM transcript")$tx_chrom)))
})
