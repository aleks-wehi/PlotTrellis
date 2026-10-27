# Phase I.1 step 5: UniProt failure modes. A network failure must never be
# cached (one offline session used to cache "no reviewed entry" for every
# gene it touched), must read differently from a real "no reviewed entry"
# answer, and must not stop a later successful lookup from caching normally.
#
# Network failure is real, not faked: .uniprot_read_lines() is wrapped so its
# rest.uniprot.org URL is rewritten to 127.0.0.1:1 (connection refused,
# instant, no internet needed), going through the real readLines() error
# path. Success is mocked with canned UniProt text. No fixture data needed.

# Captured before any test mocks it.
.real_uniprot_read_lines <- .uniprot_read_lines
.refused <- function(url) {
  suppressWarnings(.real_uniprot_read_lines(sub("https://rest.uniprot.org", "http://127.0.0.1:1", url, fixed = TRUE)))
}

# Canned UniProt responses: search -> one accession; .txt record with one
# CCDS xref and one DOMAIN; .fasta of 99 aa.
.canned <- function(search = "Q99999") {
  function(url) {
    if (grepl("/search?", url, fixed = TRUE)) return(search)
    if (grepl("\\.txt$", url)) return(c(
      "ID   TEST_MOUSE              Reviewed;          99 AA.",
      "DR   CCDS; CCDS1.1; -.",
      "DR   Ensembl; ENSMUST00000000001.1; ENSMUSP00000000001.1; ENSMUSG00000000001.1.",
      "FT   DOMAIN          1..10",
      "FT                   /note=\"Test domain\"",
      "SQ   SEQUENCE   99 AA;  11111 MW;  0000000000000000 CRC64;"))
    if (grepl("\\.fasta$", url)) return(c(">sp|Q99999|TEST_MOUSE", strrep("M", 99)))
    stop("unexpected URL in test: ", url)
  }
}

# One-transcript TxDb, CDS 300 nt (99 aa + stop), so the "ok" path can run
# end to end without fixture data.
.tiny_txdb <- function() {
  gr <- GenomicRanges::GRanges("chr1", IRanges::IRanges(rep(101, 4), rep(400, 4)), strand = "+",
    type = c("gene", "transcript", "exon", "CDS"), gene_id = "ENSMUSG00000000001",
    transcript_id = c(NA, "ENSMUST00000000001.1", "ENSMUST00000000001.1", "ENSMUST00000000001.1"),
    phase = c(NA, NA, NA, 0L))
  suppressWarnings(txdbmaker::makeTxDbFromGRanges(gr))
}

.reset_backoff <- function() options(PlotTrellis.uniprot_down_until = NULL)

test_that(".uniprot_read_lines() raises uniprot_unreachable on a refused connection", {
  .reset_backoff()
  expect_error(suppressWarnings(.uniprot_read_lines("http://127.0.0.1:1/x")),
               class = "uniprot_unreachable")
  # and restores the caller's timeout option
  old <- getOption("timeout")
  try(suppressWarnings(.uniprot_read_lines("http://127.0.0.1:1/x")), silent = TRUE)
  expect_identical(getOption("timeout"), old)
})

test_that(".uniprot_read_lines() gives up after ~10s on a host that never answers", {
  skip_on_cran()
  .reset_backoff()
  # A listening socket that never accept()s or replies: the TCP connect
  # succeeds via the kernel backlog, then the read hangs until the timeout.
  srv <- NULL
  for (port in sample(20000:60000, 5)) {
    srv <- tryCatch(serverSocket(port), error = function(e) NULL)
    if (!is.null(srv)) break
  }
  skip_if(is.null(srv), "could not open a local listening socket")
  on.exit(close(srv), add = TRUE)
  elapsed <- system.time(
    expect_error(suppressWarnings(.uniprot_read_lines(sprintf("http://127.0.0.1:%d/x", port))),
                 class = "uniprot_unreachable")
  )[["elapsed"]]
  expect_gt(elapsed, 8)
  expect_lt(elapsed, 20)
})

test_that("gene_protein_domains(): network failure is not cached, a later success is", {
  cache_dir <- tempfile()
  dir.create(cache_dir)
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  txdb <- .tiny_txdb()
  cache_file <- file.path(cache_dir, "test_genome", "Testgene.rds")
  local_mocked_bindings(.pd_gene_transcript_tags = function(...)
    data.frame(transcript_id = "ENSMUST00000000001.1", ccdsid = "CCDS1.1", stringsAsFactors = FALSE))

  # Search itself unreachable
  local_mocked_bindings(.uniprot_read_lines = .refused)
  res <- gene_protein_domains("Testgene", txdb, "unused.gtf", "org.Mm.eg.db", "test_genome", cache_dir = cache_dir)
  expect_identical(res$status, "uniprot_unreachable")
  expect_false(file.exists(cache_file))

  # Search succeeds, record fetch unreachable
  canned <- .canned()
  local_mocked_bindings(.uniprot_read_lines = function(url)
    if (grepl("/search?", url, fixed = TRUE)) canned(url) else .refused(url))
  res <- gene_protein_domains("Testgene", txdb, "unused.gtf", "org.Mm.eg.db", "test_genome", cache_dir = cache_dir)
  expect_identical(res$status, "uniprot_unreachable")
  expect_false(file.exists(cache_file))

  # Real "no reviewed entry" is a distinct status
  local_mocked_bindings(.uniprot_read_lines = .canned(search = character(0)))
  res <- gene_protein_domains("Testgene", txdb, "unused.gtf", "org.Mm.eg.db", "test_genome", cache_dir = cache_dir)
  expect_identical(res$status, "no_uniprot_entry")

  # UniProt back: full success, cached
  local_mocked_bindings(.uniprot_read_lines = .canned())
  res <- gene_protein_domains("Testgene", txdb, "unused.gtf", "org.Mm.eg.db", "test_genome", cache_dir = cache_dir)
  expect_identical(res$status, "ok")
  expect_true(file.exists(cache_file))
  expect_identical(readRDS(cache_file)$status, "ok")
})

test_that("protein-domains track shows a different message for unreachable vs no entry", {
  labels <- character(0)
  local_mocked_bindings(plotText = function(label, ...) labels <<- c(labels, label))
  ctx <- list(region_params = list(gene = "Testgene"), track_x = 0, track_w = 5, y = 0, height = 2,
              fontsize = 10, font = "sans", fontcol = "black")
  msg_for <- function(status) {
    labels <<- character(0)
    local_mocked_bindings(gene_protein_domains = function(...) list(status = status))
    render_protein_domains(list(), list(), ctx)
    labels[length(labels)]
  }
  unreachable <- msg_for("uniprot_unreachable")
  no_entry <- msg_for("no_uniprot_entry")
  expect_match(unreachable, "UniProt could not be reached", fixed = TRUE)
  expect_identical(no_entry, "No reviewed UniProt entry found for this gene")
})

test_that("uniprot_canonical_match(): network failure is not cached, a later success is", {
  cache_path <- tempfile(fileext = ".rds")
  on.exit(unlink(cache_path), add = TRUE)
  local_mocked_bindings(UNIPROT_CANONICAL_CACHE_PATH = cache_path)
  tx <- c("ENSMUST00000000001.1", "ENSMUST00000000002.1")

  local_mocked_bindings(.uniprot_read_lines = .refused)
  res <- uniprot_canonical_match("Testgene", tx, list(), orgdb = "org.Mm.eg.db")
  expect_identical(nrow(res), 2L)
  expect_true(all(is.na(res$uniprot_accession)))
  expect_false(file.exists(cache_path))

  # Search succeeds, record fetch unreachable: used to cache an
  # accession-only entry with no protein length, forever
  canned <- .canned()
  local_mocked_bindings(.uniprot_read_lines = function(url)
    if (grepl("/search?", url, fixed = TRUE)) canned(url) else .refused(url))
  res <- uniprot_canonical_match("Testgene", tx, list(), orgdb = "org.Mm.eg.db")
  expect_true(all(is.na(res$uniprot_accession)))
  expect_false(file.exists(cache_path))

  local_mocked_bindings(.uniprot_read_lines = .canned())
  res <- uniprot_canonical_match("Testgene", tx, list(), orgdb = "org.Mm.eg.db")
  expect_identical(res$uniprot_accession, rep("Q99999", 2))
  expect_identical(res$uniprot_match, c(TRUE, FALSE))
  expect_true(file.exists(cache_path))
  expect_identical(readRDS(cache_path)[["org.Mm.eg.db|Testgene"]]$accession, "Q99999")
})

test_that("after a failure, UniProt calls fail fast until the backoff window passes", {
  .reset_backoff()
  on.exit(.reset_backoff(), add = TRUE)
  expect_error(.uniprot_read_lines("http://127.0.0.1:1/x"), "could not be reached: ", class = "uniprot_unreachable")
  # inside the window: short-circuits without touching the network
  expect_error(.uniprot_read_lines("http://127.0.0.1:1/y"), "skipping until", class = "uniprot_unreachable")
  # window passed: a real attempt again
  options(PlotTrellis.uniprot_down_until = Sys.time() - 1)
  expect_error(.uniprot_read_lines("http://127.0.0.1:1/z"), "could not be reached: ", class = "uniprot_unreachable")
})
