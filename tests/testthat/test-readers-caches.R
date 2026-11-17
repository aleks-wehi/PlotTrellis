# Quant and junction readers and the per-gene caches (bugs 73, 83, 99, 100, 107, 115, 118).
# No network: UniProt reads are mocked or never reached.

# A gene with two transcripts, T1 and T2, without a real annotation.
local_two_transcripts <- function(env = parent.frame()) {
  local_mocked_bindings(gtf_transcript_structure = function(gene, anno)
    data.frame(transcript_id = c("T1.1", "T2.1"), transcript_name = c("G-201", "G-202"),
               tx_id_stripped = c("T1", "T2"), stringsAsFactors = FALSE), .env = env)
}
quant_file <- function(lines) {
  dir <- tempfile("quant-")
  dir.create(file.path(dir, "S1"), recursive = TRUE)
  f <- file.path(dir, "S1", "quant.sf")
  writeLines(c("Name\tLength\tEffectiveLength\tTPM\tNumReads", lines), f)
  f
}
one_sample <- function(path) list(list(sample_id = "S1", quant_path = path))

test_that("a nan or NA TPM for the gene stops with the file named, not a silent 0 (bug 73) (no fixture data)", {
  local_two_transcripts()
  for (bad in c("nan", "NA", "inf")) {
    f <- quant_file(c(paste0("T1.1\t1000\t900\t", bad, "\t5"), "T2.1\t1000\t900\t5\t5"))
    expect_error(filtered_transcripts("G", one_sample(f), NULL),
                 "TPM of T1.1 in S1/quant.sf is not a number", fixed = TRUE, info = bad)
  }
  # Another gene's bad value is not this gene's problem.
  f <- quant_file(c("T1.1\t1000\t900\t3\t5", "T2.1\t1000\t900\t5\t5", "T9.1\t1000\t900\tnan\t5"))
  tx <- filtered_transcripts("G", one_sample(f), NULL)
  expect_identical(attr(tx, "total_tpm"), 8)
})

test_that("a header-only quant.sf says it has no transcripts (bug 100) (no fixture data)", {
  local_two_transcripts()
  expect_error(filtered_transcripts("G", one_sample(quant_file(character(0))), NULL),
               "no transcripts in S1/quant.sf", fixed = TRUE)
  empty <- quant_file(character(0))
  file.create(empty)  # 0 bytes
  expect_error(filtered_transcripts("G", one_sample(empty), NULL), "could not read S1/quant.sf", fixed = TRUE)
})

test_that("in a downloaded script (BiocGenerics' do.call on the search path) a missing quant file says so plainly (bug 99) (no fixture data)", {
  skip_if_not_installed("BiocGenerics")
  local_two_transcripts()
  # A script defines the functions at top level, so do.call() is BiocGenerics' S4 generic there.
  as_in_script <- filtered_transcripts
  environment(as_in_script) <- list2env(list(do.call = BiocGenerics::do.call),
                                        parent = environment(filtered_transcripts))
  missing <- file.path(tempfile("quant-"), "S1", "quant.sf")
  err <- tryCatch(as_in_script("G", one_sample(missing), NULL), error = conditionMessage)
  expect_match(err, "^file not found: quant.sf")
  expect_no_match(err, "do.call", fixed = TRUE)
})

test_that("a BED12 junction spanning the whole window is not counted; a partial one is (bug 83) (no fixture data)", {
  f <- tempfile("synth-bed12-", fileext = ".bed")
  writeLines(c(
    "chr1\t2990\t7010\tinside\t12\t+\t2990\t7010\t0\t2\t10,10\t0,4010",
    "chr1\t500\t20010\tspans\t900\t+\t500\t20010\t0\t2\t10,10\t0,19500",   # intron 510-20000
    "chr1\t8990\t15010\tpartial\t40\t+\t8990\t15010\t0\t2\t10,10\t0,6010"  # intron 9000-15000
  ), f)
  j <- read_junction_file(f, "chr1", 1000L, 10000L)
  expect_identical(sort(j$nSupport), c(12L, 40L))
})

test_that("a corrupt cache file is a cache miss, not a failed track (bug 115) (no fixture data)", {
  junk <- tempfile(fileext = ".rds")
  writeLines("not an rds file", junk)
  expect_identical(.gene_display_transcripts_cache_load(junk), list())

  local_two_transcripts()
  cache_path <- tempfile(fileext = ".rds")
  writeLines("not an rds file", cache_path)
  res <- gene_display_transcripts("G", list(gtf_path = "g.gtf", tx_name_lookup_path = ""), list(), tempdir(),
                                  cache_path = cache_path)
  expect_identical(res$named, c("G-201", "G-202"))
  expect_length(readRDS(cache_path), 1)  # rewritten

  cache_dir <- tempfile("pd-")
  pd_file <- file.path(cache_dir, "mm10", "G.rds")
  dir.create(dirname(pd_file), recursive = TRUE)
  writeLines("not an rds file", pd_file)
  local_mocked_bindings(.uniprot_read_lines = function(url) character(0))
  res <- gene_protein_domains("G", NULL, "unused_tx_meta.rds", "org.Mm.eg.db", "mm10", cache_dir = cache_dir)
  expect_identical(res$status, "no_uniprot_entry")
  expect_identical(readRDS(pd_file)$status, "no_uniprot_entry")  # rewritten
  # Offline, the unreadable file is removed, not left to look cached.
  writeLines("not an rds file", pd_file)
  local_mocked_bindings(.uniprot_read_lines = function(url)
    stop(errorCondition("offline", class = "uniprot_unreachable")))
  res <- gene_protein_domains("G", NULL, "unused_tx_meta.rds", "org.Mm.eg.db", "mm10", cache_dir = cache_dir)
  expect_identical(res$status, "uniprot_unreachable")
  expect_false(file.exists(pd_file))
})

test_that("a local protein-domain error is not labelled a network issue (bug 107) (no fixture data)", {
  expect_no_match(.pd_status_text("error"), "network", fixed = TRUE)
  expect_no_match(.pd_status_text("error", panel = TRUE), "network", fixed = TRUE)
  expect_match(.pd_status_text("uniprot_unreachable", panel = TRUE), "network", fixed = TRUE)
  expect_match(.pd_status_text("uniprot_unreachable"), "could not be reached", fixed = TRUE)
})

test_that("an unsupported species says UniProt lookups are human and mouse only, not 'add one' (bug 118) (no fixture data)", {
  local_mocked_bindings(.uniprot_read_lines = function(url) stop("no network in this test"))
  for (call in list(quote(uniprot_search_url("Actb", "org.Rn.eg.db")),
                    quote(.uniprot_fetch_entry("Actb", "org.Rn.eg.db")),
                    quote(gene_protein_domains("Actb", NULL, "x.rds", "org.Rn.eg.db", "rn7", cache_dir = tempfile())))) {
    err <- tryCatch(eval(call), error = conditionMessage)
    expect_match(err, "human and mouse", fixed = TRUE)
    expect_match(err, "org.Rn.eg.db", fixed = TRUE)
    expect_no_match(err, "add one", fixed = TRUE)
    expect_no_match(err, "ORGDB_TAXON_ID", fixed = TRUE)
  }
})
