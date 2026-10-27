# Phase I.3 step 28 (stage B): per-setting effect tests for the genome
# settings (gtf_path, genome_label, orgdb). Each test runs the same consumer
# at two values of one setting and asserts a specific difference.
#
# NO FIXTURE DATA NEEDED: tiny GTFs / TxDbs are written to tempdir(), and
# build_chrom_sizes() (a UCSC download) is mocked, so these also run on
# R-universe. The orgdb render test needs org.Mm.eg.db and org.Hs.eg.db
# (Suggests) and skips without them.

# One gene with one two-exon transcript; `offset` shifts it along chr1 and
# `tx_name` names the transcript, so two GTFs differ in content, not just path.
tiny_gtf <- function(offset, tx_name) {
  dir <- tempfile("synth-gtf-"); dir.create(dir)
  path <- file.path(dir, "synth.gtf")
  g <- 'gene_id "ENSMUSG00000000001.1"; gene_name "Synth";'
  t <- sprintf('%s transcript_id "ENSMUST00000000001.1"; transcript_name "%s"; transcript_type "protein_coding";', g, tx_name)
  row <- function(type, s, e, attr) paste("chr1", "synth", type, s + offset, e + offset, ".", "+", ".", attr, sep = "\t")
  writeLines(c(row("gene", 2001, 8000, g), row("transcript", 2001, 8000, t),
               row("exon", 2001, 3000, paste(t, 'exon_number 1;')),
               row("exon", 7001, 8000, paste(t, 'exon_number 2;'))), path)
  path
}

ready <- function(gtf_path, genome_label = "synthetic") {
  labels <- character(0)
  local_mocked_bindings(build_chrom_sizes = function(genome_label, out_path) {
    labels <<- c(labels, genome_label)
    saveRDS(c(chr1 = 100000L), out_path)
  })
  suppressMessages(suppressWarnings(utils::capture.output(
    out <- ensure_genome_ready(gtf_path = gtf_path, orgdb = "org.Mm.eg.db", genome_label = genome_label)
  )))
  out$chrom_sizes_label <- labels
  out
}

test_that("gtf_path decides which annotation the TxDb and transcript names come from (no fixture data)", {
  skip_if_not_installed("org.Mm.eg.db")
  skip_if_not_installed("txdbmaker")
  a <- ready(tiny_gtf(0, "Synth-201"))
  b <- ready(tiny_gtf(10000, "Synth-202"))
  expect_false(identical(a$txdb_path, b$txdb_path))
  gene_start <- function(g) {
    txdb <- AnnotationDbi::loadDb(g$txdb_path)
    GenomicRanges::start(GenomicFeatures::genes(txdb))
  }
  expect_identical(gene_start(a), 2001L)
  expect_identical(gene_start(b), 12001L)
  expect_identical(unname(readRDS(a$tx_name_lookup_path)["ENSMUST00000000001"]), "Synth-201")
  expect_identical(unname(readRDS(b$tx_name_lookup_path)["ENSMUST00000000001"]), "Synth-202")
})

test_that("genome_label picks the chromosome-size build and scopes the protein-domain cache (no fixture data)", {
  skip_if_not_installed("org.Mm.eg.db")
  skip_if_not_installed("txdbmaker")
  # Chrom sizes are fetched for the label given (UCSC lookup, mocked here) ...
  expect_identical(ready(tiny_gtf(0, "Synth-201"), "mm10")$chrom_sizes_label, "mm10")
  expect_identical(ready(tiny_gtf(0, "Synth-201"), "mm39")$chrom_sizes_label, "mm39")
  # ... and a blank label falls back to the orgdb name rather than failing.
  blank <- ready(tiny_gtf(0, "Synth-201"), "")
  expect_identical(blank$genome_label, "org.Mm.eg.db")
  expect_identical(blank$chrom_sizes_label, "org.Mm.eg.db")

  # Protein-domain results are cached per genome_label: a cached mm10 entry
  # must not be served for mm39.
  cache_dir <- tempfile("pd-cache-")
  dir.create(file.path(cache_dir, "mm10"), recursive = TRUE)
  saveRDS(list(status = "ok", gene = "Testgene", marker = "mm10 entry"), file.path(cache_dir, "mm10", "Testgene.rds"))
  local_mocked_bindings(.uniprot_read_lines = function(url) stop(structure(
    class = c("uniprot_unreachable", "error", "condition"), list(message = "offline (test)", call = NULL))))
  hit  <- gene_protein_domains("Testgene", NULL, "unused.gtf", "org.Mm.eg.db", "mm10", cache_dir = cache_dir)
  miss <- gene_protein_domains("Testgene", NULL, "unused.gtf", "org.Mm.eg.db", "mm39", cache_dir = cache_dir)
  expect_identical(hit$marker, "mm10 entry")
  expect_null(miss$marker)
  expect_identical(miss$status, "uniprot_unreachable")
})

test_that("orgdb decides the gene symbols drawn on the genes track (no fixture data)", {
  skip_if_not_installed("org.Mm.eg.db")
  skip_if_not_installed("org.Hs.eg.db")
  skip_if_not_installed("txdbmaker")
  # A one-gene TxDb whose gene id is mouse Elf5's real Ensembl id.
  elf5 <- "ENSMUSG00000027186"
  expect_identical(unname(AnnotationDbi::mapIds(org.Mm.eg.db::org.Mm.eg.db, elf5, "SYMBOL", "ENSEMBL")), "Elf5")
  txdb <- suppressWarnings(txdbmaker::makeTxDb(
    transcripts = data.frame(tx_id = 1L, tx_name = "ENSMUST00000000001", tx_chrom = "chr2",
                             tx_strand = "+", tx_start = 2000L, tx_end = 8000L),
    splicings = data.frame(tx_id = 1L, exon_rank = 1:2, exon_start = c(2000L, 7000L), exon_end = c(3000L, 8000L)),
    genes = data.frame(tx_id = 1L, gene_id = elf5),
    chrominfo = data.frame(chrom = "chr2", length = 100000L, is_circular = FALSE)
  ))
  txdb_path <- tempfile(fileext = ".sqlite")
  AnnotationDbi::saveDb(txdb, txdb_path)
  region <- list(chrom = "chr2", start = 1000L, end = 10000L, gene = NA_character_)
  svg_for <- function(orgdb) {
    f <- tempfile(fileext = ".svg")
    on.exit(unlink(f), add = TRUE)
    dev <- NULL
    on.exit(if (!is.null(dev) && dev %in% grDevices::dev.list()) grDevices::dev.off(dev), add = TRUE, after = FALSE)
    settings <- list(txdb_path = txdb_path, orgdb_pkg = orgdb, genome_label = "synthetic", page_w = 9,
                     font = "DejaVu Sans", fontsize = 8L, chr_fmt = "dot")
    utils::capture.output(render_all(list(list(viz = "genes", height = NULL, data = NULL)), region, settings,
                                     device_fn = function(w, h) {
                                       svglite::svglite(f, width = w, height = h)
                                       dev <<- grDevices::dev.cur()
                                     }))
    grDevices::dev.off(dev)
    dev <- NULL
    paste(readLines(f), collapse = "\n")
  }
  expect_true(grepl(">Elf5</text>", svg_for("org.Mm.eg.db"), fixed = TRUE))
  # With the human OrgDb the mouse id has no symbol (today the genes track
  # reports an error in its slot), so Elf5 is not drawn.
  expect_false(grepl(">Elf5</text>", svg_for("org.Hs.eg.db"), fixed = TRUE))
})

test_that("orgdb sets the UniProt organism, uniprot_reviewed_only the Swiss-Prot filter (no fixture data)", {
  # The Gene info panel's link-out, R/gene_info_panel.R gr_extras_ui. The composer passes
  # both settings straight through; that wiring is the composer's, not tested here.
  q <- function(...) utils::URLdecode(uniprot_search_url("Pax5", ...))
  expect_match(q("org.Mm.eg.db"), "organism_id:10090", fixed = TRUE)
  expect_match(q("org.Hs.eg.db"), "organism_id:9606", fixed = TRUE)
  expect_false(grepl("reviewed", q("org.Mm.eg.db", FALSE), fixed = TRUE))
  expect_match(q("org.Mm.eg.db", TRUE), "AND (reviewed:true)", fixed = TRUE)
})
