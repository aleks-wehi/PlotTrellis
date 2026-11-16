# A GTF edited or replaced in place (same path) must not serve per-gene results computed from the
# old file: the display-transcript cache and the protein-domain cache key on the loaded genome's
# cache folder (named by the GTF's hash of path, size and modification time), which a changed GTF
# gets anew. A session still holding the old genome keeps reading and writing the old folder.

# Rewrites the GTF with different content and a later modification time, as a user replacing
# the file would (a rewrite within the same second could keep the old mtime).
replace_gtf <- function(path, text) {
  writeLines(text, path)
  Sys.setFileTime(path, Sys.time() + 5)
}

test_that("the display-transcript cache forgets a gene when its GTF changes in place", {
  gtf <- tempfile(fileext = ".gtf"); writeLines("# release 1", gtf)
  cache_dir <- tempfile("cache-")
  anno_for <- function() list(gtf_path = gtf, tx_name_lookup_path = .genome_cache_path(gtf, "_tx_names.rds", cache_dir))
  cache <- file.path(cache_dir, "gene_display_transcripts_cache.rds")
  structure_of <- function(names) function(gene, anno)
    data.frame(transcript_id = paste0(names, ".1"), transcript_name = names, tx_id_stripped = names)
  old_anno <- anno_for()
  local_mocked_bindings(gtf_transcript_structure = structure_of(c("G-201", "G-202")))
  first <- gene_display_transcripts("G", old_anno, list(), "", cache_path = cache)
  expect_identical(first$named, c("G-201", "G-202"))

  replace_gtf(gtf, "# release 2, with a renamed transcript")
  new_anno <- anno_for()
  expect_false(identical(new_anno$tx_name_lookup_path, old_anno$tx_name_lookup_path))
  local_mocked_bindings(gtf_transcript_structure = structure_of(c("G-201", "G-299")))
  second <- gene_display_transcripts("G", new_anno, list(), "", cache_path = cache)
  expect_identical(second$named, c("G-201", "G-299"))

  # Neither genome's entry is re-read: the new one hits, and a session still on the old genome
  # still gets the old names, not the new ones.
  local_mocked_bindings(gtf_transcript_structure = function(gene, anno) stop("read again"))
  expect_identical(gene_display_transcripts("G", new_anno, list(), "", cache_path = cache)$named, c("G-201", "G-299"))
  expect_identical(gene_display_transcripts("G", old_anno, list(), "", cache_path = cache)$named, c("G-201", "G-202"))
})

test_that("the protein-domain cache is per loaded genome, so a replaced GTF is looked up afresh", {
  gtf <- tempfile(fileext = ".gtf"); writeLines("# release 1", gtf)
  cache_dir <- tempfile("cache-")
  pd_dir <- file.path(cache_dir, "protein_domains")
  meta_for <- function() .genome_cache_path(gtf, "_tx_meta.rds", cache_dir)
  # A resolved result (the only kind that is cached), with UniProt canned.
  local_mocked_bindings(
    .uniprot_read_lines = function(url) {
      if (grepl("/search?", url, fixed = TRUE)) return("Q99999")
      if (grepl("\\.txt$", url)) return(c("FT   DOMAIN          1..10", "FT                   /note=\"Test domain\""))
      c(">sp|Q99999|TEST", strrep("M", 99))
    },
    .pd_gene_transcript_tags = function(...) data.frame(transcript_id = "T1.1", ccdsid = NA_character_))
  cds <- GenomicRanges::GRangesList(T1.1 = GenomicRanges::GRanges("chr1", IRanges::IRanges(101, 400), strand = "+"))
  lookup <- function(meta) gene_protein_domains("G", NULL, meta, "org.Mm.eg.db", "mm10",
                                                cache_dir = pd_dir, cds_by_tx = cds)
  old_meta <- meta_for()
  expect_identical(lookup(old_meta)$status, "ok")
  f1 <- .pd_cache_file(pd_dir, "mm10", old_meta, "G")
  expect_match(f1, file.path("mm10", basename(dirname(old_meta)), "G.rds"), fixed = TRUE)
  expect_true(file.exists(f1))
  expect_true(.pt_gi_pd_cached("G", "mm10", old_meta, pd_dir))

  replace_gtf(gtf, "# release 2, with a renamed transcript")
  new_meta <- meta_for()
  expect_false(.pt_gi_pd_cached("G", "mm10", new_meta, pd_dir))
  local_mocked_bindings(.uniprot_read_lines = function(url) stop("looked up again"))
  r2 <- lookup(new_meta)
  expect_identical(r2$status, "error")
  expect_match(r2$message, "looked up again", fixed = TRUE)
  # A session still on the old genome keeps its entry; nothing was written under the new one.
  expect_identical(lookup(old_meta)$status, "ok")
  expect_false(file.exists(.pd_cache_file(pd_dir, "mm10", new_meta, "G")))
})

test_that("placeholder lookup paths (no folder) keep the plain cache layout", {
  expect_identical(.pt_genome_id(NULL), "")
  expect_identical(.pt_genome_id("unused_tx_meta.rds"), "")
  expect_identical(.pt_genome_id("/proj/cache/genome/0123456789ab/x_tx_meta.rds"), "0123456789ab")
  expect_identical(.pd_cache_file("pd", "mm10", "unused_tx_meta.rds", "G"), file.path("pd", "mm10", "G.rds"))
})
