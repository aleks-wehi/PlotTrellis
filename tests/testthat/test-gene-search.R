# Gene search accepts a symbol in any
# case or an Ensembl gene ID (with or without version), and always resolves
# to the OrgDb's canonical symbol. Transcript IDs, Entrez IDs and aliases are
# deliberately not accepted. The resolver tests need only the OrgDb packages.

test_that("symbols resolve exactly, then case-insensitively when unambiguous", {
  expect_equal(resolve_gene_query("Pax5", "org.Mm.eg.db")$symbol, "Pax5")
  expect_equal(resolve_gene_query("PAX5", "org.Mm.eg.db")$symbol, "Pax5")
  expect_equal(resolve_gene_query(" pax5 ", "org.Mm.eg.db")$symbol, "Pax5")
  expect_equal(resolve_gene_query("ptprc", "org.Hs.eg.db")$symbol, "PTPRC")
  # "T" and "t" are both real mouse symbols: exact wins, other cases are ambiguous
  expect_equal(resolve_gene_query("T", "org.Mm.eg.db")$symbol, "T")
  amb <- resolve_gene_query("gS", "org.Mm.eg.db")
  expect_true(is.na(amb$symbol))
  expect_match(amb$error, "differ only in case")
})

test_that("Ensembl gene IDs resolve to their symbol, with or without a version", {
  for (q in c("ENSMUSG00000014030", "ENSMUSG00000014030.15", "ensmusg00000014030")) {
    r <- resolve_gene_query(q, "org.Mm.eg.db")
    expect_equal(r$symbol, "Pax5")
    expect_equal(r$ensembl, "ENSMUSG00000014030")
  }
  expect_equal(resolve_gene_query("ENSG00000081237", "org.Hs.eg.db")$symbol, "PTPRC")
  # A GENCODE PAR_Y copy is looked up as its chrX gene
  par <- resolve_gene_query("ENSG00000002586.20_PAR_Y", "org.Hs.eg.db")
  expect_equal(par$symbol, "CD99")
  expect_equal(par$ensembl, "ENSG00000002586")
})

# Which symbols an Ensembl ID has depends on the org.Mm.eg.db build: one maps
# ENSMUSG00000093931 to Amy2a1 and Amy2a3, the one R-universe checks against
# to Amy2a3 alone. So the expected symbols are read from the installed OrgDb's
# own tables (not through select(), which the resolver uses), and the several-
# symbol case runs on an ID that has several in that build.
test_that("an Ensembl ID with several symbols takes the first and names the rest", {
  map <- merge(AnnotationDbi::toTable(org.Mm.eg.db::org.Mm.egENSEMBL),
               AnnotationDbi::toTable(org.Mm.eg.db::org.Mm.egSYMBOL), by = "gene_id")
  map <- unique(map[, c("ensembl_id", "symbol")])
  takes_first_names_rest <- function(id) {
    want <- sort(map$symbol[map$ensembl_id == id], method = "radix")  # C-locale order
    r <- resolve_gene_query(id, "org.Mm.eg.db")
    expect_identical(r$symbol, want[1], label = paste(id, "symbol"))
    expect_identical(r$also, want[-1], label = paste(id, "also"))
  }
  several <- sort(unique(map$ensembl_id[duplicated(map$ensembl_id)]), method = "radix")
  expect_gt(length(several), 0)
  takes_first_names_rest(several[1])
  # The ID this test was written for, while the OrgDb still knows it
  if ("ENSMUSG00000093931" %in% map$ensembl_id) takes_first_names_rest("ENSMUSG00000093931")
  expect_length(resolve_gene_query("ENSMUSG00000014030", "org.Mm.eg.db")$also, 0)
})

test_that("unsupported or unknown queries resolve to nothing", {
  # A transcript ID, an Entrez ID, an unknown ID shape and a non-gene word
  for (q in c("ENSMUST00000000001", "18507", "notagene", "")) {
    r <- resolve_gene_query(q, "org.Mm.eg.db")
    expect_true(is.na(r$symbol), info = q)
  }
  # A well-formed Ensembl gene ID the OrgDb doesn't know says so
  r <- resolve_gene_query("ENSMUSG99999999999", "org.Mm.eg.db")
  expect_true(is.na(r$symbol))
  expect_match(r$error, "has no gene symbol")
})

test_that("gene_to_coords(): an Ensembl ID gives the same region as its symbol", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  by_sym <- gene_to_coords("Elf5", TXDBS$mouse, orgdb_pkg = "org.Mm.eg.db")
  id     <- resolve_gene_query(AnnotationDbi::select(org.Mm.eg.db::org.Mm.eg.db, "Elf5", "ENSEMBL", "SYMBOL")$ENSEMBL[1],
                               "org.Mm.eg.db")
  by_id  <- gene_to_coords(id$symbol, TXDBS$mouse, orgdb_pkg = "org.Mm.eg.db", ensembl_id = id$ensembl)
  expect_false(is.null(by_sym))
  expect_identical(by_id, by_sym)
  # A typed ID keeps its own locus when its symbol has several IDs (Amy2a1 has 3)
  amy  <- resolve_gene_query("ENSMUSG00000093931", "org.Mm.eg.db")
  want <- GenomicFeatures::genes(TXDBS$mouse, filter = list(gene_id = "ENSMUSG00000093931"))
  got  <- gene_to_coords(amy$symbol, TXDBS$mouse, orgdb_pkg = "org.Mm.eg.db", ensembl_id = amy$ensembl)
  expect_true(got$start < GenomicRanges::start(want) && got$end > GenomicRanges::end(want))
  # A symbol whose first OrgDb Ensembl ID is an alt-haplotype one still resolves (PTPRC, hg38)
  expect_false(is.null(gene_to_coords("PTPRC", TXDBS$human, orgdb_pkg = "org.Hs.eg.db")))
})

test_that("a gene the OrgDb has no symbol for opens by its GTF gene_name or gene_id (bug 258)", {
  skip_if_not_installed("org.Mm.eg.db")
  skip_if_not_installed("org.Hs.eg.db")
  mm <- suppressWarnings(pt_project(file.path(tempdir(), "pt-test-gtf-names-mm"), pt_example("mammary")))
  hs <- suppressWarnings(pt_project(file.path(tempdir(), "pt-test-gtf-names-hs"), pt_example("human")))
  opens <- function(proj, query, symbol, id) {
    fig <- pt_figure(proj, query, c("genes", "transcripts"))
    expect_identical(fig$gene, symbol, label = query)
    # The view holds the GTF's gene, and its transcripts track finds the gene's transcripts.
    g <- GenomicFeatures::genes(proj$txdb, filter = list(gene_id = id))
    expect_true(fig$region_params$start <= GenomicRanges::start(g) &&
                  fig$region_params$end >= GenomicRanges::end(g), label = query)
    tx <- gtf_transcript_structure(fig$gene, .pt_annotation(proj$render_cfg))
    expect_gt(nrow(tx), 0)
  }
  # Mouse: GTF names newer than org.Mm.eg.db, by name (any case) or ID (any version)
  opens(mm, "Gm37220", "Gm37220", "ENSMUSG00000102686")
  opens(mm, "gm37220", "Gm37220", "ENSMUSG00000102686")
  opens(mm, "ENSMUSG00000102686.1", "Gm37220", "ENSMUSG00000102686")
  opens(mm, "BC016548", "BC016548", "ENSMUSG00000070794")
  opens(mm, "ENSMUSG00000070794", "BC016548", "ENSMUSG00000070794")
  # Human: unnamed GENCODE genes (gene_name is the ID)
  for (id in c("ENSG00000259371", "ENSG00000259321", "ENSG00000259529", "ENSG00000289757"))
    opens(hs, id, id, id)
  # A gene the OrgDb has still resolves through it
  expect_identical(pt_figure(mm, "elf5", "genes")$gene, "Elf5")
})

test_that("a gene missing only from the project's GTF says so, naming the GTF (bug 258)", {
  skip_if_not_installed("org.Mm.eg.db")
  skip_if_not_installed("org.Hs.eg.db")
  mm <- suppressWarnings(pt_project(file.path(tempdir(), "pt-test-gtf-names-mm"), pt_example("mammary")))
  hs <- suppressWarnings(pt_project(file.path(tempdir(), "pt-test-gtf-names-hs"), pt_example("human")))
  expect_error(pt_figure(hs, "TP53", "genes"), "^'TP53' is not in this project's GTF \\(annotation\\.gtf\\.gz\\)$")
  expect_error(pt_figure(mm, "Trp53", "genes"), "'Trp53' is not in this project's GTF \\(gencode\\.vM25\\.example\\.gtf\\.gz\\)")
  expect_error(pt_figure(hs, "ENSG00000141510", "genes"), "'ENSG00000141510' is not in this project's GTF")
  expect_error(pt_figure(mm, "ENSMUSG99999999999", "genes"),
               "Ensembl ID 'ENSMUSG99999999999' is not in this project's GTF \\(gencode")
  # A name neither knows: as before
  expect_error(pt_figure(mm, "Foo1", "genes"), "'Foo1' not found in mm10: enter a gene symbol")
})

test_that("a GTF gene name shared by several genes, or an OrgDb symbol of another gene, is not guessed", {
  names <- structure(c(ENSMUSG99999900001 = "Zzdup1x", ENSMUSG99999900002 = "Zzdup1x", ENSMUSG99999900003 = "Pax5"),
                     gtf = "x.gtf")
  r <- resolve_gene_query("Zzdup1x", "org.Mm.eg.db", gene_names = names)
  expect_true(is.na(r$symbol))
  expect_match(r$error, "Several genes in this project's GTF \\(x.gtf\\) are named 'Zzdup1x'")
  # By ID: the gene's own ID is its symbol, so the tracks find this gene and no other
  expect_identical(resolve_gene_query("ENSMUSG99999900002", "org.Mm.eg.db", gene_names = names)$symbol,
                   "ENSMUSG99999900002")
  expect_identical(resolve_gene_query("ENSMUSG99999900003", "org.Mm.eg.db", gene_names = names)$symbol,
                   "ENSMUSG99999900003")
})

test_that("an OrgDb symbol whose Ensembl IDs the GTF has under another ID opens by the GTF's (bug 258 review)", {
  skip_if_not_installed("org.Mm.eg.db")
  # The mouse example with Elf5 under an ID org.Mm.eg.db does not have.
  ex <- dirname(pt_example("mammary"))
  dir <- tempfile("pt-test-id-drift-")
  dir.create(dir)
  gtf <- file.path(dir, "drift.gtf")
  writeLines(gsub("ENSMUSG00000027186", "ENSMUSG99999927186",
                  readLines(gzfile(file.path(ex, "gencode.vM25.example.gtf.gz"))), fixed = TRUE), gtf)
  st <- file.path(dir, "samples.txt")
  txt <- readLines(pt_example("mammary"))
  txt <- sub("^# data_root: .*", paste("# data_root:", ex), txt)
  writeLines(sub("^# gtf_path: .*", paste("# gtf_path:", gtf), txt[!grepl("^# start_", txt)]), st)
  proj <- suppressWarnings(pt_project(file.path(dir, "proj"), st))
  fig <- pt_figure(proj, "Elf5", c("genes", "transcripts"))
  expect_identical(fig$gene, "Elf5")
  g <- GenomicFeatures::genes(proj$txdb, filter = list(gene_id = "ENSMUSG99999927186"))
  expect_true(fig$region_params$start <= GenomicRanges::start(g) && fig$region_params$end >= GenomicRanges::end(g))
  expect_gt(nrow(gtf_transcript_structure("Elf5", .pt_annotation(proj$render_cfg))), 0)
})
