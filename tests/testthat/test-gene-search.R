# Phase I.1 step 4a (BUGS.md #36, #49): gene search accepts a symbol in any
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

test_that("an Ensembl ID with several symbols takes the first and names the rest", {
  r <- resolve_gene_query("ENSMUSG00000093931", "org.Mm.eg.db")
  expect_equal(r$symbol, "Amy2a1")
  expect_equal(r$also, "Amy2a3")
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
