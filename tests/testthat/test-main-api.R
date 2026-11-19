# The main screen's data answers (R/app-screen-api.R), on the bundled mammary example
# (Basal and LP; every data type; the regions of five genes on chr1, chr2, chr4, chr13 and
# chr14, mm39). Search ranking and the nearby cap run on a made-up gene index, since the
# example has only 47 genes. UniProt is mocked, so nothing here needs the network.

# The example's project and context, built once for this file.
main_ctx <- local({
  memo <- NULL
  function() {
    if (is.null(memo)) {
      proj <- tempfile("pt-main-")
      suppressWarnings(suppressMessages(pt_project(proj, pt_example("mammary"))))
      json <- file.path(normalizePath(proj, winslash = "/"), "project_settings.json")
      cfg <- .pt_read_settings(json, .pt_default_settings())
      memo <<- suppressMessages(.pt_composer_context(json, cfg))
    }
    memo
  }
})

# The JSON the page gets, read back as plain lists. A NA must never reach it: toJSON()
# writes a numeric NA as the string "NA".
as_json <- function(x) {
  txt <- as.character(jsonlite::toJSON(x, auto_unbox = TRUE, null = "null", digits = NA))
  expect_false(grepl('"NA"', txt, fixed = TRUE))
  jsonlite::fromJSON(txt, simplifyVector = FALSE)
}

# A made-up gene index: n genes on chr1, the given symbols (NA: none) and coding flags.
fake_index <- function(sym, coding, s = seq_along(sym) * 1000, e = s + 500,
                       ens = sprintf("ENSMUSG%011d", seq_along(sym))) {
  list(genes = data.frame(ens = ens, chr = "chr1", s = s, e = e, strand = "+", sym = sym,
                          display = ifelse(is.na(sym), paste0("Gm", seq_along(sym)), sym),
                          coding = coding, stringsAsFactors = FALSE),
       tx = NULL)
}

test_that("the gene index holds every example gene, named, and is built once", {
  skip_if_not_installed("org.Mm.eg.db")
  ctx <- main_ctx()
  idx <- .pt_main_gene_index(ctx)
  expect_named(idx$genes, c("ens", "chr", "s", "e", "strand", "sym", "display", "coding"))
  g <- idx$genes
  gtf <- readLines(gzfile(file.path(dirname(pt_example("mammary")), "annotation.gtf.gz")))
  gtf_genes <- sub('.*gene_name "([^"]+)".*', "\\1", grep("\tgene\t", gtf, value = TRUE))
  expect_length(gtf_genes, 47)
  expect_setequal(g$display, gtf_genes)
  # No OrgDb symbol: named from the GENCODE transcript names instead.
  expect_true(is.na(g$sym[g$display == "Gm28199"]))
  expect_identical(g$coding[match(c("Eya2", "Gm28199"), g$display)], c(TRUE, FALSE))
  expect_identical(unlist(g[g$display == "Eya2", c("chr", "s", "e", "strand")], use.names = FALSE),
                   c("chr2", "165436895", "165613647", "+"))
  expect_identical(.pt_main_gene_index(ctx), idx)
  expect_true(is.data.frame(idx$tx) && all(c("tx_name", "nexon", "tx_len", "cds_len") %in% names(idx$tx)))
})

test_that("the project gives genome, suggested genes, samples, types in order and lengths", {
  skip_if_not_installed("org.Mm.eg.db")
  ctx <- main_ctx()
  p <- .pt_main_project(ctx)
  # Nothing preset: the page opens empty.
  expect_named(p, c("genome", "suggested", "samples", "types", "chrLen"))
  expect_identical(p$genome, "Mouse \u00B7 mm39")  # the GTF is "annotation.gtf.gz": a name that says nothing is left out
  # The example's own suggested genes: its samples.txt's line, every one in its annotation.
  expect_identical(as.character(p$suggested), c("Eya2", "Lrrfip1", "Acot7", "Dzip1", "Cxcl14"))
  expect_identical(p$samples, list(list(name = "Basal", colour = "#b31e1e"),
                                   list(name = "LP", colour = "#4b7538")))
  ids <- vapply(p$types, `[[`, "", "id")
  expect_identical(ids, c("RNA", "ATAC", "H3K36me3", "IgG", "HiC", "ATAC_peaks", "H3K36me3_peaks", "junctions",
                          "transcript_expression", "genes", "protein_domains", "transcripts_all", "transcripts",
                          "transcripts_annotated", "RNA_diff", "ATAC_diff", "H3K36me3_diff", "IgG_diff",
                          "junctions_diff", "junctions_diff_tpm"))
  by_id <- setNames(p$types, ids)
  expect_identical(by_id$ATAC[c("name", "cat", "kind", "gene", "viz")],
                   list(name = "ATAC", cat = "ATAC", kind = "sample", gene = FALSE, viz = "bw"))
  expect_identical(as.character(by_id$ATAC$samples), c("Basal", "LP"))
  expect_identical(c(by_id$HiC$viz, by_id$ATAC_peaks$viz, by_id$junctions$viz), c("hic", "bed", "junctions"))
  expect_identical(by_id$transcript_expression[c("cat", "gene", "viz")],
                   list(cat = "expr", gene = TRUE, viz = "tpm_usage"))
  expect_identical(by_id$genes[c("name", "cat", "kind", "gene")],
                   list(name = "Gene annotation", cat = "ann", kind = "fixed", gene = FALSE))
  expect_null(by_id$genes$samples)
  expect_identical(by_id$ATAC_diff[c("name", "cat", "kind", "viz")],
                   list(name = "ATAC diff (stacked)", cat = "cmp", kind = "comparison", viz = "bw_diff_stacked"))
  expect_setequal(as.character(by_id$ATAC_diff$samples), c("Basal", "LP"))
  expect_identical(unlist(p$chrLen), setNames(as.numeric(ctx$CHR_LENGTHS), names(ctx$CHR_LENGTHS)))
  j <- as_json(p)
  expect_named(j$samples[[1]], c("name", "colour"))
  for (t in j$types) expect_true(all(c("id", "name", "cat", "kind", "gene") %in% names(t)))
  expect_type(j$types[[1]]$samples, "list")  # an array, even of one
  expect_identical(j$chrLen$chr2, 181755017L)
})

test_that("suggested genes: the project's list by symbol, Ensembl id or GTF name, in its order", {
  skip_if_not_installed("org.Mm.eg.db")
  ctx <- main_ctx()
  g <- .pt_main_gene_index(ctx)$genes
  named <- g[!is.na(g$sym), ]
  gtf_only <- g[is.na(g$sym), ]
  expect_gte(nrow(named), 3)
  expect_gte(nrow(gtf_only), 1)
  set <- function(x) {
    ctx$CFG$suggested_genes <- x
    ctx
  }
  # none listed: none suggested
  expect_identical(.pt_main_suggested(set("")), character(0))
  expect_identical(.pt_main_suggested(set(NULL)), character(0))
  expect_identical(.pt_main_search(set(""), ""), list())
  # a symbol in another case, an Ensembl id with a version, a gene the GTF alone names; blanks,
  # a repeat and a name the annotation does not have are dropped
  some <- set(paste(c(toupper(named$sym[2]), paste0(named$ens[1], ".7"), "", gtf_only$display[1],
                      "NoSuchGene", named$sym[2], tolower(named$sym[3])), collapse = " , "))
  want <- c(named$sym[2], named$sym[1], gtf_only$display[1], named$sym[3])
  expect_identical(.pt_main_suggested(some), want)
  expect_identical(as.character(.pt_main_project(some)$suggested), want)
  expect_identical(.pt_main_suggested_rows(some)$unknown, "NoSuchGene")
  # the empty search box lists them, each with its place, and each opens as listed
  hits <- .pt_main_search(some, "")
  expect_identical(vapply(hits, `[[`, "", "sym"), want)
  expect_identical(hits[[1]], list(sym = named$sym[2], chr = named$chr[2], s = named$s[2],
                                   e = named$e[2], ens = named$ens[2]))
  for (h in hits) expect_identical(.pt_main_locus(some, h$sym)[c("gene", "chr", "s", "e")],
                                   list(gene = h$sym, chr = h$chr, s = h$s, e = h$e), info = h$sym)
  # a name is listed once: two genes of one name (the first written is the one listed), and a
  # name the annotation does not have written twice
  twins <- named[1:3, ]
  twins$sym[2] <- twins$display[2] <- twins$sym[1]
  r <- with_mocked_bindings(
    .pt_main_suggested_rows(set(paste(twins$ens[2], twins$sym[1], twins$sym[3], "NoSuchGene", "nosuchgene",
                                      sep = ", "))),
    .pt_main_gene_index = function(ctx) list(genes = twins))
  expect_identical(r$rows$sym, c(twins$sym[1], twins$sym[3]))
  expect_identical(r$rows$ens, c(twins$ens[2], twins$ens[3]))
  expect_identical(r$unknown, "NoSuchGene")
  # typing still searches every gene, not the list
  expect_identical(vapply(.pt_main_search(some, named$sym[1]), `[[`, "", "sym")[1], named$sym[1])
  # a list written as an array in project_settings.json reads the same
  expect_identical(.pt_main_suggested(set(list(named$sym[2], named$sym[1]))), named$sym[2:1])
})

test_that("at most 8 suggested genes are offered, the first 8 listed", {
  ctx <- list(CFG = list(suggested_genes = paste(paste0("g", 12:1), collapse = ",")))
  local_mocked_bindings(.pt_main_gene_index = function(ctx) fake_index(paste0("G", 1:12), rep(TRUE, 12)))
  expect_identical(.pt_main_suggested(ctx), paste0("G", 12:5))
  expect_length(.pt_main_search(ctx, ""), 8)
})

test_that("search: case-insensitive prefix, exact match first, then coding, then shorter, at most 8", {
  skip_if_not_installed("org.Mm.eg.db")
  ctx <- main_ctx()
  local_mocked_bindings(.pt_main_gene_index = function(ctx) fake_index(
    c("Abcf", "Abcde1", "Abcd", "Abc", "Abx", "Ab1", "Ab2", "Ab3", "Ab4", "Ab5", "Ab6", NA),
    c(FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)))
  syms <- function(q) vapply(.pt_main_search(ctx, q), `[[`, "", "sym")
  expect_identical(syms("abc"), c("Abc", "Abcd", "Abcde1", "Abcf"))
  expect_identical(syms("ABCD"), c("Abcd", "Abcde1"))
  expect_length(syms("a"), 8)
  expect_identical(syms("zz"), character(0))
  # An Ensembl id prefix; a gene with no symbol is never offered.
  expect_identical(syms("ensmusg0000000000"), c("Ab1", "Ab2", "Ab3", "Ab4", "Abx", "Abcd", "Abcde1", "Abc"))
  expect_identical(syms("ENSMUSG00000000012"), character(0))
  expect_identical(.pt_main_search(ctx, "Abx")[[1]],
                   list(sym = "Abx", chr = "chr1", s = 5000, e = 5500, ens = "ENSMUSG00000000005"))
})

test_that("search on the example: prefix, Ensembl id, and the suggested genes when empty", {
  skip_if_not_installed("org.Mm.eg.db")
  ctx <- main_ctx()
  syms <- function(q) vapply(.pt_main_search(ctx, q), `[[`, "", "sym")
  expect_identical(syms("e"), c("Espn", "Eya2"))
  expect_identical(syms(" eya2 "), "Eya2")
  expect_identical(syms("ENSMUSG00000017897"), "Eya2")
  expect_identical(syms(""), .pt_main_suggested(ctx))
  expect_identical(syms(NULL), .pt_main_suggested(ctx))
  j <- as_json(.pt_main_search(ctx, "eya"))
  expect_identical(j, list(list(sym = "Eya2", chr = "chr2", s = 165436895L, e = 165613647L,
                                ens = "ENSMUSG00000017897")))
  expect_identical(as.character(jsonlite::toJSON(.pt_main_search(ctx, "zz"))), "[]")
})

test_that("locus: a gene gives the app's default window, a region with commas itself", {
  skip_if_not_installed("org.Mm.eg.db")
  ctx <- main_ctx()
  g <- .pt_main_locus(ctx, "eya2")
  expect_identical(g, list(kind = "gene", label = "Eya2", gene = "Eya2", chr = "chr2", s = 165436895L,
                           e = 165613647L, chrLen = 181755017, view = list(s = 165420000L, e = 165640000L)))
  expect_identical(.pt_main_locus(ctx, "ENSMUSG00000017897.5")$view, g$view)
  r <- .pt_main_locus(ctx, " chr2:165,430,000 - 165,620,000 ")
  expect_identical(r, list(kind = "region", label = "chr2:165,430,000-165,620,000", gene = NULL, chr = "chr2",
                           s = 165430000, e = 165620000, chrLen = 181755017,
                           view = list(s = 165430000, e = 165620000)))
  j <- as_json(r)
  expect_named(j, c("kind", "label", "gene", "chr", "s", "e", "chrLen", "view"))
  expect_null(j$gene)
  expect_identical(j$view, list(s = 165430000L, e = 165620000L))
})

test_that("locus: each error says what is wrong", {
  skip_if_not_installed("org.Mm.eg.db")
  ctx <- main_ctx()
  err <- function(q) expect_error(.pt_main_locus(ctx, q), class = "error")$message
  expect_identical(err(""), "Type a gene or a region such as chr2:165,430,000-165,620,000.")
  expect_identical(err(NULL), "Type a gene or a region such as chr2:165,430,000-165,620,000.")
  expect_identical(err("chrZ:1-10"), "No chromosome 'chrZ' in mm39.")
  expect_identical(err("chr2:500-100"), "The region's end (100) must be after its start (500).")
  expect_identical(err("chr2:0-100"), "The region must lie within chr2 (1 to 181,755,017).")
  expect_identical(err("chr2:1-999,999,999"), "The region must lie within chr2 (1 to 181,755,017).")
  expect_identical(err("Notagene"), "'Notagene' is not a gene in mm39.")
  # A real mouse gene the example's cut-down GTF does not hold.
  expect_identical(err("Cd4"), "Cd4 has no position in this project's annotation.")
  expect_identical(.pt_main_request(ctx, "locus", list(q = "chrZ:1-10")),
                   list(error = "No chromosome 'chrZ' in mm39."))
})

test_that("nearby: overlapping genes sorted by start, named from GENCODE when the OrgDb has no symbol", {
  skip_if_not_installed("org.Mm.eg.db")
  ctx <- main_ctx()
  nb <- .pt_main_nearby(ctx, "chr2", "165370000", 165440000)
  expect_identical(vapply(nb, `[[`, "", "sym"), c("Gm75470", "Gm74945", "Gm74577", "Gm75026", "Eya2"))
  expect_false(is.unsorted(vapply(nb, `[[`, 0, "s")))
  expect_identical(nb[[5]], list(sym = "Eya2", s = 165436895L, e = 165613647L, strand = "+"))
  expect_identical(.pt_main_nearby(ctx, "chr3", 1, 1e6), list())
  for (bad in list(list("chr2", NULL, 10), list(NULL, 1, 10), list("chr2", "x", 10)))
    expect_error(do.call(.pt_main_nearby, c(list(ctx), bad)), "nearby needs chr, s and e.", fixed = TRUE)
  j <- as_json(nb)
  expect_named(j[[1]], c("sym", "s", "e", "strand"))
})

test_that("nearby: over 80 genes keeps the named coding ones, then named, then the longest", {
  skip_if_not_installed("org.Mm.eg.db")
  ctx <- main_ctx()
  sym <- c(paste0("Cod", 1:30), paste0("Lnc", 1:30), rep(NA, 40))
  coding <- c(rep(TRUE, 30), rep(FALSE, 30), rep(TRUE, 40))
  s <- rev(seq_len(100)) * 1000  # start order opposite to the index order
  e <- s + c(rep(100, 60), 1:40)  # among the unnamed, the later ones are longer
  local_mocked_bindings(.pt_main_gene_index = function(ctx) fake_index(sym, coding, s, e))
  nb <- .pt_main_nearby(ctx, "chr1", 0, 2e5)
  expect_length(nb, 80)
  starts <- vapply(nb, `[[`, 0, "s")
  expect_false(is.unsorted(starts))
  kept <- vapply(nb, `[[`, "", "sym")
  expect_true(all(c(paste0("Cod", 1:30), paste0("Lnc", 1:30)) %in% kept))
  expect_setequal(setdiff(kept, sym), paste0("Gm", 81:100))
})

test_that("labels: per-sample, annotation, repeated Gene annotation and comparisons in the order given", {
  skip_if_not_installed("org.Mm.eg.db")
  ctx <- main_ctx()
  tracks <- list(list(type = "ATAC", sample = "LP"), list(type = "genes"), list(type = "ATAC_diff"),
                 list(type = "transcript_expression", sample = "Basal"), list(type = "genes"),
                 list(type = "ATAC", sample = "Basal"), list(type = "transcripts_annotated"))
  expect_identical(.pt_main_labels(ctx, c("LP", "Basal"), tracks),
                   c("bw_ATAC_LP", "genes", "bw_diff_stacked_ATAC_Basal_vs_LP", "tpm_usage_Basal", "genes",
                     "bw_ATAC_Basal", "transcripts_annotated"))
  # picked[1] is the comparison, picked[2] the reference.
  expect_identical(.pt_main_labels(ctx, list("Basal", "LP"), list(list(type = "ATAC_diff"))),
                   "bw_diff_stacked_ATAC_LP_vs_Basal")
  expect_identical(ctx$TRACK_ITEMS[["bw_diff_stacked_ATAC_LP_vs_Basal"]]$samples, c("LP", "Basal"))
  # A sample-free track ignores a sample; picked does not limit per-sample tracks.
  expect_identical(.pt_main_labels(ctx, character(0), list(list(type = "genes", sample = "LP"),
                                                            list(type = "ATAC", sample = "LP"))),
                   c("genes", "bw_ATAC_LP"))
  # The data frame jsonlite makes of a tracks array.
  df <- jsonlite::fromJSON('[{"type":"ATAC","sample":"LP"},{"type":"genes"}]')
  expect_identical(.pt_main_labels(ctx, "LP", df), c("bw_ATAC_LP", "genes"))
  expect_true(all(.pt_main_labels(ctx, c("LP", "Basal"), tracks) %in% names(ctx$TRACK_ITEMS)))
})

test_that("labels: an unknown type, a missing or unknown sample, or a comparison without two samples stop", {
  skip_if_not_installed("org.Mm.eg.db")
  ctx <- main_ctx()
  err <- function(picked, tracks) expect_error(.pt_main_labels(ctx, picked, tracks), class = "error")$message
  expect_identical(err("LP", list()), "Add at least one track.")
  expect_identical(err("LP", list(list(type = "nope"))), "Unknown track type 'nope'.")
  expect_identical(err("LP", list(list(sample = "LP"))), "Unknown track type ''.")
  expect_identical(err("LP", list(list(type = "ATAC"))), "'ATAC' is drawn per sample: say which sample.")
  expect_identical(err("LP", list(list(type = "ATAC", sample = ""))), "'ATAC' is drawn per sample: say which sample.")
  expect_identical(err("LP", list(list(type = "ATAC", sample = "ML"))), "No 'ATAC' track for sample 'ML'.")
  expect_identical(err("LP", list(list(type = "ATAC_diff"))), "'ATAC_diff' compares two samples: pick exactly two.")
  expect_identical(err(c("LP", "Basal", "ML"), list(list(type = "ATAC_diff"))),
                   "'ATAC_diff' compares two samples: pick exactly two.")
  expect_identical(err(c("LP", "LP"), list(list(type = "ATAC_diff"))), "No 'ATAC_diff' track for LP against LP.")
})

test_that("labels draw: entries built from them render an SVG figure", {
  skip_on_cran()
  skip_if_not_installed("org.Mm.eg.db")
  ctx <- main_ctx()
  labels <- .pt_main_labels(ctx, c("LP", "Basal"), list(
    list(type = "ATAC", sample = "LP"), list(type = "ATAC_diff"), list(type = "genes"), list(type = "genes")))
  cfg <- ctx$RENDER_CFG
  entries <- .pt_build_entries(labels, ctx$TRACK_ITEMS, ctx$DEFAULT_SAMPLE_COLS, cfg)
  expect_identical(vapply(entries, `[[`, "", "viz"), c("bw", "bw_diff_stacked", "genes", "genes"))
  loc <- .pt_main_locus(ctx, "Eya2")
  out <- tempfile(fileext = ".svg")
  on.exit(unlink(out), add = TRUE)
  settings <- c(.pt_settings_base(cfg), list(title = "", output_path = out))
  # a list: drawn (the error message otherwise)
  expect_type(.pt_render_svg(entries, list(chrom = loc$chr, start = loc$view$s, end = loc$view$e,
                                           gene = loc$gene), settings), "list")
  expect_true(file.exists(out))
  svg <- paste(readLines(out, warn = FALSE), collapse = "\n")
  expect_match(svg, "<svg", fixed = TRUE)
  expect_match(svg, "Eya2", fixed = TRUE)
})

test_that("geneinfo: transcripts with lengths and tags, the canonical protein's features, the UniProt link", {
  skip_if_not_installed("org.Mm.eg.db")
  ctx <- main_ctx()
  local_mocked_bindings(
    uniprot_canonical_match = function(gene_symbol, tx_ids, ...)
      data.frame(uniprot_match = tx_ids == "ENSMUST00000088132.14", uniprot_unique = TRUE),
    gene_protein_domains = function(...) list(
      status = "ok", acc = "O08575", canonical_tx = "ENSMUST00000088132.14", canonical_isoform_name = NA,
      domains = data.frame(domain_id = c("r1", "d1", "d1"), note = c("Disordered", "Eya domain", "Eya domain"),
                           aa_start = c(300, 260, 260), aa_end = c(320, 532, 532),
                           type = c("REGION", "DOMAIN", "DOMAIN"))))
  gi <- .pt_main_geneinfo(ctx, "eya2", modifyList(ctx$RENDER_CFG, list(uniprot_reviewed_only = TRUE)))
  expect_named(gi, c("gene", "chr", "s", "e", "strand", "transcripts", "protein", "uniprot_url"))
  expect_identical(gi[c("gene", "chr", "s", "e", "strand")],
                   list(gene = "Eya2", chr = "chr2", s = 165436895L, e = 165613647L, strand = "+"))
  tx <- setNames(gi$transcripts, vapply(gi$transcripts, `[[`, "", "name"))
  expect_identical(names(tx), sprintf("Eya2-%d", 201:215))
  expect_identical(tx[["Eya2-201"]][c("id", "biotype", "coding", "exons", "length", "aa")],
                   list(id = "ENSMUST00000063433.8", biotype = "protein_coding", coding = "coding",
                        exons = 16L, length = 2460L, aa = 532))
  expect_null(tx[["Eya2-203"]]$aa)
  expect_identical(tx[["Eya2-203"]]$coding, "no CDS")
  expect_identical(tx[["Eya2-209"]]$coding, "NMD")
  expect_identical(as.character(tx[["Eya2-202"]]$tags),
                   c("Ensembl canonical", "basic", "APPRIS principal_2", "CCDS", "UniProt canonical"))
  expect_identical(as.character(tx[["Eya2-201"]]$tags), c("basic", "APPRIS principal_2", "CCDS"))
  expect_identical(gi$protein, list(
    accession = "O08575", status = "ok", name = "Eya2 (Eya2-202)", length = 532,
    features = list(list(name = "Eya domain", start = 260, end = 532, source = "UniProt domain"),
                    list(name = "Disordered", start = 300, end = 320, source = "UniProt region"))))
  # The live settings' reviewed-only choice reaches the link.
  expect_match(gi$uniprot_url, "reviewed%3Atrue", fixed = TRUE)
  j <- as_json(gi)
  expect_named(j$transcripts[[1]], c("name", "id", "biotype", "coding", "exons", "length", "aa", "tags"))
  expect_type(j$transcripts[[2]]$tags, "list")
  expect_named(j$protein$features[[1]], c("name", "start", "end", "source"))
})

test_that("geneinfo: UniProt unreachable leaves the transcripts, with a protein status instead", {
  skip_if_not_installed("org.Mm.eg.db")
  ctx <- main_ctx()
  local_mocked_bindings(uniprot_canonical_match = function(...) stop("offline (test)"),
                        gene_protein_domains = function(...) list(status = "uniprot_unreachable", gene = "Eya2"))
  gi <- .pt_main_geneinfo(ctx, "Eya2", ctx$RENDER_CFG)
  expect_length(gi$transcripts, 15)
  expect_false(any(grepl("UniProt", unlist(lapply(gi$transcripts, `[[`, "tags")))))
  expect_null(gi[["protein"]])  # [[ ]]: $ would partly match protein_status
  expect_identical(gi$protein_status, .pd_status_text("uniprot_unreachable", panel = TRUE))
  expect_false(grepl("reviewed", gi$uniprot_url, fixed = TRUE))
  as_json(gi)
  expect_identical(.pt_main_request(ctx, "geneinfo", list(gene = "Cd4")),
                   list(error = "Cd4 has no position in this project's annotation."))
})

test_that("request dispatches each op and turns any error into list(error)", {
  skip_if_not_installed("org.Mm.eg.db")
  ctx <- main_ctx()
  expect_identical(.pt_main_request(ctx, "project"), .pt_main_project(ctx))
  expect_identical(.pt_main_request(ctx, "search", list(q = "e")), .pt_main_search(ctx, "e"))
  expect_identical(.pt_main_request(ctx, "locus", list(q = "Eya2")), .pt_main_locus(ctx, "Eya2"))
  expect_identical(.pt_main_request(ctx, "nearby", list(chr = "chr2", s = 165.3e6, e = 165.5e6)),
                   .pt_main_nearby(ctx, "chr2", 165.3e6, 165.5e6))
  expect_identical(.pt_main_request(ctx, "nearby", list(chr = "chr2")), list(error = "nearby needs chr, s and e."))
  expect_identical(.pt_main_request(ctx, "render"), list(error = "Unknown request 'render'."))
  expect_identical(.pt_main_request(ctx, NULL), list(error = "Unknown request ''."))
  expect_identical(as_json(.pt_main_request(ctx, "locus", list(q = ""))),
                   list(error = "Type a gene or a region such as chr2:165,430,000-165,620,000."))
})

test_that("a typed region has whole positions", {
  expect_identical(.pt_main_parse_region("chr2:1.5-100000.7"), list(chr = "chr2", s = 2, e = 100001))
  expect_identical(.pt_main_parse_region(" chr2 : 103,400,000 - 103,460,000 "),
                   list(chr = "chr2", s = 103400000, e = 103460000))
  expect_null(.pt_main_parse_region("Eya2"))
})

test_that("a gene the OrgDb has no symbol for opens by its GTF name", {
  skip_if_not_installed("org.Mm.eg.db")
  ctx <- main_ctx()
  g <- .pt_main_gene_index(ctx)$genes
  only_gtf <- g[is.na(g$sym) & g$display != g$ens & !duplicated(g$display) & !g$display %in% ctx$ALL_GENES, ]
  skip_if(nrow(only_gtf) == 0, "the example annotation has no gene named only in the GTF")
  one <- only_gtf[1, ]
  # What the shared gene lookup (pt_figure()'s too) makes of the same text.
  found <- .pt_resolve_gene_region(one$display, ctx$txdb, 0.15, ctx$ORGDB_PKG, ctx$ALL_GENES,
                                   ctx$GENOME$genome_label)
  skip_if(is.null(found$region), "the gene lookup does not find it either (no GTF gene names kept)")
  loc <- .pt_main_locus(ctx, one$display)
  expect_identical(loc$gene, found$symbol)
  expect_identical(c(loc$chr, loc$s, loc$e), c(one$chr, one$s, one$e))
})
