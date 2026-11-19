# The transcript name, transcript meta and gene name lookups (R/genome-build-tx-name-lookup.R,
# R/genome-build-tx-meta-lookup.R) are built from one read of the GTF, a chunk of lines at a
# time, keeping only the transcript and gene lines (bug 272). They must come out exactly as
# from the whole file read at once (the code before bug 272), whatever the chunk size.

lookup_builders <- list(tx_names = build_tx_name_lookup, tx_meta = build_tx_meta_lookup,
                        gene_names = build_gene_name_lookup)

# Each lookup built from `lines` (NULL: the builder's own read): the object and the saved file's md5.
build_lookups <- function(gtf, lines = NULL) {
  lapply(lookup_builders, function(build) {
    out <- tempfile(fileext = ".rds")
    on.exit(unlink(out))
    obj <- if (is.null(lines)) build(gtf, out) else build(gtf, out, lines)
    list(obj = obj, file = readRDS(out), md5 = unname(tools::md5sum(out)))
  })
}

test_that("the lookups from a chunked read are identical to those from the whole file (bug 272)", {
  ex <- system.file("extdata", c("example", "example-human"), package = "PlotTrellis")
  # the mini fixture is in the source tree only, not in a built package
  mini <- test_path("fixtures", "mini", "annotations", "gencode.vM25.annotation.gtf.gz")
  gtfs <- c(list.files(ex, "\\.gtf\\.gz$", full.names = TRUE), if (file.exists(mini)) mini)
  expect_gte(length(gtfs), 2)
  for (gtf in gtfs) {
    whole <- build_lookups(gtf, readLines(gzfile(gtf)))
    expect_gt(length(whole$tx_names$obj), 0)
    expect_gt(length(whole$gene_names$obj), 0)
    # chunk = 1 and 7 split every gene's lines across chunks; 2e5 is the default (one chunk here)
    for (chunk in c(1L, 7L, 200000L)) {
      lines <- .pt_gtf_feature_lines(gtf, c("transcript", "gene"), chunk = chunk)
      expect_identical(build_lookups(gtf, lines), whole, info = paste(basename(gtf), chunk))
    }
    # each builder reading the file itself
    expect_identical(build_lookups(gtf), whole, info = basename(gtf))
  }
})

test_that("the chunked read keeps exactly the whole file's transcript and gene lines (bug 272)", {
  row <- function(feature, attrs) paste("chr1", "src", feature, 1, 100, ".", "+", ".", attrs, sep = "\t")
  long <- strrep("x", 200000)  # longer than any read buffer
  rows <- c("##comment", row("gene", 'gene_id "G1"; gene_name "A";'),
            row("transcript", paste0('gene_id "G1"; transcript_id "T1"; note "', long, '";')),
            row("exon", 'gene_id "G1"; transcript_id "T1";'),
            row("transcript", 'gene_id "G1"; transcript_id "T2";'),
            row("gene", 'gene_id "G2";'))
  want <- function(x) x[grepl("\ttranscript\t", x, fixed = TRUE) | grepl("\tgene\t", x, fixed = TRUE)]
  for (eol in c("\n", "\r\n")) for (last_eol in c(TRUE, FALSE)) for (gz in c(TRUE, FALSE)) {
    f <- tempfile(fileext = if (gz) ".gtf.gz" else ".gtf")
    con <- if (gz) gzfile(f, "wb") else file(f, "wb")
    writeChar(paste0(paste(rows, collapse = eol), if (last_eol) eol), con, eos = NULL)
    close(con)
    whole <- suppressWarnings(readLines(gzfile(f)))  # warns on a missing last newline
    expect_identical(whole, rows)
    for (chunk in 1:4) {
      got <- suppressWarnings(.pt_gtf_feature_lines(f, c("transcript", "gene"), chunk = chunk))
      expect_identical(got, want(whole), info = paste(eol == "\n", last_eol, gz, chunk))
    }
    expect_identical(suppressWarnings(.pt_gtf_feature_lines(f, "nothing")), character(0))
    unlink(f)
  }
})
