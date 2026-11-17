# Rebuild tests/testthat/fixtures/mini/: the small, committed copy of the
# mouse test fixtures (Elf5 only) that tests/testthat/helper-fixtures.R falls
# back to when the full data/ folder is absent, so the mouse tests run on any
# machine. Same layout and file names as data/, cut down to WINDOW, so a figure
# drawn from it is identical to one drawn from the full files.
#
# Not included: Hi-C (.hic files can't be cropped with what we have), the
# human dataset, the other Salmon output files (only quant.sf is read).
#
# Needs the monorepo's data/ folder (PT_TEST_DATA_DIR overrides, as in
# helper-fixtures.R). Run from the repo root so renv activates:
#   module load R/4.5.3
#   Rscript Shiny-apps/PlotTrellis-pkg/data-raw/build-test-fixtures.R

suppressMessages({
  library(GenomicRanges)
  library(rtracklayer)
})

script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
pkg <- normalizePath(file.path(dirname(script), ".."))
data_dir <- Sys.getenv("PT_TEST_DATA_DIR", file.path(pkg, "..", "..", "data"))
out <- file.path(pkg, "tests", "testthat", "fixtures", "mini")

# chr2:102.9-104.0 Mb (mm10): Elf5 (103.41-103.45 Mb) with a wide margin, so
# zoomed-out test views stay inside. Genes overlapping it are kept whole.
WINDOW <- GRanges("chr2", IRanges(102900000, 104000000))
ATAC <- c("Ba", "LP", "P6", "P7")
LOOPS <- paste0(ATAC, "_OmniC_FitHiChIP_10kb_FDR01.bedpe")
QUANT <- c("Basal_rep1", "LP_rep1")

unlink(out, recursive = TRUE)
mouse_out <- file.path(out, "mouse_mammary_test")
dir.create(file.path(out, "annotations"), recursive = TRUE)
dir.create(file.path(mouse_out, "atac"), recursive = TRUE)
dir.create(file.path(mouse_out, "loops"), recursive = TRUE)

# ---- GTF: every line of each gene overlapping WINDOW, plus the header -------
lines <- readLines(file.path(data_dir, "annotations", "gencode.vM25.annotation.gtf.gz"))
header <- lines[startsWith(lines, "#")]
body <- lines[startsWith(lines, "chr2\t")]
cols <- do.call(rbind, strsplit(body, "\t", fixed = TRUE))[, c(3, 4, 5)]
gene_id <- sub('.*gene_id "([^"]+)".*', "\\1", body)
is_gene <- cols[, 1] == "gene"
genes <- GRanges("chr2", IRanges(as.integer(cols[is_gene, 2]), as.integer(cols[is_gene, 3])),
                 gene_id = gene_id[is_gene])
genes <- genes[overlapsAny(genes, WINDOW)]
kept <- body[gene_id %in% genes$gene_id]
con <- gzfile(file.path(out, "annotations", "gencode.vM25.annotation.gtf.gz"), "w")
writeLines(c(header, kept), con)
close(con)
tx_ids <- unique(sub('.*transcript_id "([^"]+)".*', "\\1", grep("\ttranscript\t", kept, value = TRUE)))
cat(sprintf("GTF: %d genes, %d transcripts, %d lines\n", length(genes), length(tx_ids), length(kept)))

# ---- Chromosome sizes: the full set's cache, so nothing is downloaded -------
file.copy(file.path(data_dir, "annotations", "gencode.vM25.annotation_chrom_sizes.rds"),
          file.path(out, "annotations"))

# ---- ATAC: the WINDOW bins, checked against the source over Elf5 ------------
elf5 <- GRanges("chr2", IRanges(103411688, 103450989))
for (cond in ATAC) {
  src <- file.path(data_dir, "mouse_mammary_test", "atac", paste0(cond, "_ATAC.bw"))
  dst <- file.path(mouse_out, "atac", paste0(cond, "_ATAC.bw"))
  export.bw(import(BigWigFile(src), which = WINDOW), dst)
  a <- import(src, which = elf5)
  b <- import(dst, which = elf5)
  stopifnot(identical(ranges(a), ranges(b)), isTRUE(all.equal(a$score, b$score)))
}

# ---- Loops: whole lines with either anchor in WINDOW ------------------------
# LP's source file is a 2-byte empty file; an empty one is copied as it is.
for (f in LOOPS) {
  src <- readLines(file.path(data_dir, "mouse_mammary_test", "loops", f), warn = FALSE)
  f6 <- strsplit(src, "\t", fixed = TRUE)
  in_window <- vapply(f6, function(x) {
    if (length(x) < 6) return(FALSE)
    anchor <- function(i) x[i] == "chr2" && as.numeric(x[i + 2]) > start(WINDOW) && as.numeric(x[i + 1]) < end(WINDOW)
    anchor(1) || anchor(4)
  }, logical(1))
  if (!any(nzchar(src))) {
    file.copy(file.path(data_dir, "mouse_mammary_test", "loops", f), file.path(mouse_out, "loops"))
  } else {
    writeLines(src[in_window], file.path(mouse_out, "loops", f))
  }
  cat(sprintf("Loops %s: %d of %d lines\n", f, sum(in_window), length(src)))
}

# ---- quant.sf: the rows of the GTF subset's transcripts, unchanged ----------
for (q in QUANT) {
  q_lines <- readLines(file.path(data_dir, "mouse_mammary_test", "quant", q, "quant.sf"))
  q_kept <- q_lines[c(1, which(sub("\t.*", "", q_lines) %in% tx_ids))]
  dir.create(file.path(mouse_out, "quant", q), recursive = TRUE)
  writeLines(q_kept, file.path(mouse_out, "quant", q, "quant.sf"))
}

file.copy(file.path(pkg, "data-raw", "test-fixtures-README"), file.path(out, "README"))

sizes <- file.size(list.files(out, recursive = TRUE, full.names = TRUE))
cat(sprintf("Total: %.0f KB in %d files\n", sum(sizes) / 1024, length(sizes)))
