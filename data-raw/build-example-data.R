# Rebuild the data files in inst/extdata/example/ (plottrellis_example()'s
# data, phase-I.2 step 23): the Ehf-Elf5-Cat locus of the mouse mammary
# dataset (Milevskiy et al. 2023, Cell Genomics, GEO GSE227750), Basal vs LP.
# samples.txt and README in that folder are hand-written, not built here.
#
# Needs network access (the ATAC bigWigs are cropped straight from GEO over
# HTTPS, and the chromosome sizes come from UCSC, as at a real first launch)
# and the monorepo's data/ folder (GENCODE vM25 GTF, our Salmon quant.sf
# output; PT_TEST_DATA_DIR overrides, as in tests/testthat/helper-fixtures.R).
# Run from the repo root so renv activates:
#   module load R/4.5.3
#   Rscript Shiny-apps/PlotTrellis-pkg/data-raw/build-example-data.R

suppressMessages({
  library(GenomicRanges)
  library(rtracklayer)
})

script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
pkg <- normalizePath(file.path(dirname(script), ".."))
data_dir <- Sys.getenv("PT_TEST_DATA_DIR", file.path(pkg, "..", "..", "data"))
out <- file.path(pkg, "inst", "extdata", "example")

# chr2:103.25-103.50 Mb (mm10): Ehf, Elf5, Cat and the lincRNAs/pseudogenes
# between them. No gene crosses either edge (checked below).
WINDOW <- GRanges("chr2", IRanges(103250000, 103500000))
GEO <- "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM7142nnn"
CONDITIONS <- list(
  Basal = list(atac = file.path(GEO, "GSM7142618/suppl/GSM7142618_Ba_S1_ATAC.bw"), quant = "Basal_rep1"),
  LP    = list(atac = file.path(GEO, "GSM7142621/suppl/GSM7142621_LP_S1_ATAC.bw"), quant = "LP_rep1")
)
gtf_out <- file.path(out, "gencode.vM25.example.gtf.gz")

dir.create(file.path(out, "atac"), recursive = TRUE, showWarnings = FALSE)

# ---- GTF: every line of each gene in WINDOW, plus GENCODE's header ----------
lines <- readLines(file.path(data_dir, "annotations", "gencode.vM25.annotation.gtf.gz"))
header <- lines[startsWith(lines, "#")]
body <- lines[startsWith(lines, "chr2\t")]
cols <- do.call(rbind, strsplit(body, "\t", fixed = TRUE))[, c(3, 4, 5)]
gene_id <- sub('.*gene_id "([^"]+)".*', "\\1", body)
is_gene <- cols[, 1] == "gene"
genes <- GRanges("chr2", IRanges(as.integer(cols[is_gene, 2]), as.integer(cols[is_gene, 3])),
                 gene_id = gene_id[is_gene])
genes <- genes[overlapsAny(genes, WINDOW)]
stopifnot(all(start(genes) >= start(WINDOW)), all(end(genes) <= end(WINDOW)))
kept <- body[gene_id %in% genes$gene_id]
con <- gzfile(gtf_out, "w")
writeLines(c(header, kept), con)
close(con)
tx_ids <- unique(sub('.*transcript_id "([^"]+)".*', "\\1", grep("\ttranscript\t", kept, value = TRUE)))
cat(sprintf("GTF: %d genes, %d transcripts, %d lines\n", length(genes), length(tx_ids), length(kept)))

# ---- Chromosome sizes: what ensure_genome_ready() would otherwise fetch -----
# from UCSC at first launch, under the name it looks for next to the GTF.
source(file.path(pkg, "R", "build_chrom_sizes.R"))
build_chrom_sizes("mm10", file.path(out, "gencode.vM25.example_chrom_sizes.rds"))

for (cond in names(CONDITIONS)) {
  src <- CONDITIONS[[cond]]

  # ---- ATAC: the WINDOW bins ------------------------------------------------
  # export.bw() only lists chromosomes that have data, so the file knows chr2
  # alone; a query elsewhere gets rtracklayer's "seqnames not known" warning.
  bwf <- BigWigFile(src$atac)
  gr <- import(bwf, which = WINDOW)
  bw_out <- file.path(out, "atac", paste0(cond, "_ATAC.bw"))
  export.bw(gr, bw_out)
  elf5 <- GRanges("chr2", IRanges(103411688, 103450989))
  a <- import(bwf, which = elf5)
  b <- import(bw_out, which = elf5)
  stopifnot(identical(ranges(a), ranges(b)), isTRUE(all.equal(a$score, b$score)))
  cat(sprintf("ATAC %s: %d bins, matches the GEO file over Elf5\n", cond, length(gr)))

  # ---- quant.sf: the rows of the GTF subset's transcripts, unchanged --------
  q_lines <- readLines(file.path(data_dir, "mouse_mammary_test", "quant", src$quant, "quant.sf"))
  q_name <- sub("\t.*", "", q_lines)
  q_kept <- q_lines[c(1, which(q_name %in% tx_ids))]
  stopifnot(setequal(sub("\t.*", "", q_kept[-1]), tx_ids))
  dir.create(file.path(out, "quant", src$quant), recursive = TRUE, showWarnings = FALSE)
  writeLines(q_kept, file.path(out, "quant", src$quant, "quant.sf"))
  cat(sprintf("quant.sf %s: %d transcripts\n", cond, length(q_kept) - 1))
}
