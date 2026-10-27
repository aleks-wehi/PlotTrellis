# Builds inst/extdata/example-human/, the small example dataset bundled with PlotTrellis
# (copied to a temporary folder and returned by plottrellis_example()).
#
# One ~65 kb locus on chr14 (hg38), cut from the ENCODE GM12878 vs K562 fixture in the
# monorepo's data/human_encode_test/ (see that folder's README.md for how those files
# were made). The PSME1 / EMC9 / PSME2 / RNF31 / IRF9 gene cluster was picked because
# both cell lines show ATAC signal and many splice junctions there, and PSME2 and
# PSME1 are expressed much higher in GM12878 than in K562 (a visible TPM contrast).
#
# Run from the package root (Shiny-apps/PlotTrellis-pkg/):
#   Rscript data-raw/make-example.R
# The full-size source data is looked for in $PT_TEST_DATA_DIR (the same variable the
# test suite uses), else ../../data. Not part of the installed package (.Rbuildignore).
#
# Nothing here is downloaded: the chromosome-length cache is copied from the source
# data's annotations folder (ensure_genome_ready() would otherwise fetch it from UCSC).

suppressMessages({
  library(GenomicRanges)
  library(rtracklayer)
})

SRC <- Sys.getenv("PT_TEST_DATA_DIR")
if (!nzchar(SRC)) SRC <- file.path("..", "..", "data")
stopifnot(dir.exists(SRC))
ENCODE <- file.path(SRC, "human_encode_test")
GTF_IN <- file.path(SRC, "annotations", "gencode.v46.annotation.gtf.gz")
SIZES_IN <- file.path(SRC, "annotations", "gencode.v46.annotation_chrom_sizes.rds")
stopifnot(dir.exists(ENCODE), file.exists(GTF_IN), file.exists(SIZES_IN))

OUT <- file.path("inst", "extdata", "example-human")
unlink(OUT, recursive = TRUE)
for (d in c("atac", "junctions", "quant/GM12878_rep1", "quant/K562_rep1"))
  dir.create(file.path(OUT, d), recursive = TRUE)

CHROM <- "chr14"
GENE_WINDOW <- c(24128000L, 24170000L)  # genes overlapping this are kept, with all their transcripts
PAD <- 10000L                           # flank of signal/junction data around the genes' window
CROP <- c(GENE_WINDOW[1] - PAD, GENE_WINDOW[2] + PAD)

conds <- c(GM12878 = "GM12878", K562 = "K562")

# ---- GTF subset: every line of every gene overlapping the window ----------------
lines <- readLines(gzfile(GTF_IN))
body <- lines[!startsWith(lines, "#")]
f <- strsplit(body, "\t", fixed = TRUE)
chrom <- vapply(f, `[`, "", 1); type <- vapply(f, `[`, "", 3)
st <- as.integer(vapply(f, `[`, "", 4)); en <- as.integer(vapply(f, `[`, "", 5))
gid <- sub('.*gene_id "([^"]+)".*', "\\1", vapply(f, `[`, "", 9))
keep_genes <- gid[type == "gene" & chrom == CHROM & st <= GENE_WINDOW[2] & en >= GENE_WINDOW[1]]
gtf_keep <- body[gid %in% keep_genes]
con <- gzfile(file.path(OUT, "annotation.gtf.gz"), "w")
writeLines(c(lines[startsWith(lines, "##")], gtf_keep), con)
close(con)
tx_keep <- unique(sub('.*transcript_id "([^"]+)".*', "\\1",
                      body[gid %in% keep_genes & type == "transcript"]))
cat(sprintf("GTF: %d genes, %d transcripts, %d lines\n", length(keep_genes), length(tx_keep), length(gtf_keep)))

# ---- Chromosome lengths cache (so ensure_genome_ready() never needs the network) --
# Named exactly as .genome_cache_path("annotation.gtf.gz", "_chrom_sizes.rds") derives.
file.copy(SIZES_IN, file.path(OUT, "annotation_chrom_sizes.rds"))
sizes <- readRDS(SIZES_IN)

# ---- Per-condition data -----------------------------------------------------------
which_gr <- GRanges(CHROM, IRanges(CROP[1], CROP[2]))
for (cond in names(conds)) {
  # ATAC bigWig, cropped to the locus. Only this chromosome goes in the file.
  x <- import(file.path(ENCODE, "atac", sprintf("%s_ATAC_foldchange.bigWig", cond)), which = which_gr)
  x <- x[, "score"]
  GenomeInfoDb::seqlevels(x) <- CHROM
  GenomeInfoDb::seqlengths(x) <- sizes[CHROM]
  x <- trim(x)
  export.bw(x, file.path(OUT, "atac", sprintf("%s_ATAC_foldchange.bigWig", cond)))

  # Junctions (BED12): header plus rows overlapping the window.
  jf <- file.path(ENCODE, "junctions", sprintf("%s_rep1.bed", cond))
  jl <- readLines(jf)
  jh <- jl[startsWith(jl, "#")]
  j <- strsplit(jl[!startsWith(jl, "#")], "\t", fixed = TRUE)
  jk <- vapply(j, function(r) r[1] == CHROM && as.integer(r[2]) < CROP[2] && as.integer(r[3]) > CROP[1], NA)
  writeLines(c(jh, vapply(j[jk], paste, "", collapse = "\t")),
             file.path(OUT, "junctions", sprintf("%s_rep1.bed", cond)))
  cat(sprintf("%s: %d bigWig ranges, %d junctions\n", cond, length(x), sum(jk)))

  # Salmon quant.sf: only the kept genes' transcripts, all columns, TPM untouched.
  q <- read.delim(file.path(ENCODE, "quant", sprintf("%s_rep1", cond), "quant.sf"),
                  check.names = FALSE, colClasses = "character")
  q <- q[q$Name %in% tx_keep, ]
  write.table(q, file.path(OUT, "quant", sprintf("%s_rep1", cond), "quant.sf"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("%s: %d quant rows, summed TPM %.1f\n", cond, nrow(q), sum(as.numeric(q$TPM))))
}

# ---- samples.txt (relative paths, resolved against this file's folder) -----------
writeLines(c(
  "# data_root: .",
  "# gtf_path: annotation.gtf.gz",
  "# orgdb: org.Hs.eg.db",
  "# genome_label: hg38",
  "condition\tcolour\tbw:atac\tpaired-end:junctions\tquantdir:quant",
  "GM12878\t#E41A1C\tatac/GM12878_ATAC_foldchange.bigWig\tjunctions/GM12878_rep1.bed\tquant/GM12878_rep1",
  "K562\t#377EB8\tatac/K562_ATAC_foldchange.bigWig\tjunctions/K562_rep1.bed\tquant/K562_rep1"
), file.path(OUT, "samples.txt"))

# ---- README -------------------------------------------------------------------------
writeLines(c(
  "PlotTrellis example dataset",
  "===========================",
  "",
  sprintf("A small slice of human (hg38) ENCODE data around %s:%d-%d (PSME1, EMC9, PSME2,",
          CHROM, GENE_WINDOW[1], GENE_WINDOW[2]),
  "RNF31, IRF9 and neighbours), for two cell lines: GM12878 (lymphoblastoid) and K562",
  "(myeloid leukaemia). Open it with plottrellis_example(); the tracks are called 'atac'",
  "and 'junctions' (plus transcript-usage tracks from the quant folders).",
  "",
  "Files",
  "  samples.txt           the PlotTrellis samples file (relative paths)",
  "  annotation.gtf.gz     GENCODE v46 lines for the genes overlapping the locus (complete transcripts)",
  "  annotation_chrom_sizes.rds  hg38 chromosome lengths (so no download is needed on first use)",
  "  atac/*.bigWig         ATAC-seq fold change over control, cropped to the locus +/- 10 kb",
  "  junctions/*.bed       splice junctions (BED12) overlapping the locus +/- 10 kb",
  "  quant/*/quant.sf      Salmon transcript quantification, rows for the kept genes' transcripts only",
  "                        (TPM values are as in the full run, not renormalised)",
  "",
  "Sources",
  "  ATAC-seq   ENCODE ENCFF316PVI (experiment ENCSR095QNB, GM12878) and ENCFF496GRF (ENCSR483RKN, K562),",
  "             fold change over control bigWigs",
  "  RNA-seq    ENCODE experiments ENCSR889TRN (GM12878, files ENCFF436TWZ + ENCFF286QHA) and ENCSR637VLS",
  "             (K562, file ENCFF000DWH), replicate 1 of each, quantified with Salmon 2.1.1 against GENCODE v46",
  "  Junctions  derived from the GRCh38 genome-alignment BAMs of the same replicates (ENCFF284YOU for GM12878,",
  "             ENCFF028IUE for K562) with GenomicAlignments::summarizeJunctions()",
  "  Annotation GENCODE release 46 (Ensembl 112), GRCh38",
  "",
  "Licence",
  "  ENCODE data are released without restriction on use (ENCODE data use policy, public data).",
  "  GENCODE annotation is open (EMBL-EBI terms of use). Please cite ENCODE and GENCODE if you reuse them.",
  "",
  "Rebuild: data-raw/make-example.R in the package source."
), file.path(OUT, "README"))

cat("Done. ", system(paste("du -sh", shQuote(file.path("inst", "extdata"))), intern = TRUE), "\n")
