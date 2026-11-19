# Builds inst/extdata/example/ (what pt_example() and plottrellis_example() return), also
# the paper's example: five mouse mammary genes in basal cells (Basal) and luminal
# progenitors (LP), every data type PlotTrellis draws, cropped from a full processed run of
# the public GEO data of Milevskiy et al. 2023 (GSE227750). samples.txt and README in that
# folder are hand-written, not built here.
#
# Source: one PlotTrellis-wf run (Nextflow) over replicates 1 and 2 of each assay, every
# sample downloaded from SRA by its GSM (the run's samplesheet.csv lists them; the README
# gives them). Its results/ folder already holds one file per condition (the replicates
# merged by the workflow): RNA-seq coverage, splice junctions (BED12) and the mean of the
# two replicates' Salmon TPMs, ATAC-seq coverage and narrowPeaks, H3K36me3 and IgG CUT&Tag
# coverage and H3K36me3 broadPeaks, Omni-C contacts (.hic). Reference: GENCODE vM39
# primary assembly (mm39).
#
# What it does, per gene: a window of the gene plus a flank on each side (FLANK below).
# Every file is cut to the windows with its values unchanged:
#   - bigWigs: the bins overlapping a window (rtracklayer), every chromosome's true length;
#   - junctions and peaks: the text lines overlapping a window, as they are;
#   - quant.sf: the rows of the transcripts of the genes kept in the GTF, values as they
#     are (TPM is not renormalised, so the Transcript expression bars and shares are those
#     of the full file);
#   - GTF: every line of each gene overlapping a window, plus the header;
#   - Hi-C: the .hic's contacts with both ends in a wider window per gene (HIC_WINDOWS
#     below: what the Hi-C track reads at the gene's default view), at every bin size
#     from 5 kb up (MIN_RES below; the .hic's 1 kb is left out),
#     written as a .mcool (cooler's HDF5 layout, via rhdf5: no .hic writer is installed).
#     Raw counts go in pixels/count; each normalisation of the .hic (KR, SCALE, VC,
#     VC_SQRT) goes in as a bins column holding its per-bin vector, recovered from the
#     .hic as raw / normalised = n_i * n_j. plotgardener::readCool() divides by those
#     columns the way straw does, so each normalisation reads back as from the .hic, bin
#     for bin. (Not identical in every view: for a view not on bin edges, readCool() keeps
#     the part-covered edge bins that straw leaves out, so a view can differ at its edges
#     and so in its colour scale's top, as at Cxcl14's default view: 4.45 vs 4.46.)
# Every crop is checked against the source at each gene below.
#
# Run from the repo root on a compute node (the .hic reads need a few GB):
#   module load R/4.5.3
#   PT_PAPER_SRC=<run>/results Rscript Shiny-apps/PlotTrellis-pkg/data-raw/build-example-data.R
# Not part of the installed package (.Rbuildignore). Needs network once, for the mm39
# chromosome lengths from UCSC (as at a real first launch).

options(warn = 1)
suppressMessages({
  library(GenomicRanges)
  library(rtracklayer)
  library(rhdf5)
})

# The results/ folder of that run (the one its samples.txt is in).
SRC <- Sys.getenv("PT_PAPER_SRC")
stopifnot(dir.exists(file.path(SRC, "hic", "contacts")))
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
PKG <- normalizePath(file.path(dirname(script), ".."))
OUT <- file.path(PKG, "inst", "extdata", "example")

GENES <- c("Eya2", "Lrrfip1", "Acot7", "Dzip1", "Cxcl14")
# Each side of a gene: half its length, at least 25 kb, the window's ends then rounded
# out to 10 kb (whole Hi-C bins at every size up to 10 kb). The app's default view is
# the gene plus 15% each side, well inside.
FLANK_FRACTION <- 0.5
FLANK_MIN <- 25000L
ROUND <- 10000L

# Source file of each data type, by the source's condition prefix (Ba, LP), and where
# it goes in OUT, by our condition name.
COND <- c(Basal = "Ba", LP = "LP")
FILES <- list(
  rna      = list(src = "rna/coverage/%s_RNA.bw",                  out = "rna/%s_RNA.bw"),
  junc     = list(src = "rna/junctions/%s_RNA.junctions.bed",      out = "rna/%s_junctions.bed"),
  quant    = list(src = "rna/quantification_mean/%s_RNA/quant.sf", out = "rna/%s_salmon/quant.sf"),
  atac     = list(src = "atac/coverage/%s_ATAC.bw",                out = "atac/%s_ATAC.bw"),
  atac_pk  = list(src = "atac/peaks/%s_ATAC_peaks.narrowPeak",     out = "atac/%s_ATAC_peaks.narrowPeak"),
  k36      = list(src = "cuttag/coverage/%s_H3K36me3.bw",          out = "cuttag/%s_H3K36me3.bw"),
  k36_pk   = list(src = "cuttag/peaks/%s_H3K36me3_peaks.broadPeak", out = "cuttag/%s_H3K36me3_peaks.broadPeak"),
  igg      = list(src = "cuttag/coverage/%s_IgG.bw",               out = "cuttag/%s_IgG.bw"),
  hic      = list(src = "hic/contacts/%s_HiC.hic",                 out = "hic/%s_HiC.mcool")
)
GTF_IN <- file.path(SRC, "reference", "gencode.vM39.primary_assembly.annotation.gtf.gz")
GTF_OUT <- file.path(OUT, "annotation.gtf.gz")
src_path <- function(type, cond) file.path(SRC, sprintf(FILES[[type]]$src, COND[[cond]]))
out_path <- function(type, cond) file.path(OUT, sprintf(FILES[[type]]$out, cond))

# Keep samples.txt and README; everything else is rebuilt.
old <- setdiff(list.files(OUT, recursive = TRUE), c("samples.txt", "README"))
unlink(file.path(OUT, old))
for (d in unique(dirname(unlist(lapply(names(FILES), function(t) out_path(t, names(COND)))))))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- Windows and GTF ----------------------------------------------------------------
lines <- readLines(GTF_IN)
header <- lines[startsWith(lines, "#")]
body <- lines[!startsWith(lines, "#")]
f <- strsplit(body, "\t", fixed = TRUE)
chrom <- vapply(f, `[`, "", 1)
type <- vapply(f, `[`, "", 3)
st <- as.integer(vapply(f, `[`, "", 4))
en <- as.integer(vapply(f, `[`, "", 5))
attr9 <- vapply(f, `[`, "", 9)
gid <- sub('.*gene_id "([^"]+)".*', "\\1", attr9)
gname <- sub('.*gene_name "([^"]+)".*', "\\1", attr9)
is_gene <- type == "gene"
gene_gr <- GRanges(chrom[is_gene], IRanges(st[is_gene], en[is_gene]),
                   gene_id = gid[is_gene], gene_name = gname[is_gene])

targets <- gene_gr[match(GENES, gene_gr$gene_name)]
stopifnot(!anyNA(targets$gene_id))
flank <- pmax(as.integer(ceiling(width(targets) * FLANK_FRACTION)), FLANK_MIN)
WINDOWS <- GRanges(seqnames(targets),
                   IRanges(((start(targets) - flank) %/% ROUND) * ROUND + 1L,
                           ((end(targets) + flank) %/% ROUND + 1L) * ROUND),
                   gene = GENES)
# One window per chromosome (crop_lines() below measures distance to that one window).
stopifnot(!anyDuplicated(as.character(seqnames(WINDOWS))))
print(as.data.frame(WINDOWS)[, c("seqnames", "start", "end", "width", "gene")])

# Hi-C needs more: plotgardener's Rectangle reads the file over the view widened on each
# side by the view's width times the track's height/width ratio. Covered: the app's
# default gene view (the gene plus 15% each side) with a Hi-C track up to as tall as it
# is wide, so each gene's Hi-C window is that view plus its width on each side (and at
# least the window above).
view_flank <- as.integer(ceiling(width(targets) * 0.15))
view_w <- width(targets) + 2L * view_flank
HIC_WINDOWS <- GRanges(seqnames(targets),
                       IRanges(pmin(start(WINDOWS),
                                    ((start(targets) - view_flank - view_w) %/% ROUND) * ROUND + 1L),
                               pmax(end(WINDOWS),
                                    ((end(targets) + view_flank + view_w) %/% ROUND + 1L) * ROUND)),
                       gene = GENES)
stopifnot(length(reduce(HIC_WINDOWS)) == length(HIC_WINDOWS))
print(as.data.frame(HIC_WINDOWS)[, c("seqnames", "start", "end", "width", "gene")])

keep_genes <- gene_gr$gene_id[overlapsAny(gene_gr, WINDOWS)]
kept <- body[gid %in% keep_genes]
con <- gzfile(GTF_OUT, "w")
writeLines(c(header, kept), con)
close(con)
tx_ids <- unique(sub('.*transcript_id "([^"]+)".*', "\\1", kept[type[gid %in% keep_genes] == "transcript"]))
cat(sprintf("GTF: %d genes, %d transcripts, %d lines\n", length(keep_genes), length(tx_ids), length(kept)))

# ---- Chromosome sizes: what ensure_genome_ready() would otherwise fetch at first launch,
# under the name .pt_shipped_chrom_sizes() looks for next to the GTF.
source(file.path(PKG, "R", "genome-build-chrom-sizes.R"))
sizes <- build_chrom_sizes("mm39", file.path(OUT, "annotation_chrom_sizes.rds"))

# ---- Text files: the lines overlapping a window ---------------------------------------
crop_lines <- function(src, dest, zero_based = TRUE) {
  l <- readLines(src)
  hd <- startsWith(l, "#") | startsWith(l, "track") | startsWith(l, "browser")
  x <- strsplit(l[!hd], "\t", fixed = TRUE)
  gr <- GRanges(vapply(x, `[`, "", 1),
                IRanges(as.integer(vapply(x, `[`, "", 2)) + zero_based,
                        as.integer(vapply(x, `[`, "", 3))))
  keep <- overlapsAny(gr, WINDOWS, ignore.strand = TRUE)
  # A window's chromosome with no line in it keeps its nearest line (real data, outside
  # the window): the reader says "Chromosome not found" for a chromosome a file never
  # names, where the full file draws an empty panel.
  for (ch in unique(as.character(seqnames(WINDOWS)))) {
    on <- which(as.character(seqnames(gr)) == ch)
    if (!length(on) || any(keep[on])) next
    d <- distance(gr[on], WINDOWS[as.character(seqnames(WINDOWS)) == ch], ignore.strand = TRUE)
    keep[on[which.min(d)]] <- TRUE
  }
  writeLines(c(l[hd], l[!hd][keep]), dest)
  sum(keep)
}

# ---- bigWig: the bins overlapping a window, values unchanged ------------------------
crop_bw <- function(src, dest) {
  gr <- import(BigWigFile(src), which = WINDOWS)
  gr <- gr[, "score"]
  # Every window's chromosome is named in the file, with or without data in it (as
  # for the text files above).
  used <- unique(c(as.character(unique(seqnames(gr))), as.character(seqnames(WINDOWS))))
  GenomeInfoDb::seqlevels(gr) <- used
  GenomeInfoDb::seqlengths(gr) <- sizes[used]
  export.bw(gr, dest)
  for (i in seq_along(WINDOWS)) {
    a <- import(BigWigFile(src), which = WINDOWS[i])
    b <- import(BigWigFile(dest), which = WINDOWS[i])
    stopifnot(identical(start(a), start(b)), identical(end(a), end(b)),
              isTRUE(all.equal(a$score, b$score)))
  }
  length(gr)
}

# ---- Hi-C: .hic -> .mcool of the windows' contacts ----------------------------------
NORMS <- c("KR", "SCALE", "VC", "VC_SQRT")
# Bin sizes kept: every size of the .hic from 5 kb up. 1 kb is left out (owner, 2026-11-19):
# at gene views it is the size Auto picks and it draws sparse; without it Auto picks 5 kb.
MIN_RES <- 5000L
hic_resolutions <- function(hic) {
  r <- sort(strawr::readHicBpResolutions(hic))
  r[r >= MIN_RES]
}

# The contacts of one window at one bin size: raw counts plus each normalisation's
# per-bin vector (NA where the .hic has none), recovered as raw / normalised = n_i * n_j:
# from the diagonal where it has contacts, then along rows from bins already known.
hic_window <- function(hic, chrom, s, e, res) {
  s <- (s %/% res) * res
  e <- ((e - 1) %/% res) * res
  locus <- paste0(chrom, ":", format(s, scientific = FALSE), ":", format(e, scientific = FALSE))
  raw <- strawr::straw("NONE", hic, locus, locus, "BP", res)
  bins <- seq(s, e, by = res)
  vec <- list()
  for (nm in NORMS) {
    v <- tryCatch(strawr::straw(nm, hic, locus, locus, "BP", res), error = function(err) NULL)
    if (is.null(v)) next
    m <- merge(raw, v, by = c("x", "y"), suffixes = c("", ".n"))
    m <- m[is.finite(m$counts.n) & m$counts.n > 0 & m$counts > 0, ]
    m$r <- m$counts / m$counts.n
    # Bins by their index in `bins` (names made from numbers can differ: 1e+05).
    i <- (m$x - s) %/% res + 1
    j <- (m$y - s) %/% res + 1
    stopifnot(all(i >= 1 & i <= length(bins)), all(j >= 1 & j <= length(bins)))
    n <- rep(NA_real_, length(bins))
    on_diag <- i == j
    n[i[on_diag]] <- sqrt(m$r[on_diag])
    repeat {
      ki <- !is.na(n[i]); kj <- !is.na(n[j])
      todo <- xor(ki, kj)
      if (!any(todo)) break
      known <- ifelse(ki[todo], i[todo], j[todo])
      unknown <- ifelse(ki[todo], j[todo], i[todo])
      val <- m$r[todo] / n[known]
      first <- !duplicated(unknown)
      n[unknown[first]] <- val[first]
    }
    # Every normalised contact must read back from the vector.
    back <- m$counts / (n[i] * n[j])
    stopifnot(!anyNA(back), max(abs(back / m$counts.n - 1)) < 1e-5)
    vec[[nm]] <- n
  }
  list(raw = raw, bins = bins, vec = vec)
}

write_ds <- function(file, name, x, storage) {
  h5createDataset(file, name, dims = length(x), storage.mode = storage,
                  chunk = min(length(x), 100000L), level = 9, shuffle = TRUE,
                  size = if (storage == "character") max(nchar(x)) + 1L)
  h5write(x, file, name)
}

write_mcool <- function(hic, dest) {
  chroms_all <- strawr::readHicChroms(hic)
  avail <- hic_resolutions(hic)
  used <- unique(as.character(seqnames(HIC_WINDOWS)))
  # The windows' chromosomes, in the .hic's order.
  ch <- chroms_all[chroms_all$name %in% used, ]
  ch <- ch[order(ch$index), ]
  if (file.exists(dest)) unlink(dest)
  h5createFile(dest)
  h5createGroup(dest, "resolutions")
  fid <- H5Fopen(dest)
  h5writeAttribute("HDF5::MCOOL", fid, "format")
  h5writeAttribute(2L, fid, "format-version")
  H5Fclose(fid)
  for (res in avail) {
    g <- paste0("resolutions/", as.integer(res))  # not 1e+05
    h5createGroup(dest, g)
    for (sub in c("chroms", "bins", "pixels", "indexes")) h5createGroup(dest, file.path(g, sub))
    nb <- ceiling(ch$length / res)
    offs <- c(0, cumsum(nb))
    b_chrom <- rep(seq_along(ch$name) - 1L, nb)
    b_start <- unlist(lapply(seq_along(nb), function(i) (seq_len(nb[i]) - 1) * res))
    b_end <- pmin(b_start + res, ch$length[b_chrom + 1L])
    norm_cols <- list()
    px <- list()
    for (w in seq_along(HIC_WINDOWS)) {
      ci <- match(as.character(seqnames(HIC_WINDOWS[w])), ch$name)
      hw <- hic_window(hic, ch$name[ci], start(HIC_WINDOWS[w]) - 1L, end(HIC_WINDOWS[w]), res)
      id <- function(pos) offs[ci] + pos %/% res
      for (nm in names(hw$vec)) {
        if (is.null(norm_cols[[nm]])) norm_cols[[nm]] <- rep(NA_real_, offs[length(offs)])
        norm_cols[[nm]][id(hw$bins) + 1] <- unname(hw$vec[[nm]])
      }
      if (nrow(hw$raw)) {
        b1 <- id(pmin(hw$raw$x, hw$raw$y)); b2 <- id(pmax(hw$raw$x, hw$raw$y))
        px[[w]] <- data.frame(b1 = b1, b2 = b2, count = hw$raw$counts)
      }
    }
    px <- do.call(rbind, px)
    px <- px[px$count > 0, ]
    # One zero-count contact on the file's last bin: plotgardener::readCool() reads one
    # pixel past the last bin of a view, which fails ("index exceeds HDF5-array
    # dimension") when no contact comes after it in the file, as in the last window.
    nbins <- offs[length(offs)]
    px <- rbind(px, data.frame(b1 = nbins - 1, b2 = nbins - 1, count = 0))
    px <- px[order(px$b1, px$b2), ]
    stopifnot(!anyDuplicated(px[, c("b1", "b2")]))
    is_int <- all(px$count == round(px$count))
    bin1_offset <- c(0L, cumsum(tabulate(px$b1 + 1, nbins)))
    write_ds(dest, file.path(g, "chroms/name"), ch$name, "character")
    write_ds(dest, file.path(g, "chroms/length"), as.integer(ch$length), "integer")
    write_ds(dest, file.path(g, "bins/chrom"), as.integer(b_chrom), "integer")
    write_ds(dest, file.path(g, "bins/start"), as.integer(b_start), "integer")
    write_ds(dest, file.path(g, "bins/end"), as.integer(b_end), "integer")
    for (nm in names(norm_cols)) write_ds(dest, file.path(g, "bins", nm), norm_cols[[nm]], "double")
    write_ds(dest, file.path(g, "pixels/bin1_id"), as.integer(px$b1), "integer")
    write_ds(dest, file.path(g, "pixels/bin2_id"), as.integer(px$b2), "integer")
    write_ds(dest, file.path(g, "pixels/count"), if (is_int) as.integer(px$count) else px$count,
             if (is_int) "integer" else "double")
    write_ds(dest, file.path(g, "indexes/chrom_offset"), as.integer(offs), "integer")
    write_ds(dest, file.path(g, "indexes/bin1_offset"), as.integer(bin1_offset), "integer")
    fid <- H5Fopen(dest)
    gid5 <- H5Gopen(fid, g)
    h5writeAttribute("HDF5::Cooler", gid5, "format")
    h5writeAttribute(3L, gid5, "format-version")
    h5writeAttribute("fixed", gid5, "bin-type")
    h5writeAttribute(as.integer(res), gid5, "bin-size")
    h5writeAttribute("symmetric-upper", gid5, "storage-mode")
    h5writeAttribute("mm39", gid5, "genome-assembly")
    h5writeAttribute(as.integer(nbins), gid5, "nbins")
    h5writeAttribute(length(ch$name), gid5, "nchroms")
    h5writeAttribute(nrow(px), gid5, "nnz")
    H5Gclose(gid5)
    H5Fclose(fid)
    cat(sprintf("  %s bp: %d contacts, norms %s\n", format(res, scientific = FALSE), nrow(px),
                paste(names(norm_cols), collapse = " ")))
  }
  h5closeAll()
}

# The cropped .mcool reads back the .hic's values in every window, at every bin size and
# normalisation (the pixels with a finite value; straw leaves out what readCool() drops).
check_mcool <- function(hic, cool) {
  for (res in hic_resolutions(hic)) for (w in seq_along(HIC_WINDOWS)) {
    chrom <- as.character(seqnames(HIC_WINDOWS[w]))
    s <- ((start(HIC_WINDOWS[w]) - 1L) %/% res) * res
    e <- ((end(HIC_WINDOWS[w]) - 1L) %/% res) * res
    for (nm in c("NONE", NORMS)) {
      a <- tryCatch(strawr::straw(nm, hic, paste0(chrom, ":", format(s, scientific = FALSE), ":",
                                                 format(e, scientific = FALSE)),
                                  paste0(chrom, ":", format(s, scientific = FALSE), ":",
                                         format(e, scientific = FALSE)), "BP", res),
                    error = function(err) NULL)
      if (is.null(a)) next
      # chromend e + 1: the same last bin, and never 0 bp wide when the window is one bin.
      b <- suppressWarnings(plotgardener::readCool(cool, chrom = chrom, chromstart = s, chromend = e + 1,
                                                   resolution = res, norm = nm, quiet = TRUE))
      colnames(b) <- c("x", "y", "counts")
      a <- a[is.finite(a$counts) & a$counts != 0, ]
      a <- a[order(a$x, a$y), ]; b <- b[order(b$x, b$y), ]
      stopifnot(nrow(a) == nrow(b), all(a$x == b$x), all(a$y == b$y),
                max(abs(b$counts / a$counts - 1), 0) < 1e-5)
    }
  }
}

# ---- Per condition --------------------------------------------------------------------
for (cond in names(COND)) {
  cat(cond, "\n")
  for (t in c("rna", "atac", "k36", "igg"))
    cat(sprintf("  %s: %d bins\n", t, crop_bw(src_path(t, cond), out_path(t, cond))))
  cat(sprintf("  junctions: %d\n", crop_lines(src_path("junc", cond), out_path("junc", cond))))
  for (t in c("atac_pk", "k36_pk"))
    cat(sprintf("  %s: %d peaks\n", t, crop_lines(src_path(t, cond), out_path(t, cond))))

  q <- readLines(src_path("quant", cond))
  q_name <- sub("\t.*", "", q)
  q_kept <- q[c(1, which(q_name %in% tx_ids))]
  # A GTF transcript missing from the Salmon index would have no row here either.
  writeLines(q_kept, out_path("quant", cond))
  cat(sprintf("  quant.sf: %d transcripts (%d of the GTF's have no row in the source)\n",
              length(q_kept) - 1, length(setdiff(tx_ids, q_name))))

  write_mcool(src_path("hic", cond), out_path("hic", cond))
  check_mcool(src_path("hic", cond), out_path("hic", cond))
  cat("  Hi-C: matches the .hic in every window, bin size and normalisation\n")
}

sz <- file.size(list.files(OUT, recursive = TRUE, full.names = TRUE))
cat(sprintf("Total: %.2f MB in %d files\n", sum(sz) / 1e6, length(sz)))
