# Gene symbol → genomic coordinates via GENCODE TxDb + an OrgDb package.
# Returns list(chrom, start, end) with nice proportional flanking, or NULL on lookup failure.
# `orgdb_pkg` defaults to the mouse package for callers that predate Phase
# C.0's genome config (utils/bench_slow_viz.R) -- app.R's own call site always
# passes ORGDB_PKG explicitly (the active genome, see app.R's GENOME).

nice_flank <- function(raw_start, raw_end, gene_start, gene_end) {
  # Find the largest "display unit" d (from 1/2/5 × 10^k series) such that
  # ceiling(raw_start/d)*d < gene_start  AND  floor(raw_end/d)*d > gene_end.
  # Rounding INWARD (ceiling on the start, floor on the end) rather than
  # outward guarantees s >= raw_start and e <= raw_end always -- the nice
  # coordinates can never expand the view past the configured flanking %,
  # only tidy it up within that bound. Falls back to the exact raw bounds
  # if no nice unit fits inside them.
  range_bp <- raw_end - raw_start
  mags     <- 10^(0:9)
  units    <- sort(unique(c(mags, mags * 2, mags * 5)), decreasing = TRUE)
  units    <- units[units <= range_bp]
  for (d in units) {
    s <- ceiling(raw_start / d) * d
    e <- floor(raw_end / d) * d
    if (s < gene_start && e > gene_end)
      return(list(start = as.integer(max(s, 1L)), end = as.integer(e)))
  }
  list(start = as.integer(max(raw_start, 1L)), end = as.integer(raw_end))
}

gene_to_coords <- function(gene_symbol, txdb, flank_pct = 0.15, orgdb_pkg = "org.Mm.eg.db") {
  res <- tryCatch(
    AnnotationDbi::select(
      get(orgdb_pkg, envir = asNamespace(orgdb_pkg)),
      keys    = gene_symbol,
      keytype = "SYMBOL",
      columns = "ENSEMBL"
    ),
    error = function(e) NULL
  )
  if (is.null(res) || nrow(res) == 0) return(NULL)

  ensembl_id <- res$ENSEMBL[!is.na(res$ENSEMBL)][1]
  if (length(ensembl_id) == 0 || is.na(ensembl_id)) return(NULL)

  gr <- tryCatch(
    GenomicFeatures::genes(txdb, filter = list(gene_id = ensembl_id)),
    error = function(e) NULL
  )
  if (is.null(gr) || length(gr) == 0) return(NULL)

  gene_start <- GenomicRanges::start(gr)[1]
  gene_end   <- GenomicRanges::end(gr)[1]
  flank      <- round((gene_end - gene_start) * flank_pct)

  coords <- nice_flank(
    raw_start  = gene_start - flank,
    raw_end    = gene_end   + flank,
    gene_start = gene_start,
    gene_end   = gene_end
  )

  list(
    chrom = as.character(GenomicRanges::seqnames(gr))[1],
    start = coords$start,
    end   = coords$end
  )
}
