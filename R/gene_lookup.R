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

# Ensembl gene ID, any species prefix, optional version suffix (BUGS.md #36).
# A GENCODE "_PAR_Y" copy is accepted and looked up as its chrX gene.
ENSEMBL_GENE_ID_RE <- "^ENS[A-Z]*G[0-9]{6,}(\\.[0-9]+)?(_PAR_Y)?$"

# A user's gene query -> list(symbol, ensembl, error). Accepts an OrgDb SYMBOL
# (exact, else case-insensitive when exactly one symbol matches, BUGS.md #49)
# or an Ensembl gene ID with or without a version (#36). Transcript IDs,
# Entrez IDs and aliases are deliberately not accepted (Phase I.1 step 4a).
# The app keys every track and panel on the symbol, so callers display and
# pass on `symbol`, not the raw query. An ID with several symbols takes the
# first in C-locale order and lists the rest in `also`. `symbols` is the
# OrgDb's SYMBOL keys, passed in when the caller already has them.
resolve_gene_query <- function(query, orgdb_pkg, symbols = NULL) {
  q  <- trimws(query %||% "")
  db <- get(orgdb_pkg, envir = asNamespace(orgdb_pkg))
  if (grepl(ENSEMBL_GENE_ID_RE, q, ignore.case = TRUE)) {
    id  <- toupper(sub("(\\.[0-9]+)?(_PAR_Y)?$", "", q, ignore.case = TRUE))
    res <- tryCatch(suppressMessages(AnnotationDbi::select(db, keys = id, keytype = "ENSEMBL", columns = "SYMBOL")),
                    error = function(e) NULL)
    sym <- sort(unique(as.character(res$SYMBOL[!is.na(res$SYMBOL)])), method = "radix")
    if (length(sym) == 0)
      return(list(symbol = NA_character_, ensembl = id,
                  error = paste0("Ensembl ID '", id, "' has no gene symbol in ", orgdb_pkg)))
    return(list(symbol = sym[1], ensembl = id, error = NULL, also = sym[-1]))
  }
  if (is.null(symbols)) symbols <- AnnotationDbi::keys(db, keytype = "SYMBOL")
  if (q %in% symbols) return(list(symbol = q, ensembl = NULL, error = NULL))
  hits <- symbols[toupper(symbols) == toupper(q)]
  if (length(hits) == 1) return(list(symbol = hits, ensembl = NULL, error = NULL))
  list(symbol = NA_character_, ensembl = NULL,
       error = if (length(hits) > 1) paste0("'", q, "' matches several symbols that differ only in case (",
                                            paste(hits, collapse = ", "), "); type the exact one") else NULL)
}

# `ensembl_id`, when given (the user typed an ID), picks the coordinates;
# otherwise the first of the symbol's Ensembl IDs that exists in the TxDb
# (e.g. PTPRC also maps to an alt-haplotype ID absent from a primary GTF).
gene_to_coords <- function(gene_symbol, txdb, flank_pct = 0.15, orgdb_pkg = "org.Mm.eg.db", ensembl_id = NULL) {
  if (is.null(ensembl_id)) {
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
    ensembl_id <- unique(res$ENSEMBL[!is.na(res$ENSEMBL)])
    if (length(ensembl_id) == 0) return(NULL)
  }

  gr <- tryCatch(
    GenomicFeatures::genes(txdb, filter = list(gene_id = ensembl_id)),
    error = function(e) NULL
  )
  if (is.null(gr) || length(gr) == 0) return(NULL)
  gr <- gr[match(intersect(ensembl_id, names(gr))[1], names(gr))]

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
