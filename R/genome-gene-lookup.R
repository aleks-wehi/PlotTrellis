# Gene symbol → genomic coordinates via GENCODE TxDb + an OrgDb package. Returns
# list(chrom, start, end) with nice proportional flanking, or NULL on lookup
# failure. `orgdb_pkg` is the project genome's OrgDb package; it has no default,
# so a lookup can never silently use another species'.

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

# Ensembl gene ID, any species prefix, optional version suffix.
# A GENCODE "_PAR_Y" copy is accepted and looked up as its chrX gene.
ENSEMBL_GENE_ID_RE <- "^ENS[A-Z]*G[0-9]{6,}(\\.[0-9]+)?(_PAR_Y)?$"

# The project GTF's gene names (build_gene_name_lookup()) kept beside the
# TxDb file `txdb` was opened from; NULL when there are none (a TxDb that is
# not a genome cache's).
.pt_gene_names <- function(txdb) {
  f <- tryCatch(.pt_gene_names_file(AnnotationDbi::dbfile(txdb)),
                error = function(e) "")
  if (length(f) == 1 && grepl("_gene_names\\.rds$", f) && file.exists(f))
    readRDS(f)
}

# The Ensembl gene IDs of `gene`, the symbol every track keys on, that are in
# the TxDb `txdb`: the OrgDb's; when none of those is in the project's GTF (no
# OrgDb symbol, an OrgDb symbol without Ensembl IDs, or IDs the GTF has under
# others), the GTF's, by gene_name or `gene` itself when it is a gene_id.
# character(0) when neither has it.
.pt_gene_ensembl <- function(gene, orgdb_pkg, txdb) {
  res <- tryCatch(
    suppressMessages(AnnotationDbi::select(
      get(orgdb_pkg, envir = asNamespace(orgdb_pkg)),
      keys = gene, keytype = "SYMBOL", columns = "ENSEMBL"
    )),
    error = function(e) NULL
  )
  ids <- unique(as.character(res$ENSEMBL[!is.na(res$ENSEMBL)]))
  if (length(ids)) {
    # The TxDb's own gene table; a TxDb that can't be asked keeps them all.
    known <- tryCatch(RSQLite::dbGetQuery(
      AnnotationDbi::dbconn(txdb),
      "SELECT DISTINCT gene_id FROM gene WHERE gene_id = ?",
      params = list(ids)
    )$gene_id, error = function(e) ids)
    if (any(ids %in% known)) return(ids[ids %in% known])
  }
  gtf <- .pt_gene_names(txdb)
  if (gene %in% names(gtf)) return(gene)
  as.character(names(gtf)[gtf %in% gene])
}

# A user's gene query -> list(symbol, ensembl, error). Accepts an OrgDb SYMBOL
# (exact, else case-insensitive when exactly one symbol matches)
# or an Ensembl gene ID with or without a version (#36). Transcript IDs,
# Entrez IDs and aliases are deliberately not accepted.
# The app keys every track and panel on the symbol, so callers display and
# pass on `symbol`, not the raw query. An ID with several symbols takes the
# first in C-locale order and lists the rest in `also`. `symbols` is the
# OrgDb's SYMBOL keys, passed in when the caller already has them.
# `gene_names` (.pt_gene_names()): the project GTF's, for a gene the OrgDb has
# no symbol for: found by its GTF gene_name (exact, else
# case-insensitive) or gene_id, its symbol is that gene_name, or its gene_id
# when another gene has that name (so a track looking the symbol up finds this
# gene, .pt_gene_ensembl()).
resolve_gene_query <- function(query, orgdb_pkg, symbols = NULL,
                               gene_names = NULL) {
  q  <- trimws(query %||% "")
  db <- get(orgdb_pkg, envir = asNamespace(orgdb_pkg))
  gtf_note <- if (!is.null(gene_names))
    paste0("this project's GTF (", attr(gene_names, "gtf"), ")")
  # The symbol of the GTF's gene `id`: its name, unless another gene has it
  # (in the GTF, or as an OrgDb symbol with Ensembl IDs of its own).
  gtf_symbol <- function(id) {
    nm <- unname(gene_names[id])
    if (sum(gene_names == nm) > 1) return(id)
    res <- tryCatch(
      suppressMessages(AnnotationDbi::select(
        db, keys = nm, keytype = "SYMBOL", columns = "ENSEMBL"
      )),
      error = function(e) NULL
    )
    if (any(!is.na(res$ENSEMBL))) id else nm
  }
  if (grepl(ENSEMBL_GENE_ID_RE, q, ignore.case = TRUE)) {
    id  <- toupper(sub("(\\.[0-9]+)?(_PAR_Y)?$", "", q, ignore.case = TRUE))
    res <- tryCatch(
      suppressMessages(AnnotationDbi::select(
        db, keys = id, keytype = "ENSEMBL", columns = "SYMBOL"
      )),
      error = function(e) NULL
    )
    sym <- sort(unique(as.character(res$SYMBOL[!is.na(res$SYMBOL)])),
                method = "radix")
    if (length(sym) == 0 && id %in% names(gene_names))
      return(list(symbol = gtf_symbol(id), ensembl = id, error = NULL))
    if (length(sym) == 0)
      return(list(symbol = NA_character_, ensembl = id,
                  error = if (is.null(gtf_note))
                    paste0("Ensembl ID '", id, "' has no gene symbol in ",
                           orgdb_pkg)
                  else paste0("Ensembl ID '", id, "' is not in ", gtf_note)))
    return(list(symbol = sym[1], ensembl = id, error = NULL, also = sym[-1]))
  }
  if (is.null(symbols)) symbols <- AnnotationDbi::keys(db, keytype = "SYMBOL")
  if (q %in% symbols) return(list(symbol = q, ensembl = NULL, error = NULL))
  hits <- symbols[toupper(symbols) == toupper(q)]
  if (length(hits) == 1)
    return(list(symbol = hits, ensembl = NULL, error = NULL))
  if (length(hits) == 0 && nzchar(q) && !is.null(gene_names)) {
    ids <- names(gene_names)[gene_names == q]
    if (length(ids) == 0) {
      ids <- names(gene_names)[toupper(gene_names) == toupper(q)]
      if (length(unique(gene_names[ids])) > 1) ids <- character(0)
    }
    if (length(ids) == 1)
      return(list(symbol = gtf_symbol(ids), ensembl = ids, error = NULL))
    if (length(ids) > 1)
      return(list(symbol = NA_character_, ensembl = NULL,
                  error = paste0("Several genes in ", gtf_note, " are named '",
                                 q, "' (", paste(utils::head(ids, 5),
                                                 collapse = ", "),
                                 if (length(ids) > 5) ", ...",
                                 "): enter one's Ensembl gene ID")))
  }
  list(symbol = NA_character_, ensembl = NULL,
       error = if (length(hits) > 1)
         paste0("'", q, "' matches several symbols that differ only in case (",
                paste(hits, collapse = ", "),
                "); type the exact one") else NULL)
}

# `ensembl_id`, when given (the user typed an ID), picks the coordinates;
# otherwise the first of the symbol's Ensembl IDs (.pt_gene_ensembl(): the
# OrgDb's, else the GTF's) that exists in the TxDb (e.g. PTPRC also maps to an
# alt-haplotype ID absent from a primary GTF).
gene_to_coords <- function(gene_symbol, txdb, flank_pct = 0.15, orgdb_pkg,
                           ensembl_id = NULL) {
  if (is.null(ensembl_id)) {
    ensembl_id <- .pt_gene_ensembl(gene_symbol, orgdb_pkg, txdb)
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
