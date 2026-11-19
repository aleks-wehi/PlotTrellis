# Data answers of the main screen: the project's samples and track types, gene
# search, a gene or region to a view, the genes in a view, the Gene info of one
# gene, and a render request's tracks to TRACK_ITEMS keys. Plain R lists, no
# Shiny: the page's Shiny wiring sends them as JSON (jsonlite::toJSON(x,
# auto_unbox = TRUE, null = "null", digits = NA)). Every function takes the
# composer's startup context `ctx` (.pt_composer_context()).

# Gene indexes by genome: built once per TxDb file and OrgDb, since mapping
# every gene to its symbol and every transcript to its lengths takes seconds on
# a whole genome.
.pt_main_index_memo <- new.env(parent = emptyenv())

# Biologist-facing errors: stop() without the call, which only a developer could
# read.
.pt_main_stop <- function(...) stop(..., call. = FALSE)

# list(genes, tx). genes: one row per single-strand gene of the TxDb, with ens
# (Ensembl id), chr, s, e, strand, sym (OrgDb symbol, NA when it has none),
# display (sym, else its GENCODE transcripts' gene name, else ens) and coding
# (has a transcript with a CDS). tx: transcriptLengths() with CDS lengths, for
# the Gene info table.
.pt_main_gene_index <- function(ctx) {
  path <- ctx$GENOME$txdb_path
  key <- paste(path, as.numeric(file.mtime(path)), ctx$ORGDB_PKG, sep = "|")
  hit <- .pt_main_index_memo[[key]]
  if (!is.null(hit)) return(hit)
  message("Indexing genes...")
  txdb <- ctx$txdb
  gr <- suppressMessages(
    GenomicFeatures::genes(txdb, single.strand.genes.only = TRUE)
  )
  genes <- data.frame(ens = names(gr),
                      chr = as.character(GenomicRanges::seqnames(gr)),
                      s = GenomicRanges::start(gr), e = GenomicRanges::end(gr),
                      strand = as.character(GenomicRanges::strand(gr)),
                      stringsAsFactors = FALSE)
  db <- get(ctx$ORGDB_PKG, envir = asNamespace(ctx$ORGDB_PKG))
  ids <- intersect(genes$ens, AnnotationDbi::keys(db, "ENSEMBL"))
  map <- if (length(ids)) suppressMessages(AnnotationDbi::select(
    db, keys = ids, keytype = "ENSEMBL", columns = "SYMBOL"
  )) else data.frame(ENSEMBL = character(0), SYMBOL = character(0))
  map <- map[!is.na(map$SYMBOL), ]
  # An id with several symbols takes the first in C-locale order, as
  # resolve_gene_query() does.
  map <- map[order(map$ENSEMBL, map$SYMBOL, method = "radix"), ]
  map <- map[!duplicated(map$ENSEMBL), ]
  genes$sym <- map$SYMBOL[match(genes$ens, map$ENSEMBL)]
  tl <- GenomicFeatures::transcriptLengths(txdb, with.cds_len = TRUE)
  tl$gene_id <- as.character(tl$gene_id)
  # A gene without an OrgDb symbol is named after its transcripts' GENCODE name
  # minus the "-201" suffix, else its id, so the view's gene list never shows a
  # blank.
  tx_names <- readRDS(ctx$GENOME$tx_name_lookup_path)
  nm <- unname(tx_names[sub("\\.[0-9]+$", "", tl$tx_name)])
  gtf_name <- tapply(sub("-[0-9]+$", "", nm), tl$gene_id,
                     function(x) x[!is.na(x)][1])
  genes$display <- ifelse(is.na(genes$sym), unname(gtf_name[genes$ens]),
                          genes$sym)
  genes$display[is.na(genes$display)] <- genes$ens[is.na(genes$display)]
  genes$coding <- genes$ens %in% tl$gene_id[!is.na(tl$cds_len) &
                                              tl$cds_len > 0]
  rownames(genes) <- NULL
  out <- list(genes = genes, tx = tl)
  # One index is kept: the last project's (a rebuilt TxDb or another project
  # in the same R session replaces it).
  rm(list = ls(.pt_main_index_memo), envir = .pt_main_index_memo)
  assign(key, out, envir = .pt_main_index_memo)
  out
}

# Species, build and annotation joined by middle dots ("Mouse, mm10, GENCODE
# vM25"), for the page's header.
.pt_main_genome_text <- function(ctx) {
  species <- c(org.Mm.eg.db = "Mouse", org.Hs.eg.db = "Human")[ctx$ORGDB_PKG]
  if (is.na(species)) species <- ctx$ORGDB_PKG
  gtf <- sub("\\.annotation.*$|\\.gtf.*$", "",
             basename(ctx$CFG$gtf_path %||% ""))
  gtf <- sub("^gencode\\.", "GENCODE ", gtf)
  paste(c(unname(species), ctx$GENOME$genome_label, if (nzchar(gtf)) gtf),
        collapse = " \u00B7 ")
}

# Genes offered before anything is typed: the start view's gene, then the named
# protein-coding genes nearest to it, up to 8. The package keeps no list of
# well-known genes per species, so a project without a start gene suggests none.
.pt_main_suggested <- function(ctx) {
  sym <- ctx$START_VIEW$gene
  if (length(sym) != 1 || is.na(sym) || !nzchar(sym)) return(character(0))
  g <- .pt_main_gene_index(ctx)$genes
  i <- match(sym, g$sym)
  if (is.na(i)) return(character(0))
  near <- g[g$chr == g$chr[i] & g$coding & !is.na(g$sym) & g$sym != sym, ]
  d <- pmax(0, near$s - g$e[i], g$s[i] - near$e)
  unique(c(sym, head(near$sym[order(d, near$sym, method = "radix")], 7)))
}

# One entry per track type, in TRACK_ITEMS order (per-sample types in detection
# order, then annotations, then comparisons). id is the track name pt_tracks()
# lists; cat groups the page's picker (ann, cmp, expr, else the data column's
# name); viz is the drawing module; samples lists the samples that have it (a
# comparison: every sample in one of its pairs).
.pt_main_types <- function(ctx) {
  items <- ctx$TRACK_ITEMS
  tab <- .pt_track_table(items)
  key_track <- setNames(tab$track, tab$key)
  needs_gene <- ctx$GENE_REQUIRED_TYPES %||% .pt_gene_required_types()
  out <- list()
  for (tr in unique(unname(key_track[names(items)]))) {
    rows <- tab[tab$track == tr, ]
    it <- items[[rows$key[1]]]
    kind <- rows$kind[1]
    cat_ <- if (kind == "comparison") "cmp" else if (kind == "fixed") "ann" else
      if (it$type == "tpm_usage") "expr" else it$name %||% it$type
    t <- list(id = tr, name = sub(" \u00B7 .*$", "", it$label), cat = cat_,
              kind = kind, gene = it$type %in% needs_gene, viz = it$type)
    if (kind == "sample") t$samples <- I(rows$sample)
    if (kind == "comparison")
      t$samples <- I(unique(unlist(lapply(items[rows$key], `[[`, "samples"))))
    out[[length(out) + 1]] <- t
  }
  out
}

# GET /api/project minus id: genome text, suggested genes, samples (name,
# colour, in sample order), track types and chromosome lengths. The page always
# opens empty: a project's start view (pt_preload(), a samples.txt's start
# lines) only names the gene the search box suggests first.
.pt_main_project <- function(ctx) {
  cols <- ctx$DEFAULT_SAMPLE_COLS
  list(genome = .pt_main_genome_text(ctx),
       suggested = I(.pt_main_suggested(ctx)),
       samples = lapply(names(cols), function(n)
         list(name = n, colour = unname(cols[[n]]))),
       types = .pt_main_types(ctx),
       chrLen = as.list(stats::setNames(as.numeric(ctx$CHR_LENGTHS),
                                        names(ctx$CHR_LENGTHS))))
}

# Up to 8 genes for the search box: case-insensitive symbol prefix, or an
# Ensembl id prefix when q starts with ENS; the exact match first, then
# protein-coding, then shorter symbols. An empty q gives the suggested genes.
.pt_main_search <- function(ctx, q) {
  g <- .pt_main_gene_index(ctx)$genes
  g <- g[!is.na(g$sym), ]
  q <- trimws(q %||% "")
  if (!nzchar(q)) {
    hit <- g[match(.pt_main_suggested(ctx), g$sym), ]
  } else {
    hit <- g[startsWith(toupper(g$sym), toupper(q)), ]
    if (grepl("^ENS", q, ignore.case = TRUE))
      hit <- rbind(hit, g[startsWith(g$ens, toupper(q)), ])
    hit <- hit[order(toupper(hit$sym) != toupper(q), !hit$coding,
                     nchar(hit$sym), hit$sym, method = "radix"), ]
  }
  hit <- head(hit[!duplicated(hit$ens), ], 8)
  lapply(seq_len(nrow(hit)), function(i)
    list(sym = hit$sym[i], chr = hit$chr[i], s = hit$s[i], e = hit$e[i],
         ens = hit$ens[i]))
}

# A gene query (symbol in any case, Ensembl id, or the project GTF's name for a
# gene the OrgDb has no symbol for) -> its index row as a list plus symbol, or
# list(error).
.pt_main_gene_row <- function(ctx, query) {
  q <- resolve_gene_query(query, ctx$ORGDB_PKG, ctx$ALL_GENES,
                          .pt_gene_names(ctx$txdb))
  if (is.na(q$symbol))
    return(list(error = q$error %||% paste0("'", query, "' is not a gene in ",
                                            ctx$GENOME$genome_label, ".")))
  g <- .pt_main_gene_index(ctx)$genes
  ids <- q$ensembl %||% .pt_gene_ensembl(q$symbol, ctx$ORGDB_PKG, ctx$txdb)
  i <- match(intersect(ids, g$ens)[1], g$ens)
  if (is.na(i))
    return(list(error = paste0(q$symbol,
                               " has no position in this project's annotation.")))
  c(as.list(g[i, ]), list(symbol = q$symbol))
}

# "chr2:103,400,000-103,460,000" -> list(chr, s, e), else NULL.
.pt_main_parse_region <- function(q) {
  m <- regmatches(q, regexec(
    "^\\s*([^:\\s]+)\\s*:\\s*([0-9,.]+)\\s*-\\s*([0-9,.]+)\\s*$", q,
    perl = TRUE))[[1]]
  if (length(m) != 4) return(NULL)
  # whole positions: "1.5" is read as 2
  list(chr = m[2], s = round(as.numeric(gsub(",", "", m[3]))),
       e = round(as.numeric(gsub(",", "", m[4]))))
}

# A chromosome's length, NA for one the project does not have. A chromosome
# with genes but no entry in the chromosome sizes (a scaffold) counts as
# reaching its last gene's end, so its genes can still be drawn, as on the
# classic screen.
.pt_main_chr_length <- function(ctx, chr) {
  if (!is.character(chr) || length(chr) != 1 || is.na(chr)) return(NA_real_)
  if (chr %in% names(ctx$CHR_LENGTHS))
    return(as.numeric(ctx$CHR_LENGTHS[[chr]]))
  g <- .pt_main_gene_index(ctx)$genes
  if (chr %in% g$chr) as.numeric(max(g$e[g$chr == chr])) else NA_real_
}

# GET /api/locus: a gene or a region -> list(kind, label, gene, chr, s, e,
# chrLen, view). A gene's view is the app's default window (the project's gene
# flank, rounded inward, .pt_resolve_gene_region()); a region's is the region.
# `cfg` is the live render settings (its gene flank).
.pt_main_locus <- function(ctx, q, cfg = ctx$RENDER_CFG) {
  q <- trimws(q %||% "")
  if (!nzchar(q))
    .pt_main_stop("Type a gene or a region such as chr2:103,400,000-103,460,000.")
  r <- .pt_main_parse_region(q)
  if (!is.null(r)) {
    len <- .pt_main_chr_length(ctx, r$chr)
    if (is.na(len))
      .pt_main_stop("No chromosome '", r$chr, "' in ", ctx$GENOME$genome_label,
                    ".")
    err <- .pt_region_error(r$s, r$e)
    if (!is.null(err)) .pt_main_stop(err)
    if (r$s < 1 || r$e > len)
      .pt_main_stop("The region must lie within ", r$chr, " (1 to ",
                    format(len, big.mark = ",", scientific = FALSE), ").")
    fmt <- function(x) format(x, big.mark = ",", scientific = FALSE)
    return(list(kind = "region",
                label = sprintf("%s:%s-%s", r$chr, fmt(r$s), fmt(r$e)),
                gene = NULL, chr = r$chr, s = r$s, e = r$e, chrLen = len,
                view = list(s = r$s, e = r$e)))
  }
  g <- .pt_main_gene_row(ctx, q)
  if (!is.null(g$error)) .pt_main_stop(g$error)
  res <- .pt_resolve_gene_region(q, ctx$txdb,
                                 (cfg$gene_flank_pct %||% 15) / 100,
                                 ctx$ORGDB_PKG, ctx$ALL_GENES,
                                 ctx$GENOME$genome_label)
  if (!is.null(res$error)) .pt_main_stop(res$error)
  list(kind = "gene", label = g$symbol, gene = g$symbol, chr = g$chr, s = g$s,
       e = g$e, chrLen = .pt_main_chr_length(ctx, g$chr),
       view = list(s = res$region$start, e = res$region$end))
}

# GET /api/nearby: the genes overlapping chr:s-e, sorted by start. Over 80, the
# named protein-coding ones are kept first (then other named ones, longest
# first): a 2 Mb view of a gene desert's Gm/Rik models is noise to a biologist.
.pt_main_nearby <- function(ctx, chr, s, e) {
  s <- suppressWarnings(as.numeric(s %||% NA))
  e <- suppressWarnings(as.numeric(e %||% NA))
  if (!is.character(chr) || length(chr) != 1 || length(s) != 1 ||
      length(e) != 1 || !is.finite(s) || !is.finite(e))
    .pt_main_stop("nearby needs chr, s and e.")
  g <- .pt_main_gene_index(ctx)$genes
  g <- g[g$chr == chr & g$e >= s & g$s <= e, ]
  if (nrow(g) > 80) {
    named <- !is.na(g$sym) &
      !grepl("^(Gm[0-9]|[0-9A-Z]{2,}[0-9]{5,}Rik$)", g$display)
    g <- g[order(!(g$coding & named), !named, -(g$e - g$s),
                 method = "radix"), ][1:80, ]
  }
  g <- g[order(g$s), ]
  lapply(seq_len(nrow(g)), function(i)
    list(sym = g$display[i], s = g$s[i], e = g$e[i], strand = g$strand[i]))
}

# NA -> NULL, so toJSON() writes null rather than the string "NA" (its default
# for a character or numeric NA).
.pt_main_or_null <- function(x) {
  if (length(x) != 1 || is.na(x)) NULL else x
}

# GET /api/geneinfo: the gene's transcripts (name, id, biotype, coding status,
# exons, length, aa, tags) and its canonical protein's features, with a UniProt
# search link. `cfg` is the live render settings (its lookup paths, cache
# folder and uniprot_reviewed_only). UniProt is looked up as the app's Gene info
# panel does (gene symbol and organism only, cached in the project); a lookup
# that fails leaves the transcripts as they are, without UniProt tags.
.pt_main_geneinfo <- function(ctx, gene, cfg = ctx$RENDER_CFG) {
  g <- .pt_main_gene_row(ctx, gene %||% "")
  if (!is.null(g$error)) .pt_main_stop(g$error)
  sym <- g$symbol
  orgdb <- cfg$orgdb_pkg %||% ctx$ORGDB_PKG
  cache_root <- .pt_cache_root(cfg$cache_dir %||% ctx$CACHE_DIR)
  ann <- transcript_annotations(sym, ctx$txdb, orgdb, cfg$tx_name_lookup_path,
                                cfg$tx_meta_lookup_path)
  tl <- .pt_main_gene_index(ctx)$tx
  out <- list(gene = sym, chr = g$chr, s = g$s, e = g$e, strand = g$strand)
  if (!is.null(ann) && nrow(ann)) {
    uc <- tryCatch(uniprot_canonical_match(
      sym, ann$tx_id, ctx$CDS_BY_TX, ann$ccds_id, orgdb = orgdb,
      cache_path = file.path(cache_root, "uniprot_canonical_cache.rds")
    ), error = function(e) NULL)
    if (!is.null(uc) && nrow(uc) != nrow(ann)) uc <- NULL
    i <- match(ann$tx_id, tl$tx_name)
    out$transcripts <- lapply(seq_len(nrow(ann)), function(j) {
      tags <- character(0)
      if (isTRUE(ann$ensembl_canonical[j])) tags <- c(tags, "Ensembl canonical")
      if (isTRUE(ann$gencode_basic[j])) tags <- c(tags, "basic")
      if (!is.na(ann$appris_tier[j]))
        tags <- c(tags, paste("APPRIS", ann$appris_tier[j]))
      if (!is.na(ann$ccds_id[j])) tags <- c(tags, "CCDS")
      if (!is.null(uc) && isTRUE(uc$uniprot_match[j]))
        tags <- c(tags, if (isTRUE(uc$uniprot_unique[j]))
          "UniProt canonical" else "UniProt (shared match)")
      # The TxDb's CDS includes the stop codon, unless GENCODE marks the CDS
      # end as incomplete (cds_end_NF).
      cds <- tl$cds_len[i[j]]
      aa <- if (is.na(cds) || cds == 0) NULL else
        cds %/% 3 - !grepl("cds_end_NF", ann$tag[j] %||% "")
      list(name = ann$transcript_name[j], id = ann$tx_id[j],
           biotype = .pt_main_or_null(ann$transcript_type[j]),
           coding = ann$coding_status[j],
           exons = .pt_main_or_null(tl$nexon[i[j]]),
           length = .pt_main_or_null(tl$tx_len[i[j]]),
           aa = aa, tags = I(tags))
    })
  }
  pd <- tryCatch(gene_protein_domains(
    sym, ctx$txdb, cfg$tx_meta_lookup_path, orgdb,
    cfg$genome_label %||% ctx$GENOME$genome_label,
    cache_dir = file.path(cache_root, "protein_domains"),
    cds_by_tx = ctx$CDS_BY_TX
  ), error = function(e) list(status = "error", message = conditionMessage(e)))
  if (!is.null(pd$acc)) {
    prot <- list(accession = pd$acc, status = pd$status)
    if (identical(pd$status, "ok")) {
      cds <- tl$cds_len[match(pd$canonical_tx, tl$tx_name)]
      tx_name <- if (!is.null(ann))
        ann$transcript_name[match(pd$canonical_tx, ann$tx_id)] else NA
      if (is.na(tx_name)) tx_name <- pd$canonical_tx
      iso <- pd$canonical_isoform_name %||% NA
      # The app stores no UniProt protein name: gene, isoform and transcript.
      prot$name <- paste0(sym, if (!is.na(iso)) paste0(" isoform ", iso),
                          " (", tx_name, ")")
      if (!is.na(cds)) prot$length <- cds %/% 3 - 1
      d <- pd$domains
      prot$features <- lapply(unique(d$domain_id), function(x) {
        r <- d[d$domain_id == x, ][1, ]
        list(name = r$note, start = r$aa_start, end = r$aa_end,
             source = paste("UniProt", tolower(r$type)))
      })
      prot$features <- prot$features[order(vapply(prot$features, function(f)
        as.numeric(f$start), numeric(1)))]
    }
    out$protein <- prot
  } else if (!is.null(pd$status)) {
    out$protein_status <- .pd_status_text(pd$status, panel = TRUE) %||%
      pd$status
  }
  out$uniprot_url <- tryCatch(
    uniprot_search_url(sym, orgdb, isTRUE(cfg$uniprot_reviewed_only)),
    error = function(e) NULL)
  out
}

# A render request's tracks -> names(ctx$TRACK_ITEMS) for .pt_build_entries(),
# in the order given (any type any number of times). `picked`: the picked
# sample names; a comparison type takes picked[1] as the comparison and
# picked[2] as the reference, and needs exactly two. `tracks`: a list of
# list(type, sample) (sample only for a per-sample type), or the data frame
# jsonlite makes of one.
.pt_main_labels <- function(ctx, picked, tracks) {
  picked <- as.character(unlist(picked))
  if (is.data.frame(tracks))
    tracks <- lapply(seq_len(nrow(tracks)), function(i)
      as.list(tracks[i, , drop = FALSE]))
  if (!length(tracks)) .pt_main_stop("Add at least one track.")
  tab <- .pt_track_table(ctx$TRACK_ITEMS)
  vapply(tracks, function(t) {
    type <- as.character(t$type %||% "")
    rows <- tab[tab$track == type, ]
    if (length(type) != 1 || !nrow(rows))
      .pt_main_stop("Unknown track type '", paste(type, collapse = ", "), "'.")
    kind <- rows$kind[1]
    if (kind == "fixed") return(rows$key[1])
    if (kind == "comparison") {
      if (length(picked) != 2)
        .pt_main_stop("'", type, "' compares two samples: pick exactly two.")
      k <- rows$key[rows$sample == paste(picked[1], "vs", picked[2])]
      if (!length(k))
        .pt_main_stop("No '", type, "' track for ", picked[1], " against ",
                      picked[2], ".")
      return(k[1])
    }
    sample <- as.character(unlist(t$sample))
    if (length(sample) != 1 || is.na(sample) || !nzchar(sample))
      .pt_main_stop("'", type, "' is drawn per sample: say which sample.")
    k <- rows$key[!is.na(rows$sample) & rows$sample == sample]
    if (!length(k))
      .pt_main_stop("No '", type, "' track for sample '", sample, "'.")
    k[1]
  }, character(1), USE.NAMES = FALSE)
}

# One call for the page's data requests: op in project, search, locus, nearby,
# geneinfo, with that endpoint's arguments in `args` (what the page sent) and
# `cfg`, the live render settings, which the page never supplies. Any error
# becomes list(error = <message>), the shape the page shows.
.pt_main_request <- function(ctx, op, args = list(), cfg = ctx$RENDER_CFG) {
  op <- if (is.character(op) && length(op) == 1 && !is.na(op)) op else ""
  tryCatch({
    if (!op %in% c("project", "search", "locus", "nearby", "geneinfo"))
      .pt_main_stop("Unknown request '", op, "'.")
    switch(op,
      project  = .pt_main_project(ctx),
      search   = .pt_main_search(ctx, args$q),
      locus    = .pt_main_locus(ctx, args$q, cfg),
      nearby   = .pt_main_nearby(ctx, args$chr, args$s, args$e),
      geneinfo = .pt_main_geneinfo(ctx, args$gene, cfg))
  }, error = function(e) list(error = conditionMessage(e)))
}
