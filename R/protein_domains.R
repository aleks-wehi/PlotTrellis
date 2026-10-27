# rationale: see solutions/concepts/19-relocated-code-rationale.md, "utils/protein_domains.R:1"

PROTEIN_DOMAINS_CACHE_DIR <- file.path(tools::R_user_dir("PlotTrellis", which = "cache"), "protein_domains")
PD_TYPE_COLOURS           <- c(DOMAIN = "#4C72B0", REGION = "#DD8452", MOTIF = "#55A868")

# NCBI taxon ID per orgdb this app actually supports -- UniProt's search API
# requires organism_id, and there's no generic orgdb->taxon accessor worth
# adding for just two species (mirrors GENOME_PRESETS' own "only what's
# actually supported" convention in utils/genomes.R). Extend when a new
# genome's orgdb is added.
ORGDB_TAXON_ID <- c(org.Mm.eg.db = 10090L, org.Hs.eg.db = 9606L)

.pd_strip_version <- function(x) sub("\\.[0-9]+$", "", x)

# All three fetches below go through .uniprot_read_lines() (uniprot_canonical.R):
# 10s timeout, and a network failure raises "uniprot_unreachable" instead of
# looking like an empty "no reviewed entry" answer.
.pd_uniprot_search_by_gene <- function(gene_symbol, organism_id = 10090) {
  url <- sprintf(
    "https://rest.uniprot.org/uniprotkb/search?query=gene:%s+AND+organism_id:%d+AND+reviewed:true&fields=accession&format=list",
    gene_symbol, organism_id)
  hits <- .uniprot_read_lines(url)
  if (length(hits) == 0) NA_character_ else hits[1]
}
.pd_fetch_uniprot_record <- function(acc) {
  .uniprot_read_lines(sprintf("https://rest.uniprot.org/uniprotkb/%s.txt", acc))
}
.pd_fetch_uniprot_fasta_len <- function(acc) {
  lines <- .uniprot_read_lines(sprintf("https://rest.uniprot.org/uniprotkb/%s.fasta", acc))
  if (length(lines) < 2) return(NA_integer_)
  nchar(paste(lines[-1], collapse = ""))
}
.pd_parse_canonical_isoform <- function(txt, acc) {
  alt_start <- grep("^CC\\s+-!- ALTERNATIVE PRODUCTS:", txt)
  if (length(alt_start) == 0) return(list(isoform_id = acc, name = NA_character_, n_isoforms = 1L))
  block_end <- grep("^CC   -!-", txt)
  block_end <- block_end[block_end > alt_start[1]]
  block_end <- if (length(block_end) == 0) length(txt) else block_end[1] - 1
  block <- txt[alt_start[1]:block_end]
  n_isoforms <- suppressWarnings(as.integer(sub(".*Named isoforms=(\\d+).*", "\\1",
                grep("Named isoforms=", block, value = TRUE)[1])))
  disp_idx  <- grep("Sequence=Displayed", block)
  if (length(disp_idx) == 0) return(list(isoform_id = acc, name = NA_character_, n_isoforms = n_isoforms))
  isoform_id <- sub(".*IsoId=([^;]+);.*", "\\1", block[disp_idx[1]])
  name <- sub(".*Name=([^;]+);?.*", "\\1", block[disp_idx[1] - 1])
  list(isoform_id = isoform_id, name = name, n_isoforms = n_isoforms)
}
.pd_parse_ccds_xrefs <- function(txt, acc) {
  lines <- grep("^DR   CCDS;", txt, value = TRUE)
  if (length(lines) == 0) return(data.frame(ccds_id = character(0), isoform_id = character(0)))
  ccds_id     <- sub("^DR   CCDS; ([^;]+);.*", "\\1", lines)
  has_bracket <- grepl("\\[", lines)
  isoform_id  <- ifelse(has_bracket, sub(".*\\[([^]]+)\\].*", "\\1", lines), acc)
  data.frame(ccds_id = .pd_strip_version(ccds_id), isoform_id = isoform_id)
}
.pd_parse_ft_features <- function(txt, types = c("DOMAIN", "REGION", "MOTIF")) {
  pattern <- sprintf("^FT   (%s)\\s+(\\d+)\\.\\.(\\d+)", paste(types, collapse = "|"))
  idx <- grep(pattern, txt)
  out <- data.frame(type = character(0), start = integer(0), end = integer(0), note = character(0))
  for (i in idx) {
    m <- regmatches(txt[i], regexec(pattern, txt[i]))[[1]]
    note_line <- txt[i + 1]
    note <- if (grepl("/note=", note_line)) sub('.*/note="([^"]*)".*', "\\1", note_line) else NA_character_
    out <- rbind(out, data.frame(type = m[2], start = as.integer(m[3]), end = as.integer(m[4]), note = note))
  }
  out
}
# Gene's transcript IDs + ccdsid, via the same txdb-then-precomputed-cache
# path transcript_annotations() (transcript_annotation.R) uses, instead of a
# per-call zgrep over the GTF -- reuses the tx_meta_lookup RDS that
# build_tx_meta_lookup() already built once for the whole genome.
.pd_gene_transcript_tags <- function(gene_symbol, gtf_path, txdb, orgdb) {
  ens_res <- AnnotationDbi::select(get(orgdb, envir = asNamespace(orgdb)),
    keys = gene_symbol, keytype = "SYMBOL", columns = "ENSEMBL")
  ens_id <- ens_res$ENSEMBL[!is.na(ens_res$ENSEMBL)][1]
  if (length(ens_id) == 0 || is.na(ens_id))
    return(data.frame(transcript_id = character(0), ccdsid = character(0)))

  q <- AnnotationDbi::select(txdb, keys = ens_id, keytype = "GENEID", columns = "TXNAME")
  tx_ids <- unique(q$TXNAME)
  if (length(tx_ids) == 0)
    return(data.frame(transcript_id = character(0), ccdsid = character(0)))

  meta_lookup <- readRDS(.genome_cache_path(gtf_path, "_tx_meta.rds"))
  meta <- meta_lookup[sub("\\.[0-9]+$", "", tx_ids), , drop = FALSE]
  data.frame(transcript_id = tx_ids, ccdsid = meta$ccds_id, stringsAsFactors = FALSE)
}
.pd_resolve_canonical_transcript <- function(gene_symbol, acc, txt, cds_by_tx, gtf_path, txdb, orgdb) {
  canon <- .pd_parse_canonical_isoform(txt, acc)
  ccds  <- .pd_parse_ccds_xrefs(txt, acc)
  tags  <- .pd_gene_transcript_tags(gene_symbol, gtf_path, txdb, orgdb)

  canon_ccds <- ccds$ccds_id[ccds$isoform_id == canon$isoform_id]
  id_join_tx <- tags$transcript_id[!is.na(tags$ccdsid) & .pd_strip_version(tags$ccdsid) %in% canon_ccds]
  if (length(id_join_tx) >= 1) {
    return(list(transcript_id = id_join_tx[1], method = "ccds_id_join", canonical_isoform = canon))
  }

  protein_len <- .pd_fetch_uniprot_fasta_len(acc)
  tx_ids    <- tags$transcript_id[tags$transcript_id %in% names(cds_by_tx)]
  tx_aa_len <- vapply(tx_ids, function(t) sum(GenomicRanges::width(cds_by_tx[[t]])) / 3 - 1, numeric(1))
  matched_tx <- tx_ids[!is.na(protein_len) & tx_aa_len == protein_len]
  if (length(matched_tx) == 0) {
    # Neither method resolved a transcript. This is deterministic given the
    # same GTF/UniProt/CCDS inputs (not a network hiccup), so it's returned
    # as a normal outcome rather than thrown -- an exception here would land
    # in gene_protein_domains()'s catch-all and get the same "network
    # issue? will retry" message as an actual transient fetch failure, which
    # is wrong: retrying changes nothing for this case.
    return(list(transcript_id = NA_character_, method = "unresolved", canonical_isoform = canon))
  }
  # Multiple transcripts can share the exact same CDS amino-acid length --
  # matched_tx[1] still has to pick one, but callers need to know the pick
  # was arbitrary rather than a confirmed unique match.
  method <- if (length(matched_tx) > 1) "length_match_fallback_ambiguous" else "length_match_fallback"
  list(transcript_id = matched_tx[1], method = method, canonical_isoform = canon)
}
.pd_nt_to_genomic_range <- function(cds_gr, nt_start, nt_end) {
  w <- GenomicRanges::width(cds_gr); cum_end <- cumsum(w); cum_start <- cum_end - w + 1
  segs <- list()
  for (i in seq_along(cds_gr)) {
    ov_start <- max(nt_start, cum_start[i]); ov_end <- min(nt_end, cum_end[i])
    if (ov_start > ov_end) next
    local_start <- ov_start - cum_start[i] + 1; local_end <- ov_end - cum_start[i] + 1
    ex <- cds_gr[i]
    if (as.character(GenomicRanges::strand(ex)) == "+") {
      g_start <- GenomicRanges::start(ex) + local_start - 1; g_end <- GenomicRanges::start(ex) + local_end - 1
    } else {
      g_start <- GenomicRanges::end(ex) - local_end + 1; g_end <- GenomicRanges::end(ex) - local_start + 1
    }
    segs[[length(segs) + 1]] <- GenomicRanges::GRanges(GenomicRanges::seqnames(ex), IRanges::IRanges(g_start, g_end),
                                                          strand = GenomicRanges::strand(ex))
  }
  do.call(c, segs)
}
.pd_aa_range_to_nt <- function(aa_start, aa_end) c((aa_start - 1) * 3 + 1, aa_end * 3)
.pd_map_domains_to_genomic <- function(domains, cds_gr) {
  segs_list <- lapply(seq_len(nrow(domains)), function(i) {
    d <- domains[i, ]; nt <- .pd_aa_range_to_nt(d$start, d$end); g <- .pd_nt_to_genomic_range(cds_gr, nt[1], nt[2])
    GenomicRanges::mcols(g)$domain_id <- i
    GenomicRanges::mcols(g)$type      <- d$type
    GenomicRanges::mcols(g)$note      <- d$note
    GenomicRanges::mcols(g)$aa_start  <- d$start
    GenomicRanges::mcols(g)$aa_end    <- d$end
    g
  })
  do.call(c, segs_list)
}

# `orgdb` (e.g. "org.Hs.eg.db") is required, not defaulted -- a silent mouse default here is exactly the bug that made every human-dataset gene resolve against the wrong species' UniProt record (found signing off Phase D.2, 2026-09-01). Full rationale: solutions/concepts/19-relocated-code-rationale.md, "utils/protein_domains.R:170".
gene_protein_domains <- function(gene_symbol, txdb, gtf_path, orgdb, genome_label,
                                  cache_dir = PROTEIN_DOMAINS_CACHE_DIR) {
  if (!orgdb %in% names(ORGDB_TAXON_ID))
    stop("No taxon ID mapping for orgdb '", orgdb, "' in ORGDB_TAXON_ID (utils/protein_domains.R) -- add one.")
  organism_id <- ORGDB_TAXON_ID[[orgdb]]

  # Cache is scoped by genome_label, not orgdb -- orgdb only distinguishes
  # species, not build (e.g. mm10 vs mm39 are both org.Mm.eg.db), and the
  # same case-sensitive gene_symbol string is not guaranteed unique across
  # species either (gene-symbol casing conventions differ by species but
  # aren't enforced) -- mixing species OR builds under one cache key is
  # exactly how the PAX5 bug above shipped (2026-09-01 fix scoped by orgdb;
  # 2026-09-08 GTF-consistency audit found that still collides across builds
  # sharing one orgdb, tightened to genome_label here).
  cache_file <- file.path(cache_dir, genome_label, paste0(gene_symbol, ".rds"))
  if (file.exists(cache_file)) return(readRDS(cache_file))

  result <- tryCatch({
    acc <- .pd_uniprot_search_by_gene(gene_symbol, organism_id = organism_id)
    if (is.na(acc)) return(list(status = "no_uniprot_entry", gene = gene_symbol))
    txt <- .pd_fetch_uniprot_record(acc)
    if (length(txt) == 0) return(list(status = "error", gene = gene_symbol, message = "UniProt record fetch returned nothing"))
    domains <- .pd_parse_ft_features(txt)
    if (nrow(domains) == 0) return(list(status = "no_domains", gene = gene_symbol, acc = acc))

    cds_by_tx <- GenomicFeatures::cdsBy(txdb, by = "tx", use.names = TRUE)
    resolved  <- .pd_resolve_canonical_transcript(gene_symbol, acc, txt, cds_by_tx, gtf_path, txdb, orgdb)
    if (is.na(resolved$transcript_id)) return(list(status = "no_canonical_transcript", gene = gene_symbol, acc = acc))
    cds_gr    <- cds_by_tx[[resolved$transcript_id]]
    domain_gr <- .pd_map_domains_to_genomic(domains, cds_gr)
    domain_df <- as.data.frame(domain_gr)

    list(status = "ok", gene = gene_symbol, acc = acc, canonical_tx = resolved$transcript_id,
         method = resolved$method, canonical_isoform_name = resolved$canonical_isoform$name,
         domains = domain_df)
  }, uniprot_unreachable = function(e) list(status = "uniprot_unreachable", gene = gene_symbol, message = conditionMessage(e)),
     error = function(e) list(status = "error", gene = gene_symbol, message = conditionMessage(e)))

  if (result$status %in% c("ok", "no_uniprot_entry", "no_domains", "no_canonical_transcript")) {
    dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)
    tmp <- tempfile(tmpdir = cache_dir)
    saveRDS(result, tmp)
    file.rename(tmp, cache_file)
  }
  result
}
