# rationale: see solutions/concepts/19-relocated-code-rationale.md, "utils/uniprot_canonical.R:1"

UNIPROT_CANONICAL_CACHE_PATH <- file.path(
  tools::R_user_dir("PlotTrellis", which = "cache"), "uniprot_canonical_cache.rds"
)

# rationale: see solutions/concepts/19-relocated-code-rationale.md, "utils/uniprot_canonical.R:42"
UNIPROT_CACHE_SCHEMA_VERSION <- 3L

.uniprot_canonical_cache_load <- function() {
  if (!file.exists(UNIPROT_CANONICAL_CACHE_PATH)) return(list())
  cache <- tryCatch(readRDS(UNIPROT_CANONICAL_CACHE_PATH), error = function(e) NULL)
  if (is.null(cache) || is.null(attr(cache, "schema_version")) ||
      attr(cache, "schema_version") != UNIPROT_CACHE_SCHEMA_VERSION) {
    return(list())
  }
  cache
}

.uniprot_canonical_cache_save <- function(cache) {
  attr(cache, "schema_version") <- UNIPROT_CACHE_SCHEMA_VERSION
  dir.create(dirname(UNIPROT_CANONICAL_CACHE_PATH), recursive = TRUE, showWarnings = FALSE)
  saveRDS(cache, UNIPROT_CANONICAL_CACHE_PATH)
}

# NCBI taxon ID + Ensembl-transcript-ID prefix per orgdb this app supports --
# duplicated from utils/protein_domains.R's ORGDB_TAXON_ID (same "just two
# species, no generic accessor" reasoning there) rather than shared. Keep in
# sync with protein_domains.R's copy (and utils/transcript_annotation.R's
# ORGDB_TAXON_ID copy) if any of them change.
ORGDB_TAXON_ID <- c(org.Mm.eg.db = 10090L, org.Hs.eg.db = 9606L)
ORGDB_ENSEMBL_TX_PREFIX <- c(org.Mm.eg.db = "ENSMUST", org.Hs.eg.db = "ENST")

# gene_symbol -> reviewed (Swiss-Prot) UniProt accession, or NA if no reviewed hit; spike 0.5 validated this query for Nfatc1/Klrk1.
uniprot_search_by_gene <- function(gene_symbol, organism_id = 10090) {
  url <- sprintf(
    "https://rest.uniprot.org/uniprotkb/search?query=gene:%s+AND+organism_id:%d+AND+reviewed:true&fields=accession&format=list",
    gene_symbol, organism_id
  )
  hits <- tryCatch(readLines(url, warn = FALSE), error = function(e) character(0))
  if (length(hits) == 0) NA_character_ else hits[1]
}

# Parses one DR cross-reference line (e.g. "DR   Ensembl; ENSMUST00000024747.14; ENSMUSP...; ENSMUSG.... [Q00731-3]") into id (version stripped) + isoform (NA if untagged).
.uniprot_parse_dr <- function(txt, db, id_pattern) {
  lines <- grep(sprintf("^DR   %s;", db), txt, value = TRUE)
  if (length(lines) == 0) return(data.frame(id = character(0), isoform = character(0), stringsAsFactors = FALSE))
  id <- sub(sprintf(".*(%s).*", id_pattern), "\\1", lines)
  has_iso <- grepl("\\[[A-Za-z0-9]+-[0-9]+\\]", lines)
  isoform <- rep(NA_character_, length(lines))
  isoform[has_iso] <- sub(".*\\[[A-Za-z0-9]+-([0-9]+)\\].*", "\\1", lines[has_iso])
  data.frame(id = id, isoform = isoform, stringsAsFactors = FALSE)
}

# Parses the ALTERNATIVE PRODUCTS block's Name=/IsoId= line pairs (confirmed on real records, e.g. Nfatc1/O88942's Alpha/Beta/Alpha'); empty if the entry has no such block (single-isoform gene).
.uniprot_parse_isoforms <- function(txt) {
  name_lines  <- grep("^CC\\s+Name=", txt, value = TRUE)
  isoid_lines <- grep("^CC\\s+IsoId=", txt, value = TRUE)
  if (length(name_lines) == 0 || length(name_lines) != length(isoid_lines)) {
    return(data.frame(number = character(0), name = character(0), isoid = character(0),
                       stringsAsFactors = FALSE))
  }
  name  <- sub("^CC\\s+Name=([^;]+);.*", "\\1", name_lines)
  isoid <- sub(".*IsoId=([A-Za-z0-9]+-[0-9]+);.*", "\\1", isoid_lines)
  data.frame(number = sub(".*-([0-9]+)$", "\\1", isoid), name = name, isoid = isoid,
             stringsAsFactors = FALSE)
}

# `orgdb` is required, not defaulted (2026-09-08 GTF-consistency audit) -- this used to hardcode organism_id:10090 (mouse) via uniprot_search_by_gene()'s own default AND hardcode "ENSMUST" (mouse) as the DR Ensembl id pattern. Full rationale: solutions/concepts/19-relocated-code-rationale.md, "utils/uniprot_canonical.R:123".
.uniprot_fetch_entry <- function(gene_symbol, orgdb) {
  if (!orgdb %in% names(ORGDB_TAXON_ID))
    stop("No taxon ID mapping for orgdb '", orgdb, "' in ORGDB_TAXON_ID (utils/uniprot_canonical.R) -- add one.")
  ens_tx_prefix <- ORGDB_ENSEMBL_TX_PREFIX[[orgdb]]

  empty_map <- data.frame(id = character(0), isoform = character(0), stringsAsFactors = FALSE)
  empty_result <- list(accession = NA_character_, protein_len = NA_integer_,
                        ensembl = empty_map, ccds = empty_map,
                        isoforms = data.frame(number = character(0), name = character(0), isoid = character(0)),
                        canonical_isoform = NA_character_,
                        canonical_isoform_id = NA_character_, canonical_isoform_name = NA_character_)
  acc <- uniprot_search_by_gene(gene_symbol, organism_id = ORGDB_TAXON_ID[[orgdb]])
  if (is.na(acc)) return(empty_result)

  txt <- tryCatch(readLines(sprintf("https://rest.uniprot.org/uniprotkb/%s.txt", acc), warn = FALSE),
                   error = function(e) NULL)
  if (is.null(txt)) { empty_result$accession <- acc; return(empty_result) }

  len_line <- grep("^SQ   SEQUENCE", txt, value = TRUE)
  protein_len <- if (length(len_line) == 0) NA_integer_ else
    as.integer(sub("^SQ   SEQUENCE\\s+([0-9]+) AA.*", "\\1", len_line[1]))

  isoforms <- .uniprot_parse_isoforms(txt)

  displayed <- grep("IsoId=.*Sequence=Displayed", txt, value = TRUE)
  canonical_isoform <- if (length(displayed) == 0) NA_character_ else
    sub(".*IsoId=[A-Za-z0-9]+-([0-9]+);.*", "\\1", displayed[1])
  canon_idx <- match(canonical_isoform, isoforms$number)
  # Single-isoform gene (no ALTERNATIVE PRODUCTS block): the canonical
  # sequence IS the whole entry, referenced by the bare accession, no "-N".
  canonical_isoform_id   <- if (is.na(canon_idx)) acc else isoforms$isoid[canon_idx]
  canonical_isoform_name <- if (is.na(canon_idx)) NA_character_ else isoforms$name[canon_idx]

  list(accession = acc, protein_len = protein_len,
       ensembl = .uniprot_parse_dr(txt, "Ensembl", paste0(ens_tx_prefix, "[0-9]+")),
       ccds    = .uniprot_parse_dr(txt, "CCDS", "CCDS[0-9]+\\.[0-9]+"),
       isoforms = isoforms, canonical_isoform = canonical_isoform,
       canonical_isoform_id = canonical_isoform_id, canonical_isoform_name = canonical_isoform_name)
}

# present = database cross-references this id; matched = that reference points at the canonical isoform (also TRUE when untagged, since untagged means no isoforms to distinguish); isoform = the tag if any, else NA (so a non-canonical hit can still say which isoform it is).
.uniprot_tier_lookup <- function(ids, tag_map, canonical_isoform) {
  idx <- match(ids, tag_map$id)
  present <- !is.na(idx)
  isoform <- tag_map$isoform[idx]
  matched <- rep(NA, length(ids))
  matched[present & is.na(isoform)] <- TRUE
  has_iso <- present & !is.na(isoform)
  matched[has_iso] <- !is.na(canonical_isoform) & isoform[has_iso] == canonical_isoform
  list(present = present, matched = matched, isoform = isoform)
}

TIER_ENSEMBL <- "ensembl_xref"
TIER_CCDS    <- "ccds_xref"
TIER_LENGTH  <- "length_match"
TIER_ORDER   <- c(TIER_ENSEMBL, TIER_CCDS, TIER_LENGTH)

# One row per element of tx_ids (order preserved, versioned IDs matching cds_by_tx's names, NOT tx_id_stripped); tx_ccds_id must be in the same order. Full rationale: solutions/concepts/19-relocated-code-rationale.md, "utils/uniprot_canonical.R:200".
uniprot_canonical_match <- function(gene_symbol, tx_ids, cds_by_tx, tx_ccds_id = NULL, orgdb) {
  n <- length(tx_ids)
  na_result <- data.frame(uniprot_accession = rep(NA_character_, n),
                           uniprot_match = rep(NA, n),
                           uniprot_unique = rep(NA, n),
                           uniprot_match_count = rep(NA_integer_, n),
                           uniprot_source = rep(NA_character_, n),
                           uniprot_canonical_isoform_id = rep(NA_character_, n),
                           uniprot_canonical_isoform_name = rep(NA_character_, n),
                           uniprot_row_isoform_id = rep(NA_character_, n),
                           uniprot_row_isoform_name = rep(NA_character_, n))
  if (n == 0) return(na_result)
  if (is.null(tx_ccds_id)) tx_ccds_id <- rep(NA_character_, n)

  # Cache scoped by orgdb, not just gene_symbol -- the same case-sensitive
  # symbol string isn't guaranteed unique across species (mirrors
  # utils/protein_domains.R's PAX5-bug fix: mixing species under one cache
  # key served one species' cached UniProt entry for the other).
  cache_key <- paste0(orgdb, "|", gene_symbol)
  cache <- .uniprot_canonical_cache_load()
  entry <- cache[[cache_key]]
  if (is.null(entry)) {
    entry <- .uniprot_fetch_entry(gene_symbol, orgdb)
    cache[[cache_key]] <- entry
    tryCatch(.uniprot_canonical_cache_save(cache), error = function(e) NULL)
  }
  if (is.na(entry$accession) || is.na(entry$protein_len)) return(na_result)

  tx_stripped <- sub("\\.[0-9]+$", "", tx_ids)

  ens_hit  <- .uniprot_tier_lookup(tx_stripped, entry$ensembl, entry$canonical_isoform)
  ccds_hit <- .uniprot_tier_lookup(tx_ccds_id, entry$ccds, entry$canonical_isoform)

  tx_aa_len <- vapply(tx_ids, function(t) {
    if (!t %in% names(cds_by_tx)) return(NA_real_)
    sum(GenomicRanges::width(cds_by_tx[[t]])) / 3 - 1
  }, numeric(1))
  len_matched <- !is.na(tx_aa_len) & tx_aa_len == entry$protein_len

  tier    <- ifelse(ens_hit$present, TIER_ENSEMBL, ifelse(ccds_hit$present, TIER_CCDS, TIER_LENGTH))
  matched <- ifelse(ens_hit$present, ens_hit$matched, ifelse(ccds_hit$present, ccds_hit$matched, len_matched))

  if (!is.na(entry$canonical_isoform)) {
    # Multi-isoform gene: canonical_isoform is a known, specific isoform
    # number, so a weaker tier's coincidental length match must defer to a
    # stronger tier that has already identified which transcript(s) carry
    # that isoform -- otherwise an unrelated same-length transcript could
    # dilute or contradict a curator-confirmed answer.
    best_tier <- TIER_LENGTH
    for (t in TIER_ORDER) {
      if (any(matched & (tier == t), na.rm = TRUE)) { best_tier <- t; break }
    }
    final_matched <- matched & (tier == best_tier)
  } else {
    # Single-isoform gene (no ALTERNATIVE PRODUCTS block at all): every tier
    # that cross-references a transcript is equally valid evidence for the
    # one true sequence -- confirmed on mouse Apoe, where UniProt's Ensembl
    # DR lines list 4 transcripts and its CCDS DR line (shared by many more
    # GTF transcripts via the same CCDS ID) lists a 5th route to the same
    # single protein. No tier to defer to; union them all.
    final_matched <- matched
  }
  final_matched[is.na(final_matched)] <- FALSE
  n_match <- sum(final_matched)

  row_isoform <- ifelse(ens_hit$present, ens_hit$isoform, ifelse(ccds_hit$present, ccds_hit$isoform, NA_character_))
  iso_idx <- match(row_isoform, entry$isoforms$number)
  row_isoform_id   <- ifelse(is.na(iso_idx), NA_character_, entry$isoforms$isoid[iso_idx])
  row_isoform_name <- ifelse(is.na(iso_idx), NA_character_, entry$isoforms$name[iso_idx])

  data.frame(uniprot_accession = rep(entry$accession, n),
             uniprot_match = final_matched,
             uniprot_unique = rep(n_match == 1, n),
             uniprot_match_count = rep(n_match, n),
             uniprot_source = tier,
             uniprot_canonical_isoform_id = rep(entry$canonical_isoform_id, n),
             uniprot_canonical_isoform_name = rep(entry$canonical_isoform_name, n),
             uniprot_row_isoform_id = row_isoform_id,
             uniprot_row_isoform_name = row_isoform_name)
}
