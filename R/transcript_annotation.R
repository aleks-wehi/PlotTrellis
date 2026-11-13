# Per-transcript annotation for the Gene info table (not part of any figure): coding status,
# and each canonical-transcript system's call kept separate (Ensembl canonical, APPRIS, CCDS,
# GENCODE basic), since they routinely disagree and one combined flag would hide that.

# The OrgDb and the two lookup paths are explicit arguments, so the table works
# outside the composer too (the composer's own values are not package globals).

# appris_principal_N/appris_alternative_N tags are the one GTF tag not subsumed by a plain flag -- pull out the tier so "principal_2"/"alternative_2" stay distinguishable, not collapsed into one boolean.
.appris_tier_from_tag <- function(tag) {
  # NB: regmatches(x, regexpr(...)) silently drops non-matching elements
  # instead of returning NA for them, which would desync this from `tag`'s
  # row order -- grepl + sub (both always length(tag)-preserving) avoids that.
  has_appris <- !is.na(tag) & grepl("appris_(principal|alternative)_[0-9]+", tag)
  tier <- rep(NA_character_, length(tag))
  tier[has_appris] <- sub(".*appris_(principal_[0-9]+|alternative_[0-9]+).*", "\\1", tag[has_appris])
  tier
}

# One row per annotated transcript of GENE_SYMBOL (NULL if it doesn't resolve to an Ensembl gene ID / has no transcripts).
transcript_annotations <- function(gene_symbol, txdb, orgdb_pkg, tx_name_lookup_path, tx_meta_lookup_path) {
  ens_res <- suppressMessages(AnnotationDbi::select(get(orgdb_pkg, envir = asNamespace(orgdb_pkg)),
    keys = gene_symbol, keytype = "SYMBOL", columns = "ENSEMBL"))
  ens_id <- ens_res$ENSEMBL[!is.na(ens_res$ENSEMBL)][1]
  if (length(ens_id) == 0 || is.na(ens_id)) return(NULL)

  q <- suppressMessages(AnnotationDbi::select(txdb, keys = ens_id, keytype = "GENEID", columns = "TXNAME"))
  tx_ids <- unique(q$TXNAME)
  if (length(tx_ids) == 0) return(NULL)
  tx_stripped <- sub("\\.[0-9]+$", "", tx_ids)

  name_lookup <- readRDS(tx_name_lookup_path)
  meta_lookup <- readRDS(tx_meta_lookup_path)

  tx_name <- unname(name_lookup[tx_stripped])
  tx_name[is.na(tx_name)] <- tx_ids[is.na(tx_name)]

  meta <- meta_lookup[tx_stripped, , drop = FALSE]

  out <- data.frame(
    tx_id           = tx_ids,
    transcript_name = tx_name,
    transcript_type = meta$transcript_type,
    tag             = meta$tag,
    ccds_id         = meta$ccds_id,
    stringsAsFactors = FALSE
  )
  # Three buckets, not a coding/non-coding binary: "protein_coding_CDS_not_
  # defined" (e.g. Nfatc1-207) is a real, distinct GENCODE category from a
  # genuinely non-coding transcript type -- collapsing it into "non-coding"
  # would misstate what the annotation actually says.
  out$coding_status <- ifelse(is.na(out$transcript_type), "unknown",
    ifelse(out$transcript_type == "nonsense_mediated_decay", "NMD",
    ifelse(out$transcript_type == "protein_coding", "coding",
    ifelse(out$transcript_type == "protein_coding_CDS_not_defined", "no CDS",
    "non-coding"))))

  out$ensembl_canonical <- !is.na(out$tag) & grepl("Ensembl_canonical", out$tag)
  out$gencode_basic     <- !is.na(out$tag) & grepl("(^|,)basic(,|$)", out$tag)
  out$appris_tier       <- .appris_tier_from_tag(out$tag)

  out[order(out$transcript_name), ]
}

# `orgdb` is required, not defaulted, so the UniProt link-out searches the loaded genome's species.
uniprot_search_url <- function(gene_symbol, orgdb, reviewed_only = FALSE) {
  if (!orgdb %in% names(ORGDB_TAXON_ID))
    stop("No taxon ID mapping for orgdb '", orgdb, "' in ORGDB_TAXON_ID (needed to link the gene to its UniProt search) -- add one.", call. = FALSE)
  query <- sprintf("(gene:%s) AND (organism_id:%d)", gene_symbol, ORGDB_TAXON_ID[[orgdb]])
  if (reviewed_only) query <- paste(query, "AND (reviewed:true)")
  paste0("https://www.uniprot.org/uniprotkb?query=", utils::URLencode(query, reserved = TRUE))
}
