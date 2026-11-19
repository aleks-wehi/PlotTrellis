# transcript_id -> transcript_type, tag and ccdsid, from a GENCODE- or
# Ensembl-style GTF's "transcript" lines. `tag` can appear several times per
# transcript ("Ensembl_canonical", "appris_principal_1", "CCDS"), so every one
# is kept, comma-joined; `ccdsid` (e.g. "CCDS37873.1") names the CCDS record the
# "CCDS" tag only says exists. Built by ensure_genome_ready() the first time a
# genome is used.

# `lines` as for build_tx_name_lookup().
build_tx_meta_lookup <- function(
    gtf_path, out_path,
    lines = .pt_gtf_feature_lines(gtf_path, "transcript")) {
  # Portable base-R prefilter (no zgrep dependency) -- see
  # build_tx_name_lookup()'s comment for the same tradeoff (one-time per
  # genome, not per render).
  lines <- lines[grepl("\ttranscript\t", lines, fixed = TRUE)]
  fields <- strsplit(lines, "\t")
  attrs <- vapply(fields, `[`, character(1), 9)

  tx_id <- sub('.*transcript_id "([^"]+)".*', "\\1", attrs)
  tx_id_stripped <- sub("\\.[0-9]+$", "", tx_id)

  # GENCODE uses "transcript_type", Ensembl GTFs use "transcript_biotype" for
  # the same field; sub() returns its input unchanged (not NA) on a non-match,
  # so matching only one name silently corrupted every row into garbage on a
  # real Ensembl GRCm39 GTF (149547/149547 rows). Try
  # both names; genuinely unmatched rows fall through to NA, which
  # data-transcript-meta.R/data-transcript-annotation.R already treat as "unknown".
  transcript_type <- ifelse(grepl('transcript_type "', attrs),
    sub('.*transcript_type "([^"]+)".*', "\\1", attrs),
    ifelse(grepl('transcript_biotype "', attrs),
      sub('.*transcript_biotype "([^"]+)".*', "\\1", attrs),
      NA_character_))

  tag_matches <- regmatches(attrs, gregexpr('tag "([^"]+)"', attrs))
  tag <- vapply(tag_matches, function(x) {
    if (length(x) == 0) return("")
    paste(sub('tag "([^"]+)"', "\\1", x), collapse = ",")
  }, character(1))

  # GENCODE spells the attribute "ccdsid", Ensembl GTFs "ccds_id".
  ccds_id <- ifelse(grepl('ccdsid "', attrs),
    sub('.*ccdsid "([^"]+)".*', "\\1", attrs),
    ifelse(grepl('ccds_id "', attrs),
      sub('.*ccds_id "([^"]+)".*', "\\1", attrs),
      NA_character_))

  meta <- data.frame(transcript_type = transcript_type, tag = tag,
                     ccds_id = ccds_id,
                     stringsAsFactors = FALSE, row.names = tx_id_stripped)
  .pt_save_rds(meta, out_path)
  invisible(meta)
}
