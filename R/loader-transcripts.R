load_transcripts <- function(entry, region_params, settings) {
  gene <- region_params$gene
  if (.gene_missing(gene)) {
    entry$data <- list(
      tx_all = data.frame(transcript_id = character(0), transcript_name = character(0),
                           tx_id_stripped = character(0), coding_status = character(0), is_canonical = logical(0)),
      named = character(0), n_annotated = 0L
    )
    return(entry)
  }
  .set_transcript_filter_globals(settings)
  tx_all <- gtf_transcript_structure(gene, settings$gtf_path)
  meta <- transcript_coding_meta(tx_all$tx_id_stripped, lookup_path = settings$tx_meta_lookup_path)
  tx_all$coding_status <- meta$coding_status
  tx_all$is_canonical  <- meta$is_canonical

  needs_named <- entry$viz %in% c("transcripts", "transcripts_annotated", "tpm_usage")
  named <- if (needs_named) {
    gene_display_transcripts(gene, settings$gtf_path, settings$samples, settings$data_dir)$named
  } else {
    tx_all$transcript_name
  }

  entry$data <- list(tx_all = tx_all, named = named, n_annotated = nrow(tx_all))
  entry
}
