# rationale: see solutions/concepts/19-relocated-code-rationale.md, "utils/transcript_meta.R:1"

# TX_META_LOOKUP_PATH is NOT defined here -- this file is only ever inlined
# (never sourced, see header comment above), and scripts/01_imports.R
# already defines it as a global in every assembled script (Phase C.0
# decision 6, GENOME$tx_meta_lookup_path) before this text is pasted in.

# GENCODE (vM39) biotypes confidently bucketed as non-protein-coding.
# Deliberately excludes a few biotypes whose coding status GENCODE itself
# marks as unresolved or mixed -- "TEC" literally stands for "To be
# Experimentally Confirmed", and "protein_coding_LoF"/"non_stop_decay" are
# coding-adjacent decay/loss-of-function calls, not a clean non-coding
# assertion. Those three fall through to "unknown" below rather than being
# silently folded into "non-coding".
TX_NONCODING_TYPES <- c(
  "lncRNA", "miRNA", "misc_RNA", "Mt_rRNA", "Mt_tRNA", "processed_transcript",
  "retained_intron", "ribozyme", "rRNA", "scaRNA", "scRNA", "snoRNA", "snRNA", "sRNA",
  "processed_pseudogene", "unprocessed_pseudogene", "unitary_pseudogene", "pseudogene",
  "transcribed_processed_pseudogene", "transcribed_unitary_pseudogene",
  "transcribed_unprocessed_pseudogene", "translated_unprocessed_pseudogene",
  "IG_C_gene", "IG_C_pseudogene", "IG_D_gene", "IG_D_pseudogene", "IG_J_gene",
  "IG_LV_gene", "IG_pseudogene", "IG_V_gene", "IG_V_pseudogene",
  "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_J_pseudogene", "TR_V_gene", "TR_V_pseudogene"
)

# coding_status is never NA/blank -- a missing or unbucketed transcript becomes "unknown" (a gap to flag, not a silent guess), one row per `tx_id_stripped` element in order.
transcript_coding_meta <- function(tx_id_stripped, lookup_path = TX_META_LOOKUP_PATH) {
  meta <- readRDS(lookup_path)
  m <- meta[tx_id_stripped, , drop = FALSE]

  coding_status <- ifelse(is.na(m$transcript_type), "unknown",
    ifelse(m$transcript_type == "protein_coding", "coding",
    ifelse(m$transcript_type == "nonsense_mediated_decay", "NMD",
    ifelse(m$transcript_type == "protein_coding_CDS_not_defined", "no CDS",
    ifelse(m$transcript_type %in% TX_NONCODING_TYPES, "non-coding", "unknown")))))

  is_canonical <- !is.na(m$tag) & (grepl("Ensembl_canonical", m$tag) | grepl("appris_principal", m$tag))

  data.frame(coding_status = coding_status, is_canonical = is_canonical, row.names = NULL)
}
