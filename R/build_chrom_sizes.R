# rationale: see solutions/concepts/19-relocated-code-rationale.md, "utils/build_chrom_sizes.R:1"

build_chrom_sizes <- function(genome_label, out_path) {
  si <- tryCatch(
    GenomeInfoDb::Seqinfo(genome = genome_label),
    error = function(e) stop(
      "Couldn't fetch chromosome lengths for genome_label '", genome_label,
      "' from UCSC. genome_label must be the exact UCSC genome-build id ",
      "(e.g. \"hg38\", \"mm10\", \"mm39\") for this lookup to work. ",
      "Underlying error: ", conditionMessage(e)
    )
  )
  lengths <- setNames(GenomeInfoDb::seqlengths(si), GenomeInfoDb::seqnames(si))
  lengths <- lengths[!is.na(lengths) & base::grepl("^chr", names(lengths))]
  if (length(lengths) == 0)
    stop("UCSC returned no chr-prefixed chromosome lengths for genome_label '", genome_label, "'.")

  saveRDS(lengths, out_path)
  message(sprintf("Wrote %d chromosome lengths to %s", length(lengths), out_path))
  invisible(lengths)
}
