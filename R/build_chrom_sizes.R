# Chromosome name -> length for the genome, from UCSC's genome registry (GenomeInfoDb), by
# genome_label: a GTF carries no chromosome lengths. Deliberately no fallback to a hardcoded
# species' lengths: a wrong genome_label errors here rather than silently mislabelling
# lengths. ensure_genome_ready() builds it once per genome and the app reads the cached RDS.

build_chrom_sizes <- function(genome_label, out_path) {
  si <- tryCatch(
    GenomeInfoDb::Seqinfo(genome = genome_label),
    error = function(e) stop(
      "Couldn't fetch chromosome lengths for genome_label '", genome_label,
      "' from UCSC. genome_label must be the exact UCSC genome-build id ",
      "(e.g. \"hg38\", \"mm10\", \"mm39\") for this lookup to work. ",
      "Underlying error: ", conditionMessage(e), call. = FALSE
    )
  )
  lengths <- setNames(GenomeInfoDb::seqlengths(si), GenomeInfoDb::seqnames(si))
  lengths <- lengths[!is.na(lengths) & base::grepl("^chr", names(lengths))]
  if (length(lengths) == 0)
    stop("UCSC returned no chr-prefixed chromosome lengths for genome_label '", genome_label, "'.", call. = FALSE)

  saveRDS(lengths, out_path)
  invisible(lengths)
}
