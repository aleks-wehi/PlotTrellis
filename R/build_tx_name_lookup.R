# Build: transcript_id -> transcript_name lookup, extracted from a GENCODE-
# style GTF's "transcript" feature lines (one row per transcript, not per
# exon). ensure_genome_ready() builds it the first time a genome is used;
# renders then read the cached RDS instead of re-parsing the gzipped GTF
# per gene.

build_tx_name_lookup <- function(gtf_path, out_path) {
  # Prefilter to "transcript" feature lines in base R (portable -- no zgrep
  # binary dependency, which isn't guaranteed to exist, e.g. on Windows).
  # Slower than a zgrep prefilter (~90s vs ~13s on the real mm39 GTF,
  # measured) but this only ever runs once per genome
  # (ensure_genome_ready() caches the result), never per render.
  lines <- readLines(gzfile(gtf_path))
  lines <- lines[grepl("\ttranscript\t", lines, fixed = TRUE)]
  fields <- strsplit(lines, "\t")
  attrs <- vapply(fields, `[`, character(1), 9)

  tx_id   <- sub('.*transcript_id "([^"]+)".*', "\\1", attrs)
  tx_name <- sub('.*transcript_name "([^"]+)".*', "\\1", attrs)

  lookup <- setNames(tx_name, sub("\\.[0-9]+$", "", tx_id))
  saveRDS(lookup, out_path)
  invisible(lookup)
}
