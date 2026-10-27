# Build: transcript_id -> transcript_name lookup, extracted from a GENCODE-
# style GTF's "transcript" feature lines (one row per transcript, not per
# exon -- far fewer lines than utils/transcript_filter.R's old exon-level
# zgrep). The app calls build_tx_name_lookup() directly (utils/genomes.R's
# ensure_genome_ready()) to build a new genome's cache the first time it's
# used; the app then reads the cached RDS at render time instead of
# re-parsing the gzipped GTF per gene.

build_tx_name_lookup <- function(gtf_path, out_path) {
  # Prefilter to "transcript" feature lines in base R (portable -- no zgrep
  # binary dependency, which isn't guaranteed to exist, e.g. on Windows).
  # Slower than the zgrep prefilter this replaces (~90s vs ~13s on the real
  # mm39 GTF, empirically measured) but this only ever runs once per genome
  # (ensure_genome_ready() caches the result), never per render.
  lines <- readLines(gzfile(gtf_path))
  lines <- lines[grepl("\ttranscript\t", lines, fixed = TRUE)]
  fields <- strsplit(lines, "\t")
  attrs <- vapply(fields, `[`, character(1), 9)

  tx_id   <- sub('.*transcript_id "([^"]+)".*', "\\1", attrs)
  tx_name <- sub('.*transcript_name "([^"]+)".*', "\\1", attrs)

  lookup <- setNames(tx_name, sub("\\.[0-9]+$", "", tx_id))
  saveRDS(lookup, out_path)
  message(sprintf("Wrote %d transcript_id -> transcript_name entries to %s", length(lookup), out_path))
  invisible(lookup)
}
