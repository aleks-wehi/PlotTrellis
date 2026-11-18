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

# Build: gene_id (version stripped) -> gene_name, from a GTF's "gene" lines; a
# gene without a gene_name (an Ensembl GTF leaves it out) is named by its
# gene_id. The names of the genes the OrgDb has no symbol for; its "gtf"
# attribute is the GTF's file name, for messages. Built by
# .pt_genome_files() beside the TxDb (.pt_gene_names_file()).
build_gene_name_lookup <- function(gtf_path, out_path) {
  lines <- readLines(gzfile(gtf_path))
  lines <- lines[grepl("\tgene\t", lines, fixed = TRUE)]
  attrs <- vapply(strsplit(lines, "\t"), `[`, character(1), 9)
  id <- sub("\\.[0-9]+$", "", sub('.*gene_id "([^"]+)".*', "\\1", attrs))
  name <- ifelse(grepl('gene_name "', attrs, fixed = TRUE),
                 sub('.*gene_name "([^"]+)".*', "\\1", attrs), id)
  lookup <- structure(setNames(name, id), gtf = basename(gtf_path))
  saveRDS(lookup, out_path)
  invisible(lookup)
}
