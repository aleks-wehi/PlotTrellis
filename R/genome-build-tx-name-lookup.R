# Build: transcript_id -> transcript_name lookup, extracted from a GENCODE-
# style GTF's "transcript" feature lines (one row per transcript, not per
# exon). ensure_genome_ready() builds it the first time a genome is used;
# renders then read the cached RDS instead of re-parsing the gzipped GTF
# per gene.

# `lines`: the GTF's lines, or any subset holding all its "transcript" lines
# (.pt_genome_files() reads the file once for all three lookups).
build_tx_name_lookup <- function(
    gtf_path, out_path,
    lines = .pt_gtf_feature_lines(gtf_path, "transcript")) {
  # Prefilter to "transcript" feature lines in base R (portable -- no zgrep
  # binary dependency, which isn't guaranteed to exist, e.g. on Windows).
  # Slower than a zgrep prefilter but this only ever runs once per genome
  # (ensure_genome_ready() caches the result), never per render.
  lines <- lines[grepl("\ttranscript\t", lines, fixed = TRUE)]
  fields <- strsplit(lines, "\t")
  attrs <- vapply(fields, `[`, character(1), 9)

  tx_id   <- sub('.*transcript_id "([^"]+)".*', "\\1", attrs)
  tx_name <- sub('.*transcript_name "([^"]+)".*', "\\1", attrs)

  lookup <- setNames(tx_name, sub("\\.[0-9]+$", "", tx_id))
  .pt_save_rds(lookup, out_path)
  invisible(lookup)
}

# Build: gene_id (version stripped) -> gene_name, from a GTF's "gene" lines; a
# gene without a gene_name (an Ensembl GTF leaves it out) is named by its
# gene_id. The names of the genes the OrgDb has no symbol for; its "gtf"
# attribute is the GTF's file name, for messages. Built by
# .pt_genome_files() beside the TxDb (.pt_gene_names_file()).
# `lines` as for build_tx_name_lookup(), holding all the "gene" lines.
build_gene_name_lookup <- function(
    gtf_path, out_path, lines = .pt_gtf_feature_lines(gtf_path, "gene")) {
  lines <- lines[grepl("\tgene\t", lines, fixed = TRUE)]
  attrs <- vapply(strsplit(lines, "\t"), `[`, character(1), 9)
  id <- sub("\\.[0-9]+$", "", sub('.*gene_id "([^"]+)".*', "\\1", attrs))
  name <- ifelse(grepl('gene_name "', attrs, fixed = TRUE),
                 sub('.*gene_name "([^"]+)".*', "\\1", attrs), id)
  lookup <- structure(setNames(name, id), gtf = basename(gtf_path))
  .pt_save_rds(lookup, out_path)
  invisible(lookup)
}

# The lines of a GTF (gzipped or not) whose feature (third) column is one of
# `features`, in file order. Read `chunk` lines at a time, so only the kept
# lines are ever all in memory: a whole-genome GTF is about 4 GB as R strings,
# its transcript and gene lines about 0.2 GB.
.pt_gtf_feature_lines <- function(gtf_path, features, chunk = 200000L) {
  con <- gzfile(gtf_path, "rt")
  on.exit(close(con))
  patterns <- paste0("\t", features, "\t")
  kept <- list()
  repeat {
    lines <- readLines(con, n = chunk)
    if (length(lines) == 0) break
    hit <- Reduce(`|`, lapply(patterns, grepl, x = lines, fixed = TRUE))
    kept[[length(kept) + 1]] <- lines[hit]
  }
  as.character(unlist(kept))
}

# saveRDS() to a temporary file in the same folder, then renamed into place, so
# a build stopped part-way (the R process building it killed) never leaves a
# half-written cache file at `path`.
.pt_save_rds <- function(object, path) {
  tmp <- tempfile(paste0(basename(path), "-"), tmpdir = dirname(path))
  on.exit(unlink(tmp))
  saveRDS(object, tmp)
  if (!file.rename(tmp, path))
    stop("PlotTrellis could not save ", path, ".", call. = FALSE)
}
