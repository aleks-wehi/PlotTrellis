# One transcript filter shared by Transcript expression and the transcript-structure tracks, so a
# transcript is never named in one and folded into "Other" in the other. A transcript can be
# named only where its mean TPM in a condition reaches min_display_tpm: an absolute floor, not
# a statistical test, so the same rule works for any design.

# The transcript functions below take the annotation as one argument, `anno`, from
# .pt_annotation() (boilerplate.R): the open TxDb, the OrgDb package name and the
# pre-built transcript_id -> transcript_name lookup (build_tx_name_lookup.R), one row
# per transcript, which spares gtf_transcript_structure() a ~13s re-parse of the
# gzipped GTF on every render. The minimum TPM for a transcript to be shown by name is
# the project setting min_display_tpm (default 1), an argument of
# gene_display_transcripts().

# Colours by transcript name, alphabetical, so the structure track and Transcript expression agree;
# TX_OTHER_COLOUR (not from the palette) marks the "Other" bucket.
TX_PALETTE      <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948")
TX_OTHER_COLOUR <- "grey70"

transcript_colours <- function(transcript_names, palette = TX_PALETTE) {
  levels <- sort(unique(transcript_names))
  setNames(palette[((seq_along(levels) - 1) %% length(palette)) + 1], levels)
}

strip_version <- function(x) sub("\\.[0-9]+$", "", x)

# Resolves `chrom` against the file's own observed values (`vocab`), trying only the two-way "chr"-prefix toggle (a file may say "18" where the query says "chr18"); falls back to the original chrom unresolved, rather than erroring, if neither form matches.
resolve_chrom_against_vocab <- function(chrom, vocab) {
  if (chrom %in% vocab) return(chrom)
  toggled <- if (grepl("^chr", chrom, ignore.case = TRUE)) {
    sub("^chr", "", chrom, ignore.case = TRUE)
  } else {
    paste0("chr", chrom)
  }
  if (toggled %in% vocab) toggled else chrom
}

# TRUE when `first_line` is a HiCCUPS loop file's column-name header ("#chr1 x1 x2 chr2 y1 y2
# name score ... observed ..."), as in merged_loops.bedpe and postprocessed_pixels_<res>.bedpe.
# A plain BEDPE with a commented "#chr1 x1 x2 chr2 y1 y2 ..." header has no `observed` column,
# so it stays a BEDPE.
is_hiccups_header <- function(first_line) {
  if (length(first_line) == 0) return(FALSE)
  f <- strsplit(trimws(first_line), "[ \t]+")[[1]]
  length(f) >= 6 && identical(f[1:6], c("#chr1", "x1", "x2", "chr2", "y1", "y2")) && "observed" %in% f
}

# Reads a HiCCUPS loop file (see is_hiccups_header()) as the BEDPE its loops convert to: the
# six anchor columns as they are, nSupport = the `observed` column (contacts at the loop pixel),
# found by name in the header line. Other comment lines and rows whose start is not a whole
# number are skipped; fields split on runs of spaces or tabs. Values take the types read.delim
# would give the converted file (integer where they fit, else double; a missing name is "").
read_hiccups_loops <- function(path) {
  lines <- readLines(path, warn = FALSE)
  obs <- match("observed", sub("^#", "", strsplit(trimws(lines[1]), "[ \t]+")[[1]]))
  body <- lines[-1]
  rows <- strsplit(trimws(body[!grepl("^#", body)]), "[ \t]+")
  field <- function(k) vapply(rows, `[`, character(1), k)
  rows <- rows[grepl("^[0-9]+$", field(2))]
  num <- function(k) utils::type.convert(field(k), as.is = TRUE)
  chr <- function(k) { x <- field(k); x[is.na(x)] <- ""; x }  # a short row's missing name, as read.delim fills it
  data.frame(chrom1 = chr(1), start1 = num(2), end1 = num(3),
             chrom2 = chr(4), start2 = num(5), end2 = num(6),
             nSupport = num(obs), stringsAsFactors = FALSE)
}

# Keeps a BEDPE table's cis pairs on `chrom` whose arc reaches the window, as the gap between
# the two anchors (intron_start/intron_end), not either anchor's own span.
bedpe_arcs_in_window <- function(j, chrom, start, end) {
  empty <- data.frame(intron_start = integer(0), intron_end = integer(0), nSupport = integer(0))
  j <- j[j$chrom1 == j$chrom2, ]  # Rule 1: drop trans pairs
  if (nrow(j) == 0) return(empty)
  chrom <- resolve_chrom_against_vocab(chrom, unique(j$chrom1))
  j <- j[j$chrom1 == chrom, ]
  if (nrow(j) == 0) return(empty)
  first_is_left <- j$start1 <= j$start2
  j$intron_start <- ifelse(first_is_left, j$end1, j$end2)
  j$intron_end   <- ifelse(first_is_left, j$start2, j$start1)
  # Keep an arc plotPairsArches draws (clip = FALSE): its left end in [start, end) or its right end in
  # (start, end). An arc running past a window edge so reaches the arc_clip_partial setting
  # instead of being dropped here; one spanning the whole window, or ending on an edge, draws nothing
  # and must not count toward the track's scale.
  j <- j[(j$intron_start >= start & j$intron_start < end) | (j$intron_end > start & j$intron_end < end), ]
  if (nrow(j) == 0) return(empty)
  j[, c("intron_start", "intron_end", "nSupport")]
}

# Dispatches on format: a HiCCUPS header line is a HiCCUPS loop file, otherwise exactly 12
# columns is BED12, 6 with text (a strand) in column 6 is a junction BED6, and 6-10 is BEDPE; intron_start/intron_end is the gap between the two anchors, not either anchor's own span.
read_junction_file <- function(path, chrom, start, end) {
  first_line <- readLines(path, n = 1)
  if (is_hiccups_header(first_line)) return(bedpe_arcs_in_window(read_hiccups_loops(path), chrom, start, end))
  n_cols <- length(strsplit(first_line, "\t")[[1]])
  empty <- data.frame(intron_start = integer(0), intron_end = integer(0), nSupport = integer(0))

  if (n_cols == 12) {
    j <- read.delim(path, header = FALSE, comment.char = "#",
                     col.names = c("chrom", "chromStart", "chromEnd", "name", "nSupport", "strand",
                                   "thickStart", "thickEnd", "color", "nBlocks", "blockSizes", "blockStarts"))
    chrom <- resolve_chrom_against_vocab(chrom, unique(j$chrom))
    j <- j[j$chrom == chrom & j$chromStart < end & j$chromEnd > start & j$nBlocks == 2, ]
    if (nrow(j) == 0) return(empty)
    sizes  <- strsplit(j$blockSizes, ",")
    starts <- strsplit(j$blockStarts, ",")
    j$intron_start <- j$chromStart + as.integer(vapply(sizes, `[`, character(1), 1))   # 0-based
    j$intron_end   <- j$chromStart + as.integer(vapply(starts, `[`, character(1), 2))  # 0-based
    j[, c("intron_start", "intron_end", "nSupport")]
  } else if (n_cols >= 6 && n_cols <= 10) {
    col_names <- c("chrom1", "start1", "end1", "chrom2", "start2", "end2",
                    "name", "nSupport", "strand1", "strand2")[seq_len(n_cols)]
    j <- read.delim(path, header = FALSE, comment.char = "#", col.names = col_names)
    # Junction BED6 (chrom, start, end, name or ".", count, strand): column 4 is a name, where a
    # 6-column BEDPE has its second chromosome. Told apart by column 6, text (a strand) in a BED6
    # and a number (end2) in a BEDPE, since a name can look like a chromosome ("chr1:3001-7000", "1").
    if (n_cols == 6 && is.character(j[[6]])) {
      bad <- setdiff(trimws(j[[6]]), c("+", "-", "*", "."))
      if (length(bad)) stop("6-column junction/arc file ", path, ": column 6 must be a strand (+, -, * or .) ",
                            "for a junction BED6, or a number (end2) for a BEDPE; found '", bad[1], "'", call. = FALSE)
      # Its interval is the intron: drawn as a pair whose anchors are the intron's two ends.
      j <- data.frame(chrom1 = j[[1]], start1 = j[[2]], end1 = j[[2]],
                      chrom2 = j[[1]], start2 = j[[3]], end2 = j[[3]], nSupport = j[[5]])
    }
    if (!"nSupport" %in% names(j)) j$nSupport <- 1L  # no score column: every pair counts once
    bedpe_arcs_in_window(j, chrom, start, end)
  } else {
    stop("Unrecognised junction/arc file (", n_cols, " columns): ", path,
         " -- expected 12 columns (BED12), 6 (junction BED6), 6-10 (BEDPE) or a HiCCUPS loop file", call. = FALSE)
  }
}

# Sums nSupport across replicate junction files for the SAME condition BEFORE thresholding (not filter-then-union), matching ggsashimi's multi-BAM combining convention; used by load_junctions() for the drawn arcs, independent of the transcript-display filtering below.
pool_junctions <- function(samples, chrom, start, end, min_junction_support) {
  raw_list <- lapply(samples, function(s) read_junction_file(s$junction_path, chrom, start, end))
  raw <- do.call(rbind, raw_list)
  if (nrow(raw) == 0) return(data.frame(intron_start = integer(0), intron_end = integer(0), nSupport = integer(0)))

  pooled <- aggregate(nSupport ~ intron_start + intron_end, data = raw, FUN = sum)
  pooled[pooled$nSupport >= min_junction_support, ]
}

# Returns every annotated transcript for GENE, unfiltered (condition-independent; see gene_display_transcripts() below for which become named); sourced from anno's TxDb (~0.3s indexed select vs. ~13s GTF zgrep), transcript_name from anno's id->name lookup (falls back to transcript_id).
gtf_transcript_structure <- function(gene, anno) {
  ens_res <- suppressMessages(AnnotationDbi::select(get(anno$orgdb_pkg, envir = asNamespace(anno$orgdb_pkg)),
    keys = gene, keytype = "SYMBOL", columns = "ENSEMBL"))
  ens_id <- ens_res$ENSEMBL[!is.na(ens_res$ENSEMBL)][1]
  if (length(ens_id) == 0 || is.na(ens_id)) {
    return(data.frame(transcript_id = character(0), transcript_name = character(0),
      tx_id_stripped = character(0)))
  }

  q <- suppressMessages(AnnotationDbi::select(anno$txdb, keys = ens_id, keytype = "GENEID", columns = "TXNAME"))
  tx_ids <- unique(q$TXNAME)
  tx_id_stripped <- strip_version(tx_ids)

  name_lookup <- readRDS(anno$tx_name_lookup_path)
  tx_name <- unname(name_lookup[tx_id_stripped])
  tx_name[is.na(tx_name)] <- tx_ids[is.na(tx_name)]

  tx <- data.frame(transcript_id = tx_ids, transcript_name = tx_name, stringsAsFactors = FALSE)
  tx <- tx[order(tx$transcript_name), ]
  tx$tx_id_stripped <- strip_version(tx$transcript_id)
  tx
}

# GENCODE's numeric transcript suffix ("Nfatc1-201" -> 201L), used ONLY as an ordinal fallback ranking when a project has no quant.sf data at all (see gene_display_transcripts() below) -- a weak signal (naming convention, not a real abundance ranking), but the only one available then; unparseable names return NA and sort last via na.last=TRUE at the call site.
tx_ordinal_suffix <- function(transcript_name) {
  suffix <- sub(".*-([0-9]+)$", "\\1", transcript_name)
  suffix[suffix == transcript_name] <- NA_character_
  suppressWarnings(as.integer(suffix))
}

# `samples` is a list of list(sample_id=, quant_path=) for one or more replicates of the SAME condition/gene; TPM is the mean across samples, `proportion` is its share of ALL annotated transcripts' TPM (not just the retained set), so callers can derive a meaningful "Other" residual.
filtered_transcripts <- function(gene, samples, anno) {
  tx <- gtf_transcript_structure(gene, anno)

  not_numeric <- character(0)
  quant <- do.call(rbind, lapply(samples, function(s) {
    .pt_check_file(s$quant_path)  # a missing quant file fails the track, naming the file
    d <- read.delim(s$quant_path)
    short <- file.path(basename(dirname(s$quant_path)), basename(s$quant_path))
    if (is.null(d$TPM)) stop("no TPM column in ", short, " (is it Salmon's quant.sf?)", call. = FALSE)
    # One stray word reads the whole column as text: counted as 0 TPM and flagged for
    # the caller rather than stopping, so other samples' tracks still draw.
    if (nrow(d) && !is.numeric(d$TPM)) {
      not_numeric <<- c(not_numeric, short)
      d$TPM <- 0
    }
    d$tx_stripped <- strip_version(d$Name)
    d$sample_id   <- s$sample_id
    d[d$tx_stripped %in% tx$tx_id_stripped, c("tx_stripped", "TPM", "sample_id")]
  }))
  max_tpm  <- tapply(quant$TPM, quant$tx_stripped, max)
  mean_tpm <- tapply(quant$TPM, quant$tx_stripped, mean)

  tx$max_tpm  <- as.numeric(max_tpm[tx$tx_id_stripped])
  tx$max_tpm[is.na(tx$max_tpm)] <- 0
  tx$mean_tpm <- as.numeric(mean_tpm[tx$tx_id_stripped])
  tx$mean_tpm[is.na(tx$mean_tpm)] <- 0

  # Total gene TPM = sum of mean_tpm across ALL annotated transcripts (kept
  # and not), so proportion (below) and any caller's "kept segments + Other"
  # always sum to exactly this value -- one aggregation order, not two that
  # could quietly disagree.
  total_tpm <- sum(tx$mean_tpm)
  tx$proportion <- if (total_tpm > 0) tx$mean_tpm / total_tpm else 0

  out <- tx[, c("transcript_id", "transcript_name", "max_tpm", "mean_tpm", "proportion")]
  attr(out, "total_tpm") <- total_tpm
  # 0 = none of the gene's transcripts are in any quant file (vs. found, all at 0 TPM).
  attr(out, "n_quantified") <- length(unique(quant$tx_stripped))
  attr(out, "tpm_not_numeric") <- not_numeric  # quant files whose TPM column read as text
  out
}

# Orders retained transcripts ascending by proportion, folding the smallest into "Other" until other_budget or max_segments is exceeded, whichever is stricter (a budget on the Other share, not a flat per-transcript floor, so many transcripts that are each just above a floor but together negligible still fold); `other_baseline_tpm` is TPM already excluded by the retained-set filter and always folds into Other.
truncate_by_proportion <- function(kept_df, other_baseline_tpm, total_tpm,
                                    other_budget = 0.05, max_segments = length(TX_PALETTE)) {
  if (nrow(kept_df) == 0) {
    return(list(named = kept_df, dropped = kept_df, other_tpm = other_baseline_tpm))
  }
  ord <- kept_df[order(kept_df$proportion), ]
  cum_other   <- other_baseline_tpm
  n_remaining <- nrow(ord)
  drop        <- rep(FALSE, nrow(ord))
  for (i in seq_len(nrow(ord))) {
    over_cap       <- n_remaining > max_segments
    would_be_other <- cum_other + ord$mean_tpm[i]
    within_budget  <- if (total_tpm > 0) (would_be_other / total_tpm) <= other_budget else TRUE
    if (over_cap || within_budget) {
      drop[i]     <- TRUE
      cum_other   <- would_be_other
      n_remaining <- n_remaining - 1
    } else break
  }
  named <- ord[!drop, ]
  list(named = named[order(-named$proportion), ], dropped = ord[drop, ], other_tpm = cum_other)
}

# Reconciles several samples' independently-truncated survivor sets into ONE shared named set (deciding survivors per sample first preserves isoform-switching signal a pooled ranking would wash out); the union, only if still over max_segments, is re-capped by each transcript's best proportion across samples.
join_truncated_transcripts <- function(named_list, max_segments = length(TX_PALETTE)) {
  all_named <- do.call(rbind, named_list)
  if (is.null(all_named) || nrow(all_named) == 0) return(character(0))
  best_prop   <- tapply(all_named$proportion, all_named$transcript_name, max)
  union_names <- names(best_prop)
  if (length(union_names) > max_segments) {
    union_names <- union_names[order(-best_prop)][seq_len(max_segments)]
  }
  sort(union_names)
}

.gene_display_transcripts_cache_load <- function(path) {
  if (file.exists(path)) readRDS(path) else list()
}

.gene_display_transcripts_cache_save <- function(cache, path) .pt_cache_save(cache, path)

# Gene-level (not per-run) transcript filtering, pooled over EVERY RNA-seq condition so the same gene always resolves to the same named set (falls back to tx_ordinal_suffix() ordering with no quant.sf at all); disk-cached per gene in cache_path, keyed on gene/GTF/other_budget/max_segments/samples/data_dir/min_display_tpm so different projects never share an entry.
gene_display_transcripts <- function(gene, anno, samples, data_dir, min_display_tpm = 1, cache_path,
                                      other_budget = 0.05, max_segments = length(TX_PALETTE)) {
  # Cache key derived from the samples LIST content, not a file path --
  # deparse() gives a deterministic string for a given samples list, which
  # is all a cache key needs (no cryptographic strength required). Two
  # different projects (or a project whose samples[] changes via the
  # wizard) querying the same gene symbol can't share a cache entry. gtf_path
  # is included so a project that switches GTF doesn't get the other GTF's
  # transcript structure for an already-cached gene, and min_display_tpm so
  # a changed min TPM setting applies to genes already in the cache.
  cache_key <- sprintf("%s|%s|%.6f|%d|%s|%s|%.10g", gene, anno$gtf_path, other_budget, max_segments,
                        paste(deparse(samples), collapse = ""), data_dir, min_display_tpm)
  cache <- .gene_display_transcripts_cache_load(cache_path)
  if (!is.null(cache[[cache_key]])) return(cache[[cache_key]])

  tx <- gtf_transcript_structure(gene, anno)

  pool <- load_samples_from_settings(samples, data_dir)
  pool <- Filter(function(s) !is.na(sample_data_path(s, "quantdir")), pool)

  corrupt <- FALSE
  named <- if (length(pool) == 0) {
    ord <- order(tx_ordinal_suffix(tx$transcript_name), na.last = TRUE)
    head(tx$transcript_name[ord], max_segments)
  } else {
    per_condition <- lapply(pool, function(row) {
      samples <- list(list(sample_id = row$display, quant_path = sample_data_path(row, "quantdir")))
      ft   <- filtered_transcripts(gene, samples, anno)
      if (length(attr(ft, "tpm_not_numeric"))) corrupt <<- TRUE
      kept <- ft[ft$mean_tpm >= min_display_tpm, ]
      truncate_by_proportion(kept, attr(ft, "total_tpm") - sum(kept$mean_tpm), attr(ft, "total_tpm"),
                              other_budget, max_segments)$named
    })
    join_truncated_transcripts(per_condition, max_segments)
  }

  result <- list(named = named, tx = tx)
  if (corrupt) return(result)  # not cached: the key can't tell when the file is fixed
  cache[[cache_key]] <- result
  .gene_display_transcripts_cache_save(cache, cache_path)
  result
}
