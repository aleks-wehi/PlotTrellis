# rationale: see solutions/concepts/19-relocated-code-rationale.md, "utils/transcript_filter.R:1"

# rationale: see solutions/concepts/19-relocated-code-rationale.md, "utils/transcript_filter.R:28"

# MIN_DISPLAY_TPM is NOT set here -- same mechanism as MIN_JUNCTION_SUPPORT
# above: it's a project setting (min_display_tpm, default 1) assembled into
# MIN_DISPLAY_TPM by scripts/02_params.R, which is inlined BEFORE this file
# in every assembled script (utils/assemble.R: s2 then s3). Setting it here
# would silently clobber that value.

# Pre-built transcript_id -> transcript_name lookup (utils/build_tx_name_lookup.R),
# one row per transcript instead of one per exon -- lets gtf_transcript_structure()
# below skip re-parsing the 88MB gzipped GTF (a ~13s zgrep) on every render;
# rebuild only if the GTF annotation itself changes. NOT defined here --
# this file is only ever inlined (never sourced, see header comment above),
# and scripts/01_imports.R already defines TX_NAME_LOOKUP_PATH (and
# ORGDB_PKG, used below) as globals in every assembled script (Phase C.0
# decision 6, GENOME$tx_name_lookup_path) before this text is pasted in.

# rationale: see solutions/concepts/19-relocated-code-rationale.md, "utils/transcript_filter.R:61"
TX_PALETTE      <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948")
TX_OTHER_COLOUR <- "grey70"

transcript_colours <- function(transcript_names, palette = TX_PALETTE) {
  levels <- sort(unique(transcript_names))
  setNames(palette[((seq_along(levels) - 1) %% length(palette)) + 1], levels)
}

strip_version <- function(x) sub("\\.[0-9]+$", "", x)

# Resolves `chrom` against the file's own observed values (`vocab`), trying only the two-way "chr"-prefix toggle (ARCHIVE.md bug 25); falls back to the original chrom unresolved, rather than erroring, if neither form matches. See solutions/concepts/19-relocated-code-rationale.md, "utils/transcript_filter.R:89".
resolve_chrom_against_vocab <- function(chrom, vocab) {
  if (chrom %in% vocab) return(chrom)
  toggled <- if (grepl("^chr", chrom, ignore.case = TRUE)) {
    sub("^chr", "", chrom, ignore.case = TRUE)
  } else {
    paste0("chr", chrom)
  }
  if (toggled %in% vocab) toggled else chrom
}

# Dispatches on column count: exactly 12 columns is BED12, 6-10 is BEDPE (feature 25, solutions/concepts/18-bed12-bedpe-arc-format.md); intron_start/intron_end is the gap between the two anchors, not either anchor's own span. See solutions/concepts/19-relocated-code-rationale.md, "utils/transcript_filter.R:89".
read_junction_file <- function(path, chrom, start, end) {
  first_line <- readLines(path, n = 1)
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
    if (!"nSupport" %in% names(j)) j$nSupport <- 1L  # no score column: every pair counts once
    j <- j[j$chrom1 == j$chrom2, ]  # Rule 1: drop trans pairs
    if (nrow(j) == 0) return(empty)
    chrom <- resolve_chrom_against_vocab(chrom, unique(j$chrom1))
    j <- j[j$chrom1 == chrom & pmin(j$end1, j$end2) > start & pmax(j$start1, j$start2) < end, ]
    if (nrow(j) == 0) return(empty)
    first_is_left <- j$start1 <= j$start2
    j$intron_start <- ifelse(first_is_left, j$end1, j$end2)
    j$intron_end   <- ifelse(first_is_left, j$start2, j$start1)
    j[, c("intron_start", "intron_end", "nSupport")]
  } else {
    stop("Unrecognised junction/arc file (", n_cols, " columns): ", path,
         " -- expected 12 columns (BED12) or 6-10 columns (BEDPE)")
  }
}

# Sums nSupport across replicate junction files for the SAME condition BEFORE thresholding (not filter-then-union), matching ggsashimi's multi-BAM combining convention; used by viz/sashimi.R's displayed arcs, independent of the transcript-display filtering below.
pool_junctions <- function(samples, chrom, start, end, min_junction_support = MIN_JUNCTION_SUPPORT) {
  if (is.data.frame(samples)) samples <- split(samples, seq_len(nrow(samples)))

  raw_list <- lapply(samples, function(s) read_junction_file(s$junction_path, chrom, start, end))
  raw <- do.call(rbind, raw_list)
  if (nrow(raw) == 0) return(data.frame(intron_start = integer(0), intron_end = integer(0), nSupport = integer(0)))

  pooled <- aggregate(nSupport ~ intron_start + intron_end, data = raw, FUN = sum)
  pooled[pooled$nSupport >= min_junction_support, ]
}

# Returns every annotated transcript for GENE, unfiltered (condition-independent; see gene_display_transcripts() below for which become named); sourced from the TxDb (~0.3s indexed select vs. ~13s GTF zgrep), transcript_name from TX_NAME_LOOKUP_PATH's id->name map (falls back to transcript_id); gtf_path is unused now, kept only for call-site compatibility. See solutions/concepts/19-relocated-code-rationale.md, "utils/transcript_filter.R:197".
gtf_transcript_structure <- function(gene, gtf_path) {
  ens_res <- AnnotationDbi::select(get(ORGDB_PKG, envir = asNamespace(ORGDB_PKG)),
    keys = gene, keytype = "SYMBOL", columns = "ENSEMBL")
  ens_id <- ens_res$ENSEMBL[!is.na(ens_res$ENSEMBL)][1]
  if (length(ens_id) == 0 || is.na(ens_id)) {
    return(data.frame(transcript_id = character(0), transcript_name = character(0),
      tx_id_stripped = character(0)))
  }

  q <- AnnotationDbi::select(txdb, keys = ens_id, keytype = "GENEID", columns = "TXNAME")
  tx_ids <- unique(q$TXNAME)
  tx_id_stripped <- strip_version(tx_ids)

  name_lookup <- readRDS(TX_NAME_LOOKUP_PATH)
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

# `samples` is a data.frame or list of list(sample_id=, junction_path=, quant_path=) for one or more replicates of the SAME condition/gene (junction_path accepted but unused, kept for symmetry with pool_junctions()); TPM is the mean across samples, `proportion` is its share of ALL annotated transcripts' TPM (not just the retained set), so callers can derive a meaningful "Other" residual.
filtered_transcripts <- function(gene, samples, gtf_path) {
  if (is.data.frame(samples)) samples <- split(samples, seq_len(nrow(samples)))

  tx <- gtf_transcript_structure(gene, gtf_path)

  quant <- do.call(rbind, lapply(samples, function(s) {
    d <- read.delim(s$quant_path)
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
  out
}

# Orders retained transcripts ascending by proportion, folding the smallest into "Other" until other_budget or max_segments is exceeded, whichever is stricter (BUGS.md bug 3, validated against real Gas5 data); `other_baseline_tpm` is TPM already excluded by the retained-set filter and always folds into Other. See solutions/concepts/19-relocated-code-rationale.md, "utils/transcript_filter.R:291".
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

# Reconciles several samples' independently-truncated survivor sets into ONE shared named set (deciding survivors per sample first preserves isoform-switching signal a pooled ranking would wash out); the union, only if still over max_segments, is re-capped by each transcript's best proportion across samples. See solutions/concepts/19-relocated-code-rationale.md, "utils/transcript_filter.R:333".
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

# rationale: see solutions/concepts/19-relocated-code-rationale.md, "utils/transcript_filter.R:362"
GENE_DISPLAY_TRANSCRIPTS_CACHE_PATH <- file.path(
  tools::R_user_dir("PlotTrellis", which = "cache"), "gene_display_transcripts_cache.rds"
)

.gene_display_transcripts_cache_load <- function() {
  if (file.exists(GENE_DISPLAY_TRANSCRIPTS_CACHE_PATH)) readRDS(GENE_DISPLAY_TRANSCRIPTS_CACHE_PATH) else list()
}

.gene_display_transcripts_cache_save <- function(cache) {
  dir.create(dirname(GENE_DISPLAY_TRANSCRIPTS_CACHE_PATH), recursive = TRUE, showWarnings = FALSE)
  saveRDS(cache, GENE_DISPLAY_TRANSCRIPTS_CACHE_PATH)
}

# Gene-level (not per-run) transcript filtering, pooled over EVERY RNA-seq condition so the same gene always resolves to the same named set (falls back to tx_ordinal_suffix() ordering with no quant.sf at all); disk-cached per gene, keyed on other_budget/max_segments/samples/data_dir so different projects never share an entry. See solutions/concepts/19-relocated-code-rationale.md, "utils/transcript_filter.R:362".
gene_display_transcripts <- function(gene, gtf_path, samples, data_dir,
                                      other_budget = 0.05, max_segments = length(TX_PALETTE)) {
  # Cache key derived from the samples LIST content (Phase C.3 decision 8),
  # not a file path -- deparse() gives a deterministic string for a given
  # samples list, which is all a cache key needs (no cryptographic strength
  # required). Two different projects (or a project whose samples[] changes
  # via the wizard) querying the same gene symbol still can't share a cache
  # entry, same guarantee the old samples_txt-path key gave. gtf_path is
  # included too (2026-09-08 GTF-consistency audit) -- it wasn't before, so
  # switching a project's GTF while keeping the same samples[] silently
  # served the old GTF's transcript structure for an already-cached gene.
  cache_key <- sprintf("%s|%s|%.6f|%d|%s|%s", gene, gtf_path, other_budget, max_segments,
                        paste(deparse(samples), collapse = ""), data_dir)
  cache <- .gene_display_transcripts_cache_load()
  if (!is.null(cache[[cache_key]])) return(cache[[cache_key]])

  tx <- gtf_transcript_structure(gene, gtf_path)

  pool <- load_samples_from_settings(samples, data_dir)
  pool <- Filter(function(s) !is.na(sample_data_path(s, "quantdir")), pool)

  named <- if (length(pool) == 0) {
    ord <- order(tx_ordinal_suffix(tx$transcript_name), na.last = TRUE)
    head(tx$transcript_name[ord], max_segments)
  } else {
    per_condition <- lapply(pool, function(row) {
      samples <- list(list(sample_id = row$display, junction_path = sample_data_path(row, "paired-end", name = "Junctions"), quant_path = sample_data_path(row, "quantdir")))
      ft   <- filtered_transcripts(gene, samples, gtf_path)
      kept <- ft[ft$mean_tpm >= MIN_DISPLAY_TPM, ]
      truncate_by_proportion(kept, attr(ft, "total_tpm") - sum(kept$mean_tpm), attr(ft, "total_tpm"),
                              other_budget, max_segments)$named
    })
    join_truncated_transcripts(per_condition, max_segments)
  }

  result <- list(named = named, tx = tx)
  cache[[cache_key]] <- result
  tryCatch(.gene_display_transcripts_cache_save(cache), error = function(e) NULL)
  result
}
