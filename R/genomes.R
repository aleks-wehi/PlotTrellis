# rationale: see solutions/concepts/19-relocated-code-rationale.md, "utils/genomes.R:1"; build_tx_name_lookup()/build_tx_meta_lookup()/build_chrom_sizes() are separate R/ files loaded together with this one by library(PlotTrellis), no source() needed.

# rationale: see solutions/concepts/19-relocated-code-rationale.md, "utils/genomes.R:36"
GENOME_PRESETS <- list(
  hg38 = list(
    label        = "Human (GRCh38)",
    genome_label = "hg38",
    orgdb        = "org.Hs.eg.db"
  ),
  mm10 = list(
    label        = "Mouse (GRCm38 / mm10)",
    genome_label = "mm10",
    orgdb        = "org.Mm.eg.db"
  ),
  mm39 = list(
    label        = "Mouse (GRCm39 / GENCODE vM39)",
    genome_label = "mm39",
    orgdb        = "org.Mm.eg.db"
  )
)

# Where genome caches are built: the open project's cache/genome folder (set by
# .pt_set_cache_root(), bug 53), or next to the GTF when no project is open.
GENOME_CACHE_DIR <- NULL

# Derived cache path for a GTF, used when a caller has no path yet; strips ".gtf"/".gtf.gz" (not just the last extension) so "foo.gtf.gz" -> "foo", not "foo.gtf".
# A cache already next to the GTF (built by an older version, or shipped, like
# the example's chrom sizes) is reused read-only; new ones go to GENOME_CACHE_DIR.
.genome_cache_path <- function(gtf_path, suffix) {
  base <- sub("\\.gtf(\\.gz)?$", "", basename(gtf_path))
  beside_gtf <- file.path(dirname(gtf_path), paste0(base, suffix))
  if (is.null(GENOME_CACHE_DIR) || file.exists(beside_gtf)) return(beside_gtf)
  dir.create(GENOME_CACHE_DIR, recursive = TRUE, showWarnings = FALSE)
  file.path(GENOME_CACHE_DIR, paste0(base, suffix))
}

# Any of txdb_path/tx_name_lookup_path/tx_meta_lookup_path/chrom_sizes_path that's missing or doesn't resolve to a real file gets built from gtf_path and written to a derived default path; returns the fully-resolved field list (same shape as a GENOME_PRESETS entry, plus gtf_path/orgdb/rebuilt). Full rationale: solutions/concepts/19-relocated-code-rationale.md, "utils/genomes.R:76".
ensure_genome_ready <- function(gtf_path, orgdb, genome_label,
                                 txdb_path = NULL, tx_name_lookup_path = NULL,
                                 tx_meta_lookup_path = NULL, chrom_sizes_path = NULL) {
  if (is.null(gtf_path) || !nchar(gtf_path) || !file.exists(gtf_path))
    stop("Genome GTF not found: ", if (is.null(gtf_path)) "(not set)" else gtf_path)
  if (is.null(orgdb) || !nchar(orgdb))
    stop("No OrgDb package specified for this genome (e.g. org.Mm.eg.db, org.Hs.eg.db).")
  if (!requireNamespace(orgdb, quietly = TRUE))
    stop("OrgDb package '", orgdb, "' is not installed. Install it first, e.g. ",
         "BiocManager::install('", orgdb, "').")
  # genome_label is just a display string (passed to plotgardener's
  # assembly(Genome = ...) and shown in a couple of UI messages) -- unlike
  # gtf_path/orgdb it isn't load-bearing for anything else, so a genuinely
  # minimal setup (just gtf_path + orgdb) falls back to the OrgDb package
  # name rather than failing.
  if (is.null(genome_label) || !nchar(genome_label)) genome_label <- orgdb

  rebuilt <- FALSE

  # Resolve to the derived default FIRST when the caller supplied nothing
  # (the normal case since 2026-08-25 -- no caller passes these anymore), so
  # the file.exists() check below actually looks at the path that matters.
  # Getting this order backwards means is.null(txdb_path) short-circuits the
  # `||` to TRUE before a derived path is even computed, so file.exists() is
  # only ever checked against a caller-supplied path -- silently forcing a
  # full rebuild on EVERY call when nothing is supplied, defeating caching
  # entirely (a real bug caught here before it shipped: a caller relying on
  # "same gtf_path -> reuse what's already on disk" would otherwise rebuild
  # from scratch on every single app launch).
  if (is.null(txdb_path) || !nchar(txdb_path)) txdb_path <- .genome_cache_path(gtf_path, "_txdb.sqlite")
  if (!file.exists(txdb_path)) {
    message("Building TxDb from ", gtf_path, " (this only happens once)...")
    txdb <- txdbmaker::makeTxDbFromGFF(gtf_path, format = "gtf")
    # GENCODE/Ensembl-style GTFs use bare chrom names ("1", "2", "X"); this
    # app's bigWig/Hi-C fixtures use UCSC-style ("chr1", "chr2", "chrX") --
    # must match, or gene lookups resolve to a seqname that isn't in the
    # data (see BUGS.md #12's root cause). Same recipe as r/pg/build_txdb.R.
    GenomeInfoDb::seqlevelsStyle(txdb) <- "UCSC"
    AnnotationDbi::saveDb(txdb, txdb_path)
    message("TxDb built and saved to ", txdb_path)
    rebuilt <- TRUE
  }

  # Runs on every call, cache hit or miss -- not just right after a fresh
  # build -- so a TxDb cached on disk from before either check below existed
  # (e.g. built under an older copy of this app) gets validated on load too,
  # instead of silently keeping stale/wrong data until someone manually
  # deletes the cache file (feature 34.a: gene search must stay wholly
  # derived from the GTF, not from residual cached data that can drift out
  # of sync with what the current code would produce).
  con <- RSQLite::dbConnect(RSQLite::SQLite(), txdb_path)
  # rationale: see solutions/concepts/19-relocated-code-rationale.md, "utils/genomes.R:139"
  ids  <- RSQLite::dbGetQuery(con, "SELECT DISTINCT gene_id FROM gene")$gene_id
  bare <- sub("\\.[0-9]+$", "", ids)
  versioned <- ids[bare != ids]
  if (length(versioned) > 0) {
    RSQLite::dbBegin(con)
    RSQLite::dbExecute(con, "UPDATE gene SET gene_id = ? WHERE gene_id = ?",
                        params = list(bare[bare != ids], versioned))
    RSQLite::dbCommit(con)
  }
  # rationale: see solutions/concepts/19-relocated-code-rationale.md, "utils/genomes.R:159"
  sample_ids <- RSQLite::dbGetQuery(con, "SELECT DISTINCT gene_id FROM gene LIMIT 20")$gene_id
  RSQLite::dbDisconnect(con)
  if (length(sample_ids) > 0 && !any(grepl("^ENS", sample_ids))) {
    unlink(txdb_path)
    stop(
      "This GTF's gene IDs don't look like Ensembl IDs (expected e.g. 'ENSMUSG...' ",
      "or 'ENSG...', got e.g. '", sample_ids[1], "'). This app's gene search and ",
      "annotation lookups require an Ensembl- or GENCODE-style GTF (they share the ",
      "same ID namespace) -- a RefSeq/NCBI or other symbol/accession-keyed GTF isn't ",
      "supported yet. Offending file: ", gtf_path
    )
  }

  if (is.null(tx_name_lookup_path) || !nchar(tx_name_lookup_path)) tx_name_lookup_path <- .genome_cache_path(gtf_path, "_tx_names.rds")
  if (!file.exists(tx_name_lookup_path)) {
    build_tx_name_lookup(gtf_path, tx_name_lookup_path)
    rebuilt <- TRUE
  }

  if (is.null(tx_meta_lookup_path) || !nchar(tx_meta_lookup_path)) tx_meta_lookup_path <- .genome_cache_path(gtf_path, "_tx_meta.rds")
  if (!file.exists(tx_meta_lookup_path)) {
    build_tx_meta_lookup(gtf_path, tx_meta_lookup_path)
    rebuilt <- TRUE
  }

  # Chromosome lengths: not derivable from the GTF itself (see
  # build_chrom_sizes.R's header comment), so this cache always comes from
  # UCSC, keyed by genome_label -- built once here, independent of whether
  # txdb_path above was just (re)built or already cached, so an existing
  # project's TxDb cache (built before this cache existed) picks up real
  # lengths on its next launch with no manual cache-clearing needed.
  if (is.null(chrom_sizes_path) || !nchar(chrom_sizes_path)) chrom_sizes_path <- .genome_cache_path(gtf_path, "_chrom_sizes.rds")
  if (!file.exists(chrom_sizes_path)) {
    build_chrom_sizes(genome_label, chrom_sizes_path)
    rebuilt <- TRUE
  }

  list(
    gtf_path = gtf_path, orgdb = orgdb, genome_label = genome_label,
    txdb_path = txdb_path, tx_name_lookup_path = tx_name_lookup_path,
    tx_meta_lookup_path = tx_meta_lookup_path, chrom_sizes_path = chrom_sizes_path,
    rebuilt = rebuilt
  )
}
