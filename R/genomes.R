# A project's genome is three flat settings: gtf_path, orgdb and genome_label (the exact
# UCSC build id, e.g. "mm39", since chromosome lengths are looked up by it). Everything else
# is a cache derived from gtf_path (ensure_genome_ready() below), never a setting.

# The genomes the wizard offers as a quick fill of genome_label and orgdb; mm39 is
# also a new project's default. Never read at render time. No gtf_path: the GTF is always
# the user's own file.
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

# Derived cache path for a GTF, used when a caller has no path yet; strips ".gtf"/".gtf.gz" (not just the last extension) so "foo.gtf.gz" -> "foo", not "foo.gtf".
# A cache already next to the GTF (built by an older version, or shipped, like
# the example's chrom sizes) is reused read-only; new ones go to the project's
# cache folder's genome/ (cache_dir), or next to the GTF when there is no project.
.genome_cache_path <- function(gtf_path, suffix, cache_dir = NULL) {
  base <- sub("\\.gtf(\\.gz)?$", "", basename(gtf_path))
  beside_gtf <- file.path(dirname(gtf_path), paste0(base, suffix))
  if (is.null(cache_dir) || file.exists(beside_gtf)) return(beside_gtf)
  genome_dir <- file.path(cache_dir, "genome")
  dir.create(genome_dir, recursive = TRUE, showWarnings = FALSE)
  file.path(genome_dir, paste0(base, suffix))
}

# The console line for a first-time TxDb build. A full genome annotation took 71 s (GENCODE vM25,
# mouse) to 155 s (GENCODE v46, human) for the whole first build on a compute node; the bundled
# example GTFs take about a second, so the estimate is only given for a large file.
.pt_txdb_build_message <- function(gtf_path) {
  note <- "this only happens once per GTF file"
  if (isTRUE(file.size(gtf_path) > 1e7))
    note <- paste0(note, ", and takes about 1-3 minutes, longer on a slow disk")
  paste0("Building TxDb from ", gtf_path, " (", note, ")...")
}

# The genome's caches (TxDb, transcript name and meta lookups, chromosome sizes), at paths
# derived from gtf_path (.genome_cache_path()): a file already there is reused, a missing one
# is built from gtf_path. Returns those paths with gtf_path/orgdb/genome_label, and `rebuilt`
# (whether anything was just built).
ensure_genome_ready <- function(gtf_path, orgdb, genome_label, cache_dir = NULL) {
  if (is.null(gtf_path) || !nchar(gtf_path) || !file.exists(gtf_path))
    stop("Genome GTF not found: ", if (is.null(gtf_path)) "(not set)" else gtf_path, call. = FALSE)
  if (is.null(orgdb) || !nchar(orgdb))
    stop("No OrgDb package specified for this genome (e.g. org.Mm.eg.db, org.Hs.eg.db).", call. = FALSE)
  if (!requireNamespace(orgdb, quietly = TRUE))
    stop("OrgDb package '", orgdb, "' is not installed. Install it first, e.g. ",
         "BiocManager::install('", orgdb, "').", call. = FALSE)
  # genome_label is the figure's assembly name (plotgardener's
  # assembly(Genome = ...)), shown in a couple of UI messages, and the key
  # build_chrom_sizes() looks lengths up by. A blank one falls back to the
  # OrgDb package name rather than failing here (the chromosome-length
  # lookup below then errors, naming genome_label, unless lengths are
  # already cached).
  if (is.null(genome_label) || !nchar(genome_label)) genome_label <- orgdb

  rebuilt <- FALSE

  txdb_path <- .genome_cache_path(gtf_path, "_txdb.sqlite", cache_dir)
  if (!file.exists(txdb_path)) {
    message(.pt_txdb_build_message(gtf_path))
    txdb <- suppressMessages(txdbmaker::makeTxDbFromGFF(gtf_path, format = "gtf"))
    # GENCODE/Ensembl-style GTFs use bare chrom names ("1", "2", "X"); this
    # app's bigWig/Hi-C fixtures use UCSC-style ("chr1", "chr2", "chrX") --
    # must match, or gene lookups resolve to a seqname that isn't in the
    # data.
    GenomeInfoDb::seqlevelsStyle(txdb) <- "UCSC"
    AnnotationDbi::saveDb(txdb, txdb_path)
    rebuilt <- TRUE
  }

  # Runs on every call, cache hit or miss -- not just right after a fresh
  # build -- so a TxDb cached on disk from before either check below existed
  # (e.g. built under an older copy of this app) gets validated on load too,
  # instead of silently keeping stale/wrong data until someone manually
  # deletes the cache file: gene search must stay wholly derived from the
  # GTF, not from residual cached data that can drift out of sync with what
  # the current code would produce.
  con <- RSQLite::dbConnect(RSQLite::SQLite(), txdb_path)
  # GENCODE/Ensembl gene_ids carry a version (e.g. "ENSMUSG00000027186.7"); the OrgDb's
  # ENSEMBL keys, which every gene lookup joins on, don't. Strip it in the saved SQLite
  # file (a TxDb has no API to rewrite its table); the GTF itself is never touched.
  ids  <- RSQLite::dbGetQuery(con, "SELECT DISTINCT gene_id FROM gene")$gene_id
  bare <- sub("\\.[0-9]+$", "", ids)
  versioned <- ids[bare != ids]
  if (length(versioned) > 0) {
    RSQLite::dbBegin(con)
    RSQLite::dbExecute(con, "UPDATE gene SET gene_id = ? WHERE gene_id = ?",
                        params = list(bare[bare != ids], versioned))
    RSQLite::dbCommit(con)
  }
  # Every gene lookup joins this TxDb's gene_id on the OrgDb's ENSEMBL keys, so a GTF
  # keyed otherwise (e.g. RefSeq, gene_id = symbol) stops here with a message, not with
  # a cryptic AnnotationDbi error at the first gene search.
  sample_ids <- RSQLite::dbGetQuery(con, "SELECT DISTINCT gene_id FROM gene LIMIT 20")$gene_id
  RSQLite::dbDisconnect(con)
  if (length(sample_ids) > 0 && !any(grepl("^ENS", sample_ids))) {
    unlink(txdb_path)
    stop(
      "This GTF's gene IDs don't look like Ensembl IDs (expected e.g. 'ENSMUSG...' ",
      "or 'ENSG...', got e.g. '", sample_ids[1], "'). This app's gene search and ",
      "annotation lookups require an Ensembl- or GENCODE-style GTF (they share the ",
      "same ID namespace) -- a RefSeq/NCBI or other symbol/accession-keyed GTF isn't ",
      "supported yet. Offending file: ", gtf_path, call. = FALSE
    )
  }

  tx_name_lookup_path <- .genome_cache_path(gtf_path, "_tx_names.rds", cache_dir)
  if (!file.exists(tx_name_lookup_path)) {
    build_tx_name_lookup(gtf_path, tx_name_lookup_path)
    rebuilt <- TRUE
  }

  tx_meta_lookup_path <- .genome_cache_path(gtf_path, "_tx_meta.rds", cache_dir)
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
  chrom_sizes_path <- .genome_cache_path(gtf_path, "_chrom_sizes.rds", cache_dir)
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
