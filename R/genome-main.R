# A project's genome is three flat settings: gtf_path, orgdb and genome_label
# (the exact UCSC build id, e.g. "mm39", since chromosome lengths are looked up
# by it). Everything else is a cache derived from gtf_path
# (ensure_genome_ready() below), never a setting.

# The genomes the wizard offers as a quick fill of genome_label and orgdb; mm39
# is also a new project's default. Never read at render time. No gtf_path: the
# GTF is always the user's own file.
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

# A genome cache file for a GTF, in the project's cache folder (cache_dir):
# <cache_dir>/genome/<hash>/<GTF name without .gtf/.gtf.gz><suffix>. The hash is
# of the GTF's full path, size and modification time, so a replaced or edited
# GTF gets fresh caches and two GTFs with the same name never share them.
# Nothing is read from or written beside the GTF.
.genome_cache_path <- function(gtf_path, suffix, cache_dir) {
  base <- sub("\\.gtf(\\.gz)?$", "", basename(gtf_path))
  genome_dir <- file.path(.pt_cache_root(cache_dir), "genome",
                          .pt_gtf_hash(gtf_path))
  dir.create(genome_dir, recursive = TRUE, showWarnings = FALSE)
  file.path(genome_dir, paste0(base, suffix))
}

# First 12 hex digits of the md5 of the GTF's full path, size and modification
# time (to the millisecond, written with sprintf() so options(digits)/OutDec
# can't change the key).
.pt_gtf_hash <- function(gtf_path) {
  info <- file.info(gtf_path)
  key <- tempfile(tmpdir = tempdir(check = TRUE))
  on.exit(unlink(key))
  writeLines(sprintf("%s %.0f %.3f", normalizePath(gtf_path, winslash = "/"),
                     info$size, as.numeric(info$mtime)), key)
  substr(unname(tools::md5sum(key)), 1, 12)
}

# A chromosome-sizes file shipped beside a GTF inside the PlotTrellis package
# (the examples', data-raw/build-example-data*.R), so an example works offline;
# NULL for any other GTF.
.pt_shipped_chrom_sizes <- function(gtf_path) {
  pkg <- normalizePath(system.file(package = "PlotTrellis"), winslash = "/",
                       mustWork = FALSE)
  dir <- normalizePath(dirname(gtf_path), winslash = "/", mustWork = FALSE)
  if (!nzchar(pkg) || !startsWith(dir, paste0(pkg, "/"))) return(NULL)
  shipped <- file.path(dir,
                       paste0(sub("\\.gtf(\\.gz)?$", "", basename(gtf_path)),
                              "_chrom_sizes.rds"))
  if (file.exists(shipped)) shipped
}

# The console line for a first-time TxDb build. A full genome annotation took
# 71 s (GENCODE vM25, mouse) to 155 s (GENCODE v46, human) for the whole first
# build on a compute node; the bundled example GTFs take about a second, so the
# estimate is only given for a large file.
.pt_txdb_build_message <- function(gtf_path) {
  note <- "this only happens once per GTF file"
  if (isTRUE(file.size(gtf_path) > 1e7))
    note <- paste0(note, ", and takes about 1-3 minutes, longer on a slow disk")
  paste0("Building TxDb from ", gtf_path, " (", note, ")...")
}

# An Ensembl GTF names chromosomes "1", "X", "MT"; the data files, chromosome
# lengths and region box use UCSC names ("chr1", "chrX", "chrM"). The TxDb's
# names are rewritten in its SQLite file (`txdb_path`), not with
# seqlevelsStyle<-: saveDb() does not keep that renaming, and select() and
# genes(filter =), which plotgardener's gene tracks use, ignore it anyway. A
# name with no UCSC equivalent (an unplaced scaffold, "JH584299.1") is kept, and
# so are all of them when none has one. A TxDb with any "chr" name (GENCODE) is
# left as it is.
# Several R processes may open the same project at once (a first open after a
# build, batch jobs): the rename takes the write lock first and looks again, so
# whichever process is first renames and the others find it done; a process
# kept waiting longer than `busy_ms` gets one warning, as does a file that can't
# be written (a project shared read-only), never an error.
# Returns (invisibly) "renamed", "done elsewhere", "unchanged" or "failed".
.pt_txdb_ucsc_chroms <- function(txdb_path, busy_ms = 10000) {
  # synchronous = NULL: RSQLite's connect-time PRAGMA would fail with its own
  # "database is locked" warning while another process is committing (the busy
  # handler comes after connect).
  con <- RSQLite::dbConnect(RSQLite::SQLite(), txdb_path, synchronous = NULL)
  on.exit(RSQLite::dbDisconnect(con))
  RSQLite::sqliteSetBusyHandler(con, busy_ms)
  chroms <-
    function() RSQLite::dbGetQuery(con, "SELECT chrom FROM chrominfo")$chrom
  ucsc <- function(old) tryCatch({
    si <- GenomeInfoDb::Seqinfo(old)
    GenomeInfoDb::seqlevelsStyle(si) <- "UCSC"
    GenomeInfoDb::seqnames(si)
  }, error = function(e) old)
  old <- chroms()
  if (length(old) == 0 || any(startsWith(old, "chr")) ||
      !any(ucsc(old) != old))
    return(invisible("unchanged"))
  result <- tryCatch({
    # Room for the whole rename (about a million rows on a full annotation) in
    # memory, so SQLite rewrites the file in one burst at the commit instead of
    # part-way through: other processes wait less, and AnnotationDbi's lock-free
    # readers see the old or the new file.
    RSQLite::dbExecute(con, "PRAGMA cache_size = -200000")
    # the write lock; waits for another renamer
    RSQLite::dbExecute(con, "BEGIN IMMEDIATE")
    # looked at again under the lock
    old <- chroms()
    if (length(old) == 0 || any(startsWith(old, "chr")) ||
        !any(ucsc(old) != old)) {
      RSQLite::dbExecute(con, "ROLLBACK")
      "done elsewhere"
    } else {
      new <- ucsc(old)
      changed <- new != old
      for (sql in c("UPDATE chrominfo SET chrom = ? WHERE chrom = ?",
                    "UPDATE transcript SET tx_chrom = ? WHERE tx_chrom = ?",
                    "UPDATE exon SET exon_chrom = ? WHERE exon_chrom = ?",
                    "UPDATE cds SET cds_chrom = ? WHERE cds_chrom = ?"))
        RSQLite::dbExecute(con, sql, params = list(new[changed], old[changed]))
      RSQLite::dbExecute(con, "COMMIT")
      "renamed"
    }
  }, error = function(e) {
    try(RSQLite::dbExecute(con, "ROLLBACK"), silent = TRUE)
    msg <- conditionMessage(e)
    why <- if (grepl("locked|busy", msg, ignore.case = TRUE))
      "another R process is renaming them right now; open the project again in a moment"
    else if (grepl("readonly|read-only|permission", msg, ignore.case = TRUE))
      "the file is read-only; open the project once as a user who can write it"
    else msg
    warning("PlotTrellis could not give the chromosomes in ", txdb_path,
            " their UCSC names (", why,
            "), so genes may not show in region mode.", call. = FALSE)
    "failed"
  })
  invisible(result)
}

# A TxDb built from the GTF. txdbmaker's two routine warnings for a
# GENCODE/Ensembl GTF (stop_codon rows with a phase; no genome version) are
# dropped, so options(warn = 2) does not stop a first open on them; any other
# warning still reaches the user. Newer txdbmaker words the second one "genome
# assembly is not available".
.pt_make_txdb <- function(gtf_path) {
  # txdbmaker lists every type with a phase ("features of type Selenocysteine,
  # stop_codon").
  routine <- c("features of type\\s+([A-Za-z_]+,\\s+)*stop_codon",
               "genome (version information|assembly) is not available")
  withCallingHandlers(
    suppressMessages(txdbmaker::makeTxDbFromGFF(gtf_path, format = "gtf")),
    warning = function(w) {
      if (any(vapply(routine, grepl, logical(1), conditionMessage(w))))
        invokeRestart("muffleWarning")
    })
}

# The genome's caches (TxDb, transcript name and meta lookups, chromosome
# sizes), in the project's cache folder `cache_dir` at paths derived from
# gtf_path (.genome_cache_path()): a file already there is reused, a missing one
# is built from gtf_path. Returns those paths with gtf_path/orgdb/genome_label,
# and `rebuilt` (whether anything was just built).
ensure_genome_ready <- function(gtf_path, orgdb, genome_label, cache_dir) {
  if (is.null(gtf_path) || !nchar(gtf_path) || !file.exists(gtf_path))
    stop("Genome GTF not found: ",
         if (is.null(gtf_path)) "(not set)" else gtf_path, call. = FALSE)
  if (is.null(orgdb) || !nchar(orgdb))
    stop(
      "No OrgDb package specified for this genome (e.g. org.Mm.eg.db, org.Hs.eg.db).",
      call. = FALSE
    )
  if (!requireNamespace(orgdb, quietly = TRUE))
    stop("OrgDb package '", orgdb,
         "' is not installed. Install it first, e.g. ",
         "BiocManager::install('", orgdb, "').", call. = FALSE)
  # genome_label is the figure's assembly name (plotgardener's
  # assembly(Genome = ...)), shown in a couple of UI messages, and the key
  # build_chrom_sizes() looks lengths up by. A blank one falls back to the
  # OrgDb package name rather than failing here (the chromosome-length
  # lookup below then errors, naming genome_label, unless lengths are
  # already cached).
  if (is.null(genome_label) || !nchar(genome_label)) genome_label <- orgdb

  files <- .pt_genome_files(gtf_path, cache_dir)

  # Chromosome lengths: not derivable from the GTF itself (see
  # genome-build-chrom-sizes.R's header comment), so this cache always comes from
  # UCSC, keyed by genome_label -- built once here, independent of whether
  # txdb_path above was just (re)built or already cached, so an existing
  # project's TxDb cache (built before this cache existed) picks up real
  # lengths on its next launch with no manual cache-clearing needed.
  chrom_sizes_path <- .genome_cache_path(gtf_path, "_chrom_sizes.rds",
                                         cache_dir)
  rebuilt <- files$rebuilt
  if (!file.exists(chrom_sizes_path)) {
    shipped <- .pt_shipped_chrom_sizes(gtf_path)
    if (is.null(shipped)) build_chrom_sizes(genome_label, chrom_sizes_path)
    else file.copy(shipped, chrom_sizes_path, copy.mode = FALSE)
    rebuilt <- TRUE
  }

  list(
    gtf_path = gtf_path, orgdb = orgdb, genome_label = genome_label,
    txdb_path = files$txdb_path,
    tx_name_lookup_path = files$tx_name_lookup_path,
    tx_meta_lookup_path = files$tx_meta_lookup_path,
    chrom_sizes_path = chrom_sizes_path,
    rebuilt = rebuilt
  )
}

# The genome caches built from the GTF alone (the TxDb, the transcript name
# and meta lookups and the gene names), in `cache_dir` at paths derived from
# gtf_path: a file already there is reused, a missing one is built. Shared by
# ensure_genome_ready() and a downloaded script (.pt_script_genome()). Returns
# their paths and `rebuilt` (whether anything was just built).
.pt_genome_files <- function(gtf_path, cache_dir) {
  rebuilt <- FALSE

  txdb_path <- .genome_cache_path(gtf_path, "_txdb.sqlite", cache_dir)
  if (!file.exists(txdb_path)) {
    message(.pt_txdb_build_message(gtf_path))
    txdb <- .pt_make_txdb(gtf_path)
    # Saved under a temporary name in the same folder, then renamed into place,
    # so a build stopped part-way, or a second R process building the same TxDb
    # at once, never leaves a half-written file at txdb_path. A file another
    # process put there first is kept. Only a process killed mid-save leaves its
    # temporary file behind (cleared with the cache).
    tmp <- tempfile(paste0(basename(txdb_path), "-"),
                    tmpdir = dirname(txdb_path))
    on.exit(unlink(tmp), add = TRUE)
    AnnotationDbi::saveDb(txdb, tmp)
    if (!file.exists(txdb_path) && !file.rename(tmp, txdb_path) &&
        !file.exists(txdb_path))
      stop("PlotTrellis could not save the TxDb it built to ", txdb_path, ".",
           call. = FALSE)
    unlink(tmp)
    rebuilt <- TRUE
  }

  # Runs on every call, cache hit or miss -- not just right after a fresh
  # build -- so a TxDb cached on disk from before either check below existed
  # (e.g. built under an older copy of this app) gets validated on load too,
  # instead of silently keeping stale/wrong data until someone manually
  # deletes the cache file: gene search must stay wholly derived from the
  # GTF, not from residual cached data that can drift out of sync with what
  # the current code would produce.
  con <- RSQLite::dbConnect(RSQLite::SQLite(), txdb_path, synchronous = NULL)
  # Another process may be renaming the chromosomes (.pt_txdb_ucsc_chroms()) or
  # stripping versions right now: wait for it rather than fail.
  RSQLite::sqliteSetBusyHandler(con, 10000)
  # GENCODE/Ensembl gene_ids carry a version (e.g. "ENSMUSG00000027186.7"); the
  # OrgDb's ENSEMBL keys, which every gene lookup joins on, don't. Strip it in
  # the saved SQLite file (a TxDb has no API to rewrite its table); the GTF
  # itself is never touched.
  ids  <- RSQLite::dbGetQuery(con, "SELECT DISTINCT gene_id FROM gene")$gene_id
  bare <- sub("\\.[0-9]+$", "", ids)
  versioned <- ids[bare != ids]
  if (length(versioned) > 0) {
    RSQLite::dbBegin(con)
    RSQLite::dbExecute(con, "UPDATE gene SET gene_id = ? WHERE gene_id = ?",
                        params = list(bare[bare != ids], versioned))
    RSQLite::dbCommit(con)
  }
  # Every gene lookup joins this TxDb's gene_id on the OrgDb's ENSEMBL keys, so
  # a GTF keyed otherwise (e.g. RefSeq, gene_id = symbol) stops here with a
  # message, not with a cryptic AnnotationDbi error at the first gene search.
  sample_ids <- RSQLite::dbGetQuery(
    con, "SELECT DISTINCT gene_id FROM gene LIMIT 20"
  )$gene_id
  RSQLite::dbDisconnect(con)
  if (length(sample_ids) > 0 && !any(grepl("^ENS", sample_ids))) {
    unlink(txdb_path)
    stop(
      "This GTF's gene IDs don't look like Ensembl IDs (expected e.g. 'ENSMUSG...' ",
      "or 'ENSG...', got e.g. '", sample_ids[1],
      "'). This app's gene search and ",
      "annotation lookups require an Ensembl- or GENCODE-style GTF (they share the ",
      "same ID namespace) -- a RefSeq/NCBI or other symbol/accession-keyed GTF isn't ",
      "supported yet. Offending file: ", gtf_path, call. = FALSE
    )
  }
  .pt_txdb_ucsc_chroms(txdb_path)

  tx_name_lookup_path <- .genome_cache_path(gtf_path, "_tx_names.rds",
                                            cache_dir)
  if (!file.exists(tx_name_lookup_path)) {
    build_tx_name_lookup(gtf_path, tx_name_lookup_path)
    rebuilt <- TRUE
  }

  tx_meta_lookup_path <- .genome_cache_path(gtf_path, "_tx_meta.rds", cache_dir)
  if (!file.exists(tx_meta_lookup_path)) {
    build_tx_meta_lookup(gtf_path, tx_meta_lookup_path)
    rebuilt <- TRUE
  }

  gene_names_path <- .pt_gene_names_file(txdb_path)
  if (!file.exists(gene_names_path)) {
    message("Reading the gene names in ", basename(gtf_path),
            " (once per GTF file)...")
    build_gene_name_lookup(gtf_path, gene_names_path)
    rebuilt <- TRUE
  }

  list(txdb_path = txdb_path, tx_name_lookup_path = tx_name_lookup_path,
       tx_meta_lookup_path = tx_meta_lookup_path, rebuilt = rebuilt)
}

# The genome files a downloaded script draws from: the project's cached TxDb
# and transcript lookups (settings) while they are all there, the fast route;
# where they are not (another computer, or the cache cleared), built from
# settings$gtf_path into `cache_dir` as pt_project() builds them, once, and
# reused from there on the next run; `needs`: the R packages only that build
# uses, checked before it. Returns `settings` pointing at them.
.pt_script_genome <- function(settings, cache_dir, needs = character(0)) {
  files <- c("txdb_path", "tx_name_lookup_path", "tx_meta_lookup_path")
  paths <- unlist(settings[files])
  if (length(paths) == length(files) &&
      all(file.exists(c(paths, .pt_gene_names_file(paths[[1]])))))
    return(settings)
  gtf <- settings$gtf_path %||% ""
  if (!file.exists(gtf))
    stop("The GTF '", gtf, "' is not here: set gtf_path in `settings` to ",
         "your copy of it.", call. = FALSE)
  # Said once, on the run that builds it.
  if (!file.exists(.genome_cache_path(gtf, "_txdb.sqlite", cache_dir))) {
    for (pkg in needs)
      if (!requireNamespace(pkg, quietly = TRUE))
        stop("Building the genome cache needs the R package '", pkg,
             "'. Install it with BiocManager::install(\"", pkg, "\").",
             call. = FALSE)
    message("The project's genome cache is not here: building it from the ",
            "GTF into ", cache_dir, " in the current folder.")
  }
  g <- .pt_genome_files(gtf, cache_dir)
  settings[files] <- g[files]
  settings$cache_dir <- cache_dir
  settings
}

# The gene-name lookup (build_gene_name_lookup()) kept beside a genome cache's
# TxDb file.
.pt_gene_names_file <- function(txdb_path)
  sub("_txdb\\.sqlite$", "_gene_names.rds", txdb_path)
