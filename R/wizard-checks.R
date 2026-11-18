# The setup wizard's constants and pure checks: the data types a sample's files
# can be, and the checks behind its live field status lines and its Save (also
# run by parse_samples_txt() and when a project opens). No Shiny state: each
# takes plain values and returns plain values.

# PALETTE[1] (api-figure-helpers.R); the same for every new sample, not cycled
WIZARD_DEFAULT_COLOUR <- "#1f78b4"

WIZARD_DATA_TYPES <- list(
  bw          = list(label = "BigWig coverage (ATAC / ChIP-seq / CUT&RUN / RNA-seq)",
                      ext  = "\\.(bw|bigWig)$"),
  hic         = list(label = "Hi-C contact map",
                      ext  = "\\.(hic|cool|mcool)$"),
  bed         = list(label = "BED intervals",
                      ext  = "\\.(bed|narrowPeak|broadPeak|gappedPeak)(\\.gz)?$"),
  `paired-end` = list(label = paste("Paired-interval arcs (BED12 or BED6 splice",
                                     "junctions, a STAR SJ.out.tab, or BEDPE or",
                                     "HiCCUPS loop calls)"),
                       ext  = paste0("\\.(bed|bed12|bedpe)(\\.gz)?$",
                                     "|SJ\\.out\\.tab(\\.gz)?$"))
)
# What the Add-sample dialog calls each kind of file (the type keys stay as they
# are).
WIZ_KIND_LABELS <- c(bw = "Coverage (bigWig)", hic = "Hi-C map",
                     bed = "Intervals (BED)",
                     `paired-end` = "Arcs (junctions, loops)")

# Returns paths relative to data_root (for populating an instance's path
# selectize); quantdir instead lists directories containing quant.sf, not
# extension-matched files.
wizard_scoped_choices <- function(data_root, type) {
  if (is.null(data_root) || !nchar(data_root) || !dir.exists(data_root))
    return(character(0))
  if (identical(type, "quantdir")) {
    quant_files <- list.files(data_root, pattern = "^quant\\.sf$",
                              recursive = TRUE, full.names = TRUE)
    dirs <- unique(dirname(quant_files))
    return(substring(dirs, nchar(data_root) + 2))
  }
  ext <- WIZARD_DATA_TYPES[[type]]$ext
  files <- list.files(data_root, pattern = ext, recursive = TRUE,
                      full.names = TRUE, ignore.case = TRUE)
  substring(files, nchar(data_root) + 2)
}

# Why a file chosen as paired-end can't be drawn as arcs, or NULL when it can.
# A cheap fail-fast check (is_paired_end_structural()) for the formats
# read_junction_file() reads. A peaks file (narrowPeak or similar) is pointed to
# the BED type, and so is a 3 to 5 column file, if it is BED intervals (it may be
# something else typed in by hand, such as a bedGraph).
wiz_paired_end_problem <- function(path) {
  if (is_paired_end_structural(path)) return(NULL)
  fields <- .pt_arc_first_data(path)$fields
  if (length(fields) == 0)
    return("has no data lines (only blank, comment or header lines)")
  if (.pt_is_bed_peaks(fields))
    return(paste("looks like a BED peaks file (narrowPeak or similar: a strand",
                 "in column 6); add it as Intervals (BED) instead"))
  if (length(fields) >= 3 && length(fields) <= 5)
    return(sprintf(paste("has %d columns (arcs need 6 or more); if it is BED",
                         "intervals, add it as Intervals (BED) instead"),
                   length(fields)))
  paste("is not a file the arc tracks read: BED12 or BED6 splice junctions,",
        "a STAR SJ.out.tab, a BEDPE (6 to 10 columns) or a HiCCUPS loop file")
}

# A literal "\u00b7" not surrounded by spaces breaks format_track_label()'s
# bolding regex -- reject at entry time.
wizard_name_has_bad_dot <- function(name) {
  grepl("\u00b7", name, fixed = TRUE) && !grepl("\\s\u00b7\\s", name)
}

# Same idea as clean_optional_numeric(): normalizes a blank input to NULL, not
# "", matching how an unset path is represented in the settings (CFG) and
# avoiding a spurious is_dirty()/export delta trip against NULL.
wiz_clean_optional_text <- function(v) {
  v <- trimws(v %||% "")
  if (nchar(v) == 0) NULL else v
}

# The message for a blank genome label, shared by Save/import
# (wiz_save_problems()) and opening a project (.pt_require_genome_label()).
PT_GENOME_LABEL_REQUIRED <-
  "Genome label is required: the UCSC name of the genome build, e.g. mm39, mm10 or hg38."

# Stops with PT_GENOME_LABEL_REQUIRED when a project's settings have a blank
# genome label (a hand-edited project_settings.json can: the wizard's Save and
# the samples.txt import refuse one), before the first launch of a new GTF fails
# on it with a console error.
.pt_require_genome_label <- function(genome_label, settings_json) {
  if (is.null(wiz_clean_optional_text(genome_label)))
    stop(PT_GENOME_LABEL_REQUIRED, " Set \"genome_label\" in ", settings_json,
         ", or fill it in again with plottrellis_wizard().", call. = FALSE)
}

# The checks the wizard's Save runs before writing project_settings.json, shared
# with parse_samples_txt() so the two can't drift. Returns the problems, as the
# wizard's own messages; character(0) means OK. An OrgDb package that isn't
# installed otherwise only fails at the first cold start, after Save restarts
# the app, with the hint in the console only. A genome_label that GenomeInfoDb
# can't look chromosome sizes up by (e.g. "m10") otherwise only fails at the
# first cold start, after Save; so does a blank one (the first launch of a new
# GTF stops), so both are problems.
wiz_save_problems <- function(data_root, samples, gtf_path, orgdb,
                              genome_label = NULL) {
  errs <- character(0)
  if (is.null(wiz_clean_optional_text(data_root)))
    errs <- c(errs, "Data root is required.")
  else if (!dir.exists(data_root))
    errs <- c(errs, "Data root directory does not exist.")
  if (length(samples) == 0)
    errs <- c(errs, "Add at least one sample before saving.")
  if (is.null(wiz_clean_optional_text(gtf_path)))
    errs <- c(errs, "GTF path is required.")
  orgdb <- wiz_clean_optional_text(orgdb)
  if (is.null(orgdb)) errs <- c(errs, "OrgDb package is required.")
  else if (!requireNamespace(orgdb, quietly = TRUE))
    errs <- c(errs, sprintf("OrgDb package '%s' is not installed. Install it first, e.g. BiocManager::install('%s').",
                            orgdb, orgdb))
  label <- wiz_clean_optional_text(genome_label)
  if (is.null(label))
    errs <- c(errs, PT_GENOME_LABEL_REQUIRED)
  else if (!label %in% c(GenomeInfoDb::registered_UCSC_genomes()$genome,
                         GenomeInfoDb::registered_NCBI_assemblies()$assembly))
    errs <- c(errs, sprintf("Genome label '%s' isn't a known genome build id (e.g. mm10, mm39, hg38).",
                            label))
  errs
}

# The species PlotTrellis supports, by Ensembl gene-ID prefix, with the OrgDb
# package each needs: what the wizard's live GTF and OrgDb checks compare
# against. Mouse and human only (ORGDB_TAXON_ID, the UniProt links, map only
# these two).
WIZ_SPECIES_BY_ID_PREFIX <- list(
  ENSG    = list(species = "human", orgdb = "org.Hs.eg.db"),
  ENSMUSG = list(species = "mouse", orgdb = "org.Mm.eg.db")
)

# The UCSC genome label for each build a GENCODE "##description:" or Ensembl
# "#!genome-build" header line can name, for the live genome-label check.
WIZ_LABEL_BY_BUILD <- c(GRCm39 = "mm39", GRCm38 = "mm10", GRCh38 = "hg38",
                        GRCh37 = "hg19")

# What the wizard says under the GTF path as it is typed, before Save. Reads
# only the file's first lines (gzipped or not), so it is cheap on a 1 GB file. A
# warning here never blocks Save: ensure_genome_ready() stops with the full
# message at the first launch. Returns list(state = "blank" | "missing" |
# "warn" | "ok", text, species = a WIZ_SPECIES_BY_ID_PREFIX entry or NULL,
# build = the genome build named in the header, e.g. "GRCm38", or NULL).
wiz_gtf_check <- function(path, n_lines = 2000) {
  res <- function(state, text = NULL, species = NULL, build = NULL)
    list(state = state, text = text, species = species, build = build)
  path <- wiz_clean_optional_text(path)
  if (is.null(path)) return(res("blank"))
  if (!is_real_file(path))
    return(res("missing",
               "No file at this path (it must be a path on the computer running R)."))
  if (!grepl("\\.gtf(\\.gz)?$", path, ignore.case = TRUE))
    return(res("warn",
               "Expected a .gtf or .gtf.gz file. A GFF3 (.gff, .gff3) or other annotation format is not read."))
  # A corrupt .gz warns as well as errors; neither reaches the console.
  lines <- tryCatch({
    con <- gzfile(path)
    on.exit(close(con), add = TRUE)
    suppressWarnings(readLines(con, n = n_lines, warn = FALSE))
  }, error = function(e) character(0))
  header <- lines[startsWith(lines, "#")]
  body   <- lines[!startsWith(lines, "#") & nzchar(lines)]
  # The gene_id attribute itself, not another ending in gene_id (StringTie's
  # ref_gene_id).
  id_re <- '(^|[; \t])gene_id "([^"]*)"'
  ids <- regmatches(body, regexpr(id_re, body))
  ids <- sub(id_re, "\\2", ids)
  if (length(ids) == 0)
    return(res("warn", sprintf(
      "No gene_id attribute in the first %d lines. Is this a GTF file?",
      n_lines)))
  if (!any(startsWith(ids, "ENS")))
    return(res("warn", sprintf(paste0(
      "The gene IDs here look like '%s', not Ensembl IDs (ENSG..., ENSMUSG...). This is ",
      "probably a RefSeq/NCBI GTF, which PlotTrellis cannot use: gene search joins the GTF ",
      "to the OrgDb package on Ensembl IDs. Download the GENCODE or Ensembl GTF for the ",
      "same genome build instead (see 'Which GTF?' below)."), ids[1])))
  ens_ids <- ids[startsWith(ids, "ENS")]
  # GENCODE's GRCh37-mapped ("lift37") files suffix every ID with _<n>, which no
  # version stripping removes, so no gene name ever matches.
  if (any(grepl("_[0-9]+$", ens_ids)))
    return(res("warn", sprintf(paste0(
      "Gene IDs like '%s' are from a GENCODE \"lift37\" GTF (a newer release mapped back to ",
      "GRCh37). PlotTrellis cannot match them to gene names: use the GENCODE 19 GTF instead."),
      ens_ids[1])))
  prefix  <- sub("[0-9].*$", "", ens_ids[1])
  species <- WIZ_SPECIES_BY_ID_PREFIX[[prefix]]
  # GENCODE's "##description:" line names the genome and release ("mouse genome
  # (GRCm38), version M25 (Ensembl 100)"); Ensembl's "#!genome-build" names the
  # build.
  desc <- sub("^##description: (evidence-based annotation of the )?", "",
              grep("^##description: ", header, value = TRUE))
  ens_build <- sub("^#!genome-build ", "",
                   grep("^#!genome-build ", header, value = TRUE))
  build <- unlist(regmatches(c(desc, ens_build),
                             gregexpr("GRC[mh]3[0-9]", c(desc, ens_build))))
  build <- if (length(build)) build[1] else NULL
  sp <- if (!is.null(species)) paste0(" (", species$species, ")") else ""
  src <- if (any(grepl("^##provider: GENCODE", header))) {
    paste0("a GENCODE GTF", if (length(desc)) paste0(": ", desc[1]) else sp)
  } else if (length(ens_build)) {
    paste0("an Ensembl GTF, genome build ", ens_build[1], sp)
  } else paste0("a GTF with Ensembl gene IDs", sp)
  if (is.null(species))
    return(res("warn", paste0("Looks like ", src, ", but the gene IDs (",
                              prefix, "...) are not mouse or human. ",
                              "PlotTrellis supports mouse and human only."),
               build = build))
  res("ok", paste0("Looks like ", src,
                   ". Make sure it is the release your reads were quantified with."),
      species = species, build = build)
}

# What the wizard says under the OrgDb field: whether the package is installed,
# and whether it is the one the GTF's gene IDs call for (`gtf_species`, from
# wiz_gtf_check()).
wiz_orgdb_check <- function(orgdb, gtf_species = NULL) {
  orgdb <- wiz_clean_optional_text(orgdb)
  if (is.null(orgdb)) return(list(state = "blank", text = NULL))
  if (!requireNamespace(orgdb, quietly = TRUE))
    return(list(state = "warn", text = sprintf(
      "Not installed in this R. Install it first: BiocManager::install('%s'), then come back here.",
      orgdb)))
  if (!is.null(gtf_species) && !identical(orgdb, gtf_species$orgdb))
    return(list(state = "warn", text = sprintf(
      "Installed, but the GTF's gene IDs are %s: that species' package is %s.",
      gtf_species$species, gtf_species$orgdb)))
  list(state = "ok", text = paste0("Installed",
                                   if (!is.null(gtf_species))
                                     paste0(", matches the GTF (",
                                            gtf_species$species, ")"),
                                   "."))
}

# What the wizard says under the genome label: blank stops the first launch of a
# new GTF (chromosome lengths are fetched from UCSC by this name and cached per
# GTF; ensure_genome_ready() falls back to the OrgDb name, which UCSC does not
# know); a label for another build than the GTF's header names (`gtf_build`,
# from wiz_gtf_check()) is the mm10-for-mm39 mistake; with no build in the
# header, a label of the other species (`gtf_species`) is still caught.
wiz_genome_label_check <- function(label, gtf_build = NULL,
                                   gtf_species = NULL) {
  label <- wiz_clean_optional_text(label)
  if (is.null(label))
    return(list(state = "warn", text = paste0(
      "Needed: the chromosome lengths are fetched from UCSC by this name, ",
      "and without it the first launch of a new GTF stops.")))
  expected <- if (!is.null(gtf_build) &&
                  gtf_build %in% names(WIZ_LABEL_BY_BUILD))
    WIZ_LABEL_BY_BUILD[[gtf_build]] else NULL
  if (!is.null(expected) && !identical(label, expected))
    return(list(state = "warn", text = sprintf(
      "The GTF is for %s, whose UCSC name is %s, not %s.",
      gtf_build, expected, label)))
  if (!is.null(expected))
    return(list(state = "ok",
                text = sprintf("Matches the GTF's build (%s).", gtf_build)))
  label_species <- c(mm39 = "mouse", mm10 = "mouse", hg38 = "human",
                     hg19 = "human")[label]
  if (!is.null(gtf_species) && !is.na(label_species) &&
      !identical(label_species[[1]], gtf_species$species))
    return(list(state = "warn", text = sprintf(
      "%s is a %s build, but the GTF's gene IDs are %s.",
      label, label_species[[1]], gtf_species$species)))
  list(state = "blank", text = NULL)
}
