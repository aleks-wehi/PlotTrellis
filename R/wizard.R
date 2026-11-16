# The setup wizard: a Shiny app of its own (wizard_ui()/wizard_server()) for a first-time
# setup and for editing an existing project's settings, with nothing from the composer
# loaded. Save writes project_settings.json; Save and Close then stop it with "restart",
# and run_plottrellis()'s loop starts the composer. The samples live in
# project_settings.json; a samples.txt is only imported or exported here.

WIZARD_DEFAULT_COLOUR <- "#1f78b4"  # PALETTE[1] (figure_helpers.R); the same for every new sample, not cycled

WIZARD_DATA_TYPES <- list(
  bw          = list(label = "BigWig coverage (ATAC / ChIP-seq / CUT&RUN / RNA-seq)", ext = "\\.(bw|bigWig)$"),
  hic         = list(label = "Hi-C contact map",                                     ext = "\\.(hic|cool|mcool)$"),
  bed         = list(label = "BED intervals",                                       ext = "\\.(bed|narrowPeak|broadPeak|gappedPeak)(\\.gz)?$"),
  `paired-end` = list(label = "Paired-interval arcs (BED12 or BED6 splice junctions, or BEDPE or HiCCUPS loop calls)",
                       ext  = "\\.(bed|bed12|bedpe)(\\.gz)?$")
)
# What the Add-sample dialog calls each kind of file (the type keys stay as they are).
WIZ_KIND_LABELS <- c(bw = "Coverage (bigWig)", hic = "Hi-C map", bed = "Intervals (BED)",
                     `paired-end` = "Arcs (junctions, loops)")

# Returns paths relative to data_root (for populating an instance's path selectize); quantdir instead lists directories containing quant.sf, not extension-matched files.
wizard_scoped_choices <- function(data_root, type) {
  if (is.null(data_root) || !nchar(data_root) || !dir.exists(data_root)) return(character(0))
  if (identical(type, "quantdir")) {
    quant_files <- list.files(data_root, pattern = "^quant\\.sf$", recursive = TRUE, full.names = TRUE)
    dirs <- unique(dirname(quant_files))
    return(substring(dirs, nchar(data_root) + 2))
  }
  ext <- WIZARD_DATA_TYPES[[type]]$ext
  files <- list.files(data_root, pattern = ext, recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  substring(files, nchar(data_root) + 2)
}

# A literal "\u00b7" not surrounded by spaces breaks format_track_label()'s bolding regex -- reject at entry time.
wizard_name_has_bad_dot <- function(name) {
  grepl("\u00b7", name, fixed = TRUE) && !grepl("\\s\u00b7\\s", name)
}

# Same idea as clean_optional_numeric(): normalizes a blank input to NULL, not "", matching how an unset path is represented in the settings (CFG) and avoiding a spurious is_dirty()/export delta trip against NULL.
wiz_clean_optional_text <- function(v) {
  v <- trimws(v %||% "")
  if (nchar(v) == 0) NULL else v
}

# The message for a blank genome label, shared by Save/import (wiz_save_problems()) and
# opening a project (.pt_require_genome_label()).
PT_GENOME_LABEL_REQUIRED <- "Genome label is required: the UCSC name of the genome build, e.g. mm39, mm10 or hg38."

# Stops with PT_GENOME_LABEL_REQUIRED when a project's settings have a blank genome label (a
# hand-edited project_settings.json can: the wizard's Save and the samples.txt import refuse
# one), before the first launch of a new GTF fails on it with a console error.
.pt_require_genome_label <- function(genome_label, settings_json) {
  if (is.null(wiz_clean_optional_text(genome_label)))
    stop(PT_GENOME_LABEL_REQUIRED, " Set \"genome_label\" in ", settings_json,
         ", or fill it in again with plottrellis_wizard().", call. = FALSE)
}

# The checks the wizard's Save runs before writing project_settings.json, shared with
# parse_samples_txt() so the two can't drift. Returns the problems, as
# the wizard's own messages; character(0) means OK. An OrgDb package that isn't installed
# otherwise only fails at the first cold start, after Save restarts the app, with the hint
# in the console only. A genome_label that GenomeInfoDb can't look chromosome sizes
# up by (e.g. "m10") otherwise only fails at the first cold start, after Save; so
# does a blank one (the first launch of a new GTF stops), so both are problems.
wiz_save_problems <- function(data_root, samples, gtf_path, orgdb, genome_label = NULL) {
  errs <- character(0)
  if (is.null(wiz_clean_optional_text(data_root))) errs <- c(errs, "Data root is required.")
  else if (!dir.exists(data_root)) errs <- c(errs, "Data root directory does not exist.")
  if (length(samples) == 0) errs <- c(errs, "Add at least one sample before saving.")
  if (is.null(wiz_clean_optional_text(gtf_path))) errs <- c(errs, "GTF path is required.")
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
    errs <- c(errs, sprintf("Genome label '%s' isn't a known genome build id (e.g. mm10, mm39, hg38).", label))
  errs
}

# wiz_help() (the help button beside a label) lives in help_popover.R.

# The species PlotTrellis supports, by Ensembl gene-ID prefix, with the OrgDb package each
# needs: what the wizard's live GTF and OrgDb checks compare against. Mouse and human only
# (ORGDB_TAXON_ID, the UniProt links, map only these two).
WIZ_SPECIES_BY_ID_PREFIX <- list(
  ENSG    = list(species = "human", orgdb = "org.Hs.eg.db"),
  ENSMUSG = list(species = "mouse", orgdb = "org.Mm.eg.db")
)

# The UCSC genome label for each build a GENCODE "##description:" or Ensembl "#!genome-build"
# header line can name, for the live genome-label check.
WIZ_LABEL_BY_BUILD <- c(GRCm39 = "mm39", GRCm38 = "mm10", GRCh38 = "hg38", GRCh37 = "hg19")

# What the wizard says under the GTF path as it is typed, before Save. Reads only the file's
# first lines (gzipped or not), so it is cheap on a 1 GB file. A warning here never blocks
# Save: ensure_genome_ready() stops with the full message at the first launch. Returns
# list(state = "blank" | "missing" | "warn" | "ok", text, species = a WIZ_SPECIES_BY_ID_PREFIX
# entry or NULL, build = the genome build named in the header, e.g. "GRCm38", or NULL).
wiz_gtf_check <- function(path, n_lines = 2000) {
  res <- function(state, text = NULL, species = NULL, build = NULL)
    list(state = state, text = text, species = species, build = build)
  path <- wiz_clean_optional_text(path)
  if (is.null(path)) return(res("blank"))
  if (!is_real_file(path))
    return(res("missing", "No file at this path (it must be a path on the computer running R)."))
  if (!grepl("\\.gtf(\\.gz)?$", path, ignore.case = TRUE))
    return(res("warn", "Expected a .gtf or .gtf.gz file. A GFF3 (.gff, .gff3) or other annotation format is not read."))
  # A corrupt .gz warns as well as errors; neither reaches the console.
  lines <- tryCatch({
    con <- gzfile(path)
    on.exit(close(con), add = TRUE)
    suppressWarnings(readLines(con, n = n_lines, warn = FALSE))
  }, error = function(e) character(0))
  header <- lines[startsWith(lines, "#")]
  body   <- lines[!startsWith(lines, "#") & nzchar(lines)]
  # The gene_id attribute itself, not another ending in gene_id (StringTie's ref_gene_id).
  id_re <- '(^|[; \t])gene_id "([^"]*)"'
  ids <- regmatches(body, regexpr(id_re, body))
  ids <- sub(id_re, "\\2", ids)
  if (length(ids) == 0)
    return(res("warn", sprintf("No gene_id attribute in the first %d lines. Is this a GTF file?", n_lines)))
  if (!any(startsWith(ids, "ENS")))
    return(res("warn", sprintf(paste0(
      "The gene IDs here look like '%s', not Ensembl IDs (ENSG..., ENSMUSG...). This is ",
      "probably a RefSeq/NCBI GTF, which PlotTrellis cannot use: gene search joins the GTF ",
      "to the OrgDb package on Ensembl IDs. Download the GENCODE or Ensembl GTF for the ",
      "same genome build instead (see 'Which GTF?' below)."), ids[1])))
  ens_ids <- ids[startsWith(ids, "ENS")]
  # GENCODE's GRCh37-mapped ("lift37") files suffix every ID with _<n>, which no version
  # stripping removes, so no gene name ever matches.
  if (any(grepl("_[0-9]+$", ens_ids)))
    return(res("warn", sprintf(paste0(
      "Gene IDs like '%s' are from a GENCODE \"lift37\" GTF (a newer release mapped back to ",
      "GRCh37). PlotTrellis cannot match them to gene names: use the GENCODE 19 GTF instead."), ens_ids[1])))
  prefix  <- sub("[0-9].*$", "", ens_ids[1])
  species <- WIZ_SPECIES_BY_ID_PREFIX[[prefix]]
  # GENCODE's "##description:" line names the genome and release ("mouse genome (GRCm38),
  # version M25 (Ensembl 100)"); Ensembl's "#!genome-build" names the build.
  desc <- sub("^##description: (evidence-based annotation of the )?", "", grep("^##description: ", header, value = TRUE))
  ens_build <- sub("^#!genome-build ", "", grep("^#!genome-build ", header, value = TRUE))
  build <- unlist(regmatches(c(desc, ens_build), gregexpr("GRC[mh]3[0-9]", c(desc, ens_build))))
  build <- if (length(build)) build[1] else NULL
  sp <- if (!is.null(species)) paste0(" (", species$species, ")") else ""
  src <- if (any(grepl("^##provider: GENCODE", header))) {
    paste0("a GENCODE GTF", if (length(desc)) paste0(": ", desc[1]) else sp)
  } else if (length(ens_build)) {
    paste0("an Ensembl GTF, genome build ", ens_build[1], sp)
  } else paste0("a GTF with Ensembl gene IDs", sp)
  if (is.null(species))
    return(res("warn", paste0("Looks like ", src, ", but the gene IDs (", prefix, "...) are not mouse or human. ",
                              "PlotTrellis supports mouse and human only."), build = build))
  res("ok", paste0("Looks like ", src, ". Make sure it is the release your reads were quantified with."),
      species = species, build = build)
}

# What the wizard says under the OrgDb field: whether the package is installed, and whether
# it is the one the GTF's gene IDs call for (`gtf_species`, from wiz_gtf_check()).
wiz_orgdb_check <- function(orgdb, gtf_species = NULL) {
  orgdb <- wiz_clean_optional_text(orgdb)
  if (is.null(orgdb)) return(list(state = "blank", text = NULL))
  if (!requireNamespace(orgdb, quietly = TRUE))
    return(list(state = "warn", text = sprintf(
      "Not installed in this R. Install it first: BiocManager::install('%s'), then come back here.", orgdb)))
  if (!is.null(gtf_species) && !identical(orgdb, gtf_species$orgdb))
    return(list(state = "warn", text = sprintf(
      "Installed, but the GTF's gene IDs are %s: that species' package is %s.", gtf_species$species, gtf_species$orgdb)))
  list(state = "ok", text = paste0("Installed", if (!is.null(gtf_species)) paste0(", matches the GTF (", gtf_species$species, ")"), "."))
}

# What the wizard says under the genome label: blank stops the first launch of a new GTF
# (chromosome lengths are fetched from UCSC by this name and cached per GTF;
# ensure_genome_ready() falls back to the OrgDb name, which UCSC does not know); a label for
# another build than the GTF's header names (`gtf_build`, from wiz_gtf_check()) is the
# mm10-for-mm39 mistake; with no build in the header, a label of the other species
# (`gtf_species`) is still caught.
wiz_genome_label_check <- function(label, gtf_build = NULL, gtf_species = NULL) {
  label <- wiz_clean_optional_text(label)
  if (is.null(label))
    return(list(state = "warn", text = paste0("Needed: the chromosome lengths are fetched from UCSC by this name, ",
                                              "and without it the first launch of a new GTF stops.")))
  expected <- if (!is.null(gtf_build) && gtf_build %in% names(WIZ_LABEL_BY_BUILD)) WIZ_LABEL_BY_BUILD[[gtf_build]] else NULL
  if (!is.null(expected) && !identical(label, expected))
    return(list(state = "warn", text = sprintf("The GTF is for %s, whose UCSC name is %s, not %s.", gtf_build, expected, label)))
  if (!is.null(expected)) return(list(state = "ok", text = sprintf("Matches the GTF's build (%s).", gtf_build)))
  label_species <- c(mm39 = "mouse", mm10 = "mouse", hg38 = "human", hg19 = "human")[label]
  if (!is.null(gtf_species) && !is.na(label_species) && !identical(label_species[[1]], gtf_species$species))
    return(list(state = "warn", text = sprintf("%s is a %s build, but the GTF's gene IDs are %s.",
                                               label, label_species[[1]], gtf_species$species)))
  list(state = "blank", text = NULL)
}

# The messages a failed Save or import shows: one bullet each, never one joined line.
# `kind` "error" is a red alert; "warn" an amber note (an import that succeeded with warnings).
wiz_error_list <- function(heading, msgs, kind = "error") {
  div(class = if (identical(kind, "warn")) "wiz-warn wiz-warn-box" else "wiz-error wiz-error-box",
      role = if (identical(kind, "warn")) "status" else "alert",
      if (!is.null(heading)) tags$b(heading),
      tags$ul(style = "margin:4px 0 0 0; padding-left:20px;", lapply(msgs, tags$li)))
}

# A colour swatch with no shiny.js binding (type="color"): oninput reports changes, and the
# inline script reports the shown value on render, so a re-rendered swatch (the
# composer's Revert) never leaves input[[id]] at a stale colour.
.pt_colour_input <- function(id, label, title, value) {
  tagList(tags$label(label, title = title), br(),
          tags$input(type = "color", id = id, value = value,
                     style = "width:60px; height:34px; padding:2px; border:1px solid #ccc; border-radius:4px; cursor:pointer; display:block;",
                     oninput = sprintf("Shiny.setInputValue('%s', this.value, {priority: 'event'});", id)),
          tags$script(sprintf("Shiny.setInputValue('%s', %s);", id, jsonlite::toJSON(value, auto_unbox = TRUE))))
}

# Core settings: what identifies and locates the project's data. Set only in the wizard;
# every other .pt_default_settings() key is a display setting, adjustable live in the
# composer.
.PT_CORE_SETTINGS <- c("gtf_path", "orgdb", "genome_label", "data_dir", "samples")

# Inputs for every display setting, ids paste0(prefix, key), initial values from cfg: the
# composer's Figure settings panel.
# Returns list(display =, diff =): the general settings and the 2-sample comparison ones.
.pt_display_settings_ui <- function(cfg, prefix) {
  list(
    display = tagList(
      .pt_settings_group("Layout & display",
      fluidRow(
        column(3, selectizeInput(paste0(prefix, "font"), "Font:", choices = c("DejaVu Sans", "Helvetica", "Arial"), selected = cfg$font %||% "DejaVu Sans")),
        column(3, numericInput(paste0(prefix, "fontsize"), "Font size:", value = cfg$fontsize %||% 8, min = 4, max = 24)),
        column(3, numericInput(paste0(prefix, "page_w"),
                                wiz_help("Page width (in):", "Width of the exported figure, in inches. Track heights below adapt to whatever's left after labels and margins."),
                                value = cfg$page_w %||% 9, min = 3, max = 30)),
        column(3, selectizeInput(paste0(prefix, "chr_fmt"),
                               wiz_help("Chromosome format:", "How chromosome names are displayed on the figure's axis."),
                               choices = c("Dot (chr.19)" = "dot", "Space (chr 19)" = "space", "None (chr19)" = "plain"),
                               selected = cfg$chr_fmt %||% "dot")),
        column(3, checkboxInput(paste0(prefix, "show_track_axes"),
                                 wiz_help("Show track value axes", "Draws a numeric value axis on coverage (bigWig) and junction tracks. On by default."),
                                 value = cfg$show_track_axes %||% TRUE))
      ),
      fluidRow(
        column(3, checkboxInput(paste0(prefix, "show_title"),
                                 wiz_help("Show title", "Draws the title (the Title field above) over the figure; with no title text, nothing is drawn. On by default."),
                                 value = cfg$show_title %||% TRUE)),
        column(3, checkboxInput(paste0(prefix, "show_scale_bar"),
                                 wiz_help("Show scale bar", "Draws the genomic distance scale bar on the first coverage-type track. On by default."),
                                 value = cfg$show_scale_bar %||% TRUE))
      )),

      .pt_settings_group("Coverage tracks (BigWig: ATAC / ChIP-seq / CUT&RUN / RNA-seq)",
      fluidRow(
        column(3, numericInput(paste0(prefix, "atac_height"),
                                wiz_help("Track height (in):", "Height in inches of each BigWig coverage track (ATAC, ChIP-seq, CUT&RUN, RNA-seq)."),
                                value = cfg$atac_height %||% 0.85, min = 0.3, max = 8, step = 0.1))
      )),

      .pt_settings_group("Hi-C",
      fluidRow(
        column(3, selectizeInput(paste0(prefix, "hic_style"),
                               wiz_help("Display style:", "How the Hi-C contact matrix is drawn: a full rectangle, or a triangle (upper-left half)."),
                               choices = c("Rectangle" = "rectangle", "Triangle" = "triangle"), selected = cfg$hic_style %||% "rectangle")),
        column(3, textInput(paste0(prefix, "hic_norm"),
                             wiz_help("Normalisation:", "Hi-C matrix normalisation method understood by the .hic/.cool/.mcool file, e.g. KR, VC, SCALE. Leave as KR unless your data needs a different one."),
                             value = cfg$hic_norm %||% "KR")),
        column(3, numericInput(paste0(prefix, "zrange_cap"),
                                wiz_help("Colour scale cap:", "Caps the Hi-C colour scale at this contact-count value. Leave blank to scale automatically per plot."),
                                value = clean_optional_numeric(cfg$zrange_cap) %||% NA, min = 1, step = 10)),
        column(3, numericInput(paste0(prefix, "hic_height"),
                                wiz_help("Track height (in):", "Height in inches of each Hi-C track. Leave blank to size automatically based on the display style."),
                                value = clean_optional_numeric(cfg$hic_height) %||% NA, min = 0.5, max = 15, step = 0.5))
      ),
      fluidRow(
        column(3, selectizeInput(paste0(prefix, "hic_res_mode"),
                               wiz_help("Resolution mode:", "Auto picks the coarsest resolution that still gives at least the \"Max bins per view\" number of bins across the view; Finest always uses the highest-resolution matrix available (slower for wide views); Fixed uses the one bin size you type (a Hi-C track whose file lacks that bin size shows an error). plotgardener draws at most 1,000 bins across a view, so for a wider view it uses a coarser size instead (\"too many Hi-C pixels\"), in every mode."),
                               choices = c("Auto" = "auto", "Finest" = "finest", "Fixed" = "fixed"), selected = cfg$hic_res_mode %||% "auto")),
        column(3, numericInput(paste0(prefix, "hic_max_bins_per_view"),
                                wiz_help("Max bins per view:", "Auto resolution mode picks the coarsest resolution that still gives at least this many bins across the view (the finest if none does). Despite the name it is a minimum, not an upper limit."),
                                value = cfg$hic_max_bins_per_view %||% 500, min = 50, max = 5000, step = 50)),
        column(3, conditionalPanel(paste0("input.", prefix, "hic_res_mode === 'fixed'"),
                 numericInput(paste0(prefix, "hic_resolution"),
                                wiz_help("Bin size (bp):", "The Hi-C bin size, in base pairs, used for every Hi-C track (e.g. 10000). A track whose file has no bin of exactly this size shows an error naming the closest sizes the file has."),
                                value = clean_optional_numeric(cfg$hic_resolution) %||% NA, min = 0, step = 1000)))
      )),

      .pt_settings_group("Junction / arc tracks",
      fluidRow(
        column(4, checkboxInput(paste0(prefix, "arc_clip_partial"),
                                 wiz_help("Clip partial arcs", "For an arc that runs past the edge of the plotted region (one anchor inside it, the other outside): on, the arc is left out; off, it is drawn, cut off at the edge. Off by default. Long-range BEDPE-sourced pairs (e.g. Hi-C loop calls) run past the edge more often than short RNA-seq splice junctions."),
                                 value = cfg$arc_clip_partial %||% FALSE)),
        column(3, numericInput(paste0(prefix, "min_junction_support"),
                                wiz_help("Min arc support:", "Minimum pooled nSupport for a paired-end arc (BED12 or BED6 splice junction, or BEDPE or HiCCUPS loop call) to be drawn. Default 0 (show everything) -- a BEDPE score column's scale is caller-specific (e.g. a HiChIP loop caller's continuous significance score, not a read count), so raise this only once you know your own data's scale."),
                                value = cfg$min_junction_support %||% 0, min = 0, step = 1)),
        column(3, numericInput(paste0(prefix, "junction_height"),
                                wiz_help("Junction track height (in):", "Height in inches of a junction track's arcs, below its heading. A mirrored comparison draws two samples, so its arcs are twice this tall."),
                                value = cfg$junction_height %||% 0.85, min = 0.3, max = 8, step = 0.1))
      )),

      .pt_settings_group("Gene view",
      fluidRow(
        column(3, numericInput(paste0(prefix, "gene_flank_pct"),
                                wiz_help("Flanking (%):", "Percent of the gene's length added to each side when jumping to a gene, so nearby regulatory regions stay in view."),
                                value = cfg$gene_flank_pct %||% 15, min = 0, max = 200, step = 5)),
        column(3, numericInput(paste0(prefix, "min_display_tpm"),
                                wiz_help("Min display TPM:", "Per-condition mean-TPM floor a transcript must clear to be individually named/displayed; below it, a transcript folds into 'Other'. Default 1."),
                                value = cfg$min_display_tpm %||% 1, min = 0, step = 0.5)),
        column(4, checkboxInput(paste0(prefix, "uniprot_reviewed_only"),
                                 wiz_help("UniProt: reviewed only", "Restrict the Gene info panel's UniProt link-out to manually reviewed (Swiss-Prot) entries only."),
                                 value = cfg$uniprot_reviewed_only %||% FALSE))
      ))
    ),
    diff = tagList(
      fluidRow(
        column(3, .pt_colour_input(paste0(prefix, "diff_colour_lo"), "Loss colour (coverage diff):", "Colour for coverage lost in the comparison condition, in the coverage difference track.", cfg$diff_colour_lo %||% "#0033CC")),
        column(3, .pt_colour_input(paste0(prefix, "diff_colour_hi"), "Gain colour (coverage diff):", "Colour for coverage gained in the comparison condition, in the coverage difference track.", cfg$diff_colour_hi %||% "#CC0000")),
        column(3, numericInput(paste0(prefix, "junction_opacity_cap"),
                                wiz_help("Junctions opacity cap (log2FC):", "log2 fold-change value at which a junction arc reaches full opacity in a comparison track. Smaller values make modest changes more visible."),
                                value = cfg$junction_opacity_cap %||% 2, min = 0.1, step = 0.5)),
        column(3, numericInput(paste0(prefix, "diff_shared_alpha"),
                                wiz_help("Shared coverage tint:", "How strongly the coverage both conditions share is tinted in the coverage difference track: 0 = white, 1 = the full reference colour."),
                                value = cfg$diff_shared_alpha %||% 0.15, min = 0, max = 1, step = 0.05))
      ),
      fluidRow(
        column(3, .pt_colour_input(paste0(prefix, "reference_colour"), "Reference colour:", "Fixed colour always used for the reference condition in a 2-sample comparison track, regardless of that condition's own colour.", cfg$reference_colour %||% "#4D4D4D")),
        column(3, .pt_colour_input(paste0(prefix, "comparison_colour"), "Comparison colour:", "Fixed colour always used for the comparison condition in a 2-sample comparison track, regardless of that condition's own colour.", cfg$comparison_colour %||% "#56B4E9"))
      )
    )
  )
}

# Reads the display settings back from .pt_display_settings_ui()'s inputs. An input that
# isn't there (NULL, e.g. not rendered yet) or is blank falls back to cfg, then to the
# default, so it never silently resets a value. The optional numbers (blank = auto) only
# fall back when the input is missing. names() of the result are the display-setting keys.
.pt_collect_display_settings <- function(input, prefix, cfg) {
  dflt <- .pt_default_settings()
  get <- function(k) input[[paste0(prefix, k)]] %||% cfg[[k]] %||% dflt[[k]]
  opt <- function(k) {
    v <- input[[paste0(prefix, k)]]
    clean_optional_numeric(if (is.null(v)) cfg[[k]] else v)
  }
  list(
    font        = trimws(get("font")),
    fontsize    = get("fontsize"),
    page_w      = get("page_w"),
    chr_fmt     = get("chr_fmt"),
    show_track_axes = get("show_track_axes"),
    show_title      = get("show_title"),
    show_scale_bar  = get("show_scale_bar"),
    hic_norm    = trimws(get("hic_norm")),
    hic_style   = get("hic_style"),
    zrange_cap  = opt("zrange_cap"),
    atac_height = get("atac_height"),
    hic_height  = opt("hic_height"),
    hic_res_mode = get("hic_res_mode"),
    hic_max_bins_per_view = get("hic_max_bins_per_view"),
    hic_resolution = opt("hic_resolution"),
    diff_colour_lo    = trimws(get("diff_colour_lo")),
    diff_colour_hi    = trimws(get("diff_colour_hi")),
    junction_opacity_cap = get("junction_opacity_cap"),
    min_junction_support = get("min_junction_support"),
    reference_colour  = trimws(get("reference_colour")),
    comparison_colour = trimws(get("comparison_colour")),
    gene_flank_pct    = get("gene_flank_pct"),
    min_display_tpm   = get("min_display_tpm"),
    uniprot_reviewed_only = get("uniprot_reviewed_only"),
    arc_clip_partial  = get("arc_clip_partial"),
    junction_height   = get("junction_height"),
    diff_shared_alpha = get("diff_shared_alpha")
  )
}

# Writes project_settings.json and, beside it, the project_samples.txt copy the project can
# be rebuilt from, with absolute paths since the importer resolves relative ones
# against the file's own folder. Shared by the wizard's Save and the composer's Save as
# project default. Returns whether the samples.txt copy was written; a failure there
# doesn't block.
.pt_write_project_files <- function(settings, settings_json) {
  # null = "null": jsonlite's toJSON() default serializes a NULL list element
  # (e.g. zrange_cap/hic_height left blank for "auto") as {} rather than the
  # JSON literal null. On reload, fromJSON() then reads that back as an empty
  # list, not NULL, and clean_optional_numeric()'s is.finite() crashes on it. digits = NA: every
  # writer of project_settings.json keeps full precision (toJSON() rounds to 4 decimals).
  # Through a temporary file renamed into place, so a failed write (a full disk) never
  # leaves a truncated settings file; the error names the reason (e.g. "Permission
  # denied"), which writeLines() only gives as a warning.
  tmp <- tempfile("project_settings-", tmpdir = dirname(settings_json), fileext = ".json")
  why <- character(0)
  ok <- tryCatch(withCallingHandlers({
    writeLines(jsonlite::toJSON(settings, auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA), tmp)
    file.rename(tmp, settings_json)
  }, warning = function(w) { why <<- c(why, conditionMessage(w)); invokeRestart("muffleWarning") }),
  error = function(e) { why <<- c(why, conditionMessage(e)); FALSE })
  if (!isTRUE(ok)) {
    unlink(tmp)
    stop(paste(unique(why), collapse = "; "), call. = FALSE)
  }
  .pt_write_samples_copy(settings, settings_json)
}

# The project_samples.txt beside settings_json, from `settings`; whether it was written. A
# failure (error or warning) only gives FALSE.
.pt_write_samples_copy <- function(settings, settings_json) {
  tryCatch(withCallingHandlers({
    flat <- settings[setdiff(names(settings), c("data_dir", "samples"))]
    if (length(flat$gtf_path) == 1 && nzchar(flat$gtf_path)) flat$gtf_path <- normalizePath(flat$gtf_path, mustWork = FALSE)
    txt <- export_samples_txt(c(list(data_root = normalizePath(settings$data_dir, mustWork = FALSE),
                                     samples = settings$samples), flat),
                              defaults = .pt_default_settings())
    # The start view as start lines after the other "# key:" lines, so the copy restores it too.
    if (is.list(settings$start_view))
      txt <- append(txt, tryCatch(.pt_start_view_lines(settings$start_view, settings), error = function(e) character(0)),
                    after = match(FALSE, startsWith(txt, "#"), nomatch = length(txt) + 1) - 1)
    writeLines(txt, file.path(dirname(settings_json), "project_samples.txt"))
    TRUE
  }, warning = function(w) stop(conditionMessage(w), call. = FALSE)), error = function(e) FALSE)
}

# The start view (pt_preload(), or a samples.txt's start_gene/start_tracks) of the saved
# settings `old`, if the wizard's new settings `new` still have everything it refers to: the
# same GTF and OrgDb (its gene is looked up in them), each of its samples, and each of its
# tracks (same key, i.e. same sample and data name; a changed file path or data folder keeps
# it). Otherwise NULL: the app then opens blank rather than on a view the user changed
# underneath it. Other samples can come and go.
.pt_kept_start_view <- function(old, new) {
  sv <- old$start_view
  if (!is.list(sv)) return(NULL)
  path <- function(p) if (is.character(p) && length(p) == 1 && nzchar(p)) normalizePath(p, mustWork = FALSE) else ""
  if (!identical(path(old$gtf_path), path(new$gtf_path)) || !identical(old$orgdb, new$orgdb)) return(NULL)
  samples <- load_samples_from_settings(new$samples, new$data_dir)
  conds <- vapply(samples, function(s) s$display, character(1))
  keys <- names(c(build_track_items(samples), build_comparison_track_items(samples)))
  if (!all(as.character(unlist(sv$conditions)) %in% conds) ||
      !all(.pt_renamed_junctions_key(as.character(unlist(sv$tracks))) %in% keys)) return(NULL)
  sv
}

# One numbered section card of the setup form (Genome, Data folder, Samples).
wiz_section <- function(number, title, intro, ...) {
  div(class = "wiz-section",
      div(class = "wiz-section-head",
          tags$span(class = "wiz-step", number),
          tags$h3(title)),
      if (!is.null(intro)) p(class = "wiz-section-intro", intro),
      ...)
}

# The "Which GTF?" guide under the GTF field: what works and why, where to get it, and the
# OrgDb package and genome label that go with each. Every claim here is checked against
# ensure_genome_ready() (Ensembl IDs required, UCSC chromosome names applied, lengths fetched
# by genome_label); a plain Ensembl GTF was run end to end before the Ensembl claims were written.
wiz_gtf_guide <- function() {
  row <- function(...) tags$tr(lapply(list(...), tags$td))
  tags$details(class = "wiz-guide",
    tags$summary("Which GTF? Where do I get it?"),
    tags$ul(
      tags$li(tags$b("GENCODE or Ensembl, not RefSeq. "),
              "PlotTrellis finds a gene by looking its name up in the OrgDb package, which gives the ",
              "gene's Ensembl ID (ENSG..., ENSMUSG...), then finds that ID in the GTF. GENCODE and ",
              "Ensembl GTFs use Ensembl IDs; RefSeq (NCBI) GTFs use other IDs, so no gene would be ",
              "found. GENCODE gives the most: an Ensembl GTF draws every track, but has no APPRIS labels ",
              "for the Gene info table. Mouse and human are supported."),
      tags$li(tags$b("The release your quantification used. "),
              "Transcripts in your quant.sf files are matched to the GTF by ID, so use the GTF from ",
              "the same GENCODE or Ensembl release as the transcript sequences Salmon was run on; a ",
              "transcript the GTF does not have is left out of the expression tracks. For GENCODE, ",
              "Salmon must have been run with --gencode, or no transcript name matches. A different ",
              "genome build (mm10 vs mm39) puts every gene at the wrong coordinates."),
      tags$li(tags$b("Keep it gzipped "), "(.gtf.gz). The first launch builds an index from it, ",
              "about 1 to 3 minutes for a whole genome, once per GTF file.")
    ),
    p(tags$b("Where to download it: "),
      tags$a(href = "https://www.gencodegenes.org/human/", target = "_blank", "GENCODE human"), " or ",
      tags$a(href = "https://www.gencodegenes.org/mouse/", target = "_blank", "GENCODE mouse"),
      ", the \"Comprehensive gene annotation\" GTF (region CHR); earlier releases under \"Releases\" on the same page. ",
      "Ensembl: ", tags$code("https://ftp.ensembl.org/pub/release-<N>/gtf/<species>/"), ", the file ",
      tags$code("<Species>.<build>.<N>.gtf.gz"), "; human GRCh37 is under ", tags$code("pub/grch37/"), " instead."),
    tags$table(class = "table table-condensed wiz-guide-table",
      tags$thead(tags$tr(tags$th("Species, build"), tags$th("GENCODE releases"), tags$th("Ensembl releases"), tags$th("OrgDb package"), tags$th("Genome label"))),
      tags$tbody(
        row("Mouse, GRCm39", "M26 and later", "103 and later", tags$code("org.Mm.eg.db"), tags$code("mm39")),
        row("Mouse, GRCm38", "M25 and earlier", "68 to 102", tags$code("org.Mm.eg.db"), tags$code("mm10")),
        row("Human, GRCh38", "20 and later", "76 and later", tags$code("org.Hs.eg.db"), tags$code("hg38")),
        row("Human, GRCh37", "3c to 19 (not the \"lift37\" files)", "55 to 75, or the GRCh37 archive", tags$code("org.Hs.eg.db"), tags$code("hg19"))
      )
    ),
    p(class = "text-muted", style = "margin-bottom:0;",
      "The build must be the one your reads were aligned to. GENCODE's \"lift37\" GTFs carry IDs such as ",
      tags$code("ENSG00000290825.1_1"), " that match no gene name: use GENCODE 19 instead. ",
      "Install the OrgDb package once in R: ", tags$code("BiocManager::install(\"org.Mm.eg.db\")"),
      " (mouse) or ", tags$code("BiocManager::install(\"org.Hs.eg.db\")"), " (human).")
  )
}

# Shared Genome/Data folder/Samples/Save form for every entry path (import, manual, edit-mode re-entry; import only prefills it); `cfg` supplies initial values (CFG when editing, DEFAULT_SETTINGS when fresh), `is_edit` governs what else differs (data-root's blank default, Save's label, the extra Close button).
wiz_shared_form <- function(cfg, is_edit) {
  tagList(
    wiz_section(1, "Genome",
      "The gene annotation (GTF) PlotTrellis draws genes and transcripts from, and the R package it looks gene names up in.",
      fluidRow(
        column(5, selectizeInput("wiz_genome_preset",
                    wiz_help("Species and genome build",
                             "Fills in the OrgDb package and genome label for a common genome; for anything else type them yourself below. The GTF is always your own file."),
                    choices = c("Choose, or fill the fields below yourself" = "",
                                setNames(names(GENOME_PRESETS), vapply(GENOME_PRESETS, function(g) g$label, character(1)))),
                    width = "100%"))
      ),
      textInput("wiz_gtf_path",
                wiz_help("GTF file (path on the computer running R)",
                         "Your genome's gene annotation, from GENCODE or Ensembl (mouse or human), gzipped or not. It must be the same genome build and release your reads were aligned or quantified with. RefSeq/NCBI GTFs do not work here. Open 'Which GTF?' below for where to download it."),
                value = cfg$gtf_path %||% "", width = "100%",
                placeholder = "e.g. /data/reference/gencode.vM32.annotation.gtf.gz"),
      uiOutput("wiz_gtf_status"),
      wiz_gtf_guide(),
      fluidRow(
        column(6, textInput("wiz_orgdb",
                            wiz_help("OrgDb package",
                                     "The Bioconductor annotation package for your species' gene names: org.Mm.eg.db for mouse, org.Hs.eg.db for human. It must already be installed in this R (BiocManager::install('org.Mm.eg.db')). Gene search uses it to turn a gene name into the GTF's gene ID."),
                            value = cfg$orgdb %||% "", width = "100%", placeholder = "org.Mm.eg.db"),
               uiOutput("wiz_orgdb_status")),
        column(6, textInput("wiz_genome_label",
                            wiz_help("Genome label (UCSC name)",
                                     "The genome build's UCSC name: mm39, mm10, hg38 or hg19. PlotTrellis fetches the chromosome lengths from UCSC by this name (a GTF does not contain them) and prints it on the figure. Needed: without it the first launch of a new GTF stops. It must be the build of the GTF (mm10 is GRCm38, mm39 is GRCm39)."),
                            value = cfg$genome_label %||% "", width = "100%", placeholder = "mm39"),
               uiOutput("wiz_genome_label_status"))
      )
    ),
    wiz_section(2, "Data folder",
      "The folder your data files are in (bigWig, Hi-C, BED, junction or loop files, salmon quant folders), on the computer running R. Each sample's files are then picked from it, so one folder for everything is easiest.",
      textInput("wiz_data_root",
                wiz_help("Data folder (full path)",
                         "Files in any subfolder of it are found. If your data is spread over several places, a folder of symbolic links to them works too."),
                value = if (is_edit) cfg$data_dir else "", width = "100%",
                placeholder = "e.g. /data/smithlab/tcell_2026"),
      uiOutput("wiz_data_root_status")
    ),
    wiz_section(3, "Samples",
      "A sample is one condition you want to look at (Control, KO, Basal, LP...), with its files. Add one per condition; the figure then draws them side by side or compares two of them.",
      actionButton("wiz_add_sample", "+ Add a sample", class = "pt-once btn-primary"),
      uiOutput("wiz_samples_ui")
    ),
    p(class = "text-muted",
      "Fonts, track heights and other display settings are changed in the app's Figure settings panel, not here."),
    uiOutput("wiz_completion_error"),
    div(class = "wiz-actions",
      if (is_edit) tagList(
        actionButton("wiz_save", "Save and return to the app", class = "pt-once btn-success btn-lg"),
        actionButton("wiz_close", "Close without saving", class = "pt-once btn-default btn-lg")
      ) else actionButton("wiz_save", "Save and open the app", class = "pt-once btn-success btn-lg"),
      # Export only, never import, and never triggered by Save: kept apart from the primary
      # buttons and captioned so it cannot be mistaken for a way to load samples.
      div(class = "wiz-export",
          actionButton("wiz_export_samples", "Export these settings as a samples.txt", class = "pt-once btn-link btn-sm"),
          div(class = "text-muted", "A text copy you can keep, share, or import next time. Export only."))
    )
  )
}

# The wizard's page for the project whose settings file is (or will be) `SETTINGS_JSON`.
wizard_ui <- function(SETTINGS_JSON) {
  is_edit <- file.exists(SETTINGS_JSON)
  fluidPage(class = "wiz-page",
    tags$head(tags$title(if (is_edit) "PlotTrellis: project settings" else "PlotTrellis: set up your project")),
    .pt_theme_css(),
    .pt_once_js(),
    .pt_help_assets(),
    tags$style(HTML("
      .wiz-page { max-width: 1100px; margin: 0 auto; padding: 28px 24px 60px; font-size: 15px; }
      .wiz-page h1 { font-size: 34px; margin: 0 0 6px; }
      .wiz-page h3 { font-size: 20px; margin: 0; }
      .wiz-lead { font-size: 17px; color: var(--pt-muted); max-width: 820px; margin-bottom: 26px; }
      .wiz-page label { font-size: 15px; }
      /* Landing: what a project needs, then the two ways in. */
      .wiz-needs { display: flex; gap: 14px; flex-wrap: wrap; margin-bottom: 26px; }
      .wiz-need { flex: 1 1 240px; background: var(--pt-surface); border: 1px solid var(--pt-border); border-radius: 12px; padding: 14px 16px; }
      .wiz-need b { display: block; color: var(--pt-heading); margin-bottom: 4px; }
      .wiz-need, .wiz-choice p { color: var(--pt-text); }
      .wiz-choices { display: flex; gap: 16px; flex-wrap: wrap; }
      .wiz-choice { flex: 1 1 300px; background: var(--pt-surface); border: 1px solid var(--pt-border); border-radius: 12px; padding: 20px 20px 16px; box-shadow: var(--pt-shadow); }
      .wiz-choice h3 { margin-bottom: 6px; }
      .wiz-choice p { min-height: 44px; color: var(--pt-muted); }
      .wiz-choice .btn { width: 100%; }
      /* Numbered section cards. */
      .wiz-section { background: var(--pt-surface); border: 1px solid var(--pt-border); border-radius: 12px; box-shadow: var(--pt-shadow); padding: 20px 22px 14px; margin-bottom: 18px; }
      .wiz-section-head { display: flex; align-items: center; gap: 12px; margin-bottom: 6px; }
      .wiz-step { display: inline-flex; align-items: center; justify-content: center; width: 30px; height: 30px; border-radius: 50%; background: var(--pt-primary); color: var(--pt-on-primary); font-weight: 700; font-size: 15px; flex: none; }
      .wiz-section-intro { color: var(--pt-muted); margin-bottom: 14px; }
      .wiz-guide { margin: -4px 0 16px; }
      .wiz-guide > summary { display: list-item; cursor: pointer; font-weight: 600; color: var(--pt-link); }
      .wiz-guide > summary:hover { color: var(--pt-link-hover); }
      .wiz-guide ul { margin-top: 10px; }
      .wiz-guide li { margin-bottom: 6px; }
      .wiz-guide-table { margin-bottom: 8px; }
      .wiz-guide-table th { color: var(--pt-muted); font-weight: 600; }
      .wiz-guide-table td, .wiz-guide-table th { border-top-color: var(--pt-border) !important; }
      /* Field status lines and errors. */
      .wiz-status { margin: -8px 0 12px; font-size: 14px; }
      .wiz-row { padding: 6px 0; border-bottom: 1px solid var(--pt-border); }
      .wiz-row-error { background: var(--pt-danger-soft); }
      .wiz-error { color: var(--pt-danger); margin-bottom: 8px; }
      .wiz-error-box { background: var(--pt-danger-soft); border: 1px solid var(--pt-danger); border-radius: 10px; padding: 10px 14px; }
      .wiz-warn-box { background: #FBF3E3; border: 1px solid var(--pt-warn-edge); border-radius: 10px; padding: 10px 14px; margin-bottom: 8px; }
      .wiz-ok    { color: #2e7d32; }
      .wiz-warn  { color: var(--pt-accent-strong); }
      .wiz-empty { color: var(--pt-muted); padding: 14px 0 6px; }
      .wiz-samples { margin-top: 12px; }
      .wiz-actions { display: flex; align-items: center; gap: 12px; flex-wrap: wrap; margin-top: 8px; }
      .wiz-export { margin-left: auto; text-align: right; font-size: 12px; }
      .wiz-export .btn.btn-link { padding: 0 0 2px; font-weight: 600; box-shadow: none; border: 0; background: transparent; color: var(--pt-link); margin: 0; text-decoration: underline; }
      .wiz-page code { background: var(--pt-code-bg); color: var(--pt-code-text); }
      .wiz-back { display: inline-block; margin-bottom: 14px; font-weight: 600; }
      /* Add-sample modal (modalDialog size='l') is Bootstrap's fixed 900px
         .modal-lg -- widen it so the instance-row selectize boxes have more
         room. #shiny-modal is the one live modal at a time, so this only
         affects a modal-lg while it's showing. */
      #shiny-modal .modal-lg { width: 1100px; max-width: 95vw; }
      #shiny-modal .modal-body { font-size: 15px; }
      .wiz-file-kinds .btn { margin-right: 6px; margin-bottom: 8px; }
      .wiz-kind { font-weight: 600; color: var(--pt-heading); padding-top: 7px; }
      /* selectize's single-select box doesn't truncate long selected paths by
         default -- they overflow past the visible input. Clip with ellipsis
         instead. */
      .selectize-input, .selectize-input > input, .selectize-input > .item {
        overflow: hidden;
        text-overflow: ellipsis;
        white-space: nowrap;
      }
    ")),
    if (is_edit) tagList(
      h1("Project settings"),
      p(class = "wiz-lead",
        "Change your project's genome, data folder or samples. Save restarts the app with the new settings; Close keeps the old ones.")
    ) else tagList(
      h1("Welcome to PlotTrellis"),
      p(class = "wiz-lead",
        "Let's set up your project: where your data is, which genome annotation to use, and what your samples are. ",
        "It takes a few minutes, and everything can be changed later from the app.")
    ),
    uiOutput("wiz_body")
  )
}
# Builds the wizard's server function for one project: `SETTINGS_JSON` is the settings file it
# writes on Save, `CFG` the project's current settings (shown for editing when that file
# exists), `DEFAULT_SETTINGS` the defaults a first-time setup starts from.
wizard_server <- function(SETTINGS_JSON, CFG, DEFAULT_SETTINGS) function(input, output, session) {
  .pt_once_server(input, session)  # repeat-click gate for the "pt-once" buttons (R/once_gate.R)

  is_edit <- file.exists(SETTINGS_JSON)
  cfg     <- if (is_edit) CFG else DEFAULT_SETTINGS

  wiz_samples     <- reactiveVal(if (is_edit) cfg$samples else list())  # committed samples: list(condition=,colour=,data=list(list(type=,name=,path=)))
  modal_instances <- reactiveVal(list())  # bw/hic/bed/paired-end rows in the currently-open add-sample modal
  modal_row_id    <- reactiveVal(0)
  # NULL = the modal is adding a brand-new sample; a string = the condition
  # name of the sample currently being edited, so wiz_confirm_sample can
  # replace it in place instead of appending a duplicate. A name, not a row
  # number: a Remove handled while the modal opens would shift the rows.
  # Set by wiz_add_sample/wiz_edit_sample below, cleared on
  # confirm -- a re-click of either always resets it, so a stale value can't
  # survive a Cancel.
  modal_edit_name <- reactiveVal(NULL)
  # Row of the sample being edited in wiz_samples(), NULL when adding, NA if it was removed meanwhile.
  modal_edit_index <- function() {
    nm <- modal_edit_name()
    if (is.null(nm)) NULL else match(nm, vapply(wiz_samples(), function(s) s$condition, character(1)))
  }
  # Per-row/per-field error messages from the last failed wiz_confirm_sample
  # attempt, so a validation error can highlight the specific offending
  # row/file instead of only listing errors in one combined banner. Keyed by
  # row id (character) for modal_row_errors, by "quantdir" for
  # modal_field_errors. Both reset to list() on a genuinely fresh modal open.
  modal_row_errors   <- reactiveVal(list())
  modal_field_errors <- reactiveVal(list())
  # Bumped on every failed Save: re-showing the modal replaces
  # uiOutput("wiz_instance_rows"), and shiny.js repaints it from its cached
  # last render, so without a fresh render the rows snap back to stale values.
  modal_fail_count   <- reactiveVal(0)
  # NULL = landing screen (choose a path); "import" / "manual" = which entry
  # path was chosen; "edit" = re-entered from the running app to edit an
  # existing project, which skips the landing screen entirely --
  # there's no "start fresh vs. import" choice to make when a project already
  # exists. Every non-NULL mode renders the same wiz_shared_form() below an
  # optional import box -- see wiz_shared_form()'s comment.
  # A first-time wizard that finds the project_samples.txt an
  # earlier Save left beside project_settings.json (settings since reset or
  # deleted) opens straight on Import, prefilled, so the project rebuilds in
  # one click. Its own name, so it never overwrites a user's samples.txt.
  mirror_txt <- file.path(dirname(SETTINGS_JSON), "project_samples.txt")
  has_mirror <- !is_edit && file.exists(mirror_txt)
  wiz_mode <- reactiveVal(if (is_edit) "edit" else if (has_mirror) "import" else NULL)
  # Display settings imported from a samples.txt (the wizard has no inputs for them):
  # Save writes cfg's values with these on top.
  import_visual <- reactiveVal(list())
  imported      <- reactiveVal(NULL)  # the lines of the samples.txt last imported in this wizard session

  observeEvent(input$wiz_mode_import, wiz_mode("import"))
  observeEvent(input$wiz_mode_manual, wiz_mode("manual"))
  observeEvent(input$wiz_mode_back, { wiz_mode(NULL); save_failed(FALSE); write_error(NULL) }, ignoreInit = TRUE)

  # Snapshot of every field this form can change, taken once at edit-mode
  # entry -- the baseline is_dirty() compares live values against, for the
  # Close-without-saving confirm guard. Not used
  # (and not needed) in fresh-setup mode: there's no prior configuration to
  # accidentally discard.
  original_snapshot <- if (is_edit) {
    list(
      data_root = cfg$data_dir,
      samples   = cfg$samples,
      flat      = cfg[c("font", "fontsize", "page_w", "chr_fmt", "show_track_axes", "hic_norm", "hic_style",
                         "zrange_cap", "atac_height", "hic_height", "hic_res_mode",
                         "hic_max_bins_per_view", "hic_resolution", "diff_colour_lo", "diff_colour_hi",
                         "junction_opacity_cap", "min_junction_support", "reference_colour", "comparison_colour",
                         "gene_flank_pct", "min_display_tpm", "uniprot_reviewed_only", "arc_clip_partial",
                         "show_title", "show_scale_bar", "junction_height", "diff_shared_alpha",
                         "gtf_path", "orgdb", "genome_label")]
    )
  } else NULL

  output$wiz_body <- renderUI({
    mode <- wiz_mode()
    if (is.null(mode)) {
      need <- function(title, ...) div(class = "wiz-need", tags$b(title), ...)
      return(tagList(
        h3("What a project needs", style = "margin-bottom:12px;"),
        div(class = "wiz-needs",
          need("Your data files",
               "Coverage (bigWig), Hi-C maps, BED intervals, splice-junction or loop files, salmon quant folders. ",
               "In one folder on the computer running R."),
          need("A GTF",
               "The gene annotation from GENCODE or Ensembl, for the same genome build and release your data was quantified with. ",
               "Not RefSeq. The next screen says where to get it."),
          need("The species' annotation package",
               "org.Mm.eg.db (mouse) or org.Hs.eg.db (human), installed in R. Gene search uses it to find genes by name.")
        ),
        h3("How do you want to start?", style = "margin-bottom:12px;"),
        div(class = "wiz-choices",
          div(class = "wiz-choice",
              h3("Fill in a form"),
              p("Type the genome and data folder, then add each sample and its files from a list. Best for a handful of samples."),
              actionButton("wiz_mode_manual", "Set up by hand", class = "btn-primary btn-lg")),
          div(class = "wiz-choice",
              h3("Import a samples.txt"),
              p("A tab-separated table listing your samples, written by you, a colleague or a chat AI (the vignette shows how). It prefills the form."),
              actionButton("wiz_mode_import", "Import a samples.txt", class = "btn-default btn-lg"))
        )
      ))
    }
    tagList(
      if (!identical(mode, "edit"))
        actionLink("wiz_mode_back", "\u2190 Back to the start", class = "wiz-back"),
      if (identical(mode, "import")) div(class = "wiz-section",
        div(class = "wiz-section-head", tags$h3("Import a samples.txt")),
        if (has_mirror)
          p(class = "wiz-section-intro", "A project_samples.txt saved with your last project was found. Click Import to restore it.")
        else
          p(class = "wiz-section-intro", "The file fills in the three sections below; check them and add anything missing before saving."),
        fluidRow(
          column(9, textInput("wiz_import_path", "Path to the samples.txt (on the computer running R)",
                              value = if (has_mirror) mirror_txt else "", width = "100%",
                              placeholder = "e.g. /data/smithlab/tcell_2026/samples.txt")),
          column(3, style = "margin-top: 25px;", actionButton("wiz_import_load", "Import", class = "pt-once btn-primary"))
        ),
        uiOutput("wiz_import_status")
      ),
      wiz_shared_form(cfg, is_edit)
    )
  })

  # Quick-fill the genome identity fields from a bundled GENOME_PRESETS
  # entry -- pure convenience, nothing is stored under the preset's id; Save
  # always writes the plain field values. GENOME_PRESETS carries no
  # gtf_path, so the GTF path field is left alone: always the user's own
  # file, and a path already entered survives a preset pick.
  observeEvent(input$wiz_genome_preset, {
    preset <- GENOME_PRESETS[[input$wiz_genome_preset]]
    if (is.null(preset)) return()
    updateTextInput(session, "wiz_orgdb", value = preset$orgdb)
    updateTextInput(session, "wiz_genome_label", value = preset$genome_label)
  }, ignoreInit = TRUE)

  data_root <- reactive({
    dr <- trimws(input$wiz_data_root %||% "")
    if (nchar(dr) == 0) NULL else dr
  })

  output$wiz_data_root_status <- renderUI({
    dr <- data_root()
    if (is.null(dr)) return(NULL)
    if (dir.exists(dr)) div(class = "wiz-status wiz-ok", "\u2713 Folder found.")
    else div(class = "wiz-status wiz-error", "\u2717 No folder at this path (it must be a path on the computer running R).")
  })

  # Live checks under the GTF and OrgDb fields (wiz_gtf_check(), wiz_orgdb_check()); advice,
  # never a block: Save's own checks are wiz_save_problems().
  # The typed path is debounced, so the file is read once per pause, not per keystroke.
  gtf_path_settled <- debounce(reactive(input$wiz_gtf_path), 400)
  gtf_check <- reactive(wiz_gtf_check(gtf_path_settled()))
  status_line <- function(chk) {
    if (identical(chk$state, "blank")) return(NULL)
    mark <- switch(chk$state, ok = "\u2713 ", warn = "\u26a0 ", "\u2717 ")
    div(class = paste("wiz-status", switch(chk$state, ok = "wiz-ok", warn = "wiz-warn", "wiz-error")),
        paste0(mark, chk$text))
  }
  output$wiz_gtf_status   <- renderUI(status_line(gtf_check()))
  output$wiz_orgdb_status <- renderUI(status_line(wiz_orgdb_check(input$wiz_orgdb, gtf_check()$species)))
  output$wiz_genome_label_status <- renderUI(status_line(wiz_genome_label_check(input$wiz_genome_label, gtf_check()$build, gtf_check()$species)))

  # ── Add-sample modal ──────────────────────────────────────────────────────

  # `error` NULL means a fresh modal open (fields blank/default); non-NULL means a re-display after a failed wiz_confirm_sample validation (fields seeded from current inputs instead, since Shiny never clears an unmounted input, so `error`'s presence is the only way to tell the two apart).
  add_sample_modal <- function(error = NULL) {
    dr    <- data_root()
    fresh <- is.null(error)
    seed  <- function(input_id, default) if (fresh) default else (isolate(input[[input_id]]) %||% default)
    field_errors <- if (fresh) list() else modal_field_errors()
    # Editing an existing sample (modal_edit_index() set) prefills name/
    # colour/quantdir from that sample. Adding a brand-new sample prefills
    # only the quantdir NAME (not path): scans every already-committed sample for a quantdir
    # entry, most recent match wins, rather than only the very last sample --
    # so a naming convention set early on isn't lost just because the most
    # recent sample happens to have no quantdir at all. Nothing
    # sample-specific (condition name, colour, actual paths) is prefilled
    # into a new sample.
    edit_idx       <- modal_edit_index()
    editing_sample <- if (!is.null(edit_idx) && !is.na(edit_idx)) wiz_samples()[[edit_idx]] else NULL
    carry_quantdir <- if (is.null(editing_sample)) {
      quantdirs <- Filter(Negate(is.null), lapply(wiz_samples(), function(s) Find(function(d) identical(d$type, "quantdir"), s$data)))
      if (length(quantdirs) > 0) quantdirs[[length(quantdirs)]] else NULL
    } else NULL
    quantdir_inst  <- if (!is.null(editing_sample)) Find(function(d) identical(d$type, "quantdir"), editing_sample$data) else carry_quantdir
    default_name    <- editing_sample$condition %||% ""
    default_colour  <- editing_sample$colour %||% WIZARD_DEFAULT_COLOUR
    default_qname   <- quantdir_inst$name %||% "Quant"
    default_qpath   <- if (!is.null(editing_sample)) (quantdir_inst$path %||% "") else ""
    # Explicit "" choice: with only real file paths in `choices`,
    # selected="" matches no <option>, and
    # shiny's underlying native <select> then falls back to auto-selecting
    # the first real option before selectize.js upgrades it -- so a field the
    # user never touched silently pointed at whatever file happened to be
    # first. A real blank option gives selected="" something to bind to.
    quantdir_choices <- c("", wizard_scoped_choices(dr, "quantdir"))
    modalDialog(
      title = if (is.null(edit_idx)) "Add a sample" else "Edit sample", size = "l", easyClose = FALSE,
      if (!is.null(error)) wiz_error_list("Please fix before adding:", error),
      fluidRow(
        column(6, textInput("wiz_sample_name",
                            wiz_help("Condition", "The name this sample's tracks carry in the figure, e.g. Control, KO, Basal, LP. One per sample; it must be unique."),
                            value = seed("wiz_sample_name", default_name), placeholder = "e.g. Control")),
        column(6, tags$label(wiz_help("Colour", "The colour of this sample's coverage and junction tracks. Pick one per condition that stays readable next to the others.")), br(),
               # shiny.js has no input binding for type="color" (.pt_colour_input()
               # works around the same gap the same way) -- without the oninput
               # handler below, input$wiz_sample_colour never updates from its initial value.
               tags$input(type = "color", id = "wiz_sample_colour", value = seed("wiz_sample_colour", default_colour),
                          style = "width:60px; height:34px; padding:2px; border:1px solid #ccc; border-radius:4px; cursor:pointer; display:block;",
                          oninput = "Shiny.setInputValue('wiz_sample_colour', this.value, {priority: 'event'});"),
               # Report the shown colour as soon as the modal renders:
               # oninput alone leaves input$wiz_sample_colour at whatever the
               # previous modal set, so an untouched swatch saved that colour.
               tags$script(sprintf("Shiny.setInputValue('wiz_sample_colour', %s);",
                                   jsonlite::toJSON(seed("wiz_sample_colour", default_colour), auto_unbox = TRUE))))
      ),
      hr(),
      tags$label(wiz_help("Files for this sample",
                          "Add one row per file. Give each a track name (e.g. ATAC, H3K27ac, HiC) and use the same track name for the same kind of data in every sample: that is how the app lines them up across conditions. A row left without a path is ignored.")),
      p(class = "text-muted", style = "margin-top:-2px;",
        "Pick the kind of file, then its path from the list (every matching file under the data folder is listed; typing a path works too)."),
      div(class = "wiz-file-kinds",
        actionButton("wiz_add_bw",         paste0("+ ", WIZ_KIND_LABELS[["bw"]]),         class = "pt-once btn-sm btn-default"),
        actionButton("wiz_add_hic",        paste0("+ ", WIZ_KIND_LABELS[["hic"]]),        class = "pt-once btn-sm btn-default"),
        actionButton("wiz_add_bed",        paste0("+ ", WIZ_KIND_LABELS[["bed"]]),        class = "pt-once btn-sm btn-default"),
        actionButton("wiz_add_paired_end", paste0("+ ", WIZ_KIND_LABELS[["paired-end"]]), class = "pt-once btn-sm btn-default")
      ),
      uiOutput("wiz_instance_rows"),
      hr(),
      tags$label(wiz_help("Transcript quantification (optional)",
                          "A salmon output folder for this sample (the one holding quant.sf). It gives the transcript expression and transcript count tracks. One per sample.")),
      fluidRow(class = if (is.null(field_errors$quantdir)) "wiz-row" else "wiz-row wiz-row-error",
        column(2, div(class = "wiz-kind", "Salmon folder")),
        column(4, textInput("wiz_quantdir_name", NULL, value = seed("wiz_quantdir_name", default_qname), placeholder = "track name")),
        column(6, selectizeInput("wiz_quantdir_path", NULL,
                                  choices = unique(c(quantdir_choices, seed("wiz_quantdir_path", default_qpath))),
                                  selected = seed("wiz_quantdir_path", default_qpath),
                                  options = list(create = TRUE, placeholder = "folder with quant.sf, relative to the data folder"))),
        if (!is.null(field_errors$quantdir)) column(12, tags$div(class = "wiz-error", style = "font-size:13px; margin-top:4px;", field_errors$quantdir))
      ),
      footer = tagList(modalButton("Cancel"), actionButton("wiz_confirm_sample", if (is.null(edit_idx)) "Add sample" else "Save changes", class = "pt-once btn-primary"))
    )
  }

  # `samples` is a list to prefill rows from: keep_path=TRUE for one sample being edited, FALSE for the deduped union (first-seen order) of every committed sample's type:name pairs, so a pair missing from the latest sample isn't lost; quantdir isn't included here, see add_sample_modal().
  wiz_prefill_rows <- function(samples, keep_path) {
    rows <- list()
    seen <- character(0)
    for (s in samples) {
      for (d in s$data) {
        if (identical(d$type, "quantdir")) next
        key <- paste0(d$type, ":", d$name)  # ":" isn't valid in a WIZARD_DATA_TYPES key, so this can't collide across types
        if (key %in% seen) next
        seen <- c(seen, key)
        rid <- modal_row_id() + 1
        modal_row_id(rid)
        rows[[length(rows) + 1]] <- list(id = rid, type = d$type, name = d$name,
                                          path = if (keep_path) d$path else "")
      }
    }
    rows
  }

  observeEvent(input$wiz_add_sample, {
    modal_edit_name(NULL)
    modal_instances(wiz_prefill_rows(wiz_samples(), keep_path = FALSE))
    modal_row_errors(list())
    modal_field_errors(list())
    # modal_row_id is NOT reset here: its ids are
    # input-element ids (wiz_inst_name_<id>/wiz_inst_path_<id>) that live
    # for the whole Shiny session even after their row is removed from the
    # DOM -- Shiny never clears input$wiz_inst_name_1 just because the
    # element disappeared. Resetting to 0 on every "+ Add sample" click meant
    # a fresh row on sample 2 reused id 1 and silently inherited sample 1's
    # already-committed values. Keeping the counter monotonic for the whole
    # session means an id is retired for good once used, so this can't recur.
    showModal(add_sample_modal())
  })

  observeEvent(input$wiz_edit_sample, {
    # By condition name, like Remove below: a row number can be stale after a quick Remove.
    idx <- match(input$wiz_edit_sample, vapply(wiz_samples(), function(s) s$condition, character(1)))
    if (is.na(idx)) return()
    modal_edit_name(input$wiz_edit_sample)
    modal_instances(wiz_prefill_rows(list(wiz_samples()[[idx]]), keep_path = TRUE))
    modal_row_errors(list())
    modal_field_errors(list())
    showModal(add_sample_modal())
  }, ignoreInit = TRUE)

  add_modal_row <- function(type) {
    rid <- modal_row_id() + 1
    modal_row_id(rid)
    modal_instances(c(modal_instances(), list(list(id = rid, type = type))))
  }
  observeEvent(input$wiz_add_bw,         add_modal_row("bw"))
  observeEvent(input$wiz_add_hic,        add_modal_row("hic"))
  observeEvent(input$wiz_add_bed,        add_modal_row("bed"))
  observeEvent(input$wiz_add_paired_end, add_modal_row("paired-end"))

  observeEvent(input$wiz_remove_row, {
    modal_instances(Filter(function(r) r$id != input$wiz_remove_row, modal_instances()))
  }, ignoreInit = TRUE)

  output$wiz_instance_rows <- renderUI({
    rows <- modal_instances()
    dr <- data_root()
    errs <- modal_row_errors()
    modal_fail_count()  # dependency only: re-render from current inputs after every failed Save
    if (length(rows) == 0) return(p(class = "wiz-empty", "No files yet. Use the buttons above to add one."))
    tagList(lapply(rows, function(r) {
      name_id <- paste0("wiz_inst_name_", r$id)
      path_id <- paste0("wiz_inst_path_", r$id)
      row_err <- errs[[as.character(r$id)]]
      fluidRow(class = if (is.null(row_err)) "wiz-row" else "wiz-row wiz-row-error",
        column(2, div(class = "wiz-kind", WIZ_KIND_LABELS[[r$type]])),
        # r$name/r$path: a row added via wiz_prefill_rows()
        # (edit-sample prefill, or a new sample carrying over the previous
        # one's instance layout) carries its prefill value here, since its
        # input id is freshly minted and isolate(input[[...]]) is NULL until
        # the user actually touches it. A row added via "+ <type> instance"
        # has no r$name/r$path (NULL), so this falls through to "" exactly
        # as before.
        column(4, textInput(name_id, NULL, value = isolate(input[[name_id]]) %||% r$name %||% "", placeholder = "track name, e.g. ATAC")),
        # "" prepended (see add_sample_modal()): otherwise
        # selected="" matches no real <option> and the browser silently
        # falls back to the first matching file.
        # The current value is also added to choices: a path typed via
        # create = TRUE that isn't in the scoped file list would otherwise be
        # dropped on re-render, since selected must be one of the choices.
        column(5, selectizeInput(path_id, NULL,
                                  choices = unique(c("", isolate(input[[path_id]]) %||% r$path %||% "",
                                                     wizard_scoped_choices(dr, r$type))),
                                  selected = isolate(input[[path_id]]) %||% r$path %||% "",
                                  options = list(create = TRUE, placeholder = "file, relative to the data folder"))),
        column(1, actionButton(paste0("wiz_rm_btn_", r$id), "\u2715", class = "btn-xs", title = "Remove this row",
                                onclick = sprintf("Shiny.setInputValue('wiz_remove_row', %d, {priority:'event'})", r$id))),
        if (!is.null(row_err)) column(12, tags$div(class = "wiz-error", style = "font-size:13px; margin-top:4px;", row_err))
      )
    }))
  })

  observeEvent(input$wiz_confirm_sample, {
    edit_idx <- modal_edit_index()
    if (identical(edit_idx, NA_integer_)) {  # the sample being edited was removed meanwhile
      showNotification(sprintf("Sample '%s' was removed, so these changes weren't saved.", modal_edit_name()), type = "warning")
      modal_edit_name(NULL)
      removeModal()
      return()
    }
    name <- trimws(input$wiz_sample_name %||% "")
    existing_names <- vapply(wiz_samples(), function(s) s$condition, character(1))
    # Editing a sample must not collide with its OWN pre-existing name --
    # excluding it here (rather than just allowing name == old name as a
    # special case) also lets the uniqueness check below stay a plain %in%.
    if (!is.null(edit_idx)) existing_names <- existing_names[-edit_idx]
    errs <- character(0)
    row_errs   <- list()  # row id (character) -> message, for wiz_instance_rows highlighting
    field_errs <- list()  # "quantdir" -> message, for add_sample_modal() highlighting

    if (nchar(name) == 0) errs <- c(errs, "Condition is required.")
    else if (name %in% existing_names) errs <- c(errs, sprintf("A condition named '%s' already exists.", name))

    data <- list()
    for (r in modal_instances()) {
      nm  <- trimws(input[[paste0("wiz_inst_name_", r$id)]] %||% "")
      pth <- trimws(input[[paste0("wiz_inst_path_", r$id)]] %||% "")
      if (nchar(pth) == 0) next  # blank instance = not declared, silently dropped
      if (nchar(nm) == 0) {
        msg <- sprintf("A %s instance has a path but no name.", r$type)
        errs <- c(errs, msg); row_errs[[as.character(r$id)]] <- msg
        next
      }
      if (wizard_name_has_bad_dot(nm)) {
        msg <- sprintf("Instance name '%s' contains a bare '\u00b7' (surround it with spaces, or remove it). This breaks track-label formatting.", nm)
        errs <- c(errs, msg); row_errs[[as.character(r$id)]] <- msg
        next
      }
      if (any(vapply(data, function(d) d$type == r$type && d$name == nm, logical(1)))) {
        msg <- sprintf("Duplicate %s instance name '%s'.", r$type, nm)
        errs <- c(errs, msg); row_errs[[as.character(r$id)]] <- msg
        next
      }
      # is_real_file(), not file.exists(): file.exists() is TRUE for
      # directories too, which would otherwise be accepted here only to fail
      # (or, for paired-end, crash) at first render.
      if (!is.null(data_root()) && !is_real_file(file.path(data_root(), pth))) {
        msg <- sprintf("Path for %s:%s does not resolve to a file under the data root: %s", r$type, nm, pth)
        errs <- c(errs, msg); row_errs[[as.character(r$id)]] <- msg
      } else if (identical(r$type, "paired-end") && !is.null(data_root()) &&
                 !is_paired_end_structural(file.path(data_root(), pth))) {
        # Cheap fail-fast check (is_paired_end_structural()): must be BED12
        # (12 cols), BED6 or BEDPE (6-10 cols) or a HiCCUPS loop file, the
        # formats read_junction_file() reads.
        msg <- sprintf("Paired-end file for %s is not BED12 (12 columns), BED6 or BEDPE (6-10 columns), or a HiCCUPS loop file: %s", nm, pth)
        errs <- c(errs, msg); row_errs[[as.character(r$id)]] <- msg
      }
      data[[length(data) + 1]] <- list(type = r$type, name = nm, path = pth)
    }

    quantdir_path <- trimws(input$wiz_quantdir_path %||% "")
    if (nchar(quantdir_path) > 0) {
      quantdir_name <- trimws(input$wiz_quantdir_name %||% "")
      if (nchar(quantdir_name) == 0) quantdir_name <- "Quant"
      if (!is.null(data_root()) && !is_real_file(file.path(data_root(), quantdir_path, "quant.sf"))) {
        msg <- sprintf("No quant.sf found under: %s", quantdir_path)
        errs <- c(errs, msg); field_errs$quantdir <- msg
      }
      data[[length(data) + 1]] <- list(type = "quantdir", name = quantdir_name, path = quantdir_path)
    }

    if (length(errs) > 0) {
      modal_row_errors(row_errs)
      modal_field_errors(field_errs)
      modal_fail_count(modal_fail_count() + 1)
      showModal(add_sample_modal(error = errs))
      return()
    }

    colour <- input$wiz_sample_colour %||% WIZARD_DEFAULT_COLOUR
    new_sample <- list(condition = name, colour = colour, data = data)
    if (is.null(edit_idx)) {
      wiz_samples(c(wiz_samples(), list(new_sample)))
    } else {
      samples <- wiz_samples()
      samples[[edit_idx]] <- new_sample
      wiz_samples(samples)
    }
    modal_edit_name(NULL)
    removeModal()
  })

  # ── Committed-sample list + edit/remove ───────────────────────────────────
  output$wiz_samples_ui <- renderUI({
    samples <- wiz_samples()
    if (length(samples) == 0) return(p(class = "wiz-empty", "No samples yet. Add one per condition; you need at least one to save."))
    div(class = "wiz-samples", lapply(seq_along(samples), function(i) {
      s <- samples[[i]]
      inst_labels <- if (length(s$data) > 0) vapply(s$data, function(d) paste0(d$type, ":", d$name), character(1)) else character(0)
      fluidRow(class = "wiz-row",
        column(1, tags$span(style = sprintf("display:inline-block;width:14px;height:14px;background:%s;border-radius:2px;", s$colour %||% "#ccc"))),
        column(3, tags$b(s$condition)),
        column(5, paste(inst_labels, collapse = ", ")),
        column(3,
               actionButton(paste0("wiz_edit_sample_btn_", i), "Edit", class = "pt-once btn-xs btn-default",
                             onclick = sprintf("Shiny.setInputValue('wiz_edit_sample', %s, {priority:'event'})",
                                               jsonlite::toJSON(s$condition, auto_unbox = TRUE))),
               actionButton(paste0("wiz_rm_sample_btn_", i), "Remove", class = "pt-once btn-xs btn-default",
                             onclick = sprintf("Shiny.setInputValue('wiz_remove_sample', %s, {priority:'event'})",
                                               jsonlite::toJSON(s$condition, auto_unbox = TRUE))))
      )
    }))
  })

  # By condition name (unique in the wizard), not row number: a second quick Remove sends a row
  # number from the table as it was before the first removal, which then names another sample.
  observeEvent(input$wiz_remove_sample, {
    wiz_samples(Filter(function(s) !identical(s$condition, input$wiz_remove_sample), wiz_samples()))
  }, ignoreInit = TRUE)

  # ── Import from samples.txt ───────────────────────────────────────────────
  output$wiz_import_status <- renderUI(NULL)

  observeEvent(input$wiz_import_load, {
    path <- trimws(input$wiz_import_path %||% "")
    if (nchar(path) == 0) {
      output$wiz_import_status <- renderUI(div(class = "wiz-error", "Enter a path to samples.txt first."))
      return()
    }
    if (!is_real_file(path)) {
      output$wiz_import_status <- renderUI(div(class = "wiz-error", sprintf("No such file: %s", path)))
      return()
    }
    res <- tryCatch(
      # defaults = DEFAULT_SETTINGS: lets any `# <key>: value` line in the
      # file be parsed back into `res$overrides`, applied below: the 3 genome
      # fields (gtf_path/orgdb/genome_label) fill their inputs, the display
      # settings are kept for Save. A line for a key that isn't a setting
      # (e.g. a cache path like `# txdb_path: ...`) is ignored.
      import_samples_txt(path, defaults = DEFAULT_SETTINGS),
      error = function(e) list(errors = conditionMessage(e))
    )
    if (length(res$errors) > 0) {
      output$wiz_import_status <- renderUI(
        wiz_error_list("Nothing imported. The file has these problems:", res$errors)
      )
      return()
    }
    if (!is.null(res$data_root)) updateTextInput(session, "wiz_data_root", value = res$data_root)
    wiz_samples(res$samples)
    imported(readLines(path, warn = FALSE))

    # Apply any `# <key>: value` overrides: the genome fields onto their
    # inputs, the display settings into import_visual() for Save. Numbers for the
    # NULL-default ("auto") keys are converted the way parse_samples_txt() does.
    coerced <- .pt_coerce_auto_overrides(res$overrides, DEFAULT_SETTINGS)
    ov <- coerced$overrides
    res$warnings <- c(res$layout_warning, res$warnings, coerced$ignored)
    visual_keys <- setdiff(names(DEFAULT_SETTINGS), .PT_CORE_SETTINGS)
    import_visual(ov[intersect(names(ov), visual_keys)])
    if (!is.null(ov$gtf_path))            updateTextInput(session, "wiz_gtf_path", value = ov$gtf_path)
    if (!is.null(ov$orgdb))               updateTextInput(session, "wiz_orgdb", value = ov$orgdb)
    if (!is.null(ov$genome_label))        updateTextInput(session, "wiz_genome_label", value = ov$genome_label)

    status_lines <- list(div(class = "wiz-ok", sprintf("Imported %d sample(s).", length(res$samples))))
    if (length(res$warnings) > 0)
      status_lines <- c(status_lines, list(wiz_error_list("Imported, with warnings:", res$warnings, kind = "warn")))
    output$wiz_import_status <- renderUI(tagList(status_lines))
  })

  # ── Export samples.txt ────────────────────────────────────────────────────

  # A separate action from Save, never triggered automatically.
  wiz_collect_flat_settings <- function() {
    c(list(
      gtf_path            = wiz_clean_optional_text(input$wiz_gtf_path),
      orgdb               = wiz_clean_optional_text(input$wiz_orgdb),
      genome_label        = wiz_clean_optional_text(input$wiz_genome_label)
    ), wiz_display_settings())
  }
  # The saved display settings with any imported ones on top (no wizard inputs for them).
  wiz_display_settings <- function() {
    x <- cfg
    iv <- import_visual()
    x[names(iv)] <- iv
    .pt_collect_display_settings(list(), "", x)
  }

  observeEvent(input$wiz_export_samples, {
    dr <- data_root()
    showModal(modalDialog(
      title = "Export settings to txt", easyClose = TRUE,
      textInput("wiz_export_path", "Write to (full path):",
                value = if (is.null(dr)) "" else file.path(dr, "samples.txt"), width = "100%"),
      uiOutput("wiz_export_status"),
      footer = tagList(modalButton("Cancel"), actionButton("wiz_confirm_export", "Export", class = "pt-once btn-primary"))
    ))
  })

  output$wiz_export_status <- renderUI(NULL)

  observeEvent(input$wiz_confirm_export, {
    path <- trimws(input$wiz_export_path %||% "")
    if (nchar(path) == 0) {
      output$wiz_export_status <- renderUI(div(class = "wiz-error", "Enter a path first."))
      return()
    }
    settings <- c(
      list(data_root = data_root(), samples = wiz_samples()),
      wiz_collect_flat_settings()
    )
    ok <- tryCatch({
      writeLines(export_samples_txt(settings, defaults = DEFAULT_SETTINGS), path)
      TRUE
    }, error = function(e) {
      output$wiz_export_status <- renderUI(div(class = "wiz-error", sprintf("Failed to write: %s", conditionMessage(e))))
      FALSE
    })
    if (ok) output$wiz_export_status <- renderUI(div(class = "wiz-ok", sprintf("Wrote %s.", path)))
  })

  # ── Dirty guard (edit mode only) ──────────────────────────────────────────
  is_dirty <- reactive({
    if (!is_edit) return(FALSE)
    if (!identical(data_root(), original_snapshot$data_root)) return(TRUE)
    if (!identical(wiz_samples(), original_snapshot$samples)) return(TRUE)
    flat <- wiz_collect_flat_settings()
    any(vapply(names(flat), function(k) settings_value_differs(flat[[k]], original_snapshot$flat[[k]]), logical(1)))
  })

  wiz_do_close <- function() {
    showModal(modalDialog(
      title = "Closing", easyClose = FALSE, footer = NULL,
      p("Returning to the app without saving."),
      p(strong("You can close this tab now."), " The app reopens in a new tab from RStudio in a few seconds. ",
        "If nothing opens, go back to RStudio: the app's address is in the Console.")
    ))
    message("PlotTrellis: settings unchanged, returning to the app.")
    later::later(function() stopApp(returnValue = "restart"), 0.2)
  }

  observeEvent(input$wiz_close, {
    if (is_dirty()) {
      showModal(modalDialog(
        title = "Discard changes?", easyClose = FALSE,
        p("You have unsaved changes. Closing now will discard them and return to the app with its previous settings."),
        footer = tagList(modalButton("Keep editing"), actionButton("wiz_confirm_close", "Discard and close", class = "pt-once btn-danger"))
      ))
    } else {
      wiz_do_close()
    }
  })

  observeEvent(input$wiz_confirm_close, {
    removeModal()
    wiz_do_close()
  }, ignoreInit = TRUE)

  # ── Save ───────────────────────────────────────────────────────────────────
  # After a Save that failed, the "Not saved yet" list re-checks as the fields change, so it
  # never shows a problem already fixed (e.g. by an Import), and goes once all pass.
  save_failed <- reactiveVal(FALSE)
  # The last Save's write failure (a folder that can't be written), cleared by the next Save.
  write_error <- reactiveVal(NULL)
  save_problems <- reactive({
    input$wiz_save  # each Save click re-checks too (see the Save observer)
    wiz_save_problems(data_root(), wiz_samples(), input$wiz_gtf_path, input$wiz_orgdb, input$wiz_genome_label)
  })
  output$wiz_completion_error <- renderUI({
    if (!is.null(write_error())) return(wiz_error_list("Not saved:", list(write_error())))
    if (!save_failed()) return(NULL)
    errs <- save_problems()
    if (length(errs) > 0) wiz_error_list("Not saved yet. Please fix:", errs)
  })

  observeEvent(input$wiz_save, {
    write_error(NULL)
    dr <- data_root()
    # Checked afresh, not from save_problems(): a folder created or an OrgDb installed since
    # changes the answer without changing any field.
    errs <- wiz_save_problems(dr, wiz_samples(), input$wiz_gtf_path, input$wiz_orgdb, input$wiz_genome_label)
    if (length(errs) > 0) {
      save_failed(TRUE)
      return()
    }

    settings <- c(
      list(
        data_dir        = dr,
        samples         = wiz_samples()
      ),
      wiz_collect_flat_settings()
    )
    # Kept across edits: pt_project() uses it to tell this project's samples.txt
    # from a different one pointed at the same folder; and the start view, while
    # what it refers to is unchanged (.pt_kept_start_view()). After an Import the imported
    # file restates the project: the start view is the one its start lines give, as in
    # pt_project(), or none. Lines that don't fit the samples don't block the Save.
    if (file.exists(SETTINGS_JSON)) {
      old <- tryCatch(jsonlite::read_json(SETTINGS_JSON), error = function(e) list())
      settings$samples_txt_source <- old$samples_txt_source
      if (is.null(imported())) settings$start_view <- .pt_kept_start_view(old, settings)
    }
    start_problem <- NULL
    if (!is.null(imported()))
      settings$start_view <- tryCatch(withCallingHandlers(.pt_txt_start_view(imported(), settings),
        warning = function(w) { start_problem <<- c(start_problem, conditionMessage(w)); invokeRestart("muffleWarning") }),
        error = function(e) { start_problem <<- c(start_problem, conditionMessage(e)); NULL })
    if (!is.null(start_problem))
      message("PlotTrellis: the imported samples.txt's start lines: ", paste(start_problem, collapse = " "))
    # Always keep a plain-text project_samples.txt beside it, the
    # human-readable record the project can be rebuilt from (mirror_txt).
    # A folder that can't be written (permissions, a full disk) says so here instead of
    # ending the session.
    txt_ok <- tryCatch(.pt_write_project_files(settings, SETTINGS_JSON), error = function(e) {
      write_error(paste0("Couldn't write ", SETTINGS_JSON, " (", conditionMessage(e), "). ",
                         "Check that you can write to that folder, then Save again."))
      NULL
    })
    if (is.null(txt_ok)) return()

    save_failed(FALSE)
    showModal(modalDialog(
      title = "Settings saved", easyClose = FALSE, footer = NULL,
      p("project_settings.json has been written."),
      p(if (txt_ok) paste0("A copy of your samples is in ", mirror_txt, "; it can rebuild the project if the settings ever break.")
        else "Couldn't write project_samples.txt beside it (use Export samples.txt to save a copy yourself)."),
      if (!is.null(start_problem))
        p(class = "wiz-warn", if (is.null(settings$start_view))
            "The imported samples.txt's # start_gene: / # start_tracks: lines don't fit its samples, so the app opens with no gene or tracks picked. "
          else "The imported samples.txt's # start_ lines gave this note. ",
          paste(start_problem, collapse = " ")),
      p(strong("You can close this tab now."), " The app reopens with your project in a new tab from RStudio in a few seconds. ",
        "If nothing opens, go back to RStudio: the app's address is in the Console.")
    ))
    message("PlotTrellis: settings saved, launching the app.")
    later::later(function() stopApp(returnValue = "restart"), 0.2)
  })
}
