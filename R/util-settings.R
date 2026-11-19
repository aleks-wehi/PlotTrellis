# Project settings shared by the setup wizard and the composer: which settings
# are core (set only in the wizard) and which are display settings, the display
# settings' inputs (the composer's Figure settings panel) and reading them back,
# and writing a project's project_settings.json with its project_samples.txt
# copy.

# A colour swatch with no shiny.js binding (type="color"): oninput reports
# changes, and the inline script reports the shown value on render, so a
# re-rendered swatch (the composer's Revert) never leaves input[[id]] at a stale
# colour.
.pt_colour_input <- function(id, label, title, value) {
  tagList(tags$label(label, title = title), br(),
          tags$input(type = "color", id = id, value = value,
                     style = "width:60px; height:34px; padding:2px; border:1px solid #ccc; border-radius:4px; cursor:pointer; display:block;",
                     oninput = sprintf("Shiny.setInputValue('%s', this.value, {priority: 'event'});",
                                       id)),
          tags$script(sprintf("Shiny.setInputValue('%s', %s);", id,
                              jsonlite::toJSON(value, auto_unbox = TRUE))))
}

# Core settings: what identifies and locates the project's data, set only in
# the wizard, and the genes the experienced screen's gene box suggests, set only
# by a samples.txt (`# suggested_genes:`). Every other .pt_default_settings()
# key is a display setting, adjustable live in the app.
.PT_CORE_SETTINGS <- c("gtf_path", "orgdb", "genome_label", "data_dir",
                       "samples", "suggested_genes")

# Inputs for every display setting, ids paste0(prefix, key), initial values from
# cfg: the composer's Figure settings panel. Returns list(display =, diff =):
# the general settings and the 2-sample comparison ones.
.pt_display_settings_ui <- function(cfg, prefix) {
  list(
    display = tagList(
      .pt_settings_group("Layout & display",
      fluidRow(
        column(3, selectizeInput(paste0(prefix, "font"), "Font:",
                                 choices = c("DejaVu Sans", "Helvetica",
                                             "Arial"),
                                 selected = cfg$font %||% "DejaVu Sans")),
        column(3, numericInput(paste0(prefix, "fontsize"), "Font size:",
                               value = cfg$fontsize %||% 8, min = 4, max = 24)),
        column(3, numericInput(paste0(prefix, "page_w"),
                                wiz_help("Page width (in):",
                                         "Width of the exported figure, in inches."),
                                value = cfg$page_w %||% 9, min = 3, max = 30)),
        column(3, selectizeInput(paste0(prefix, "chr_fmt"),
                               wiz_help("Chromosome format:",
                                        "How chromosome names are written on the axis."),
                               choices = c("Dot (chr.19)" = "dot",
                                           "Space (chr 19)" = "space",
                                           "None (chr19)" = "plain"),
                               selected = cfg$chr_fmt %||% "dot")),
        column(3, checkboxInput(paste0(prefix, "show_track_axes"),
                                 wiz_help("Show track value axes",
                                          "Draw a value axis on coverage and junction tracks."),
                                 value = cfg$show_track_axes %||% TRUE))
      ),
      fluidRow(
        column(3, checkboxInput(paste0(prefix, "show_title"),
                                 wiz_help("Show title",
                                          "Draw the title over the figure."),
                                 value = cfg$show_title %||% TRUE)),
        column(3, checkboxInput(paste0(prefix, "show_scale_bar"),
                                 wiz_help("Show scale bar",
                                          "Draw a genomic distance scale bar."),
                                 value = cfg$show_scale_bar %||% TRUE))
      )),

      .pt_settings_group(
        "Coverage tracks (BigWig: ATAC / ChIP-seq / CUT&RUN / RNA-seq)",
      fluidRow(
        column(3, numericInput(paste0(prefix, "atac_height"),
                                wiz_help("Track height (in):",
                                         "Height in inches of each coverage track."),
                                value = cfg$atac_height %||% 0.85, min = 0.3,
                                max = 8, step = 0.1))
      )),

      .pt_settings_group("Hi-C",
      fluidRow(
        column(3, selectizeInput(paste0(prefix, "hic_style"),
                               wiz_help("Display style:",
                                        "How the Hi-C contact matrix is drawn: a full rectangle, or a triangle (upper-left half)."),
                               choices = c("Rectangle" = "rectangle",
                                           "Triangle" = "triangle"),
                               selected = cfg$hic_style %||% "rectangle")),
        # The usual names (HIC_NORMS); another can be typed in. Which ones a
        # file has is only known per file, so the list is not limited to them.
        column(3, selectizeInput(paste0(prefix, "hic_norm"),
                             wiz_help("Normalisation:",
                                      "How contact counts are balanced: a method your files have (KR, SCALE, VC, VC_SQRT; BALANCE in a cooler-balanced .cool or .mcool), another typed in, or NONE for raw counts, which every file has."),
                             choices = unique(c(HIC_NORMS,
                                                cfg$hic_norm %||% "KR")),
                             selected = cfg$hic_norm %||% "KR",
                             options = list(create = TRUE))),
        column(3, selectizeInput(paste0(prefix, "hic_colour_scale"),
                               wiz_help("Colour scale:",
                                        "Linear (as Juicebox) shows the strongest contacts; Log (as HiGlass) also shows the faint ones."),
                               choices = c("Linear" = "linear", "Log" = "log"),
                               selected = cfg$hic_colour_scale %||% "linear")),
        # One setting: a choice, or a number typed in (as Normalisation takes
        # a typed name). A saved number is listed so it shows as selected.
        column(3, selectizeInput(paste0(prefix, "zrange_cap"),
                               wiz_help("Colour max:",
                                        "The value that gets the darkest colour, from the contacts drawn: Auto is the 95th percentile on Linear and Max on Log, or type a number."),
                               choices = .pt_zrange_cap_choices(cfg$zrange_cap),
                               selected = .pt_zrange_cap_value(cfg$zrange_cap),
                               options = list(
                                 create = TRUE,
                                 createFilter = "^[0-9]*[.]?[0-9]+$"
                               )))
      ),
      fluidRow(
        column(3, selectizeInput(paste0(prefix, "hic_res_mode"),
                               wiz_help("Resolution:",
                                        "The bin size: Auto fits it to the view, Finest uses your file's smallest, Fixed the one you type; a view draws at most 1,000 bins across, so a wider one uses a coarser bin."),
                               choices = c("Auto" = "auto", "Finest" = "finest",
                                           "Fixed" = "fixed"),
                               selected = cfg$hic_res_mode %||% "auto")),
        # One column, whichever of the two the mode uses. The Auto setting's key
        # is kept from when it was named "Max bins per view" (saved settings and
        # samples.txt files use it).
        column(3, conditionalPanel(paste0("input.", prefix,
                                          "hic_res_mode === 'auto'"),
                 numericInput(paste0(prefix, "hic_max_bins_per_view"),
                                wiz_help("Min bins across view:",
                                         "Auto uses the coarsest bin size that still gives at least this many bins across the view."),
                                value = cfg$hic_max_bins_per_view %||% 500,
                                min = 10, max = HIC_MAX_BINS, step = 10)),
               conditionalPanel(paste0("input.", prefix,
                                       "hic_res_mode === 'fixed'"),
                 numericInput(paste0(prefix, "hic_resolution"),
                                wiz_help("Bin size (bp):",
                                         "The bin size for every Hi-C track (e.g. 10000); a file without it shows an error naming the sizes it has."),
                                value = clean_optional_numeric(
                                  cfg$hic_resolution
                                ) %||% NA,
                                min = 0, step = 1000))),
        column(3, numericInput(paste0(prefix, "hic_height"),
                                wiz_help("Track height (in):",
                                         "Height in inches of each Hi-C track. Leave blank to size automatically based on the display style."),
                                value = clean_optional_numeric(
                                  cfg$hic_height
                                ) %||% NA,
                                min = 0.5, max = 15, step = 0.5))
      )),

      .pt_settings_group("Junction / arc tracks",
      fluidRow(
        column(4, checkboxInput(paste0(prefix, "arc_clip_partial"),
                                 wiz_help("Clip partial arcs",
                                          "Leave out arcs that run past the edge of the region, instead of cutting them off."),
                                 value = cfg$arc_clip_partial %||% FALSE)),
        column(3, numericInput(paste0(prefix, "min_junction_support"),
                                wiz_help("Min arc support:",
                                         "Hide arcs whose read count (or loop score) is below this; 0 shows all."),
                                value = cfg$min_junction_support %||% 0,
                                min = 0, step = 1)),
        column(3, numericInput(paste0(prefix, "junction_height"),
                                wiz_help("Junction track height (in):",
                                         "Height in inches of a junction track's arcs (each half, when mirrored)."),
                                value = cfg$junction_height %||% 0.85,
                                min = 0.3, max = 8, step = 0.1))
      )),

      .pt_settings_group("Gene view",
      fluidRow(
        column(3, numericInput(paste0(prefix, "gene_flank_pct"),
                                wiz_help("Flanking (%):",
                                         "Extra region shown each side of a gene, as % of its length."),
                                value = cfg$gene_flank_pct %||% 15, min = 0,
                                max = 200, step = 5)),
        column(3, numericInput(paste0(prefix, "min_display_tpm"),
                                wiz_help("Min display TPM:",
                                         "Transcripts below this mean TPM in every condition are folded into 'Other'."),
                                value = cfg$min_display_tpm %||% 1, min = 0,
                                step = 0.5)),
        column(4, checkboxInput(paste0(prefix, "uniprot_reviewed_only"),
                                 wiz_help("UniProt: reviewed only",
                                          "Link only to manually reviewed (Swiss-Prot) UniProt entries."),
                                 value = cfg$uniprot_reviewed_only %||% FALSE))
      ))
    ),
    diff = tagList(
      fluidRow(
        column(3, .pt_colour_input(paste0(prefix, "diff_colour_lo"),
                                   "Loss colour (coverage diff):",
                                   "Colour of coverage lost in the comparison.",
                                   cfg$diff_colour_lo %||% "#0033CC")),
        column(3, .pt_colour_input(paste0(prefix, "diff_colour_hi"),
                                   "Gain colour (coverage diff):",
                                   "Colour of coverage gained in the comparison.",
                                   cfg$diff_colour_hi %||% "#CC0000")),
        column(3, numericInput(paste0(prefix, "junction_opacity_cap"),
                                wiz_help("Junctions opacity cap (log2FC):",
                                         "log2 fold change at which a comparison arc is fully opaque."),
                                value = cfg$junction_opacity_cap %||% 2,
                                min = 0.1, step = 0.5)),
        column(3, numericInput(paste0(prefix, "diff_shared_alpha"),
                                wiz_help("Shared coverage tint:",
                                         "Tint of the shared coverage, from 0 (white) to 1 (full reference colour)."),
                                value = cfg$diff_shared_alpha %||% 0.15,
                                min = 0, max = 1, step = 0.05))
      ),
      fluidRow(
        column(3, .pt_colour_input(paste0(prefix, "reference_colour"),
                                   "Reference colour:",
                                   "Colour of the reference condition in comparison tracks.",
                                   cfg$reference_colour %||% "#4D4D4D")),
        column(3, .pt_colour_input(paste0(prefix, "comparison_colour"),
                                   "Comparison colour:",
                                   "Colour of the comparison condition in comparison tracks.",
                                   cfg$comparison_colour %||% "#56B4E9"))
      )
    )
  )
}

# The Colour max list: Auto, the statistics, and a saved number.
.pt_zrange_cap_choices <- function(cap) {
  cap <- .pt_zrange_cap_value(cap)
  c(c("Auto" = "auto", "Max" = "max", "95th percentile" = "p95",
      "90th percentile" = "p90"),
    if (!cap %in% c("auto", "max", names(HIC_TOP_QUANTILES)))
      stats::setNames(cap, cap))
}

# A saved Colour max as the list's value: "auto" for NULL (or a value that is
# not a Colour max), else the statistic or the number as text.
.pt_zrange_cap_value <- function(cap) {
  cap <- .pt_hic_zrange_cap(cap, bad = NULL)
  if (is.null(cap)) "auto" else
    if (is.numeric(cap)) format(cap, scientific = FALSE, digits = 15) else cap
}

# Reads the display settings back from .pt_display_settings_ui()'s inputs. An
# input that isn't there (NULL, e.g. not rendered yet) or is blank falls back to
# cfg, then to the default, so it never silently resets a value. The optional
# numbers (blank = auto) only fall back when the input is missing. names() of
# the result are the display-setting keys.
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
    # A cleared Normalisation box sends "", which no file has: the saved one,
    # else KR, instead.
    hic_norm    = Find(nzchar, trimws(c(get("hic_norm"), cfg$hic_norm %||% "",
                                        "KR"))),
    hic_style   = get("hic_style"),
    hic_colour_scale = get("hic_colour_scale"),
    zrange_cap  = .pt_hic_zrange_cap(
      (function(v) if (is.null(v)) cfg$zrange_cap else v)(
        input[[paste0(prefix, "zrange_cap")]]
      ), bad = NULL
    ),
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

# Writes project_settings.json and, beside it, the project_samples.txt copy the
# project can be rebuilt from, with absolute paths since the importer resolves
# relative ones against the file's own folder. Shared by the wizard's Save and
# the composer's Save as project default. Returns whether the samples.txt copy
# was written; a failure there doesn't block.
.pt_write_project_files <- function(settings, settings_json) {
  # null = "null": jsonlite's toJSON() default serializes a NULL list element
  # (e.g. zrange_cap/hic_height left blank for "auto") as {} rather than the
  # JSON literal null. On reload, fromJSON() then reads that back as an empty
  # list, not NULL, and clean_optional_numeric()'s is.finite() crashes on it.
  # digits = NA: every writer of project_settings.json keeps full precision
  # (toJSON() rounds to 4 decimals). Through a temporary file renamed into
  # place, so a failed write (a full disk) never leaves a truncated settings
  # file; the error names the reason (e.g. "Permission denied"), which
  # writeLines() only gives as a warning.
  tmp <- tempfile("project_settings-", tmpdir = dirname(settings_json),
                  fileext = ".json")
  why <- character(0)
  ok <- tryCatch(withCallingHandlers({
    writeLines(jsonlite::toJSON(settings, auto_unbox = TRUE, pretty = TRUE,
                                null = "null", digits = NA), tmp)
    file.rename(tmp, settings_json)
  }, warning = function(w) {
    why <<- c(why, conditionMessage(w)); invokeRestart("muffleWarning")
  }),
  error = function(e) { why <<- c(why, conditionMessage(e)); FALSE })
  if (!isTRUE(ok)) {
    unlink(tmp)
    stop(paste(unique(why), collapse = "; "), call. = FALSE)
  }
  .pt_write_samples_copy(settings, settings_json)
}

# The project_samples.txt beside settings_json, from `settings`; whether it was
# written. A failure (error or warning) only gives FALSE.
.pt_write_samples_copy <- function(settings, settings_json) {
  tryCatch(withCallingHandlers({
    flat <- settings[setdiff(names(settings), c("data_dir", "samples"))]
    if (length(flat$gtf_path) == 1 && nzchar(flat$gtf_path))
      flat$gtf_path <- normalizePath(flat$gtf_path, mustWork = FALSE)
    txt <- export_samples_txt(
      c(list(data_root = normalizePath(settings$data_dir, mustWork = FALSE),
             samples = settings$samples), flat),
      defaults = .pt_default_settings()
    )
    # The start view as start lines after the other "# key:" lines, so the copy
    # restores it too.
    if (is.list(settings$start_view))
      txt <- append(txt, tryCatch(.pt_start_view_lines(settings$start_view,
                                                       settings),
                                  error = function(e) character(0)),
                    after = match(FALSE, startsWith(txt, "#"),
                                  nomatch = length(txt) + 1) - 1)
    writeLines(txt, file.path(dirname(settings_json), "project_samples.txt"))
    TRUE
  }, warning = function(w) stop(conditionMessage(w), call. = FALSE)),
  error = function(e) FALSE)
}

# The start view (pt_preload(), or a samples.txt's start_gene/start_tracks) of
# the saved settings `old`, if the wizard's new settings `new` still have
# everything it refers to: the same GTF and OrgDb (its gene is looked up in
# them), each of its samples, and each of its tracks (same key, i.e. same sample
# and data name; a changed file path or data folder keeps it). Otherwise NULL:
# the app then opens blank rather than on a view the user changed underneath it.
# Other samples can come and go.
.pt_kept_start_view <- function(old, new) {
  sv <- old$start_view
  if (!is.list(sv)) return(NULL)
  path <- function(p)
    if (is.character(p) && length(p) == 1 && nzchar(p))
      normalizePath(p, mustWork = FALSE) else ""
  if (!identical(path(old$gtf_path), path(new$gtf_path)) ||
      !identical(old$orgdb, new$orgdb))
    return(NULL)
  samples <- load_samples_from_settings(new$samples, new$data_dir)
  conds <- vapply(samples, function(s) s$display, character(1))
  keys <- names(c(build_track_items(samples),
                  build_comparison_track_items(samples)))
  if (!all(as.character(unlist(sv$conditions)) %in% conds) ||
      !all(.pt_renamed_junctions_key(as.character(unlist(sv$tracks))) %in%
             keys))
    return(NULL)
  sv
}
