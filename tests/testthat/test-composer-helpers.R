# Unit tests of the composer's small helpers (label formatting, chromosome sorting, number
# formatting and parsing, track pools, the Generate button, panels, the inlined SVG). They
# were written against the helpers inside .pt_composer_app(); `h` maps each old closure
# name to its top-level function (util-format.R, app-utils-ui.R, app-tracks.R,
# api-figure-helpers.R), with the app's startup context filled in where the old closure read it.
h <- function(name, cp) {
  ctx <- get("ctx", envir = environment(cp$app$serverFuncSource()))
  switch(name,
    start_condition     = function(i) .pt_start_condition(i, ctx$START_VIEW),
    track_labels        = function(keys) .pt_track_labels(keys, ctx$TRACK_ITEMS),
    sort_keys_by_label  = function(keys) .pt_sort_keys_by_label(keys, ctx$TRACK_ITEMS),
    default_figure_keys = function(conditions) .pt_default_figure_keys(conditions, ctx$START_VIEW),
    track_split         = function(conditions) .pt_track_split(conditions, ctx$TRACK_ITEMS, ctx$START_VIEW),
    condition_choices   = function() .pt_condition_choices(ctx$AVAILABLE_CONDITIONS),
    build_entries       = function(labels, col_map, settings) .pt_build_entries(labels, ctx$TRACK_ITEMS, col_map, settings),
    get(paste0(".pt_", name)))
}

need_fixture_data <- function()
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")

# A composer on the mouse fixture (samples Ba, LP, P6, P7), with an optional start view
# written before it is built.
helper_composer <- function(start_view = NULL) {
  proj <- tempfile("pt-helpers-"); dir.create(proj)
  seed_project_genome(proj)
  json <- file.path(proj, "project_settings.json")
  jsonlite::write_json(list(data_dir = FIXTURES$mouse$data_dir, samples = FIXTURES$mouse$samples,
                            gtf_path = GENOMES$mouse$gtf_path, orgdb = GENOMES$mouse$orgdb,
                            genome_label = GENOMES$mouse$genome_label),
                       json, auto_unbox = TRUE, null = "null")
  if (!is.null(start_view)) do.call(.pt_write_start_view, c(list(proj), start_view))
  list(proj = proj, json = json,
       app = suppressMessages(.pt_composer_app(json, .pt_read_settings(json, .pt_default_settings()))))
}

# Built once per file, on first use; removed when the tests end.
.helpers_cp <- new.env()
cached_composer <- function(name, start_view = NULL) {
  if (is.null(.helpers_cp[[name]])) {
    cp <- helper_composer(start_view)
    .helpers_cp[[name]] <- cp
    withr::defer(unlink(cp$proj, recursive = TRUE), envir = testthat::teardown_env())
  }
  .helpers_cp[[name]]
}
plain_cp  <- function() cached_composer("plain")
single_cp <- function() cached_composer("single", list(view_mode = "single", conditions = "P6", gene = "Elf5",
                                                       keys = c("bw_atac_P6", "genes", "hic_hic_P6")))

# ---- format_track_label ------------------------------------------------------------------

test_that("format_track_label bolds what precedes the first ' \u00B7 '", {
  need_fixture_data()
  f <- h("format_track_label", plain_cp())
  expect_identical(f("atac \u00B7 Ba"), "<strong>atac</strong> \u00B7 Ba")
  expect_identical(f("a \u00B7 b \u00B7 c"), "<strong>a</strong> \u00B7 b \u00B7 c")
  expect_identical(f("Gene annotation"), "<strong>Gene annotation</strong>")
  expect_identical(f("x \u00B7 "), "<strong>x</strong> \u00B7 ")
  expect_identical(f(c("atac \u00B7 Ba", "Transcripts")), c("<strong>atac</strong> \u00B7 Ba", "<strong>Transcripts</strong>"))
})

test_that("format_track_label leaves labels it cannot split alone", {
  need_fixture_data()
  f <- h("format_track_label", plain_cp())
  expect_identical(f(""), "")
  expect_identical(f(" \u00B7 x"), " \u00B7 x")       # nothing before the separator but its space
  expect_identical(f("atac\u00B7Ba"), "atac\u00B7Ba") # a middle dot without spaces is not the separator
  expect_identical(f(character(0)), character(0))
  expect_identical(f(NA_character_), NA_character_)
})

# ---- generate_button_ui and collapsible_panel ------------------------------------------

test_that("generate_button_ui is an action button when enabled, a disabled plain button otherwise", {
  need_fixture_data()
  f <- h("generate_button_ui", plain_cp())
  expect_identical(as.character(f("generate", TRUE)), paste0(
    '<button class="btn btn-default action-button btn-success btn-lg generate-btn" id="generate" type="button">',
    '<span class="action-label">Generate figure \u2192</span></button>'))
  # Disabled: no id, so Shiny never binds it, whatever id is given.
  expect_identical(as.character(f("anything", FALSE)),
                   '<button class="btn btn-success btn-lg generate-btn" disabled="disabled">Generate figure \u2192</button>')
})

test_that("collapsible_panel builds an open panel by default", {
  need_fixture_data()
  f <- h("collapsible_panel", plain_cp())
  expect_identical(as.character(f("sec", "Title", "body")), paste0(
    '<div class="panel panel-default plottrellis-section ">\n',
    '  <div class="panel-heading">\n',
    '    <h4 class="panel-title">\n',
    '      <a href="#sec-body" data-toggle="collapse">\n',
    '        <span class="collapse-arrow">\u25BE</span>\n',
    '         Title\n',
    '      </a>\n',
    '    </h4>\n',
    '  </div>\n',
    '  <div id="sec-body" class="collapse in">\n',
    '    <div class="panel-body">body</div>\n',
    '  </div>\n',
    '</div>'))
})

test_that("collapsible_panel collapsed, with a wrapper id and style", {
  need_fixture_data()
  f <- h("collapsible_panel", plain_cp())
  html <- as.character(f("out", "Output", tags$p("a"), tags$p("b"), wrapper_id = "w",
                         wrapper_style = "display:none;", collapsed = TRUE))
  expect_match(html, '^<div id="w" class="panel panel-default plottrellis-section collapsed-panel" style="display:none;">')
  expect_match(html, '<span class="collapse-arrow">\u25B8</span>', fixed = TRUE)
  expect_match(html, '<div id="out-body" class="collapse ">', fixed = TRUE)
  expect_match(html, '<div class="panel-body">\n      <p>a</p>\n      <p>b</p>\n    </div>', fixed = TRUE)
})

# ---- start view helpers ----------------------------------------------------------------

test_that("start_condition is '' with no start view, else the view's i-th condition or ''", {
  need_fixture_data()
  expect_identical(h("start_condition", plain_cp())(1), "")
  expect_identical(h("start_condition", plain_cp())(2), "")
  expect_identical(h("start_condition", single_cp())(1), "P6")
  expect_identical(h("start_condition", single_cp())(2), "")
})

test_that("default_figure_keys is Gene annotation unless the conditions are the start view's", {
  need_fixture_data()
  for (conds in list("Ba", c("Ba", "LP"), NULL, character(0)))
    expect_identical(h("default_figure_keys", plain_cp())(conds), "genes", info = deparse(conds))
  f <- h("default_figure_keys", single_cp())
  expect_identical(f("P6"), c("bw_atac_P6", "genes", "hic_hic_P6"))
  expect_identical(f("Ba"), "genes")
  expect_identical(f(c("P6", "Ba")), "genes")
})

# ---- sort_chrs -------------------------------------------------------------------------

test_that("sort_chrs puts numbered chromosomes in numeric order, then X, Y, M, then the rest", {
  need_fixture_data()
  f <- h("sort_chrs", plain_cp())
  expect_identical(f(c("chrM", "chr10", "chr2", "chrY", "chr1", "chrX", "chrUn_x", "chr1_random",
                       "scaffold_1", "chr2_alt", "MT", "1")),
                   c("chr1", "chr2", "chr10", "chrX", "chrY", "chrM",
                     "1", "MT", "chr1_random", "chr2_alt", "chrUn_x", "scaffold_1"))  # C collation
})

test_that("sort_chrs edge cases: empty, one, missing sex chromosomes, zero padding, repeats", {
  need_fixture_data()
  f <- h("sort_chrs", plain_cp())
  expect_identical(f(character(0)), character(0))
  expect_identical(f("chrX"), "chrX")
  expect_identical(f(c("chr3", "chrM", "chr1")), c("chr1", "chr3", "chrM"))
  expect_identical(f(c("chr2", "chr01", "chr1")), c("chr01", "chr1", "chr2"))  # chr01 is number 1
  expect_identical(f(c("chr2", "chr2", "chr1")), c("chr1", "chr2", "chr2"))
  # Repeated X/Y/M collapse to one each.
  expect_identical(f(c("chrX", "chrX", "chr1")), c("chr1", "chrX"))
  # Ensembl-style names (no "chr") are not main chromosomes.
  expect_identical(f(c("2", "10", "X")), c("10", "2", "X"))
})

# ---- fmt_num, format_bp, format_coord --------------------------------------------------

test_that("fmt_num rounds to the given digits and drops trailing zeros", {
  need_fixture_data()
  f <- h("fmt_num", plain_cp())
  expect_identical(f(1), "1")
  expect_identical(f(1.5), "1.5")
  expect_identical(f(2.10), "2.1")
  expect_identical(f(1.234), "1.23")
  expect_identical(f(1.999), "2")
  expect_identical(f(1.005), "1")   # 1.005 is just below 1.005 in binary
  expect_identical(f(0), "0")
  expect_identical(f(1e-9), "0")
  expect_identical(f(-1.5), "-1.5")
  expect_identical(f(1234567.891), "1234567.89")
  expect_identical(f(1.2345, 3), "1.234")
  expect_identical(f(c(1.5, 2)), c("1.5", "2"))
  expect_identical(f(NA), " NA")
})

test_that("fmt_num with digits = 0 strips an integer's own zeros", {
  need_fixture_data()
  f <- h("fmt_num", plain_cp())
  # current behaviour, possibly a bug: with no decimal point, the trailing-zero strip eats
  # the integer's zeros (no caller in the app passes digits = 0).
  expect_identical(f(100, 0), "1")
  expect_identical(f(0, 0), "")
  expect_identical(f(1.5, 0), "2")
})

test_that("format_bp switches unit at 1 kb and 1 Mb, two decimals, no space", {
  need_fixture_data()
  f <- h("format_bp", plain_cp())
  expect_identical(f(0), "0bp")
  expect_identical(f(999), "999bp")
  expect_identical(f(999.5), "999.5bp")
  expect_identical(f(1000), "1Kb")
  expect_identical(f(1500), "1.5Kb")
  expect_identical(f(45000), "45Kb")
  expect_identical(f(1e6), "1Mb")
  expect_identical(f(182113224), "182.11Mb")
  expect_error(f(NA), "missing value")
})

test_that("format_bp just under 1 Mb rounds up to 1000Kb", {
  need_fixture_data()
  # current behaviour, possibly a bug: the unit is picked before rounding, so 999,999 bp
  # reads "1000Kb" rather than "1Mb".
  expect_identical(h("format_bp", plain_cp())(999999), "1000Kb")
})

test_that("format_coord uses three decimals and a spaced bp unit", {
  need_fixture_data()
  f <- h("format_coord", plain_cp())
  expect_identical(f(1L), "1 bp")
  expect_identical(f(999), "999 bp")
  expect_identical(f(1000), "1Kb")
  expect_identical(f(16299), "16.299Kb")
  expect_identical(f(999999), "999.999Kb")
  expect_identical(f(1e6), "1Mb")
  expect_identical(f(1234567), "1.235Mb")
  expect_identical(f(182113224), "182.113Mb")
  expect_identical(f(-5), "-5 bp")
})

# ---- parse_bp ----------------------------------------------------------------------------

test_that("parse_bp accepts plain digits and kb/Mb in any case, with or without a space", {
  need_fixture_data()
  f <- h("parse_bp", plain_cp())
  ok <- c("17000000" = 17000000L, "17Mb" = 17000000L, "17.56Mb" = 17560000L, "15kb" = 15000L,
          "150KB" = 150000L, "1.5 Mb" = 1500000L, "12 kb" = 12000L, "  42  " = 42L, "0" = 0L,
          ".5kb" = 500L, "17mB" = 17000000L, "1.5 MB" = 1500000L, "1.2345678Mb" = 1234568L,
          "15.5kb" = 15500L, "007" = 7L)
  for (s in names(ok)) expect_identical(f(s), ok[[s]], info = s)
})

test_that("parse_bp gives NA for anything else", {
  need_fixture_data()
  f <- h("parse_bp", plain_cp())
  for (s in c("1,000", "-5", "abc", "", "   ", "1e6", "1.5", "5 gb", "5bp", "1.  5kb", "Mb", "1.5.5kb",
              "5 k b", "1 000"))
    expect_identical(f(s), NA_integer_, info = s)
  expect_identical(f(NA), NA_integer_)
})

test_that("parse_bp past the integer range warns and gives NA; NULL and empty error", {
  need_fixture_data()
  f <- h("parse_bp", plain_cp())
  expect_warning(v <- f("3000000000"), "integer range")
  expect_identical(v, NA_integer_)
  expect_warning(v <- f("3000Mb"), "integer range")
  expect_identical(v, NA_integer_)
  expect_error(f(NULL), "argument is of length zero")
  expect_error(f(character(0)), "argument is of length zero")
})

# ---- build_entries -------------------------------------------------------------------------

test_that("build_entries turns picker keys into render entries and drops unknown keys", {
  need_fixture_data()
  f <- h("build_entries", plain_cp())
  e <- f(c("genes", "nope", "bw_atac_Ba", "genes"), c(Ba = "#E41A1C"), list())
  expect_length(e, 3)
  expect_identical(e[[1]], list(viz = "genes", height = NULL, data = NULL))
  expect_identical(e[[3]], e[[1]])
  expect_identical(e[[2]][c("viz", "sample", "colour", "name")],
                   list(viz = "bw", sample = "Ba", colour = "#E41A1C", name = "atac"))
  expect_identical(basename(e[[2]]$path), "Ba_ATAC.bw")
  expect_identical(f(character(0), list(), list()), list())
})

test_that("build_entries gives comparison tracks the role colours, defaulting when unset", {
  need_fixture_data()
  f <- h("build_entries", plain_cp())
  key <- "bw_diff_stacked_atac_Ba_vs_LP"
  e <- f(key, list(), list(reference_colour = "#111111", comparison_colour = "#222222"))[[1]]
  expect_identical(e[c("viz", "control_label", "treatment_label", "control_colour", "treatment_colour")],
                   list(viz = "bw_diff_stacked", control_label = "Ba", treatment_label = "LP",
                        control_colour = "#111111", treatment_colour = "#222222"))
  e <- f(key, list(), list())[[1]]
  expect_identical(c(e$control_colour, e$treatment_colour), c("#4D4D4D", "#56B4E9"))
})

# ---- track_labels and sort_keys_by_label -------------------------------------------------

test_that("track_labels gives each key its formatted label as HTML, named by key", {
  need_fixture_data()
  f <- h("track_labels", plain_cp())
  l <- f(c("bw_atac_Ba", "genes"))
  expect_identical(names(l), c("bw_atac_Ba", "genes"))
  expect_s3_class(l$genes, "html")
  expect_identical(as.character(l$bw_atac_Ba), "<strong>atac</strong> \u00B7 Ba")
  expect_identical(as.character(l$genes), "<strong>Gene annotation</strong>")
  expect_identical(as.character(f("tpm_usage_LP")[[1]]), "<strong>Transcript expression</strong> \u00B7 LP")
  expect_identical(f(character(0)), setNames(list(), character(0)))
  # An unknown key keeps its name, with an empty label.
  expect_identical(as.character(f("nope")$nope), "")
})

test_that("sort_keys_by_label groups by assay, then sample, ignoring case", {
  need_fixture_data()
  f <- h("sort_keys_by_label", plain_cp())
  expect_identical(f(c("tpm_usage_Ba", "bw_atac_P6", "genes", "bw_atac_Ba", "hic_hic_LP", "transcripts",
                       "transcripts_all", "protein_domains", "transcripts_annotated", "junctions_Junctions_P6")),
                   c("bw_atac_Ba", "bw_atac_P6", "protein_domains", "genes", "hic_hic_LP",
                     "junctions_Junctions_P6", "tpm_usage_Ba", "transcripts_all", "transcripts",
                     "transcripts_annotated"))
  # Comparison labels sort by their whole sample part.
  expect_identical(f(c("bw_diff_stacked_atac_P6_vs_Ba", "bw_diff_stacked_atac_Ba_vs_LP", "bw_atac_LP")),
                   c("bw_atac_LP", "bw_diff_stacked_atac_P6_vs_Ba", "bw_diff_stacked_atac_Ba_vs_LP"))
})

test_that("sort_keys_by_label edge cases: empty, one, ties keep their order, unknown errors", {
  need_fixture_data()
  f <- h("sort_keys_by_label", plain_cp())
  expect_identical(f(character(0)), character(0))
  expect_identical(f("genes"), "genes")
  expect_identical(f(c("bw_atac_P6", "bw_atac_Ba", "bw_atac_P6")), c("bw_atac_Ba", "bw_atac_P6", "bw_atac_P6"))
  expect_identical(f(c("genes", "genes")), c("genes", "genes"))
  expect_error(f("nope"), "argument is of length zero")
})

# ---- track_split and condition_choices ---------------------------------------------------

test_that("track_split for one condition: its tracks and the sample-free ones; Figure holds Gene annotation", {
  need_fixture_data()
  s <- h("track_split", plain_cp())("Ba")
  expect_identical(names(s), c("available", "figure"))
  expect_identical(names(s$available),
                   c("bw_atac_Ba", "protein_domains", "genes", "hic_hic_Ba", "junctions_Junctions_Ba",
                     "tpm_usage_Ba", "transcripts_all", "transcripts", "transcripts_annotated"))
  expect_identical(names(s$figure), "genes")
  expect_identical(as.character(s$figure$genes), "<strong>Gene annotation</strong>")
})

test_that("track_split for a pair adds that ordered pair's comparisons only", {
  need_fixture_data()
  f <- h("track_split", plain_cp())
  s <- f(c("Ba", "P6"))
  expect_identical(names(s$available),
                   c("bw_atac_Ba", "bw_atac_P6", "bw_diff_stacked_atac_Ba_vs_P6", "protein_domains", "genes",
                     "hic_hic_Ba", "hic_hic_P6", "junctions_Junctions_Ba", "junctions_Junctions_P6",
                     "junctions_mirrored_diff_Junctions_Ba_vs_P6", "tpm_usage_Ba", "transcripts_all",
                     "transcripts", "transcripts_annotated"))
  expect_identical(grep("_vs_", names(f(c("P6", "Ba"))$available), value = TRUE),
                   c("bw_diff_stacked_atac_P6_vs_Ba", "junctions_mirrored_diff_Junctions_P6_vs_Ba"))
  # LP has no junctions: its pair with Ba has the coverage comparison only.
  expect_identical(grep("_vs_", names(f(c("Ba", "LP"))$available), value = TRUE), "bw_diff_stacked_atac_Ba_vs_LP")
})

test_that("track_split with no, unknown or three conditions", {
  need_fixture_data()
  f <- h("track_split", plain_cp())
  free <- c("protein_domains", "genes", "transcripts_all", "transcripts", "transcripts_annotated")
  expect_identical(names(f(character(0))$available), free)
  expect_identical(names(f(character(0))$figure), "genes")
  expect_identical(names(f("Nobody")$available), free)
  three <- names(f(c("Ba", "LP", "P6"))$available)
  expect_false(any(grepl("_vs_", three)))  # comparisons are for exactly two
  expect_true(all(c("bw_atac_Ba", "bw_atac_LP", "bw_atac_P6") %in% three))
})

test_that("track_split takes Figure tracks from a matching start view; Gene annotation stays in Available", {
  need_fixture_data()
  f <- h("track_split", single_cp())
  s <- f("P6")
  expect_identical(names(s$figure), c("bw_atac_P6", "genes", "hic_hic_P6"))
  expect_identical(names(s$available),
                   c("protein_domains", "genes", "junctions_Junctions_P6", "transcripts_all",
                     "transcripts", "transcripts_annotated"))
  expect_identical(names(f("Ba")$figure), "genes")
})

test_that("condition_choices is a blank placeholder then every sample", {
  need_fixture_data()
  expect_identical(h("condition_choices", plain_cp())(),
                   c("(select)" = "", Ba = "Ba", LP = "LP", P6 = "P6", P7 = "P7"))
})

# ---- inline_svg ------------------------------------------------------------------------------

test_that("inline_svg sizes the first <svg> to its box and wraps it in a capped, centred div", {
  need_fixture_data()
  f <- h("inline_svg", plain_cp())
  svg <- tempfile(fileext = ".svg")
  on.exit(unlink(svg), add = TRUE)
  writeLines(c("<svg width='1'>", "<g/>", "</svg>"), svg)
  expect_identical(as.character(f(svg)), paste0(
    '<div style="width:100%; max-width:1400px; margin:0 auto; overflow-x:auto;">',
    "<svg style='width:100%;height:auto;' width='1'>\n<g/>\n</svg></div>"))
  writeLines(c("<?xml version='1.0'?>", "<svg a='1'><svg b='2'/></svg>"), svg)
  expect_identical(as.character(f(svg)), paste0(
    '<div style="width:100%; max-width:1400px; margin:0 auto; overflow-x:auto;">',
    "<?xml version='1.0'?>\n<svg style='width:100%;height:auto;' a='1'><svg b='2'/></svg></div>"))
})

test_that("inline_svg of an empty file is an empty box; a missing file errors", {
  need_fixture_data()
  f <- h("inline_svg", plain_cp())
  empty <- tempfile(fileext = ".svg")
  on.exit(unlink(empty), add = TRUE)
  file.create(empty)
  expect_identical(as.character(f(empty)),
                   '<div style="width:100%; max-width:1400px; margin:0 auto; overflow-x:auto;"></div>')
  expect_error(suppressWarnings(f(tempfile(fileext = ".svg"))), "cannot open")
})
