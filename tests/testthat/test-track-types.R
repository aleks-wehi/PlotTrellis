# Every list of track types the app, the code API and the fixture helpers use is derived
# from the modules' own info_<key>() declarations (see the viz-module-contract in
# R/render_all.R). The literals below pin those lists' values, and the derived
# values must stay equal to them (in order where the
# order decides something: detection order sets pt_tracks()' order and which of two
# same-named types keeps the plain name). Only the track types below are pinned, so a new
# track type needs no change here; it must not change what these give for the old ones.
EXISTING <- c("bw", "hic", "bed", "junctions", "tpm_usage", "genes", "protein_domains", "transcripts_all",
              "transcripts", "transcripts_annotated", "bw_diff_stacked", "junctions_mirrored_diff",
              "junctions_mirrored_diff_tpm_scaled")

# The declarations of the pinned types, in detection order.
infos <- function() {
  all <- .pt_viz_infos()
  all[names(all) %in% EXISTING]
}
keys_where <- function(pred) names(Filter(pred, infos()))
# The keys and the API names of the pinned types' items (added types' items left out).
pinned_keys <- function(items) names(Filter(function(it) it$type %in% EXISTING, items))
pinned_table <- function(items) {
  tab <- .pt_track_table(items)
  tab[vapply(tab$key, function(k) items[[k]]$type, "") %in% EXISTING, ]
}

test_that("the declared track types, by kind, in detection order", {
  expect_identical(keys_where(function(i) i$kind == "sample"), c("bw", "hic", "bed", "junctions", "tpm_usage"))
  expect_identical(keys_where(function(i) i$kind == "fixed"),
                   c("genes", "protein_domains", "transcripts_all", "transcripts", "transcripts_annotated"))
  expect_identical(keys_where(function(i) i$kind == "comparison"),
                   c("bw_diff_stacked", "junctions_mirrored_diff", "junctions_mirrored_diff_tpm_scaled"))
  # Every module is offered but the genome axis, which render_all() adds itself.
  expect_identical(setdiff(.pt_viz_keys(), names(.pt_viz_infos())), "genome_axis")
  expect_false(.pt_viz_info("genome_axis")$offered)
  expect_identical(vapply(infos(), `[[`, "", "data_type")[c("bw", "hic", "bed", "junctions", "tpm_usage",
                                                            "bw_diff_stacked", "junctions_mirrored_diff",
                                                            "junctions_mirrored_diff_tpm_scaled")],
                   c(bw = "bw", hic = "hic", bed = "bed", junctions = "paired-end", tpm_usage = "quantdir",
                     bw_diff_stacked = "bw", junctions_mirrored_diff = "paired-end",
                     junctions_mirrored_diff_tpm_scaled = "paired-end"))
})

test_that("scale-bar, row-panel and gene-required types are the hand-kept lists", {
  expect_setequal(keys_where(function(i) i$scale_bar), c("bw", "hic", "bw_diff_stacked"))  # SCALE_BAR_TYPES
  expect_identical(keys_where(function(i) i$row_panel), "tpm_usage")  # ROW_PANEL_TYPES
  expect_setequal(.pt_gene_required_types(infos()),  # GENE_REQUIRED_TYPES
                  c("transcripts_all", "transcripts", "transcripts_annotated", "protein_domains", "tpm_usage",
                    "junctions_mirrored_diff_tpm_scaled"))
  # A type without a declaration takes neither.
  expect_false(.pt_viz_info("genome_axis")$scale_bar)
  expect_false(.pt_viz_info("genome_axis")$row_panel)
})

test_that("sample-free tracks keep their labels and order (FIXED_ITEMS)", {
  fixed <- Filter(function(i) i$kind == "fixed", infos())
  expect_identical(lapply(fixed, `[[`, "label"), list(
    genes                 = "Gene annotation",
    protein_domains       = "Canonical protein annotation",
    transcripts_all       = "Transcripts",
    transcripts           = "Transcripts (filtered)",
    transcripts_annotated = "Transcripts annotated (filtered)"
  ))
})

test_that("the app's filter buttons and type-to-filter map are the hand-kept ones", {
  map <- .pt_type_filters(infos())
  expect_setequal(names(map), EXISTING)
  # COMPARISON_VIZ_TYPES_JS and ANNOTATION_VIZ_TYPES_JS; every other type was its own filter.
  expect_setequal(names(Filter(function(f) f == "comparison", map)),
                  c("bw_diff_stacked", "junctions_mirrored_diff", "junctions_mirrored_diff_tpm_scaled"))
  expect_setequal(names(Filter(function(f) f == "annotation", map)),
                  c("genes", "transcripts", "transcripts_all", "transcripts_annotated", "protein_domains"))
  other <- Filter(function(f) !f %in% c("comparison", "annotation"), map)
  expect_identical(unlist(other), setNames(names(other), names(other)))
  # The buttons, written out by hand.
  old <- tagList(
    conditionalPanel("input.view_mode === 'two'",
      tags$button("Comparisons", class = "btn btn-lg btn-default filter-btn", `data-filter` = "comparison")
    ),
    tags$button("Annotations", class = "btn btn-lg btn-default filter-btn", `data-filter` = "annotation"),
    tags$button("BED",         class = "btn btn-lg btn-default filter-btn", `data-filter` = "bed"),
    tags$button("Coverage",    class = "btn btn-lg btn-default filter-btn", `data-filter` = "bw"),
    tags$button("HiC",         class = "btn btn-lg btn-default filter-btn", `data-filter` = "hic"),
    tags$button("Junctions",   class = "btn btn-lg btn-default filter-btn", `data-filter` = "junctions"),
    tags$button("TPMs",        class = "btn btn-lg btn-default filter-btn", `data-filter` = "tpm_usage")
  )
  expect_identical(as.character(tagList(.pt_filter_buttons(infos()))), as.character(old))
})

# Two samples with one of every data type (named as in the human test data), and two whose
# data of four types share one name. Paths need not exist.
fake_rows <- function() lapply(c("A", "B"), function(s) list(display = s, data = list(
  list(type = "bw", name = "atac", path = paste0(s, "_atac.bw")),
  list(type = "bw", name = "h3k27ac", path = paste0(s, "_h3k27ac.bw")),
  list(type = "hic", name = "hic", path = paste0(s, ".hic")),
  list(type = "bed", name = "peaks", path = paste0(s, "_peaks.bed")),
  list(type = "paired-end", name = "Junctions", path = paste0(s, "_junctions.bed")),
  list(type = "quantdir", name = "Quant", path = paste0(s, "_quant.sf")))))
clash_rows <- function() lapply(c("A", "B"), function(s) list(display = s, data = list(
  list(type = "bed", name = "atac", path = paste0(s, "_atac.bed")),
  list(type = "paired-end", name = "atac", path = paste0(s, "_atac_junctions.bed")),
  list(type = "bw", name = "atac", path = paste0(s, "_atac.bw")),
  list(type = "hic", name = "atac", path = paste0(s, "_atac.hic")))))

test_that("detection offers the same track items, in the same order", {
  items <- build_track_items(fake_rows())
  expect_identical(pinned_keys(items), c(
    "bw_atac_A", "bw_h3k27ac_A", "hic_hic_A", "bed_peaks_A", "junctions_Junctions_A", "tpm_usage_A",
    "bw_atac_B", "bw_h3k27ac_B", "hic_hic_B", "bed_peaks_B", "junctions_Junctions_B", "tpm_usage_B",
    "genes", "protein_domains", "transcripts_all", "transcripts", "transcripts_annotated"))
  expect_identical(unname(vapply(items[pinned_keys(items)[1:6]], `[[`, "", "label")), c(
    "atac \u00B7 A", "h3k27ac \u00B7 A", "hic \u00B7 A", "peaks \u00B7 A", "Junctions \u00B7 A",
    "Transcript expression \u00B7 A"))
  expect_identical(items$bw_atac_A, list(key = "bw_atac_A", label = "atac \u00B7 A", type = "bw", sample = "A",
    file = "A_atac.bw", colour_key = "A", name = "atac", broken = TRUE))
  expect_identical(items$junctions_Junctions_A, list(key = "junctions_Junctions_A", label = "Junctions \u00B7 A",
    type = "junctions", sample = "A", file = NA_character_, colour_key = "A", name = "Junctions", broken = TRUE,
    junction_path = "A_junctions.bed"))
  # Field order aside (quant_path now comes first), as before.
  expect_identical(items$tpm_usage_A, list(key = "tpm_usage_A", label = "Transcript expression \u00B7 A",
    type = "tpm_usage", sample = "A", file = NA_character_, colour_key = "A", quant_path = "A_quant.sf",
    broken = TRUE))
  expect_identical(items$genes, list(key = "genes", label = "Gene annotation", type = "genes",
    sample = NA_character_, file = NA_character_, colour_key = NA_character_))

  cmp <- build_comparison_track_items(fake_rows())
  expect_identical(pinned_keys(cmp), c(
    "bw_diff_stacked_atac_A_vs_B", "bw_diff_stacked_h3k27ac_A_vs_B", "junctions_mirrored_diff_Junctions_A_vs_B",
    "junctions_mirrored_diff_tpm_scaled_Junctions_A_vs_B", "bw_diff_stacked_atac_B_vs_A",
    "bw_diff_stacked_h3k27ac_B_vs_A", "junctions_mirrored_diff_Junctions_B_vs_A",
    "junctions_mirrored_diff_tpm_scaled_Junctions_B_vs_A"))
  expect_identical(unname(vapply(cmp[pinned_keys(cmp)[1:4]], `[[`, "", "label")), c(
    "atac diff (stacked) \u00B7 B (Cmp) vs A (Ref)", "h3k27ac diff (stacked) \u00B7 B (Cmp) vs A (Ref)",
    "Junctions mirrored diff \u00B7 B (Cmp) vs A (Ref)", "Junctions mirrored diff (TPM-scaled) \u00B7 B (Cmp) vs A (Ref)"))
  expect_identical(cmp$junctions_mirrored_diff_tpm_scaled_Junctions_A_vs_B, list(
    key = "junctions_mirrored_diff_tpm_scaled_Junctions_A_vs_B",
    label = "Junctions mirrored diff (TPM-scaled) \u00B7 B (Cmp) vs A (Ref)",
    type = "junctions_mirrored_diff_tpm_scaled", sample = NA_character_, file = NA_character_,
    colour_key = NA_character_, name = "Junctions", samples = c("A", "B"), control_label = "A",
    treatment_label = "B", control_path = "A_junctions.bed", treatment_path = "B_junctions.bed",
    control_quant = "A_quant.sf", treatment_quant = "B_quant.sf"))
  # The TPM-scaled comparison needs a quantification on both sides.
  no_quant <- lapply(fake_rows(), function(r) { r$data <- r$data[1:5]; r })
  expect_false(any(grepl("tpm_scaled", names(build_comparison_track_items(no_quant)))))
})

test_that("API track names are the same, the comparison suffixes and transcript_expression included", {
  rows <- fake_rows()
  tab <- pinned_table(c(build_track_items(rows), build_comparison_track_items(rows)))
  expect_identical(tab$track, c(
    "genes", "protein_domains", "transcripts_all", "transcripts", "transcripts_annotated",
    "atac", "h3k27ac", "hic", "peaks", "Junctions", "transcript_expression",
    "atac", "h3k27ac", "hic", "peaks", "Junctions", "transcript_expression",
    "atac_diff", "h3k27ac_diff", "Junctions_diff", "Junctions_diff_tpm",
    "atac_diff", "h3k27ac_diff", "Junctions_diff", "Junctions_diff_tpm"))
  # One name, four data types: the first detected keeps it.
  rows <- clash_rows()
  tab <- pinned_table(c(build_track_items(rows), build_comparison_track_items(rows)))
  expect_identical(tab$track, c(
    "genes", "protein_domains", "transcripts_all", "transcripts", "transcripts_annotated",
    "atac", "atac_hic", "atac_bed", "atac_junctions", "atac", "atac_hic", "atac_bed", "atac_junctions",
    "atac_diff", "atac_diff_junctions_mirrored_diff", "atac_diff", "atac_diff_junctions_mirrored_diff"))
})

test_that("figure entries are the same for one track of every type", {
  rows <- fake_rows()
  items <- c(build_track_items(rows), build_comparison_track_items(rows))
  keys <- c("bw_atac_A", "hic_hic_A", "bed_peaks_A", "junctions_Junctions_A", "tpm_usage_A", "genes",
            "protein_domains", "transcripts_all", "transcripts", "transcripts_annotated",
            "bw_diff_stacked_atac_A_vs_B", "junctions_mirrored_diff_Junctions_A_vs_B",
            "junctions_mirrored_diff_tpm_scaled_Junctions_A_vs_B")
  settings <- list(hic_style = "triangle", reference_colour = "#111111", comparison_colour = "#222222")
  entries <- .pt_build_entries(keys, items, c(A = "#E41A1C", B = "#377EB8"), settings)
  fixed <- function(viz) list(viz = viz, height = NULL, data = NULL)
  cmp <- function(viz, name, ext) list(viz = viz, name = name, control_label = "A", treatment_label = "B",
    control_path = paste0("A", ext), treatment_path = paste0("B", ext), control_colour = "#111111",
    treatment_colour = "#222222", height = NULL, data = NULL)
  expect_identical(entries, list(
    list(viz = "bw", height = NULL, data = NULL, sample = "A", colour = "#E41A1C", name = "atac", path = "A_atac.bw"),
    list(viz = "hic", height = NULL, data = NULL, sample = "A", colour = "#E41A1C", name = "hic", path = "A.hic",
         hic_style = "triangle"),
    list(viz = "bed", height = NULL, data = NULL, sample = "A", colour = "#E41A1C", name = "peaks",
         path = "A_peaks.bed"),
    list(viz = "junctions", height = NULL, data = NULL, sample = "A", colour = "#E41A1C", name = "Junctions",
         path = "A_junctions.bed"),
    list(viz = "tpm_usage", height = NULL, data = NULL, sample = "A", colour = "#E41A1C",
         quant_path = "A_quant.sf"),
    fixed("genes"), fixed("protein_domains"), fixed("transcripts_all"), fixed("transcripts"),
    fixed("transcripts_annotated"),
    cmp("bw_diff_stacked", "atac", "_atac.bw"),
    cmp("junctions_mirrored_diff", "Junctions", "_junctions.bed"),
    c(cmp("junctions_mirrored_diff_tpm_scaled", "Junctions", "_junctions.bed"),
      list(control_quant = "A_quant.sf", treatment_quant = "B_quant.sf"))
  ))
  # No hic_style setting: none in the entry, as before.
  expect_false("hic_style" %in% names(.pt_build_entries("hic_hic_A", items, c(A = "#E41A1C"), list())[[1]]))
})

test_that("the test datasets offer the same tracks, in the same order", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  mouse <- test_rows("mouse")
  expect_identical(pinned_keys(build_track_items(mouse)), c(
    "bw_atac_Ba", "hic_hic_Ba", "junctions_Junctions_Ba", "tpm_usage_Ba", "bw_atac_LP", "hic_hic_LP",
    "tpm_usage_LP", "bw_atac_P6", "hic_hic_P6", "junctions_Junctions_P6", "bw_atac_P7", "hic_hic_P7",
    "junctions_Junctions_P7", "genes", "protein_domains", "transcripts_all", "transcripts",
    "transcripts_annotated"))
  expect_identical(pinned_keys(build_comparison_track_items(mouse)), c(
    "bw_diff_stacked_atac_Ba_vs_LP", "bw_diff_stacked_atac_Ba_vs_P6", "junctions_mirrored_diff_Junctions_Ba_vs_P6",
    "bw_diff_stacked_atac_Ba_vs_P7", "junctions_mirrored_diff_Junctions_Ba_vs_P7", "bw_diff_stacked_atac_LP_vs_Ba",
    "bw_diff_stacked_atac_LP_vs_P6", "bw_diff_stacked_atac_LP_vs_P7", "bw_diff_stacked_atac_P6_vs_Ba",
    "junctions_mirrored_diff_Junctions_P6_vs_Ba", "bw_diff_stacked_atac_P6_vs_LP", "bw_diff_stacked_atac_P6_vs_P7",
    "junctions_mirrored_diff_Junctions_P6_vs_P7", "bw_diff_stacked_atac_P7_vs_Ba",
    "junctions_mirrored_diff_Junctions_P7_vs_Ba", "bw_diff_stacked_atac_P7_vs_LP", "bw_diff_stacked_atac_P7_vs_P6",
    "junctions_mirrored_diff_Junctions_P7_vs_P6"))
  human <- test_rows("human")
  expect_identical(pinned_keys(build_track_items(human)), c(
    "bw_atac_GM12878", "bw_h3k27ac_GM12878", "hic_hic_GM12878", "bed_peaks_GM12878", "junctions_Junctions_GM12878",
    "tpm_usage_GM12878", "bw_atac_K562", "bw_h3k27ac_K562", "hic_hic_K562", "bed_peaks_K562",
    "junctions_Junctions_K562", "tpm_usage_K562", "genes", "protein_domains", "transcripts_all", "transcripts",
    "transcripts_annotated"))
  expect_identical(pinned_keys(build_comparison_track_items(human)), c(
    "bw_diff_stacked_atac_GM12878_vs_K562", "bw_diff_stacked_h3k27ac_GM12878_vs_K562",
    "junctions_mirrored_diff_Junctions_GM12878_vs_K562", "junctions_mirrored_diff_tpm_scaled_Junctions_GM12878_vs_K562",
    "bw_diff_stacked_atac_K562_vs_GM12878", "bw_diff_stacked_h3k27ac_K562_vs_GM12878",
    "junctions_mirrored_diff_Junctions_K562_vs_GM12878", "junctions_mirrored_diff_tpm_scaled_Junctions_K562_vs_GM12878"))
})

# pt_tracks() of a project, the pinned types' rows only.
pinned_tracks <- function(proj) {
  tab <- pt_tracks(proj)[vapply(proj$tracks$key, function(k) proj$items[[k]]$type, "") %in% EXISTING, ]
  rownames(tab) <- NULL
  tab
}

test_that("pt_tracks() of both bundled examples is the same", {
  skip_if_not_installed("org.Mm.eg.db")
  skip_if_not_installed("org.Hs.eg.db")
  mammary <- pinned_tracks(suppressWarnings(pt_project(file.path(tempdir(), "pt-test-script-mammary"),
                                                       pt_example("mammary"))))
  expect_identical(mammary, data.frame(
    track = c("genes", "protein_domains", "transcripts_all", "transcripts", "transcripts_annotated", "atac",
              "transcript_expression", "atac", "transcript_expression", "atac_diff", "atac_diff"),
    sample = c(NA, NA, NA, NA, NA, "Basal", "Basal", "LP", "LP", "LP vs Basal", "Basal vs LP"),
    kind = c(rep("fixed", 5), rep("sample", 4), rep("comparison", 2)),
    label = c("Gene annotation", "Canonical protein annotation", "Transcripts", "Transcripts (filtered)",
              "Transcripts annotated (filtered)", "atac \u00B7 Basal", "Transcript expression \u00B7 Basal",
              "atac \u00B7 LP", "Transcript expression \u00B7 LP", "atac diff (stacked) \u00B7 LP (Cmp) vs Basal (Ref)",
              "atac diff (stacked) \u00B7 Basal (Cmp) vs LP (Ref)"),
    stringsAsFactors = FALSE))
  human <- pinned_tracks(suppressWarnings(pt_project(file.path(tempdir(), "pt-test-example-human"), pt_example("human"))))
  expect_identical(human, data.frame(
    track = c("genes", "protein_domains", "transcripts_all", "transcripts", "transcripts_annotated", "atac",
              "junctions", "transcript_expression", "atac", "junctions", "transcript_expression", "atac_diff",
              "junctions_diff", "junctions_diff_tpm", "atac_diff", "junctions_diff", "junctions_diff_tpm"),
    sample = c(NA, NA, NA, NA, NA, "GM12878", "GM12878", "GM12878", "K562", "K562", "K562", "K562 vs GM12878",
               "K562 vs GM12878", "K562 vs GM12878", "GM12878 vs K562", "GM12878 vs K562", "GM12878 vs K562"),
    kind = c(rep("fixed", 5), rep("sample", 6), rep("comparison", 6)),
    label = c("Gene annotation", "Canonical protein annotation", "Transcripts", "Transcripts (filtered)",
              "Transcripts annotated (filtered)", "atac \u00B7 GM12878", "junctions \u00B7 GM12878",
              "Transcript expression \u00B7 GM12878", "atac \u00B7 K562", "junctions \u00B7 K562",
              "Transcript expression \u00B7 K562", "atac diff (stacked) \u00B7 K562 (Cmp) vs GM12878 (Ref)",
              "junctions mirrored diff \u00B7 K562 (Cmp) vs GM12878 (Ref)",
              "junctions mirrored diff (TPM-scaled) \u00B7 K562 (Cmp) vs GM12878 (Ref)",
              "atac diff (stacked) \u00B7 GM12878 (Cmp) vs K562 (Ref)",
              "junctions mirrored diff \u00B7 GM12878 (Cmp) vs K562 (Ref)",
              "junctions mirrored diff (TPM-scaled) \u00B7 GM12878 (Cmp) vs K562 (Ref)"),
    stringsAsFactors = FALSE))
})

# ---- Registered wrongly: loud, not silent --------------------------------------

# .pt_module_fn() with extra modules, defined here: `fns` maps a function name to it.
with_modules <- function(fns, env = parent.frame()) {
  real <- .pt_module_fn
  local_mocked_bindings(.pt_module_fn = function(prefix, key) {
    fn <- fns[[paste0(prefix, "_", key)]]
    if (is.null(fn)) real(prefix, key) else fn
  }, .env = env)
}

test_that("a malformed declaration stops with an error naming the module and the field", {
  info_of <- function(info) {
    with_modules(list(info_toy = function() info))
    .pt_viz_info("toy")
  }
  expect_error(info_of(list(kind = "per-sample", data_type = "bw")), "info_toy() field 'kind'", fixed = TRUE)
  expect_error(info_of(list(data_type = "bw")), "info_toy() field 'kind'", fixed = TRUE)
  expect_error(info_of(list(kind = "sample")), "info_toy() field 'data_type'", fixed = TRUE)
  expect_error(info_of(list(kind = "comparison", data_type = "bigwig")), "info_toy() field 'data_type'", fixed = TRUE)
  expect_error(info_of(list(kind = "sample", data_type = "bw", scalebar = TRUE)),
               "info_toy() has an unknown field 'scalebar'", fixed = TRUE)
  expect_error(info_of(list(kind = "sample", data_type = "bw", scale_bar = "yes")), "info_toy() field 'scale_bar'",
               fixed = TRUE)
  expect_error(info_of(list(kind = "fixed", filter = "Annotations")), "info_toy() field 'filter'", fixed = TRUE)
  expect_error(info_of(list(kind = "fixed", order = "first")), "info_toy() field 'order'", fixed = TRUE)
  # A sample-free track needs no data type; the defaults are filled in.
  ok <- info_of(list(kind = "fixed"))
  expect_identical(ok[c("label", "api", "scale_bar", "row_panel", "needs_gene")],
                   list(label = "toy", api = "toy", scale_bar = FALSE, row_panel = FALSE, needs_gene = FALSE))
})

test_that("a declaration without both module functions, or two labels for one filter button, is an error", {
  two <- list(a = .pt_viz_info("bw"), b = .pt_viz_info("hic"))
  two$b$filter <- c(bw = "Signal")
  expect_error(.pt_filter_buttons(two), "filter button 'bw' has two labels: Coverage, Signal", fixed = TRUE)
  expect_error(.pt_check_filter_labels(two), "filter button 'bw' has two labels: Coverage, Signal", fixed = TRUE)
  # Every module but bed has both functions: bed's declaration is then the only orphan.
  keys <- setdiff(.pt_viz_keys(), "bed")
  local_mocked_bindings(.pt_viz_keys = function() keys)
  expect_error(.pt_viz_infos(), "info_bed() is defined but settings_bed() and render_bed() are not both", fixed = TRUE)
})

test_that("two sample-free types under one name get distinct track names", {
  items <- list(a = list(key = "a", label = "A", type = "undeclared_a", sample = NA_character_, name = "x"),
                b = list(key = "b", label = "B", type = "undeclared_b", sample = NA_character_, name = "x"))
  expect_identical(.pt_track_table(items)$track, c("x", "x_undeclared_b"))
})

test_that("a track type with its own detection and its own path field gets its path as `path`", {
  local_mocked_bindings(
    .pt_viz_info = function(key) list(key = key, path_field = "bw_path", entry_settings = character(0)),
    .pt_module_fn = function(prefix, key) if (prefix == "detect") function(sample_row, info) list())
  item <- list(key = "one", label = "x", type = "toy", sample = "A", file = NA_character_, colour_key = "A",
               name = "atac", bw_path = "a.bw", extra = 1, broken = FALSE)
  entry <- .pt_build_entries("one", list(one = item), list(A = "#000000"), list())[[1]]
  expect_identical(entry[c("sample", "colour", "name", "extra", "path")],
                   list(sample = "A", colour = "#000000", name = "atac", extra = 1, path = "a.bw"))
  expect_false("bw_path" %in% names(entry))
  # Its item without a path gives an entry without one (as Transcript expression, which reads a folder).
  item$bw_path <- NA_character_
  expect_false("path" %in% names(.pt_build_entries("one", list(one = item), list(A = "#000000"), list())[[1]]))
  # Several paths arrive as they are, and a field named like one of the entry's own is not copied.
  item <- utils::modifyList(item, list(bw_path = c("a.bw", "b.bw"), colour = "#ffffff", height = 9, viz = "other"))
  entry <- .pt_build_entries("one", list(one = item), list(A = "#000000"), list())[[1]]
  expect_identical(entry$path, c("a.bw", "b.bw"))
  expect_identical(entry[c("viz", "colour")], list(viz = "toy", colour = "#000000"))
  expect_null(entry$height)
  expect_false(anyDuplicated(names(entry)) > 0)
  # The entry's path comes from the path field only: a field named `path` is not copied, and
  # an empty path field counts as none.
  item <- utils::modifyList(item, list(bw_path = NA_character_, path = "stale.bw"))
  expect_false("path" %in% names(.pt_build_entries("one", list(one = item), list(A = "#000000"), list())[[1]]))
  item$bw_path <- character(0)
  expect_false("path" %in% names(.pt_build_entries("one", list(one = item), list(A = "#000000"), list())[[1]]))
  item$bw_path <- "a.bw"
  expect_identical(.pt_build_entries("one", list(one = item), list(A = "#000000"), list())[[1]]$path, "a.bw")
})

test_that("an entry setting comes from the settings, never from a track item's field of that name", {
  local_mocked_bindings(
    .pt_viz_info = function(key) list(key = key, path_field = "file", entry_settings = "toy_style"),
    .pt_module_fn = function(prefix, key) NULL)
  item <- list(key = "one", label = "x", type = "toy", sample = "A", file = "a.bw", colour_key = "A",
               name = "atac", toy_style = "from the item", broken = FALSE)
  build <- function(settings) .pt_build_entries("one", list(one = item), list(A = "#000000"), settings)[[1]]
  expect_identical(build(list(toy_style = "square"))$toy_style, "square")
  expect_false("toy_style" %in% names(build(list())))
})

test_that("tracks of one type share a value range per data name", {
  key <- function(name) .pt_shared_range_key("bw_range_", list(name = name))
  expect_identical(key("ATAC seq/1"), "bw_range_ATAC_seq_1")
  expect_false(identical(key("atac"), key("h3k27ac")))
})

test_that("a track type can't take the old junction key's prefix", {
  real <- .pt_viz_keys()
  keys <- c(real, "sashimi_plot")
  local_mocked_bindings(.pt_viz_keys = function() keys)
  expect_error(.pt_viz_infos(), "can't be 'sashimi' or start with 'sashimi_' (sashimi_plot)", fixed = TRUE)
  keys <- c(real, "sashimi")
  expect_error(.pt_viz_infos(), "can't be 'sashimi' or start with 'sashimi_' (sashimi)", fixed = TRUE)
  # A key that only begins with the same letters is an ordinary key.
  keys <- c(real, "sashimi2")
  said <- tryCatch({ .pt_viz_infos(); "" }, error = conditionMessage)
  expect_false(grepl("old name of the junctions track", said, fixed = TRUE), info = said)
})

test_that("a per-sample track with no file keeps an NA path, so its module can say which file is missing", {
  items <- list(one = list(key = "one", label = "x", type = "hic", sample = "A", file = NA_character_,
                           colour_key = "A", name = "hic", broken = TRUE))
  entry <- .pt_build_entries("one", items, list(A = "#000000"), list(hic_style = "rectangle"))[[1]]
  expect_true("path" %in% names(entry))
  expect_identical(entry$path, NA_character_)
})

test_that("a track item of an undeclared type keeps its data and gets a name", {
  items <- list(
    one = list(key = "one", label = "x", type = "undeclared", sample = "A", file = "a.bw", colour_key = "A",
               name = "atac", broken = FALSE),
    cmp = list(key = "cmp", label = "x", type = "undeclared_pair", sample = NA_character_, file = NA_character_,
               colour_key = NA_character_, name = "atac", samples = c("A", "B"), control_label = "A",
               treatment_label = "B", control_path = "a.bw", treatment_path = "b.bw"))
  entry <- .pt_build_entries("one", items, c(A = "#E41A1C"), list())[[1]]
  expect_identical(entry[c("viz", "sample", "colour", "name", "path")],
                   list(viz = "undeclared", sample = "A", colour = "#E41A1C", name = "atac", path = "a.bw"))
  expect_identical(.pt_track_table(items)$track, c("atac", "atac_undeclared_pair"))
})

test_that("the fixture helpers stop on a type they do not know", {
  row <- fake_rows()[[1]]
  expect_error(build_entry("genome_axis", row), "not a per-sample or sample-free track type")
  expect_error(build_entry("bw_diff_stacked", row), "not a per-sample or sample-free track type")
  expect_error(build_comparison_entry("bw", row, row), "not a comparison track type")
  expect_error(build_comparison_entry("nope", row, row), "not a comparison track type")
  expect_identical(build_comparison_entry("junctions_mirrored_diff", row, row)$control_path, "A_junctions.bed")
})

# ---- Every module draws ----------------------------------------------------------

test_that("every track type draws a one-track figure without 'track failed'", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  # The tracks the app would offer for GM12878 and for K562 against it, and their entries,
  # built as the app builds them: a new track type is covered without a fixture of its own.
  rows <- unname(test_rows("human")[c("GM12878", "K562")])
  settings <- test_settings("human")
  items <- c(build_track_items(rows), build_comparison_track_items(rows))
  colours <- setNames(lapply(rows, function(r) r$colour %||% "#1B9E77"), vapply(rows, `[[`, "", "display"))
  for (key in .pt_viz_keys()) {
    of_type <- Filter(function(it) identical(it$type, key), items)
    if (.pt_viz_info(key)$offered) expect_true(length(of_type) > 0, info = key)
    # Not offered: the genome axis, which render_all() adds to every figure.
    entries <- if (length(of_type)) .pt_build_entries(of_type[[1]]$key, items, colours, settings) else list()
    # The human Hi-C data has contacts only in a wider window (see test-modules.R).
    region <- test_region("human", flank_pct = if (key == "hic") 25 else 0.15)
    svg <- tempfile(fileext = ".svg")
    msgs <- character(0)
    withCallingHandlers(
      suppressWarnings(render_all(entries, region, settings,
                                  device_fn = function(w, h) svglite::svglite(svg, width = w, height = h))),
      message = function(m) {
        msgs <<- c(msgs, conditionMessage(m))
        invokeRestart("muffleMessage")
      })
    grDevices::dev.off()
    expect_false(any(grepl("track failed", msgs)), info = key)
    expect_false(any(grepl("track failed", readLines(svg, warn = FALSE))), info = key)
  }
})
