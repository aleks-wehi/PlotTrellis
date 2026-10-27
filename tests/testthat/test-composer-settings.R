# Phase I.3 step 28: the composer's Figure settings panel. Every non-core setting is
# adjustable live (s_<key> inputs), reaches render_all() and the downloaded R script,
# and is persisted only by Save as project default. Drives the real composer server
# headlessly (pattern of test-gene-search-app.R).

# One non-default value per non-core setting. Must list exactly the non-core keys of
# .pt_default_settings(), so a newly added setting fails here until it is wired up.
SETTING_TEST_VALUES <- list(
  font = "Helvetica", fontsize = 14, page_w = 12, chr_fmt = "space", show_track_axes = FALSE,
  hic_norm = "VC", hic_style = "triangle", zrange_cap = 50, atac_height = 1.6, hic_height = 3,
  hic_res_mode = "finest", hic_max_bins_per_view = 800, gene_flank_pct = 50,
  min_junction_support = 3, min_display_tpm = 5, uniprot_reviewed_only = TRUE,
  arc_clip_partial = TRUE, diff_colour_lo = "#111111", diff_colour_hi = "#222222",
  junction_opacity_cap = 0.5, reference_colour = "#333333", comparison_colour = "#444444",
  show_title = FALSE, show_scale_bar = FALSE, junction_height = 1.4, diff_shared_alpha = 0.4
)
# Settings that never reach render_all(): the flank is applied at gene lookup (into the
# region), and reviewed-only only changes the Gene info link-out. Checked separately.
NOT_RENDER_SETTINGS <- c("gene_flank_pct", "uniprot_reviewed_only")

# A fresh project (json + its folder) and the composer built for it.
new_composer <- function(extra = list()) {
  proj <- tempfile("pt-composer-"); dir.create(proj)
  json <- file.path(proj, "project_settings.json")
  jsonlite::write_json(c(list(data_dir = FIXTURES$mouse$data_dir, samples = FIXTURES$mouse$samples,
                              gtf_path = GENOMES$mouse$gtf_path, orgdb = GENOMES$mouse$orgdb,
                              genome_label = GENOMES$mouse$genome_label), extra),
                       json, auto_unbox = TRUE, null = "null")
  list(proj = proj, json = json,
       app = suppressMessages(.pt_composer_app(json, .pt_read_settings(json, .pt_default_settings()))))
}

# Inputs for a "genes"-only coordinates figure on the mouse fixture.
COORD_FIGURE <- list(view_mode = "compose", region_mode = "coords", region_chr = "chr2",
                     region_range = c(103300000L, 103500000L), figure_tracks = "genes")

test_that("the non-core settings list matches the panel, the collector and the test values", {
  core <- .PT_CORE_SETTINGS
  expect_true(all(core %in% names(.pt_default_settings())))
  non_core <- setdiff(names(.pt_default_settings()), core)
  expect_setequal(names(SETTING_TEST_VALUES), non_core)
  expect_setequal(names(.pt_collect_display_settings(list(), "s_", .pt_default_settings())), non_core)
  # Every non-core key has an s_<key> input in the panel, and no core key does.
  ui <- .pt_display_settings_ui(.pt_default_settings(), "s_")
  html <- as.character(tagList(ui$display, ui$diff))
  for (k in non_core) expect_true(grepl(sprintf('id="s_%s"', k), html, fixed = TRUE), info = k)
  for (k in core) expect_false(grepl(sprintf('id="s_%s"', k), html, fixed = TRUE), info = k)
})

test_that("every composer setting input reaches render_all() and the R script", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  cp <- new_composer()
  on.exit(unlink(cp$proj, recursive = TRUE), add = TRUE)
  seen <- NULL
  local_mocked_bindings(render_all = function(entries, region_params, settings, device_fn) {
    seen <<- settings
    invisible(NULL)
  })
  shiny::testServer(cp$app, {
    # Nothing set: the live settings are the saved ones.
    for (k in names(SETTING_TEST_VALUES)) expect_identical(live_cfg()[[k]], saved_cfg()[[k]], info = k)
    expect_false(settings_dirty())
    do.call(session$setInputs, COORD_FIGURE)
    do.call(session$setInputs, setNames(SETTING_TEST_VALUES, paste0("s_", names(SETTING_TEST_VALUES))))
    expect_true(settings_dirty())
    session$setInputs(generate = 1)
    expect_null(state$plot_error)
    for (k in setdiff(names(SETTING_TEST_VALUES), c(NOT_RENDER_SETTINGS, "zrange_cap"))) {
      v <- SETTING_TEST_VALUES[[k]]
      expect_equal(seen[[k]], v, info = k)
      # The downloaded script carries the same value (deparsed into its settings list).
      expect_true(grepl(paste0(k, " = ", deparse(if (k == "fontsize") as.integer(v) else v)),
                        state$rcode, fixed = TRUE), info = k)
    }
    expect_equal(seen$zrange, c(0, 50))
    # Saved values are untouched by live changes.
    expect_equal(saved_cfg()$fontsize, 8L)
  })
})

test_that("a composer font size change changes the rendered figure", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  cp <- new_composer()
  on.exit(unlink(cp$proj, recursive = TRUE), add = TRUE)
  shiny::testServer(cp$app, {
    do.call(session$setInputs, COORD_FIGURE)
    session$setInputs(generate = 1)
    svg8 <- paste(readLines(state$svg_path, warn = FALSE), collapse = "\n")
    session$setInputs(s_fontsize = 14, generate = 2)
    expect_null(state$plot_error)
    svg14 <- paste(readLines(state$svg_path, warn = FALSE), collapse = "\n")
    expect_match(state$rcode, "fontsize = 14L", fixed = TRUE)
    # The auto title is drawn at fontsize + 2.
    expect_false(grepl("font-size: 16", svg8, fixed = TRUE))
    expect_true(grepl("font-size: 16", svg14, fixed = TRUE))
  })
})

test_that("a composer flank change re-resolves the current gene's region", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  cp <- new_composer()
  on.exit(unlink(cp$proj, recursive = TRUE), add = TRUE)
  shiny::testServer(cp$app, {
    session$setInputs(gene_name = "ENSMUSG00000093931")  # Amy2a1 by one of its 3 IDs
    session$setInputs(gene_name = "Amy2a1")              # the echo
    r15 <- state$gene_region
    session$setInputs(s_gene_flank_pct = 50)
    r50 <- state$gene_region
    expect_gt(r50$end - r50$start, r15$end - r15$start)
    # Still the typed ID's gene, not whichever ID the symbol resolves to first.
    expect_equal((r50$start + r50$end) / 2, (r15$start + r15$end) / 2, tolerance = 1)
  })
})

test_that("zoom and Apply to figure rebuild the figure with the live settings", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  cp <- new_composer()
  on.exit(unlink(cp$proj, recursive = TRUE), add = TRUE)
  items <- get("TRACK_ITEMS", envir = environment(cp$app$serverFuncSource()))
  hic_key <- names(items)[vapply(items, function(x) identical(x$type, "hic"), logical(1))][1]
  expect_false(is.na(hic_key))
  calls <- list()
  local_mocked_bindings(render_all = function(entries, region_params, settings, device_fn) {
    calls[[length(calls) + 1]] <<- list(entries = entries, region = region_params, settings = settings)
    invisible(NULL)
  })
  hic_style <- function(call) Filter(function(e) identical(e$viz, "hic"), call$entries)[[1]]$hic_style
  shiny::testServer(cp$app, {
    fig <- utils::modifyList(COORD_FIGURE, list(figure_tracks = c("genes", hic_key)))
    do.call(session$setInputs, fig)
    session$setInputs(plot_title = "My title", generate = 1)
    first <- calls[[length(calls)]]
    expect_equal(hic_style(first), "rectangle")
    session$setInputs(s_hic_style = "triangle", s_fontsize = 14, zoom_in_15x = 1)
    zoomed <- calls[[length(calls)]]
    expect_equal(hic_style(zoomed), "triangle")  # baked into the entry, rebuilt on zoom
    expect_equal(zoomed$settings$fontsize, 14L)
    expect_equal(zoomed$settings$title, "My title")
    expect_lt(zoomed$region$end - zoomed$region$start, first$region$end - first$region$start)
    expect_match(state$rcode, "fontsize = 14L", fixed = TRUE)
    session$setInputs(s_fontsize = 10, apply_settings = 1)
    applied <- calls[[length(calls)]]
    expect_equal(applied$settings$fontsize, 10L)
    expect_equal(applied$settings$title, "My title")
    expect_equal(applied$region[c("chrom", "start", "end")], zoomed$region[c("chrom", "start", "end")])
  })
})

test_that("zoom keeps a custom title in the rendered figure", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  cp <- new_composer()
  on.exit(unlink(cp$proj, recursive = TRUE), add = TRUE)
  shiny::testServer(cp$app, {
    do.call(session$setInputs, COORD_FIGURE)
    session$setInputs(plot_title = "My title", generate = 1)
    expect_true(any(grepl("My title", readLines(state$svg_path, warn = FALSE), fixed = TRUE)))
    session$setInputs(zoom_in_15x = 1)
    expect_null(state$plot_error)
    expect_true(any(grepl("My title", readLines(state$svg_path, warn = FALSE), fixed = TRUE)))
  })
})

test_that("live changes never touch the project files; Save writes only the non-core keys", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  cp <- new_composer(list(zrange_cap = NULL, page_w = 10))
  on.exit(unlink(cp$proj, recursive = TRUE), add = TRUE)
  mirror <- file.path(cp$proj, "project_samples.txt")
  before <- readLines(cp$json)
  # An unsaved change, then the session ends: nothing is written.
  shiny::testServer(cp$app, session$setInputs(s_fontsize = 14))
  expect_identical(readLines(cp$json), before)
  expect_false(file.exists(mirror))

  shiny::testServer(cp$app, {
    do.call(session$setInputs, setNames(SETTING_TEST_VALUES, paste0("s_", names(SETTING_TEST_VALUES))))
    expect_true(settings_dirty())
    expect_false(is.null(unsaved_settings_note()))
    session$setInputs(save_settings = 1)
    expect_false(settings_dirty())
    expect_null(unsaved_settings_note())
    expect_equal(saved_cfg()$fontsize, 14)
  })
  old  <- jsonlite::fromJSON(paste(before, collapse = "\n"), simplifyDataFrame = FALSE)
  json <- jsonlite::fromJSON(cp$json, simplifyDataFrame = FALSE)
  for (k in names(SETTING_TEST_VALUES)) expect_equal(json[[k]], SETTING_TEST_VALUES[[k]], info = k)
  for (k in .PT_CORE_SETTINGS) expect_identical(json[[k]], old[[k]], info = k)
  # The saved file reads back as the same settings, and the samples.txt copy carries them.
  cfg <- .pt_read_settings(cp$json, .pt_default_settings())
  for (k in names(SETTING_TEST_VALUES)) expect_equal(cfg[[k]], SETTING_TEST_VALUES[[k]], info = k)
  expect_true(file.exists(mirror))
  expect_true("# fontsize: 14" %in% readLines(mirror))
  res <- import_samples_txt(mirror, defaults = .pt_default_settings())
  expect_length(res$errors, 0)
  expect_equal(length(res$samples), length(old$samples))

  # Saving "auto" back writes null, not a dropped key or {}.
  shiny::testServer(cp$app, {
    session$setInputs(s_zrange_cap = NA, save_settings = 1)
    expect_null(saved_cfg()$zrange_cap)
  })
  json <- jsonlite::fromJSON(cp$json, simplifyDataFrame = FALSE)
  expect_true("zrange_cap" %in% names(json))
  expect_null(json$zrange_cap)
})

test_that("Revert to saved re-renders the panel from the saved settings", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  cp <- new_composer()
  on.exit(unlink(cp$proj, recursive = TRUE), add = TRUE)
  before <- readLines(cp$json)
  shiny::testServer(cp$app, {
    session$setInputs(s_fontsize = 14, s_reference_colour = "#333333")
    expect_true(settings_dirty())
    session$setInputs(revert_settings = 1)
    html <- output$display_settings_ui$html
    # The re-rendered inputs show the saved values (font size 8, the default role colour).
    expect_match(html, '<input id="s_fontsize" type="number"[^>]* value="8"')
    expect_match(html, "Shiny.setInputValue('s_reference_colour', \"#4D4D4D\")", fixed = TRUE)
    # What the browser then reports back:
    session$setInputs(s_fontsize = 8, s_reference_colour = "#4D4D4D")
    expect_false(settings_dirty())
    expect_equal(live_cfg()$fontsize, 8)
  })
  expect_identical(readLines(cp$json), before)
})
