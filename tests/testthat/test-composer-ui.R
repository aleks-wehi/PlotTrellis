# The composer page's builders (R/composer-ui.R) on a hand-made startup context
# (helper-composer-fakes.R): plain HTML, no genome, no app.

# The <option>s of one select in `html`, as c(value = selected?).
select_options <- function(html, id) {
  block <- regmatches(html, regexpr(sprintf('<select[^>]*id="%s"[^>]*>.*?</select>', id), html))
  opts <- regmatches(block, gregexpr('<option value="[^"]*"( selected)?>', block))[[1]]
  setNames(grepl("selected", opts, fixed = TRUE), sub('<option value="([^"]*)".*', "\\1", opts))
}

# The radio choice of `id` that is checked, or character(0).
checked_radio <- function(html, id) {
  hits <- regmatches(html, gregexpr(sprintf('name="%s" value="[^"]+" checked="checked"', id), html))[[1]]
  sub('.*value="([^"]+)".*', "\\1", hits)
}

rank_ids <- function(html) sub('data-rank-id="([^"]+)"', "\\1",
                               regmatches(html, gregexpr('data-rank-id="[^"]+"', html))[[1]])

test_that(".pt_track_js_maps maps each track key to its type and its broken flag", {
  ctx <- fake_ctx()
  m <- .pt_track_js_maps(ctx$TRACK_ITEMS, ctx$GENE_REQUIRED_TYPES)
  expect_named(m, c("TRACK_TYPE_JS", "TRACK_BROKEN_JS", "TRACK_FILTER_JS", "GENE_REQUIRED_TYPES_JS"))
  types <- jsonlite::fromJSON(m$TRACK_TYPE_JS)
  expect_identical(names(types), names(ctx$TRACK_ITEMS))
  expect_identical(unlist(types), vapply(ctx$TRACK_ITEMS, `[[`, "", "type"))
  broken <- unlist(jsonlite::fromJSON(m$TRACK_BROKEN_JS))
  # Only B's ATAC track is broken; an item with no `broken` field counts as not broken.
  expect_identical(names(broken)[broken], "bw_atac_B")
  expect_identical(broken[["genes"]], FALSE)
  expect_match(m$TRACK_BROKEN_JS, '"bw_atac_B":true', fixed = TRUE)
})

test_that(".pt_track_js_maps gives the module filters and the gene-required types as JSON", {
  m <- .pt_track_js_maps(fake_track_items(), "tpm_usage")
  expect_identical(as.character(m$TRACK_FILTER_JS), as.character(jsonlite::toJSON(.pt_type_filters(), auto_unbox = TRUE)))
  expect_identical(jsonlite::fromJSON(m$TRACK_FILTER_JS)$bw_diff_stacked, "comparison")
  # One gene-required type stays an array (not unboxed).
  expect_identical(as.character(m$GENE_REQUIRED_TYPES_JS), '["tpm_usage"]')
  expect_identical(as.character(.pt_track_js_maps(fake_track_items(), character(0))$GENE_REQUIRED_TYPES_JS), "[]")
})

test_that(".pt_track_js_maps of no items gives empty JSON arrays", {
  m <- .pt_track_js_maps(list(), character(0))
  # An empty (unnamed) list is written as [], not {}; the page only looks keys up in it.
  expect_identical(as.character(m$TRACK_TYPE_JS), "[]")
  expect_identical(as.character(m$TRACK_BROKEN_JS), "[]")
})

test_that("the page head carries the theme, the composer's stylesheet and script, and the JS maps", {
  ctx <- fake_ctx()
  # tags$head() contents render into the page head.
  html <- htmltools::renderTags(.pt_composer_head(ctx))$head
  m <- .pt_track_js_maps(ctx$TRACK_ITEMS, ctx$GENE_REQUIRED_TYPES)
  expect_match(html, paste0("var TRACK_TYPE_MAP = ", m$TRACK_TYPE_JS, ";"), fixed = TRUE)
  expect_match(html, paste0("var TRACK_BROKEN_MAP = ", m$TRACK_BROKEN_JS, ";"), fixed = TRUE)
  expect_match(html, paste0("var TRACK_FILTER_MAP = ", m$TRACK_FILTER_JS, ";"), fixed = TRUE)
  expect_match(html, paste0("var GENE_REQUIRED_TYPES = ", m$GENE_REQUIRED_TYPES_JS, ";"), fixed = TRUE)
  css <- readLines(system.file("www", "composer.css", package = "PlotTrellis"), warn = FALSE)
  js  <- readLines(system.file("www", "composer.js", package = "PlotTrellis"), warn = FALSE)
  expect_match(html, trimws(css[nzchar(trimws(css))][1]), fixed = TRUE)
  expect_match(html, "function applyTrackFilter(blockEl)", fixed = TRUE)
  expect_true(length(js) > 0)
  # The maps are defined before the page script that reads them.
  expect_lt(regexpr("var TRACK_TYPE_MAP", html, fixed = TRUE), regexpr("function applyTrackFilter", html, fixed = TRUE))
})

test_that("the Compose section starts blank with no start view", {
  html <- as.character(.pt_compose_panel_ui(fake_ctx()))
  expect_identical(checked_radio(html, "view_mode"), character(0))
  expect_identical(checked_radio(html, "region_mode"), character(0))
  for (id in c("sel_condition", "sel_reference", "sel_comparison"))
    expect_identical(select_options(html, id), setNames(c(TRUE, FALSE, FALSE), c("", "A", "B")), info = id)
  expect_match(html, 'data-display-if="input.view_mode === &#39;single&#39;"', fixed = TRUE)
  expect_match(html, 'data-display-if="input.view_mode === &#39;two&#39;"', fixed = TRUE)
})

test_that("a Two-sample start view preselects the mode, Reference, Comparison and Gene name", {
  html <- as.character(.pt_compose_panel_ui(fake_ctx(list(view_mode = "two", conditions = c("B", "A"), gene = "Pax5"))))
  expect_identical(checked_radio(html, "view_mode"), "two")
  expect_identical(checked_radio(html, "region_mode"), "gene")
  # Every condition is offered, in project order, with the start view's one selected.
  expect_identical(select_options(html, "sel_reference"), setNames(c(FALSE, FALSE, TRUE), c("", "A", "B")))
  expect_identical(names(which(select_options(html, "sel_comparison"))), "A")
  expect_identical(names(which(select_options(html, "sel_condition"))), "")
})

test_that("a One-sample start view preselects its condition only", {
  html <- as.character(.pt_compose_panel_ui(fake_ctx(list(view_mode = "single", conditions = "B"))))
  expect_identical(checked_radio(html, "view_mode"), "single")
  expect_identical(names(which(select_options(html, "sel_condition"))), "B")
  expect_identical(names(which(select_options(html, "sel_reference"))), "")
  expect_identical(names(which(select_options(html, "sel_comparison"))), "")
  # No gene in the start view: no region mode chosen.
  expect_identical(checked_radio(html, "region_mode"), character(0))
})

test_that("the track picker starts with Multi's pool and Gene annotation, whatever the start view", {
  ctx <- fake_ctx(list(view_mode = "compose", tracks = c("bw_atac_A", "genes")))
  html <- as.character(.pt_compose_panel_ui(ctx))
  # The server replaces both lists at startup; the markup holds Multi's defaults.
  expect_identical(rank_ids(html), c(ctx$COMPOSE_AVAILABLE_KEYS, "genes"))
  expect_match(html, "Available tracks (double-click)", fixed = TRUE)
  expect_match(html, "Figure tracks (drag to reorder)", fixed = TRUE)
  expect_match(html, '<div id="tracks-block" data-available-id="available_tracks">', fixed = TRUE)
})

test_that("filter buttons cover only the track types the project has, All first", {
  html <- as.character(.pt_compose_panel_ui(fake_ctx()))
  filters <- sub('data-filter="([^"]+)"', "\\1", regmatches(html, gregexpr('data-filter="[^"]+"', html))[[1]])
  expect_identical(filters[1], "all")
  expect_setequal(filters, c("all", "comparison", "annotation", "bw", "hic", "tpm_usage"))
  expect_false(any(c("bed", "junctions") %in% filters))
  # Comparisons only in Two-sample mode.
  expect_match(html, "data-display-if=\"input.view_mode === &#39;two&#39;\"[^>]*>\\s*<button[^>]*data-filter=\"comparison\"")
  items <- fake_track_items()[c("genes", "bw_atac_A")]
  html <- as.character(.pt_compose_panel_ui(fake_ctx(items = items)))
  filters <- sub('data-filter="([^"]+)"', "\\1", regmatches(html, gregexpr('data-filter="[^"]+"', html))[[1]])
  expect_identical(filters, c("all", "annotation", "bw"))
})

test_that("the Compose section holds the region, settings and Generate outputs", {
  html <- as.character(.pt_compose_panel_ui(fake_ctx()))
  for (id in c("gene_result_ui", "chr_select_ui", "range_slider_ui", "start_num_ui", "end_num_ui", "region_error_ui",
               "tracks_hint", "display_settings_ui", "apply_settings_ui", "settings_dirty_ui", "generate_ui",
               "generate_warning_ui", "save_settings", "revert_settings", "gene_name", "compose-body", "figure_settings-body"))
    expect_match(html, sprintf('id="%s"', id), fixed = TRUE, info = id)
  expect_match(html, '<div id="compose-body" class="collapse in">', fixed = TRUE)
  expect_match(html, '<div id="figure_settings-body" class="collapse ">', fixed = TRUE)
  expect_match(html, 'id="save_settings" type="button"', fixed = TRUE)
  expect_match(html, "btn-default pt-once", fixed = TRUE)
  expect_match(html, "&quot;placeholder&quot;:&quot;e.g. Pax5&quot;|\"placeholder\":\"e.g. Pax5\"")
  expect_match(html, '"create":true', fixed = TRUE)
})

test_that("the Output section is hidden until the first render and holds every output", {
  html <- as.character(.pt_output_panel_ui())
  expect_match(html, '^<div id="output-section" class="panel panel-default plottrellis-section " style="display:none;">')
  for (id in c("render_error_ui", "zoom_ui", "svg_display_ui", "dl_svg", "dl_rcode", "rcode_output", "rcode-body"))
    expect_match(html, sprintf('id="%s"', id), fixed = TRUE, info = id)
  expect_match(html, '<div id="rcode-body" class="collapse ">', fixed = TRUE)
  expect_match(html, "Download SVG", fixed = TRUE)
  expect_match(html, "Download R script", fixed = TRUE)
  expect_match(html, "An R script that reproduces this figure.", fixed = TRUE)
})

test_that("the composer page puts the head, the project actions, Compose and Output in one tab", {
  r <- htmltools::renderTags(.pt_composer_ui(fake_ctx()))
  expect_match(r$head, "var TRACK_TYPE_MAP = ", fixed = TRUE)
  expect_match(r$html, '<small class="plottrellis-version">v', fixed = TRUE)
  expect_match(r$html, paste0("v", utils::packageVersion("PlotTrellis")), fixed = TRUE)
  expect_match(r$html, 'id="main_tabs"', fixed = TRUE)
  expect_match(r$head, "<title>PlotTrellis</title>", fixed = TRUE)
  pos <- vapply(c('id="topbar-actions"', 'id="compose-body"', 'id="output-section"'),
                function(p) regexpr(p, r$html, fixed = TRUE)[[1]], 1L)
  expect_true(all(pos > 0))
  expect_identical(unname(order(pos)), 1:3)
  # One page tab only, with an empty title (the other tabs are the Gene info panel's own).
  expect_length(regmatches(r$html, gregexpr('data-toggle="tab" data-bs-toggle="tab" data-value=""', r$html, fixed = TRUE))[[1]], 1)
  expect_match(r$html, 'data-value="Gene info"', fixed = TRUE)
})
