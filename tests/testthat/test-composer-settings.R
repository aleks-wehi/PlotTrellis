# The Figure settings panel through the app's server (.pt_main_server(),
# R/app-screen-server.R). Every non-core setting is adjustable live (s_<key> inputs),
# reaches render_all() and the downloaded R script, and is persisted only by Save as
# project default. Drives the real server headlessly on the mouse fixture, with requests
# sent as the page sends them (pattern of test-main-server.R); the block on its own is
# test-composer-settings-block.R's.

# One non-default value per non-core setting. Must list exactly the non-core keys of
# .pt_default_settings(), so a newly added setting fails here until it is wired up.
SETTING_TEST_VALUES <- list(
  font = "Helvetica", fontsize = 14, page_w = 12, chr_fmt = "space", show_track_axes = FALSE,
  hic_norm = "VC", hic_style = "triangle", hic_colour_scale = "log", zrange_cap = 50, atac_height = 1.6, hic_height = 3,
  hic_res_mode = "finest", hic_max_bins_per_view = 800, hic_resolution = 25000, gene_flank_pct = 50,
  min_junction_support = 3, min_display_tpm = 5, uniprot_reviewed_only = TRUE,
  arc_clip_partial = TRUE, diff_colour_lo = "#111111", diff_colour_hi = "#222222",
  junction_opacity_cap = 0.5, reference_colour = "#333333", comparison_colour = "#444444",
  show_title = FALSE, show_scale_bar = FALSE, junction_height = 1.4, diff_shared_alpha = 0.4
)
# Settings that never reach render_all(): the flank is applied at gene lookup (into the
# gene's view), and reviewed-only only changes the Gene info UniProt link. Checked
# separately (below, and test-main-server.R / test-main-api.R).
NOT_RENDER_SETTINGS <- c("gene_flank_pct", "uniprot_reviewed_only")

# A fresh project (json + its folder), its startup context and the app built on it.
new_composer <- function(extra = list()) {
  proj <- tempfile("pt-composer-"); dir.create(proj)
  seed_project_genome(proj)
  json <- file.path(proj, "project_settings.json")
  jsonlite::write_json(c(list(data_dir = FIXTURES$mouse$data_dir, samples = FIXTURES$mouse$samples,
                              gtf_path = GENOMES$mouse$gtf_path, orgdb = GENOMES$mouse$orgdb,
                              genome_label = GENOMES$mouse$genome_label), extra),
                       json, auto_unbox = TRUE, null = "null")
  ctx <- suppressMessages(.pt_composer_context(json, .pt_read_settings(json, .pt_default_settings())))
  list(proj = proj, json = json, ctx = ctx, app = suppressMessages(.pt_main_app(ctx, json)))
}

# A pt_rpc input value: the requests as the JSON array the page sends.
rpc <- function(...) as.character(jsonlite::toJSON(list(...), auto_unbox = TRUE, digits = NA))

# The reply to request `id` (a string; the last one with that id), its JSON read back.
reply_to <- function(rec, id) {
  r <- Filter(function(m) identical(m$type, "pt_rpc_reply") && identical(m$message$id, id), rec$custom)
  expect_gt(length(r), 0)
  jsonlite::fromJSON(r[[length(r)]]$message$json, simplifyVector = FALSE)
}

# Sends one request inside testServer() and returns its reply.
ask <- function(session, rec, id, op, args = list()) {
  session$setInputs(pt_rpc = rpc(list(id = id, op = op, args = args)))
  reply_to(rec, id)
}

# A render request for a "genes"-only region figure on the mouse fixture.
COORD_FIGURE <- list(chr = "chr2", s = 103300000, e = 103500000, picked = list(),
                     tracks = list(list(type = "genes")))

# render_all() stand-in that records its arguments and draws an empty figure, so the
# server still has an SVG file to send.
recording_render_all <- function(calls) function(entries, region_params, settings, device_fn) {
  calls$all[[length(calls$all) + 1]] <- list(entries = entries, region = region_params, settings = settings)
  device_fn(1, 1)
  invisible(NULL)
}

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

test_that("every setting input reaches render_all() and the R script", {
  skip_if_not(PT_TEST_MOUSE_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  cp <- new_composer()
  on.exit(unlink(cp$proj, recursive = TRUE), add = TRUE)
  calls <- new.env(); calls$all <- list()
  local_mocked_bindings(render_all = recording_render_all(calls))
  rec <- recording_session()
  shiny::testServer(cp$app, {
    session$flushReact()
    # Nothing set: the live settings are the saved ones.
    for (k in names(SETTING_TEST_VALUES)) expect_identical(fs$live_cfg()[[k]], fs$saved_cfg()[[k]], info = k)
    expect_false(fs$settings_dirty())
    do.call(session$setInputs, setNames(SETTING_TEST_VALUES, paste0("s_", names(SETTING_TEST_VALUES))))
    expect_true(fs$settings_dirty())
    fig <- ask(session, rec, "r1", "render", COORD_FIGURE)
    expect_null(fig$error)
    expect_length(calls$all, 1)
    seen <- calls$all[[1]]$settings
    for (k in setdiff(names(SETTING_TEST_VALUES), c(NOT_RENDER_SETTINGS, "zrange_cap"))) {
      v <- SETTING_TEST_VALUES[[k]]
      expect_equal(seen[[k]], v, info = k)
      # The downloaded script carries the same value (deparsed into its settings list).
      # hic_style is carried by each Hi-C track instead (checked with a Hi-C figure below).
      if (k != "hic_style")
        expect_true(grepl(paste0(k, " = ", deparse(if (k == "fontsize") as.integer(v) else v)),
                          fig$script, fixed = TRUE), info = k)
    }
    expect_equal(seen$zrange_cap, 50)
    # Saved values are untouched by live changes.
    expect_equal(fs$saved_cfg()$fontsize, 8L)
  }, session = rec$session)
})

test_that("a Hi-C track takes the live Hi-C style in its entry, named once in the script", {
  skip_if_not(PT_TEST_MOUSE_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  cp <- new_composer()
  on.exit(unlink(cp$proj, recursive = TRUE), add = TRUE)
  tab <- .pt_track_table(cp$ctx$TRACK_ITEMS)
  hic <- tab[vapply(cp$ctx$TRACK_ITEMS[tab$key], function(x) identical(x$type, "hic"), logical(1)) &
               tab$kind == "sample", ][1, ]
  expect_false(is.na(hic$key))
  figure <- COORD_FIGURE
  figure$tracks <- list(list(type = "genes"), list(type = hic$track, sample = hic$sample))
  calls <- new.env(); calls$all <- list()
  local_mocked_bindings(render_all = recording_render_all(calls))
  hic_style <- function(call) Filter(function(e) identical(e$viz, "hic"), call$entries)[[1]]$hic_style
  rec <- recording_session()
  shiny::testServer(cp$app, {
    session$flushReact()
    expect_null(ask(session, rec, "r1", "render", figure)$error)
    expect_equal(hic_style(calls$all[[1]]), "rectangle")
    session$setInputs(s_hic_style = "triangle", s_fontsize = 14)
    fig <- ask(session, rec, "r2", "render", figure)
    expect_null(fig$error)
    expect_equal(hic_style(calls$all[[2]]), "triangle")
    # The script has it once, on the Hi-C track, not in the settings block.
    expect_identical(lengths(regmatches(fig$script, gregexpr("hic_style = \"triangle\"", fig$script))), 1L)
    expect_equal(calls$all[[2]]$settings$fontsize, 14L)
    expect_match(fig$script, "fontsize = 14L", fixed = TRUE)
  }, session = rec$session)
})

test_that("a font size change changes the rendered figure", {
  skip_if_not(PT_TEST_MOUSE_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  cp <- new_composer()
  on.exit(unlink(cp$proj, recursive = TRUE), add = TRUE)
  rec <- recording_session()
  shiny::testServer(cp$app, {
    session$flushReact()
    fig8 <- ask(session, rec, "r1", "render", COORD_FIGURE)
    expect_null(fig8$error)
    session$setInputs(s_fontsize = 14)
    fig14 <- ask(session, rec, "r2", "render", COORD_FIGURE)
    expect_null(fig14$error)
    expect_match(fig14$script, "fontsize = 14L", fixed = TRUE)
    expect_false(grepl("font-size: 16", fig8$svg, fixed = TRUE))
    expect_true(grepl("font-size: 16", fig14$svg, fixed = TRUE))
  }, session = rec$session)
})

test_that("a flank change widens the current gene's view, around the typed ID's own gene", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  cp <- new_composer()
  on.exit(unlink(cp$proj, recursive = TRUE), add = TRUE)
  rec <- recording_session()
  shiny::testServer(cp$app, {
    session$flushReact()
    # Amy2a1 (Milton's OrgDb) by one of its 3 IDs.
    r15 <- ask(session, rec, "l1", "locus", list(q = "ENSMUSG00000093931"))
    expect_null(r15$error)
    session$setInputs(s_gene_flank_pct = 50)
    r50 <- ask(session, rec, "l2", "locus", list(q = "ENSMUSG00000093931"))
    expect_gt(r50$view$e - r50$view$s, r15$view$e - r15$view$s)
    # Still the typed ID's gene, not whichever ID the symbol resolves to first.
    expect_equal((r50$view$s + r50$view$e) / 2, (r15$view$s + r15$view$e) / 2, tolerance = 1)
    expect_identical(r50[c("chr", "s", "e")], r15[c("chr", "s", "e")])
  }, session = rec$session)
})

test_that("the app draws with the project's cache folder, and Edit project settings stops it with that word", {
  skip_if_not(PT_TEST_MOUSE_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  cp <- new_composer()
  on.exit(unlink(cp$proj, recursive = TRUE), add = TRUE)
  expect_identical(cp$ctx$RENDER_CFG$cache_dir, file.path(cp$proj, "cache"))
  # The stop is scheduled a moment later (so the page can show its message first): catch
  # what is scheduled instead of stopping anything.
  scheduled <- list()
  local_mocked_bindings(later = function(func, ...) {
    scheduled[[length(scheduled) + 1]] <<- func
    invisible(NULL)
  }, .package = "later")
  shiny::testServer(cp$app, {
    session$flushReact()  # the observers exist; the click below is then a change
    session$setInputs(confirm_edit_settings = 1)
  })
  stops <- unlist(lapply(scheduled, function(f) paste(deparse(body(f)), collapse = " ")))
  expect_identical(grep("stopApp", stops, value = TRUE), "stopApp(returnValue = \"edit_settings\")")
})

test_that("live changes never touch the project files; Save writes only the non-core keys", {
  skip_if_not(PT_TEST_MOUSE_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
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
    expect_true(fs$settings_dirty())
    expect_false(is.null(unsaved_settings_note()))
    session$setInputs(save_settings = 1)
    expect_false(fs$settings_dirty())
    expect_null(unsaved_settings_note())
    expect_equal(fs$saved_cfg()$fontsize, 14)
  })
  old  <- jsonlite::fromJSON(paste(before, collapse = "\n"), simplifyDataFrame = FALSE)
  json <- jsonlite::fromJSON(cp$json, simplifyDataFrame = FALSE)
  for (k in names(SETTING_TEST_VALUES)) expect_equal(json[[k]], SETTING_TEST_VALUES[[k]], info = k)
  # Core keys as they were; one the file lacked (suggested_genes here) is written as its default.
  for (k in .PT_CORE_SETTINGS)
    expect_identical(json[[k]], old[[k]] %||% .pt_default_settings()[[k]], info = k)
  expect_null(old$suggested_genes)
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
    expect_null(fs$saved_cfg()$zrange_cap)
  })
  json <- jsonlite::fromJSON(cp$json, simplifyDataFrame = FALSE)
  expect_true("zrange_cap" %in% names(json))
  expect_null(json$zrange_cap)
})

test_that("Revert to saved re-renders the panel from the saved settings", {
  skip_if_not(PT_TEST_MOUSE_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  cp <- new_composer()
  on.exit(unlink(cp$proj, recursive = TRUE), add = TRUE)
  before <- readLines(cp$json)
  shiny::testServer(cp$app, {
    session$setInputs(s_fontsize = 14, s_reference_colour = "#333333")
    expect_true(fs$settings_dirty())
    session$setInputs(revert_settings = 1)
    html <- output$display_settings_ui$html
    # The re-rendered inputs show the saved values (font size 8, the default role colour).
    expect_match(html, '<input id="s_fontsize" type="number"[^>]* value="8"')
    expect_match(html, "Shiny.setInputValue('s_reference_colour', \"#4D4D4D\")", fixed = TRUE)
    # What the browser then reports back:
    session$setInputs(s_fontsize = 8, s_reference_colour = "#4D4D4D")
    expect_false(fs$settings_dirty())
    expect_equal(fs$live_cfg()$fontsize, 8)
  })
  expect_identical(readLines(cp$json), before)
})
