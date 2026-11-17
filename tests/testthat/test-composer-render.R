# The composer's rendering (R/composer-render.R) on its own: .pt_render_svg() and
# assemble_download_script() are replaced by recorders, the track items are hand-made
# (helper-composer-fakes.R), and the region, mode and live settings are reactive values the
# test sets. No genome, no figure is drawn.

RENDER_TEST_CFG <- c(.pt_default_settings(),
                     list(txdb_path = "/no/such/txdb.sqlite", cache_dir = "/no/such/cache", orgdb_pkg = "org.Fake.eg.db",
                          tx_name_lookup_path = "/no/such/tx_names.rds", tx_meta_lookup_path = "/no/such/tx_meta.rds"))

render_app <- function(worker = NULL, ctx = fake_ctx(), rec = NULL) function(input, output, session) {
  rec <- rec  # a local, so the test code sees it
  state <- reactiveValues(svg_path = NULL, rcode = NULL, plot_error = NULL, last_region = NULL,
                          last_labels = NULL, last_title = "", rendering = FALSE)
  cfg  <- reactiveVal(RENDER_TEST_CFG)
  rr   <- reactiveVal(NULL)   # resolved_region()
  mode <- reactiveVal(TRUE)   # has_mode()
  .pt_render_server(input, output, session, state, ctx, worker, live_cfg = reactive(cfg()),
                    resolved_region = reactive(rr()), has_mode = reactive(mode()))
}

# A stand-in for the session's background render process: records each render and keeps its
# `done` callback for the test to call; counts cancels.
fake_worker <- function() {
  w <- new.env()
  w$renders <- list(); w$cancels <- 0
  w$render <- function(entries, region_params, settings, done)
    w$renders[[length(w$renders) + 1]] <- list(entries = entries, region_params = region_params,
                                               settings = settings, done = done)
  w$cancel <- function() w$cancels <- w$cancels + 1
  w
}

# Runs `code` against the render server (it sees session, output, state, cfg, rr, mode and
# `rec`). `rec` holds the in-process renders, the script assemblies, the notifications and the
# messages to the page; `rec$render_result` / `rec$script_error` make the next draw or script fail.
drive_render <- function(code, worker = NULL, ctx = fake_ctx()) {
  code <- rlang::enquo(code)
  rec <- recording_session()
  rec$renders <- list(); rec$scripts <- list(); rec$notes <- list()
  rec$render_result <- NULL; rec$script_error <- NULL
  local_mocked_bindings(
    .pt_render_svg = function(entries, region_params, settings) {
      rec$renders[[length(rec$renders) + 1]] <- list(entries = entries, region_params = region_params, settings = settings)
      writeLines("<svg width='1'></svg>", settings$output_path)
      rec$render_result
    },
    assemble_download_script = function(entries, region_params, settings) {
      rec$scripts[[length(rec$scripts) + 1]] <- list(entries = entries, region_params = region_params, settings = settings)
      if (!is.null(rec$script_error)) stop(rec$script_error)
      paste("script", length(rec$scripts))
    })
  local_mocked_bindings(showNotification = function(ui, action = NULL, duration = 5, closeButton = TRUE, id = NULL,
                                                    type = "default", session = NULL)
    rec$notes[[length(rec$notes) + 1]] <- list(ui = ui, type = type))
  rlang::inject(shiny::testServer(render_app(worker, ctx, rec), !!code, session = rec$session))
  rec
}

REGION <- list(rgn = list(chrom = "chr1", start = 1000000, end = 2000000), gene_name = "Pax5")

# A Generate of Gene annotation and A's ATAC track over REGION, titled "My title".
generate <- quote({
  rr(REGION)
  session$setInputs(view_mode = "single", figure_tracks = c("genes", "bw_atac_A"), plot_title = "  My title  ")
  session$setInputs(generate = 1)
})

custom_types <- function(rec) vapply(rec$custom, `[[`, "", "type")

# The region each in-process render drew, as c(start, end).
drawn <- function(rec) lapply(rec$renders, function(r) c(r$region_params$start, r$region_params$end))

test_that("Generate draws the picked tracks over the resolved region with the live settings", {
  rec <- drive_render(!!generate)
  expect_length(rec$renders, 1)
  r <- rec$renders[[1]]
  expect_identical(r$region_params, list(chrom = "chr1", start = 1000000L, end = 2000000L, gene = "Pax5"))
  ctx <- fake_ctx()
  expect_identical(r$entries, .pt_build_entries(c("genes", "bw_atac_A"), ctx$TRACK_ITEMS, ctx$DEFAULT_SAMPLE_COLS, RENDER_TEST_CFG))
  expect_identical(vapply(r$entries, `[[`, "", "viz"), c("genes", "bw"))
  # The settings: .pt_settings_base() of the live settings, the trimmed title, a new .svg file.
  expect_identical(r$settings[names(.pt_settings_base(RENDER_TEST_CFG))], .pt_settings_base(RENDER_TEST_CFG))
  expect_identical(r$settings$title, "My title")
  expect_match(r$settings$output_path, "\\.svg$")
  # The script is assembled from exactly what was drawn.
  expect_identical(rec$scripts[[1]], r)
})

test_that("a finished render puts the figure, script, region, tracks, title and gene in state", {
  rec <- drive_render({
    !!generate
    out <- rec$renders[[1]]$settings$output_path
    expect_identical(state$svg_path, out)
    expect_true(file.exists(out))
    expect_identical(state$rcode, "script 1")
    expect_identical(state$last_region, list(chrom = "chr1", start = 1000000L, end = 2000000L, gene = "Pax5"))
    expect_identical(state$last_labels, c("genes", "bw_atac_A"))
    expect_identical(state$last_title, "My title")
    expect_identical(state$rendered_gene, "Pax5")
    expect_null(state$plot_error)
    expect_false(state$rendering)
  })
  expect_identical(custom_types(rec), c("show_output", "hide_spinner"))
})

test_that("Generate does nothing without a mode, tracks or a region, or while a render runs", {
  rec <- drive_render({
    rr(REGION)
    session$setInputs(view_mode = "single", figure_tracks = character(0), generate = 1)  # no tracks
    session$setInputs(figure_tracks = "genes")
    rr(NULL); session$setInputs(generate = 2)                                           # no region
    rr(REGION); mode(FALSE); session$setInputs(generate = 3)                            # no mode
    mode(TRUE); state$rendering <- TRUE; session$setInputs(generate = 4)                # busy
  })
  expect_length(rec$renders, 0)
  expect_length(rec$custom, 0)
})

test_that("Two samples needs both conditions before Generate draws", {
  rec <- drive_render({
    rr(REGION)
    session$setInputs(view_mode = "two", sel_reference = "A", sel_comparison = "", figure_tracks = "genes")
    session$setInputs(generate = 1)
    expect_length(rec$renders, 0)
    session$setInputs(sel_comparison = "B", generate = 2)
    expect_length(rec$renders, 1)
  })
})

test_that("a blank title is an empty string, and a coordinates figure has no gene", {
  rec <- drive_render({
    rr(list(rgn = list(chrom = "chr2", start = 10, end = 5000), gene_name = NA_character_))
    session$setInputs(view_mode = "compose", figure_tracks = "genes", generate = 1)
    expect_identical(state$last_title, "")
    expect_identical(state$rendered_gene, NA_character_)
  })
  expect_identical(rec$renders[[1]]$settings$title, "")
  expect_identical(rec$renders[[1]]$region_params, list(chrom = "chr2", start = 10L, end = 5000L, gene = NA_character_))
})

test_that("a failed draw keeps the previous figure, shows the error and deletes the new file", {
  rec <- drive_render({
    !!generate
    first <- state$svg_path
    rec$render_result <- "Couldn't draw: boom"
    session$setInputs(generate = 2)
    expect_identical(state$plot_error, "Couldn't draw: boom")
    expect_identical(state$svg_path, first)
    expect_true(file.exists(first))
    expect_false(file.exists(rec$renders[[2]]$settings$output_path))
    expect_identical(state$rcode, "script 1")
    expect_false(state$rendering)
  })
  # The failed draw's script is never assembled; the page is still told the render ended.
  expect_length(rec$scripts, 1)
  expect_identical(tail(custom_types(rec), 2), c("show_output", "hide_spinner"))
})

test_that("a script that cannot be assembled is reported as the render's error", {
  rec <- drive_render({
    rec$script_error <- "no script for you"
    !!generate
    expect_identical(state$plot_error, "no script for you")
    expect_null(state$svg_path)
    expect_false(file.exists(rec$renders[[1]]$settings$output_path))
  })
})

test_that("a new figure deletes the previous figure's file and clears its error", {
  rec <- drive_render({
    rec$render_result <- "first failed"
    !!generate
    expect_identical(state$plot_error, "first failed")
    rec$render_result <- NULL
    session$setInputs(generate = 2)
    first <- state$svg_path
    session$setInputs(generate = 3)
    expect_false(file.exists(first))
    expect_true(file.exists(state$svg_path))
    expect_null(state$plot_error)
  })
})

test_that("each zoom button recentres the last region by its factor", {
  expected <- list(zoom_in_15x  = c(1166667L, 1833333L), zoom_in_3x   = c(1333333L, 1666667L),
                   zoom_in_10x  = c(1450000L, 1550000L), zoom_out_15x = c(750000L, 2250000L),
                   zoom_out_3x  = c(1L, 3000000L),       zoom_out_10x = c(1L, 6500000L))
  for (btn in names(expected)) {
    rec <- drive_render({
      !!generate
      do.call(session$setInputs, setNames(list(1), btn))
    })
    expect_identical(drawn(rec)[[2]], expected[[btn]], info = btn)
  }
})

test_that("a zoom keeps the last figure's tracks, gene and title, and moves the region inputs", {
  rec <- drive_render({
    !!generate
    session$setInputs(figure_tracks = "genes", plot_title = "Other")
    rr(list(rgn = list(chrom = "chr2", start = 1, end = 10), gene_name = "Elf5"))
    session$setInputs(zoom_in_3x = 1)
    expect_identical(state$last_region, list(chrom = "chr1", start = 1333333L, end = 1666667L, gene = "Pax5"))
    expect_identical(state$last_labels, c("genes", "bw_atac_A"))
    expect_identical(state$last_title, "My title")
  })
  expect_identical(vapply(rec$renders[[2]]$entries, `[[`, "", "viz"), c("genes", "bw"))
  expect_identical(rec$renders[[2]]$settings$title, "My title")
  sent <- setNames(lapply(rec$input, function(m) m$message$value), vapply(rec$input, `[[`, "", "id"))
  expect_identical(as.numeric(trimws(sent$region_range)), c(1333333, 1666667))
  expect_identical(sent$region_start_text, "1.333Mb")
  expect_identical(sent$region_end_text, "1.667Mb")
})

test_that("zoom stops at the chromosome's ends and never shows less than 1,000 bp", {
  ctx <- fake_ctx()
  rec <- drive_render({
    rr(list(rgn = list(chrom = "chr2", start = 4990000, end = 5000000), gene_name = NA_character_))
    session$setInputs(view_mode = "compose", figure_tracks = "genes", generate = 1)
    session$setInputs(zoom_out_10x = 1)   # chr2 is 5Mb long
    rr(list(rgn = list(chrom = "chr1", start = 1000, end = 2000), gene_name = NA_character_))
    session$setInputs(generate = 2)
    session$setInputs(zoom_in_10x = 1)    # half width at least 500
    session$setInputs(zoom_in_3x = 1)
  })
  expect_identical(drawn(rec), list(c(4990000L, 5000000L), c(4945000L, 5000000L),
                                    c(1000L, 2000L), c(1000L, 2000L), c(1000L, 2000L)))
})

test_that("zoom does nothing before the first figure or while a render runs", {
  rec <- drive_render({
    session$setInputs(zoom_in_3x = 1, zoom_out_3x = 1, apply_settings = 1)
    !!generate
    state$rendering <- TRUE
    session$setInputs(zoom_in_10x = 1)
  })
  expect_length(rec$renders, 1)
})

test_that("zoom on a chromosome missing from the chromosome lengths uses the 200 Mb default", {
  # Bug 110 (fixed on main 2026-11-17): a contig the chromosome sizes do not list no longer
  # stops the zoom observer; the new region is clipped to the 2e8 default.
  rec <- drive_render({
    rr(list(rgn = list(chrom = "chrUn_x", start = 1000, end = 5000), gene_name = "Gm1"))
    session$setInputs(view_mode = "compose", figure_tracks = "genes", generate = 1)
    session$setInputs(zoom_in_3x = 1)
    expect_identical(state$last_region[c("chrom", "start", "end")], list(chrom = "chrUn_x", start = 2333L, end = 3667L))
  })
  expect_length(rec$renders, 2)
})

test_that("Apply to figure appears after the first figure and redraws it with the live settings", {
  rec <- drive_render({
    expect_error(output$apply_settings_ui, class = "shiny.silent.error")
    !!generate
    html <- ui_html(output$apply_settings_ui)
    expect_match(html, 'id="apply_settings"', fixed = TRUE)
    expect_match(html, "btn-primary generate-btn", fixed = TRUE)
    cfg(modifyList(RENDER_TEST_CFG, list(fontsize = 12, page_w = 6)))
    session$setInputs(apply_settings = 1)
  })
  expect_identical(drawn(rec), list(c(1000000L, 2000000L), c(1000000L, 2000000L)))
  expect_identical(rec$renders[[2]]$settings$fontsize, 12L)
  expect_identical(rec$renders[[2]]$settings$page_w, 6)
  expect_identical(rec$renders[[2]]$region_params$gene, "Pax5")
})

test_that("with a render process, the figure arrives when the process calls back", {
  w <- fake_worker()
  rec <- drive_render({
    !!generate
    expect_true(state$rendering)
    expect_null(state$svg_path)
    # A second Generate while it runs is ignored.
    session$setInputs(generate = 2)
    expect_length(w$renders, 1)
    job <- w$renders[[1]]
    writeLines("<svg></svg>", job$settings$output_path)
    job$done(NULL)
    session$flushReact()
    expect_false(state$rendering)
    expect_identical(state$svg_path, job$settings$output_path)
    expect_identical(state$rcode, "script 1")
  }, worker = w)
  expect_length(rec$renders, 0)   # never drawn in-process
  expect_identical(custom_types(rec), c("show_output", "hide_spinner"))
})

test_that("with a render process, an error it reports becomes the render's error", {
  w <- fake_worker()
  rec <- drive_render({
    !!generate
    w$renders[[1]]$done("The render process stopped.")
    session$flushReact()
    expect_identical(state$plot_error, "The render process stopped.")
    expect_false(state$rendering)
  }, worker = w)
  expect_length(rec$scripts, 0)
})

test_that("Cancel stops a running render, keeps the previous figure and says so", {
  w <- fake_worker()
  rec <- drive_render({
    !!generate
    session$setInputs(cancel_render = 1)
    expect_false(state$rendering)
    expect_null(state$svg_path)
  }, worker = w)
  expect_identical(w$cancels, 1)
  expect_identical(rec$notes, list(list(ui = "Rendering cancelled. The figure shown, if any, is the previous one.",
                                        type = "warning")))
  expect_identical(custom_types(rec), "hide_spinner")
})

test_that("Cancel with nothing running, or with no render process, only hides the spinner", {
  w <- fake_worker()
  rec <- drive_render(session$setInputs(cancel_render = 1), worker = w)
  expect_identical(w$cancels, 0)
  expect_length(rec$notes, 0)
  expect_identical(custom_types(rec), "hide_spinner")
  rec <- drive_render({
    state$rendering <- TRUE
    session$setInputs(cancel_render = 1)
    expect_true(state$rendering)
  })
  expect_length(rec$notes, 0)
  expect_identical(custom_types(rec), "hide_spinner")
})
