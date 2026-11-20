# Characterisation of the app server's shared pieces, as wired into the main screen's server
# (.pt_main_server(), R/app-screen-server.R): the Figure settings block's live and saved
# settings, the live settings reaching a gene's view and a render, the repeat-click gate,
# and the unsaved-settings warning of the project dialogs. The block's own checks with no
# genome are test-composer-settings-block.R's; the page's data answers are test-main-api.R's.

need_fixture_data <- function()
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")

# A fresh project (json + its folder) on the mouse fixture, and its startup context.
char_project <- function(extra = list()) {
  proj <- tempfile("pt-char-"); dir.create(proj)
  seed_project_genome(proj)
  json <- file.path(proj, "project_settings.json")
  jsonlite::write_json(c(list(data_dir = FIXTURES$mouse$data_dir, samples = FIXTURES$mouse$samples,
                              gtf_path = GENOMES$mouse$gtf_path, orgdb = GENOMES$mouse$orgdb,
                              genome_label = GENOMES$mouse$genome_label), extra),
                       json, auto_unbox = TRUE, null = "null")
  list(proj = proj, json = json,
       ctx = suppressMessages(.pt_composer_context(json, .pt_read_settings(json, .pt_default_settings()))))
}

# One project shared by the tests that never write its files (each testServer() call is
# its own session, so they share no state). Built on first use.
.char_shared <- new.env()
shared_project <- function() {
  if (is.null(.char_shared$cp)) {
    .char_shared$cp <- char_project()
    withr::defer(unlink(.char_shared$cp$proj, recursive = TRUE), envir = testthat::teardown_env())
  }
  .char_shared$cp
}

# The main screen's app for a project.
char_app <- function(cp) suppressMessages(.pt_main_app(cp$ctx, cp$json))

# A pt_rpc input value: the requests as the JSON array the page sends.
char_rpc <- function(...) as.character(jsonlite::toJSON(list(...), auto_unbox = TRUE, digits = NA))

# The custom messages of one type a recording session got, in the order sent.
char_sent <- function(rec, type)
  lapply(Filter(function(m) identical(m$type, type), rec$custom), `[[`, "message")

# The reply to request `id` (a string), its JSON read back; an error if there is none.
char_reply <- function(rec, id) {
  r <- Filter(function(m) identical(m$id, id), char_sent(rec, "pt_rpc_reply"))
  if (!length(r)) stop("no reply to request '", id, "'", call. = FALSE)
  jsonlite::fromJSON(r[[length(r)]]$json, simplifyVector = FALSE)
}

# The Figure settings block alone, on the default settings (fontsize 8L, page_w 9, font
# "DejaVu Sans", show_title TRUE). Nothing here saves, so the settings file is never read.
block_app <- function(input, output, session)
  fs <- .pt_figure_settings_server(input = input, output = output, render_cfg = .pt_default_settings(),
                                   settings_json = file.path(tempdir(), "pt-char-no-such-project.json"))

# ---- Figure settings ------------------------------------------------------------------

test_that("display_keys are the non-core settings, and live_cfg starts as saved_cfg", {
  need_fixture_data()
  cp <- shared_project()
  shiny::testServer(char_app(cp), {
    session$flushReact()
    expect_identical(fs$display_keys, setdiff(names(.pt_default_settings()), .PT_CORE_SETTINGS))
    expect_identical(fs$live_cfg()[fs$display_keys], fs$saved_cfg()[fs$display_keys])
    expect_false(fs$settings_dirty())
    expect_null(output$settings_dirty_ui)
    expect_identical(fs$settings_ui_gen(), 0)
    expect_identical(fs$saved_cfg()$cache_dir, file.path(cp$proj, "cache"))
  })
})

test_that("settings_dirty compares numbers by value, everything else exactly", {
  shiny::testServer(block_app, {
    session$setInputs(s_fontsize = 8, s_page_w = 9L)  # saved 8L and 9
    expect_false(fs$settings_dirty())
    session$setInputs(s_page_w = 9 + 1e-10)
    expect_false(fs$settings_dirty())
    session$setInputs(s_page_w = 9.01)
    expect_true(fs$settings_dirty())
    expect_identical(ui_html(output$settings_dirty_ui),
                     '<span class="text-warning">Unsaved changes (this session only).</span>')
    session$setInputs(s_page_w = 9, s_zrange_cap = NA, s_hic_height = NA)  # empty optional boxes are "auto"
    expect_false(fs$settings_dirty())
    session$setInputs(s_font = " DejaVu Sans ")  # trimmed
    expect_false(fs$settings_dirty())
    session$setInputs(s_fontsize = "8")  # a string is not the number
    expect_true(fs$settings_dirty())
    session$setInputs(s_fontsize = 8, s_show_title = FALSE)
    expect_true(fs$settings_dirty())
    session$setInputs(s_show_title = TRUE)
    expect_false(fs$settings_dirty())
    expect_null(output$settings_dirty_ui)
  })
})

test_that("the settings panel holds every s_ input and the comparison group", {
  shiny::testServer(block_app, {
    html <- ui_html(output$display_settings_ui)
    for (k in fs$display_keys) expect_match(html, sprintf('id="s_%s"', k), fixed = TRUE, info = k)
    expect_match(html, "Comparison (diff) tracks", fixed = TRUE)
  })
})

# ---- Live settings in the page's requests -----------------------------------------------

test_that("the gene flank comes from s_gene_flank_pct, else the saved setting", {
  need_fixture_data()
  cp <- char_project(list(gene_flank_pct = 50))
  on.exit(unlink(cp$proj, recursive = TRUE), add = TRUE)
  rec <- recording_session()
  shiny::testServer(char_app(cp), {
    session$flushReact()
    n <- 0
    view <- function() {
      n <<- n + 1
      id <- paste0("q", n)
      session$setInputs(pt_rpc = char_rpc(list(id = id, op = "locus", args = list(q = "Elf5"))))
      char_reply(rec, id)$view
    }
    saved50 <- view()
    expect_identical(saved50, list(s = 103400000L, e = 103460000L))
    session$setInputs(s_gene_flank_pct = 15)
    expect_identical(view(), list(s = 103410000L, e = 103455000L))
    for (typed in list(NA, "abc")) {  # an empty box, or junk: the saved value
      session$setInputs(s_gene_flank_pct = typed)
      expect_identical(view(), saved50, info = deparse(typed))
    }
  }, session = rec$session)
})

test_that("a render uses the live settings, comparison role colours included", {
  need_fixture_data()
  cp <- shared_project()
  rec <- recording_session()
  seen <- list()
  local_mocked_bindings(.pt_render_svg = function(entries, region_params, settings) {
    seen[[length(seen) + 1]] <<- list(entries = entries, settings = settings)
    writeLines("<svg width='10' height='10'><rect/></svg>", settings$output_path)
    list(hic_tops = NULL)
  })
  # picked[1] is the comparison, picked[2] the reference.
  args <- list(chr = "chr2", s = 103300000, e = 103500000, picked = list("LP", "Ba"),
               tracks = list(list(type = "atac_diff"), list(type = "genes")))
  shiny::testServer(char_app(cp), {
    session$flushReact()
    session$setInputs(pt_rpc = char_rpc(list(id = "r1", op = "render", args = args)))
    expect_named(char_reply(rec, "r1"), c("svg", "script", "key", "job"))
    expect_length(seen, 1)
    e <- seen[[1]]$entries[[1]]
    expect_identical(e[c("viz", "control_label", "treatment_label", "control_colour", "treatment_colour")],
                     list(viz = "bw_diff_stacked", control_label = "Ba", treatment_label = "LP",
                          control_colour = "#4D4D4D", treatment_colour = "#56B4E9"))
    session$setInputs(s_reference_colour = "#123456", s_fontsize = 14)
    session$setInputs(pt_rpc = char_rpc(list(id = "r2", op = "render", args = args)))
    expect_named(char_reply(rec, "r2"), c("svg", "script", "key", "job"))
    expect_length(seen, 2)
    expect_identical(seen[[2]]$entries[[1]]$control_colour, "#123456")
    expect_equal(seen[[2]]$settings$fontsize, 14)
  }, session = rec$session)
})

# ---- Small server pieces ----------------------------------------------------------------

test_that("a pt-once click is released after the flush that handled it", {
  need_fixture_data()
  rec <- recording_session()
  shiny::testServer(char_app(shared_project()), {
    session$flushReact()
    session$setInputs(pt_once_click = "save_settings")
    got <- char_sent(rec, "pt_once_ready")
    expect_identical(got[[length(got)]], "save_settings")
  }, session = rec$session)
})

test_that("Edit project settings and Clear all caches warn about unsaved Figure settings only when dirty", {
  need_fixture_data()
  modals <- list()
  local_mocked_bindings(showModal = function(ui, ...) modals[[length(modals) + 1]] <<- as.character(ui))
  shiny::testServer(char_app(shared_project()), {
    session$flushReact()
    session$setInputs(open_settings = 1)
    expect_false(grepl("not saved and will be lost", modals[[1]], fixed = TRUE))
    session$setInputs(s_fontsize = 14, open_settings = 2)
    expect_match(modals[[2]], "Your Figure settings changes are not saved and will be lost.", fixed = TRUE)
    session$setInputs(clear_caches_btn = 1)
    expect_match(modals[[3]], "not saved and will be lost", fixed = TRUE)
  })
})
