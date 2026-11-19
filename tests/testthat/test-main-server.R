# The main screen's server side (R/app-screen-server.R) and the screen switch
# (.pt_screen(), R/app-screen-ui.R), on the bundled mammary example (Basal and LP; chr2,
# mm10). setup-screen.R selects the classic screen for the whole test session, so each
# test here that needs the new screen asks for it locally. The data answers themselves
# are test-main-api.R's.

# The example's project and context, built once for this file.
main_srv <- local({
  memo <- NULL
  function() {
    if (is.null(memo)) {
      proj <- tempfile("pt-main-srv-")
      suppressWarnings(suppressMessages(pt_project(proj, pt_example("mammary"))))
      json <- file.path(normalizePath(proj, winslash = "/"), "project_settings.json")
      cfg <- .pt_read_settings(json, .pt_default_settings())
      memo <<- list(json = json, cfg = cfg,
                    ctx = suppressMessages(.pt_composer_context(json, cfg)))
    }
    memo
  }
})

# A render request as the page sends it, with any field replaced or dropped (NULL).
render_args <- function(...) {
  args <- list(gene = "Elf5", chr = "chr2", s = 103400000, e = 103460000,
               picked = list("LP", "Basal"),
               tracks = list(list(type = "atac", sample = "LP"), list(type = "genes")))
  over <- list(...)
  args[names(over)] <- over
  args
}

# The message a job stops with.
job_error <- function(ctx, args, title = "")
  expect_error(.pt_main_render_job(ctx, args, ctx$RENDER_CFG, title), class = "error")$message

# The custom messages of one type a recording session got, in the order sent.
sent <- function(rec, type)
  lapply(Filter(function(m) identical(m$type, type), rec$custom), `[[`, "message")

# The page's answers so far, as list(id, reply), the reply's JSON read back.
replies <- function(rec)
  lapply(sent(rec, "pt_rpc_reply"), function(m)
    list(id = m$id, reply = jsonlite::fromJSON(m$json, simplifyVector = FALSE)))

# A pt_rpc input value: the requests as the JSON array the page sends.
rpc <- function(...) as.character(jsonlite::toJSON(list(...), auto_unbox = TRUE, digits = NA))

DISPLAY_KEYS <- setdiff(names(.pt_default_settings()), .PT_CORE_SETTINGS)

# ---- .pt_screen() ------------------------------------------------------------------------

test_that("a test run drives the classic screen unless a test asks for the new one", {
  # setup-screen.R: what the older composer tests rely on.
  expect_identical(getOption("PlotTrellis.screen"), "classic")
  expect_identical(.pt_screen(), "classic")
})

test_that(".pt_screen(): new by default, the option, else the environment variable", {
  withr::local_options(PlotTrellis.screen = NULL)
  withr::local_envvar(PLOTTRELLIS_SCREEN = NA)
  expect_identical(.pt_screen(), "new")
  withr::local_envvar(PLOTTRELLIS_SCREEN = "classic")
  expect_identical(.pt_screen(), "classic")
  withr::local_options(PlotTrellis.screen = "new")
  expect_identical(.pt_screen(), "new")
  withr::local_envvar(PLOTTRELLIS_SCREEN = "new")
  withr::local_options(PlotTrellis.screen = "classic")
  expect_identical(.pt_screen(), "classic")
})

test_that(".pt_screen(): any value but \"classic\" means the new screen", {
  withr::local_envvar(PLOTTRELLIS_SCREEN = NA)
  for (v in list("Classic", "old", "", NA, TRUE, c("classic", "new"), 1))
    expect_identical(withr::with_options(list(PlotTrellis.screen = v), .pt_screen()), "new")
  withr::local_options(PlotTrellis.screen = NULL)
  withr::local_envvar(PLOTTRELLIS_SCREEN = "CLASSIC")
  expect_identical(.pt_screen(), "new")
})

# ---- .pt_composer_app() picks the page ---------------------------------------------------

test_that(".pt_composer_app() gives the new page by default and the classic one with the option", {
  m <- main_srv()
  page_html <- function(app) {
    ui <- environment(app$httpHandler)$ui
    htmltools::renderTags(if (is.function(ui)) ui(NULL) else ui)$html
  }
  withr::local_envvar(PLOTTRELLIS_SCREEN = NA)
  new_app <- withr::with_options(list(PlotTrellis.screen = NULL),
                                 suppressMessages(.pt_composer_app(m$json, m$cfg)))
  expect_s3_class(new_app, "shiny.appobj")
  html <- page_html(new_app)
  expect_match(html, 'id="ptApp"', fixed = TRUE)
  expect_match(html, 'id="ptSettings"', fixed = TRUE)
  expect_false(grepl('id="view_mode"', html, fixed = TRUE))

  classic <- withr::with_options(list(PlotTrellis.screen = "classic"),
                                 suppressMessages(.pt_composer_app(m$json, m$cfg)))
  expect_s3_class(classic, "shiny.appobj")
  html <- page_html(classic)
  expect_match(html, 'id="view_mode"', fixed = TRUE)
  expect_false(grepl('id="ptApp"', html, fixed = TRUE))
})

# ---- .pt_main_settings_key() -------------------------------------------------------------

test_that("the settings key changes with any display setting or the title, and only with those", {
  cfg <- .pt_default_settings()
  key <- function(cfg, title = "") .pt_main_settings_key(cfg, DISPLAY_KEYS, title)
  k0 <- key(cfg)
  expect_type(k0, "character")
  expect_length(k0, 1)
  expect_identical(key(cfg), k0)
  # The title, trimmed; NULL is no title.
  expect_identical(key(cfg, "  "), k0)
  expect_identical(key(cfg, NULL), k0)
  expect_false(identical(key(cfg, "Elf5"), k0))
  expect_identical(key(cfg, " Elf5 "), key(cfg, "Elf5"))
  # One display setting.
  for (k in c("fontsize", "page_w", "show_title", "reference_colour", "uniprot_reviewed_only")) {
    x <- cfg
    x[[k]] <- if (is.logical(x[[k]])) !x[[k]] else if (is.numeric(x[[k]])) x[[k]] + 1 else "#010203"
    expect_false(identical(key(x), k0), info = k)
  }
  # A core setting is not a display setting.
  x <- cfg
  x$data_dir <- "/elsewhere"
  expect_identical(key(x), k0)
  # A NULL ("auto") setting: no error, and not the same key as a set value.
  expect_null(cfg$zrange_cap)
  x <- cfg
  x["zrange_cap"] <- list(5)
  expect_false(identical(key(x), k0))
  x["zrange_cap"] <- list(NULL)
  expect_identical(key(x), k0)
  # Two auto settings, one set at a time, give different keys.
  a <- cfg; a["zrange_cap"] <- list(5)
  b <- cfg; b["hic_height"] <- list(5)
  expect_false(identical(key(a), key(b)))
})

# ---- .pt_main_render_job() ---------------------------------------------------------------

test_that("a render job draws what the classic screen would for the same choice", {
  skip_if_not_installed("org.Mm.eg.db")
  ctx <- main_srv()$ctx
  cfg <- ctx$RENDER_CFG
  cfg$fontsize <- 11L
  args <- render_args(s = 103400000.4, e = 103460000.6)
  job <- .pt_main_render_job(ctx, args, cfg, "  My title ")
  expect_named(job, c("entries", "region_params", "settings"))
  labels <- .pt_main_labels(ctx, args$picked, args$tracks)
  expect_identical(labels, c("bw_atac_LP", "genes"))
  expect_identical(job$entries, .pt_build_entries(labels, ctx$TRACK_ITEMS, ctx$DEFAULT_SAMPLE_COLS, cfg))
  expect_identical(job$region_params,
                   list(chrom = "chr2", start = 103400000L, end = 103460001L, gene = "Elf5"))
  expect_identical(job$settings$title, "My title")
  expect_identical(job$settings$fontsize, 11L)
  expect_identical(job$settings[names(.pt_settings_base(cfg))], .pt_settings_base(cfg))
  expect_match(job$settings$output_path, "\\.svg$")
  expect_false(file.exists(job$settings$output_path))
  # Each job draws to its own file.
  expect_false(identical(.pt_main_render_job(ctx, args, cfg)$settings$output_path,
                         job$settings$output_path))
})

test_that("a render job clamps the region to the chromosome and gives a region no gene", {
  skip_if_not_installed("org.Mm.eg.db")
  ctx <- main_srv()$ctx
  len <- ctx$CHR_LENGTHS[["chr2"]]
  rp <- function(...) .pt_main_render_job(ctx, render_args(...), ctx$RENDER_CFG)$region_params
  r <- rp(s = -500, e = len + 1e6)
  expect_identical(r[c("start", "end")], list(start = 1L, end = as.integer(len)))
  # Numbers sent as text are read as numbers.
  expect_identical(rp(s = "1000", e = "2000")[c("start", "end")], list(start = 1000L, end = 2000L))
  # A region: no gene, an empty gene or a non-text gene is NA.
  for (g in list(NULL, "", NA, 5, list("Elf5", "Ehf"))) {
    args <- render_args()
    args["gene"] <- list(g)
    expect_identical(.pt_main_render_job(ctx, args, ctx$RENDER_CFG)$region_params$gene, NA_character_)
  }
  # No title is an empty title.
  expect_identical(.pt_main_render_job(ctx, render_args(), ctx$RENDER_CFG, NULL)$settings$title, "")
})

test_that("a render job that cannot be drawn says why", {
  skip_if_not_installed("org.Mm.eg.db")
  ctx <- main_srv()$ctx
  expect_identical(job_error(ctx, render_args(tracks = list())), "Add at least one track.")
  expect_identical(job_error(ctx, render_args(tracks = NULL)), "Add at least one track.")
  expect_identical(job_error(ctx, render_args(chr = "chrZ")), "No chromosome 'chrZ'.")
  expect_identical(job_error(ctx, render_args(chr = NULL)), "No chromosome ''.")
  expect_identical(job_error(ctx, render_args(chr = list("chr1", "chr2"))), "No chromosome 'chr1chr2'.")
  needs <- "The region needs a start and an end."
  expect_identical(job_error(ctx, render_args(s = NULL)), needs)
  expect_identical(job_error(ctx, render_args(e = NULL)), needs)
  expect_identical(job_error(ctx, render_args(s = "start")), needs)
  expect_identical(job_error(ctx, render_args(e = list(1, 2))), needs)
  # .pt_region_error()'s cases, word for word.
  for (se in list(c(1000, 1000), c(2000, 1000), c(1000, Inf), c(1, 3e9)))
    expect_identical(job_error(ctx, render_args(s = se[1], e = se[2])), .pt_region_error(se[1], se[2]))
  expect_match(job_error(ctx, render_args(s = 2000, e = 1000)), "must be after its start", fixed = TRUE)
  expect_match(job_error(ctx, render_args(s = 1, e = 3e9)), "must be below", fixed = TRUE)
  # .pt_main_labels()'s errors pass through, before the region is read.
  expect_identical(job_error(ctx, render_args(tracks = list(list(type = "nope")), chr = "chrZ")),
                   "Unknown track type 'nope'.")
  expect_identical(job_error(ctx, render_args(tracks = list(list(type = "atac")))),
                   "'atac' is drawn per sample: say which sample.")
  expect_identical(job_error(ctx, render_args(picked = list("LP"), tracks = list(list(type = "atac_diff")))),
                   "'atac_diff' compares two samples: pick exactly two.")
})

test_that("a render job takes a gene only as one of the project's own, under its symbol", {
  skip_if_not_installed("org.Mm.eg.db")
  ctx <- main_srv()$ctx
  gene_of <- function(g) .pt_main_render_job(ctx, render_args(gene = g), ctx$RENDER_CFG)$region_params$gene
  expect_identical(gene_of("elf5"), "Elf5")
  # Text that is no gene never reaches the render (a track module builds a cache file's
  # path from the gene's name): not a path, not a made-up name.
  for (g in c("../../../somewhere/else", "/tmp/x", "Elf5/../x", "NoSuchGene123"))
    expect_identical(gene_of(g), NA_character_, info = g)
  # Nor does a symbol the organism has but this project's annotation does not (the example
  # annotates part of chr2 only): some symbols contain a slash.
  outside <- setdiff(ctx$ALL_GENES, .pt_main_gene_index(ctx)$genes$sym)
  expect_gt(length(outside), 0)
  for (g in c(outside[1], grep("/", outside, value = TRUE, fixed = TRUE)[1]))
    if (!is.na(g)) expect_identical(gene_of(g), NA_character_, info = g)
})

test_that("a render job refuses what the classic screen will not generate", {
  skip_if_not_installed("org.Mm.eg.db")
  ctx <- main_srv()$ctx
  tab <- .pt_track_table(ctx$TRACK_ITEMS)
  # A gene-only track for a region (no gene).
  needs <- vapply(ctx$TRACK_ITEMS[tab$key], function(it) it$type %in% ctx$GENE_REQUIRED_TYPES, logical(1))
  type <- tab$track[needs & tab$kind == "fixed"][1]
  expect_false(is.na(type))
  one <- list(list(type = type))
  expect_match(job_error(ctx, render_args(gene = NULL, tracks = one)), "^Remove before generating: .*needs a gene")
  expect_named(.pt_main_render_job(ctx, render_args(tracks = one), ctx$RENDER_CFG),
               c("entries", "region_params", "settings"))
  # A track whose file is missing.
  broken <- ctx
  broken$TRACK_ITEMS[["bw_atac_LP"]]$broken <- TRUE
  expect_identical(job_error(broken, render_args()),
                   paste0("Remove before generating: ", ctx$TRACK_ITEMS[["bw_atac_LP"]]$label,
                          " (its file is missing)."))
})

test_that("a chromosome with genes but no chromosome size is drawn up to its last gene", {
  skip_if_not_installed("org.Mm.eg.db")
  ctx <- main_srv()$ctx
  bare <- ctx
  bare$CHR_LENGTHS <- ctx$CHR_LENGTHS[setdiff(names(ctx$CHR_LENGTHS), "chr2")]
  g <- .pt_main_gene_index(ctx)$genes
  last <- max(g$e[g$chr == "chr2"])
  expect_identical(.pt_main_chr_length(bare, "chr2"), as.numeric(last))
  expect_identical(.pt_main_chr_length(ctx, "chr2"), as.numeric(ctx$CHR_LENGTHS[["chr2"]]))
  for (x in list("chrZ", NA_character_, NULL, c("chr1", "chr2"), 2))
    expect_identical(.pt_main_chr_length(bare, x), NA_real_)
  r <- .pt_main_render_job(bare, render_args(s = 103400000, e = 3e8), ctx$RENDER_CFG)$region_params
  expect_identical(r$end, as.integer(last))
})

test_that("a gene's view follows the live gene flank, not the one at launch", {
  skip_if_not_installed("org.Mm.eg.db")
  ctx <- main_srv()$ctx
  wide <- ctx$RENDER_CFG
  wide$gene_flank_pct <- 200
  a <- .pt_main_locus(ctx, "Elf5")$view
  b <- .pt_main_locus(ctx, "Elf5", wide)$view
  expect_lt(b$s, a$s)
  expect_gt(b$e, a$e)
  expect_identical(.pt_main_request(ctx, "locus", list(q = "Elf5"), wide)$view, b)
  # The page cannot supply settings of its own.
  expect_identical(.pt_main_request(ctx, "locus", list(q = "Elf5", cfg = wide))$view, a)
})

# ---- .pt_main_closer() -------------------------------------------------------------------

# A closer whose timers are held (run them with fire()) and whose stop is counted.
held_closer <- function() {
  timers <- list()
  stops <- 0L
  closer <- .pt_main_closer(grace = 8, stop = function() stops <<- stops + 1L,
                            later = function(f, delay) timers[[length(timers) + 1L]] <<- f)
  c(closer, list(fire = function() { for (f in timers) f(); timers <<- list() },
                 pending = function() length(timers), stops = function() stops))
}

test_that("the app stops once, after the last page is closed by its user", {
  k <- held_closer()
  k$opened()
  k$ended(by_user = TRUE)
  expect_identical(k$pending(), 1L)
  expect_identical(k$stops(), 0L)   # not before the grace period is over
  k$fire()
  expect_identical(k$stops(), 1L)
  # Two pages closed together: still one stop.
  k <- held_closer()
  k$opened(); k$opened()
  k$ended(TRUE); k$ended(TRUE)
  k$fire()
  expect_identical(k$stops(), 1L)
})

test_that("a page still open, a reload or a lost connection does not stop the app", {
  k <- held_closer()
  k$opened(); k$opened()
  k$ended(TRUE)                     # one of two tabs closed
  k$fire()
  expect_identical(k$stops(), 0L)
  k <- held_closer()
  k$opened()
  k$ended(TRUE); k$opened()         # a reload: the new page arrives within the grace period
  k$fire()
  expect_identical(k$stops(), 0L)
  k <- held_closer()
  k$opened(); k$opened()
  k$ended(TRUE)
  k$opened(); k$ended(TRUE)         # the other order: the new page first, then the old one ends
  k$fire()
  expect_identical(k$stops(), 0L)
  k <- held_closer()
  k$opened()
  k$ended(by_user = FALSE)          # cut off, not closed: no timer at all
  expect_identical(k$pending(), 0L)
  expect_identical(k$stops(), 0L)
})

test_that("a timer left from an app that has already stopped does not stop the next one", {
  k <- held_closer()
  k$opened()
  k$ended(TRUE)
  k$stopped()                       # Edit project settings, Reset: the app stopped by itself
  k$fire()
  expect_identical(k$stops(), 0L)
})

# ---- .pt_main_figure() -------------------------------------------------------------------

test_that("the figure is the SVG's text and its script; the temporary file is removed", {
  skip_on_cran()
  skip_if_not_installed("org.Mm.eg.db")
  ctx <- main_srv()$ctx
  job <- .pt_main_render_job(ctx, render_args(), ctx$RENDER_CFG, "Fig")
  path <- job$settings$output_path
  on.exit(unlink(path), add = TRUE)
  res <- .pt_render_svg(job$entries, job$region_params, job$settings)
  expect_type(res, "list")          # drawn (the error message otherwise)
  expect_true(file.exists(path))
  fig <- .pt_main_figure(job, res$hic_tops)
  expect_named(fig, c("svg", "script"))
  expect_false(file.exists(path))
  expect_type(fig$svg, "character")
  expect_length(fig$svg, 1)
  expect_match(fig$svg, "^<svg|<svg ")
  expect_match(fig$svg, "Elf5", fixed = TRUE)
  expect_type(fig$script, "character")
  expect_length(fig$script, 1)
  expect_match(fig$script, "render_all", fixed = TRUE)
  expect_match(fig$script, "chr2", fixed = TRUE)
  expect_false(grepl(basename(path), fig$script, fixed = TRUE))
  expect_false(grepl(path, fig$script, fixed = TRUE))
})

# ---- the server --------------------------------------------------------------------------
# Custom messages are captured by recording_session() (helper-composer-fakes.R): a
# MockShinySession whose sendCustomMessage is replaced by one that records each
# list(type, message), passed to testServer() as its `session`.

test_that("one pt_rpc value with several requests gets one reply each, with its id", {
  skip_if_not_installed("org.Mm.eg.db")
  m <- main_srv()
  rec <- recording_session()
  shiny::testServer(.pt_main_app(m$ctx, m$json), {
    session$flushReact()
    session$setInputs(pt_rpc = rpc(list(id = 1, op = "project"),
                                   list(id = "two", op = "search", args = list(q = "el")),
                                   list(id = 3, op = "no_such_op"),
                                   list(id = 4, op = "cancel")))
    r <- replies(rec)
    expect_identical(lapply(r, `[[`, "id"), list(1L, "two", 3L, 4L))
    p <- r[[1]]$reply
    expect_identical(p$id, basename(m$ctx$PROJECT_DIR))
    expect_identical(p$version, as.character(utils::packageVersion("PlotTrellis")))
    expect_true(all(c("genome", "samples", "types", "chrLen") %in% names(p)))
    expect_null(r[[2]]$reply$error)
    expect_type(r[[3]]$reply$error, "character")
    # Nothing drawing: cancel just answers.
    expect_identical(r[[4]]$reply, list(ok = TRUE))
    expect_length(sent(rec, "pt_rpc_reply"), 4)
  }, session = rec$session)
})

test_that("malformed JSON and odd requests do not stop the session", {
  skip_if_not_installed("org.Mm.eg.db")
  m <- main_srv()
  rec <- recording_session()
  shiny::testServer(.pt_main_app(m$ctx, m$json), {
    session$flushReact()
    session$setInputs(pt_rpc = "[{\"id\": 1, \"op\": ")
    expect_length(replies(rec), 0)
    session$setInputs(pt_rpc = "[]")
    expect_length(replies(rec), 0)
    # No op, then the session still answers.
    session$setInputs(pt_rpc = rpc(list(id = 5), list(id = 6, op = "cancel")))
    r <- replies(rec)
    expect_identical(lapply(r, `[[`, "id"), list(5L, 6L))
    expect_type(r[[1]]$reply$error, "character")
    expect_identical(r[[2]]$reply, list(ok = TRUE))
  }, session = rec$session)
})

test_that("a pt_rpc value that is a JSON object, not an array, does not stop the session", {
  skip_if_not_installed("org.Mm.eg.db")
  m <- main_srv()
  rec <- recording_session()
  shiny::testServer(.pt_main_app(m$ctx, m$json), {
    session$flushReact()
    expect_no_error(session$setInputs(pt_rpc = "{\"id\": 1, \"op\": \"cancel\"}"))
    session$setInputs(pt_rpc = rpc(list(id = 2, op = "cancel")))
    expect_identical(replies(rec)[[length(replies(rec))]], list(id = 2L, reply = list(ok = TRUE)))
  }, session = rec$session)
})

test_that("geneinfo is answered with the live settings, and its error as {error}", {
  skip_if_not_installed("org.Mm.eg.db")
  m <- main_srv()
  rec <- recording_session()
  local_mocked_bindings(.pt_main_geneinfo = function(ctx, gene, cfg) {
    if (identical(gene, "bad")) stop("no such gene")
    list(gene = gene, reviewed = cfg$uniprot_reviewed_only)
  })
  shiny::testServer(.pt_main_app(m$ctx, m$json), {
    session$flushReact()
    session$setInputs(s_uniprot_reviewed_only = TRUE)
    session$setInputs(pt_rpc = rpc(list(id = 1, op = "geneinfo", args = list(gene = "Elf5")),
                                   list(id = 2, op = "geneinfo", args = list(gene = "bad"))))
    r <- replies(rec)
    expect_identical(r[[1]], list(id = 1L, reply = list(gene = "Elf5", reviewed = TRUE)))
    expect_identical(r[[2]], list(id = 2L, reply = list(error = "no such gene")))
  }, session = rec$session)
})

test_that("render in-process answers with the SVG and script; one that cannot be drawn with {error}", {
  skip_on_cran()
  skip_if_not_installed("org.Mm.eg.db")
  m <- main_srv()
  rec <- recording_session()
  shiny::testServer(.pt_main_app(m$ctx, m$json), {
    session$flushReact()
    session$setInputs(plot_title = "Title here")
    session$setInputs(pt_rpc = rpc(list(id = "r1", op = "render", args = render_args()),
                                   list(id = "r2", op = "render", args = render_args(tracks = list())),
                                   list(id = "r3", op = "render", args = render_args(chr = "chrZ"))))
    r <- replies(rec)
    expect_identical(lapply(r, `[[`, "id"), list("r1", "r2", "r3"))
    fig <- r[[1]]$reply
    expect_named(fig, c("svg", "script", "key"))
    # The key is the settings the figure was drawn with (the title included), as the
    # "pt_settings" message gives it: the page files the figure under it.
    expect_identical(fig$key, .pt_main_settings_key(m$cfg, DISPLAY_KEYS, "Title here"))
    expect_identical(fig$key, sent(rec, "pt_settings")[[length(sent(rec, "pt_settings"))]]$key)
    expect_match(fig$svg, "<svg", fixed = TRUE)
    expect_match(fig$svg, "Title here", fixed = TRUE)
    expect_match(fig$script, "render_all", fixed = TRUE)
    expect_identical(r[[2]]$reply, list(error = "Add at least one track."))
    expect_identical(r[[3]]$reply, list(error = "No chromosome 'chrZ'."))
  }, session = rec$session)
})

test_that("a render that fails answers {error} and removes its file", {
  skip_if_not_installed("org.Mm.eg.db")
  m <- main_srv()
  rec <- recording_session()
  paths <- character(0)
  local_mocked_bindings(.pt_render_svg = function(entries, region_params, settings) {
    paths <<- c(paths, settings$output_path)
    writeLines("<svg/>", settings$output_path)
    "the figure failed"
  })
  shiny::testServer(.pt_main_app(m$ctx, m$json), {
    session$flushReact()
    session$setInputs(pt_rpc = rpc(list(id = 1, op = "render", args = render_args()),
                                   list(id = 2, op = "cancel")))
    r <- replies(rec)
    expect_identical(r[[1]], list(id = 1L, reply = list(error = "the figure failed")))
    expect_identical(r[[2]], list(id = 2L, reply = list(ok = TRUE)))
    expect_length(paths, 1)
    expect_false(file.exists(paths))
  }, session = rec$session)
})

test_that("pt_settings is sent with a key, and a new key after a setting or the title changes", {
  skip_if_not_installed("org.Mm.eg.db")
  m <- main_srv()
  rec <- recording_session()
  shiny::testServer(.pt_main_app(m$ctx, m$json), {
    session$flushReact()
    keys <- function() vapply(sent(rec, "pt_settings"), `[[`, "", "key")
    k <- keys()
    expect_gte(length(k), 1)
    k0 <- k[length(k)]
    expect_identical(k0, .pt_main_settings_key(m$ctx$RENDER_CFG, DISPLAY_KEYS, ""))
    session$setInputs(s_fontsize = m$ctx$RENDER_CFG$fontsize + 3L)
    k1 <- keys()[length(keys())]
    expect_false(identical(k1, k0))
    session$setInputs(plot_title = "A title")
    k2 <- keys()[length(keys())]
    expect_false(k2 %in% c(k0, k1))
    # Back to the earlier settings and title: the earlier key.
    session$setInputs(s_fontsize = m$ctx$RENDER_CFG$fontsize, plot_title = "")
    expect_identical(keys()[length(keys())], k0)
  }, session = rec$session)
})
