# The main screen's server side (R/app-screen-server.R), on the bundled mammary example
# (Basal and LP; five genes, mm39). The data answers themselves are test-main-api.R's.

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
  args <- list(gene = "Eya2", chr = "chr2", s = 165430000, e = 165620000,
               picked = list("LP", "Basal"),
               tracks = list(list(type = "ATAC", sample = "LP"), list(type = "genes")))
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

# ---- .pt_composer_app() builds the page --------------------------------------------------

test_that(".pt_composer_app() gives the main screen's page", {
  m <- main_srv()
  app <- suppressMessages(.pt_composer_app(m$json, m$cfg))
  expect_s3_class(app, "shiny.appobj")
  ui <- environment(app$httpHandler)$ui
  html <- htmltools::renderTags(if (is.function(ui)) ui(NULL) else ui)$html
  expect_match(html, 'id="ptApp"', fixed = TRUE)
  expect_match(html, 'id="ptSettings"', fixed = TRUE)
  # nothing of the earlier screen: its mode switch, its Close app button
  expect_false(grepl('id="view_mode"', html, fixed = TRUE))
  expect_false(grepl("close_app_btn", html, fixed = TRUE))
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
  expect_false(identical(key(cfg, "Eya2"), k0))
  expect_identical(key(cfg, " Eya2 "), key(cfg, "Eya2"))
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

test_that("a render job is the choice's entries, its region and the live settings", {
  skip_if_not_installed("org.Mm.eg.db")
  ctx <- main_srv()$ctx
  cfg <- ctx$RENDER_CFG
  cfg$fontsize <- 11L
  args <- render_args(s = 165430000.4, e = 165620000.6)
  job <- .pt_main_render_job(ctx, args, cfg, "  My title ")
  expect_named(job, c("entries", "region_params", "settings"))
  labels <- .pt_main_labels(ctx, args$picked, args$tracks)
  expect_identical(labels, c("bw_ATAC_LP", "genes"))
  expect_identical(job$entries, .pt_build_entries(labels, ctx$TRACK_ITEMS, ctx$DEFAULT_SAMPLE_COLS, cfg))
  expect_identical(job$region_params,
                   list(chrom = "chr2", start = 165430000L, end = 165620001L, gene = "Eya2"))
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
  for (g in list(NULL, "", NA, 5, list("Eya2", "Acot7"))) {
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
  expect_identical(job_error(ctx, render_args(tracks = list(list(type = "ATAC")))),
                   "'ATAC' is drawn per sample: say which sample.")
  expect_identical(job_error(ctx, render_args(picked = list("LP"), tracks = list(list(type = "ATAC_diff")))),
                   "'ATAC_diff' compares two samples: pick exactly two.")
})

test_that("a render job takes a gene only as one of the project's own, under its symbol", {
  skip_if_not_installed("org.Mm.eg.db")
  ctx <- main_srv()$ctx
  gene_of <- function(g) .pt_main_render_job(ctx, render_args(gene = g), ctx$RENDER_CFG)$region_params$gene
  expect_identical(gene_of("eya2"), "Eya2")
  # Text that is no gene never reaches the render (a track module builds a cache file's
  # path from the gene's name): not a path, not a made-up name.
  for (g in c("../../../somewhere/else", "/tmp/x", "Eya2/../x", "NoSuchGene123"))
    expect_identical(gene_of(g), NA_character_, info = g)
  # Nor does a symbol the organism has but this project's annotation does not (the example
  # annotates five genes' regions only): some symbols contain a slash.
  outside <- setdiff(ctx$ALL_GENES, .pt_main_gene_index(ctx)$genes$sym)
  expect_gt(length(outside), 0)
  for (g in c(outside[1], grep("/", outside, value = TRUE, fixed = TRUE)[1]))
    if (!is.na(g)) expect_identical(gene_of(g), NA_character_, info = g)
})

test_that("a render job refuses a track whose file is missing, and a gene-only track for a region", {
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
  broken$TRACK_ITEMS[["bw_ATAC_LP"]]$broken <- TRUE
  expect_identical(job_error(broken, render_args()),
                   paste0("Remove before generating: ", ctx$TRACK_ITEMS[["bw_ATAC_LP"]]$label,
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
  r <- .pt_main_render_job(bare, render_args(s = 165430000, e = 3e8), ctx$RENDER_CFG)$region_params
  expect_identical(r$end, as.integer(last))
})

test_that("a gene's view follows the live gene flank, not the one at launch", {
  skip_if_not_installed("org.Mm.eg.db")
  ctx <- main_srv()$ctx
  wide <- ctx$RENDER_CFG
  wide$gene_flank_pct <- 200
  a <- .pt_main_locus(ctx, "Eya2")$view
  b <- .pt_main_locus(ctx, "Eya2", wide)$view
  expect_lt(b$s, a$s)
  expect_gt(b$e, a$e)
  expect_identical(.pt_main_request(ctx, "locus", list(q = "Eya2"), wide)$view, b)
  # The page cannot supply settings of its own.
  expect_identical(.pt_main_request(ctx, "locus", list(q = "Eya2", cfg = wide))$view, a)
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
  expect_match(fig$svg, "Eya2", fixed = TRUE)
  expect_type(fig$script, "character")
  expect_length(fig$script, 1)
  expect_match(fig$script, "render_all", fixed = TRUE)
  expect_match(fig$script, "chr2", fixed = TRUE)
  expect_false(grepl(basename(path), fig$script, fixed = TRUE))
  expect_false(grepl(path, fig$script, fixed = TRUE))
})

# ---- the server --------------------------------------------------------------------------
# Custom messages are captured by recording_session() (helper-session.R): a
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
    session$setInputs(pt_rpc = rpc(list(id = 1, op = "geneinfo", args = list(gene = "Eya2")),
                                   list(id = 2, op = "geneinfo", args = list(gene = "bad"))))
    r <- replies(rec)
    expect_identical(r[[1]], list(id = 1L, reply = list(gene = "Eya2", reviewed = TRUE)))
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
    expect_named(fig, c("svg", "script", "key", "job"))
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

# ---- a drawn figure as a PNG or JPEG file ------------------------------------------------

# Width and height in pixels, from a PNG's header and from a JPEG's frame header.
png_size <- function(raw) {
  stopifnot(identical(raw[1:8], as.raw(c(0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a))))
  c(sum(as.integer(raw[17:20]) * 256^(3:0)), sum(as.integer(raw[21:24]) * 256^(3:0)))
}
jpeg_size <- function(raw) {
  stopifnot(identical(raw[1:2], as.raw(c(0xff, 0xd8))))
  i <- 3
  repeat {
    stopifnot(raw[i] == as.raw(0xff))
    two <- function(at) as.integer(raw[at]) * 256 + as.integer(raw[at + 1])
    if (as.integer(raw[i + 1]) %in% c(0xc0, 0xc1, 0xc2)) return(c(two(i + 7), two(i + 5)))
    i <- i + 2 + two(i + 2)
  }
}

test_that("a drawn figure comes back as a PNG and a JPEG of its own size at 300 dpi, as it was drawn", {
  skip_on_cran()
  skip_if_not_installed("org.Mm.eg.db")
  m <- main_srv()
  rec <- recording_session()
  shiny::testServer(.pt_main_app(m$ctx, m$json), {
    session$flushReact()
    session$setInputs(pt_rpc = rpc(list(id = "r1", op = "render", args = render_args())))
    fig <- replies(rec)[[1]]$reply
    expect_identical(fig$job, "j1")
    inches <- as.numeric(regmatches(fig$svg, regexec("<svg[^>]* width='([0-9.]+)pt' height='([0-9.]+)pt'",
                                                     fig$svg))[[1]][2:3]) / 72
    expect_identical(inches[1], 9)
    # The page width changes after the figure was drawn: its files are still the figure as drawn.
    session$setInputs(s_page_w = 5)
    expect_identical(fs$live_cfg()$page_w, 5)
    session$setInputs(pt_rpc = rpc(list(id = "p1", op = "raster", args = list(job = "j1", format = "png")),
                                   list(id = "p2", op = "raster", args = list(job = "j1", format = "jpeg"))))
    r <- replies(rec)
    expect_identical(lapply(r, `[[`, "id"), list("r1", "p1", "p2"))
    expect_named(r[[2]]$reply, c("format", "data"))
    expect_identical(c(r[[2]]$reply$format, r[[3]]$reply$format), c("png", "jpeg"))
    png <- jsonlite::base64_dec(r[[2]]$reply$data)
    jpg <- jsonlite::base64_dec(r[[3]]$reply$data)
    expect_equal(png_size(png), round(inches * 300))
    expect_equal(jpeg_size(jpg), round(inches * 300))
    expect_gt(length(png), 20000)
    expect_gt(length(jpg), 20000)
  }, session = rec$session)
})

test_that("a PNG or JPEG is drawn from the figure's own job, into a file that is removed; what cannot be drawn says why", {
  skip_if_not_installed("org.Mm.eg.db")
  m <- main_srv()
  rec <- recording_session()
  paths <- character(0)
  fail <- FALSE
  # Stands in for the drawing: the file says which job it was drawn from (its title) and as what.
  local_mocked_bindings(.pt_render_svg = function(entries, region_params, settings) {
    paths <<- c(paths, settings$output_path)
    writeLines(paste(settings$title, tools::file_ext(settings$output_path)), settings$output_path)
    if (fail) "the figure failed" else list(hic_tops = NULL)   # a failed draw leaves its file half written
  })
  drawn <- function(reply) rawToChar(jsonlite::base64_dec(reply$data))
  shiny::testServer(.pt_main_app(m$ctx, m$json), {
    session$flushReact()
    session$setInputs(plot_title = "First")
    session$setInputs(pt_rpc = rpc(list(id = 1, op = "render", args = render_args())))
    session$setInputs(plot_title = "Second")
    session$setInputs(pt_rpc = rpc(list(id = 2, op = "render", args = render_args())))
    r <- replies(rec)
    expect_identical(list(r[[1]]$reply$job, r[[2]]$reply$job), list("j1", "j2"))
    session$setInputs(pt_rpc = rpc(
      list(id = 3, op = "raster", args = list(job = "j1", format = "png")),
      list(id = 4, op = "raster", args = list(job = "j2", format = "jpeg")),
      list(id = 5, op = "raster", args = list(job = "j9", format = "png")),
      list(id = 6, op = "raster", args = list(job = "j1", format = "gif")),
      list(id = 7, op = "raster", args = list(job = list("j1", "j2"), format = "jpeg")),
      list(id = 8, op = "raster", args = list(format = "png")),
      list(id = 9, op = "raster", args = list(job = "j1"))))
    r <- replies(rec)
    expect_identical(drawn(r[[3]]$reply), "First png\n")
    expect_identical(drawn(r[[4]]$reply), "Second jpg\n")
    expect_identical(r[[5]]$reply, list(error = "Generate this figure again to save it as PNG."))
    expect_identical(r[[6]]$reply, list(error = "Unknown file format."))
    expect_identical(r[[7]]$reply, list(error = "Generate this figure again to save it as JPEG."))
    expect_identical(r[[8]]$reply, list(error = "Generate this figure again to save it as PNG."))
    expect_identical(r[[9]]$reply, list(error = "Unknown file format."))
    # two renders and two files drawn, each into its own temporary file, none left
    expect_length(unique(paths), 4)
    expect_false(any(file.exists(paths)))
    fail <<- TRUE
    session$setInputs(pt_rpc = rpc(list(id = 10, op = "raster", args = list(job = "j1", format = "png"))))
    expect_identical(replies(rec)[[10]]$reply, list(error = "the figure failed"))
    expect_length(unique(paths), 5)
    expect_false(any(file.exists(paths)))
  }, session = rec$session)
})

test_that("a PNG or JPEG is 300 dpi, or as much as fits in 32767 pixels a side", {
  expect_identical(.pt_raster_dpi(9, 12), 300)
  expect_identical(.pt_raster_dpi(9, 109.2), 300)
  expect_identical(.pt_raster_dpi(9, 110), 297)
  expect_identical(.pt_raster_dpi(400, 9), 81)
  expect_lte(110 * .pt_raster_dpi(9, 110), 32767)
})

test_that("the server keeps the jobs of the last figures only: an older figure says to generate it again", {
  skip_if_not_installed("org.Mm.eg.db")
  m <- main_srv()
  rec <- recording_session()
  local_mocked_bindings(.pt_render_svg = function(entries, region_params, settings) {
    writeLines(paste(settings$title, tools::file_ext(settings$output_path)), settings$output_path)
    list(hic_tops = NULL)
  })
  drawn <- function(reply) rawToChar(jsonlite::base64_dec(reply$data))
  shiny::testServer(.pt_main_server(m$ctx, m$json, jobs_kept = 3L), {
    session$flushReact()
    for (i in 1:3) {
      session$setInputs(plot_title = paste("Figure", i))
      session$setInputs(pt_rpc = rpc(list(id = i, op = "render", args = render_args())))
    }
    session$setInputs(pt_rpc = rpc(list(id = 4, op = "raster", args = list(job = "j1", format = "png"))))
    expect_identical(drawn(replies(rec)[[4]]$reply), "Figure 1 png\n")
    session$setInputs(pt_rpc = rpc(list(id = 5, op = "render", args = render_args()),
                                   list(id = 6, op = "raster", args = list(job = "j1", format = "png")),
                                   list(id = 7, op = "raster", args = list(job = "j2", format = "png")),
                                   list(id = 8, op = "raster", args = list(job = "j4", format = "jpeg"))))
    r <- replies(rec)
    expect_identical(r[[5]]$reply$job, "j4")
    expect_identical(r[[6]]$reply, list(error = "Generate this figure again to save it as PNG."))
    expect_identical(drawn(r[[7]]$reply), "Figure 2 png\n")
    expect_identical(drawn(r[[8]]$reply), "Figure 3 jpg\n")
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

# ---- the render process ------------------------------------------------------------------
# testServer() draws in-process unless it is given a pool. This one hands out a render process
# that only records what it is asked: renders[[i]]$done(result) is the call the real one
# (R/viz-render-worker.R) makes when its drawing ends.
fake_pool <- function() {
  w <- new.env()
  w$renders <- list()
  w$cancels <- 0L
  w$back <- 0L
  worker <- list(
    render = function(entries, region_params, settings, done)
      w$renders[[length(w$renders) + 1L]] <- list(settings = settings, done = done),
    cancel = function() w$cancels <- w$cancels + 1L)
  w$pool <- list(take = function() worker, give_back = function(x) w$back <- w$back + 1L)
  w
}

test_that("with a render process the figure is answered when the process calls back, one drawing at a time", {
  skip_if_not_installed("org.Mm.eg.db")
  m <- main_srv()
  rec <- recording_session()
  w <- fake_pool()
  local_mocked_bindings(.pt_render_svg = function(...) stop("drawn in-process"))
  shiny::testServer(.pt_main_server(m$ctx, m$json, w$pool), {
    session$flushReact()
    session$setInputs(pt_rpc = rpc(list(id = "a", op = "render", args = render_args())))
    expect_length(w$renders, 1)
    expect_length(replies(rec), 0)   # no answer before the process calls back
    # While it draws, another Generate is refused and nothing more is sent to the process.
    session$setInputs(pt_rpc = rpc(list(id = "b", op = "render", args = render_args())))
    expect_identical(replies(rec), list(list(id = "b", reply = list(error = "A figure is already being drawn."))))
    expect_length(w$renders, 1)
    job <- w$renders[[1]]
    writeLines("<svg>drawn by the process</svg>", job$settings$output_path)
    job$done(list(hic_tops = NULL))
    r <- replies(rec)
    expect_identical(r[[2]]$id, "a")
    expect_identical(r[[2]]$reply$svg, "<svg>drawn by the process</svg>")
    expect_identical(r[[2]]$reply$job, "j1")
    expect_match(r[[2]]$reply$script, "Eya2", fixed = TRUE)
    expect_false(file.exists(job$settings$output_path))
    # A PNG of it is drawn by the process too, and holds off a Generate the same way.
    session$setInputs(pt_rpc = rpc(list(id = "p", op = "raster", args = list(job = "j1", format = "png")),
                                   list(id = "c", op = "render", args = render_args())))
    expect_length(w$renders, 2)
    expect_identical(replies(rec)[[3]], list(id = "c", reply = list(error = "A figure is already being drawn.")))
    png <- w$renders[[2]]
    expect_identical(tools::file_ext(png$settings$output_path), "png")
    writeBin(as.raw(1:9), png$settings$output_path)
    png$done(list(hic_tops = NULL))
    r <- replies(rec)
    expect_identical(r[[4]]$id, "p")
    expect_identical(jsonlite::base64_dec(r[[4]]$reply$data), as.raw(1:9))
    expect_false(file.exists(png$settings$output_path))
    # What the process reports as an error is the answer, and the next Generate is taken.
    session$setInputs(pt_rpc = rpc(list(id = "d", op = "render", args = render_args())))
    expect_length(w$renders, 3)
    w$renders[[3]]$done("The render process stopped.")
    expect_identical(replies(rec)[[5]], list(id = "d", reply = list(error = "The render process stopped.")))
    session$setInputs(pt_rpc = rpc(list(id = "e", op = "render", args = render_args())))
    expect_length(w$renders, 4)
    expect_identical(w$cancels, 0L)
  }, session = rec$session)
  expect_identical(w$back, 1L)   # the session's end gives the process back to the pool
})

test_that("Cancel stops the render process and answers the drawing as cancelled, a figure or a PNG", {
  skip_if_not_installed("org.Mm.eg.db")
  m <- main_srv()
  rec <- recording_session()
  w <- fake_pool()
  shiny::testServer(.pt_main_server(m$ctx, m$json, w$pool), {
    session$flushReact()
    cancelled <- list(error = "Cancelled.", cancelled = TRUE)
    # Nothing drawing: the process is left alone.
    session$setInputs(pt_rpc = rpc(list(id = 0, op = "cancel")))
    expect_identical(w$cancels, 0L)
    session$setInputs(pt_rpc = rpc(list(id = "a", op = "render", args = render_args())))
    session$setInputs(pt_rpc = rpc(list(id = 1, op = "cancel")))
    expect_identical(w$cancels, 1L)
    r <- replies(rec)
    expect_identical(r[[2]], list(id = "a", reply = cancelled))
    expect_identical(r[[3]], list(id = 1L, reply = list(ok = TRUE)))
    # The next Generate is taken at once; its figure, then a PNG of it cancelled.
    session$setInputs(pt_rpc = rpc(list(id = "b", op = "render", args = render_args())))
    expect_length(w$renders, 2)
    writeLines("<svg/>", w$renders[[2]]$settings$output_path)
    w$renders[[2]]$done(list(hic_tops = NULL))
    job <- replies(rec)[[4]]$reply$job
    session$setInputs(pt_rpc = rpc(list(id = "p", op = "raster", args = list(job = job, format = "jpeg"))))
    session$setInputs(pt_rpc = rpc(list(id = 2, op = "cancel")))
    expect_identical(w$cancels, 2L)
    expect_identical(replies(rec)[[5]], list(id = "p", reply = cancelled))
    # The figure is still remembered: its PNG can be asked for again.
    session$setInputs(pt_rpc = rpc(list(id = "q", op = "raster", args = list(job = job, format = "png"))))
    expect_length(w$renders, 4)
  }, session = rec$session)
})

test_that("the script is told the Hi-C colour tops the drawing reports; a script that cannot be written is the answer's error", {
  skip_if_not_installed("org.Mm.eg.db")
  m <- main_srv()
  rec <- recording_session()
  tops <- list()
  paths <- character(0)
  fail <- FALSE
  local_mocked_bindings(
    .pt_render_svg = function(entries, region_params, settings) {
      paths <<- c(paths, settings$output_path)
      writeLines("<svg/>", settings$output_path)
      list(hic_tops = c(hic_zrange_HiC = 41.5))
    },
    assemble_download_script = function(entries, region_params, settings, hic_tops = NULL) {
      tops[[length(tops) + 1L]] <<- hic_tops
      if (fail) stop("no script for you")
      "the script"
    })
  shiny::testServer(.pt_main_app(m$ctx, m$json), {
    session$flushReact()
    session$setInputs(pt_rpc = rpc(list(id = 1, op = "render", args = render_args())))
    expect_identical(replies(rec)[[1]]$reply$script, "the script")
    fail <<- TRUE
    session$setInputs(pt_rpc = rpc(list(id = 2, op = "render", args = render_args())))
    expect_identical(replies(rec)[[2]]$reply, list(error = "no script for you"))
    # Not remembered as a figure, and its file is gone.
    session$setInputs(pt_rpc = rpc(list(id = 3, op = "raster", args = list(job = "j2", format = "png"))))
    expect_match(replies(rec)[[3]]$reply$error, "Generate this figure again", fixed = TRUE)
  }, session = rec$session)
  expect_identical(tops, list(c(hic_zrange_HiC = 41.5), c(hic_zrange_HiC = 41.5)))
  expect_false(any(file.exists(paths)))
})
