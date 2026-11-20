# The main screen's app and server. The page (R/app-screen-ui.R, inst/www/main/)
# keeps its own state; this side answers its questions and draws its figures.
# Questions arrive as one input, pt_rpc: a JSON array of {id, op, args}. Each
# is answered with a "pt_rpc_reply" message {id, json}. The data answers are
# R/app-screen-api.R; "render" draws through the session's render worker
# (.pt_build_entries() -> .pt_render_svg()), with the Figure settings panel's
# live settings.

# An answer as the JSON text the page parses. Written here, not left to
# Shiny's own message encoding, so one number format holds for every answer
# (digits = NA: coordinates in full).
.pt_main_json <- function(x)
  as.character(jsonlite::toJSON(x, auto_unbox = TRUE, null = "null",
                                digits = NA))

# What makes two figures of the same choice differ: every display setting and
# the title, as one string. The page keeps it with each drawn figure, so a
# settings change empties the canvas and going back to the earlier settings
# brings the earlier figure back.
.pt_main_settings_key <- function(cfg, display_keys, title) {
  .pt_main_json(list(settings = unname(lapply(cfg[display_keys], function(v)
    if (is.null(v)) NA else v)), title = trimws(title %||% "")))
}

# A render request from the page ({gene, chr, s, e, picked, tracks}) as what
# the render path takes: entries, region_params and settings. Stops with a
# message for the page when the request cannot be drawn. `cfg` is the live
# render settings, `title` the plot title.
.pt_main_render_job <- function(ctx, args, cfg, title = "") {
  picked <- as.character(unlist(args$picked))
  tracks <- args$tracks
  if (!length(tracks)) stop("Add at least one track.", call. = FALSE)
  labels <- .pt_main_labels(ctx, picked, tracks)

  chrom <- as.character(unlist(args$chr %||% ""))
  len <- if (length(chrom) == 1L) .pt_main_chr_length(ctx, chrom) else NA
  if (is.na(len))
    stop("No chromosome '", paste(chrom, collapse = ""), "'.", call. = FALSE)
  s <- suppressWarnings(as.numeric(args$s %||% NA))
  e <- suppressWarnings(as.numeric(args$e %||% NA))
  if (length(s) != 1L || length(e) != 1L || is.na(s) || is.na(e))
    stop("The region needs a start and an end.", call. = FALSE)
  err <- .pt_region_error(s, e)
  if (!is.null(err)) stop(err, call. = FALSE)
  # Clamped to the chromosome, as the page's zoom and move can overshoot its
  # ends by a little; a region with nothing left inside it is an error.
  s <- max(1, round(s))
  e <- min(len, round(e))
  if (s >= e)
    stop("That region is outside ", chrom, " (1 to ",
         format(len, big.mark = ",", scientific = FALSE), ").", call. = FALSE)

  # A region has no gene: NA.
  # A gene is taken only as one this project's annotation has, under its
  # symbol (.pt_main_gene_row()): the page's text is never trusted as a name
  # (a track module may build a cache file's path from it).
  gene <- args$gene
  gene <- if (is.character(gene) && length(gene) == 1L && !is.na(gene) &&
              nzchar(gene)) {
    row <- tryCatch(.pt_main_gene_row(ctx, gene),
                    error = function(e) list(error = conditionMessage(e)))
    if (is.null(row$error)) row$symbol else NA_character_
  } else NA_character_

  # Never generated: a track whose file is missing, and a gene-only track for a
  # region.
  bad <- Filter(function(k) {
    it <- ctx$TRACK_ITEMS[[k]]
    isTRUE(it$broken) ||
      (is.na(gene) && it$type %in%
         (ctx$GENE_REQUIRED_TYPES %||% .pt_gene_required_types()))
  }, unique(labels))
  if (length(bad))
    stop("Remove before generating: ",
         paste(vapply(bad, function(k) ctx$TRACK_ITEMS[[k]]$label %||% k,
                      character(1)), collapse = ", "),
         if (is.na(gene)) " (needs a gene, or its file is missing)."
         else " (its file is missing).", call. = FALSE)

  list(
    entries = .pt_build_entries(labels, ctx$TRACK_ITEMS,
                                ctx$DEFAULT_SAMPLE_COLS, cfg),
    region_params = list(chrom = chrom, start = as.integer(s),
                         end = as.integer(e), gene = gene),
    settings = c(.pt_settings_base(cfg),
                 list(title = trimws(title %||% ""),
                      output_path = tempfile(fileext = ".svg")))
  )
}

# The drawn figure as the page takes it: the SVG's text and the R script that
# reproduces it. The SVG file is the render's temporary one and is removed.
# `hic_tops`: the Hi-C colour tops the render found (.pt_render_svg()'s result),
# which the script names.
.pt_main_figure <- function(job, hic_tops = NULL) {
  path <- job$settings$output_path
  on.exit(unlink(path))
  list(
    svg = paste(readLines(path, warn = FALSE, encoding = "UTF-8"),
                collapse = "\n"),
    script = paste(assemble_download_script(job$entries, job$region_params,
                                            job$settings, hic_tops = hic_tops),
                   collapse = "\n")
  )
}

# Stopping the app (it has no Close app button): the app stops, with the value
# "closed", shortly after the last open page is closed by its user. A page that was only cut off (a dropped
# connection, a sleeping laptop) does not count, so the app is still there to
# reload; nor does a reload, whose new session arrives within `grace` seconds.
# `stop` and `later` are arguments so a test can stand in for them.
.pt_main_closer <- function(grace = 8,
                            stop = function() stopApp(returnValue = "closed"),
                            later = later::later) {
  n <- 0L
  over <- FALSE
  list(
    opened = function() n <<- n + 1L,
    ended = function(by_user) {
      n <<- n - 1L
      if (isTRUE(by_user))
        later(function() if (!over && n <= 0L) {
          over <<- TRUE
          stop()
        }, grace)
      invisible()
    },
    # the app has stopped for another reason (Edit project settings, Reset):
    # a timer still pending must not stop the app that follows it
    stopped = function() over <<- TRUE
  )
}

.pt_main_server <- function(ctx, SETTINGS_JSON, RENDER_POOL = NULL,
                            closer = NULL, jobs_kept = 500L) {
  function(input, output, session) {
    # The page says when its user closes or reloads it: the pt_bye input, sent
    # over the connection just before it closes, and a request to this URL as
    # well. The session's end alone cannot tell a closed page from a lost
    # connection.
    leaving <- FALSE
    bye_url <- tryCatch(session$registerDataObj("bye", NULL, function(data, req) {
      leaving <<- TRUE
      httpResponse(204L, "text/plain", "")
    }), error = function(e) NULL)
    observeEvent(input$pt_bye, leaving <<- TRUE)
    if (!is.null(closer)) {
      closer$opened()
      session$onSessionEnded(function() closer$ended(leaving))
    }

    # repeat-click gate for the "pt-once" buttons (R/app-once-gate.R)
    .pt_once_server(input, session)

    # This session's background render process (R/viz-render-worker.R), from
    # RENDER_POOL (the launcher's, .pt_render_pool()) when there is one. Without
    # a pool shiny::testServer() has none and draws in-process.
    render_worker <- if (!is.null(RENDER_POOL)) RENDER_POOL$take()
      else if (!inherits(session, "MockShinySession")) .pt_render_worker()
    session$onSessionEnded(function() {
      if (is.null(render_worker)) return()
      if (is.null(RENDER_POOL)) render_worker$close()
      else RENDER_POOL$give_back(render_worker)
    })

    # Figure settings: the panel, saved and live settings, Save and Revert
    # (R/app-settings.R). Edit project settings and its dialogs:
    # R/app-project-actions.R.
    fs <- .pt_figure_settings_server(input, output, ctx$RENDER_CFG,
                                     SETTINGS_JSON)
    unsaved_settings_note <- function()
      .pt_unsaved_settings_note(fs$settings_dirty())
    .pt_project_actions_server(input, session, unsaved_settings_note,
                               ctx$PROJECT_DIR, SETTINGS_JSON, ctx$CACHE_DIR)

    reply <- function(id, x)
      session$sendCustomMessage("pt_rpc_reply",
                                list(id = id, json = .pt_main_json(x)))

    # The page learns when the settings change (and what they are, as a key).
    observe({
      session$sendCustomMessage("pt_settings", list(
        key = .pt_main_settings_key(fs$live_cfg(), fs$display_keys,
                                    input$plot_title)))
    })

    # The render in flight (its request id), or NULL. One at a time: the page
    # asks for the next only when the last one has answered.
    drawing <- NULL
    # The figures drawn in this session, as the jobs they were drawn from, by
    # the name the page was given for each (the reply's `job`): what a PNG or
    # JPEG of a figure is drawn from, so it has the settings the figure was
    # drawn with. A job is a small list; the last `jobs_kept` are kept, far more
    # than the page holds (40 figures).
    jobs <- list()
    n_jobs <- 0L
    remember <- function(job) {
      n_jobs <<- n_jobs + 1L
      jobs[[paste0("j", n_jobs)]] <<- job
      if (length(jobs) > jobs_kept) jobs <<- utils::tail(jobs, jobs_kept)
      paste0("j", n_jobs)
    }
    start_render <- function(id, args) {
      if (!is.null(drawing))
        return(reply(id, list(error = "A figure is already being drawn.")))
      job <- tryCatch(
        .pt_main_render_job(ctx, args, fs$live_cfg(), input$plot_title),
        error = function(e) conditionMessage(e))
      if (is.character(job)) return(reply(id, list(error = job)))
      # The settings this figure is drawn with, as the key the page files it
      # under: a setting changed just before Generate may not have reached the
      # page as a "pt_settings" message yet.
      key <- .pt_main_settings_key(fs$live_cfg(), fs$display_keys,
                                   input$plot_title)
      drawing <<- id
      # Called from the app's event loop when the background render ends (or
      # at once, in-process), with .pt_render_svg()'s result: a list, or the
      # error message. A cancelled render never gets here: cancel() below has
      # answered it, and the worker drops its callback.
      done <- function(res) {
        drawing <<- NULL
        out <- if (is.character(res)) {
          unlink(job$settings$output_path)
          list(error = res)
        } else tryCatch(c(.pt_main_figure(job, res$hic_tops),
                          list(key = key, job = remember(job))),
                        error = function(e) list(error = conditionMessage(e)))
        reply(id, out)
      }
      if (is.null(render_worker))
        done(.pt_render_svg(job$entries, job$region_params, job$settings))
      else render_worker$render(job$entries, job$region_params, job$settings,
                                done)
    }
    # A drawn figure again, as a PNG or JPEG file (args: the figure's `job`
    # name and the format): drawn from its remembered job by the render path,
    # answered as the file's bytes in base64. It takes the place of a render
    # while it draws: one at a time, and a "cancel" request stops it too.
    start_raster <- function(id, args) {
      format <- args$format
      if (!is.character(format) || length(format) != 1L ||
          !format %in% c("png", "jpeg"))
        return(reply(id, list(error = "Unknown file format.")))
      name <- args$job
      job <- if (is.character(name) && length(name) == 1L && !is.na(name))
        jobs[name][[1]]
      if (is.null(job))
        return(reply(id, list(error = paste0(
          "Generate this figure again to save it as ", toupper(format), "."))))
      if (!is.null(drawing))
        return(reply(id, list(error = "A figure is being drawn: try again when it is done.")))
      path <- tempfile(fileext = if (format == "png") ".png" else ".jpg")
      job$settings$output_path <- path
      drawing <<- id
      done <- function(res) {
        drawing <<- NULL
        on.exit(unlink(path))
        reply(id, if (is.character(res)) list(error = res) else tryCatch(
          list(format = format, data = jsonlite::base64_enc(
            readBin(path, "raw", file.size(path)))),
          error = function(e) list(error = conditionMessage(e))))
      }
      if (is.null(render_worker))
        done(.pt_render_svg(job$entries, job$region_params, job$settings))
      else render_worker$render(job$entries, job$region_params, job$settings,
                                done)
    }
    cancel_render <- function() {
      if (is.null(drawing)) return()
      id <- drawing
      drawing <<- NULL
      if (!is.null(render_worker)) render_worker$cancel()
      reply(id, list(error = "Cancelled.", cancelled = TRUE))
    }

    observeEvent(input$pt_rpc, {
      batch <- tryCatch(jsonlite::fromJSON(input$pt_rpc, simplifyVector = FALSE),
                        error = function(e) list())
      # The page sends an array of requests. Anything else must not end the
      # session (an error in an observer does): one request sent bare is taken
      # as an array of one, and what is not a request is skipped.
      if (!is.list(batch)) batch <- list()
      if (!is.null(names(batch))) batch <- list(batch)
      for (rq in batch) {
        if (!is.list(rq) || is.null(rq$id)) next
        id <- rq$id
        op <- rq$op %||% ""
        args <- rq$args %||% list()
        if (identical(op, "render")) start_render(id, args)
        else if (identical(op, "raster")) start_raster(id, args)
        else if (identical(op, "cancel")) {
          cancel_render()
          reply(id, list(ok = TRUE))
        } else {
          # the gene flank (locus) and the lookups (geneinfo) are the live
          # settings'
          res <- .pt_main_request(ctx, op, args, fs$live_cfg())
          if (identical(op, "project") && is.null(res$error)) {
            res$id <- basename(ctx$PROJECT_DIR)
            res$version <- as.character(utils::packageVersion("PlotTrellis"))
            res$bye <- bye_url
          }
          reply(id, res)
        }
      }
    })
  }
}

# The app with the main screen, from the startup context
# (.pt_composer_context()).
# options(PlotTrellis.stop_on_close = FALSE) keeps the app running after its
# last page is closed (browser-driven checks open and close many pages).
.pt_main_app <- function(ctx, SETTINGS_JSON, RENDER_POOL = NULL) {
  # The gene index behind the search box and the nearby genes takes seconds on
  # a whole genome: built here, with the rest of the launch, not at the first
  # keystroke (it is kept for the R session, .pt_main_gene_index()).
  .pt_main_gene_index(ctx)
  # A suggested gene this annotation does not have is left out: said once here,
  # at launch (immediate: runApp() would hold it until the app closes).
  unknown <- .pt_main_suggested_rows(ctx)$unknown
  if (length(unknown))
    warning("Not suggested, as this project's annotation has no such gene: ",
            paste(unknown, collapse = ", "), " (# suggested_genes:).",
            call. = FALSE, immediate. = TRUE)
  closer <- if (isTRUE(getOption("PlotTrellis.stop_on_close", TRUE)))
    .pt_main_closer()
  shinyApp(ui = .pt_main_ui(ctx),
           server = .pt_main_server(ctx, SETTINGS_JSON, RENDER_POOL, closer),
           onStart = if (!is.null(closer)) function() onStop(closer$stopped))
}
