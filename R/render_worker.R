# The composer's figure renders (Generate, zoom, Apply to figure) run in a
# background R process, one per browser session, kept alive between renders so
# its loaded packages and opened annotation stay warm. The app starts one as it
# launches (loading the packages takes seconds) and a session that ends hands
# its process to the next one (a page reload). The app's own R process stays
# free to react while a figure draws, so the page's Cancel button can stop a
# long render: it kills the background process and starts a fresh one, and the
# previous figure stays on the page. The first render after a Cancel waits for
# the fresh process to load its packages.

# Draws one figure to settings$output_path as a tidied SVG file. Returns NULL,
# or the error message when the whole figure failed (a single failing track is
# drawn as its own in-figure message by render_all(), not reported here). Run in
# the render process.
.pt_render_svg <- function(entries, region_params, settings) {
  dev_before <- dev.cur()
  tryCatch({
    # .pt_quiet(): plotgardener prints a "text[text2]"-style line per grob,
    # hundreds per figure.
    .pt_quiet(render_all(entries, region_params, settings,
                         device_fn = function(w, h) svglite::svglite(
                           settings$output_path, width = w, height = h
                         )))
    if (dev.cur() != dev_before) dev.off()
    .pt_tidy_svg(settings$output_path, unlock_text = FALSE)
    NULL
  }, error = function(e) {
    if (dev.cur() != dev_before) try(dev.off(), silent = TRUE)
    conditionMessage(e)
  })
}

# How the render process loads this package: from the same source tree when this
# session loaded it from source (pkgload::load_all(), as the tests do), else the
# installed copy.
.pt_worker_source <- function() {
  path <- getNamespaceInfo(asNamespace("PlotTrellis"), "path")
  if (length(list.files(file.path(path, "R"), pattern = "\\.R$")) > 0 &&
      requireNamespace("pkgload", quietly = TRUE)) path else ""
}

# Run in the render process before each render: loads this package once. It is
# sent to the process detached from this package (environment set to the base
# one), since calling it through this package's namespace would load the
# installed copy before the source tree. pkgload is a Suggests dependency: `src`
# is non-empty only when .pt_worker_source() found it installed, so an installed
# PlotTrellis never needs it.
.pt_worker_load <- function(src) {
  if (!isNamespaceLoaded("PlotTrellis")) {
    if (nzchar(src) && requireNamespace("pkgload", quietly = TRUE))
      suppressMessages(pkgload::load_all(src, quiet = TRUE, helpers = FALSE,
                                         export_all = FALSE))
    else loadNamespace("PlotTrellis")
  }
  invisible(NULL)
}

# The render process for one browser session. render() hands a figure to it and
# returns at once; done(err) is called later, from the app's event loop, with
# NULL or an error message. cancel() kills a running render (done() is then
# never called) and starts a fresh process. close() ends the process (the
# session ended). What a render prints (a failed track's line, notes, warnings)
# is shown on this process's console when it ends, or when it is cancelled, as
# far as it got. draw is the function run in the process for each render (NULL:
# .pt_render_svg()); the tests pass a slow one. If the process cannot be started
# or cannot load the package, figures are drawn in the app's own process
# instead, as before there was a render process (no Cancel then).
.pt_render_worker <- function(poll_secs = 0.1, draw = NULL) {
  src     <- .pt_worker_source()
  load    <- .pt_worker_load
  environment(load) <- baseenv()
  rs      <- NULL
  # NULL (idle), "start" (process starting), "load", or "render"
  job     <- NULL
  # a render asked for before the process was ready
  pending <- NULL
  done    <- NULL
  # the running render's SVG file, removed if it is cancelled
  out     <- NULL
  # the file the running render's printed output and messages go to
  log     <- NULL
  closed  <- FALSE
  # the process could not start or load the package: draw in-process
  broken  <- FALSE
  # a poll is scheduled; one at a time
  polling <- FALSE

  schedule <- function() {
    if (polling || closed) return(invisible())
    polling <<- TRUE
    later::later(function() {
      polling <<- FALSE
      # An error here would reach the app's event loop and stop the app: report
      # it as the render's error instead, and start over with a fresh process.
      tryCatch(poll(), error = function(e) {
        message("PlotTrellis: the render process handler failed: ",
                conditionMessage(e))
        stop_process()
        start()
        finish(paste("The figure could not be drawn:", conditionMessage(e)))
      })
    }, poll_secs)
  }

  # Not waiting long for the process: the app carries on while it starts.
  # Without the user's .Rprofile (a project's renv activation there can take
  # seconds); the library paths are passed on as this process has them.
  start <- function() {
    if (closed || broken) return(invisible())
    rs <<- tryCatch(callr::r_session$new(
      options = callr::r_session_options(user_profile = FALSE), wait = FALSE
    ), error = function(e) NULL)
    if (is.null(rs)) return(give_up())
    job <<- "start"
    # Usually up in well under a second: then the package load starts now, not
    # at the first poll (which, at app launch, only runs once the app is
    # serving).
    if (rs$poll_process(2000) == "ready")
      tryCatch(poll(), error = function(e) give_up())
    schedule()
  }

  stop_process <- function() {
    if (!is.null(rs)) {
      try(rs$kill(), silent = TRUE)
      try(rs$close(), silent = TRUE)
    }
    rs <<- NULL
    job <<- NULL
  }

  in_process <- function(r) (if (is.null(draw)) .pt_render_svg else draw)(
    r$entries, r$region_params, r$settings
  )

  # The process could not be started or set up: draw in this process from now
  # on.
  give_up <- function() {
    stop_process()
    broken <<- TRUE
    if (!is.null(pending)) {
      r <- pending
      pending <<- NULL
      r$done(in_process(r))
    }
  }

  dispatch <- function(r) {
    done <<- r$done
    out  <<- r$settings$output_path
    log  <<- tempfile("pt-render-", fileext = ".log")
    rs$call(function(load, src, draw, log, entries, region_params, settings) {
      load(src)
      # Looked up by name in the process's loaded namespace (the source tree's
      # or the installed one).
      if (is.null(draw))
        draw <- get(".pt_render_svg", envir = asNamespace("PlotTrellis"))
      # Printed output, messages and warnings, in order, to `log`; written out
      # at each message or warning, so a killed render leaves what it had said
      # so far.
      con <- file(log, open = "w")
      sink(con)
      on.exit({ sink(); close(con) })
      withCallingHandlers(
        draw(entries, region_params, settings),
        message = function(m) {
          cat(conditionMessage(m), file = con)
          flush(con)
          invokeRestart("muffleMessage")
        },
        warning = function(w) {
          cat("Warning: ", conditionMessage(w), "\n", sep = "", file = con)
          flush(con)
          invokeRestart("muffleWarning")
        }
      )
    }, list(load, src, draw, log, r$entries, r$region_params, r$settings))
    job <<- "render"
  }

  # Shows on this process's console what the render process printed for the last
  # render.
  relay <- function() {
    if (is.null(log)) return(invisible())
    # No file when the render was cancelled before it began.
    txt <- if (file.exists(log))
      tryCatch(
        readChar(log, file.size(log), useBytes = TRUE),
        error = function(e) ""
      ) else ""
    unlink(log)
    log <<- NULL
    if (length(txt) == 1 && nzchar(txt))
      message(txt, appendLF = !endsWith(txt, "\n"))
  }

  finish <- function(err) {
    relay()
    cb <- done
    done <<- NULL
    # the file is the caller's now: a later Cancel must not remove it
    out <<- NULL
    if (!is.null(cb)) cb(err)
  }

  poll <- function() {
    if (closed || is.null(rs) || is.null(job)) return(invisible())
    if (rs$poll_process(0) == "ready") {
      res <- tryCatch(rs$read(),
                      error = function(e) list(code = 501L, error = e))
      if (!is.null(res)) {
        # code 501: the process died (out of memory, a crash).
        crashed <- isTRUE(res$code == 501) || !rs$is_alive()
        err <- if (!is.null(res$error))
          paste("The figure could not be drawn:",
                conditionMessage(res$error$parent %||% res$error))
        else if (crashed)
          "The figure could not be drawn: the render process stopped."
        was <- job
        job <<- NULL
        if (was %in% c("start", "load") && (crashed || !is.null(err)))
          return(give_up())
        if (identical(was, "start")) {
          rs$call(load, list(src))
          job <<- "load"
        } else if (identical(was, "load")) {
          if (!is.null(pending)) { r <- pending; pending <<- NULL; dispatch(r) }
        } else {
          if (crashed) { stop_process(); start() }
          finish(if (is.null(err)) res$result else err)
          # A render asked for while this one ran (the app asks for one at a
          # time).
          if (!is.null(pending) && is.null(job) && !is.null(rs)) {
            r <- pending; pending <<- NULL
            dispatch(r)
          }
        }
      }
    }
    if (!is.null(job)) schedule()
    invisible()
  }

  start()
  list(
    render = function(entries, region_params, settings, done) {
      r <- list(entries = entries, region_params = region_params,
                settings = settings, done = done)
      if (broken) return(done(in_process(r)))
      # it died while idle
      if (is.null(rs) || !rs$is_alive()) { stop_process(); start() }
      if (broken) return(done(in_process(r)))
      if (!is.null(job))
        # starting, loading, or drawing another figure
        pending <<- r
      else { dispatch(r); schedule() }
    },
    cancel = function() {
      pending <<- NULL
      done <<- NULL
      stop_process()
      relay()
      if (!is.null(out)) unlink(out)
      start()
    },
    rendering = function() identical(job, "render") || !is.null(pending),
    close = function() {
      closed <<- TRUE
      pending <<- NULL
      done <<- NULL
      stop_process()
      if (!is.null(out)) unlink(out)
      if (!is.null(log)) unlink(log)
    }
  )
}

# The app's render processes: take() hands a session the process started at
# launch, or the one an earlier session gave back, else a new one; give_back()
# keeps a session's process for the next session unless a render is still
# running in it (at most one is kept).
.pt_render_pool <- function() {
  kept <- .pt_render_worker()
  # a session can end after the app has stopped: its process is then closed
  closed <- FALSE
  list(
    take = function() {
      w <- kept
      kept <<- NULL
      if (is.null(w)) .pt_render_worker() else w
    },
    give_back = function(w) {
      if (!closed && is.null(kept) && !w$rendering()) kept <<- w else w$close()
    },
    close = function() {
      closed <<- TRUE
      if (!is.null(kept)) kept$close()
      kept <<- NULL
    }
  )
}
