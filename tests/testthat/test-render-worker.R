# The composer draws its figures in a background R process (.pt_render_worker()), so a
# render can be cancelled. These check that the process draws the same figure the app
# drew in its own process, and that Cancel stops a slow render quickly and leaves a
# working process behind.

# Runs the app's event loop until done() has been called or `secs` have passed.
wait_for <- function(got, secs) {
  t0 <- Sys.time()
  while (!got() && as.numeric(Sys.time() - t0, units = "secs") < secs) later::run_now(0.1)
  got()
}

test_that("a figure drawn in the render process is the one drawn in-process (mouse dataset)", {
  skip_on_cran()
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  rows <- test_rows("mouse")
  specs <- list(single_track("genes", rows[["Ba"]]), single_track("transcripts", rows[["Ba"]]),
                single_track("bw", rows[["Ba"]]), single_track("hic", rows[["Ba"]]))
  here  <- c(test_settings("mouse"), list(output_path = tempfile(fileext = ".svg")))
  there <- c(test_settings("mouse"), list(output_path = tempfile(fileext = ".svg")))
  on.exit(unlink(c(here$output_path, there$output_path)), add = TRUE)
  expect_null(.pt_render_svg(specs, test_region("mouse"), here))

  w <- .pt_render_worker()
  on.exit(w$close(), add = TRUE)
  res <- list()
  w$render(specs, test_region("mouse"), there, function(err) res$err <<- list(err))
  expect_true(wait_for(function() length(res) > 0, 120))
  expect_null(res$err[[1]])
  expect_identical(readLines(there$output_path), readLines(here$output_path))
})

test_that("Cancel stops a slow render within seconds and the next render still runs", {
  skip_on_cran()
  draw <- function(entries, region_params, settings) {
    if (isTRUE(settings$slow)) { writeLines("started", settings$started); Sys.sleep(300) }
    writeLines("drawn", settings$output_path)
    NULL
  }
  w <- .pt_render_worker(draw = draw)
  on.exit(w$close(), add = TRUE)
  slow <- list(slow = TRUE, output_path = tempfile(fileext = ".svg"), started = tempfile())
  fast <- list(output_path = tempfile(fileext = ".svg"))
  on.exit(unlink(c(slow$output_path, slow$started, fast$output_path)), add = TRUE)

  called <- FALSE
  w$render(list(), list(), slow, function(err) called <<- TRUE)
  # Cancelled while it draws, not while the process is still loading the package.
  expect_true(wait_for(function() file.exists(slow$started), 180))
  t0 <- Sys.time()
  w$cancel()
  expect_lt(as.numeric(Sys.time() - t0, units = "secs"), 5)
  wait_for(function() called, 5)
  expect_false(called)  # a cancelled render never reports back
  expect_false(file.exists(slow$output_path))

  res <- list()
  w$render(list(), list(), fast, function(err) res$err <<- list(err))
  expect_true(wait_for(function() length(res) > 0, 120))
  expect_null(res$err[[1]])
  expect_identical(readLines(fast$output_path), "drawn")
})

test_that("a render process that dies reports an error and is replaced", {
  skip_on_cran()
  draw <- function(entries, region_params, settings) {
    if (isTRUE(settings$crash)) quit(save = "no", status = 3)
    writeLines("drawn", settings$output_path)
    NULL
  }
  w <- .pt_render_worker(draw = draw)
  on.exit(w$close(), add = TRUE)
  res <- list()
  w$render(list(), list(), list(crash = TRUE, output_path = tempfile()), function(err) res$a <<- list(err))
  expect_true(wait_for(function() length(res) > 0, 120))
  expect_match(res$a[[1]], "could not be drawn")
  out <- tempfile()
  on.exit(unlink(out), add = TRUE)
  w$render(list(), list(), list(output_path = out), function(err) res$b <<- list(err))
  expect_true(wait_for(function() !is.null(res$b), 120))
  expect_null(res$b[[1]])
  expect_identical(readLines(out), "drawn")
})

test_that("a render process that cannot start leaves figures drawn in the app's own process", {
  skip_on_cran()
  local_mocked_bindings(r_session_options = function(...) stop("no R here (test)"), .package = "callr")
  w <- .pt_render_worker(draw = function(entries, region_params, settings) { writeLines(as.character(Sys.getpid()), settings$output_path); NULL })
  on.exit(w$close(), add = TRUE)
  out <- tempfile()
  on.exit(unlink(out), add = TRUE)
  res <- list()
  w$render(list(), list(), list(output_path = out), function(err) res$a <<- list(err))
  expect_length(res, 1)  # done() was called straight away, before render() returned
  expect_null(res$a[[1]])
  expect_identical(readLines(out), as.character(Sys.getpid()))
})

test_that("the app's pool hands a session its kept process and keeps one given back", {
  skip_on_cran()
  pool <- .pt_render_pool()
  on.exit(pool$close(), add = TRUE)
  a <- pool$take()
  b <- pool$take()  # none kept now: a new one
  expect_false(identical(a, b))
  pool$give_back(a)
  expect_identical(pool$take(), a)  # a page reload gets the same, already loaded, process
  pool$give_back(a)
  pool$give_back(b)  # one is kept already: this one is closed
  out <- tempfile()
  on.exit(unlink(out), add = TRUE)
  res <- list()
  w <- pool$take()
  expect_identical(w, a)
  w$render(list(), list(), list(output_path = out), function(err) res$a <<- list(err))
  expect_true(wait_for(function() length(res) > 0, 120))
})

test_that("a render process that exits cleanly mid-render reports an error, not a figure", {
  skip_on_cran()
  w <- .pt_render_worker(draw = function(entries, region_params, settings) quit(save = "no", status = 0))
  on.exit(w$close(), add = TRUE)
  res <- list()
  w$render(list(), list(), list(output_path = tempfile()), function(err) res$a <<- list(err))
  expect_true(wait_for(function() length(res) > 0, 120))
  expect_match(res$a[[1]], "render process stopped")
})

# Collects the messages shown on this process's console while f() runs.
console_of <- function(f) {
  msgs <- character(0)
  withCallingHandlers(f(), message = function(m) {
    msgs <<- c(msgs, conditionMessage(m))
    invokeRestart("muffleMessage")
  })
  paste(msgs, collapse = "")
}

test_that("what a render prints reaches the app's console, also from a cancelled render", {
  skip_on_cran()
  draw <- function(entries, region_params, settings) {
    cat("printed line\n")
    message("PlotTrellis: hic (A) track failed: no file")
    warning("could not write the cache file")
    if (isTRUE(settings$slow)) { writeLines("started", settings$started); Sys.sleep(300) }
    writeLines("drawn", settings$output_path)
    NULL
  }
  w <- .pt_render_worker(draw = draw)
  on.exit(w$close(), add = TRUE)
  fast <- list(output_path = tempfile(fileext = ".svg"))
  slow <- list(slow = TRUE, output_path = tempfile(fileext = ".svg"), started = tempfile())
  on.exit(unlink(c(fast$output_path, slow$output_path, slow$started)), add = TRUE)

  res <- list()
  w$render(list(), list(), fast, function(err) res$err <<- list(err))
  shown <- console_of(function() expect_true(wait_for(function() length(res) > 0, 120)))
  expect_null(res$err[[1]])
  expect_identical(shown, "printed line\nPlotTrellis: hic (A) track failed: no file\nWarning: could not write the cache file\n")

  # Cancelled at once, before the idle process has begun the render: nothing to show, and no
  # warning about a log file not yet written.
  w$render(list(), list(), slow, function(err) NULL)
  expect_no_warning(expect_identical(console_of(w$cancel), ""))

  w$render(list(), list(), slow, function(err) NULL)
  expect_true(wait_for(function() file.exists(slow$started), 120))
  logs_before <- list.files(tempdir(), "^pt-render-")
  expect_identical(console_of(w$cancel), "printed line\nPlotTrellis: hic (A) track failed: no file\nWarning: could not write the cache file\n")
  expect_length(intersect(list.files(tempdir(), "^pt-render-"), logs_before), 0)  # its log is removed
})

test_that("a failed track's line from a real render reaches the app's console", {
  skip_on_cran()
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  rows <- test_rows("mouse")
  specs <- list(single_track("hic", rows[["Ba"]], extra = list(path = tempfile(fileext = ".hic"))))
  settings <- c(test_settings("mouse"), list(output_path = tempfile(fileext = ".svg")))
  on.exit(unlink(settings$output_path), add = TRUE)
  w <- .pt_render_worker()
  on.exit(w$close(), add = TRUE)
  res <- list()
  w$render(specs, test_region("mouse"), settings, function(err) res$err <<- list(err))
  shown <- console_of(function() expect_true(wait_for(function() length(res) > 0, 180)))
  expect_null(res$err[[1]])
  expect_match(shown, "PlotTrellis: hic (Ba) track failed:", fixed = TRUE)
})
