# Phase G.1 step 1 verify: assemble_download_script()'s discovery/pasting
# mechanism has no automated coverage otherwise (the 3 bugs fixed in
# 205d546 -- backtick-quoting `%||%`, all.names=TRUE for dot-prefixed
# helpers, environment(sys.function()) instead of a hardcoded namespace --
# were all found by manually running a generated script, not a test). This
# doesn't re-do that full end-to-end execution check (real UniProt/data
# access, a real installed-with-srcref build); it's a cheap regression guard
# against the exact 3 bugs recurring.

test_that("assemble_download_script() produces syntactically valid, self-contained R code", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "no /vast test data on this machine")

  region_params <- test_region("mouse")
  settings <- test_settings("mouse")
  entries <- list(list(viz = "genes", height = NULL, data = NULL))

  script <- assemble_download_script(entries, region_params, settings)

  parsed <- NULL
  expect_no_error(parsed <- parse(text = script))

  # Re-deparse the parsed expressions rather than grepping the raw pasted
  # source text: every name in the script is backtick-quoted (so
  # ".pt_dispatch_fn <- function" is really "`.pt_dispatch_fn` <- function"),
  # and srcref-recovered functions keep their original comments (one of
  # which *names* the old asNamespace("PlotTrellis") pattern to explain why
  # it's no longer used) -- deparse() strips both quirks, leaving only
  # actual code to check.
  code_only <- paste(vapply(parsed, function(e) paste(deparse(e), collapse = "\n"), character(1)), collapse = "\n")

  # Regression guard for the backtick-quoting fix: a bare `%||% <- ...`
  # (unquoted infix-operator LHS) is a parse error, so this also implicitly
  # covers itself via the parse() check above -- checked explicitly too so a
  # future change that silently drops %||% entirely (rather than
  # mis-quoting it) still fails this test.
  expect_true(grepl("%||%", code_only, fixed = TRUE))

  # Regression guard for the all.names=TRUE fix: dot-prefixed internal
  # helpers must be discovered and pasted in, not silently dropped.
  expect_true(grepl(".pt_dispatch_fn", code_only, fixed = TRUE))

  # Regression guard for the environment(sys.function()) fix: the assembled
  # script's actual CODE (not a comment describing the old pattern) must not
  # call asNamespace("PlotTrellis"), which doesn't exist once pasted
  # standalone with no PlotTrellis package installed.
  expect_false(grepl('asNamespace("PlotTrellis")', code_only, fixed = TRUE))
})
