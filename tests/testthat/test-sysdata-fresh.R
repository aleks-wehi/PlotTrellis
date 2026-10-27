# R/sysdata.rda's PT_SOURCE_TEXT is a build-time copy of every top-level
# function's source text (data-raw/build-source-text.R), which downloaded R
# scripts are assembled from. Nothing rebuilds it automatically, so an edit
# to R/ without a rebuild silently ships the old text. Runs from a source
# tree only (the sync gate's load_all): data-raw/ isn't in the built package.
test_that("sysdata.rda's PT_SOURCE_TEXT matches the functions in R/", {
  pkg <- normalizePath(test_path("..", ".."))
  script <- file.path(pkg, "data-raw", "build-source-text.R")
  skip_if_not(file.exists(script), "data-raw/ not available (built package)")

  # Everything in the build script except its final save().
  exprs <- Filter(function(e) !(is.call(e) && identical(e[[1]], as.name("save"))), as.list(parse(script)))
  env <- new.env()
  old_wd <- setwd(pkg)
  on.exit(setwd(old_wd), add = TRUE)
  for (e in exprs) eval(e, env)

  fresh <- env$PT_SOURCE_TEXT
  stale <- union(setdiff(names(fresh), names(PT_SOURCE_TEXT)), setdiff(names(PT_SOURCE_TEXT), names(fresh)))
  common <- intersect(names(fresh), names(PT_SOURCE_TEXT))
  stale <- c(stale, common[!vapply(common, function(n) identical(fresh[[n]], PT_SOURCE_TEXT[[n]]), logical(1))])
  expect_true(length(stale) == 0, info = sprintf(
    "stale in R/sysdata.rda: %s -- rebuild with `Rscript --vanilla data-raw/build-source-text.R` in a clean `git archive` export",
    paste(stale, collapse = ", ")
  ))
})
