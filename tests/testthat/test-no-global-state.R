# No package function writes to the user's global environment or attaches a package,
# and the TxDb and its CDS index are opened or built once per genome per session
# (.pt_annotation(), R/boilerplate.R).

# Renders to a throwaway SVG; returns the messages render_all() gave.
render_messages <- function(items, region, settings) {
  f <- tempfile(fileext = ".svg")
  on.exit(unlink(f), add = TRUE)
  dev <- NULL
  on.exit(if (!is.null(dev) && dev %in% grDevices::dev.list()) grDevices::dev.off(dev), add = TRUE, after = FALSE)
  msgs <- character(0)
  withCallingHandlers(
    render_all(items, region, settings, device_fn = function(w, h) {
      svglite::svglite(f, width = w, height = h)
      dev <<- grDevices::dev.cur()
    }),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  msgs
}

# UniProt is never contacted: every read fails as if offline, which is never cached.
offline_uniprot <- function(url) stop(errorCondition("offline in this test", class = "uniprot_unreachable"))

test_that("a render of every transcript-type track leaves globalenv and search() as they were (human dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  rows <- test_rows("human")
  user <- list(txdb = "my own txdb", MIN_DISPLAY_TPM = 99)
  for (n in names(user)) assign(n, user[[n]], envir = globalenv())
  on.exit(rm(list = names(user), envir = globalenv()), add = TRUE)
  local_mocked_bindings(.uniprot_read_lines = offline_uniprot)
  items <- list(
    single_track("transcripts_all", rows[["GM12878"]]),
    single_track("transcripts", rows[["GM12878"]]),
    single_track("transcripts_annotated", rows[["GM12878"]]),
    single_track("tpm_usage", rows[["GM12878"]]),
    single_track("tpm_usage", rows[["K562"]]),
    pair_track("junctions_mirrored_diff_tpm_scaled", rows[["GM12878"]], rows[["K562"]]),
    single_track("protein_domains", rows[["GM12878"]])
  )
  globals  <- function() setdiff(ls(globalenv(), all.names = TRUE), ".Random.seed")
  before   <- globals()
  attached <- search()
  msgs <- render_messages(items, test_region("human"), test_settings("human"))
  expect_false(any(grepl("track failed", msgs)), info = paste(msgs, collapse = "\n"))
  expect_identical(globals(), before)
  expect_identical(mget(names(user), envir = globalenv()), user)
  expect_identical(search(), attached)
})

# Canned UniProt answers: a reviewed hit with one domain, so the protein-domains lookup
# goes on to need the CDS index, and a protein of `aa` amino acids.
canned_uniprot <- function(aa = 99) function(url) {
  if (grepl("/search?", url, fixed = TRUE)) return("Q99999")
  if (grepl("\\.txt$", url)) return(c(
    "ID   TEST_PROTEIN            Reviewed;          99 AA.",
    "FT   DOMAIN          1..10",
    "FT                   /note=\"Test domain\""))
  if (grepl("\\.fasta$", url)) return(c(">sp|Q99999|TEST_PROTEIN", strrep("M", aa)))
  stop("unexpected URL in test: ", url)
}

test_that("a render opens the TxDb at most once and two genes' protein tracks build the CDS index once (mouse dataset)", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  rows <- test_rows("mouse")
  cache_dir <- tempfile("pt-cache-")
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)
  settings <- test_settings("mouse", overrides = list(cache_dir = cache_dir))
  local_mocked_bindings(.uniprot_read_lines = canned_uniprot())
  calls <- new.env()
  calls$loadDb <- 0L
  calls$cdsBy <- 0L
  real_load_db <- AnnotationDbi::loadDb
  real_cds_by  <- GenomicFeatures::cdsBy
  local_mocked_bindings(loadDb = function(...) { calls$loadDb <- calls$loadDb + 1L; real_load_db(...) },
                        .package = "AnnotationDbi")
  local_mocked_bindings(cdsBy = function(...) { calls$cdsBy <- calls$cdsBy + 1L; real_cds_by(...) },
                        .package = "GenomicFeatures")
  rm(list = ls(.pt_txdb_memo), envir = .pt_txdb_memo)  # nothing open yet this session
  items <- list(single_track("transcripts", rows[["Ba"]]), single_track("tpm_usage", rows[["Ba"]]),
                single_track("protein_domains", rows[["Ba"]]))
  for (gene in c("Elf5", "Ehf")) {  # neither in the (empty) protein-domain cache
    before <- calls$loadDb
    msgs <- render_messages(items, test_region("mouse", gene), settings)
    expect_false(any(grepl("track failed", msgs)), info = paste(msgs, collapse = "\n"))
    expect_lte(calls$loadDb - before, 1L)
  }
  expect_identical(calls$loadDb, 1L)
  expect_identical(calls$cdsBy, 1L)
})

test_that("an annotation file replaced on disk is opened again, and one left alone is not", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")
  path <- tempfile(fileext = ".sqlite")
  on.exit(unlink(path), add = TRUE)
  n_tx <- function(txdb) length(GenomicFeatures::transcripts(txdb))
  file.copy(GENOMES$mouse$txdb_path, path)
  first <- .pt_annotation(list(txdb_path = path))$txdb
  expect_identical(n_tx(first), n_tx(TXDBS$mouse))
  expect_identical(.pt_annotation(list(txdb_path = path))$txdb, first)
  file.copy(GENOMES$human$txdb_path, path, overwrite = TRUE)
  Sys.setFileTime(path, Sys.time() + 60)
  expect_identical(n_tx(.pt_annotation(list(txdb_path = path))$txdb), n_tx(TXDBS$human))
  expect_false(n_tx(TXDBS$human) == n_tx(TXDBS$mouse))
})

test_that("a project, a figure and a saved figure attach nothing and leave nothing behind, in a fresh R session", {
  skip_on_cran()
  skip_if_not_installed("org.Mm.eg.db")
  # The package as this test session has it: the source tree when there is one (load_all),
  # else the installed package. Without the test helpers: they open annotation databases
  # themselves, which would attach the packages this test says the package leaves alone.
  src <- normalizePath(test_path("..", ".."))
  load <- if (file.exists(file.path(src, "DESCRIPTION")) && requireNamespace("pkgload", quietly = TRUE))
    sprintf("suppressMessages(pkgload::load_all(%s, quiet = TRUE, helpers = FALSE))", encodeString(src, quote = "\""))
  else "library(PlotTrellis)"
  run_dir <- tempfile("pt-fresh-session-")
  dir.create(run_dir)
  on.exit(unlink(run_dir, recursive = TRUE), add = TRUE)
  writeLines(c(
    load,
    "attached <- search()",
    "objects <- ls(globalenv(), all.names = TRUE)",
    "proj <- pt_project(file.path(getwd(), 'proj'), pt_example('mammary'))",
    "fig <- pt_figure(proj, 'Elf5', c('genes', 'transcripts', 'atac@Basal', 'transcript_expression@Basal'))",
    "pt_save(fig, 'figure.svg')",
    "cat('ATTACHED:', setdiff(search(), attached), '\n')",
    "cat('OBJECTS:', setdiff(ls(globalenv(), all.names = TRUE), c(objects, 'objects', 'proj', 'fig', '.Random.seed')), '\n')"
  ), file.path(run_dir, "session.R"))
  old <- setwd(run_dir)
  on.exit(setwd(old), add = TRUE)
  out <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"), c("--vanilla", "session.R"), stdout = TRUE, stderr = TRUE,
    timeout = 600,
    env = c(paste0("R_LIBS=", paste(.libPaths(), collapse = .Platform$path.sep)),
            paste0("R_USER_CACHE_DIR=", file.path(run_dir, "user-cache")))))
  expect_null(attr(out, "status"), info = paste(out, collapse = "\n"))
  expect_true(file.exists(file.path(run_dir, "figure.svg")))
  expect_identical(trimws(grep("^ATTACHED:", out, value = TRUE)), "ATTACHED:")
  expect_identical(trimws(grep("^OBJECTS:", out, value = TRUE)), "OBJECTS:")
  expect_false(any(grepl("Attaching package|masked from|Loading required package", out)),
               info = paste(out, collapse = "\n"))
})

test_that("two projects in one session keep their own caches, and a figure's script names its own project's", {
  skip_if_not_installed("org.Hs.eg.db")
  # Real path: pt_project() stores its folder normalized, and macOS's tempdir is under a
  # symlink (/var -> /private/var), so the script names the resolved path.
  root <- tempfile("pt-two-")
  dir.create(root)
  root <- normalizePath(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  a_dir <- file.path(root, "A")
  b_dir <- file.path(root, "B")
  A <- suppressMessages(suppressWarnings(pt_project(a_dir, pt_example("human"))))
  B <- suppressMessages(suppressWarnings(pt_project(b_dir, pt_example("human"))))
  b_files <- function() list.files(file.path(b_dir, "cache"), recursive = TRUE, all.files = TRUE)
  b_before <- b_files()
  # UniProt answers with a protein as long as one of PSME2's coding transcripts, so the
  # protein-domains lookup succeeds and is cached.
  anno <- .pt_annotation(A$render_cfg)
  tx <- .pd_gene_transcript_tags("PSME2", A$genome$tx_meta_lookup_path, anno$txdb, A$genome$orgdb)$transcript_id
  cds <- anno$cds_by_tx()
  aa <- sum(GenomicRanges::width(cds[[intersect(tx, names(cds))[1]]])) / 3 - 1
  local_mocked_bindings(.uniprot_read_lines = canned_uniprot(aa))
  fig <- pt_figure(A, "PSME2", c("genes", "transcripts", "protein_domains"))
  suppressMessages(pt_save(fig, file.path(root, "A.svg")))
  expect_setequal(b_files(), b_before)
  expect_true(file.exists(file.path(a_dir, "cache", "gene_display_transcripts_cache.rds")))
  expect_true(file.exists(file.path(a_dir, "cache", "protein_domains", A$genome$genome_label, "PSME2.rds")))
  pt_script(fig, file.path(root, "A.R"))
  script <- paste(readLines(file.path(root, "A.R")), collapse = "\n")
  expect_true(grepl(file.path(a_dir, "cache"), script, fixed = TRUE))
  expect_false(grepl(b_dir, script, fixed = TRUE))
})

# ---- Static checks over every function in the package namespace (no data needed) ----

ns_functions <- function() {
  ns <- asNamespace("PlotTrellis")
  Filter(is.function, mget(ls(ns, all.names = TRUE), envir = ns))
}

# Every call in a function (formals and body, nested functions included), as a list of calls.
calls_in <- function(f) {
  out <- list()
  walk <- function(x) {
    if (is.call(x)) out[[length(out) + 1]] <<- x
    if (is.call(x) || is.pairlist(x) || is.expression(x))
      for (i in seq_along(x)) if (!identical(x[[i]], quote(expr = ))) walk(x[[i]])
  }
  walk(formals(f))
  walk(body(f))
  out
}

# A call's function name: `f(...)`, `pkg::f(...)` and `pkg:::f(...)` all give "f".
call_name <- function(cl) {
  f <- cl[[1]]
  if (is.name(f)) return(as.character(f))
  if (is.call(f) && as.character(f[[1]]) %in% c("::", ":::")) return(as.character(f[[3]]))
  ""
}

test_that("no package function reads a variable or function that is bound nowhere", {
  ns  <- asNamespace("PlotTrellis")
  imp <- parent.env(ns)
  known <- function(n) exists(n, envir = ns, inherits = FALSE) || exists(n, envir = imp, inherits = FALSE) ||
    exists(n, envir = baseenv(), inherits = FALSE)
  fns <- ns_functions()
  free <- as.character(unlist(lapply(names(fns), function(nm) {
    g <- codetools::findGlobals(fns[[nm]], merge = TRUE)
    miss <- g[!vapply(g, known, logical(1))]
    if (length(miss)) paste0(nm, ": ", paste(miss, collapse = ", "))
  })))
  expect_identical(free, character(0), info = paste(free, collapse = "; "))
})

# The three checks below are tripwires for the ways this package once leaked state, and for
# their obvious variants. They are strict on purpose (a function may not name the global
# environment at all, even to read from it) and they are not a proof: code that builds an
# environment or a function name at run time gets past them.

test_that("no package function sets environment variables, attaches packages, unlocks bindings or changes folder", {
  banned <- c("Sys.setenv", "Sys.unsetenv", "library", "require", "attach", "attachNamespace", "unlockBinding",
              "setwd")
  # A banned function called by name, or looked up by name or string and then called.
  lookups <- c("do.call", "match.fun", "get", "get0", "getExportedValue", "getFromNamespace")
  banned_in <- function(cl) {
    nm <- call_name(cl)
    if (nm %in% banned) return(nm)
    if (nm %in% lookups) {
      args <- Filter(function(a) is.name(a) || is.character(a), as.list(cl)[-1])
      return(intersect(vapply(args, as.character, ""), banned))
    }
    character(0)
  }
  fns <- ns_functions()
  hits <- as.character(unlist(lapply(names(fns), function(nm) {
    used <- unique(unlist(lapply(calls_in(fns[[nm]]), banned_in)))
    if (length(used)) paste0(nm, ": ", paste(used, collapse = ", "))
  })))
  expect_identical(hits, character(0), info = paste(hits, collapse = "; "))
})

test_that("no package function names the global environment, another frame's, or one by its search-path position", {
  # Environments a package function has no business holding: the global one, its caller's
  # (the global one, when the caller is the console), and one picked by position.
  elsewhere <- c("globalenv", "parent.frame", "parent.env", "eval.parent", "as.environment", "sys.frame",
                 "sys.frames", "sys.function", "topenv", "pos.to.env")
  # .pt_module_fn() and the functions that list modules read environment(sys.function()):
  # their own enclosing environment, which is the namespace (or a downloaded script's).
  own_env <- function(cl) call_name(cl) == "environment" && length(cl) == 2 && identical(cl[[2]], quote(sys.function()))
  fns <- ns_functions()
  hits <- as.character(unlist(lapply(names(fns), function(nm) {
    calls <- calls_in(fns[[nm]])
    allowed <- sum(vapply(calls, own_env, logical(1)))
    names_used <- vapply(calls, call_name, "")
    used <- intersect(names_used, elsewhere)
    if (sum(names_used == "sys.function") <= allowed) used <- setdiff(used, "sys.function")
    if (any(vapply(calls, function(cl) any(vapply(as.list(cl), identical, logical(1), quote(.GlobalEnv))), logical(1))))
      used <- c(used, ".GlobalEnv")
    # assign(x, value, pos): a third argument given by position, or one named pos.
    by_pos <- any(vapply(calls, function(cl) {
      if (call_name(cl) != "assign") return(FALSE)
      args <- as.list(cl)[-1]
      nms <- names(args) %||% rep("", length(args))
      "pos" %in% nms || sum(!nzchar(nms)) >= 3
    }, logical(1)))
    if (by_pos) used <- c(used, "assign(pos)")
    if (length(used)) paste0(nm, ": ", paste(used, collapse = ", "))
  })))
  expect_identical(hits, character(0), info = paste(hits, collapse = "; "))
})
