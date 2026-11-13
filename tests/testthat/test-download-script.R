# assemble_download_script()'s discovery/pasting mechanism: past bugs here
# (backtick-quoting `%||%`, all.names=TRUE for dot-prefixed helpers,
# environment(sys.function()) instead of a hardcoded namespace) were all found
# by manually running a generated script. The checks below parse one-track
# scripts for every track type with fixture data, read their data block back,
# and run one script end to end with Rscript.

test_that("assemble_download_script() produces syntactically valid, self-contained R code", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "no /vast test data on this machine")

  region_params <- test_region("mouse")
  settings <- test_settings("mouse")
  entries <- list(list(viz = "genes", height = NULL, data = NULL))

  script <- assemble_download_script(entries, region_params, settings)

  parsed <- NULL
  expect_no_error(parsed <- parse(text = script))

  # Re-deparse the parsed expressions rather than grepping the raw pasted
  # source text: pasted functions keep their original comments (one of
  # which *names* asNamespace("PlotTrellis") to explain why it isn't
  # used); deparse() drops the comments, leaving only
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
  # script's actual CODE (not a comment naming it) must not
  # call asNamespace("PlotTrellis"), which doesn't exist once pasted
  # standalone with no PlotTrellis package installed.
  expect_false(grepl('asNamespace("PlotTrellis")', code_only, fixed = TRUE))
})

test_that("the script opens with a header naming the version, then the editable block", {
  script <- assemble_download_script(list(list(viz = "genes", height = NULL, data = NULL)),
                                     list(chrom = "chr2", start = 1000L, end = 5000L, gene = "Elf5"),
                                     list(orgdb_pkg = "org.Mm.eg.db", genome_label = "mm10"))
  lines <- strsplit(script, "\n")[[1]]
  header <- lines[seq_len(match(FALSE, startsWith(lines, "#")) - 1)]
  expect_true(length(header) >= 10 && length(header) <= 15)
  expect_true(any(grepl(paste("PlotTrellis", utils::packageVersion("PlotTrellis")), header, fixed = TRUE)))
  expect_true(any(grepl("Elf5, chr2:1,000-5,000 (mm10)", header, fixed = TRUE)))
  expect_lt(match("entries <- list(", lines), grep("<- function", lines)[1])
  expect_lt(grep("^OUTPUT <- ", lines), match("entries <- list(", lines))
  # The block says which track types the script can draw.
  expect_true(any(grepl("# (genes, the genome axis): a track of another type", lines, fixed = TRUE)))
})

test_that("the header names tracks as the app labels them, one line whatever a sample is called", {
  entries <- list(
    list(viz = "bw", sample = "Basal\ncells", colour = "#000000", name = "atac", path = "a.bw"),
    list(viz = "tpm_usage", sample = "K562", colour = "#000000", quant_path = "q"),
    list(viz = "junctions_mirrored_diff", name = "junctions", control_label = "A", treatment_label = "B",
         control_path = "a.bed", treatment_path = "b.bed", control_colour = "#000000", treatment_colour = "#111111"),
    list(viz = "genes")
  )
  script <- assemble_download_script(entries, list(chrom = "chr2", start = 1000L, end = 5000L, gene = "Elf5"),
                                     list(orgdb_pkg = "org.Mm.eg.db", genome_label = "mm10"))
  lines <- strsplit(script, "\n")[[1]]
  expect_identical(grep("^# Tracks: ", lines, value = TRUE),
                   "# Tracks: atac (Basal cells), Transcript expression (K562), junctions mirrored diff (B vs A), Gene annotation")
  expect_no_error(parse(text = script))
})

# Every name exported to the search path once `pkgs` are attached in an Rscript
# session: R's default packages, the packages themselves and what they Depend on.
visible_names <- function(pkgs) {
  todo <- c("methods", "datasets", "utils", "grDevices", "graphics", "stats", pkgs)
  seen <- character(0)
  while (length(todo)) {
    p <- todo[1]
    todo <- todo[-1]
    if (p %in% seen) next
    seen <- c(seen, p)
    dep <- utils::packageDescription(p, fields = "Depends")
    if (!is.na(dep)) todo <- c(todo, setdiff(trimws(sub("\\(.*", "", strsplit(dep, ",")[[1]])), c("R", "")))
  }
  c(ls(baseenv(), all.names = TRUE), unlist(lapply(seen, getNamespaceExports)))
}

# One track of each type, on the human fixture: single-sample (or sample-free)
# when the fixture helpers can build it, else a GM12878-vs-K562 comparison.
one_track_entries <- function() {
  rows <- test_rows("human")
  keys <- .pt_viz_keys()
  entries <- lapply(keys, function(k) tryCatch(single_track(k, rows[["GM12878"]], settings = test_settings("human")),
    error = function(e) tryCatch(pair_track(k, rows[["GM12878"]], rows[["K562"]]), error = function(e) NULL)))
  Filter(Negate(is.null), setNames(entries, keys))
}

test_that("every track type's script uses only what it defines, loads or gets from base R", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "no /vast test data on this machine")
  region_params <- test_region("human")
  settings <- test_settings("human")
  entries <- one_track_entries()
  expect_gt(length(entries), 10)

  for (key in names(entries)) {
    script <- assemble_download_script(list(entries[[key]]), region_params, settings)
    exprs <- parse(text = script, keep.source = FALSE)
    attached <- regmatches(script, gregexpr("(?m)^\\s*library\\(\\K[A-Za-z0-9.]+(?=\\))", script, perl = TRUE))[[1]]
    checked <- eval(Find(function(e) is.call(e) && identical(e[[1]], as.name("for")), exprs)[[3]], baseenv())
    # The whole script as one function body: findGlobals() then reports exactly
    # the names it uses but never assigns.
    whole <- function() NULL
    body(whole) <- as.call(c(as.name("{"), as.list(exprs)))
    used <- codetools::findGlobals(whole, merge = TRUE)
    # Names a pasted function writes into the global environment with assign("<name>", ...).
    assigned <- regmatches(script, gregexpr('assign\\("\\K[A-Za-z0-9._]+(?=")', script, perl = TRUE))[[1]]
    missing <- setdiff(used, c(visible_names(attached), assigned))
    expect_identical(missing, character(0), info = key)
    # ...and attaches nothing it doesn't use.
    for (p in attached) expect_true(any(used %in% getNamespaceExports(p)), info = paste(key, p))

    pd <- utils::getParseData(parse(text = script, keep.source = TRUE))
    colons <- unique(pd$text[pd$token == "SYMBOL_PACKAGE"])
    is_base <- vapply(colons, function(p) identical(utils::packageDescription(p, fields = "Priority"), "base"), logical(1))
    expect_identical(setdiff(colons[!is_base], checked), character(0), info = key)
    expect_true(settings$orgdb_pkg %in% checked, info = key)
    # The figure's own track type is pasted; render_all() comes first, the call last.
    expect_true(grepl(paste0("\nrender_", key, " <- function"), script, fixed = TRUE), info = key)
    expect_lt(regexpr("\nrender_all <- function", script, fixed = TRUE),
              regexpr(paste0("\nsettings_", key, " <- function"), script, fixed = TRUE))
  }
})

test_that("a script stops with an install hint when a package it needs is missing", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "no /vast test data on this machine")
  settings <- test_settings("mouse", overrides = list(orgdb_pkg = "org.Xx.eg.db"))
  script <- assemble_download_script(list(list(viz = "genes")), test_region("mouse"), settings)
  check <- Find(function(e) is.call(e) && identical(e[[1]], as.name("for")), parse(text = script, keep.source = FALSE))
  expect_error(eval(check, new.env()),
               "needs the R package 'org.Xx.eg.db'. Install it with BiocManager::install(\"org.Xx.eg.db\")", fixed = TRUE)
})

test_that("a gzipped BED track needs R.utils: imported by the package and checked by its script", {
  expect_match(utils::packageDescription("PlotTrellis")$Imports, "R.utils", fixed = TRUE)
  skip_if_not(PT_TEST_DATA_AVAILABLE, "no /vast test data on this machine")
  script <- assemble_download_script(list(list(viz = "bed", name = "peaks", path = "peaks.bed.gz")),
                                     test_region("human"), test_settings("human"))
  check <- Find(function(e) is.call(e) && identical(e[[1]], as.name("for")), parse(text = script, keep.source = FALSE))
  expect_true("R.utils" %in% eval(check[[3]], baseenv()))
  # A plain .bed file or another track type doesn't ask for it.
  plain <- assemble_download_script(list(list(viz = "bed", name = "peaks", path = "peaks.bed")),
                                    test_region("human"), test_settings("human"))
  needed <- function(s) eval(Find(function(e) is.call(e) && identical(e[[1]], as.name("for")),
                                    parse(text = s, keep.source = FALSE))[[3]], baseenv())
  expect_false("R.utils" %in% needed(plain))
  other <- assemble_download_script(list(list(viz = "genes")), test_region("human"), test_settings("human"))
  expect_false("R.utils" %in% needed(other))
})

test_that("load_bed() reads a gzipped BED file's chromosomes", {
  f <- tempfile(fileext = ".bed.gz")
  con <- gzfile(f, "w")
  writeLines(c("chr1\t100\t200", "chr1\t300\t400"), con)
  close(con)
  out <- load_bed(list(path = f), list(chrom = "chr1"))
  expect_false(out$data$no_chrom)
  expect_equal(out$data$chrom, "chr1")
})

# The data block's three objects, read back from the script text.
data_block <- function(script) {
  exprs <- parse(text = script, keep.source = FALSE)
  env <- new.env()
  for (e in exprs) {
    if (is.call(e) && identical(e[[1]], as.name("<-")) &&
        as.character(e[[2]]) %in% c("OUTPUT", "entries", "region_params", "settings")) eval(e, env)
  }
  as.list(env)
}

render_svg <- function(entries, region_params, settings) {
  f <- tempfile(fileext = ".svg")
  # A Hi-C track with no contacts in the region warns; both renders do alike.
  suppressWarnings(.pt_quiet(render_all(entries, region_params, settings,
                                        device_fn = function(w, h) svglite::svglite(f, width = w, height = h))))
  grDevices::dev.off()
  readBin(f, "raw", file.size(f))
}

test_that("the editable block reads back as the figure's tracks, region and settings, and draws the same", {
  skip_if_not(PT_TEST_DATA_AVAILABLE, "no /vast test data on this machine")
  skip_on_cran()
  region_params <- test_region("human")
  settings <- test_settings("human")
  entries <- unname(one_track_entries())

  back <- data_block(assemble_download_script(entries, region_params, settings))
  # NULL fields are left out of a track, and nothing else.
  expect_identical(back$entries, lapply(entries, function(e) e[!vapply(e, is.null, logical(1))]))
  expect_identical(back$region_params, region_params)
  # Every setting but the app's output_path and hic_style (each Hi-C track carries its own),
  # value for value; only the order differs.
  expected <- settings[!names(settings) %in% c("output_path", "hic_style")]
  expect_setequal(names(back$settings), names(expected))
  expect_identical(back$settings[names(expected)], expected)

  expect_identical(render_svg(back$entries, back$region_params, back$settings),
                   render_svg(entries, region_params, settings))
})

test_that("a comparison track carries its name once", {
  items <- list(cmp = list(key = "cmp", label = "x", type = "bw_diff_stacked", name = "atac",
                           samples = c("A", "B"), control_label = "A", treatment_label = "B",
                           control_path = "a.bw", treatment_path = "b.bw"))
  entry <- .pt_build_entries("cmp", items, character(0), list())[[1]]
  expect_identical(sum(names(entry) == "name"), 1L)
  expect_identical(entry$name, "atac")
})

test_that("a setting that is drawn from each track's own copy is not offered in the settings block", {
  script <- assemble_download_script(
    list(list(viz = "hic", sample = "S", colour = "#000000", name = "hic", path = "a.hic", hic_style = "square")),
    list(chrom = "chr2", start = 1000L, end = 5000L, gene = NA_character_),
    list(orgdb_pkg = "org.Mm.eg.db", genome_label = "mm10", hic_style = "rectangle", hic_norm = "KR"))
  back <- data_block(script)
  expect_identical(back$entries[[1]]$hic_style, "square")
  expect_false("hic_style" %in% names(back$settings))
  expect_identical(back$settings$hic_norm, "KR")
})

test_that("the bundled example's script runs alone with Rscript and draws pt_save()'s SVG", {
  skip_on_cran()
  skip_if_not_installed("org.Mm.eg.db")
  proj <- suppressWarnings(pt_project(file.path(tempdir(), "pt-test-script-mammary"), pt_example("mammary")))
  # A coverage track first (it draws the scale bar) and two Transcript expression tracks (they
  # fold into one panel): both come from the track types' declarations, which the script
  # must carry.
  fig <- pt_figure(proj, "Elf5", c("atac@Basal", "genes", "transcript_expression@Basal", "transcript_expression@LP"))
  api_svg <- tempfile(fileext = ".svg")
  pt_save(fig, api_svg)

  run_dir <- tempfile("pt-script-run-")
  dir.create(run_dir)
  script <- pt_script(fig, file.path(run_dir, "figure.R"))
  head <- readLines(script, n = 15)
  expect_true(any(grepl(paste("PlotTrellis", utils::packageVersion("PlotTrellis")), head, fixed = TRUE)))
  expect_true(any(grepl("Elf5", head, fixed = TRUE)))

  old <- setwd(run_dir)
  on.exit(setwd(old), add = TRUE)
  out <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", "figure.R"), stdout = TRUE, stderr = TRUE,
                 env = paste0("R_LIBS=", paste(.libPaths(), collapse = .Platform$path.sep)))
  expect_null(attr(out, "status"))
  svg <- list.files(run_dir, pattern = "\\.svg$")
  expect_length(svg, 1)
  expect_identical(readBin(svg, "raw", file.size(svg)), readBin(api_svg, "raw", file.size(api_svg)))
  # A quiet run: no per-object plotgardener lines, no "null device", one closing line.
  expect_false(any(grepl("^\\w+\\[\\w+\\]\\s*$|null device", out)))
  expect_identical(sum(grepl("^Figure written to", out)), 1L)
})
