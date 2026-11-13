# pt_figure()'s track grammar (.pt_resolve_tracks()) on a hand-made track
# table: pure, no data needed. The end-to-end example-data tests are below.

tab <- data.frame(
  track  = c("genes", "transcripts", "atac", "atac", "atac", "junctions", "junctions_diff", "junctions_diff"),
  sample = c(NA, NA, "A", "B", "C", "A", "B vs A", "C vs A"),
  kind   = c("fixed", "fixed", rep("sample", 4), "comparison", "comparison"),
  label  = "",
  key    = c("genes", "transcripts", "atac_A", "atac_B", "atac_C", "junc_A", "diff_A_vs_B", "diff_A_vs_C"),
  stringsAsFactors = FALSE
)

test_that("explicit name@sample picks exactly that track, in the order given", {
  expect_equal(.pt_resolve_tracks(tab, c("genes", "atac@B", "atac@A", "junctions_diff@C vs A"), NULL),
               c("genes", "atac_B", "atac_A", "diff_A_vs_C"))
  # Spacing around @ and vs doesn't matter.
  expect_equal(.pt_resolve_tracks(tab, "junctions_diff @ B  vs  A", NULL), "diff_A_vs_B")
})

test_that("parallel expands bare per-sample names track-major, skipping samples without it", {
  expect_message(keys <- .pt_resolve_tracks(tab, c("genes", "atac", "junctions"), c("A", "B")),
                 "no 'junctions' for B")
  expect_equal(keys, c("genes", "atac_A", "atac_B", "junc_A"))
  # An explicit @ still means just that one when running in parallel.
  expect_equal(.pt_resolve_tracks(tab, c("atac@C", "atac"), c("A", "B")), c("atac_C", "atac_A", "atac_B"))
})

test_that("a bare comparison uses the 2 paralleld samples, the first as reference", {
  expect_equal(.pt_resolve_tracks(tab, "junctions_diff", c("A", "C")), "diff_A_vs_C")
  expect_error(.pt_resolve_tracks(tab, "junctions_diff", c("A", "B", "C")), "compares 2 samples")
})

test_that("mistakes stop with a message saying what to write instead", {
  expect_error(.pt_resolve_tracks(tab, "atac", NULL), "atac@A.*parallel = TRUE")
  expect_error(.pt_resolve_tracks(tab, "ATAC", NULL), "Unknown track 'ATAC'.*atac")
  expect_error(.pt_resolve_tracks(tab, "genes@A", NULL), "isn't per-sample")
  expect_error(.pt_resolve_tracks(tab, "atac@D", NULL), "No 'atac' track for 'D'")
  expect_error(.pt_resolve_tracks(tab, "junctions_diff@A vs B", NULL), "No 'junctions_diff' track for 'A vs B'")
})

test_that("transcript_counts, the old name of the quantification track, finds it even when another type took the new name", {
  tab2 <- data.frame(track = c("transcript_expression", "transcript_expression_tpm_usage"),
                     sample = c("A", "A"), kind = "sample", label = "",
                     key = c("bw_transcript_expression_A", "tpm_usage_A"), stringsAsFactors = FALSE)
  expect_equal(.pt_resolve_tracks(tab2, "transcript_counts@A", NULL), "tpm_usage_A")
  expect_equal(.pt_resolve_tracks(tab2, "transcript_expression@A", NULL), "bw_transcript_expression_A")
})

# End to end on the shipped example data (inst/extdata/example-human), so these run
# everywhere, R-universe included. Needs org.Hs.eg.db (Suggests). The data is read
# in place; the project (settings, caches) lives in EX_DIR.
EX_DIR <- file.path(tempdir(), "pt-test-example-human")
test_that("the example goes from pt_example(human) to a saved figure", {
  skip_if_not_installed("org.Hs.eg.db")
  # First open builds the example's TxDb; txdbmaker warns about the GTF's
  # stop_codon phase and missing genome version, both harmless.
  attached <- search()
  proj <- suppressWarnings(pt_project(EX_DIR, pt_example("human")))
  expect_s3_class(proj, "pt_project")
  expect_setequal(unique(pt_tracks(proj)$track),
                  c("genes", "protein_domains", "transcripts_all", "transcripts", "transcripts_annotated",
                    "atac", "junctions", "transcript_expression", "atac_diff", "junctions_diff", "junctions_diff_tpm"))
  fig <- pt_figure(proj, "psme2", c("genes", "atac", "junctions_diff"), parallel = TRUE)
  expect_equal(fig$gene, "PSME2")
  expect_equal(fig$keys, c("genes", "bw_atac_GM12878", "bw_atac_K562", "junctions_mirrored_diff_junctions_GM12878_vs_K562"))
  png <- pt_save(fig, tempfile(fileext = ".png"), dpi = 50)
  # rtracklayer and plotgardener cannot read bigWig on Windows, so the
  # atac tracks draw empty there and the PNG is small.
  if (.Platform$OS.type != "windows") expect_gt(file.size(png), 10000)
  expect_true(file.exists(png))
  # The project's OrgDb is used through its namespace, never attached to the search path.
  expect_false(any(grepl("^package:org\\.", setdiff(search(), attached))))
  expect_gt(length(readLines(pt_script(fig, tempfile(fileext = ".R")))), 100)
  # Preloaded for the app, and every preloaded track exists.
  sv <- .pt_project(proj$dir)$cfg$start_view
  expect_equal(sv$gene, "PSME2")
  expect_true(all(unlist(sv$tracks) %in% names(proj$items)))
})

test_that("transcript_counts, the old name of transcript_expression, still draws the same figure", {
  skip_if_not_installed("org.Hs.eg.db")
  proj <- suppressWarnings(pt_project(EX_DIR, pt_example("human")))
  tr <- pt_tracks(proj)
  expect_true("transcript_expression" %in% tr$track)
  expect_false("transcript_counts" %in% tr$track)           # listed under the new name only
  expect_match(tr$label[tr$track == "transcript_expression"][1], "^Transcript expression")
  new <- pt_figure(proj, "PSME2", c("genes", "transcript_expression@GM12878"))
  old <- expect_silent(pt_figure(proj, "PSME2", c("genes", "transcript_counts@GM12878")))   # no warning or message
  expect_identical(old$keys, new$keys)
  expect_identical(old$entries, new$entries)
  expect_identical(pt_figure(proj, "PSME2", "transcript_counts", parallel = TRUE)$keys,
                   pt_figure(proj, "PSME2", "transcript_expression", parallel = TRUE)$keys)
})

test_that("no 'Transcript counts' is left in text a user sees", {
  # The track's label, as the composer and the downloaded script show it.
  expect_identical(.pt_viz_info("tpm_usage")$label, "Transcript expression")
  pkg <- testthat::test_path("..", "..")
  skip_if_not(file.exists(file.path(pkg, "DESCRIPTION")) && dir.exists(file.path(pkg, "R")), "package source not available (installed package)")
  files <- c(list.files(file.path(pkg, "R"), "\\.R$", full.names = TRUE),
             list.files(file.path(pkg, "man"), "\\.Rd$", full.names = TRUE),
             list.files(file.path(pkg, "vignettes"), "\\.Rmd$", full.names = TRUE),
             list.files(file.path(pkg, "inst", "extdata"), "samples\\.txt$", recursive = TRUE, full.names = TRUE))
  expect_gt(length(files), 20)
  for (f in files) {
    hits <- grep("transcript counts", readLines(f, warn = FALSE), ignore.case = TRUE, value = TRUE)
    expect_identical(hits, character(0), info = basename(f))
  }
  # What users read or copy (examples, vignette, READMEs) also doesn't use the old code name or
  # the older "transcript-usage"; the alias is only documented in the manual page.
  copied <- c(list.files(file.path(pkg, "vignettes"), "\\.Rmd$", full.names = TRUE),
              list.files(file.path(pkg, "inst", "extdata"), "(samples\\.txt|README)$", recursive = TRUE, full.names = TRUE))
  expect_gt(length(copied), 3)
  for (f in copied) {
    hits <- grep("transcript[ _-](counts|usage)", readLines(f, warn = FALSE), ignore.case = TRUE, value = TRUE)
    expect_identical(hits, character(0), info = basename(f))
  }
})

test_that("pt_preload() stores the figure as the app's start view", {
  skip_if_not_installed("org.Hs.eg.db")
  proj <- pt_project(EX_DIR, pt_example("human"))
  json <- file.path(proj$dir, "project_settings.json")
  before <- readLines(json)
  on.exit(writeLines(before, json))
  pt_preload(proj, pt_figure(proj, "PSME2", c("genes", "atac@K562", "junctions_diff@K562 vs GM12878")))
  sv <- jsonlite::read_json(json)$start_view
  expect_equal(sv$view_mode, "two")
  expect_equal(unlist(sv$conditions), c("GM12878", "K562"))   # reference first
  expect_equal(sv$gene, "PSME2")
  expect_equal(unlist(sv$tracks), c("genes", "bw_atac_K562", "junctions_mirrored_diff_junctions_GM12878_vs_K562"))
  # The app reads it back through the normal settings reader.
  expect_equal(.pt_project(proj$dir)$cfg$start_view$gene, "PSME2")
})

test_that("pt_project() re-imports a samples.txt only when it is newer than the project", {
  skip_if_not_installed("org.Hs.eg.db")
  src <- pt_example("human")
  st  <- file.path(tempfile(), "samples.txt")
  dir.create(dirname(st))
  txt <- readLines(src)
  txt <- sub("^# data_root: .*", paste("# data_root:", dirname(pt_example("human"))), txt)
  txt <- sub("^# gtf_path: .*", paste("# gtf_path:", file.path(dirname(pt_example("human")), "annotation.gtf.gz")), txt)
  writeLines(txt, st)
  dir <- file.path(dirname(st), "proj")
  suppressWarnings(pt_project(dir, st))   # builds its own TxDb: harmless txdbmaker warnings
  json <- file.path(dir, "project_settings.json")
  Sys.setFileTime(st, Sys.time() - 60)
  expect_silent(stamp <- file.mtime(json))
  pt_project(dir, st)
  expect_equal(file.mtime(json), stamp)                        # older samples.txt: opened as is
  Sys.setFileTime(st, Sys.time() + 60)
  expect_message(pt_project(dir, st), "Re-importing")
})

test_that("a project keeps every cache in its own folder, none next to the data", {
  skip_if_not_installed("org.Hs.eg.db")
  data <- tempfile("pt-data-")
  dir.create(data)
  file.copy(list.files(system.file("extdata", "example-human", package = "PlotTrellis"), full.names = TRUE),
            data, recursive = TRUE)
  before <- list.files(data, recursive = TRUE)
  proj_dir <- tempfile("pt-proj-")
  proj <- suppressWarnings(pt_project(proj_dir, file.path(data, "samples.txt")))
  # The data folder is untouched; the shipped chrom-sizes cache there is reused, not copied.
  expect_setequal(list.files(data, recursive = TRUE), before)
  expect_same_path(proj$genome$txdb_path, file.path(proj_dir, "cache", "genome", "annotation_txdb.sqlite"))
  expect_true(file.exists(proj$genome$txdb_path))
  expect_same_path(proj$render_cfg$cache_dir, file.path(proj_dir, "cache"))
  # A gene render's transcript-filter cache lands there too.
  pt_save(pt_figure(proj, "PSME2", c("genes", "transcripts")), tempfile(fileext = ".png"), dpi = 30)
  expect_true(file.exists(file.path(proj_dir, "cache", "gene_display_transcripts_cache.rds")))
})

test_that("a relative project_dir gives absolute cache paths, so a script or pt_save() works after setwd()", {
  skip_if_not_installed("org.Hs.eg.db")
  base <- tempfile("pt-rel-")
  dir.create(base)
  old <- setwd(base)
  on.exit(setwd(old), add = TRUE)
  proj <- suppressWarnings(pt_project("rel-project", pt_example("human")))
  fig <- pt_figure(proj, "PSME2", c("genes", "transcripts"))
  paths <- unlist(proj$render_cfg[c("txdb_path", "cache_dir", "tx_name_lookup_path", "tx_meta_lookup_path")])
  # Absolute: normalising again changes nothing (any OS), and none is the relative form.
  expect_identical(unname(paths), normalizePath(paths, winslash = "/", mustWork = FALSE))
  expect_false(any(startsWith(paths, "rel-project")))
  script <- readLines(pt_script(fig, "fig.R"))
  in_script <- sub('.*= "(.*)",?$', "\\1", grep("(txdb_path|cache_dir) = ", script, value = TRUE))
  expect_gt(length(in_script), 0)
  expect_false(any(startsWith(in_script, "rel-project")))
  setwd(tempdir())   # elsewhere: the project's files are still found
  expect_true(file.exists(pt_save(fig, file.path(base, "fig.png"), dpi = 30)))
})

test_that("one folder holds one project: a different samples.txt there is an error", {
  skip_if_not_installed("org.Hs.eg.db")
  mk <- function() {
    st <- file.path(tempfile(), "samples.txt")
    dir.create(dirname(st))
    txt <- readLines(pt_example("human"))
    txt <- sub("^# data_root: .*", paste("# data_root:", dirname(pt_example("human"))), txt)
    writeLines(sub("^# gtf_path: .*", paste("# gtf_path:", file.path(dirname(pt_example("human")), "annotation.gtf.gz")), txt), st)
    st
  }
  a <- mk(); b <- mk()
  dir <- tempfile("pt-one-")
  suppressWarnings(pt_project(dir, a))
  expect_error(pt_project(dir, b), "already holds a project, made from .*project_dir")
  expect_message(suppressWarnings(pt_project(dir, b, overwrite = TRUE)), "Re-importing")
  expect_error(pt_project(dir, a), "made from .*samples.txt")
})

test_that("bare pt_project() with no project says how to make one", {
  withr_home <- Sys.getenv("HOME")
  on.exit(Sys.setenv(HOME = withr_home))
  Sys.setenv(HOME = tempfile("pt-home-"))
  expect_error(pt_project(), "No PlotTrellis project in .*samples_txt = \"path/to/samples.txt\"")
})

test_that("a sample name in the wrong case gets a did-you-mean", {
  expect_error(.pt_resolve_tracks(tab, "atac@a", NULL), "Did you mean 'A'\\?")
})
