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
S <- c("A", "B", "C")   # samples.txt order
res <- function(tracks, ...) .pt_resolve_tracks(tab, tracks, S, ...)

test_that("with neither samples nor comparison, each per-sample track is drawn for every sample", {
  expect_message(keys <- res(c("genes", "atac", "junctions")), "no 'junctions' for B")
  expect_equal(keys, c("genes", "atac_A", "atac_B", "atac_C", "junc_A"))
})

test_that("order = \"sample\" groups a run of per-sample tracks by sample; fixed tracks stay put", {
  keys <- suppressMessages(res(c("genes", "atac", "junctions", "transcripts"), order = "sample"))
  expect_equal(keys, c("genes", "atac_A", "junc_A", "atac_B", "atac_C", "transcripts"))
  expect_equal(suppressMessages(res(c("atac", "junctions"), order = "track")), c("atac_A", "atac_B", "atac_C", "junc_A"))
  expect_error(res("atac", order = "samples"), "order must be \"track\".*or \"sample\"")
})

test_that("samples of length 1 picks that sample for every per-sample track", {
  expect_equal(res(c("atac", "genes", "junctions"), samples = "A"), c("atac_A", "genes", "junc_A"))
  expect_error(res(c("atac", "junctions"), samples = "B"), "No 'junctions' track for 'B'. Samples that have it: A")
  expect_error(res("atac", samples = NA), "Track 1 \\('atac'\\) needs a sample, but samples\\[1\\] is NA")
})

test_that("samples as long as tracks gives track i sample i, NA for fixed tracks", {
  expect_equal(res(c("genes", "atac", "atac", "junctions"), samples = c(NA, "C", "A", "A")),
               c("genes", "atac_C", "atac_A", "junc_A"))
  expect_error(res(c("genes", "atac"), samples = c("A", "B")), "samples\\[1\\] is 'A', but track 1 \\('genes'\\) has no sample: put NA there")
  expect_error(res(c("genes", "atac"), samples = c(NA, NA)), "Track 2 \\('atac'\\) needs a sample, but samples\\[2\\] is NA")
})

test_that("any other length of samples is an error saying what lengths work", {
  expect_error(res(c("genes", "atac", "atac"), samples = c("A", "B")), "samples has 2 entries for 3 tracks: give 1 .* or 3")
  expect_error(res("atac", samples = c("A", "B")), "samples has 2 entries for 1 track: give 1")
})

test_that("comparison repeats per-sample tracks over its samples and gives comparison tracks their 2", {
  # comparison = c("C", "A") is C against the reference A ("C vs A"); per-sample tracks in that order.
  expect_equal(res(c("atac", "junctions_diff", "genes"), comparison = c("C", "A")), c("atac_C", "atac_A", "diff_A_vs_C", "genes"))
  expect_equal(res("atac", comparison = c("C", "A")), c("atac_C", "atac_A"))
  expect_equal(suppressMessages(res(c("junctions", "atac", "junctions_diff"), comparison = c("B", "A"), order = "sample")),
               c("atac_B", "junc_A", "atac_A", "diff_A_vs_B"))
  expect_error(res("atac", comparison = c("A", "A")), "comparison lists 'A' twice")
  expect_error(res("atac", comparison = c("A", NA)), "comparison can't contain NA")
})

test_that("a comparison track needs comparison with exactly 2 samples, the second the reference", {
  expect_error(res("junctions_diff"), "'junctions_diff' compares 2 samples: set comparison = c\\(\"B\", \"A\"\\) \\(reference second\\)")
  expect_error(res("junctions_diff", comparison = c("A", "B", "C")), "comparison needs exactly 2 \\(reference second\\); it has 3")
  expect_error(res(c("atac", "junctions_diff"), samples = "A"), "which samples = can't give: drop samples and set comparison = c\\(\"B\", \"A\"\\) \\(reference second\\)")
  expect_error(res("junctions_diff", comparison = c("B", "C")), "No 'junctions_diff' track for B vs C. Pairs that have it: B vs A, C vs A")
  # The pair the heading writes, "B vs A", is comparison = c("B", "A"); the other way round is another pair.
  expect_equal(res("junctions_diff", comparison = c("B", "A")), "diff_A_vs_B")
  expect_error(res("junctions_diff", comparison = c("A", "B")), "No 'junctions_diff' track for A vs B")
})

test_that("samples and comparison together are an error", {
  expect_error(res("atac", samples = "A", comparison = c("A", "B")), "Set samples or comparison, not both")
})

test_that("sample names match exactly; a near miss lists the names with their numbers", {
  expect_error(res("atac", samples = "a"), "Unknown sample 'a' \\(in samples\\).*exactly as in samples.txt.*1 A, 2 B, 3 C")
  expect_error(res("atac", comparison = c("A", "Bb")), "Unknown sample 'Bb' \\(in comparison\\)")
  expect_error(res("atac", samples = TRUE), "samples must be sample names or numbers.*parallel = TRUE is gone")
})

test_that("tracks and samples also take numbers (pt_tracks() ids, samples.txt rows)", {
  expect_equal(res(c(3, 1), samples = 2), c("atac_B", "genes"))
  expect_equal(res(c(1, 3, 3), samples = c(NA, 3, 1)), c("genes", "atac_C", "atac_A"))
  expect_equal(res(c(3, 5), comparison = c(3, 1)), c("atac_C", "atac_A", "diff_A_vs_C"))
  expect_error(res(9), "There is no track number 9. This project's tracks: 1 genes, 2 transcripts, 3 atac, 4 junctions, 5 junctions_diff")
  expect_error(res(1.5), "There is no track number 1.5")
  expect_error(res("atac", samples = 4), "There is no sample number 4 \\(in samples\\).*1 A, 2 B, 3 C")
})

test_that("the old name@sample form and unknown names stop with a message saying what to write", {
  # A suggested call is one that draws exactly the same tracks; with none, it says how.
  e <- expect_error(res(c("atac@B", "atac@A", "genes")), fixed = TRUE,
                    'The same figure: tracks = c("atac", "atac", "genes"), samples = c("B", "A", NA).')
  expect_false(grepl("start_", conditionMessage(e), fixed = TRUE))   # no samples.txt start lines any more
  expect_error(res(c("genes", "atac@B", "atac@A", "junctions_diff@B vs A")), fixed = TRUE,
               'The same figure: tracks = c("genes", "atac", "junctions_diff"), comparison = c("B", "A").')
  tab_jb <- rbind(tab, data.frame(track = "junctions", sample = "B", kind = "sample", label = "", key = "junc_B"))
  expect_error(.pt_resolve_tracks(tab_jb, c("genes", "atac@B", "junctions@B", "atac@A", "junctions@A", "junctions_diff@B vs A"), S),
               fixed = TRUE, 'comparison = c("B", "A"), order = "sample".')
  # Per-sample tracks reference first, then the comparison: no comparison = call draws that order.
  expect_error(res(c("genes", "atac@A", "atac@B", "junctions_diff@B vs A")), fixed = TRUE,
               'comparison = c(other, reference) for comparison tracks')
  expect_error(res(c("atac@C", "junctions_diff@B vs A")), "form is gone.*one per track.*See \\?pt_figure")
  expect_error(res(c("atac", "junctions@A")), "form is gone.*See \\?pt_figure")
  expect_error(res("atac@D"), "form is gone.*See \\?pt_figure")
  expect_error(res("ATAC"), "Unknown track 'ATAC'.*3 atac")
  expect_error(res(NA_character_), "tracks must be track names or numbers")
  expect_error(res(character(0)), "tracks is empty")
  expect_error(res(numeric(0)), "tracks is empty")
  # Numbers in quotes are names; the error says so.
  expect_error(res("1"), "'1' is in quotes, so it is read as a name: write 1 without quotes")
  expect_error(res("atac", samples = "2"), "'2' is in quotes, so it is read as a name")
  # Only a TRUE/FALSE samples mentions the old parallel.
  expect_error(res("atac", samples = factor("A")), "must be sample names or numbers, e.g. samples = \"A\" or samples = 1. This")
})

test_that("the name@sample hint quotes sample names with quotes or backslashes safely", {
  tabq <- data.frame(track = c("genes", "atac", "atac"), sample = c(NA, 'A"x', 'B\\q'),
                     kind = c("fixed", "sample", "sample"), label = "",
                     key = c("genes", "atac_1", "atac_2"), stringsAsFactors = FALSE)
  smp <- c('A"x', 'B\\q', 'x", stop("ran"), "')
  expect_error(.pt_resolve_tracks(tabq, c('atac@A"x', 'atac@B\\q', "genes"), smp), fixed = TRUE,
               'samples = c("A\\"x", "B\\\\q", NA)')
  expect_error(.pt_resolve_tracks(tabq, 'atac@x", stop("ran"), "', smp), "form is gone")
})

test_that("transcript_counts, the old name of the quantification track, finds it even when another type took the new name", {
  tab2 <- data.frame(track = c("transcript_expression", "transcript_expression_tpm_usage"),
                     sample = c("A", "A"), kind = "sample", label = "",
                     key = c("bw_transcript_expression_A", "tpm_usage_A"), stringsAsFactors = FALSE)
  expect_equal(.pt_resolve_tracks(tab2, "transcript_counts", "A"), "tpm_usage_A")
  expect_equal(.pt_resolve_tracks(tab2, "transcript_expression", "A"), "bw_transcript_expression_A")
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
  fig <- pt_figure(proj, "psme2", c("genes", "atac", "junctions_diff"), comparison = c("K562", "GM12878"))
  expect_equal(fig$gene, "PSME2")
  expect_equal(fig$keys, c("genes", "bw_atac_K562", "bw_atac_GM12878", "junctions_mirrored_diff_junctions_GM12878_vs_K562"))
  png <- pt_save(fig, tempfile(fileext = ".png"), dpi = 50)
  # rtracklayer and plotgardener cannot read bigWig on Windows, so the
  # atac tracks draw empty there and the PNG is small.
  if (.Platform$OS.type != "windows") expect_gt(file.size(png), 10000)
  expect_true(file.exists(png))
  # The project's OrgDb is used through its namespace, never attached to the search path.
  expect_false(any(grepl("^package:org\\.", setdiff(search(), attached))))
  expect_gt(length(readLines(pt_script(fig, tempfile(fileext = ".R")))), 100)
})

test_that("pt_figure() says what to write instead of parallel, and rejects unknown arguments", {
  skip_if_not_installed("org.Hs.eg.db")
  proj <- suppressWarnings(pt_project(EX_DIR, pt_example("human")))
  expect_error(pt_figure(proj, "PSME2", "atac", parallel = TRUE), "parallel is gone: leave samples and comparison unset")
  expect_error(pt_figure(proj, "PSME2", "atac", colour = "red"), "no argument 'colour'")
  expect_error(pt_figure(proj, "PSME2", "atac", NULL, NULL, NULL, "", "track", "extra"), "takes 8 arguments; this call gives more")
})

test_that("transcript_counts, the old name of transcript_expression, still draws the same figure", {
  skip_if_not_installed("org.Hs.eg.db")
  proj <- suppressWarnings(pt_project(EX_DIR, pt_example("human")))
  tr <- pt_tracks(proj)
  expect_true("transcript_expression" %in% tr$track)
  expect_false("transcript_counts" %in% tr$track)           # listed under the new name only
  expect_match(tr$label[tr$track == "transcript_expression"][1], "^Transcript expression")
  new <- pt_figure(proj, "PSME2", c("genes", "transcript_expression"), samples = "GM12878")
  old <- expect_silent(pt_figure(proj, "PSME2", c("genes", "transcript_counts"), samples = "GM12878"))   # no warning or message
  expect_identical(old$keys, new$keys)
  expect_identical(old$entries, new$entries)
  expect_identical(pt_figure(proj, "PSME2", "transcript_counts")$keys,
                   pt_figure(proj, "PSME2", "transcript_expression")$keys)
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

test_that("pt_tracks(), pt_samples() and print(proj) give the numbers pt_figure() takes", {
  skip_if_not_installed("org.Hs.eg.db")
  proj <- suppressWarnings(pt_project(EX_DIR, pt_example("human")))
  tr <- pt_tracks(proj)
  expect_named(tr, c("id", "track", "kind", "label", "samples"))
  expect_equal(tr$id, seq_len(nrow(tr)))
  expect_false(anyDuplicated(tr$track) > 0)                    # one row per track name
  expect_equal(tr$samples[tr$track == "atac"], "GM12878, K562")
  expect_equal(tr$samples[tr$track == "genes"], "")
  expect_equal(tr$label[tr$track == "atac"], "atac")           # the label without its sample
  smp <- pt_samples(proj)
  expect_equal(smp$name, c("GM12878", "K562"))
  expect_equal(smp$id, 1:2)
  expect_match(smp$colour, "^#")
  out <- capture.output(print(proj))
  expect_true(any(grepl("^ +1 GM12878$", out)) && any(grepl("^ +2 K562$", out)))
  expect_true(any(grepl(sprintf("^ +%d atac +GM12878, K562$", tr$id[tr$track == "atac"]), out)))
  # Numbers and names give the same figure.
  ids <- pt_figure(proj, "PSME2", tr$id[match(c("genes", "atac"), tr$track)], samples = 2)
  expect_identical(ids$keys, pt_figure(proj, "PSME2", c("genes", "atac"), samples = "K562")$keys)
  expect_equal(ids$tracks, c("Gene annotation", "atac \u00b7 K562"))   # print(fig) shows the names
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
  data <- dirname(example_copy("human"))
  # Caches an older version left beside the GTF are not reused: a TxDb there that isn't one
  # would break the project if it were read.
  writeLines("not a TxDb", file.path(data, "annotation_txdb.sqlite"))
  before <- tools::md5sum(list.files(data, recursive = TRUE, full.names = TRUE))
  local_offline_chrom_sizes()
  proj_dir <- tempfile("pt-proj-")
  proj <- suppressWarnings(pt_project(proj_dir, file.path(data, "samples.txt")))
  # The data folder is untouched, file for file.
  expect_identical(tools::md5sum(list.files(data, recursive = TRUE, full.names = TRUE)), before)
  genome_files <- unlist(proj$genome[c("txdb_path", "tx_name_lookup_path", "tx_meta_lookup_path", "chrom_sizes_path")])
  expect_true(all(file.exists(genome_files)))
  expect_true(all(startsWith(normalizePath(genome_files), normalizePath(file.path(proj_dir, "cache", "genome")))))
  expect_identical(basename(proj$genome$txdb_path), "annotation_txdb.sqlite")
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

test_that("a project made from a bundled example reopens after PlotTrellis moves (bug 263)", {
  skip_if_not_installed("org.Hs.eg.db")
  dir <- tempfile("pt-moved-")
  suppressWarnings(pt_project(dir, pt_example("human")))
  json <- file.path(dir, "project_settings.json")
  made_from <- function(src) {
    s <- jsonlite::read_json(json)
    s$samples_txt_source <- src
    s$data_dir <- "/old/R/4.4/library/PlotTrellis/extdata/example-human"
    writeLines(jsonlite::toJSON(s, auto_unbox = TRUE, pretty = TRUE, null = "null", digits = NA), json)
  }
  # As the same example in another library recorded it (an R upgrade or reinstall).
  made_from("/old/R/4.4/library/PlotTrellis/extdata/example-human/samples.txt")
  expect_message(proj <- suppressWarnings(pt_project(dir, pt_example("human"))), "Re-importing")
  s <- jsonlite::read_json(json)
  expect_identical(s$samples_txt_source, normalizePath(pt_example("human")))
  # Its data is read from this install now, not the old one.
  expect_identical(normalizePath(s$data_dir), normalizePath(dirname(pt_example("human"))))
  expect_s3_class(pt_figure(proj, "PSME2", "genes"), "pt_figure")
  # A source tree of the package (a clone of its repository) too.
  made_from("/home/me/src/PlotTrellis/inst/extdata/example-human/samples.txt")
  expect_message(suppressWarnings(pt_project(dir, pt_example("human"))), "Re-importing")
  # Opened again from the same place: as before, no re-import.
  expect_no_message(suppressWarnings(pt_project(dir, pt_example("human"))), message = "Re-importing")
  # The other example, or a samples.txt outside a PlotTrellis install, is still another project.
  made_from("/old/R/4.4/library/PlotTrellis/extdata/example/samples.txt")
  expect_error(pt_project(dir, pt_example("human")), "already holds a project")
  made_from("/home/me/extdata/example-human/samples.txt")
  expect_error(pt_project(dir, pt_example("human")), "already holds a project")
})

test_that("bare pt_project() with no project says how to make one", {
  withr_home <- Sys.getenv("HOME")
  on.exit(Sys.setenv(HOME = withr_home))
  Sys.setenv(HOME = tempfile("pt-home-"))
  expect_error(pt_project(), "No PlotTrellis project in .*samples_txt = \"path/to/samples.txt\"")
})

test_that("a sample name in the wrong case is not matched; the error lists names and numbers", {
  expect_error(.pt_resolve_tracks(tab, "atac", S, samples = "a"), "Unknown sample 'a'.*1 A, 2 B, 3 C")
})

# A human-example samples.txt in a new temp folder, its "# start_" lines replaced by `start`
# (lines put first); returns it and the project folder beside it.
tab_samples_txt <- function(start = character(0)) {
  st <- file.path(tempfile(), "samples.txt")
  dir.create(dirname(st))
  txt <- readLines(pt_example("human"))
  txt <- sub("^# data_root: .*", paste("# data_root:", dirname(pt_example("human"))), txt)
  txt <- sub("^# gtf_path: .*", paste("# gtf_path:", file.path(dirname(pt_example("human")), "annotation.gtf.gz")), txt)
  writeLines(c(start, txt[!grepl("^# start_", txt)]), st)
  list(st = st, dir = file.path(dirname(st), "proj"))
}

test_that("a project kept with saveRDS() and read back still finds genes (bug 122)", {
  skip_if_not_installed("org.Hs.eg.db")
  proj <- suppressWarnings(pt_project(EX_DIR, pt_example("human")))
  f <- tempfile(fileext = ".rds")
  saveRDS(proj, f)
  # readRDS() gives a TxDb with a dead database connection; with the open ones closed and
  # forgotten, as in a new session, pt_figure() has to reopen it.
  .pt_close_txdbs()
  expect_equal(pt_figure(readRDS(f), "PSME2", "genes")$gene, "PSME2")
})

test_that("region = takes the chromosome with or without chr (bug 202)", {
  skip_if_not_installed("org.Hs.eg.db")
  proj <- suppressWarnings(pt_project(EX_DIR, pt_example("human")))
  chrom <- function(r) pt_figure(proj, region = r, tracks = "genes")$region_params$chrom
  expect_equal(chrom("14:24,140,000-24,160,000"), "chr14")
  expect_equal(chrom("chr14:24,140,000-24,160,000"), "chr14")
  expect_equal(chrom("zz9:1-1000"), "zz9")   # no such chromosome: as typed
})

test_that("knit_print() of a figure outside knitr writes nothing in the working directory (bug 131)", {
  skip_if_not_installed("org.Hs.eg.db")
  skip_if_not_installed("knitr")
  proj <- suppressWarnings(pt_project(EX_DIR, pt_example("human")))
  wd <- tempfile("pt-wd-")
  dir.create(wd)
  old <- setwd(wd)
  on.exit(setwd(old))
  out <- knitr::knit_print(pt_figure(proj, "PSME2", "genes"))
  expect_match(as.character(out), "\\.png")
  expect_length(list.files(wd, recursive = TRUE, all.files = TRUE, no.. = TRUE), 0)
})

test_that("a samples.txt edited before Save as project default is still re-imported (bug 220)", {
  skip_if_not_installed("org.Hs.eg.db")
  p <- tab_samples_txt()
  Sys.setFileTime(p$st, Sys.time() - 120)
  suppressWarnings(pt_project(p$dir, p$st))
  json <- file.path(p$dir, "project_settings.json")
  # Edited after the import, then the app's Save as project default rewrites the settings,
  # which are now newer than the edit.
  writeLines(sub("#377EB8", "#00AA00", readLines(p$st), fixed = TRUE), p$st)
  Sys.setFileTime(p$st, Sys.time() - 60)
  .pt_write_project_files(.pt_read_settings(json, .pt_default_settings()), json)
  expect_message(proj <- pt_project(p$dir, p$st), "Re-importing")
  expect_equal(proj$colours[["K562"]], "#00AA00")
  # Not edited since: opened as is.
  stamp <- file.mtime(json)
  pt_project(p$dir, p$st)
  expect_equal(file.mtime(json), stamp)
})

test_that("a samples.txt's '# start_' lines are ignored, with one message", {
  skip_if_not_installed("org.Hs.eg.db")
  # Lines the app once opened on, one of them naming a sample the project doesn't have.
  p <- tab_samples_txt(c("# start_gene: PSME2", "# start_tracks: genes, atac", "# start_samples: NA, Zzz",
                         "# start_comparison: K562, GM12878", "# start_order: sample"))
  msgs <- character(0)
  warns <- character(0)
  withCallingHandlers(pt_project(p$dir, p$st),
                      message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") },
                      warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
  ignoring <- grep("Ignoring the '# start_' lines of", msgs, fixed = TRUE, value = TRUE)
  expect_length(ignoring, 1)
  expect_true(startsWith(ignoring, "Ignoring the '# start_' lines of 'samples.txt'"))
  expect_length(grep("start", warns), 0)
  settings <- jsonlite::read_json(file.path(p$dir, "project_settings.json"))
  expect_false("start_view" %in% names(settings))
  expect_false(any(grepl("^# start_", readLines(file.path(p$dir, "project_samples.txt")))))
  expect_s3_class(pt_figure(pt_project(p$dir), "PSME2", "genes"), "pt_figure")
})
