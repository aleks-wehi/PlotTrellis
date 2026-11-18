# The composer's startup context (R/app-context.R). The settings checks stop before
# the genome loads (ensure_genome_ready() is replaced by a stop), so they need no genome;
# the full context is built from the bundled mammary example, as test-filter-buttons.R does.

# Settings for a context call that never gets as far as the genome.
pre_genome_cfg <- function(...) modifyList(list(gtf_path = "/no/such.gtf", orgdb = "org.Mm.eg.db",
                                                genome_label = "mm10", data_dir = "", samples = list()),
                                           list(...))

# Runs .pt_composer_context() up to ensure_genome_ready(), which stops with "reached the
# genome" after recording what it was given.
until_genome <- function(json, cfg) {
  seen <- NULL
  local_mocked_bindings(ensure_genome_ready = function(gtf_path, orgdb, genome_label, cache_dir) {
    seen <<- list(gtf_path = gtf_path, orgdb = orgdb, genome_label = genome_label, cache_dir = cache_dir)
    stop("reached the genome")
  })
  expect_error(.pt_composer_context(json, cfg), "reached the genome", fixed = TRUE)
  seen
}

# The bundled example's project and its context, built once for this file.
example_ctx <- local({
  memo <- NULL
  function() {
    if (is.null(memo)) {
      proj <- tempfile("pt-context-")
      suppressWarnings(suppressMessages(pt_project(proj, pt_example("mammary"))))
      json <- file.path(normalizePath(proj, winslash = "/"), "project_settings.json")
      cfg <- .pt_read_settings(json, .pt_default_settings())
      memo <<- list(json = json, cfg = cfg, ctx = suppressMessages(.pt_composer_context(json, cfg)))
    }
    memo
  }
})

test_that("the genome comes from the settings, with the caches in the project's cache folder", {
  dir <- tempfile("pt-ctx-"); dir.create(dir); on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  json <- file.path(dir, "project_settings.json")
  seen <- until_genome(json, pre_genome_cfg())
  expect_identical(seen, list(gtf_path = "/no/such.gtf", orgdb = "org.Mm.eg.db", genome_label = "mm10",
                              cache_dir = file.path(dir, "cache")))
})

test_that("an old genome_id that is a preset is migrated to flat fields and written back", {
  dir <- tempfile("pt-ctx-"); dir.create(dir); on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  json <- file.path(dir, "project_settings.json")
  cfg <- pre_genome_cfg(genome_id = "hg38", orgdb = NULL, genome_label = NULL)
  expect_message(seen <- until_genome(json, cfg),
                 "Migrated project_settings.json's old genome_id ('Human (GRCh38)') to flat genome fields.", fixed = TRUE)
  expect_identical(seen[c("orgdb", "genome_label")], list(orgdb = "org.Hs.eg.db", genome_label = "hg38"))
  written <- jsonlite::fromJSON(json, simplifyDataFrame = FALSE)
  expect_null(written$genome_id)
  expect_identical(written$orgdb, "org.Hs.eg.db")
  expect_identical(written$genome_label, "hg38")
  expect_identical(written$gtf_path, "/no/such.gtf")
})

test_that("an old genome_id that is not a preset stops with the exact message", {
  dir <- tempfile("pt-ctx-"); dir.create(dir); on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  json <- file.path(dir, "project_settings.json")
  expect_error(.pt_composer_context(json, pre_genome_cfg(genome_id = "zz9", orgdb = NULL)),
               paste0("project_settings.json has an old-format genome_id ('zz9') that isn't a known built-in preset, ",
                      "so it can't be auto-migrated. Re-open Settings and fill in the Genome section (GTF path, OrgDb, ",
                      "genome label) by hand."), fixed = TRUE)
  expect_false(file.exists(json))
})

test_that("a genome_id beside an orgdb is left alone", {
  dir <- tempfile("pt-ctx-"); dir.create(dir); on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  json <- file.path(dir, "project_settings.json")
  expect_no_message(seen <- until_genome(json, pre_genome_cfg(genome_id = "hg38")))
  expect_identical(seen$orgdb, "org.Mm.eg.db")
  expect_false(file.exists(json))
})

test_that("a blank genome label stops before the genome loads", {
  dir <- tempfile("pt-ctx-"); dir.create(dir); on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  json <- file.path(dir, "project_settings.json")
  for (lbl in list(NULL, "", "  ")) {
    local_mocked_bindings(ensure_genome_ready = function(...) stop("reached the genome"))
    expect_error(.pt_composer_context(json, pre_genome_cfg(genome_label = lbl)),
                 paste0(PT_GENOME_LABEL_REQUIRED, ' Set "genome_label" in ', json), fixed = TRUE)
  }
})

test_that("the example's context holds every startup field, built from its settings", {
  ex <- example_ctx(); ctx <- ex$ctx
  expect_named(ctx, c("SETTINGS_JSON", "CFG", "PROJECT_DIR", "CACHE_DIR", "GENOME", "ORGDB_PKG", "txdb", "CDS_BY_TX",
                      "CHR_LENGTHS", "SAMPLES_DF", "AVAILABLE_CONDITIONS", "TRACK_ITEMS", "ALL_GENES", "RENDER_CFG",
                      "DEFAULT_SAMPLE_COLS", "START_VIEW", "COMPOSE_AVAILABLE_KEYS", "GENE_REQUIRED_TYPES"))
  expect_identical(ctx$SETTINGS_JSON, ex$json)
  expect_identical(ctx$PROJECT_DIR, dirname(ex$json))
  expect_identical(ctx$CACHE_DIR, file.path(dirname(ex$json), "cache"))
  expect_identical(ctx$ORGDB_PKG, ex$cfg$orgdb)
  expect_identical(ctx$RENDER_CFG$orgdb_pkg, ctx$ORGDB_PKG)
  expect_identical(ctx$RENDER_CFG$genome_label, ex$cfg$genome_label)
  expect_s4_class(ctx$txdb, "TxDb")
})

test_that("the example's conditions, colours and track items agree", {
  ctx <- example_ctx()$ctx
  expect_identical(ctx$AVAILABLE_CONDITIONS, vapply(ctx$SAMPLES_DF, `[[`, "", "display"))
  expect_true(length(ctx$AVAILABLE_CONDITIONS) >= 2)
  expect_identical(names(ctx$DEFAULT_SAMPLE_COLS), ctx$AVAILABLE_CONDITIONS)
  expect_identical(ctx$TRACK_ITEMS, c(build_track_items(ctx$SAMPLES_DF), build_comparison_track_items(ctx$SAMPLES_DF)))
  expect_true("genes" %in% names(ctx$TRACK_ITEMS))
  expect_identical(ctx$COMPOSE_AVAILABLE_KEYS, .pt_compose_available_keys(ctx$TRACK_ITEMS))
  expect_identical(ctx$GENE_REQUIRED_TYPES, .pt_gene_required_types())
})

test_that("the example's chromosome lengths and gene list come from the genome", {
  ctx <- example_ctx()$ctx
  expect_true(is.numeric(ctx$CHR_LENGTHS))
  expect_false(is.null(names(ctx$CHR_LENGTHS)))
  expect_true(all(ctx$CHR_LENGTHS > 0))
  expect_identical(ctx$CHR_LENGTHS, readRDS(ctx$GENOME$chrom_sizes_path))
  expect_type(ctx$ALL_GENES, "character")
  expect_false(is.unsorted(ctx$ALL_GENES))
  expect_true(length(ctx$ALL_GENES) > 1000)
})

test_that("the start view keeps known tracks in order, repeats of Gene annotation only", {
  ex <- example_ctx()
  keys <- setdiff(names(ex$ctx$TRACK_ITEMS), "genes")[1:2]
  cfg <- ex$cfg
  cfg$start_view <- list(view_mode = list("compose"),
                         tracks = list("genes", keys[1], "no_such_track", keys[1], "genes", keys[2], "sashimi_gone"))
  ctx <- suppressMessages(.pt_composer_context(ex$json, cfg))
  # Lists are unlisted; unknown keys (and an old junctions key still unknown once renamed)
  # dropped; a repeated track kept once, a repeated Gene annotation kept.
  expect_identical(ctx$START_VIEW$view_mode, "compose")
  expect_identical(ctx$START_VIEW$tracks, c("genes", keys[1], "genes", keys[2]))
})

test_that("no start view, or one that is not a list, gives NULL", {
  ex <- example_ctx()
  cfg <- ex$cfg
  cfg$start_view <- NULL
  expect_null(suppressMessages(.pt_composer_context(ex$json, cfg))$START_VIEW)
  cfg$start_view <- "compose"
  expect_null(suppressMessages(.pt_composer_context(ex$json, cfg))$START_VIEW)
})
