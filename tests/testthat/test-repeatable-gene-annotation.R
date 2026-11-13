# Gene annotation ("genes") can be in a figure many times, every instance the same
# track; Available tracks always holds exactly one. The drag/double-click behaviour
# is client-side (checked in a headless browser, outside this suite); these check the
# server side: the pools the composer fills, the figure API and the start view.

test_that("pt_figure() resolves repeated Gene annotations in place", {
  tab <- data.frame(track = c("genes", "atac", "atac"), sample = c(NA, "A", "B"),
                    kind = c("fixed", "sample", "sample"), label = "",
                    key = c("genes", "atac_A", "atac_B"), stringsAsFactors = FALSE)
  expect_equal(.pt_resolve_tracks(tab, c("genes", "atac@A", "genes", "atac@B", "genes"), NULL),
               c("genes", "atac_A", "genes", "atac_B", "genes"))
})

test_that("a Multi start view keeps repeated Gene annotations and their order", {
  dir <- tempfile("pt-sv-"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  writeLines("{}", file.path(dir, "project_settings.json"))
  items <- list(genes = list(sample = NA), a_A = list(sample = "A"), a_B = list(sample = "B"),
                a_C = list(sample = "C"), cmp = list(sample = NA, samples = c("A", "B")))
  expect_warning(.pt_start_view(dir, items, c("genes", "a_A", "genes", "cmp", "a_B", "a_C", "genes"),
                                "X", c("A", "B", "C")), "no comparison tracks")
  sv <- jsonlite::read_json(file.path(dir, "project_settings.json"))$start_view
  expect_equal(sv$view_mode, "compose")
  expect_equal(unlist(sv$tracks), c("genes", "a_A", "genes", "a_B", "a_C", "genes"))
})

# The mammary example (Basal, LP), its project in a temp folder.
REPEAT_DIR <- file.path(tempdir(), "pt-test-repeat-genes")
repeat_project <- function() suppressWarnings(pt_project(REPEAT_DIR, pt_example("mammary")))

test_that("several Gene annotations reach the entries, the script and the start view", {
  proj <- repeat_project()
  fig <- pt_figure(proj, "Elf5", c("genes", "atac@Basal", "genes", "atac@LP", "genes"))
  expect_equal(fig$keys, c("genes", "bw_atac_Basal", "genes", "bw_atac_LP", "genes"))
  expect_equal(vapply(fig$entries, `[[`, character(1), "viz"), c("genes", "bw", "genes", "bw", "genes"))
  script <- readLines(pt_script(fig, tempfile(fileext = ".R")))
  expect_equal(sum(grepl("viz = \"genes\"", script, fixed = TRUE)), 3)
  json <- file.path(proj$dir, "project_settings.json")
  before <- readLines(json)
  on.exit(writeLines(before, json), add = TRUE)
  pt_preload(proj, fig)
  expect_equal(unlist(jsonlite::read_json(json)$start_view$tracks), fig$keys)
})

# What the composer server puts in Available/Figure tracks (update_rank_list() calls,
# by list), for each mode.
composer_lists <- function(json, inputs) {
  sent <- list()
  local_mocked_bindings(update_rank_list = function(css_id, text = NULL, labels = NULL, session = NULL)
    sent[[css_id]] <<- names(labels))
  app <- suppressMessages(.pt_composer_app(json, .pt_read_settings(json, .pt_default_settings())))
  shiny::testServer(app, do.call(session$setInputs, inputs))
  sent
}

test_that("Available tracks holds exactly one Gene annotation in every mode", {
  proj <- repeat_project()
  json <- file.path(proj$dir, "project_settings.json")
  before <- readLines(json)
  on.exit(writeLines(before, json), add = TRUE)
  # Without the example's preloaded start view: a new project's start.
  cfg <- jsonlite::read_json(json)
  cfg$start_view <- NULL
  writeLines(jsonlite::toJSON(cfg, auto_unbox = TRUE, pretty = TRUE, null = "null"), json)
  for (inputs in list(list(view_mode = "compose"), list(view_mode = "single", sel_condition = "Basal"),
                      list(view_mode = "two", sel_reference = "Basal", sel_comparison = "LP"))) {
    l <- composer_lists(json, inputs)
    expect_equal(sum(l$available_tracks == "genes"), 1, label = inputs$view_mode)
    expect_equal(l$figure_tracks, "genes", label = inputs$view_mode)
  }
})

test_that("a start view with several Gene annotations restores them, Available keeps one", {
  proj <- repeat_project()
  json <- file.path(proj$dir, "project_settings.json")
  before <- readLines(json)
  on.exit(writeLines(before, json), add = TRUE)
  pt_preload(proj, pt_figure(proj, "Elf5", c("genes", "atac@Basal", "genes", "atac@LP", "genes")))
  l <- composer_lists(json, list(view_mode = "two", sel_reference = "Basal", sel_comparison = "LP"))
  expect_equal(l$figure_tracks, c("genes", "bw_atac_Basal", "genes", "bw_atac_LP", "genes"))
  expect_equal(sum(l$available_tracks == "genes"), 1)
  expect_false(any(c("bw_atac_Basal", "bw_atac_LP") %in% l$available_tracks))
  # Only Gene annotation repeats: any other track repeated in the start view (a
  # pt_figure() can repeat one) is restored once.
  pt_preload(proj, pt_figure(proj, "Elf5", c("genes", "atac", "atac@Basal", "genes", "protein_domains",
                                             "protein_domains"), parallel = TRUE))
  l <- composer_lists(json, list(view_mode = "two", sel_reference = "Basal", sel_comparison = "LP"))
  expect_equal(l$figure_tracks, c("genes", "bw_atac_Basal", "bw_atac_LP", "genes", "protein_domains"))
  expect_false(any(c("bw_atac_Basal", "protein_domains") %in% l$available_tracks))
  # Multi (written directly: the example has 2 samples, so pt_preload() picks Two).
  .pt_write_start_view(proj$dir, "compose", character(0), "Elf5", c("bw_atac_LP", "genes", "bw_atac_Basal", "genes"))
  l <- composer_lists(json, list(view_mode = "compose"))
  expect_equal(l$figure_tracks, c("bw_atac_LP", "genes", "bw_atac_Basal", "genes"))
  expect_equal(sum(l$available_tracks == "genes"), 1)
  expect_false(any(c("bw_atac_Basal", "bw_atac_LP") %in% l$available_tracks))
})
