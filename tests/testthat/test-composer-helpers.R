# .pt_build_entries() (R/api-figure-helpers.R): track keys to render entries, on the mouse
# fixture's track items (samples Ba, LP, P6, P7), built as the startup context builds them
# (.pt_composer_context(), test-composer-context.R).

need_fixture_data <- function()
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")

# The mouse fixture's track items, as ctx$TRACK_ITEMS holds them.
mouse_track_items <- function() {
  rows <- load_samples_from_settings(FIXTURES$mouse$samples, FIXTURES$mouse$data_dir)
  c(build_track_items(rows), build_comparison_track_items(rows))
}

test_that("build_entries turns picker keys into render entries and drops unknown keys", {
  need_fixture_data()
  items <- mouse_track_items()
  f <- function(labels, col_map, settings) .pt_build_entries(labels, items, col_map, settings)
  e <- f(c("genes", "nope", "bw_atac_Ba", "genes"), c(Ba = "#E41A1C"), list())
  expect_length(e, 3)
  expect_identical(e[[1]], list(viz = "genes", height = NULL, data = NULL))
  expect_identical(e[[3]], e[[1]])
  expect_identical(e[[2]][c("viz", "sample", "colour", "name")],
                   list(viz = "bw", sample = "Ba", colour = "#E41A1C", name = "atac"))
  expect_identical(basename(e[[2]]$path), "Ba_ATAC.bw")
  expect_identical(f(character(0), list(), list()), list())
})

test_that("build_entries gives comparison tracks the role colours, defaulting when unset", {
  need_fixture_data()
  items <- mouse_track_items()
  f <- function(labels, col_map, settings) .pt_build_entries(labels, items, col_map, settings)
  key <- "bw_diff_stacked_atac_Ba_vs_LP"
  e <- f(key, list(), list(reference_colour = "#111111", comparison_colour = "#222222"))[[1]]
  expect_identical(e[c("viz", "control_label", "treatment_label", "control_colour", "treatment_colour")],
                   list(viz = "bw_diff_stacked", control_label = "Ba", treatment_label = "LP",
                        control_colour = "#111111", treatment_colour = "#222222"))
  e <- f(key, list(), list())[[1]]
  expect_identical(c(e$control_colour, e$treatment_colour), c("#4D4D4D", "#56B4E9"))
})
