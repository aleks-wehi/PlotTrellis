# Gene annotation ("genes") can be in a figure many times, every instance the same
# track. These check the figure API; the app's side is test-main-api.R's (labels).

test_that("pt_figure() resolves repeated Gene annotations in place", {
  tab <- data.frame(track = c("genes", "atac", "atac"), sample = c(NA, "A", "B"),
                    kind = c("fixed", "sample", "sample"), label = "",
                    key = c("genes", "atac_A", "atac_B"), stringsAsFactors = FALSE)
  expect_equal(.pt_resolve_tracks(tab, c("genes", "atac", "genes", "atac", "genes"), c("A", "B"),
                                  samples = c(NA, "A", NA, "B", NA)),
               c("genes", "atac_A", "genes", "atac_B", "genes"))
  # Repeated over every sample, each Gene annotation stays where it is written.
  expect_equal(.pt_resolve_tracks(tab, c("genes", "atac", "genes"), c("A", "B"), order = "sample"),
               c("genes", "atac_A", "atac_B", "genes"))
})

# The mammary example (Basal, LP), its project in a temp folder.
REPEAT_DIR <- file.path(tempdir(), "pt-test-repeat-genes")
repeat_project <- function() suppressWarnings(pt_project(REPEAT_DIR, pt_example("mammary")))

test_that("several Gene annotations reach the entries and the script", {
  proj <- repeat_project()
  fig <- pt_figure(proj, "Eya2", c("genes", "ATAC", "genes", "ATAC", "genes"), samples = c(NA, "Basal", NA, "LP", NA))
  expect_equal(fig$keys, c("genes", "bw_ATAC_Basal", "genes", "bw_ATAC_LP", "genes"))
  expect_equal(vapply(fig$entries, `[[`, character(1), "viz"), c("genes", "bw", "genes", "bw", "genes"))
  script <- readLines(pt_script(fig, tempfile(fileext = ".R")))
  expect_equal(sum(grepl("viz = \"genes\"", script, fixed = TRUE)), 3)
})
