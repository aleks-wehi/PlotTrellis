# The junction track's viz key was once "sashimi". Projects, track names and figures
# saved with that key still open.

test_that("old junction keys map to the new ones, and nothing else changes", {
  expect_equal(.pt_renamed_junctions_key(c("sashimi", "sashimi_junctions_GM12878", "genes", "bw_atac_K562",
                                           "junctions_junctions_K562", "junctions_mirrored_diff_junctions_GM12878_vs_K562")),
               c("junctions", "junctions_junctions_GM12878", "genes", "bw_atac_K562",
                 "junctions_junctions_K562", "junctions_mirrored_diff_junctions_GM12878_vs_K562"))
  expect_null(.pt_renamed_junctions_key(NULL))
})

test_that("a track name with the old _sashimi suffix still resolves", {
  tab <- data.frame(track = c("rnaseq", "rnaseq_junctions"), sample = c("A", "A"), kind = "sample",
                    label = "", key = c("bw_rnaseq_A", "junctions_rnaseq_A"), stringsAsFactors = FALSE)
  expect_equal(.pt_resolve_tracks(tab, c("rnaseq", "rnaseq_sashimi", "rnaseq_junctions"), "A"),
               c("bw_rnaseq_A", "junctions_rnaseq_A", "junctions_rnaseq_A"))
  expect_error(.pt_resolve_tracks(tab, "atac_sashimi", "A"), "Unknown track 'atac_sashimi'")
})

test_that("a pt_figure saved with the old viz key renders and exports as the new one", {
  skip_if_not_installed("org.Hs.eg.db")
  proj <- suppressWarnings(pt_project(file.path(tempdir(), "pt-test-renamed-key"), pt_example("human")))
  fig <- pt_figure(proj, "PSME2", c("genes", "junctions"), samples = "GM12878")
  expect_equal(fig$entries[[2]]$viz, "junctions")
  old <- fig
  old$entries[[2]]$viz <- "sashimi"
  svgs <- c(pt_save(fig, tempfile(fileext = ".svg")), pt_save(old, tempfile(fileext = ".svg")))
  expect_identical(readLines(svgs[2]), readLines(svgs[1]))
  expect_false(any(grepl("track failed", readLines(svgs[2]), fixed = TRUE)))
  scripts <- c(pt_script(fig, tempfile(fileext = ".R")), pt_script(old, tempfile(fileext = ".R")))
  expect_identical(readLines(scripts[2]), readLines(scripts[1]))
})

# A project as the version before the rename made it: junction items typed and keyed "sashimi".
old_project <- function(proj) {
  old <- function(x) sub("^junctions(_|$)", "sashimi\\1", x)
  is_jn <- vapply(proj$items, function(it) identical(it$type, "junctions"), logical(1))
  proj$items[is_jn] <- lapply(proj$items[is_jn], function(it) modifyList(it, list(key = old(it$key), type = "sashimi")))
  names(proj$items)[is_jn] <- old(names(proj$items)[is_jn])
  proj$tracks$key[proj$tracks$key %in% names(is_jn)[is_jn]] <- old(proj$tracks$key[proj$tracks$key %in% names(is_jn)[is_jn]])
  proj
}

test_that("a pt_project kept from before the rename draws its junction tracks", {
  skip_if_not_installed("org.Hs.eg.db")
  proj <- suppressWarnings(pt_project(file.path(tempdir(), "pt-test-renamed-key"), pt_example("human")))
  old <- old_project(proj)
  expect_true("sashimi_junctions_GM12878" %in% names(old$items))
  region <- "chr14:24,140,000-24,150,000"
  fig <- pt_figure(proj, region = region, tracks = "junctions", samples = "GM12878")
  old_fig <- pt_figure(old, region = region, tracks = "junctions", samples = "GM12878")
  expect_identical(old_fig$entries, fig$entries)
  expect_false(is.null(old_fig$entries[[1]]$path))
  expect_identical(old_fig$keys, fig$keys)
})

