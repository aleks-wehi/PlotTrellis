# Bug 106: an inverted or empty region stops with a clear message, not "argument is of length
# zero" from deep inside a track.

test_that(".pt_region_error() accepts start < end and explains anything else", {
  expect_null(.pt_region_error(100, 200))
  expect_match(.pt_region_error(200, 100), "end \\(100\\) must be after its start \\(200\\)")
  expect_match(.pt_region_error(100, 100), "must be after its start")
  expect_match(.pt_region_error(NA, 100), "needs a start and an end")
  expect_match(.pt_region_error(integer(0), 100), "needs a start and an end")
  expect_match(.pt_region_error(100, 3e9), "must be below 2,147,483,647")
})

test_that("render_all() stops on an inverted or empty region before drawing anything", {
  for (r in list(list(start = 2000L, end = 1000L), list(start = 1000L, end = 1000L)))
    expect_error(render_all(list(), c(list(chrom = "chr1", gene = NA), r), list()), "must be after its start")
})

test_that("pt_figure() rejects an inverted region and a gene flank of -50% or less (example project)", {
  skip_if_not_installed("org.Mm.eg.db")
  proj <- suppressMessages(suppressWarnings(pt_project(tempfile("pt-region-"), pt_example("mammary"))))
  expect_error(pt_figure(proj, region = "chr2:103,500,000-103,400,000", tracks = "genes"),
               "end \\(103,400,000\\) must be after its start \\(103,500,000\\)")
  expect_error(pt_figure(proj, region = "chr2:103,400,000-103,400,000", tracks = "genes"),
               "must be after its start")
  expect_error(pt_figure(proj, region = "chr2:100-3,000,100,000", tracks = "genes"), "must be below")
  expect_s3_class(pt_figure(proj, region = "chr2:103,400,000-103,500,000", tracks = "genes"), "pt_figure")
  # Bug 265: two regions stop, instead of the first being used.
  expect_error(pt_figure(proj, region = c("chr2:103,400,000-103,500,000", "chr2:1-2"), tracks = "genes"),
               "region must be one string")
  expect_error(pt_figure(proj, region = 103400000, tracks = "genes"), "region must be one string")
  # The gene lookup (also the app's gene box) says why, instead of handing on a backwards region.
  res <- .pt_resolve_gene_region("Elf5", proj$txdb, -0.6, proj$genome$orgdb, proj$genes, proj$genome$genome_label)
  expect_null(res$region)
  expect_match(res$error, "flank of -60% leaves nothing of Elf5 to show: use a larger flank")
  expect_false(is.null(.pt_resolve_gene_region("Elf5", proj$txdb, -0.4, proj$genome$orgdb, proj$genes,
                                               proj$genome$genome_label)$region))
})
