# The paper's bundled example (pt_example("paper"), inst/extdata/example-paper/): five
# genes, every data type. Uses only the installed example data, so these run on
# R-universe too. Needs org.Mm.eg.db (Suggests).

PAPER_GENES <- c("Eya2", "Lrrfip1", "Acot7", "Dzip1", "Cxcl14")
PAPER_DIR <- file.path(tempdir(), "pt-test-example-paper")

paper_project <- function() {
  # txdbmaker's two notices for any GENCODE GTF are expected; anything else still warns.
  withCallingHandlers(
    pt_project(PAPER_DIR, pt_example("paper")),
    warning = function(w) {
      if (grepl("\"phase\" metadata column|genome version information", conditionMessage(w)))
        invokeRestart("muffleWarning")
    })
}

test_that("pt_example(\"paper\") gives the installed paper example, small enough to ship", {
  txt <- pt_example("paper")
  expect_same_path(dirname(txt), system.file("extdata", "example-paper", package = "PlotTrellis"))
  files <- list.files(dirname(txt), recursive = TRUE, full.names = TRUE)
  # Under 5 MB in all and 1 MB per file: R-universe builds the package with it.
  expect_lt(sum(file.size(files)), 5e6)
  expect_true(all(file.size(files) < 1e6), info = basename(files[file.size(files) >= 1e6]))
  res <- import_samples_txt(txt, defaults = .pt_default_settings())
  expect_length(res$errors, 0)
  expect_length(res$warnings, 0)
})

test_that("the paper example opens offline with every data type and its start view", {
  skip_if_not_installed("org.Mm.eg.db")
  # The shipped chromosome sizes must be used, not fetched from UCSC.
  local_mocked_bindings(build_chrom_sizes = function(...) stop("tried to download chromosome sizes"))
  proj <- paper_project()
  expect_equal(pt_samples(proj)$name, c("Basal", "LP"))
  expect_setequal(unique(pt_tracks(proj)$track),
                  c("genes", "protein_domains", "transcripts_all", "transcripts", "transcripts_annotated",
                    "RNA", "junctions", "transcript_expression", "ATAC", "ATAC_peaks", "H3K36me3",
                    "H3K36me3_peaks", "IgG", "HiC", "RNA_diff", "ATAC_diff", "H3K36me3_diff", "IgG_diff",
                    "junctions_diff", "junctions_diff_tpm"))
  expect_equal(.pt_project(proj$dir)$cfg$start_view$gene, PAPER_GENES[1])
  # Each gene opens by name, on its window's chromosome.
  chroms <- c(Eya2 = "chr2", Lrrfip1 = "chr1", Acot7 = "chr4", Dzip1 = "chr14", Cxcl14 = "chr13")
  for (g in PAPER_GENES)
    expect_equal(pt_figure(proj, g, "genes")$region_params$chrom, chroms[[g]], info = g)
})

# Every track type at each gene, both samples and the comparison: no track may fail (a
# failed track is drawn as a message in the figure and said on the console). Protein
# domains come from UniProt over the internet, so they are tried below only when online.
# About 2 minutes, so not on CRAN-mode builds (R-universe); bigWig can't be read on Windows.
test_that("each of the paper's genes draws every track type with no error", {
  skip_if_not_installed("org.Mm.eg.db")
  skip_on_cran()
  skip_on_os("windows")
  proj <- paper_project()
  tracks <- setdiff(unique(pt_tracks(proj)$track), "protein_domains")
  for (g in PAPER_GENES) {
    fig <- pt_figure(proj, g, tracks, comparison = c("Basal", "LP"))
    msgs <- character(0)
    # plotgardener's "Not enough data within range to plot." (a sparse track) and "No
    # range data to plot." (a peaks track with no peak in view) are said with the full
    # data too; every other warning is still reported.
    svg <- withCallingHandlers(pt_save(fig, tempfile(fileext = ".svg")),
                               message = function(m) {
                                 msgs <<- c(msgs, conditionMessage(m))
                                 invokeRestart("muffleMessage")
                               },
                               warning = function(w) {
                                 if (grepl("Not enough data within range|No range data to plot", conditionMessage(w)))
                                   invokeRestart("muffleWarning")
                               })
    expect_identical(grep("track failed", msgs, value = TRUE), character(0), info = g)
    # A crop must still name each window's chromosome in every file: the full data
    # draws an empty panel there, not this message.
    expect_false(any(grepl("Chromosome not found", readLines(svg, warn = FALSE), fixed = TRUE)), info = g)
  }
})

test_that("the paper's genes draw their protein domains (online)", {
  skip_if_not_installed("org.Mm.eg.db")
  skip_on_cran()
  skip_if_offline()
  proj <- paper_project()
  for (g in PAPER_GENES) {
    msgs <- character(0)
    withCallingHandlers(pt_save(pt_figure(proj, g, c("protein_domains", "genes")),
                                tempfile(fileext = ".png"), dpi = 30),
                        message = function(m) {
                          msgs <<- c(msgs, conditionMessage(m))
                          invokeRestart("muffleMessage")
                        })
    expect_identical(grep("track failed", msgs, value = TRUE), character(0), info = g)
  }
})
