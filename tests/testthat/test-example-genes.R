# The bundled mouse example (pt_example(), inst/extdata/example/; also the paper's): five
# genes, every data type. Uses only the installed example data, so these run on
# R-universe too. Needs org.Mm.eg.db (Suggests).

EXAMPLE_GENES <- c("Eya2", "Lrrfip1", "Acot7", "Dzip1", "Cxcl14")
EXAMPLE_DIR <- file.path(tempdir(), "pt-test-example-genes")

example_project <- function() {
  # txdbmaker's two notices for any GENCODE GTF are expected; anything else still warns.
  withCallingHandlers(
    pt_project(EXAMPLE_DIR, pt_example()),
    warning = function(w) {
      if (grepl("\"phase\" metadata column|genome version information", conditionMessage(w)))
        invokeRestart("muffleWarning")
    })
}

test_that("pt_example() gives the installed five-gene example, small enough to ship", {
  txt <- pt_example()
  expect_same_path(dirname(txt), system.file("extdata", "example", package = "PlotTrellis"))
  files <- list.files(dirname(txt), recursive = TRUE, full.names = TRUE)
  # Under 5 MB in all and 1 MB per file: R-universe builds the package with it.
  expect_lt(sum(file.size(files)), 5e6)
  expect_true(all(file.size(files) < 1e6), info = basename(files[file.size(files) >= 1e6]))
  res <- import_samples_txt(txt, defaults = .pt_default_settings())
  expect_length(res$errors, 0)
  expect_length(res$warnings, 0)
})

test_that("the example opens offline with every data type", {
  skip_if_not_installed("org.Mm.eg.db")
  # The shipped chromosome sizes must be used, not fetched from UCSC.
  local_mocked_bindings(build_chrom_sizes = function(...) stop("tried to download chromosome sizes"))
  proj <- example_project()
  expect_equal(pt_samples(proj)$name, c("Basal", "LP"))
  expect_setequal(unique(pt_tracks(proj)$track),
                  c("genes", "protein_domains", "transcripts_all", "transcripts", "transcripts_annotated",
                    "RNA", "junctions", "transcript_expression", "ATAC", "ATAC_peaks", "H3K36me3",
                    "H3K36me3_peaks", "IgG", "HiC", "RNA_diff", "ATAC_diff", "H3K36me3_diff", "IgG_diff",
                    "junctions_diff", "junctions_diff_tpm"))
  # Each gene opens by name, on its window's chromosome.
  chroms <- c(Eya2 = "chr2", Lrrfip1 = "chr1", Acot7 = "chr4", Dzip1 = "chr14", Cxcl14 = "chr13")
  for (g in EXAMPLE_GENES)
    expect_equal(pt_figure(proj, g, "genes")$region_params$chrom, chroms[[g]], info = g)
})

# Every track type at each gene, both samples and the comparison: no track may fail (a
# failed track is drawn as a message in the figure and said on the console). Protein
# domains come from UniProt over the internet, so they are tried below only when online.
# About 2 minutes, so not on CRAN-mode builds (R-universe); bigWig can't be read on Windows.
test_that("each of the example's genes draws every track type with no error", {
  skip_if_not_installed("org.Mm.eg.db")
  skip_on_cran()
  skip_on_os("windows")
  proj <- example_project()
  tracks <- setdiff(unique(pt_tracks(proj)$track), "protein_domains")
  for (g in EXAMPLE_GENES) {
    fig <- pt_figure(proj, g, tracks, comparison = c("LP", "Basal"))
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

test_that("the example's genes draw their protein domains (online)", {
  skip_if_not_installed("org.Mm.eg.db")
  skip_on_cran()
  skip_if_offline()
  proj <- example_project()
  for (g in EXAMPLE_GENES) {
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

# The bundled Hi-C has the .hic's bin sizes from 5 kb up (1 kb left out, owner 2026-11-19),
# each with the .hic's normalisations, so Auto (500 min bins) draws every gene at 5 kb, and
# Fixed 1 kb says which sizes there are. Reads the .mcool only (no bigWig): runs everywhere.
test_that("the example Hi-C has 5 kb and coarser bins; Auto draws the genes at 5 kb; Fixed 1 kb says the closest", {
  skip_if_not_installed("org.Mm.eg.db")
  proj <- example_project()
  hic <- file.path(dirname(pt_example()), "hic", c("Basal_HiC.mcool", "LP_HiC.mcool"))
  sizes <- c(5000, 10000, 25000, 50000, 1e5, 2.5e5, 5e5, 1e6, 2.5e6)
  for (f in hic) {
    expect_equal(as.numeric(plotgardener::readCoolBpResolutions(f)), sizes, info = basename(f))
    for (r in sizes)
      expect_setequal(plotgardener::readCoolNorms(f, resolution = r), c("KR", "SCALE", "VC", "VC_SQRT", "NONE"))
  }
  for (g in EXAMPLE_GENES) {
    fig <- pt_figure(proj, g, "genes")
    expect_identical(fig$settings$hic_res_mode, "auto")
    expect_equal(fig$settings$hic_max_bins_per_view, 500)
    for (f in hic)
      expect_equal(load_hic(list(path = f), fig$region_params, fig$settings)$data$resolution, 5000,
                   info = paste(g, basename(f)))
  }
  fixed <- modifyList(fig$settings, list(hic_res_mode = "fixed", hic_resolution = 1000))
  expect_error(load_hic(list(path = hic[1]), fig$region_params, fixed),
               "this file has no 1,000 bp bins; closest: 5,000, 10,000, 25,000 bp.", fixed = TRUE)
})
