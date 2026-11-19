# import_samples_txt() warns when a samples.txt has one column per sample instead of one per
# data type (a real user's layout, W75). No fixture data needed.

write_samples <- function(lines) {
  path <- tempfile(fileext = ".txt")
  writeLines(lines, path)
  path
}

# A real user's layout, cut down: each sample's data in columns of its own, filled on its row only.
one_column_per_sample <- c(
  "condition\tbw:Grans ATAC\tbw:Grans CUT&RUN\tbw:pDC ATAC\thic:Grans HiC\thic:pDC HiC",
  "Grans\tatac/Grans.bw\tcutrun/Grans.bw\t\thic/Grans.hic\t",
  "pDC\t\t\tatac/pDC.bw\t\thic/pDC.hic"
)

test_that("one column per sample warns once, naming the columns and the merged layout (no fixture data)", {
  res <- import_samples_txt(write_samples(one_column_per_sample))
  expect_length(res$errors, 0)
  expect_length(res$layout_warning, 1)
  msg <- res$layout_warning
  for (col in c("bw:Grans ATAC", "bw:pDC ATAC", "hic:Grans HiC", "hic:pDC HiC"))
    expect_match(msg, col, fixed = TRUE)
  # The suggested layout: the columns' shared name without the sample's name, once per data type.
  expect_match(msg, "a row per sample: bw:ATAC, hic:HiC.", fixed = TRUE)
  expect_match(msg, 'See vignette("PlotTrellis"), section "Your own data in code (writing a samples.txt)".', fixed = TRUE)
  # Grans' CUT&RUN has no pDC counterpart: an extra data type, not a per-sample column.
  expect_false(grepl("CUT&RUN", msg, fixed = TRUE))
  # A warning, not an error: the samples still import.
  expect_length(res$samples, 2)
})

test_that("a normal file, a one-sample file and the bundled examples don't warn (no fixture data)", {
  normal <- c("condition\tbw:ATAC\thic:HiC\tpaired-end:Loops",
              "Grans\tatac/Grans.bw\thic/Grans.hic\tloops/Grans.bedpe",
              "pDC\tatac/pDC.bw\thic/pDC.hic\t")
  expect_length(import_samples_txt(write_samples(normal))$layout_warning, 0)
  one_sample <- c("condition\tbw:Grans ATAC\thic:Grans HiC", "Grans\tatac/Grans.bw\thic/Grans.hic")
  expect_length(import_samples_txt(write_samples(one_sample))$layout_warning, 0)
  for (ex in c("example", "example-human")) {
    txt <- system.file("extdata", ex, "samples.txt", package = "PlotTrellis")
    expect_length(import_samples_txt(txt)$layout_warning, 0)
  }
})

test_that("different data on different samples is a normal file and doesn't warn (no fixture data)", {
  extras <- c("condition\tbw:ATAC\tbw:CUT&RUN\tbw:H3K27ac\tpaired-end:Loops\tpaired-end:Junctions",
              "A\ta.bw\tcutrun.bw\t\tloops.bedpe\t",
              "B\tb.bw\t\th3k27ac.bw\t\tjunctions.bed")
  expect_length(import_samples_txt(write_samples(extras))$layout_warning, 0)
  # Names that only differ by the sample's name but sit in one row each still warn, whatever the case.
  mixed <- c("condition\tbw:grans_atac\tbw:PDC-ATAC", "Grans\tg.bw\t", "pDC\t\tp.bw")
  expect_match(import_samples_txt(write_samples(mixed))$layout_warning, "bw:atac", fixed = TRUE)
  colons <- c("condition\tbw:Grans:ATAC\tbw:pDC:ATAC", "Grans\tg.bw\t", "pDC\t\tp.bw")
  expect_match(import_samples_txt(write_samples(colons))$layout_warning, "bw:ATAC.", fixed = TRUE)
})

test_that("parse_samples_txt() (and so pt_project()) gives the layout warning once, on a successful import", {
  skip_if_not_installed("org.Mm.eg.db")
  # Written to a temporary folder, reading the installed example's data by absolute path: the
  # installed example must not change.
  dir <- dirname(pt_example("mammary"))
  shipped <- list.files(dir, recursive = TRUE, all.files = TRUE)
  per_sample <- c(
    paste("# data_root:", dir), paste("# gtf_path:", file.path(dir, "annotation.gtf.gz")),
    "# orgdb: org.Mm.eg.db", "# genome_label: mm39",
    "condition\tbw:Basal atac\tbw:LP atac",
    "Basal\tatac/Basal_ATAC.bw\t",
    "LP\t\tatac/LP_ATAC.bw")
  path <- tempfile("per-sample-", fileext = ".txt")
  writeLines(per_sample, path)
  warns <- character(0)
  settings <- withCallingHandlers(
    parse_samples_txt(path, project_dir = tempfile("pt-layout-")),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_length(settings$samples, 2)
  expect_identical(sum(grepl("are each used by one sample only", warns, fixed = TRUE)), 1L)
  expect_true(any(grepl("a row per sample: bw:atac.", warns, fixed = TRUE)))
  expect_identical(list.files(dir, recursive = TRUE, all.files = TRUE), shipped)
})
