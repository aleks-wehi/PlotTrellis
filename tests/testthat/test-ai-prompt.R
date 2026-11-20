# The AI prompt for writing a samples.txt (inst/samples_txt_ai_prompt.md, shown in the
# vignettes) must keep up with what the importer reads, so it cannot drift.

prompt_text <- function() {
  f <- system.file("samples_txt_ai_prompt.md", package = "PlotTrellis")
  expect_true(nzchar(f))
  paste(readLines(f, warn = FALSE), collapse = "\n")
}

test_that("the AI prompt names every data type the importer accepts, as a type:name column", {
  txt <- prompt_text()
  for (type in DATA_TYPE_VOCAB) {
    expect_true(grepl(paste0("\n  ", type, " "), txt, fixed = TRUE), info = paste("type listed:", type))
  }
  # The fixed columns, spelled as the importer spells them.
  expect_match(txt, "condition", fixed = TRUE)
  expect_match(txt, "colour", fixed = TRUE)
})

test_that("the AI prompt names every file extension the app looks for", {
  txt <- prompt_text()
  # Every extension in the wizard's file patterns, e.g. "\\.(bw|bigWig)$" gives bw and bigWig.
  exts <- unique(unlist(lapply(WIZARD_DATA_TYPES, function(d)
    strsplit(sub("^\\\\\\.\\(([^)]*)\\).*$", "\\1", d$ext), "|", fixed = TRUE)[[1]])))
  expect_gt(length(exts), 8)
  # as ".bw" or as a word ("BED12", "narrowPeak"), in any case
  for (e in exts) expect_true(grepl(paste0("\\b", e, "\\b"), txt, ignore.case = TRUE), info = paste("extension:", e))
  expect_match(txt, ".gz", fixed = TRUE)          # bed and paired-end may be gzipped
  expect_match(txt, "quant.sf", fixed = TRUE)     # quantdir is the folder holding it
})

test_that("the AI prompt asks about R and the data first, and ends by opening the project", {
  txt <- prompt_text()
  expect_match(txt, "How comfortable are you with R", fixed = TRUE)
  expect_match(txt, "Do you know where your data files are", fixed = TRUE)
  expect_match(txt, "ask the two questions.", fixed = TRUE)   # closing line
  for (h in c("\n6. Show the whole table", "\n7. Check before writing", "\n8. Check their setup",
              "\n9. Write the file and open it"))
    expect_true(grepl(h, txt, fixed = TRUE), info = h)
  expect_match(txt, 'run_plottrellis(project_dir = "my-project",', fixed = TRUE)
  expect_false(grepl("[^\\x01-\\x7F]", txt, perl = TRUE))  # project rule: ASCII only, no em dash
})

test_that("the AI assistant vignette shows the prompt from the shipped file", {
  pkg <- testthat::test_path("..", "..")
  skip_if_not(dir.exists(file.path(pkg, "vignettes")), "vignette sources not available (installed package)")
  for (v in "samples-txt-ai-prompt.Rmd") {
    src <- readLines(file.path(pkg, "vignettes", v), warn = FALSE)
    expect_true(any(grepl("system.file(\"samples_txt_ai_prompt.md\"", src, fixed = TRUE)), info = v)
  }
})
