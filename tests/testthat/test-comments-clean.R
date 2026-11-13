# Comments say what the code does now, not where it came from: no bug, feature or work-item
# numbers, project phases, dates, or names of files outside the package. Checks the R/ source
# (when the source tree is there) and the source text pasted into downloaded scripts
# (PT_SOURCE_TEXT, with each function's "comment" attribute), so a reader of either never meets
# a reference they can't follow.

COMMENT_HISTORY_PATTERNS <- c(
  "\\b[Bb]ug [0-9]",
  "\\b[Ff]eature [0-9]",
  "\\b[Pp]hase[- ][A-Z]\\b",                # "Phase C.3", "phase-I.3"
  "\\b[A-Z]\\.[0-9]+[a-z]? [Ss]tep [0-9]+", # "I.3 step 28"; a plain "step 2" is fine
  "\\bW[0-9]{2}\\b",
  "BUGS\\.md",
  "FEATURE_IDEAS",
  "solutions/",
  "utils/",
  "\\bmvp[0-9]",
  "20[0-9]{2}-[0-9]{2}-[0-9]{2}",
  "\u2014"                              # em dash
)

# "<where>:<line>: <text>" for every line of `lines` matching a pattern.
comment_history_hits <- function(lines, where) {
  hit <- Reduce(`|`, lapply(COMMENT_HISTORY_PATTERNS, grepl, x = lines, perl = TRUE))
  sprintf("%s:%d: %s", where, which(hit), trimws(lines[hit]))
}

test_that("R/ files carry no history references", {
  r_dir <- test_path("..", "..", "R")
  skip_if_not(dir.exists(r_dir), "R/ not available (built package)")
  files <- list.files(r_dir, pattern = "\\.R$", full.names = TRUE)
  hits <- unlist(lapply(files, function(f)
    comment_history_hits(readLines(f, warn = FALSE, encoding = "UTF-8"), file.path("R", basename(f)))))
  expect_true(length(hits) == 0,
              info = paste0(length(hits), " line(s):\n", paste(hits, collapse = "\n")))
})

test_that("the source text pasted into downloaded scripts carries no history references", {
  hits <- unlist(lapply(names(PT_SOURCE_TEXT), function(nm) {
    src <- PT_SOURCE_TEXT[[nm]]
    c(comment_history_hits(strsplit(src, "\n", fixed = TRUE)[[1]], sprintf("PT_SOURCE_TEXT[[\"%s\"]]", nm)),
      if (!is.null(attr(src, "comment")))
        comment_history_hits(strsplit(attr(src, "comment"), "\n", fixed = TRUE)[[1]],
                             sprintf("PT_SOURCE_TEXT[[\"%s\"]] comment", nm)))
  }))
  expect_true(length(hits) == 0,
              info = paste0(length(hits), " line(s):\n", paste(hits, collapse = "\n")))
})
