# Paths compared as files, not strings. On Windows the same file can be
# spelled several ways (tempdir()'s backslashes and RUNNER~1 short names vs
# normalizePath()'s long names, dirname() turning "\" into "/"), so the
# R-universe Windows checks failed on paths that pointed at the same place.
norm_path <- function(p) normalizePath(p, winslash = "/", mustWork = FALSE)

expect_same_path <- function(object, expected) {
  expect_identical(norm_path(object), norm_path(expected))
}
