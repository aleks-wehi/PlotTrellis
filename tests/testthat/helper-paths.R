# Paths compared as files, not strings. On Windows the same file can be
# spelled several ways (tempdir()'s backslashes and RUNNER~1 short names vs
# normalizePath()'s long names, dirname() turning "\" into "/", a drive letter
# added to "/x"), and on macOS tempdir() is under /var, a symlink to
# /private/var. normalizePath() resolves only paths that exist, so the deepest
# existing folder is resolved and the rest appended.
norm_path <- function(p) {
  vapply(p, function(x) {
    rest <- character()
    while (!file.exists(x) && dirname(x) != x) {
      rest <- c(basename(x), rest)
      x <- dirname(x)
    }
    paste(c(sub("/$", "", normalizePath(x, winslash = "/", mustWork = FALSE)), rest), collapse = "/")
  }, "", USE.NAMES = FALSE)
}

expect_same_path <- function(object, expected) {
  expect_identical(norm_path(object), norm_path(expected))
}
