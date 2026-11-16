# R CMD check warns when a code file holds a non-ASCII character outside comments (R-universe
# builds on Windows and Linux with different native encodings). Non-ASCII text in code is
# written escaped, "\u00B7" for a middle dot, never the character itself. This runs the same check R CMD check
# does, so a new one cannot slip in unnoticed.

test_that("R/ code has no non-ASCII characters outside comments", {
  pkg_dir <- test_path("..", "..")
  skip_if_not(dir.exists(file.path(pkg_dir, "R")), "R/ not available (built package)")
  check <- get0(".check_package_ASCII_code", envir = asNamespace("tools"))
  skip_if(is.null(check), "tools:::.check_package_ASCII_code() not in this R")
  bad <- check(pkg_dir)
  expect_true(length(bad) == 0, info = paste("non-ASCII code in:", paste(bad, collapse = ", ")))
})
