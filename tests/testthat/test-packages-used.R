# R CMD check's own dependency scan of the code: every package reached with `::`,
# requireNamespace() or loadNamespace() is declared in DESCRIPTION, and the code never calls its
# own namespace with `:::`. Run here so a new dependency cannot slip in unnoticed before a check.

test_that("R/ code uses only declared packages and no own-namespace :::", {
  pkg_dir <- test_path("..", "..")
  skip_if_not(file.exists(file.path(pkg_dir, "DESCRIPTION")) && dir.exists(file.path(pkg_dir, "R")),
              "source tree not available (built package)")
  check <- get0(".check_packages_used", envir = asNamespace("tools"))
  skip_if(is.null(check), "tools:::.check_packages_used() not in this R")
  findings <- capture.output(print(check(dir = pkg_dir)))
  expect_true(length(findings) == 0, info = paste(findings, collapse = "\n"))
})
