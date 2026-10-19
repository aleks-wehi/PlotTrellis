# Decompresses golden-image snapshots (tests/testthat/_snaps_gz/**/*.svg.gz)
# into the plain .svg vdiffr::expect_doppelganger() actually reads, under
# tests/testthat/_snaps/. The largest Hi-C fixtures stay gzipped in git
# (~12-13x smaller) to keep blobs under GitHub's push size limits; this
# restores the plain .svg on disk before any test runs so vdiffr never has to
# know about the .gz extension.
#
# The .gz files live in a sibling _snaps_gz/ directory, not inside _snaps/
# itself: testthat::test_dir() prunes any file under _snaps/<test-file>/ that
# a test run didn't touch as an "unused snapshot", and it doesn't recognise
# .svg.gz as a snapshot format, so it silently deleted these 6 files after
# every run when they lived alongside the plain .svg (BUGS.md #45).
#
# To update one of these snapshots: run the normal vdiffr snapshot-update
# workflow (produces a plain .svg under _snaps/), then replace the
# corresponding _snaps_gz/ file with a gzip'd copy and delete the plain .svg,
# e.g.:
#   gzip -9 -c tests/testthat/_snaps/foo/bar.svg > tests/testthat/_snaps_gz/foo/bar.svg.gz
#   rm tests/testthat/_snaps/foo/bar.svg

# Gated on PT_TEST_DATA_AVAILABLE (helper-fixtures.R, sourced first -- "f" <
# "g") -- every expect_doppelganger() call is behind
# skip_if_not(PT_TEST_DATA_AVAILABLE) too, so without fixture data none of
# these references would ever be touched by a test in the run. Decompressing
# them anyway left all 45 sitting in _snaps/ untouched, and testthat's own
# unused-snapshot pruning (see above) turns that into an R CMD check FAILURE
# on any machine without /vast data -- R-universe's build farm, CRAN, any
# fresh clone -- every single time.
if (PT_TEST_DATA_AVAILABLE) {
  snaps_gz_dir <- testthat::test_path("_snaps_gz")
  snap_gz_files <- list.files(snaps_gz_dir, pattern = "\\.svg\\.gz$", recursive = TRUE, full.names = TRUE)
  for (gz_path in snap_gz_files) {
    rel_path <- sub("\\.gz$", "", substring(gz_path, nchar(snaps_gz_dir) + 2))
    svg_path <- testthat::test_path("_snaps", rel_path)
    dir.create(dirname(svg_path), recursive = TRUE, showWarnings = FALSE)
    con <- gzfile(gz_path, "rt")
    writeLines(readLines(con, warn = FALSE), svg_path)
    close(con)
  }
}
