# The render-time check (helper-fixtures.R) must not fail on one slow draw of a
# busy node, but must still fail on a render that is slow every time.

test_that("a time under the ceiling is kept and never re-timed", {
  called <- FALSE
  expect_equal(perf_seconds_to_judge(0.15, 0.10, function() { called <<- TRUE; 9 }), 0.15)
  expect_false(called)
  expect_equal(perf_seconds_to_judge(0.15, 0.10), 0.15)   # no retime given
})

test_that("a slow first draw is judged by the faster second one", {
  expect_equal(perf_seconds_to_judge(4.43, 0.10, function() 0.11), 0.11)
})

test_that("a render slow both times still counts as slow", {
  expect_equal(perf_seconds_to_judge(4.43, 0.10, function() 5), 4.43)
})

test_that("a failing re-draw leaves the first time standing", {
  expect_equal(perf_seconds_to_judge(4.43, 0.10, function() stop("boom")), 4.43)
})
