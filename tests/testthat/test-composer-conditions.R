# The composer's condition selection (R/app-conditions.R) on its own: hand-made track
# items (helper-composer-fakes.R), no genome. What the server sends to the track picker
# (update_rank_list()) and to the condition selects (updateSelectizeInput()) is recorded.

conditions_app <- function(ctx = fake_ctx()) function(input, output, session) {
  cs <- .pt_conditions_server(input, output, session, ctx)
}

# Runs `code` against the conditions server (it sees session, output, cs and `sent`);
# returns the rank-list updates (one list(id, keys, labels) each) and the select updates
# (one list(id, selected, choices) each).
drive_conditions <- function(code, ctx = fake_ctx()) {
  code <- rlang::enquo(code)
  sent <- new.env(); sent$ranks <- list(); sent$selects <- list()
  local_mocked_bindings(update_rank_list = function(css_id, text = NULL, labels = NULL, session = NULL)
    sent$ranks[[length(sent$ranks) + 1]] <- list(id = css_id, keys = names(labels), labels = labels))
  local_mocked_bindings(updateSelectizeInput = function(session, inputId, label = NULL, choices = NULL,
                                                        selected = NULL, ...)
    sent$selects[[length(sent$selects) + 1]] <- list(id = inputId, selected = selected, choices = choices))
  rlang::inject(shiny::testServer(conditions_app(ctx), !!code))
  sent
}

# The keys last sent to one rank list.
last_keys <- function(sent, id) {
  hits <- Filter(function(r) identical(r$id, id), sent$ranks)
  if (length(hits) == 0) return("<none sent>")
  k <- hits[[length(hits)]]$keys
  if (is.null(k)) character(0) else k
}

test_that("at startup, with no mode, both track lists are emptied", {
  s <- drive_conditions({
    session$flushReact()
    expect_identical(cs$figure_pool_key(), "not_ready")
  })
  expect_identical(last_keys(s, "available_tracks"), character(0))
  expect_identical(last_keys(s, "figure_tracks"), character(0))
  expect_length(s$ranks, 2)
})

test_that("Multi offers every track except comparisons, with Gene annotation in Figure tracks", {
  s <- drive_conditions({
    session$setInputs(view_mode = "compose")
    expect_identical(cs$figure_pool_key(), "compose")
  })
  expect_identical(last_keys(s, "available_tracks"),
                   c("bw_atac_A", "bw_atac_B", "genes", "hic_hic_A", "tpm_usage_Quant_A", "transcripts"))
  expect_identical(last_keys(s, "figure_tracks"), "genes")
  figure <- Filter(function(r) identical(r$id, "figure_tracks"), s$ranks)
  expect_identical(as.character(figure[[length(figure)]]$labels$genes), "<strong>Gene annotation</strong>")
})

test_that("Multi with a Multi start view starts from its tracks, Gene annotation kept in Available", {
  ctx <- fake_ctx(start_view = list(view_mode = "compose", tracks = c("bw_atac_A", "genes", "genes")))
  s <- drive_conditions(session$setInputs(view_mode = "compose"), ctx)
  expect_identical(last_keys(s, "figure_tracks"), c("bw_atac_A", "genes", "genes"))
  expect_identical(last_keys(s, "available_tracks"),
                   c("bw_atac_B", "genes", "hic_hic_A", "tpm_usage_Quant_A", "transcripts"))
})

test_that("a start view of another mode does not change Multi's start", {
  ctx <- fake_ctx(start_view = list(view_mode = "single", conditions = "A", tracks = "bw_atac_A"))
  s <- drive_conditions(session$setInputs(view_mode = "compose"), ctx)
  expect_identical(last_keys(s, "figure_tracks"), "genes")
})

test_that("One sample waits for its condition, then offers that sample's tracks", {
  s <- drive_conditions({
    session$setInputs(view_mode = "single")
    expect_identical(cs$figure_pool_key(), "not_ready")
    session$setInputs(sel_condition = "")
    expect_identical(cs$figure_pool_key(), "not_ready")
    session$setInputs(sel_condition = "A")
    expect_identical(cs$figure_pool_key(), "single::A")
  })
  expect_identical(last_keys(s, "available_tracks"),
                   c("bw_atac_A", "genes", "hic_hic_A", "tpm_usage_Quant_A", "transcripts"))
  expect_identical(last_keys(s, "figure_tracks"), "genes")
})

test_that("Two samples wait for both, then offer that pair's comparison", {
  s <- drive_conditions({
    # The page's boxes start blank (""), not NULL.
    session$setInputs(view_mode = "two", sel_reference = "B", sel_comparison = "")
    expect_identical(cs$figure_pool_key(), "not_ready")
    session$setInputs(sel_comparison = "A")
    expect_identical(cs$figure_pool_key(), "two::B|A")
  })
  expect_true("bw_diff_stacked_atac_B_vs_A" %in% last_keys(s, "available_tracks"))
  expect_false("bw_diff_stacked_atac_A_vs_B" %in% last_keys(s, "available_tracks"))
})

test_that("Two samples with the Comparison box never set builds the Reference's pool", {
  s <- drive_conditions({
    session$setInputs(view_mode = "two", sel_reference = "B")
    # current behaviour, possibly a bug: a NULL Comparison drops out of c(ref, cmp), so the
    # pool is built for the Reference alone. The page always sends "" for an empty box, so
    # this does not happen in the app; has_conditions_ready() stays FALSE meanwhile.
    expect_identical(cs$figure_pool_key(), "two::B")
    expect_false(cs$has_conditions_ready())
  })
  expect_identical(last_keys(s, "available_tracks"), c("bw_atac_B", "genes", "transcripts"))
})

test_that("a change that leaves the mode and conditions as they were does not rebuild the lists", {
  s <- drive_conditions({
    session$setInputs(view_mode = "single", sel_condition = "A")
    # Two-sample boxes change while in One sample: same pool, nothing sent.
    session$setInputs(sel_reference = "B")
    session$setInputs(sel_comparison = "A")
    session$setInputs(sel_condition = "B")
    expect_identical(cs$figure_pool_key(), "single::B")
  })
  # One rebuild for A (the first flush), one for B: two lists each.
  expect_identical(vapply(s$ranks, `[[`, "", "id"),
                   c("available_tracks", "figure_tracks", "available_tracks", "figure_tracks"))
  expect_true("bw_atac_A" %in% s$ranks[[1]]$keys)
  expect_true("bw_atac_B" %in% s$ranks[[3]]$keys)
})

test_that("Reference and Comparison: picking the other box's condition clears the other box", {
  s <- drive_conditions({
    session$flushReact()  # the observer skips its first run (ignoreInit), the blank start
    session$setInputs(sel_reference = "A")
    session$setInputs(sel_comparison = "B")
    expect_identical(c(cs$sel_prev$control, cs$sel_prev$treatment), c("A", "B"))
    # Reference just changed to Comparison's condition: Comparison is cleared.
    session$setInputs(sel_reference = "B")
    expect_identical(c(cs$sel_prev$control, cs$sel_prev$treatment), c("B", ""))
  })
  last <- tail(s$selects, 2)
  expect_identical(vapply(last, `[[`, "", "id"), c("sel_reference", "sel_comparison"))
  expect_identical(vapply(last, `[[`, "", "selected"), c("B", ""))
  # Both boxes always keep every condition.
  for (u in s$selects) expect_identical(u$choices, c("(select)" = "", A = "A", B = "B"))
})

test_that("Comparison changed to Reference's condition clears Reference", {
  s <- drive_conditions({
    session$flushReact()  # the observer skips its first run (ignoreInit), the blank start
    session$setInputs(sel_reference = "A")
    session$setInputs(sel_comparison = "A")
    expect_identical(c(cs$sel_prev$control, cs$sel_prev$treatment), c("", "A"))
  })
  expect_identical(vapply(tail(s$selects, 2), `[[`, "", "selected"), c("", "A"))
})

test_that("both boxes set to one condition at once keep Reference and clear Comparison", {
  s <- drive_conditions({
    session$flushReact()  # the observer skips its first run (ignoreInit), the blank start
    session$setInputs(sel_reference = "A", sel_comparison = "A")
    expect_identical(c(cs$sel_prev$control, cs$sel_prev$treatment), c("A", ""))
  })
})

test_that("distinct or blank picks are echoed back unchanged", {
  s <- drive_conditions({
    session$flushReact()  # the observer skips its first run (ignoreInit), the blank start
    session$setInputs(sel_reference = "", sel_comparison = "")
    expect_identical(c(cs$sel_prev$control, cs$sel_prev$treatment), c("", ""))
    session$setInputs(sel_reference = "B", sel_comparison = "A")
    expect_identical(c(cs$sel_prev$control, cs$sel_prev$treatment), c("B", "A"))
  })
  expect_identical(vapply(tail(s$selects, 2), `[[`, "", "selected"), c("B", "A"))
})

test_that("selected_conditions follows the mode", {
  drive_conditions({
    expect_null(cs$selected_conditions())
    session$setInputs(sel_condition = "A", sel_reference = "B", sel_comparison = "A")
    expect_identical(cs$selected_conditions(), "A")
    session$setInputs(view_mode = "two")
    expect_identical(cs$selected_conditions(), c("B", "A"))
    # Multi has no condition boxes; it still reads the One-sample box.
    session$setInputs(view_mode = "compose")
    expect_identical(cs$selected_conditions(), "A")
  })
})

test_that("has_conditions_ready: Multi always, One and Two once their boxes are filled", {
  drive_conditions({
    expect_false(cs$has_conditions_ready())
    session$setInputs(view_mode = "compose")
    expect_true(cs$has_conditions_ready())
    session$setInputs(view_mode = "single")
    expect_false(cs$has_conditions_ready())
    session$setInputs(sel_condition = "A")
    expect_true(cs$has_conditions_ready())
    session$setInputs(view_mode = "two", sel_reference = "A")
    expect_false(cs$has_conditions_ready())
    session$setInputs(sel_comparison = "B")
    expect_true(cs$has_conditions_ready())
    session$setInputs(view_mode = "something else")
    expect_false(cs$has_conditions_ready())
  })
})

test_that("the tracks hint says what to pick first, and goes once the pick is made", {
  drive_conditions({
    hint <- function() ui_html(output$tracks_hint)
    expect_identical(hint(), '<style>#tracks-block { --pt-tracks-hint: "Pick the number of samples above to start."; }</style>')
    session$setInputs(view_mode = "single")
    expect_match(hint(), '"Pick a condition above to see its tracks."', fixed = TRUE)
    session$setInputs(view_mode = "two")
    expect_match(hint(), '"Pick a reference and a comparison above to see their tracks."', fixed = TRUE)
    session$setInputs(view_mode = "compose")
    expect_null(output$tracks_hint)
  })
})
