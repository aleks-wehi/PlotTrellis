# Two quick "Remove" clicks in the wizard's sample list: the second click comes from the table as
# it was before the first removal. Each must remove the sample on the row clicked. Drives the real
# wizard_server() headlessly, no fixture data needed.

test_that("two quick Removes on different rows remove those two samples (no fixture data)", {
  json <- file.path(tempfile("pt-proj-"), "project_settings.json")
  defaults <- .pt_default_settings()
  shiny::testServer(wizard_server(json, defaults, defaults), {
    wiz_samples(lapply(c("A's \"x\"", "B", "C"), function(cond)
      list(condition = cond, colour = "#1B9E77", data = list(list(type = "bw", name = "ATAC", path = "a.bw")))))
    # The value each row's Remove button sends, read off the rendered list (one render, as the
    # browser has it before either click is handled).
    html <- as.character(output$wiz_samples_ui$html)
    for (e in list(c("&#39;", "'"), c("&quot;", "\""), c("&lt;", "<"), c("&gt;", ">"), c("&amp;", "&")))
      html <- gsub(e[1], e[2], html, fixed = TRUE)
    sends <- regmatches(html, gregexpr("setInputValue\\('wiz_remove_sample', .*?, \\{priority", html, perl = TRUE))[[1]]
    sends <- lapply(sub("^setInputValue\\('wiz_remove_sample', (.*), \\{priority$", "\\1", sends), jsonlite::fromJSON)
    expect_length(sends, 3)
    session$setInputs(wiz_remove_sample = sends[[1]])  # row A (a name that needs escaping)
    session$setInputs(wiz_remove_sample = sends[[2]])  # row B, from the same (now stale) table
    expect_identical(vapply(wiz_samples(), `[[`, "", "condition"), "C")
    # A repeat of an already-handled Remove changes nothing.
    session$setInputs(wiz_remove_sample = sends[[2]])
    expect_identical(vapply(wiz_samples(), `[[`, "", "condition"), "C")
  })
})

test_that("Edit after a quick Remove opens the sample on the row clicked (no fixture data)", {
  json <- file.path(tempfile("pt-proj-"), "project_settings.json")
  defaults <- .pt_default_settings()
  shiny::testServer(wizard_server(json, defaults, defaults), {
    wiz_samples(lapply(c("A", "B", "C"), function(cond)
      list(condition = cond, colour = "#1B9E77", data = list(list(type = "bw", name = cond, path = "a.bw")))))
    html <- gsub("&quot;", "\"", gsub("&#39;", "'", as.character(output$wiz_samples_ui$html), fixed = TRUE), fixed = TRUE)
    edits <- regmatches(html, gregexpr("setInputValue\\('wiz_edit_sample', .*?, \\{priority", html, perl = TRUE))[[1]]
    edits <- lapply(sub("^setInputValue\\('wiz_edit_sample', (.*), \\{priority$", "\\1", edits), jsonlite::fromJSON)
    expect_length(edits, 3)
    session$setInputs(wiz_remove_sample = "A")
    session$setInputs(wiz_edit_sample = edits[[3]])  # row C, from the table before A was removed
    expect_identical(modal_edit_index(), 2L)
    expect_identical(modal_instances()[[1]]$name, "C")
  })
})

test_that("a Remove handled while an Edit is open doesn't make Save overwrite another sample (no fixture data)", {
  json <- file.path(tempfile("pt-proj-"), "project_settings.json")
  defaults <- .pt_default_settings()
  conditions <- function(samples) vapply(samples, `[[`, "", "condition")
  shiny::testServer(wizard_server(json, defaults, defaults), {
    invisible(output$wiz_samples_ui)  # a first flush, so the events below aren't skipped as initial values
    wiz_samples(lapply(c("A", "B", "C"), function(cond) list(condition = cond, colour = "#1B9E77", data = list())))
    session$setInputs(wiz_edit_sample = "B")
    session$setInputs(wiz_remove_sample = "A")  # clicked before the Edit modal appeared
    session$setInputs(wiz_sample_name = "B2", wiz_confirm_sample = 1)
    expect_identical(conditions(wiz_samples()), c("B2", "C"))
    # Saving under the unchanged name isn't a clash with itself.
    session$setInputs(wiz_edit_sample = "C")
    session$setInputs(wiz_remove_sample = "B2")
    session$setInputs(wiz_sample_name = "C", wiz_confirm_sample = 2)
    expect_identical(conditions(wiz_samples()), "C")
    expect_null(modal_edit_name())
    # The sample being edited is removed meanwhile: Save closes without adding it back.
    session$setInputs(wiz_edit_sample = "C")
    session$setInputs(wiz_remove_sample = "C")
    session$setInputs(wiz_sample_name = "C", wiz_confirm_sample = 3)
    expect_length(wiz_samples(), 0)
  })
})
