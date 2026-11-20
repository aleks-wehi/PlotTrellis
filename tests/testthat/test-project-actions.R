# The app's project actions (R/app-project-actions.R) on their own: no genome, no
# fixture data, a temp project folder. The delayed stopApp() is caught, not run.

# A server made of nothing but the project actions, for one project folder.
project_actions_app <- function(proj, dirty = FALSE) {
  json <- file.path(proj, "project_settings.json")
  function(input, output, session) {
    .pt_project_actions_server(input, session, function() .pt_unsaved_settings_note(dirty),
                               project_dir = proj, settings_json = json, cache_dir = file.path(proj, "cache"))
  }
}

new_project <- function() {
  proj <- tempfile("pt-actions-"); dir.create(file.path(proj, "cache"), recursive = TRUE)
  writeLines("{}", file.path(proj, "project_settings.json"))
  writeLines("cached", file.path(proj, "cache", "x.rds"))
  writeLines("keep", file.path(proj, "keep.txt"))
  proj
}

# Runs `steps` against the project-actions server, returning the dialogs shown and the
# stopApp() calls scheduled (as deparsed text).
drive <- function(proj, dirty = FALSE, steps) {
  shown <- list(); scheduled <- list()
  local_mocked_bindings(showModal = function(ui, ...) shown[[length(shown) + 1]] <<- as.character(ui))
  local_mocked_bindings(later = function(func, ...) { scheduled[[length(scheduled) + 1]] <<- func; invisible(NULL) },
                        .package = "later")
  shiny::testServer(project_actions_app(proj, dirty), { session$flushReact(); steps(session) })
  list(shown = shown,
       stops = grep("stopApp", unlist(lapply(scheduled, function(f) paste(deparse(body(f)), collapse = " "))), value = TRUE))
}

test_that("the app's server wires the project actions to its own paths, in order", {
  # Guards the positional call in app-screen-server.R without fixture data: the call's arguments, by name.
  app_r <- testthat::test_path("..", "..", "R", "app-screen-server.R")
  skip_if_not(file.exists(app_r), "source tree not available (installed package)")
  src <- paste(trimws(readLines(app_r, warn = FALSE)), collapse = " ")
  call <- regmatches(src, gregexpr("\\.pt_project_actions_server\\([^)]*\\)", src))[[1]]
  expect_length(call, 1)
  expect_identical(call, paste(".pt_project_actions_server(input, session, unsaved_settings_note,",
                               "ctx$PROJECT_DIR, SETTINGS_JSON, ctx$CACHE_DIR)"))
})

test_that("Edit project settings confirms, names the project folder, and stops with 'edit_settings'", {
  proj <- new_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  r <- drive(proj, dirty = TRUE, function(session) {
    session$setInputs(open_settings = 1)
    session$setInputs(confirm_edit_settings = 1)
  })
  expect_length(r$shown, 2)
  expect_match(r$shown[[1]], proj, fixed = TRUE)
  expect_match(r$shown[[1]], "Your Figure settings changes are not saved", fixed = TRUE)
  expect_match(r$shown[[1]], '<details id="settings_advanced">', fixed = TRUE)
  expect_identical(r$stops, 'stopApp(returnValue = "edit_settings")')
  expect_true(file.exists(file.path(proj, "project_settings.json")))
})

test_that("Reset deletes project_settings.json and nothing else, then restarts", {
  proj <- new_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  r <- drive(proj, steps = function(session) {
    session$setInputs(reset_project_btn = 1)
    session$setInputs(confirm_reset_project = 1)
  })
  expect_false(file.exists(file.path(proj, "project_settings.json")))
  expect_true(file.exists(file.path(proj, "cache", "x.rds")))
  expect_true(file.exists(file.path(proj, "keep.txt")))
  expect_identical(r$stops, 'stopApp(returnValue = "restart")')
})

test_that("Clear all caches removes the cache folder only, then restarts", {
  proj <- new_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  r <- drive(proj, steps = function(session) {
    session$setInputs(clear_caches_btn = 1)
    session$setInputs(confirm_clear_caches = 1)
  })
  expect_match(r$shown[[1]], file.path(proj, "cache"), fixed = TRUE)
  expect_false(dir.exists(file.path(proj, "cache")))
  expect_true(file.exists(file.path(proj, "project_settings.json")))
  expect_true(file.exists(file.path(proj, "keep.txt")))
  expect_identical(r$stops, 'stopApp(returnValue = "restart")')
})

test_that("Edit project settings and Clear all caches show no unsaved note when the Figure settings are saved", {
  proj <- new_project(); on.exit(unlink(proj, recursive = TRUE), add = TRUE)
  r <- drive(proj, dirty = FALSE, function(session) {
    session$setInputs(open_settings = 1)
    session$setInputs(clear_caches_btn = 1)
  })
  expect_length(r$shown, 2)
  expect_false(grepl("not saved", r$shown[[1]], fixed = TRUE))
  expect_false(grepl("not saved", r$shown[[2]], fixed = TRUE))
  expect_length(r$stops, 0)
})
