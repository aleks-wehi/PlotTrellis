# The wizard's Close without saving and its dirty guard (.pt_wizard_dirty_server(),
# R/wizard-dirty.R) on their own: plain reactives for the fields, no genome. The dialogs are
# caught, not shown, and the delayed stopApp() is caught, not run.

edit_cfg <- function() {
  cfg <- utils::modifyList(.pt_default_settings(), list(
    data_dir = "/data/proj", gtf_path = "/x/genes.gtf.gz", orgdb = "org.Mm.eg.db", genome_label = "mm10"))
  cfg$samples <- list(list(condition = "Ba", colour = "#1B9E77", data = list()))  # not merged by modifyList()
  cfg
}

# The settings wiz_collect_flat_settings() gives for `cfg` unchanged.
flat_of <- function(cfg) c(cfg[c("gtf_path", "orgdb", "genome_label")],
                           .pt_collect_display_settings(list(), "", cfg))

# A server made of nothing but the dirty block, over fields that start as `cfg` has them.
dirty_app <- function(cfg, is_edit = TRUE) {
  function(input, output, session) {
    data_root   <- reactiveVal(cfg$data_dir)
    wiz_samples <- reactiveVal(cfg$samples)
    flat        <- reactiveVal(flat_of(cfg))
    res <- .pt_wizard_dirty_server(input, cfg, is_edit, data_root, wiz_samples, function() flat())
  }
}

# Catches the dialogs shown and the stopApp() calls scheduled (as deparsed text), for the rest
# of the calling test.
catch_close <- function(env = parent.frame()) {
  got <- new.env(); got$shown <- character(0); got$stops <- character(0); got$removed <- 0L
  local_mocked_bindings(showModal = function(ui, ...) got$shown <- c(got$shown, as.character(ui)),
                        removeModal = function(...) got$removed <- got$removed + 1L, .env = env)
  local_mocked_bindings(later = function(func, ...) {
    got$stops <- c(got$stops, grep("stopApp", deparse(body(func)), value = TRUE))
    invisible(NULL)
  }, .package = "later", .env = env)
  got
}

test_that("a project opened unchanged is not dirty; each kind of change makes it dirty", {
  cfg <- edit_cfg()
  shiny::testServer(dirty_app(cfg), {
    expect_false(res$is_dirty())
    data_root("/data/other"); expect_true(res$is_dirty())
    data_root(cfg$data_dir);  expect_false(res$is_dirty())
    wiz_samples(list()); expect_true(res$is_dirty())
    wiz_samples(cfg$samples); expect_false(res$is_dirty())
    flat(utils::modifyList(flat_of(cfg), list(genome_label = "mm39"))); expect_true(res$is_dirty())
    # A number equal in value is no change.
    flat(utils::modifyList(flat_of(cfg), list(fontsize = as.integer(cfg$fontsize)))); expect_false(res$is_dirty())
    flat(utils::modifyList(flat_of(cfg), list(fontsize = cfg$fontsize + 1))); expect_true(res$is_dirty())
  })
})

test_that("a first-time setup is never dirty", {
  cfg <- .pt_default_settings()
  shiny::testServer(dirty_app(cfg, is_edit = FALSE), {
    data_root("/data/anything")
    expect_false(res$is_dirty())
  })
})

test_that("Close with nothing changed closes at once and restarts into the app", {
  got <- catch_close()
  expect_message(shiny::testServer(dirty_app(edit_cfg()), session$setInputs(wiz_close = 1)),
                 "settings unchanged", fixed = TRUE)
  expect_length(got$shown, 1)
  expect_match(got$shown[[1]], "Returning to the app without saving.", fixed = TRUE)
  expect_identical(trimws(got$stops), 'stopApp(returnValue = "restart")')
})

test_that("Close with changes asks first; Discard and close then closes", {
  got <- catch_close()
  shiny::testServer(dirty_app(edit_cfg()), {
    session$flushReact()  # Discard and close ignores its initial value
    data_root("/data/other")
    session$setInputs(wiz_close = 1)
    expect_length(got$stops, 0)
    expect_message(session$setInputs(wiz_confirm_close = 1), "settings unchanged", fixed = TRUE)
  })
  expect_match(got$shown[[1]], "Discard changes?", fixed = TRUE)
  expect_match(got$shown[[1]], 'id="wiz_confirm_close"', fixed = TRUE)
  expect_match(got$shown[[2]], "Returning to the app without saving.", fixed = TRUE)
  expect_identical(got$removed, 1L)
  expect_identical(trimws(got$stops), 'stopApp(returnValue = "restart")')
})
