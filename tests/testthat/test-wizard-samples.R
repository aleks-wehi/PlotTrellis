# The wizard's Samples section (.pt_wizard_samples_server(), R/wizard-samples.R) on its own:
# the add-sample dialog, its file rows and checks, and the committed-sample list. A temp data
# folder with empty files, plain reactiveVals, no genome. The dialogs are caught, not shown.

new_data_root <- function() {
  root <- tempfile("pt-samples-")
  dir.create(file.path(root, "atac"), recursive = TRUE)
  dir.create(file.path(root, "salmon", "ctl"), recursive = TRUE)
  file.create(file.path(root, "atac", c("ctl.bw", "ko.bw")))
  file.create(file.path(root, "salmon", "ctl", "quant.sf"))
  writeLines("chr1\t1\t2\tx", file.path(root, "atac", "peaks.bed"))  # not a paired-end file
  root
}

sample <- function(condition, ...) list(condition = condition, colour = "#1B9E77", data = list(...))
bw <- function(name, path) list(type = "bw", name = name, path = path)

# A server made of nothing but the Samples block over `root`, starting with `samples`, with
# the dialog's state as wizard_server() keeps it.
samples_app <- function(root, samples = list()) {
  function(input, output, session) {
    wiz_samples <- reactiveVal(samples)
    data_root   <- reactive(root)
    modal_instances    <- reactiveVal(list())
    modal_row_id       <- reactiveVal(0)
    modal_edit_name    <- reactiveVal(NULL)
    modal_edit_index   <- function() {
      nm <- modal_edit_name()
      if (is.null(nm)) NULL else match(nm, vapply(wiz_samples(), function(s) s$condition, character(1)))
    }
    modal_row_errors   <- reactiveVal(list())
    modal_field_errors <- reactiveVal(list())
    modal_fail_count   <- reactiveVal(0)
    res <- .pt_wizard_samples_server(input, output, wiz_samples, data_root, modal_instances, modal_row_id,
                                     modal_edit_name, modal_edit_index, modal_row_errors, modal_field_errors,
                                     modal_fail_count)
  }
}

# Catches the dialogs, notifications and removeModal() calls the code under test makes, for
# the rest of the calling test.
catch_dialogs <- function(env = parent.frame()) {
  got <- new.env()
  got$shown <- character(0); got$notes <- character(0); got$removed <- 0L
  local_mocked_bindings(showModal = function(ui, ...) got$shown <- c(got$shown, as.character(ui)),
                        removeModal = function(...) got$removed <- got$removed + 1L,
                        showNotification = function(ui, ...) got$notes <- c(got$notes, as.character(ui)),
                        .env = env)
  got
}

conditions <- function(samples) vapply(samples, `[[`, "", "condition")

test_that("prefill rows: the deduped union of type:name pairs, no quantdir, ids after last_id", {
  samples <- list(sample("A", bw("ATAC", "a.bw"), list(type = "quantdir", name = "Quant", path = "q")),
                  sample("B", bw("ATAC", "b.bw"), list(type = "hic", name = "HiC", path = "b.hic")))
  rows <- .pt_wizard_prefill_rows(samples, keep_path = FALSE, last_id = 4)
  expect_identical(vapply(rows, `[[`, 0, "id"), c(5, 6))
  expect_identical(vapply(rows, `[[`, "", "type"), c("bw", "hic"))
  expect_identical(vapply(rows, `[[`, "", "path"), c("", ""))
  kept <- .pt_wizard_prefill_rows(samples[2], keep_path = TRUE, last_id = 0)
  expect_identical(vapply(kept, `[[`, "", "path"), c("b.bw", "b.hic"))
  expect_length(.pt_wizard_prefill_rows(list(), keep_path = FALSE, last_id = 0), 0)
})

test_that("Add a sample: a bw row and a salmon folder become a committed sample", {
  root <- new_data_root(); on.exit(unlink(root, recursive = TRUE), add = TRUE)
  got <- catch_dialogs()
  shiny::testServer(samples_app(root), {
    session$flushReact()
    session$setInputs(wiz_add_sample = 1)
    session$setInputs(wiz_add_bw = 1)
    expect_identical(modal_instances()[[1]][c("id", "type")], list(id = 1, type = "bw"))
    session$setInputs(wiz_sample_name = " Ctl ", wiz_sample_colour = "#123456",
                      wiz_inst_name_1 = "ATAC", wiz_inst_path_1 = "atac/ctl.bw",
                      wiz_quantdir_name = "", wiz_quantdir_path = "salmon/ctl")
    session$setInputs(wiz_confirm_sample = 1)
    expect_identical(wiz_samples(), list(list(condition = "Ctl", colour = "#123456", data = list(
      list(type = "bw", name = "ATAC", path = "atac/ctl.bw"),
      list(type = "quantdir", name = "Quant", path = "salmon/ctl")))))
    expect_null(modal_edit_name())
  })
  expect_length(got$shown, 1)
  expect_match(got$shown[[1]], "Add a sample", fixed = TRUE)
  expect_match(got$shown[[1]], "salmon/ctl", fixed = TRUE)  # the salmon folder is offered
  expect_identical(got$removed, 1L)
})

test_that("a second sample's dialog prefills the first's rows with new ids and blank paths", {
  root <- new_data_root(); on.exit(unlink(root, recursive = TRUE), add = TRUE)
  first <- sample("Ctl", bw("ATAC", "atac/ctl.bw"), list(type = "quantdir", name = "Salmon", path = "salmon/ctl"))
  got <- catch_dialogs()
  shiny::testServer(samples_app(root, list(first)), {
    session$flushReact()
    session$setInputs(wiz_add_sample = 1)
    expect_identical(modal_instances(), list(list(id = 1, type = "bw", name = "ATAC", path = "")))
    session$setInputs(wiz_add_sample = 2)  # opened again: ids never reused
    expect_identical(modal_instances()[[1]]$id, 2)
    html <- paste(as.character(output$wiz_instance_rows$html), collapse = "")
    expect_match(html, 'id="wiz_inst_name_2"', fixed = TRUE)
  })
  # The quantdir's track name carries over, not its path.
  expect_match(got$shown[[1]], 'id="wiz_quantdir_name"[^>]*value="Salmon"')
})

test_that("each bad row is named in the re-shown dialog and nothing is committed", {
  root <- new_data_root(); on.exit(unlink(root, recursive = TRUE), add = TRUE)
  got <- catch_dialogs()
  shiny::testServer(samples_app(root, list(sample("Ctl"))), {
    session$flushReact()
    session$setInputs(wiz_add_sample = 1)
    for (i in 1:5) session$setInputs(wiz_add_bw = i)
    session$setInputs(wiz_add_paired_end = 1)
    session$setInputs(wiz_sample_name = "Ctl",
                      wiz_inst_name_1 = "",          wiz_inst_path_1 = "atac/ctl.bw",
                      wiz_inst_name_2 = "a·b",  wiz_inst_path_2 = "atac/ctl.bw",
                      wiz_inst_name_3 = "ATAC",      wiz_inst_path_3 = "atac/ctl.bw",
                      wiz_inst_name_4 = "ATAC",      wiz_inst_path_4 = "atac/ko.bw",
                      wiz_inst_name_5 = "Gone",      wiz_inst_path_5 = "atac/none.bw",
                      wiz_inst_name_6 = "Loops",     wiz_inst_path_6 = "atac/peaks.bed",
                      wiz_quantdir_path = "atac")
    session$setInputs(wiz_confirm_sample = 1)
    expect_identical(conditions(wiz_samples()), "Ctl")
    expect_identical(modal_field_errors(), list(quantdir = "No quant.sf found under: atac"))
    rows <- paste(as.character(output$wiz_instance_rows$html), collapse = "")
    expect_identical(lengths(regmatches(rows, gregexpr("wiz-row wiz-row-error", rows, fixed = TRUE))), 5L)
  })
  expect_length(got$shown, 2)
  err <- got$shown[[2]]
  for (msg in c("A condition named 'Ctl' already exists.",
                "A bw instance has a path but no name.",
                "contains a bare",
                "Duplicate bw instance name 'ATAC'.",
                "Path for bw:Gone does not resolve to a file under the data root: atac/none.bw",
                "Paired-end file for Loops (atac/peaks.bed) has 4 columns",
                "No quant.sf found under: atac"))
    expect_match(err, msg, fixed = TRUE)
  expect_identical(got$removed, 0L)
})

test_that("a blank condition is required; a row with no path is dropped, not an error", {
  root <- new_data_root(); on.exit(unlink(root, recursive = TRUE), add = TRUE)
  got <- catch_dialogs()
  shiny::testServer(samples_app(root), {
    session$flushReact()
    session$setInputs(wiz_add_sample = 1, wiz_add_hic = 1)
    session$setInputs(wiz_sample_name = "", wiz_confirm_sample = 1)
    session$setInputs(wiz_sample_name = "KO", wiz_confirm_sample = 2)
    expect_identical(wiz_samples(), list(list(condition = "KO", colour = WIZARD_DEFAULT_COLOUR, data = list())))
  })
  expect_match(got$shown[[2]], "Condition is required.", fixed = TRUE)
})

test_that("Edit opens the sample's own rows and saves it in place, under its own name", {
  root <- new_data_root(); on.exit(unlink(root, recursive = TRUE), add = TRUE)
  start <- list(sample("A", bw("ATAC", "atac/ctl.bw")), sample("B", bw("ATAC", "atac/ko.bw")))
  got <- catch_dialogs()
  shiny::testServer(samples_app(root, start), {
    session$flushReact()
    session$setInputs(wiz_edit_sample = "missing")  # no such sample: nothing opens
    session$setInputs(wiz_edit_sample = "A")
    expect_identical(modal_edit_name(), "A")
    expect_identical(modal_edit_index(), 1L)
    expect_identical(modal_instances(), list(list(id = 1, type = "bw", name = "ATAC", path = "atac/ctl.bw")))
    session$setInputs(wiz_sample_name = "A", wiz_sample_colour = "#000000",
                      wiz_inst_name_1 = "ATAC", wiz_inst_path_1 = "atac/ko.bw", wiz_confirm_sample = 1)
    expect_identical(conditions(wiz_samples()), c("A", "B"))
    expect_identical(wiz_samples()[[1]]$data[[1]]$path, "atac/ko.bw")
    expect_identical(wiz_samples()[[1]]$colour, "#000000")
  })
  expect_length(got$shown, 1)
  expect_match(got$shown[[1]], "Edit sample", fixed = TRUE)
  expect_match(got$shown[[1]], "Save changes", fixed = TRUE)
})

test_that("a sample removed while its Edit dialog is open is not saved back", {
  root <- new_data_root(); on.exit(unlink(root, recursive = TRUE), add = TRUE)
  got <- catch_dialogs()
  shiny::testServer(samples_app(root, list(sample("A"), sample("B"))), {
    session$flushReact()
    session$setInputs(wiz_edit_sample = "A")
    session$setInputs(wiz_remove_sample = "A")
    expect_identical(modal_edit_index(), NA_integer_)
    session$setInputs(wiz_sample_name = "A", wiz_confirm_sample = 1)
    expect_identical(conditions(wiz_samples()), "B")
    expect_null(modal_edit_name())
  })
  expect_match(got$notes, "Sample 'A' was removed", fixed = TRUE)
  expect_identical(got$removed, 1L)
})

test_that("a row's remove button takes that row out of the dialog", {
  root <- new_data_root(); on.exit(unlink(root, recursive = TRUE), add = TRUE)
  catch_dialogs()
  shiny::testServer(samples_app(root), {
    session$flushReact()
    session$setInputs(wiz_add_sample = 1)
    session$setInputs(wiz_add_bw = 1); session$setInputs(wiz_add_bed = 1)
    session$setInputs(wiz_remove_row = 1)
    expect_identical(vapply(modal_instances(), `[[`, "", "type"), "bed")
  })
})

test_that("the sample list shows nothing and the dialog's rows say when they are empty", {
  root <- new_data_root(); on.exit(unlink(root, recursive = TRUE), add = TRUE)
  catch_dialogs()
  shiny::testServer(samples_app(root), {
    session$flushReact()
    expect_identical(paste(as.character(output$wiz_samples_ui$html), collapse = ""), "")
    expect_match(paste(as.character(output$wiz_instance_rows$html)), "No files yet", fixed = TRUE)
    wiz_samples(list(sample("Ctl", bw("ATAC", "atac/ctl.bw"))))
    session$flushReact()
    html <- paste(as.character(output$wiz_samples_ui$html), collapse = "")
    expect_match(html, "bw:ATAC", fixed = TRUE)
    expect_match(html, 'id="wiz_rm_sample_btn_1"', fixed = TRUE)
  })
})
