# The composer's Figure settings block (R/composer-settings.R) on its own: render settings
# made from the defaults, a temporary project_settings.json, no genome. Notifications are
# recorded. (test-composer-settings.R drives the same block through the whole app.)

BLOCK_RENDER_CFG <- c(.pt_default_settings(),
                      list(orgdb_pkg = "org.Fake.eg.db", txdb_path = "/no/such/txdb.sqlite"))

# A project folder with a project_settings.json holding core keys only; returns its path.
new_settings_json <- function() {
  proj <- tempfile("pt-settings-block-"); dir.create(proj)
  json <- file.path(proj, "project_settings.json")
  jsonlite::write_json(list(gtf_path = "/data/genome.gtf", orgdb = "org.Mm.eg.db", genome_label = "mm10",
                            data_dir = proj, samples = list()), json, auto_unbox = TRUE)
  json
}

settings_app <- function(json, render_cfg = BLOCK_RENDER_CFG) function(input, output, session) {
  # By name: the block takes no session.
  fs <- .pt_figure_settings_server(input = input, output = output, render_cfg = render_cfg, settings_json = json)
}

# Runs `code` against the settings block (it sees session, output and fs); returns the
# notifications shown, each list(ui, type, duration).
drive_settings <- function(json, code) {
  code <- rlang::enquo(code)
  rec <- new.env(); rec$notes <- list()
  local_mocked_bindings(showNotification = function(ui, action = NULL, duration = 5, closeButton = TRUE, id = NULL,
                                                    type = "default", session = NULL)
    rec$notes[[length(rec$notes) + 1]] <- list(ui = ui, type = type, duration = duration))
  rlang::inject(shiny::testServer(settings_app(json), !!code))
  rec$notes
}

test_that("the block starts with live = saved = the render settings, and nothing unsaved", {
  json <- new_settings_json(); on.exit(unlink(dirname(json), recursive = TRUE), add = TRUE)
  drive_settings(json, {
    expect_identical(fs$saved_cfg(), BLOCK_RENDER_CFG)
    for (k in fs$display_keys) expect_identical(fs$live_cfg()[[k]], BLOCK_RENDER_CFG[[k]], info = k)
    expect_false(fs$settings_dirty())
    expect_identical(fs$settings_ui_gen(), 0)
    expect_null(output$settings_dirty_ui)
  })
})

test_that("display_keys are every setting except the core ones", {
  json <- new_settings_json(); on.exit(unlink(dirname(json), recursive = TRUE), add = TRUE)
  drive_settings(json, {
    expect_identical(fs$display_keys, setdiff(names(.pt_default_settings()), .PT_CORE_SETTINGS))
    expect_false(any(.PT_CORE_SETTINGS %in% fs$display_keys))
    expect_true(all(c("fontsize", "zrange_cap", "gene_flank_pct") %in% fs$display_keys))
  })
})

test_that("a panel input changes the live settings only, and marks them unsaved", {
  json <- new_settings_json(); on.exit(unlink(dirname(json), recursive = TRUE), add = TRUE)
  drive_settings(json, {
    session$setInputs(s_fontsize = 12, s_font = "  Helvetica ")
    expect_identical(fs$live_cfg()$fontsize, 12)
    expect_identical(fs$live_cfg()$font, "Helvetica")
    expect_identical(fs$saved_cfg()$fontsize, 8L)
    expect_true(fs$settings_dirty())
    expect_identical(ui_html(output$settings_dirty_ui),
                     '<span class="text-warning">Unsaved changes (this session only).</span>')
  })
})

test_that("a value equal to the saved one is not unsaved, whatever its numeric type", {
  json <- new_settings_json(); on.exit(unlink(dirname(json), recursive = TRUE), add = TRUE)
  drive_settings(json, {
    session$setInputs(s_fontsize = 8, s_page_w = 9L, s_zrange_cap = "")
    expect_false(fs$settings_dirty())
    # A blank optional number stays an explicit NULL ("auto") in the live settings.
    expect_true("zrange_cap" %in% names(fs$live_cfg()))
    expect_null(fs$live_cfg()$zrange_cap)
    session$setInputs(s_zrange_cap = 50)
    expect_true(fs$settings_dirty())
  })
})

test_that("the live settings keep the render-only keys, and core keys ignore panel inputs", {
  json <- new_settings_json(); on.exit(unlink(dirname(json), recursive = TRUE), add = TRUE)
  drive_settings(json, {
    session$setInputs(s_gtf_path = "/other.gtf", s_orgdb_pkg = "org.Other.eg.db", s_fontsize = 10)
    expect_identical(fs$live_cfg()$gtf_path, "")
    expect_identical(fs$live_cfg()$orgdb_pkg, "org.Fake.eg.db")
    expect_identical(fs$live_cfg()$txdb_path, "/no/such/txdb.sqlite")
    expect_identical(setdiff(names(fs$live_cfg()), names(BLOCK_RENDER_CFG)), character(0))
  })
})

test_that("Save writes the display settings into the file, keeping its core keys", {
  json <- new_settings_json(); on.exit(unlink(dirname(json), recursive = TRUE), add = TRUE)
  notes <- drive_settings(json, {
    session$setInputs(s_fontsize = 12, s_hic_norm = "VC", s_zrange_cap = "")
    session$setInputs(save_settings = 1)
    expect_identical(fs$saved_cfg()$fontsize, 12)
    expect_false(fs$settings_dirty())
    expect_null(output$settings_dirty_ui)
  })
  expect_identical(notes, list(list(ui = "Saved as the project default.", type = "message", duration = 5)))
  saved <- jsonlite::fromJSON(json, simplifyDataFrame = FALSE)
  expect_identical(saved$fontsize, 12L)
  expect_identical(saved$hic_norm, "VC")
  expect_null(saved$zrange_cap)
  expect_true("zrange_cap" %in% names(saved))
  # Core keys come from the file, not from the render settings (whose gtf_path is "").
  expect_identical(saved$gtf_path, "/data/genome.gtf")
  expect_identical(saved$genome_label, "mm10")
  # Render-only keys are never written.
  expect_false(any(c("orgdb_pkg", "txdb_path") %in% names(saved)))
  # project_samples.txt is refreshed beside it.
  expect_true(file.exists(file.path(dirname(json), "project_samples.txt")))
})

test_that("Save with the settings file gone shows the exact error and saves nothing", {
  json <- new_settings_json(); on.exit(unlink(dirname(json), recursive = TRUE), add = TRUE)
  unlink(json)
  notes <- drive_settings(json, {
    session$setInputs(s_fontsize = 12, save_settings = 1)
    expect_identical(fs$saved_cfg()$fontsize, 8L)
    expect_true(fs$settings_dirty())
  })
  expect_identical(notes, list(list(ui = paste0("Couldn't save the settings: ", json, " no longer exists; restart the app."),
                                    type = "error", duration = NULL)))
  expect_false(file.exists(json))
})

test_that("Save with an unreadable settings file shows the read error and leaves the file as it was", {
  json <- new_settings_json(); on.exit(unlink(dirname(json), recursive = TRUE), add = TRUE)
  writeLines("{not json", json)
  notes <- drive_settings(json, {
    session$setInputs(s_fontsize = 12, save_settings = 1)
    expect_identical(fs$saved_cfg()$fontsize, 8L)
  })
  expect_length(notes, 1)
  expect_identical(notes[[1]]$type, "error")
  expect_null(notes[[1]]$duration)
  expect_match(notes[[1]]$ui, paste0("Couldn't save the settings: Can't read ", json, " ("), fixed = TRUE)
  expect_identical(readLines(json), "{not json")
})

test_that("the panel shows the saved values, and only Revert re-renders it", {
  json <- new_settings_json(); on.exit(unlink(dirname(json), recursive = TRUE), add = TRUE)
  drive_settings(json, {
    html <- ui_html(output$display_settings_ui)
    expect_match(html, '<input id="s_fontsize" type="number" class="[^"]*" value="8"')
    expect_match(html, "Comparison (diff) tracks", fixed = TRUE)
    session$setInputs(s_fontsize = 12, save_settings = 1)
    # Save does not re-render the panel (the browser already shows 12).
    expect_match(ui_html(output$display_settings_ui), '<input id="s_fontsize" type="number" class="[^"]*" value="8"')
    session$setInputs(revert_settings = 1)
    expect_identical(fs$settings_ui_gen(), 1)
    expect_match(ui_html(output$display_settings_ui), '<input id="s_fontsize" type="number" class="[^"]*" value="12"')
  })
})

test_that("Revert re-renders from the saved settings, dropping unsaved panel values", {
  json <- new_settings_json(); on.exit(unlink(dirname(json), recursive = TRUE), add = TRUE)
  drive_settings(json, {
    session$setInputs(s_fontsize = 14)
    session$setInputs(revert_settings = 1)
    session$setInputs(revert_settings = 2)
    expect_identical(fs$settings_ui_gen(), 2)
    expect_match(ui_html(output$display_settings_ui), '<input id="s_fontsize" type="number" class="[^"]*" value="8"')
    # The live value follows once the browser reports the re-rendered input.
    expect_true(fs$settings_dirty())
    session$setInputs(s_fontsize = 8)
    expect_false(fs$settings_dirty())
  })
})
