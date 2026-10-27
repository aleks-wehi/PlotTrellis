# Tier 1 entry points (phase-I.3 step 25): the exported ways to start
# PlotTrellis. run_plottrellis() is the usual one; the others open one piece on
# its own. The project-settings defaults/reader and the composer app builder
# they share are in R/app.R.

# The NULL-default ("auto") settings, zrange_cap and hic_height, are numbers, but the
# importer can't infer that from a NULL default and returns their text; without this the
# composer's clean_optional_numeric() would silently read "50" as auto. A value that isn't a
# finite number is dropped (so the default stays) and named in `ignored`. Shared by
# parse_samples_txt() and the wizard's Import.
.pt_coerce_auto_overrides <- function(overrides, defaults) {
  ignored <- character(0)
  for (k in names(overrides)) if (k %in% names(defaults) && is.null(defaults[[k]])) {
    n <- suppressWarnings(as.numeric(overrides[[k]]))
    if (!is.finite(n)) ignored <- c(ignored, paste0("Ignoring '# ", k, ": ", overrides[[k]], "' (not a number); using auto."))
    overrides[[k]] <- if (is.finite(n)) n else NULL
  }
  list(overrides = overrides, ignored = ignored)
}

#' Create project settings from a samples.txt
#'
#' Reads a `samples.txt` (the file the setup wizard imports and exports), checks it the same
#' way the wizard's Save does, and writes `project_settings.json` into `project_dir`, so
#' [plottrellis_composer()] or [run_plottrellis()] opens straight into the composer, no wizard.
#' Relative `# data_root:` and `# gtf_path:` values are resolved against the folder the
#' samples.txt is in; absolute ones are used as written. Settings the file doesn't mention
#' take PlotTrellis's defaults.
#'
#' @param samples_txt Path to a samples.txt file.
#' @param project_dir Directory to write `project_settings.json` into (created if needed).
#'   Same default as [run_plottrellis()].
#' @param overwrite If `FALSE` (the default), stop instead of replacing an existing
#'   `project_settings.json` in `project_dir`.
#' @return Invisibly, the settings list that was written. Stops, writing nothing, if
#'   settings already exist and `overwrite = FALSE`, if the file has structural errors (all
#'   listed in the message), or if a required field is missing (a data root that exists, at
#'   least one sample, `gtf_path`, `orgdb`). Sample files that don't resolve give one warning
#'   and are saved anyway, like the wizard's import.
#' @seealso [run_plottrellis()], [plottrellis_composer()]
#' @examples
#' if (interactive()) {
#'   project <- tempfile()
#'   parse_samples_txt("path/to/samples.txt", project_dir = project)
#'   plottrellis_composer(project_dir = project)
#' }
#' @export
parse_samples_txt <- function(samples_txt,
                              project_dir = .pt_default_dir(),
                              overwrite = FALSE) {
  # Order: the no-clobber guard first (so a parse error can't mask it), then
  # the file, the importer's errors, the wizard's 4 checks. Nothing is
  # written or created on any error.
  settings_json <- file.path(project_dir, "project_settings.json")
  if (file.exists(settings_json) && !isTRUE(overwrite))
    stop("'", settings_json, "' already exists. Pass overwrite = TRUE to replace it, ",
         "or use a different project_dir.", call. = FALSE)
  if (!is_real_file(samples_txt)) stop("No such file: ", samples_txt, call. = FALSE)

  defaults <- .pt_default_settings()
  res <- import_samples_txt(samples_txt, defaults = defaults)
  if (length(res$errors) > 0)
    stop("Can't import '", samples_txt, "':\n", paste0("- ", res$errors, collapse = "\n"), call. = FALSE)

  coerced <- .pt_coerce_auto_overrides(res$overrides, defaults)
  for (msg in coerced$ignored) warning(msg, call. = FALSE)
  res$overrides <- coerced$overrides
  # The wizard shows these defaults in pre-filled fields; here nothing would, so a human
  # samples.txt without '# orgdb:' would silently get the mouse OrgDb.
  defaulted <- setdiff(c("orgdb", "genome_label"), names(res$overrides))
  if (length(defaulted) > 0)
    message("Not in '", basename(samples_txt), "', using the default: ",
            paste0(defaulted, " = ", vapply(defaulted, function(k) defaults[[k]], ""), collapse = ", "))

  settings <- modifyList(defaults, res$overrides)
  settings$data_dir <- res$data_root %||% ""
  # Set directly, not through modifyList(): see .pt_read_settings()'s comment
  # (an unnamed list merged into DEFAULT_SETTINGS$samples = list() is dropped).
  settings$samples <- res$samples

  problems <- wiz_save_problems(settings$data_dir, settings$samples, settings$gtf_path, settings$orgdb,
                                settings$genome_label)
  if (length(problems) > 0)
    stop("Can't create project settings from '", samples_txt, "':\n",
         paste0("- ", problems, collapse = "\n"), call. = FALSE)

  dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)
  # Same toJSON() call as the wizard's Save (null = "null", see its comment).
  writeLines(jsonlite::toJSON(settings, auto_unbox = TRUE, pretty = TRUE, null = "null"), settings_json)
  # Refresh the wizard Save's project_samples.txt copy too (feature 43): Reset and the
  # broken-settings error both point at it, so a stale one would restore the old samples.
  # As in Save, a failure here doesn't block.
  tryCatch(writeLines(export_samples_txt(c(list(data_root = settings$data_dir), settings[names(settings) != "data_dir"]),
                                         defaults = defaults),
                      file.path(project_dir, "project_samples.txt")),
           error = function(e) warning("Couldn't write project_samples.txt: ", conditionMessage(e), call. = FALSE))
  # After the write, so "Saved" is true when it's shown.
  if (length(res$warnings) > 0)
    warning(paste(c("Saved, but some sample files were not found:", res$warnings), collapse = "\n"),
            call. = FALSE)
  invisible(settings)
}

#' Open the PlotTrellis setup wizard on its own
#'
#' Opens only the setup wizard for `project_dir`: first-time setup if it has no
#' `project_settings.json`, otherwise the current settings, pre-filled for editing. Save
#' writes `project_settings.json` and stops the app; it does not continue into the composer
#' (use [run_plottrellis()] for that, or call [plottrellis_composer()] next). Blocks until
#' the wizard is saved or closed.
#'
#' @inheritParams run_plottrellis
#' @return Invisibly, the value the wizard stopped with (`"restart"` after Save or Close).
#'   Check `file.exists(file.path(project_dir, "project_settings.json"))` to see whether
#'   settings were saved.
#' @seealso [run_plottrellis()], [parse_samples_txt()]
#' @examples
#' if (interactive()) {
#'   plottrellis_wizard(project_dir = tempfile())
#' }
#' @export
plottrellis_wizard <- function(project_dir = .pt_default_dir(), ...) {
  p <- .pt_project(project_dir)
  .set_wizard_globals(p$settings_json, p$cfg, p$defaults)
  invisible(shiny::runApp(shinyApp(ui = wizard_ui(), server = wizard_server), ...))
}

#' Open the PlotTrellis composer on its own
#'
#' Opens the composer for a project that already has a `project_settings.json` (made by
#' [parse_samples_txt()] or [plottrellis_wizard()]). Unlike [run_plottrellis()] it never
#' opens the wizard: with no settings it stops with an error. Its Edit project settings,
#' Reset and Clear caches buttons stop the app rather than relaunching it. Blocks until
#' the app is stopped.
#'
#' @inheritParams run_plottrellis
#' @return Invisibly, the value the app stopped with: `"closed"` (Close app),
#'   `"edit_settings"` (Edit project settings) or `"restart"` (Reset project / Clear caches).
#' @seealso [run_plottrellis()], [parse_samples_txt()]
#' @examples
#' if (interactive()) {
#'   plottrellis_composer(project_dir = "my_project")
#' }
#' @export
plottrellis_composer <- function(project_dir = .pt_default_dir(), ...) {
  if (inherits(project_dir, "pt_project")) project_dir <- project_dir$dir
  # Checked before .pt_project() so a failed call creates no directory.
  if (!file.exists(file.path(project_dir, "project_settings.json")))
    stop("No project_settings.json in '", project_dir, "'. Create one with ",
         "parse_samples_txt(\"path/to/samples.txt\", project_dir = ...) or plottrellis_wizard(), ",
         "or launch with run_plottrellis(), which opens the wizard first when needed.", call. = FALSE)
  p <- .pt_project(project_dir)
  invisible(shiny::runApp(.pt_composer_app(p$settings_json, p$cfg), ...))
}

#' Launch PlotTrellis
#'
#' The usual way to start PlotTrellis: a 2-page Shiny app (Compose | Output) for building
#' plotgardener figures without writing code, and downloading a standalone R script that
#' reproduces each one. Opens the setup wizard when `project_dir` has no
#' `project_settings.json` yet (or after Edit project settings), otherwise the composer,
#' and switches between the two by itself after Save, Reset or Clear caches, so a first
#' run needs no re-launch. Blocks until Close app is used.
#'
#' `PlotTrellisApp()` is this function's original name, kept as an alias with the same
#' behaviour.
#'
#' @param project_dir Directory holding this project's `project_settings.json` and its
#'   cache/downloaded-script state, or a [pt_project()]. Defaults to a `PlotTrellis` folder
#'   in your home directory (created if needed); pass another folder per project.
#' @param samples_txt Optional path to a `samples.txt` to create (or refresh) the project in
#'   `project_dir` from, e.g. [pt_example()]; the same rules as [pt_project()]. Without it the
#'   app opens the project already there, or the setup wizard if there is none.
#' @param ... Passed through to [shiny::runApp()] (e.g. `launch.browser`, `port`, `host`).
#' @return Invisibly, the value the app finally stopped with (`"closed"` from Close app).
#' @seealso [plottrellis_wizard()], [plottrellis_composer()], [parse_samples_txt()]
#' @examples
#' if (interactive()) {
#'   run_plottrellis(project_dir = tempfile())
#' }
#' @export
run_plottrellis <- function(project_dir = .pt_default_dir(), samples_txt = NULL, ...) {
  if (inherits(project_dir, "pt_project")) project_dir <- project_dir$dir
  # Same rules as pt_project(): created or refreshed from samples_txt, one project per folder.
  if (!is.null(samples_txt)) .pt_import(project_dir, samples_txt)
  # Loops so the wizard-to-app transition needs no manual re-launch: the
  # settings file and WIZARD_EDIT_MODE are re-read each pass since a completed
  # wizard Save (or the main app's "Edit project settings" restart) changes
  # both on disk/env between iterations. Mirrors the restart contract
  # wizard.R already relies on (stopApp(returnValue = "restart"/"edit_settings")).
  # Each pass calls the exported launchers, so the standalone and looped paths
  # are the same code.
  repeat {
    # Phase C.1 step 1: no project_settings.json yet -> guided first-run wizard
    # instead of the normal app. Phase C.3 decision 1 extends this branch to also
    # cover re-entering the wizard to EDIT an existing settings file: the main
    # app's "Edit project settings" button sets this env var and restarts --
    # same process-restart requirement as first-run, see utils/wizard.R's header
    # comment for why. wizard_server() itself decides prefill-vs-blank-landing
    # purely from file.exists(SETTINGS_JSON), not this flag -- this flag only
    # controls which of the two cases routes here at all.
    WIZARD_EDIT_MODE <- identical(Sys.getenv("PLOTTRELLIS_EDIT_SETTINGS", "0"), "1")
    result <- if (!file.exists(file.path(project_dir, "project_settings.json")) || WIZARD_EDIT_MODE) {
      plottrellis_wizard(project_dir, ...)
    } else {
      plottrellis_composer(project_dir, ...)
    }
    if (!identical(result, "restart") && !identical(result, "edit_settings")) return(invisible(result))
  }
}

#' @rdname run_plottrellis
#' @export
PlotTrellisApp <- function(project_dir = .pt_default_dir(), ...) {
  run_plottrellis(project_dir = project_dir, ...)
}
