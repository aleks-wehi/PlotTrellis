# The app's project actions: the dialogs behind the page's Edit project
# settings button, Reset project configuration and Clear all caches included.
# A top-level function that takes the session objects and only the state it
# reads, keeps every input id (the page, the browser drivers and the tests
# address them by id), and tests on its own with no genome.

# The warning shown in the Edit project settings and Clear caches confirms
# while the Figure settings have unsaved changes; NULL (nothing shown)
# otherwise.
.pt_unsaved_settings_note <- function(dirty) {
  if (isTRUE(dirty))
    p(class = "text-warning",
      strong("Your Figure settings changes are not saved and will be lost."),
      " Use Save as project default first to keep them.")
}

# `unsaved_settings_note` is the server's zero-argument function returning
# that warning; the three paths are the project's.
.pt_project_actions_server <- function(input, session, unsaved_settings_note,
                                       project_dir, settings_json, cache_dir) {

  # ── Settings (re-enter wizard) ──────────────────────────────────────────────
  # Editing settings stops this app and restarts into the wizard, a separate app
  # (run_plottrellis()'s loop), since everything here is built from the settings
  # the wizard edits -- confirm first so this isn't a silent surprise
  # mid-session, rather than restarting immediately on click.
  # The dialog also shows the project folder (everything PlotTrellis writes is
  # in it) and, in a collapsed Advanced section, Clear all caches and Reset
  # project configuration.
  observeEvent(input$open_settings, {
    showModal(modalDialog(
      title = "Project settings", easyClose = FALSE,
      p("Project folder: ", tags$code(project_dir)),
      # Light or dark: the page's script (inst/www/main/main.js) sets this switch
      # to the mode on screen when the dialog shows, and acts on it.
      div(class = "pt-setrow", span("Appearance"),
          div(class = "flip", id = "modeFlip", role = "switch", tabindex = "0",
              `aria-checked` = "false", `aria-label` = "Dark mode",
              span("Light"), span("Dark"))),
      p("Edit settings restarts the app into the setup wizard, pre-filled with your current settings. Any in-progress figure/selection in this session will be lost."),
      unsaved_settings_note(),
      tags$details(id = "settings_advanced",
        tags$summary(class = "btn btn-default btn-sm", "Advanced"),
        p("Everything this project keeps (settings and caches) is in its project folder. Your data files are never changed."),
        actionButton("clear_caches_btn", "Clear all caches...",
                     class = "btn-warning btn-sm pt-once"),
        actionButton("reset_project_btn", "Reset project configuration...",
                     class = "btn-danger btn-sm pt-once")
      ),
      footer = tagList(modalButton("Close"),
                       actionButton("confirm_edit_settings", "Edit settings",
                                    class = "btn-primary pt-once"))
    ))
  })

  observeEvent(input$confirm_edit_settings, {
    # Same "show feedback, then delay stopApp()" pattern as
    # confirm_reset_project below -- an immediate stopApp() can kill the
    # process before the modal's removal/replacement actually
    # reaches the browser.
    showModal(modalDialog(
      title = "Restarting into the setup wizard", easyClose = FALSE,
      footer = NULL,
      p(strong("You can close this tab now."),
        " The setup wizard opens in a new tab from RStudio in a few seconds. ",
        "If nothing opens, go back to RStudio: the app's address is in the Console.")
    ))
    later::later(function() stopApp(returnValue = "edit_settings"), 0.2)
  }, ignoreInit = TRUE)

  observeEvent(input$reset_project_btn, {
    showModal(modalDialog(
      title = "Reset project configuration?", easyClose = FALSE,
      p("This deletes project_settings.json and restarts the app into first-time setup, which offers to restore it from project_samples.txt (kept). Your own samples.txt is not touched."),
      footer = tagList(modalButton("Cancel"),
                       actionButton("confirm_reset_project",
                                    "Delete and restart",
                                    class = "btn-danger pt-once"))
    ))
  })

  observeEvent(input$confirm_reset_project, {
    # file.remove() gives its reason ("Permission denied") only as a warning
    why <- NULL
    if (file.exists(settings_json))
      withCallingHandlers(file.remove(settings_json), warning = function(w) {
        why <<- sub("^cannot remove file .*, reason '(.*)'$", "\\1",
                    conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    # A folder that can't be written keeps the file: say so, and don't restart
    # into it.
    if (file.exists(settings_json)) {
      showModal(modalDialog(
        title = "Nothing was reset", easyClose = TRUE,
        footer = modalButton("Close"),
        p("Couldn't delete ", settings_json,
          if (!is.null(why)) paste0(" (", why, ")"), ". ",
          "Check that you can write to ", dirname(settings_json),
          ", then try again.")
      ))
      return()
    }
    showModal(modalDialog(
      title = "Project configuration reset", easyClose = FALSE, footer = NULL,
      p("project_settings.json has been deleted."),
      p("Restarting into first-time setup."),
      p(strong("You can close this tab now."),
        " The app reopens in a new tab from RStudio in a few seconds. ",
        "If nothing opens, go back to RStudio: the app's address is in the Console.")
    ))
    later::later(function() stopApp(returnValue = "restart"), 0.2)
  }, ignoreInit = TRUE)

  # ── Clear all caches ────────────────────────────────────────────────────────
  # Deletes the project's cache/ folder and nothing else (never files beside the
  # GTF or in data folders); project_settings.json is left alone, so the restart
  # goes straight back into the app, and each cache rebuilds on next use.
  observeEvent(input$clear_caches_btn, {
    showModal(modalDialog(
      title = "Clear all caches?", easyClose = FALSE,
      p("This deletes this project's cache folder, ", tags$code(cache_dir),
        ":"),
      tags$ul(
        tags$li("Genome files built from your GTF (rebuilding takes 1-3 minutes for a full annotation)"),
        tags$li("Chromosome lengths (fetched again from UCSC, which needs an internet connection)"),
        tags$li("Transcripts shown per gene, UniProt canonical isoforms, protein domains")
      ),
      p("Project settings, your samples.txt and your data files are not touched. Everything above rebuilds automatically the next time it's needed."),
      unsaved_settings_note(),
      footer = tagList(modalButton("Cancel"),
                       actionButton("confirm_clear_caches",
                                    "Delete and restart",
                                    class = "btn-warning pt-once"))
    ))
  })

  observeEvent(input$confirm_clear_caches, {
    # Close the open genome databases first: on a network drive a deleted file
    # still open stays behind as a .nfs file.
    .pt_close_txdbs()
    unlink(cache_dir, recursive = TRUE)
    # unlink() doesn't say what it left. A .nfs file is a network drive's
    # stand-in for a file still open in this R session, and goes when R exits:
    # not a failure.
    left <- if (dir.exists(cache_dir))
      list.files(cache_dir, recursive = TRUE, all.files = TRUE, no.. = TRUE,
                 include.dirs = TRUE) else character(0)
    left <- left[!grepl("(^|/)\\.nfs", left)]
    left <- left[!vapply(file.path(cache_dir, left), function(f) {
      isTRUE(file.info(f)$isdir) &&
        length(list.files(f, all.files = TRUE, no.. = TRUE,
                          recursive = TRUE)) > 0 &&
        all(grepl("(^|/)\\.nfs", list.files(f, all.files = TRUE, no.. = TRUE,
                                            recursive = TRUE)))
    }, logical(1))]

    showModal(modalDialog(
      title = if (length(left)) "Caches partly cleared" else "Caches cleared",
      easyClose = FALSE, footer = NULL,
      if (length(left))
        p(length(left), " item", if (length(left) > 1) "s", " in ", cache_dir,
          " could not be deleted (",
          paste(utils::head(left, 3), collapse = ", "),
          if (length(left) > 3) ", ...", "); ",
          "check you can write to that folder, or delete ",
          if (length(left) > 1) "them" else "it", " by hand.")
      else p("All caches have been deleted."),
      p("Restarting with your current project settings. Rebuilding the genome files can take 1-3 minutes for a full annotation."),
      p(strong("You can close this tab now."),
        " The app reopens in a new tab from RStudio in a few seconds. ",
        "If nothing opens, go back to RStudio: the app's address is in the Console.")
    ))
    later::later(function() stopApp(returnValue = "restart"), 0.2)
  }, ignoreInit = TRUE)
  invisible(NULL)
}
