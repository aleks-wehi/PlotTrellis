# The setup wizard's Samples section: one block of wizard_server(). The
# add-sample dialog (add and edit), its file rows and their checks, and the
# committed-sample list with its Edit and Remove buttons. `wiz_samples` is the
# committed samples, `data_root` the Data folder reactive, and the modal_*
# arguments the dialog's state (see wizard_server()). Returns
# list(add_sample_modal =, wiz_prefill_rows =, add_modal_row =).

# `samples` is a list to prefill rows from: keep_path=TRUE for one sample being
# edited, FALSE for the deduped union (first-seen order) of every committed
# sample's type:name pairs, so a pair missing from the latest sample isn't lost;
# quantdir isn't included here, see add_sample_modal(). Row ids continue from
# `last_id`, the highest id given out so far.
.pt_wizard_prefill_rows <- function(samples, keep_path, last_id) {
  rows <- list()
  seen <- character(0)
  for (s in samples) {
    for (d in s$data) {
      if (identical(d$type, "quantdir")) next
      # ":" isn't valid in a WIZARD_DATA_TYPES key, so this can't collide across
      # types
      key <- paste0(d$type, ":", d$name)
      if (key %in% seen) next
      seen <- c(seen, key)
      rid <- last_id + length(rows) + 1
      rows[[length(rows) + 1]] <- list(id = rid, type = d$type, name = d$name,
                                        path = if (keep_path) d$path else "")
    }
  }
  rows
}

# The add-sample dialog, its rows and checks, and the committed-sample list.
.pt_wizard_samples_server <- function(input, output, wiz_samples, data_root,
                                      modal_instances, modal_row_id,
                                      modal_edit_name, modal_edit_index,
                                      modal_row_errors, modal_field_errors,
                                      modal_fail_count) {
  # ── Add-sample modal ────────────────────────────────────────────────────────

  # `error` NULL means a fresh modal open (fields blank/default); non-NULL means
  # a re-display after a failed wiz_confirm_sample validation (fields seeded
  # from current inputs instead, since Shiny never clears an unmounted input, so
  # `error`'s presence is the only way to tell the two apart).
  add_sample_modal <- function(error = NULL) {
    dr    <- data_root()
    fresh <- is.null(error)
    seed  <- function(input_id, default)
      if (fresh) default else (isolate(input[[input_id]]) %||% default)
    field_errors <- if (fresh) list() else modal_field_errors()
    # Editing an existing sample (modal_edit_index() set) prefills name/
    # colour/quantdir from that sample. Adding a brand-new sample prefills only
    # the quantdir NAME (not path): scans every already-committed sample for a
    # quantdir entry, most recent match wins, rather than only the very last
    # sample -- so a naming convention set early on isn't lost just because the
    # most recent sample happens to have no quantdir at all. Nothing
    # sample-specific (condition name, colour, actual paths) is prefilled into a
    # new sample.
    edit_idx       <- modal_edit_index()
    editing_sample <- if (!is.null(edit_idx) && !is.na(edit_idx))
      wiz_samples()[[edit_idx]] else NULL
    carry_quantdir <- if (is.null(editing_sample)) {
      quantdirs <- Filter(Negate(is.null), lapply(wiz_samples(), function(s)
        Find(function(d) identical(d$type, "quantdir"), s$data)))
      if (length(quantdirs) > 0) quantdirs[[length(quantdirs)]] else NULL
    } else NULL
    quantdir_inst  <- if (!is.null(editing_sample))
      Find(function(d) identical(d$type, "quantdir"),
           editing_sample$data) else carry_quantdir
    default_name    <- editing_sample$condition %||% ""
    default_colour  <- editing_sample$colour %||% WIZARD_DEFAULT_COLOUR
    default_qname   <- quantdir_inst$name %||% "Quant"
    default_qpath   <- if (!is.null(editing_sample))
      (quantdir_inst$path %||% "") else ""
    # Explicit "" choice: with only real file paths in `choices`,
    # selected="" matches no <option>, and
    # shiny's underlying native <select> then falls back to auto-selecting
    # the first real option before selectize.js upgrades it -- so a field the
    # user never touched silently pointed at whatever file happened to be
    # first. A real blank option gives selected="" something to bind to.
    quantdir_choices <- c("", wizard_scoped_choices(dr, "quantdir"))
    modalDialog(
      title = if (is.null(edit_idx)) "Add a sample" else "Edit sample",
      size = "l", easyClose = FALSE,
      if (!is.null(error)) wiz_error_list("Please fix before adding:", error),
      fluidRow(
        column(6, textInput("wiz_sample_name",
                            wiz_help("Condition",
                                     "The name this sample carries in the figure; it must be unique."),
                            value = seed("wiz_sample_name", default_name),
                            placeholder = "e.g. Control")),
        column(6, tags$label(wiz_help("Colour",
                                      "The colour of this sample's coverage and junction tracks.")),
               br(),
               # shiny.js has no input binding for type="color"
               # (.pt_colour_input() works around the same gap the same way) --
               # without the oninput handler below, input$wiz_sample_colour
               # never updates from its initial value.
               tags$input(type = "color", id = "wiz_sample_colour",
                          value = seed("wiz_sample_colour", default_colour),
                          style = "width:60px; height:34px; padding:2px; border:1px solid #ccc; border-radius:4px; cursor:pointer; display:block;",
                          oninput = "Shiny.setInputValue('wiz_sample_colour', this.value, {priority: 'event'});"),
               # Report the shown colour as soon as the modal renders:
               # oninput alone leaves input$wiz_sample_colour at whatever the
               # previous modal set, so an untouched swatch saved that colour.
               tags$script(sprintf(
                 "Shiny.setInputValue('wiz_sample_colour', %s);",
                 jsonlite::toJSON(seed("wiz_sample_colour", default_colour),
                                  auto_unbox = TRUE)
               )))
      ),
      hr(),
      tags$label(wiz_help("Files for this sample",
                          "One row per file; use the same track name for the same data in every sample.")),
      p(class = "text-muted", style = "margin-top:-2px;",
        "Pick the kind of file, then its path from the list (every matching file under the data folder is listed; typing a path works too)."),
      div(class = "wiz-file-kinds",
        actionButton("wiz_add_bw",
                     paste0("+ ", WIZ_KIND_LABELS[["bw"]]),
                     class = "pt-once btn-sm btn-default"),
        actionButton("wiz_add_hic",
                     paste0("+ ", WIZ_KIND_LABELS[["hic"]]),
                     class = "pt-once btn-sm btn-default"),
        actionButton("wiz_add_bed",
                     paste0("+ ", WIZ_KIND_LABELS[["bed"]]),
                     class = "pt-once btn-sm btn-default"),
        actionButton("wiz_add_paired_end",
                     paste0("+ ", WIZ_KIND_LABELS[["paired-end"]]),
                     class = "pt-once btn-sm btn-default")
      ),
      uiOutput("wiz_instance_rows"),
      hr(),
      tags$label(wiz_help("Transcript quantification (optional)",
                          "This sample's salmon folder (the one holding quant.sf).")),
      fluidRow(class = if (is.null(field_errors$quantdir))
                         "wiz-row" else "wiz-row wiz-row-error",
        column(2, div(class = "wiz-kind", "Salmon folder")),
        column(4, textInput("wiz_quantdir_name", NULL,
                            value = seed("wiz_quantdir_name", default_qname),
                            placeholder = "track name")),
        column(6, selectizeInput(
          "wiz_quantdir_path", NULL,
          choices = unique(c(quantdir_choices,
                             seed("wiz_quantdir_path", default_qpath))),
          selected = seed("wiz_quantdir_path", default_qpath),
          options = list(create = TRUE,
                         placeholder = "folder with quant.sf, relative to the data folder")
        )),
        if (!is.null(field_errors$quantdir))
          column(12, tags$div(class = "wiz-error",
                              style = "font-size:13px; margin-top:4px;",
                              field_errors$quantdir))
      ),
      footer = tagList(modalButton("Cancel"),
                       actionButton("wiz_confirm_sample",
                                    if (is.null(edit_idx))
                                      "Add sample" else "Save changes",
                                    class = "pt-once btn-primary"))
    )
  }

  # The rows to prefill the dialog with (.pt_wizard_prefill_rows()), each with a
  # new row id.
  wiz_prefill_rows <- function(samples, keep_path) {
    rows <- .pt_wizard_prefill_rows(samples, keep_path, modal_row_id())
    modal_row_id(modal_row_id() + length(rows))
    rows
  }

  observeEvent(input$wiz_add_sample, {
    modal_edit_name(NULL)
    modal_instances(wiz_prefill_rows(wiz_samples(), keep_path = FALSE))
    modal_row_errors(list())
    modal_field_errors(list())
    # modal_row_id is NOT reset here: its ids are
    # input-element ids (wiz_inst_name_<id>/wiz_inst_path_<id>) that live
    # for the whole Shiny session even after their row is removed from the
    # DOM -- Shiny never clears input$wiz_inst_name_1 just because the
    # element disappeared. Resetting to 0 on every "+ Add sample" click meant
    # a fresh row on sample 2 reused id 1 and silently inherited sample 1's
    # already-committed values. Keeping the counter monotonic for the whole
    # session means an id is retired for good once used, so this can't recur.
    showModal(add_sample_modal())
  })

  observeEvent(input$wiz_edit_sample, {
    # By condition name, like Remove below: a row number can be stale after a
    # quick Remove.
    idx <- match(input$wiz_edit_sample,
                 vapply(wiz_samples(), function(s) s$condition, character(1)))
    if (is.na(idx)) return()
    modal_edit_name(input$wiz_edit_sample)
    modal_instances(wiz_prefill_rows(list(wiz_samples()[[idx]]),
                                     keep_path = TRUE))
    modal_row_errors(list())
    modal_field_errors(list())
    showModal(add_sample_modal())
  }, ignoreInit = TRUE)

  add_modal_row <- function(type) {
    rid <- modal_row_id() + 1
    modal_row_id(rid)
    modal_instances(c(modal_instances(), list(list(id = rid, type = type))))
  }
  observeEvent(input$wiz_add_bw,         add_modal_row("bw"))
  observeEvent(input$wiz_add_hic,        add_modal_row("hic"))
  observeEvent(input$wiz_add_bed,        add_modal_row("bed"))
  observeEvent(input$wiz_add_paired_end, add_modal_row("paired-end"))

  observeEvent(input$wiz_remove_row, {
    modal_instances(Filter(function(r) r$id != input$wiz_remove_row,
                           modal_instances()))
  }, ignoreInit = TRUE)

  output$wiz_instance_rows <- renderUI({
    rows <- modal_instances()
    dr <- data_root()
    errs <- modal_row_errors()
    # dependency only: re-render from current inputs after every failed Save
    modal_fail_count()
    if (length(rows) == 0)
      return(p(class = "wiz-empty",
               "No files yet. Use the buttons above to add one."))
    tagList(lapply(rows, function(r) {
      name_id <- paste0("wiz_inst_name_", r$id)
      path_id <- paste0("wiz_inst_path_", r$id)
      row_err <- errs[[as.character(r$id)]]
      fluidRow(class = if (is.null(row_err))
                         "wiz-row" else "wiz-row wiz-row-error",
        column(2, div(class = "wiz-kind", WIZ_KIND_LABELS[[r$type]])),
        # r$name/r$path: a row added via wiz_prefill_rows()
        # (edit-sample prefill, or a new sample carrying over the previous
        # one's instance layout) carries its prefill value here, since its
        # input id is freshly minted and isolate(input[[...]]) is NULL until
        # the user actually touches it. A row added via "+ <type> instance"
        # has no r$name/r$path (NULL), so this falls through to "" exactly
        # as before.
        column(4, textInput(name_id, NULL,
                            value = isolate(input[[name_id]]) %||%
                                    r$name %||% "",
                            placeholder = "track name, e.g. ATAC")),
        # "" prepended (see add_sample_modal()): otherwise
        # selected="" matches no real <option> and the browser silently
        # falls back to the first matching file.
        # The current value is also added to choices: a path typed via
        # create = TRUE that isn't in the scoped file list would otherwise be
        # dropped on re-render, since selected must be one of the choices.
        column(5, selectizeInput(
          path_id, NULL,
          choices = unique(c("", isolate(input[[path_id]]) %||% r$path %||% "",
                             wizard_scoped_choices(dr, r$type))),
          selected = isolate(input[[path_id]]) %||% r$path %||% "",
          options = list(create = TRUE,
                         placeholder = "file, relative to the data folder")
        )),
        column(1, actionButton(paste0("wiz_rm_btn_", r$id), "\u2715",
                               class = "btn-xs", title = "Remove this row",
                               onclick = sprintf(
                                 "Shiny.setInputValue('wiz_remove_row', %d, {priority:'event'})",
                                 r$id
                               ))),
        if (!is.null(row_err))
          column(12, tags$div(class = "wiz-error",
                              style = "font-size:13px; margin-top:4px;",
                              row_err))
      )
    }))
  })

  observeEvent(input$wiz_confirm_sample, {
    edit_idx <- modal_edit_index()
    # the sample being edited was removed meanwhile
    if (identical(edit_idx, NA_integer_)) {
      showNotification(
        sprintf("Sample '%s' was removed, so these changes weren't saved.",
                modal_edit_name()),
        type = "warning"
      )
      modal_edit_name(NULL)
      removeModal()
      return()
    }
    name <- trimws(input$wiz_sample_name %||% "")
    existing_names <- vapply(wiz_samples(), function(s) s$condition,
                             character(1))
    # Editing a sample must not collide with its OWN pre-existing name --
    # excluding it here (rather than just allowing name == old name as a
    # special case) also lets the uniqueness check below stay a plain %in%.
    if (!is.null(edit_idx)) existing_names <- existing_names[-edit_idx]
    errs <- character(0)
    # row id (character) -> message, for wiz_instance_rows highlighting
    row_errs   <- list()
    # "quantdir" -> message, for add_sample_modal() highlighting
    field_errs <- list()

    if (nchar(name) == 0) errs <- c(errs, "Condition is required.")
    else if (name %in% existing_names)
      errs <- c(errs, sprintf("A condition named '%s' already exists.", name))

    data <- list()
    for (r in modal_instances()) {
      nm  <- trimws(input[[paste0("wiz_inst_name_", r$id)]] %||% "")
      pth <- trimws(input[[paste0("wiz_inst_path_", r$id)]] %||% "")
      # blank instance = not declared, silently dropped
      if (nchar(pth) == 0) next
      if (nchar(nm) == 0) {
        msg <- sprintf("A %s instance has a path but no name.", r$type)
        errs <- c(errs, msg); row_errs[[as.character(r$id)]] <- msg
        next
      }
      if (wizard_name_has_bad_dot(nm)) {
        msg <- sprintf("Instance name '%s' contains a bare '\u00b7' (surround it with spaces, or remove it). This breaks track-label formatting.",
                       nm)
        errs <- c(errs, msg); row_errs[[as.character(r$id)]] <- msg
        next
      }
      if (any(vapply(data, function(d) d$type == r$type && d$name == nm,
                     logical(1)))) {
        msg <- sprintf("Duplicate %s instance name '%s'.", r$type, nm)
        errs <- c(errs, msg); row_errs[[as.character(r$id)]] <- msg
        next
      }
      # is_real_file(), not file.exists(): file.exists() is TRUE for
      # directories too, which would otherwise be accepted here only to fail
      # (or, for paired-end, crash) at first render.
      if (!is.null(data_root()) && !is_real_file(file.path(data_root(), pth))) {
        msg <- sprintf(
          "Path for %s:%s does not resolve to a file under the data root: %s",
          r$type, nm, pth
        )
        errs <- c(errs, msg); row_errs[[as.character(r$id)]] <- msg
      } else if (identical(r$type, "paired-end") && !is.null(data_root()) &&
                 !is.null(problem <- wiz_paired_end_problem(
                   file.path(data_root(), pth)))) {
        msg <- sprintf("Paired-end file for %s (%s) %s.", nm, pth, problem)
        errs <- c(errs, msg); row_errs[[as.character(r$id)]] <- msg
      }
      data[[length(data) + 1]] <- list(type = r$type, name = nm, path = pth)
    }

    quantdir_path <- trimws(input$wiz_quantdir_path %||% "")
    if (nchar(quantdir_path) > 0) {
      quantdir_name <- trimws(input$wiz_quantdir_name %||% "")
      if (nchar(quantdir_name) == 0) quantdir_name <- "Quant"
      if (!is.null(data_root()) &&
          !is_real_file(file.path(data_root(), quantdir_path, "quant.sf"))) {
        msg <- sprintf("No quant.sf found under: %s", quantdir_path)
        errs <- c(errs, msg); field_errs$quantdir <- msg
      }
      data[[length(data) + 1]] <- list(type = "quantdir", name = quantdir_name,
                                       path = quantdir_path)
    }

    if (length(errs) > 0) {
      modal_row_errors(row_errs)
      modal_field_errors(field_errs)
      modal_fail_count(modal_fail_count() + 1)
      showModal(add_sample_modal(error = errs))
      return()
    }

    colour <- input$wiz_sample_colour %||% WIZARD_DEFAULT_COLOUR
    new_sample <- list(condition = name, colour = colour, data = data)
    if (is.null(edit_idx)) {
      wiz_samples(c(wiz_samples(), list(new_sample)))
    } else {
      samples <- wiz_samples()
      samples[[edit_idx]] <- new_sample
      wiz_samples(samples)
    }
    modal_edit_name(NULL)
    removeModal()
  })

  # ── Committed-sample list + edit/remove ─────────────────────────────────────
  output$wiz_samples_ui <- renderUI({
    samples <- wiz_samples()
    if (length(samples) == 0)
      return(p(class = "wiz-empty",
               "No samples yet. Add one per condition; you need at least one to save."))
    div(class = "wiz-samples", lapply(seq_along(samples), function(i) {
      s <- samples[[i]]
      inst_labels <- if (length(s$data) > 0)
        vapply(s$data, function(d) paste0(d$type, ":", d$name),
               character(1)) else character(0)
      fluidRow(class = "wiz-row",
        column(1, tags$span(style = sprintf("display:inline-block;width:14px;height:14px;background:%s;border-radius:2px;",
                                            s$colour %||% "#ccc"))),
        column(3, tags$b(s$condition)),
        column(5, paste(inst_labels, collapse = ", ")),
        column(3,
               actionButton(paste0("wiz_edit_sample_btn_", i), "Edit",
                            class = "pt-once btn-xs btn-default",
                            onclick = sprintf(
                              "Shiny.setInputValue('wiz_edit_sample', %s, {priority:'event'})",
                              jsonlite::toJSON(s$condition, auto_unbox = TRUE)
                            )),
               actionButton(paste0("wiz_rm_sample_btn_", i), "Remove",
                            class = "pt-once btn-xs btn-default",
                            onclick = sprintf(
                              "Shiny.setInputValue('wiz_remove_sample', %s, {priority:'event'})",
                              jsonlite::toJSON(s$condition, auto_unbox = TRUE)
                            )))
      )
    }))
  })

  # By condition name (unique in the wizard), not row number: a second quick
  # Remove sends a row number from the table as it was before the first removal,
  # which then names another sample.
  observeEvent(input$wiz_remove_sample, {
    wiz_samples(Filter(function(s) !identical(s$condition,
                                              input$wiz_remove_sample),
                       wiz_samples()))
  }, ignoreInit = TRUE)

  list(add_sample_modal = add_sample_modal, wiz_prefill_rows = wiz_prefill_rows,
       add_modal_row = add_modal_row)
}
