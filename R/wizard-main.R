# The setup wizard: a Shiny app of its own (wizard_ui()/wizard_server()) for a
# first-time setup and for editing an existing project's settings, with nothing
# from the composer loaded. Save writes project_settings.json; Save and Close
# then stop it with "restart", and run_plottrellis()'s loop starts the composer.
# The samples live in project_settings.json; a samples.txt is only imported or
# exported here.

# Builds the wizard's server function for one project: `SETTINGS_JSON` is the
# settings file it writes on Save, `CFG` the project's current settings (shown
# for editing when that file exists), `DEFAULT_SETTINGS` the defaults a
# first-time setup starts from.
wizard_server <- function(SETTINGS_JSON, CFG,
                          DEFAULT_SETTINGS) function(input, output, session) {
  # repeat-click gate for the "pt-once" buttons (R/app-once-gate.R)
  .pt_once_server(input, session)

  is_edit <- file.exists(SETTINGS_JSON)
  cfg     <- if (is_edit) CFG else DEFAULT_SETTINGS

  # committed samples:
  # list(condition=,colour=,data=list(list(type=,name=,path=)))
  wiz_samples     <- reactiveVal(if (is_edit) cfg$samples else list())
  # bw/hic/bed/paired-end rows in the currently-open add-sample modal
  modal_instances <- reactiveVal(list())
  modal_row_id    <- reactiveVal(0)
  # NULL = the modal is adding a brand-new sample; a string = the condition
  # name of the sample currently being edited, so wiz_confirm_sample can
  # replace it in place instead of appending a duplicate. A name, not a row
  # number: a Remove handled while the modal opens would shift the rows.
  # Set by wiz_add_sample/wiz_edit_sample (R/wizard-samples.R), cleared on
  # confirm -- a re-click of either always resets it, so a stale value can't
  # survive a Cancel.
  modal_edit_name <- reactiveVal(NULL)
  # Row of the sample being edited in wiz_samples(), NULL when adding, NA if it
  # was removed meanwhile.
  modal_edit_index <- function() {
    nm <- modal_edit_name()
    if (is.null(nm)) NULL else match(nm, vapply(wiz_samples(),
                                                function(s) s$condition,
                                                character(1)))
  }
  # Per-row/per-field error messages from the last failed wiz_confirm_sample
  # attempt, so a validation error can highlight the specific offending
  # row/file instead of only listing errors in one combined banner. Keyed by
  # row id (character) for modal_row_errors, by "quantdir" for
  # modal_field_errors. Both reset to list() on a genuinely fresh modal open.
  modal_row_errors   <- reactiveVal(list())
  modal_field_errors <- reactiveVal(list())
  # Bumped on every failed Save: re-showing the modal replaces
  # uiOutput("wiz_instance_rows"), and shiny.js repaints it from its cached
  # last render, so without a fresh render the rows snap back to stale values.
  modal_fail_count   <- reactiveVal(0)
  # NULL = landing screen (choose a path); "import" / "manual" = which entry
  # path was chosen; "edit" = re-entered from the running app to edit an
  # existing project, which skips the landing screen entirely --
  # there's no "start fresh vs. import" choice to make when a project already
  # exists. Every non-NULL mode renders the same wiz_shared_form() below an
  # optional import box -- see wiz_shared_form()'s comment.
  # A first-time wizard that finds the project_samples.txt an
  # earlier Save left beside project_settings.json (settings since reset or
  # deleted) opens straight on Import, prefilled, so the project rebuilds in
  # one click. Its own name, so it never overwrites a user's samples.txt.
  mirror_txt <- file.path(dirname(SETTINGS_JSON), "project_samples.txt")
  has_mirror <- !is_edit && file.exists(mirror_txt)
  wiz_mode <- reactiveVal(
    if (is_edit) "edit" else if (has_mirror) "import" else NULL
  )
  # Display settings imported from a samples.txt (the wizard has no inputs for
  # them): Save writes cfg's values with these on top.
  import_visual <- reactiveVal(list())
  # the lines of the samples.txt last imported in this wizard session
  imported      <- reactiveVal(NULL)
  # After a Save that failed, the "Not saved yet" list re-checks as the fields
  # change, so it never shows a problem already fixed (e.g. by an Import), and
  # goes once all pass.
  save_failed <- reactiveVal(FALSE)
  # The last Save's write failure (a folder that can't be written), cleared by
  # the next Save.
  write_error <- reactiveVal(NULL)

  # Each block lives in R/wizard-<block>.R. Keep this order: it is the order the
  # observers were created in before the split. What a block returns keeps its
  # name here, where the tests read it.
  .pt_wizard_mode_server(input, output, wiz_mode, save_failed, write_error,
                         cfg, is_edit, has_mirror, mirror_txt)
  data_root <- .pt_wizard_data_folder_server(input, output)$data_root
  gtf_check <- .pt_wizard_genome_server(input, output, session)$gtf_check
  samples_block <- .pt_wizard_samples_server(
    input, output, wiz_samples, data_root, modal_instances, modal_row_id,
    modal_edit_name, modal_edit_index, modal_row_errors, modal_field_errors,
    modal_fail_count
  )
  add_sample_modal <- samples_block$add_sample_modal
  wiz_prefill_rows <- samples_block$wiz_prefill_rows
  add_modal_row    <- samples_block$add_modal_row
  .pt_wizard_import_server(input, output, session, wiz_samples, imported,
                           import_visual, DEFAULT_SETTINGS)
  export_block <- .pt_wizard_export_server(
    input, output, data_root, wiz_samples, import_visual, cfg, DEFAULT_SETTINGS
  )
  wiz_collect_flat_settings <- export_block$wiz_collect_flat_settings
  wiz_display_settings      <- export_block$wiz_display_settings
  is_dirty <- .pt_wizard_dirty_server(input, cfg, is_edit, data_root,
                                      wiz_samples,
                                      wiz_collect_flat_settings)$is_dirty
  save_problems <- .pt_wizard_save_server(
    input, output, SETTINGS_JSON, mirror_txt, data_root, wiz_samples, imported,
    wiz_collect_flat_settings, save_failed, write_error
  )$save_problems
}
