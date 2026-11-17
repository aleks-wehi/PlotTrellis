# The setup wizard's Close without saving (edit mode only), and the dirty guard
# that confirms it when a field differs from the project as opened: one block of
# wizard_server(). Returns list(is_dirty =).

# Closes the wizard unchanged: says so, then stops the app with "restart" so
# run_plottrellis()'s loop opens the composer on the old settings.
.pt_wizard_do_close <- function() {
  showModal(modalDialog(
    title = "Closing", easyClose = FALSE, footer = NULL,
    p("Returning to the app without saving."),
    p(strong("You can close this tab now."),
      " The app reopens in a new tab from RStudio in a few seconds. ",
      "If nothing opens, go back to RStudio: the app's address is in the Console.")
  ))
  message("PlotTrellis: settings unchanged, returning to the app.")
  later::later(function() stopApp(returnValue = "restart"), 0.2)
}

# The dirty guard: is_dirty() compares the form with the project as opened;
# Close asks before discarding changes.
.pt_wizard_dirty_server <- function(input, cfg, is_edit, data_root, wiz_samples,
                                    wiz_collect_flat_settings) {
  # Snapshot of every field this form can change, taken once at edit-mode
  # entry -- the baseline is_dirty() compares live values against, for the
  # Close-without-saving confirm guard. Not used
  # (and not needed) in fresh-setup mode: there's no prior configuration to
  # accidentally discard.
  original_snapshot <- if (is_edit) {
    list(
      data_root = cfg$data_dir,
      samples   = cfg$samples,
      flat      = cfg[c("font", "fontsize", "page_w", "chr_fmt",
                         "show_track_axes", "hic_norm", "hic_style",
                         "hic_colour_scale", "zrange_cap", "atac_height",
                         "hic_height",
                         "hic_res_mode", "hic_max_bins_per_view",
                         "hic_resolution", "diff_colour_lo", "diff_colour_hi",
                         "junction_opacity_cap", "min_junction_support",
                         "reference_colour", "comparison_colour",
                         "gene_flank_pct", "min_display_tpm",
                         "uniprot_reviewed_only", "arc_clip_partial",
                         "show_title", "show_scale_bar", "junction_height",
                         "diff_shared_alpha", "gtf_path", "orgdb",
                         "genome_label")]
    )
  } else NULL

  is_dirty <- reactive({
    if (!is_edit) return(FALSE)
    if (!identical(data_root(), original_snapshot$data_root)) return(TRUE)
    if (!identical(wiz_samples(), original_snapshot$samples)) return(TRUE)
    flat <- wiz_collect_flat_settings()
    any(vapply(names(flat),
               function(k) settings_value_differs(flat[[k]],
                                                  original_snapshot$flat[[k]]),
               logical(1)))
  })

  observeEvent(input$wiz_close, {
    if (is_dirty()) {
      showModal(modalDialog(
        title = "Discard changes?", easyClose = FALSE,
        p("You have unsaved changes. Closing now will discard them and return to the app with its previous settings."),
        footer = tagList(modalButton("Keep editing"),
                         actionButton("wiz_confirm_close", "Discard and close",
                                      class = "pt-once btn-danger"))
      ))
    } else {
      .pt_wizard_do_close()
    }
  })

  observeEvent(input$wiz_confirm_close, {
    removeModal()
    .pt_wizard_do_close()
  }, ignoreInit = TRUE)

  list(is_dirty = is_dirty)
}
