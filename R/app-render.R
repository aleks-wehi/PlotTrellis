# The composer's rendering (R/app-main.R): Generate, Cancel, the zoom buttons and
# Apply to figure all build entries from the picker and the live settings and
# draw through the session's render worker (NULL draws in-process), writing the
# result into `state`.
.pt_render_server <- function(input, output, session, state, ctx,
                              render_worker, live_cfg, resolved_region,
                              has_mode) {
  # Render helper (Generate, zoom and Apply to figure): .pt_render_svg() draws
  # the figure in the background render process and the page updates when it is
  # done; the displayed/downloadable R script is assemble_download_script() on
  # the same entries/region_params/settings.
  do_render <- function(region_params, entries, settings_base, title_str = "",
                        labels = NULL) {
    settings <- c(settings_base, list(
      title       = title_str,
      output_path = tempfile(fileext = ".svg")
    ))
    # The previous figure's error, if any, stays until this one is drawn:
    # clearing it now would end the error output in a silent req() failure,
    # which the page would take for an error and hide the loading overlay (and
    # its Cancel button) mid-render.
    state$rendering <- TRUE
    # Called from the app's event loop when the background render ends, outside
    # any reactive context: isolate() lets it read state. `res` is
    # .pt_render_svg()'s result: a list (the Hi-C colour tops the script
    # names) or the error message.
    done <- function(res) isolate({
      state$rendering <- FALSE
      err <- if (is.character(res)) res
      if (is.null(err)) err <- tryCatch({
        rcode <- assemble_download_script(entries, region_params, settings,
                                          hic_tops = res$hic_tops)
        NULL
      }, error = function(e) conditionMessage(e))
      if (is.null(err)) {
        if (!is.null(state$svg_path) && file.exists(state$svg_path))
          unlink(state$svg_path)
        state$plot_error       <- NULL
        state$svg_path         <- settings$output_path
        state$rcode             <- rcode
        state$last_region      <- region_params
        state$last_labels      <- labels
        state$last_title       <- title_str
        # Gene / protein info follows the figure
        state$rendered_gene    <- region_params$gene
      } else {
        unlink(settings$output_path)
        state$plot_error <- err
      }
      session$sendCustomMessage("show_output", list())
      session$sendCustomMessage("hide_spinner", list())
    })
    if (is.null(render_worker))
      done(.pt_render_svg(entries, region_params, settings))
    else render_worker$render(entries, region_params, settings, done)
  }

  # Cancel (the loading overlay's button): stop the render, keep the previous
  # figure.
  observeEvent(input$cancel_render, {
    if (isTRUE(state$rendering) && !is.null(render_worker)) {
      render_worker$cancel()
      state$rendering <- FALSE
      showNotification(
        "Rendering cancelled. The figure shown, if any, is the previous one.",
        type = "warning"
      )
    }
    session$sendCustomMessage("hide_spinner", list())
  })

  # Shared by Generate, zoom and Apply to figure: build_entries() +
  # settings_base + do_render(), from a region/labels/gene-name/title already
  # resolved by the caller.
  render_from_selection <- function(rgn, labels, gene_name, title_str) {
    cfg     <- live_cfg()
    entries <- .pt_build_entries(labels, ctx$TRACK_ITEMS,
                                 ctx$DEFAULT_SAMPLE_COLS, cfg)

    region_params <- list(
      chrom = rgn$chrom,
      start = as.integer(rgn$start),
      end   = as.integer(rgn$end),
      gene  = gene_name
    )
    settings_base <- .pt_settings_base(cfg)

    do_render(region_params, entries, settings_base, title_str, labels)
  }

  # ── Generate: build specs from UI then render ───────────────────────────────
  observeEvent(input$generate, {
    req(!isTRUE(state$rendering))
    req(has_mode())
    labels <- input$figure_tracks %||% character(0)
    req(length(labels) > 0)

    rr <- resolved_region()
    req(!is.null(rr))
    rgn       <- rr$rgn
    gene_name <- rr$gene_name

    if (identical(input$view_mode, "two")) {
      req(nchar(input$sel_reference %||% "") > 0,
          nchar(input$sel_comparison %||% "") > 0)
    }

    render_from_selection(rgn, labels, gene_name,
                          trimws(input$plot_title %||% ""))
  })

  # ── Zoom: recentre and rescale the last rendered region ─────────────────────
  # Rebuilds entries and settings from live_cfg(), keeping the
  # last figure's tracks, gene and title.
  zoom_to <- function(factor) {
    req(!isTRUE(state$rendering))
    req(!is.null(state$last_region), !is.null(state$last_labels))
    r         <- state$last_region
    mid       <- (r$start + r$end) / 2
    hw        <- max((r$end - r$start) / 2 / factor, 500L)
    clen      <- .pt_chr_length(ctx$CHR_LENGTHS, r$chrom, 2e8)
    new_start <- as.integer(max(1L,   round(mid - hw)))
    new_end   <- as.integer(min(clen, round(mid + hw)))
    new_rgn   <- list(chrom = r$chrom, start = new_start, end = new_end)
    updateSliderInput(session, "region_range",
                      value = c(new_start, new_end))
    updateTextInput(session, "region_start_text",
                    value = .pt_format_coord(new_start))
    updateTextInput(session, "region_end_text",
                    value = .pt_format_coord(new_end))
    render_from_selection(new_rgn, state$last_labels, r$gene, state$last_title)
  }
  observeEvent(input$zoom_in_15x,  zoom_to(1.5))
  observeEvent(input$zoom_in_3x,   zoom_to(3))
  observeEvent(input$zoom_in_10x,  zoom_to(10))
  observeEvent(input$zoom_out_15x, zoom_to(1 / 1.5))
  observeEvent(input$zoom_out_3x,  zoom_to(1 / 3))
  observeEvent(input$zoom_out_10x, zoom_to(0.1))

  # Apply to figure: re-render the last figure, same region, with the live
  # settings. A generate-btn, so the spinner/double-click gate covers it; shown
  # only once a figure exists.
  output$apply_settings_ui <- renderUI({
    req(!is.null(state$last_region))
    actionButton("apply_settings", "Apply to figure",
                 class = "btn-primary generate-btn")
  })
  observeEvent(input$apply_settings, zoom_to(1))

  invisible(NULL)
}
