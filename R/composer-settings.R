# The composer's Figure settings block (R/app.R): the saved and live settings,
# the panel, the unsaved marker, Save as project default and Revert.
# `render_cfg` is the project's render settings at launch, `settings_json` its
# project_settings.json path. Returns the reactives the other blocks read.
.pt_figure_settings_server <- function(input, output, render_cfg,
                                       settings_json) {
  # ── Figure settings ─────────────────────────────────────────────────────────
  # saved_cfg: the project's saved settings (render_cfg, then whatever Save as
  # project default writes). live_cfg: saved_cfg with the composer panel's
  # s_<key> inputs on top. Every render reads live_cfg(); core keys (genome,
  # data root, samples) and the derived paths always come from saved_cfg.
  # Assigned with x[names(v)] <- v, not modifyList(), which would drop the
  # NULL ("auto") values of zrange_cap/hic_height.
  saved_cfg <- reactiveVal(render_cfg)
  # bumped to re-render the panel from saved_cfg()
  settings_ui_gen <- reactiveVal(0)
  display_keys <- setdiff(names(.pt_default_settings()), .PT_CORE_SETTINGS)
  live_cfg <- reactive({
    x <- saved_cfg()
    v <- .pt_collect_display_settings(input, "s_", x)
    x[names(v)] <- v
    x
  })
  settings_dirty <- reactive({
    l <- live_cfg(); sv <- saved_cfg()
    any(vapply(display_keys,
               function(k) settings_value_differs(l[[k]], sv[[k]]),
               logical(1)))
  })

  output$display_settings_ui <- renderUI({
    settings_ui_gen()
    ui <- .pt_display_settings_ui(isolate(saved_cfg()), "s_")
    tagList(ui$display, .pt_settings_group("Comparison (diff) tracks", ui$diff))
  })
  # Rendered while the panel is still collapsed, so every s_ input reports
  # its value from the start.
  outputOptions(output, "display_settings_ui", suspendWhenHidden = FALSE)

  output$settings_dirty_ui <- renderUI({
    if (settings_dirty())
      tags$span(class = "text-warning", "Unsaved changes (this session only).")
  })

  # Save as project default: writes only the non-core keys into the saved
  # project_settings.json (re-read, so core keys and samples stay exactly as
  # saved), and refreshes project_samples.txt the way the wizard's Save does.
  observeEvent(input$save_settings, {
    live <- live_cfg()
    ok <- tryCatch({
      # Never write a project from defaults (e.g. reset from another tab
      # meanwhile).
      if (!file.exists(settings_json))
        stop(settings_json, " no longer exists; restart the app.",
             call. = FALSE)
      x <- .pt_read_settings(settings_json, .pt_default_settings())
      x[display_keys] <- live[display_keys]
      .pt_write_project_files(x, settings_json)
      TRUE
    }, error = function(e) {
      showNotification(paste("Couldn't save the settings:",
                             conditionMessage(e)),
                       type = "error", duration = NULL)
      FALSE
    })
    if (ok) {
      saved_cfg(live)
      showNotification("Saved as the project default.", type = "message")
    }
  })

  # Revert to saved: re-render the panel from saved_cfg(); the browser then
  # reports the re-rendered inputs' values.
  observeEvent(input$revert_settings, settings_ui_gen(settings_ui_gen() + 1))

  list(saved_cfg = saved_cfg, live_cfg = live_cfg,
       settings_dirty = settings_dirty,
       settings_ui_gen = settings_ui_gen, display_keys = display_keys)
}
