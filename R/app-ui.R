# The JS maps the page script reads: track key -> type, key -> broken, type ->
# filter button, and the gene-required types, as JSON strings.
.pt_track_js_maps <- function(items, gene_required_types) {
  # JS map from track key → type, embedded in page for client-side filter
  TRACK_TYPE_JS <- jsonlite::toJSON(
    setNames(lapply(items, function(x) x$type), names(items)),
    auto_unbox = TRUE
  )
  # JS map from track key -> broken: a declared file that is missing (the
  # `broken` field viz-track-detection.R sets on per-sample items, checked once at
  # startup) gets red-boxed in the track picker, the same look as the
  # gene-required flag. isTRUE() turns the common case (an item with no `broken`
  # field at all, e.g. genes, transcripts, comparison tracks) from NULL into
  # FALSE, rather than needing every module to declare it.
  TRACK_BROKEN_JS <- jsonlite::toJSON(
    setNames(lapply(items, function(x) isTRUE(x$broken)), names(items)),
    auto_unbox = TRUE
  )
  # JS map from track type -> the id of the filter button that lists it, from
  # each module's info_<key>() `filter` (.pt_type_filters()): e.g. the
  # comparison types share "comparison" ("Comparisons", Two-sample mode only),
  # and the gene annotation and transcript-structure tracks share "annotation"
  # ("Annotations").
  TRACK_FILTER_JS <- jsonlite::toJSON(.pt_type_filters(), auto_unbox = TRUE)
  GENE_REQUIRED_TYPES_JS <- jsonlite::toJSON(gene_required_types)
  list(TRACK_TYPE_JS = TRACK_TYPE_JS, TRACK_BROKEN_JS = TRACK_BROKEN_JS,
       TRACK_FILTER_JS = TRACK_FILTER_JS,
       GENE_REQUIRED_TYPES_JS = GENE_REQUIRED_TYPES_JS)
}

# The composer page's head: the shared theme, the once gate and help assets, the
# composer's own stylesheet and page script (inst/www/composer.css and
# composer.js), and the JS maps the script reads, which depend on the project's
# tracks.
.pt_composer_head <- function(ctx) {
  maps <- .pt_track_js_maps(ctx$TRACK_ITEMS, ctx$GENE_REQUIRED_TYPES)
  tagList(
    tags$head(.pt_theme_css(), .pt_once_js(), .pt_help_assets(),
              includeCSS(system.file("www", "composer.css",
                                     package = "PlotTrellis",
                                     mustWork = TRUE))),
    tags$head(tags$script(HTML(paste0("
    var TRACK_TYPE_MAP = ", maps$TRACK_TYPE_JS, ";
    var TRACK_BROKEN_MAP = ", maps$TRACK_BROKEN_JS, ";
    var TRACK_FILTER_MAP = ", maps$TRACK_FILTER_JS, ";
    var GENE_REQUIRED_TYPES = ", maps$GENE_REQUIRED_TYPES_JS, ";"))),
              includeScript(system.file("www", "composer.js",
                                        package = "PlotTrellis",
                                        mustWork = TRUE)))
  )
}

# The composer page (R/app-main.R's `ui`): a navbarPage with the head assets and one
# tab holding the project actions, the Compose section and the Output section.
# Built once per launch from the startup context.
.pt_composer_ui <- function(ctx) {
  navbarPage(
    # shared with the setup wizard's title bar (.pt_title_bar())
    title = .pt_brand(),
    windowTitle = "PlotTrellis",
    id    = "main_tabs",

    # In `header`, not loose children: navbarPage treats every other child as a
    # nav panel and bslib warns about each one. The tags$head() contents are
    # hoisted to <head> either way.
    header = .pt_composer_head(ctx),

    tabPanel("",

      # Edit settings / close buttons, fixed to the top-right navbar area via
      # CSS. Clear all caches and Reset project configuration are in the Edit
      # project settings dialog's collapsed Advanced section.
      .pt_project_actions_ui(),
      .pt_compose_panel_ui(ctx),
      .pt_output_panel_ui()
    ) # /tabPanel
  )
}

# The Compose section: mode radio, condition selects, region (gene or
# coordinates), the track picker and the Figure settings panel, then Generate.
.pt_compose_panel_ui <- function(ctx) {
  # ── Compose: one 3-way mode radio ───────────────────────────────────────────
  # "single"/"two" (One/Two) pick 1-2 conditions from dropdowns and narrow the
  # track pool to them; "compose" (Multi) offers every sample's tracks, with the
  # track picker as its only sample selector.
  .pt_collapsible_panel("compose", "Compose",

    h4("Number of samples"),
    radioButtons("view_mode", NULL,
                 choices  = c("One" = "single",
                              "Two" = "two",
                              "Multi" = "compose"),
                 selected = ctx$START_VIEW$view_mode %||% character(0),
                 inline = TRUE),
    hr(),

    # Condition selection -- single/two-sample modes only; track-compose has
    # no dropdown, the track picker below is the only sample selector there.
    # Placed above Genomic region so Reference/Comparison (or Condition) is
    # chosen before searching for a gene. All start unselected:
    # region_mode/gene_name and the track picker stay locked until these
    # resolve too, not just until a mode is chosen (see
    # has_conditions_ready() in server()).
    conditionalPanel("input.view_mode === 'single'",
      fluidRow(column(4,
        selectizeInput("sel_condition", "Condition:",
                    choices = .pt_condition_choices(ctx$AVAILABLE_CONDITIONS),
                    selected = if (identical(ctx$START_VIEW$view_mode,
                                             "single"))
                      .pt_start_condition(1, ctx$START_VIEW) else "")
      ))
    ),
    conditionalPanel("input.view_mode === 'two'",
      fluidRow(
        column(4,
          selectizeInput("sel_reference", "Reference:",
                      choices = .pt_condition_choices(ctx$AVAILABLE_CONDITIONS),
                      selected = if (identical(ctx$START_VIEW$view_mode, "two"))
                        .pt_start_condition(1, ctx$START_VIEW) else "")
        ),
        column(4,
          selectizeInput("sel_comparison", "Comparison:",
                      choices = .pt_condition_choices(ctx$AVAILABLE_CONDITIONS),
                      selected = if (identical(ctx$START_VIEW$view_mode, "two"))
                        .pt_start_condition(2, ctx$START_VIEW) else "")
        )
      )
    ),

    # Region (full width)
    fluidRow(
      column(12,
        h4("Genomic region"),
        radioButtons("region_mode", NULL,
                     choices  = c("Gene name" = "gene",
                                  "Coordinates" = "coords"),
                     selected = if (is.null(ctx$START_VIEW$gene))
                       character(0) else "gene",
                     inline = TRUE),
        hr(style = "margin:8px 0;"),

        conditionalPanel("input.region_mode === 'gene'",
          fluidRow(
            column(4,
              selectizeInput("gene_name", "Gene symbol or Ensembl ID:",
                             choices = NULL,
                             # maxOptions caps how many server-side search
                             # matches selectize will ever render, not a "load
                             # more on scroll" limit -- a small cap silently
                             # truncates most typed queries well before the user
                             # could scroll to a match. create = TRUE lets an
                             # Ensembl ID (not in the symbol list) be entered;
                             # .pt_resolve_gene_region() maps it to its symbol.
                             options = list(placeholder = "e.g. Pax5",
                                            maxOptions = 1000, create = TRUE))
            ),
            column(6, br(), uiOutput("gene_result_ui"))
          )
        ),

        # Coordinates panel, shown whenever region_mode is 'coords' --
        # selectable in every sample-count mode; tracks that need a gene symbol
        # can't draw from coordinates and are red-boxed in the picker.
        conditionalPanel("input.region_mode === 'coords'",
          fluidRow(
            column(3, uiOutput("chr_select_ui")),
            column(9, uiOutput("range_slider_ui"))
          ),
          fluidRow(
            column(3, uiOutput("start_num_ui")),
            column(3, uiOutput("end_num_ui")),
            column(6, br(), uiOutput("region_error_ui"))
          )
        )
      )
    ),

    hr(),

    # Track assignment (full width)
    fluidRow(
      column(12,
        h4("Track assignment"),

        div(id = "tracks-block", `data-available-id` = "available_tracks",
          # what to pick first, while the lists are empty
          uiOutput("tracks_hint"),
          div(class = "track-filters",
            tags$button("All",
                        class = "btn btn-lg btn-default filter-btn active",
                        `data-filter` = "all"),
            # One button per filter group the modules declare (info_<key>()'s
            # `filter`). "Comparisons" shows in Two-sample mode only: one-sample
            # has none, and track-compose's pool (ctx$COMPOSE_AVAILABLE_KEYS)
            # excludes every comparison item outright, so there is nothing for
            # it to select there either. Only the types this project has tracks
            # of: an empty filter would show nothing.
            .pt_filter_buttons(Filter(
              function(info) info$key %in%
                vapply(ctx$TRACK_ITEMS, `[[`, "", "type"),
              .pt_viz_infos()
            ))
          ),

          bucket_list(
            header      = NULL,
            group_name  = "tracks",
            orientation = "horizontal",
            add_rank_list(
              text     = "Available tracks (double-click)",
              labels   = .pt_track_labels(ctx$COMPOSE_AVAILABLE_KEYS,
                                          ctx$TRACK_ITEMS),
              input_id = "available_tracks"
            ),
            add_rank_list(
              text     = "Figure tracks (drag to reorder)",
              labels   = list(genes = HTML(.pt_format_track_label(
                ctx$TRACK_ITEMS[["genes"]]$label
              ))),
              input_id = "figure_tracks"
            )
          )
        )
      )
    ),

    # Every non-core setting, live for this session. Changes apply
    # on Generate, zoom or Apply to figure (no automatic re-render).
    .pt_collapsible_panel("figure_settings", "Figure settings",
      collapsed = TRUE,
      .pt_plot_title_input(),
      p(class = "text-muted",
        "Changes apply on the next Generate, zoom or Apply to figure, for this session only. Save as project default keeps them for next time."),
      uiOutput("display_settings_ui"),
      tags$div(style = "margin-top:10px; display:flex; gap:8px; align-items:center; flex-wrap:wrap;",
        uiOutput("apply_settings_ui", inline = TRUE),
        actionButton("save_settings", "Save as project default",
                     class = "btn-default pt-once"),
        actionButton("revert_settings", "Revert to saved",
                     class = "btn-default pt-once"),
        uiOutput("settings_dirty_ui", inline = TRUE)
      )
    ),
    tags$div(style = "height:8px;"),
    fluidRow(
      column(4, uiOutput("generate_ui")),
      column(8, uiOutput("generate_warning_ui"))
    ),
    tags$div(style = "height:16px;")
  ) # /collapsible_panel: compose
}

# The Output section: the figure, its errors, zoom, downloads, Gene / protein
# info and the R script, hidden until the first render.
.pt_output_panel_ui <- function() {
  # ── Output section (hidden until first render) ──────────────────────────────
  # Full-width, no inner scrollbox -- the figure renders at its natural size
  # and the whole page scrolls to reveal more of it (a max-height/overflow-y
  # box makes tall figures feel clipped/hidden).
  .pt_collapsible_panel("output", "Output",
    uiOutput("render_error_ui"),
    uiOutput("zoom_ui"),
    uiOutput("svg_display_ui"),
    tags$div(style = "margin-top:10px;",
      downloadButton("dl_svg", "Download SVG", class = "btn-lg pt-once"),
      tags$span(style = "margin-left:6px;",
        `data-toggle` = "tooltip", `data-placement` = "top",
        title = "An R script that reproduces this figure.",
        downloadButton("dl_rcode", "Download R script",
                       class = "btn-lg pt-once")
      )
    ),
    .pt_gene_info_ui(.pt_collapsible_panel),
    tags$div(style = "margin-top:12px;",
      .pt_collapsible_panel("rcode", "R script", collapsed = TRUE,
        div(style = "max-height: 60vh; overflow-y: auto; font-size: var(--bio-font-code);",
          verbatimTextOutput("rcode_output"))
      )
    ),
    wrapper_id = "output-section", wrapper_style = "display:none;"
  )
}
