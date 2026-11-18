# The composer's drag-and-drop track lists (Available tracks, Figure tracks):
# two lists in one drag group, on SortableJS (inst/www/track-lists/, with its
# licence). The same names, markup, inputs and update message as the sortable
# package's bucket_list(), add_rank_list() and update_rank_list(), which they
# replace, so the composer's script and styles (R/app-main.R, plottrellis.css) work
# on them unchanged. Only what the composer uses is here.
#
# Inputs, per list: `input$<input_id>`, the items' ids (each item's name in
# `labels`, else its text) in list order, sent on load and after every drag or
# update; and `input$<group_name>`, every list's ids by input id.

# One list of a bucket_list(): its title, items (a named list of labels; names
# are the ids) and input id, which is also the list element's id.
add_rank_list <- function(text, labels = NULL, input_id) {
  list(text = text, labels = labels, input_id = input_id)
}

# The lists given (add_rank_list()), side by side, sharing one drag group.
bucket_list <- function(header = NULL, ..., group_name,
                        orientation = c("horizontal", "vertical")) {
  orientation <- match.arg(orientation)
  lists <- list(...)
  ids <- vapply(lists, `[[`, character(1), "input_id")
  boxes <- lapply(seq_along(lists), function(i) {
    l <- lists[[i]]
    tags$div(class = paste0("rank-list-container default-sortable column_", i),
             id = paste0("rank-list-", l$input_id),
             if (!is.null(l$text) && nchar(l$text) > 0)
               tags$p(class = "rank-list-title", l$text),
             tags$div(class = "rank-list", id = l$input_id,
                      .pt_rank_list_items(l$labels)))
  })
  tagList(
    tags$div(class = "bucket-list-container default-sortable",
             id = paste0("bucket-list-", group_name),
             if (!is.null(header)) tags$p(class = "bucket-list-header", header),
             tags$div(
               class = paste0("default-sortable bucket-list bucket-list-",
                              orientation),
               boxes
             )),
    tags$script(HTML(sprintf("ptTrackLists.create(%s, %s, %s);",
                             jsonlite::toJSON(group_name, auto_unbox = TRUE),
                             jsonlite::toJSON(ids), jsonlite::toJSON(ids)))),
    .pt_track_lists_dependency()
  )
}

# Replaces a list's items (labels, as for add_rank_list()) and/or its title; the
# list then reports its new ids as `input$<css_id>`. An empty `labels` empties
# the list.
update_rank_list <- function(css_id, text = NULL, labels = NULL,
                             session = shiny::getDefaultReactiveDomain()) {
  input_id <- paste0("rank-list-", css_id)
  if (!is.null(labels) && length(labels) > 0)
    labels <- as.character(tagList(.pt_rank_list_items(labels)))
  message <- list(id = input_id, text = text, labels = labels)
  session$sendInputMessage(input_id,
                           message[!vapply(message, is.null, logical(1))])
}

# One item per label; a label's name (if any) is its id.
.pt_rank_list_items <- function(labels) {
  ids <- names(labels) %||% rep("", length(labels))
  ids[is.na(ids)] <- ""
  mapply(function(label, id) tags$div(class = "rank-list-item",
                                      `data-rank-id` = if (nzchar(id)) id,
                                      label),
         labels, ids, SIMPLIFY = FALSE, USE.NAMES = FALSE)
}

# SortableJS, then the lists' script (Shiny input bindings and Sortable set-up)
# and styles, versioned with the package so a browser never keeps an older copy.
.pt_track_lists_dependency <- function() {
  list(htmltools::htmlDependency("SortableJS", "1.15.6",
                                 src = "www/track-lists",
                                 package = "PlotTrellis",
                                 script = "Sortable.min.js", all_files = FALSE),
       htmltools::htmlDependency(
         "pt-track-lists", as.character(utils::packageVersion("PlotTrellis")),
         src = "www/track-lists",
         package = "PlotTrellis",
         script = "track-lists.js",
         stylesheet = "track-lists.css",
         all_files = FALSE
       ))
}

# The lists' input values arrive as JSON arrays: as character vectors (a
# bucket's, one per list).
.pt_register_track_list_inputs <- function() {
  as_ids <- function(x) if (length(x) == 0) character(0) else unlist(x)
  shiny::registerInputHandler(
    "sortablejs.rank_list",
    function(val, shinysession, name) as_ids(val),
    force = TRUE
  )
  shiny::registerInputHandler(
    "sortablejs.bucket_list",
    function(val, shinysession, name) lapply(val, as_ids),
    force = TRUE
  )
}

.onLoad <- function(libname, pkgname) {
  .pt_register_track_list_inputs()
}
