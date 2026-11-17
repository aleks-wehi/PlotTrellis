# The composer's track lists (R/track-lists.R): the markup, update message and input values the
# composer's script (R/app.R), plottrellis.css and the browser checks rely on.

composer_lists <- function() {
  bucket_list(
    header = NULL, group_name = "tracks", orientation = "horizontal",
    add_rank_list(text = "Available tracks (double-click)",
                  labels = list(bw_a = HTML("<b>A</b>"), genes = "Gene annotation", "plain"),
                  input_id = "available_tracks"),
    add_rank_list(text = "Figure tracks (drag to reorder)", labels = list(genes = HTML("G")),
                  input_id = "figure_tracks"))
}

test_that("bucket_list() writes the two lists, their items and the Sortable set-up", {
  r <- htmltools::renderTags(composer_lists())
  html <- gsub(">\\s+<", "><", r$html)
  expect_match(html, paste0(
    '<div class="bucket-list-container default-sortable" id="bucket-list-tracks">',
    '<div class="default-sortable bucket-list bucket-list-horizontal">',
    '<div class="rank-list-container default-sortable column_1" id="rank-list-available_tracks">',
    '<p class="rank-list-title">Available tracks (double-click)</p>',
    '<div class="rank-list" id="available_tracks">',
    '<div class="rank-list-item" data-rank-id="bw_a"><b>A</b></div>',
    '<div class="rank-list-item" data-rank-id="genes">Gene annotation</div>',
    '<div class="rank-list-item">plain</div></div></div>',
    '<div class="rank-list-container default-sortable column_2" id="rank-list-figure_tracks">',
    '<p class="rank-list-title">Figure tracks (drag to reorder)</p>',
    '<div class="rank-list" id="figure_tracks"><div class="rank-list-item" data-rank-id="genes">G</div>',
    '</div></div></div></div>'), fixed = TRUE)
  expect_match(html, 'ptTrackLists.create("tracks", ["available_tracks","figure_tracks"], ["available_tracks","figure_tracks"]);',
               fixed = TRUE)
  deps <- Filter(function(d) d$name %in% c("SortableJS", "pt-track-lists"), r$dependencies)
  expect_identical(vapply(deps, `[[`, "", "name"), c("SortableJS", "pt-track-lists"))  # Sortable first
  expect_identical(deps[[2]]$version, as.character(utils::packageVersion("PlotTrellis")))
  files <- file.path(deps[[1]]$src$file, c("Sortable.min.js", "track-lists.js", "track-lists.css", "LICENSE"))
  expect_true(all(file.exists(files)))
  expect_identical(c(deps[[1]]$script, deps[[2]]$script, deps[[2]]$stylesheet), basename(files[1:3]))
  expect_match(readLines(files[1], n = 1), "Sortable 1.15.6 - MIT", fixed = TRUE)
})

test_that("update_rank_list() sends the list's new items, an empty list, or a new title", {
  sent <- list()
  session <- list(sendInputMessage = function(id, message) sent[[length(sent) + 1]] <<- list(id = id, message = message))
  update_rank_list("figure_tracks", labels = list(genes = HTML("G"), "plain"), session = session)
  update_rank_list("available_tracks", labels = list(), session = session)
  update_rank_list("figure_tracks", text = "Figure", session = session)
  expect_identical(sent[[1]], list(id = "rank-list-figure_tracks", message = list(
    id = "rank-list-figure_tracks",
    labels = '<div class="rank-list-item" data-rank-id="genes">G</div>\n<div class="rank-list-item">plain</div>')))
  expect_identical(sent[[2]], list(id = "rank-list-available_tracks",
                                   message = list(id = "rank-list-available_tracks", labels = list())))
  expect_identical(sent[[3]]$message, list(id = "rank-list-figure_tracks", text = "Figure"))
})

test_that("the lists' input values arrive as character vectors (handlers registered on load)", {
  rank <- shiny:::inputHandlers$get("sortablejs.rank_list")
  bucket <- shiny:::inputHandlers$get("sortablejs.bucket_list")
  expect_identical(rank(list("genes", "bw_a"), NULL, "figure_tracks"), c("genes", "bw_a"))
  expect_identical(rank(list(), NULL, "figure_tracks"), character(0))
  expect_identical(bucket(list(available_tracks = list("bw_a"), figure_tracks = list()), NULL, "tracks"),
                   list(available_tracks = "bw_a", figure_tracks = character(0)))
})
