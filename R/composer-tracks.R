# The composer's track pool (R/app.R): which track items a mode and its chosen
# conditions offer, how they are labelled and ordered, and the preloaded start
# view's picks. `items` is the project's TRACK_ITEMS (build_track_items() plus
# build_comparison_track_items()), `start_view` the settings' start_view list or
# NULL.

.pt_start_condition <- function(i, start_view) {
  v <- start_view$conditions
  if (length(v) >= i) v[[i]] else ""
}

.pt_track_labels <- function(keys, items) {
  setNames(lapply(keys,
                  function(k) HTML(.pt_format_track_label(items[[k]]$label))),
           keys)
}

# Alphabetical (case-insensitive) by items' label, split on " · ", so tracks
# group by assay name, then by sample; never applied to figure_keys, whose order
# is deliberately curated.
.pt_sort_keys_by_label <- function(keys, items) {
  label_parts <- function(label) {
    label <- tolower(label)
    idx <- regexpr(" \u00B7 ", label, fixed = TRUE)
    if (idx == -1) return(c(label, ""))
    c(substr(label, 1, idx - 1),
      substr(label, idx + attr(idx, "match.length"), nchar(label)))
  }
  parts <- vapply(keys, function(k) label_parts(items[[k]]$label), character(2))
  keys[order(parts[1, ], parts[2, ])]
}

# Every sample-count mode starts with just the gene-annotation track (or the
# preloaded start view's tracks), a plain, predictable starting point; anything
# else goes in from Available tracks explicitly.
.pt_default_figure_keys <- function(conditions, start_view)
  if (!is.null(start_view) && identical(conditions, start_view$conditions))
    start_view$tracks else "genes"

# Item pool + default figure/available split for 1-2 conditions, filtered by
# each item's own `sample` field; fixed_item()s are found by scanning items
# itself, never a hardcoded type list, so a new sample-free track type shows up
# here without an edit.
.pt_track_split <- function(conditions, items, start_view) {
  cond_items <- Filter(function(sp) isTRUE(sp$sample %in% conditions), items)
  comparison_keys <- character(0)
  if (length(conditions) > 1) {
    comparison_keys <- names(items)[vapply(items, function(sp) {
      !is.null(sp$samples) && length(sp$samples) == 2 &&
        identical(sp$samples, conditions)
    }, logical(1))]
  }
  # is.null(sp$samples) is required alongside is.na(sp$sample) here -- every
  # comparison item (detect_pair_generic()) also sets sample = NA (it uses
  # `samples`, plural, instead), so without this every comparison track of
  # EVERY condition pair (not just the one selected above) would land in
  # "fixed" and appear whatever the conditions.
  fixed_keys  <- names(Filter(
    function(sp) isTRUE(is.na(sp$sample)) && is.null(sp$samples),
    items
  ))
  pool_keys   <- c(names(cond_items), comparison_keys, fixed_keys)
  figure_keys <- .pt_default_figure_keys(conditions, start_view)
  figure_keys <- figure_keys[figure_keys %in% pool_keys]
  # Gene annotation stays in Available whatever Figure tracks holds: a figure
  # can hold it many times (the composer JS keeps exactly one in Available).
  avail_keys  <- .pt_sort_keys_by_label(union(setdiff(pool_keys, figure_keys),
                                              intersect("genes", pool_keys)),
                                        items)
  list(available = .pt_track_labels(avail_keys, items),
       figure = .pt_track_labels(figure_keys, items))
}

# Multi's "Available tracks" pool: items minus EVERY two-sample comparison item:
# a comparison is a Reference-vs-Comparison concept that Two-sample mode offers
# for its own pair, and offering every pair's again here would only be clutter.
# Gene annotation is in it too, as in every mode (Figure tracks also starts with
# one).
# Sorted alphabetically by track label (.pt_sort_keys_by_label()) so the "All"
# filter shows a stable, name-grouped order (ATAC vs RNA-seq vs Junctions, etc.)
# rather than items' own sample-then-type insertion order.
.pt_compose_available_keys <- function(items)
  .pt_sort_keys_by_label(Filter(function(k) {
    sp <- items[[k]]
    is.null(sp$samples) || length(sp$samples) != 2
  }, names(items)), items)

# Choices for the Condition/Reference/Comparison selects: every
# available_conditions entry is always offered plus a blank placeholder;
# collisions (picking the value already sat in the other box) are resolved by
# clearing the other box server-side, not by removing options here.
.pt_condition_choices <- function(available_conditions) {
  c(setNames("", "(select)"),
    setNames(available_conditions, available_conditions))
}
