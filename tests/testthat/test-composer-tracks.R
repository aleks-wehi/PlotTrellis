# The composer's track pool (R/composer-tracks.R) on hand-made track items
# (helper-composer-fakes.R): samples A and B, no genome, no data files.

SORTED_ALL <- c("bw_atac_A", "bw_atac_B", "bw_diff_stacked_atac_B_vs_A", "bw_diff_stacked_atac_A_vs_B",
                "genes", "hic_hic_A", "tpm_usage_Quant_A", "transcripts")

test_that(".pt_start_condition reads the start view's i-th condition, else blank", {
  sv <- list(conditions = c("A", "B"))
  expect_identical(.pt_start_condition(1, sv), "A")
  expect_identical(.pt_start_condition(2, sv), "B")
  expect_identical(.pt_start_condition(3, sv), "")
  expect_identical(.pt_start_condition(1, NULL), "")
  expect_identical(.pt_start_condition(1, list(conditions = character(0))), "")
  expect_identical(.pt_start_condition(2, list(conditions = list("A", "B"))), "B")
  expect_identical(.pt_start_condition(1, list(view_mode = "compose")), "")
})

test_that(".pt_track_labels names each label by its key and bolds the track name", {
  items <- fake_track_items()
  labs <- .pt_track_labels(c("bw_atac_A", "genes"), items)
  expect_identical(names(labs), c("bw_atac_A", "genes"))
  expect_s3_class(labs$bw_atac_A, "html")
  expect_identical(as.character(labs$bw_atac_A), "<strong>atac</strong> \u00B7 A")
  expect_identical(as.character(labs$genes), "<strong>Gene annotation</strong>")
  # A repeated key (Gene annotation can be in a figure several times) gives a repeated entry.
  expect_identical(names(.pt_track_labels(c("genes", "genes"), items)), c("genes", "genes"))
})

test_that(".pt_track_labels of no keys is an empty list, and an unknown key an empty label", {
  items <- fake_track_items()
  expect_length(.pt_track_labels(character(0), items), 0)
  # An unknown key is not an error: its label is empty.
  lab <- .pt_track_labels("no_such_key", items)
  expect_identical(names(lab), "no_such_key")
  expect_identical(as.character(lab[[1]]), "")
})

test_that(".pt_sort_keys_by_label sorts by track name, then by sample, ignoring case", {
  items <- fake_track_items()
  expect_identical(.pt_sort_keys_by_label(rev(names(items)), items), SORTED_ALL)
  expect_identical(.pt_sort_keys_by_label(names(items), items), SORTED_ALL)
  ci <- list(b = list(label = "beta \u00B7 x"), a = list(label = "Alpha \u00B7 y"),
             c = list(label = "alpha \u00B7 X"), d = list(label = "ALPHA"))
  # "alpha" with no sample first (its second part is blank), then samples x, y.
  expect_identical(.pt_sort_keys_by_label(names(ci), ci), c("d", "c", "a", "b"))
})

test_that(".pt_sort_keys_by_label keeps the given order for equal labels, and takes no keys", {
  ci <- list(p = list(label = "same \u00B7 s"), q = list(label = "Same \u00B7 S"))
  expect_identical(.pt_sort_keys_by_label(c("p", "q"), ci), c("p", "q"))
  expect_identical(.pt_sort_keys_by_label(c("q", "p"), ci), c("q", "p"))
  expect_identical(.pt_sort_keys_by_label(character(0), ci), character(0))
})

test_that(".pt_default_figure_keys is Gene annotation unless the start view's conditions match", {
  sv <- list(conditions = c("A", "B"), tracks = c("bw_atac_A", "genes"))
  expect_identical(.pt_default_figure_keys("A", NULL), "genes")
  expect_identical(.pt_default_figure_keys(c("A", "B"), sv), c("bw_atac_A", "genes"))
  expect_identical(.pt_default_figure_keys(c("B", "A"), sv), "genes")
  expect_identical(.pt_default_figure_keys("A", sv), "genes")
  # identical(): a list of the same conditions is not a match.
  expect_identical(.pt_default_figure_keys(list("A", "B"), sv), "genes")
  # A start view with no tracks gives NULL.
  expect_null(.pt_default_figure_keys("A", list(conditions = "A")))
})

test_that(".pt_track_split for one condition offers its tracks and the sample-free ones", {
  items <- fake_track_items()
  s <- .pt_track_split("A", items, NULL)
  expect_identical(names(s$available), c("bw_atac_A", "genes", "hic_hic_A", "tpm_usage_Quant_A", "transcripts"))
  expect_identical(names(s$figure), "genes")
  expect_identical(as.character(s$figure$genes), "<strong>Gene annotation</strong>")
  s <- .pt_track_split("B", items, NULL)
  expect_identical(names(s$available), c("bw_atac_B", "genes", "transcripts"))
})

test_that(".pt_track_split for two conditions adds only that pair's comparison, in that order", {
  items <- fake_track_items()
  s <- .pt_track_split(c("A", "B"), items, NULL)
  expect_identical(names(s$available), c("bw_atac_A", "bw_atac_B", "bw_diff_stacked_atac_A_vs_B", "genes",
                                         "hic_hic_A", "tpm_usage_Quant_A", "transcripts"))
  s <- .pt_track_split(c("B", "A"), items, NULL)
  expect_true("bw_diff_stacked_atac_B_vs_A" %in% names(s$available))
  expect_false("bw_diff_stacked_atac_A_vs_B" %in% names(s$available))
})

test_that(".pt_track_split with no or unknown conditions offers only the sample-free tracks", {
  items <- fake_track_items()
  for (conds in list(character(0), "Z", c("Y", "Z"))) {
    s <- .pt_track_split(conds, items, NULL)
    expect_identical(names(s$available), c("genes", "transcripts"), info = paste(conds))
    expect_identical(names(s$figure), "genes", info = paste(conds))
  }
})

test_that(".pt_track_split with an NA condition offers every comparison", {
  # current behaviour, possibly a bug: comparison items have sample = NA, so an NA
  # condition matches every one of them (NA %in% NA). The conditions block never passes
  # an NA condition (it waits for both boxes to be filled), so the app is not affected.
  s <- .pt_track_split(NA_character_, fake_track_items(), NULL)
  expect_identical(names(s$available), c("bw_diff_stacked_atac_B_vs_A", "bw_diff_stacked_atac_A_vs_B",
                                         "genes", "transcripts"))
})

test_that(".pt_track_split starts from a matching start view, dropping tracks outside the pool", {
  items <- fake_track_items()
  sv <- list(conditions = "A", tracks = c("bw_atac_A", "genes", "bw_atac_B", "genes", "no_such_key"))
  s <- .pt_track_split("A", items, sv)
  # Order and repeated Gene annotation kept; B's track and the unknown key dropped.
  expect_identical(names(s$figure), c("bw_atac_A", "genes", "genes"))
  # Gene annotation stays in Available, once.
  expect_identical(names(s$available), c("genes", "hic_hic_A", "tpm_usage_Quant_A", "transcripts"))
  # Other conditions start from Gene annotation alone.
  expect_identical(names(.pt_track_split("B", items, sv)$figure), "genes")
})

test_that(".pt_track_split with no Gene annotation item leaves Figure tracks empty", {
  items <- fake_track_items()
  items$genes <- NULL
  s <- .pt_track_split("A", items, NULL)
  expect_length(s$figure, 0)
  expect_false("genes" %in% names(s$available))
})

test_that(".pt_compose_available_keys drops every two-sample comparison and sorts by label", {
  items <- fake_track_items()
  expect_identical(.pt_compose_available_keys(items),
                   c("bw_atac_A", "bw_atac_B", "genes", "hic_hic_A", "tpm_usage_Quant_A", "transcripts"))
  # An item with a `samples` field of another length is not a comparison.
  items$odd <- list(label = "zzz", samples = c("A", "B", "C"))
  expect_true("odd" %in% .pt_compose_available_keys(items))
  # No items: NULL (Filter() over no names), not character(0).
  expect_null(.pt_compose_available_keys(list()))
})

test_that(".pt_condition_choices puts a blank '(select)' first, then every condition", {
  expect_identical(.pt_condition_choices(c("A", "B")), c("(select)" = "", A = "A", B = "B"))
  expect_identical(.pt_condition_choices(character(0)), c("(select)" = ""))
  expect_identical(.pt_condition_choices(NULL), c("(select)" = ""))
})
