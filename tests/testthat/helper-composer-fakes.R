# Hand-made stand-ins for the composer's startup context (.pt_composer_context()), so the
# composer's lifted UI builders and server blocks can be tested on their own, with no
# genome and no data files. The items have the shape track_detection.R gives them.

# Two samples, A and B: an ATAC coverage track each (B's file missing, so `broken`), a Hi-C
# track and a Transcript expression track for A, the A-vs-B and B-vs-A coverage
# comparisons, and two sample-free tracks (Gene annotation, Transcripts (filtered)).
fake_track_items <- function() {
  per_sample <- function(type, name, sample, broken = FALSE, label = name)
    list(key = paste0(type, "_", name, "_", sample), label = paste0(label, " \u00B7 ", sample),
         type = type, sample = sample, file = paste0("/no/such/", sample, "_", name),
         colour_key = sample, name = name, broken = broken)
  pair <- function(ref, cmp)
    list(key = paste0("bw_diff_stacked_atac_", ref, "_vs_", cmp),
         label = paste0("atac diff (stacked) \u00B7 ", cmp, " (Cmp) vs ", ref, " (Ref)"),
         type = "bw_diff_stacked", sample = NA_character_, file = NA_character_,
         colour_key = NA_character_, name = "atac", samples = c(ref, cmp),
         control_label = ref, treatment_label = cmp,
         control_path = paste0("/no/such/", ref, ".bw"), treatment_path = paste0("/no/such/", cmp, ".bw"))
  fixed <- function(key, label)
    list(key = key, label = label, type = key, sample = NA_character_, file = NA_character_,
         colour_key = NA_character_)
  items <- list(
    per_sample("bw", "atac", "A"),
    per_sample("hic", "hic", "A"),
    per_sample("tpm_usage", "Quant", "A", label = "Transcript expression"),
    per_sample("bw", "atac", "B", broken = TRUE),
    pair("A", "B"),
    pair("B", "A"),
    fixed("genes", "Gene annotation"),
    fixed("transcripts", "Transcripts (filtered)")
  )
  setNames(items, vapply(items, `[[`, "", "key"))
}

# A startup context with the fields the UI builders and server blocks read. `start_view`
# is the settings' start view (NULL: none).
fake_ctx <- function(start_view = NULL, items = fake_track_items()) {
  list(
    TRACK_ITEMS = items,
    AVAILABLE_CONDITIONS = c("A", "B"),
    START_VIEW = start_view,
    COMPOSE_AVAILABLE_KEYS = .pt_compose_available_keys(items),
    GENE_REQUIRED_TYPES = .pt_gene_required_types(),
    DEFAULT_SAMPLE_COLS = c(A = "#E41A1C", B = "#377EB8"),
    CHR_LENGTHS = c(chr1 = 1e7, chr2 = 5e6),
    ALL_GENES = c("Elf5", "Pax5"),
    txdb = "fake-txdb",
    ORGDB_PKG = "org.Fake.eg.db",
    GENOME = list(genome_label = "mm10")
  )
}

# A mock session that records what the server sends the page: custom messages (in
# `$custom`, as list(type, message)) and input updates (in `$input`, as list(id, message)).
recording_session <- function() {
  rec <- new.env()
  rec$custom <- list()
  rec$input <- list()
  session <- shiny::MockShinySession$new()
  session$sendCustomMessage <- function(type, message)
    rec$custom[[length(rec$custom) + 1]] <- list(type = type, message = message)
  session$sendInputMessage <- function(inputId, message)
    rec$input[[length(rec$input) + 1]] <- list(id = inputId, message = message)
  rec$session <- session
  rec
}

# The HTML of a renderUI() output read inside testServer(), as one string.
ui_html <- function(x) paste(as.character(x$html), collapse = "")
