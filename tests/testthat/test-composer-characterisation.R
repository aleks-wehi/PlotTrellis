# Characterisation of the composer server: what each observer, reactive and output does
# now, driven only through input ids, output ids, `state` and the server's named reactives
# (saved_cfg, live_cfg, settings_dirty, chr_len_r, region_valid, resolved_region,
# selected_conditions, has_conditions_ready, has_mode, invalid_figure_tracks,
# figure_pool_key, sel_prev, settings_ui_gen, display_keys). These names and ids stay when
# the server is split into top-level functions, so the same tests check that the split
# changes no behaviour. Nothing draws: render_all() is mocked to write a tiny SVG.
#
# A test marked "current behaviour, possibly a bug" pins what the code does today, not
# what it should do.

need_fixture_data <- function()
  skip_if_not(PT_TEST_DATA_AVAILABLE, "PlotTrellis test fixture data not available on this machine")

# A fresh project (json + its folder) on the mouse fixture, and the composer built for it.
# The same shape as test-composer-settings.R's new_composer(), which lives in that file.
char_composer <- function(extra = list()) {
  proj <- tempfile("pt-char-"); dir.create(proj)
  seed_project_genome(proj)
  json <- file.path(proj, "project_settings.json")
  jsonlite::write_json(c(list(data_dir = FIXTURES$mouse$data_dir, samples = FIXTURES$mouse$samples,
                              gtf_path = GENOMES$mouse$gtf_path, orgdb = GENOMES$mouse$orgdb,
                              genome_label = GENOMES$mouse$genome_label), extra),
                       json, auto_unbox = TRUE, null = "null")
  list(proj = proj, json = json,
       app = suppressMessages(.pt_composer_app(json, .pt_read_settings(json, .pt_default_settings()))))
}

# One composer shared by the tests that never write its project files (each testServer()
# call is its own session, so they share no state). Built on first use.
.char_shared <- new.env()
shared_composer <- function() {
  if (is.null(.char_shared$cp)) {
    .char_shared$cp <- char_composer()
    withr::defer(unlink(.char_shared$cp$proj, recursive = TRUE), envir = testthat::teardown_env())
  }
  .char_shared$cp
}

# Records what the server does: renders (render_all() calls), track-list updates, notifications
# and modals (mocked in the calling test), and, once connect() has hooked the session,
# its custom messages and input updates. Set rec$fail to a message to make the next renders fail.
composer_recorder <- function(env = parent.frame()) {
  rec <- new.env()
  rec$renders <- list(); rec$ranks <- list(); rec$notes <- list(); rec$modals <- list()
  rec$msgs <- list(); rec$inputs <- list(); rec$fail <- NULL
  local_mocked_bindings(
    render_all = function(entries, region_params, settings, device_fn) {
      rec$renders[[length(rec$renders) + 1]] <- list(entries = entries, region = region_params, settings = settings)
      if (!is.null(rec$fail)) stop(rec$fail, call. = FALSE)
      writeLines("<svg width='10' height='10'><rect/></svg>", settings$output_path)
      invisible(NULL)
    },
    update_rank_list = function(css_id, text = NULL, labels = NULL, session = NULL) {
      rec$ranks[[length(rec$ranks) + 1]] <- list(id = css_id, keys = names(labels),
                                                 labels = unname(vapply(labels, as.character, character(1))))
    },
    showNotification = function(ui, ..., type = "default") {
      rec$notes[[length(rec$notes) + 1]] <- list(ui = ui, type = type, args = list(...))
    },
    showModal = function(ui, ...) rec$modals[[length(rec$modals) + 1]] <- as.character(ui),
    .env = env
  )
  rec
}

# Hooks the test session's custom messages and input updates (update*Input()) into `rec`,
# then runs the first flush, as a browser connecting does: until then no observer has run,
# and an observer with ignoreInit = TRUE would take the test's first input change for its
# start. Call it first in every testServer() block.
connect <- function(session, rec) {
  unlockBinding("sendCustomMessage", session)
  session$sendCustomMessage <- function(type, message)
    rec$msgs[[length(rec$msgs) + 1]] <- list(type = type, message = message)
  unlockBinding("sendInputMessage", session)
  session$sendInputMessage <- function(inputId, message)
    rec$inputs[[length(rec$inputs) + 1]] <- list(id = inputId, message = message)
  session$flushReact()
  invisible()
}

last_msg <- function(rec, type) {
  m <- Filter(function(x) identical(x$type, type), rec$msgs)
  if (length(m)) m[[length(m)]]$message else NULL
}
msg_types <- function(rec) vapply(rec$msgs, `[[`, character(1), "type")
# Every value sent to one input by update*Input(), in order.
sent_values <- function(rec, id)
  lapply(Filter(function(x) identical(x$id, id), rec$inputs), function(x) x$message$value)
last_sent <- function(rec, id) { v <- sent_values(rec, id); if (length(v)) v[[length(v)]] else NULL }
# The last update_rank_list() call for one list.
last_rank <- function(rec, id) {
  r <- Filter(function(x) identical(x$id, id), rec$ranks)
  if (length(r)) r[[length(r)]] else NULL
}
last_render <- function(rec) rec$renders[[length(rec$renders)]]
html_of <- function(x) if (is.null(x)) NULL else as.character(x$html)

# A "genes"-only coordinates figure on chr2.
COORDS <- list(view_mode = "compose", region_mode = "coords", region_chr = "chr2",
               region_range = c(103300000L, 103500000L), figure_tracks = "genes")
PARSE_START_ERR <- "Cannot parse start. Accepts: 17000000, 17Mb, 15.5kb (case-insensitive)"
PARSE_END_ERR   <- "Cannot parse end. Accepts: 17000000, 17Mb, 15.5kb (case-insensitive)"
CHR2_LEN <- 182113224
CHRM_LEN <- 16299
COMPOSE_KEYS <- c("bw_atac_Ba", "bw_atac_LP", "bw_atac_P6", "bw_atac_P7", "protein_domains", "genes",
                  "hic_hic_Ba", "hic_hic_LP", "hic_hic_P6", "hic_hic_P7", "junctions_Junctions_Ba",
                  "junctions_Junctions_P6", "junctions_Junctions_P7", "tpm_usage_Ba", "tpm_usage_LP",
                  "transcripts_all", "transcripts", "transcripts_annotated")
SINGLE_BA_KEYS <- c("bw_atac_Ba", "protein_domains", "genes", "hic_hic_Ba", "junctions_Junctions_Ba",
                    "tpm_usage_Ba", "transcripts_all", "transcripts", "transcripts_annotated")

# ---- Region: chromosome, slider and text boxes ------------------------------------------

test_that("chr_select_ui lists main chromosomes first, then contigs, and selects chr1", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    html <- html_of(output$chr_select_ui)
    vals <- regmatches(html, gregexpr('<option value="[^"]*"', html))[[1]]
    vals <- sub('^<option value="', "", sub('"$', "", vals))
    expect_identical(vals[1:22], c(paste0("chr", 1:19), "chrX", "chrY", "chrM"))
    # Every other contig after them, in sort() order (testthat collates in C).
    expect_true(all(grepl("_", vals[-(1:22)])))
    expect_identical(vals[-(1:22)], sort(vals[-(1:22)]))
    expect_match(html, '<option value="chr1" selected>', fixed = TRUE)
    expect_match(html, 'id="region_chr"', fixed = TRUE)
  })
})

test_that("chr_len_r is the chosen chromosome's length, and waits for a chromosome", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    expect_error(chr_len_r(), class = "shiny.silent.error")
    session$setInputs(region_chr = "chr2")
    expect_equal(chr_len_r(), CHR2_LEN)
    session$setInputs(region_chr = "chrM")
    expect_equal(chr_len_r(), CHRM_LEN)
  })
})

test_that("an unknown chromosome gets the default length (bug 110)", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    # Bug 110 (fixed on main 2026-11-17): a chromosome the chrom sizes lack gets the
    # default length (200 Mb for the slider, no upper bound for the validator).
    session$setInputs(region_chr = "chrUnknown", region_range = c(1L, 1000L), region_mode = "coords")
    expect_identical(chr_len_r(), 200000000L)
    expect_true(region_valid()$valid)
    expect_null(output$region_error_ui)
  })
})

test_that("range_slider_ui spans the chromosome in 1 kb steps; the text boxes show its ends", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    session$setInputs(region_chr = "chr2")
    slider <- html_of(output$range_slider_ui)
    expect_match(slider, 'id="region_range"', fixed = TRUE)
    expect_match(slider, 'data-min="1" data-max="182113224" data-from="1" data-to="182113224" data-step="1000"', fixed = TRUE)
    expect_match(slider, 'data-prettify-separator=","', fixed = TRUE)
    expect_match(html_of(output$start_num_ui), 'id="region_start_text"[^>]*value="1 bp"')
    expect_match(html_of(output$end_num_ui), 'id="region_end_text"[^>]*value="182.113Mb"')
    session$setInputs(region_chr = "chrM")
    expect_match(html_of(output$range_slider_ui), 'data-max="16299"', fixed = TRUE)
    expect_match(html_of(output$end_num_ui), 'value="16.299Kb"', fixed = TRUE)
    # The start box never depends on the chromosome.
    expect_match(html_of(output$start_num_ui), 'value="1 bp"', fixed = TRUE)
  })
})

test_that("moving the slider writes both text boxes in coordinate format", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    session$setInputs(region_chr = "chr2", region_range = c(1000000L, 2500000L))
    expect_identical(sent_values(rec, "region_start_text"), list("1Mb"))
    expect_identical(sent_values(rec, "region_end_text"), list("2.5Mb"))
    session$setInputs(region_range = c(999L, 1234567L))
    expect_identical(last_sent(rec, "region_start_text"), "999 bp")
    expect_identical(last_sent(rec, "region_end_text"), "1.235Mb")
    expect_length(sent_values(rec, "region_range"), 0)  # no echo back to the slider
  })
})

test_that("a parsed start or end text moves only its own end of the slider", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    session$setInputs(region_chr = "chr2", region_range = c(1000000L, 2500000L))
    session$setInputs(region_start_text = "1.5 Mb")
    expect_equal(as.numeric(last_sent(rec, "region_range")), c(1500000, 2500000))
    expect_null(state$start_parse_err)
    session$setInputs(region_end_text = "12 kb")
    # The test session never moves the slider itself, so the start is the input's own.
    expect_equal(as.numeric(last_sent(rec, "region_range")), c(1000000, 12000))
    expect_null(state$end_parse_err)
    session$setInputs(region_start_text = "17000000")
    expect_equal(as.numeric(last_sent(rec, "region_range")), c(17000000, 2500000))
  })
})

test_that("a text equal to the slider's value sends no slider update", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    session$setInputs(region_chr = "chr2", region_range = c(1000000L, 2500000L))
    session$setInputs(region_start_text = "1Mb", region_end_text = "2500kb")
    expect_length(sent_values(rec, "region_range"), 0)
  })
})

test_that("start/end beyond the chromosome or crossed over are passed to the slider unchecked", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    session$setInputs(region_chr = "chr2", region_range = c(1000000L, 2500000L))
    session$setInputs(region_end_text = "200Mb")
    expect_equal(as.numeric(last_sent(rec, "region_range")), c(1000000, 200000000))
    expect_null(state$end_parse_err)
    session$setInputs(region_start_text = "5Mb")  # past the slider's end
    expect_equal(as.numeric(last_sent(rec, "region_range")), c(5000000, 2500000))
    expect_null(state$start_parse_err)
    session$setInputs(region_start_text = "0")
    expect_equal(as.numeric(last_sent(rec, "region_range")), c(0, 2500000))
  })
})

test_that("unparseable start or end text sets its parse error and leaves the slider", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    session$setInputs(region_chr = "chr2", region_range = c(1000000L, 2500000L))
    for (bad in c("abc", "1,000", "-5", "1e6", "1.5", "5 gb")) {
      session$setInputs(region_start_text = bad)
      expect_identical(state$start_parse_err, PARSE_START_ERR, info = bad)
      session$setInputs(region_end_text = bad)
      expect_identical(state$end_parse_err, PARSE_END_ERR, info = bad)
    }
    expect_length(sent_values(rec, "region_range"), 0)
    # A good value clears only its own box's error.
    session$setInputs(region_start_text = "2Mb")
    expect_null(state$start_parse_err)
    expect_identical(state$end_parse_err, PARSE_END_ERR)
  })
})

test_that("a blank text is ignored and keeps an earlier parse error", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    session$setInputs(region_chr = "chr2", region_range = c(1000000L, 2500000L))
    session$setInputs(region_start_text = "abc")
    session$setInputs(region_start_text = "")
    expect_identical(state$start_parse_err, PARSE_START_ERR)
    session$setInputs(region_start_text = "   ")
    expect_identical(state$start_parse_err, PARSE_START_ERR)
    expect_length(sent_values(rec, "region_range"), 0)
  })
})

test_that("text typed before the slider exists is ignored, even when unparseable", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    session$setInputs(region_start_text = "abc", region_end_text = "xyz")
    expect_null(state$start_parse_err)
    expect_null(state$end_parse_err)
    expect_length(sent_values(rec, "region_range"), 0)
  })
})

test_that("region_valid accepts a range inside the chromosome and names each fault", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    expect_error(region_valid(), class = "shiny.silent.error")
    session$setInputs(region_chr = "chr2")
    expect_error(region_valid(), class = "shiny.silent.error")  # no range yet
    session$setInputs(region_range = c(1000000L, 2500000L))
    expect_identical(region_valid(), list(valid = TRUE, chrom = "chr2", start = 1000000L, end = 2500000L))
    session$setInputs(region_range = c(1L, as.integer(CHR2_LEN)))
    expect_true(region_valid()$valid)
    session$setInputs(region_range = c(5000L, 4000L))
    expect_identical(region_valid(), list(valid = FALSE, msg = "End must be > start"))
    session$setInputs(region_range = c(5000L, 5000L))
    expect_identical(region_valid(), list(valid = FALSE, msg = "End must be > start"))
    session$setInputs(region_range = c(0L, 4000L))
    expect_identical(region_valid(), list(valid = FALSE, msg = "Start must be \u2265 1"))
    session$setInputs(region_range = c(1L, as.integer(CHR2_LEN) + 1L))
    expect_identical(region_valid(), list(valid = FALSE, msg = "End exceeds chromosome length"))
    # Doubles from the slider come back as integers.
    session$setInputs(region_range = c(1500.7, 3000.2))
    expect_identical(region_valid()[c("start", "end")], list(start = 1500L, end = 3000L))
  })
})

test_that("region_error_ui shows only in Coordinates mode; parse errors win over a bad range", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    session$setInputs(region_chr = "chr2", region_range = c(5000L, 4000L))
    expect_null(output$region_error_ui)  # no region mode yet
    session$setInputs(region_mode = "gene")
    expect_null(output$region_error_ui)
    session$setInputs(region_mode = "coords")
    expect_identical(html_of(output$region_error_ui), '<div class="alert alert-danger">End must be &gt; start</div>')
    session$setInputs(region_start_text = "abc")
    expect_identical(html_of(output$region_error_ui),
                     paste0('<div class="alert alert-warning">', PARSE_START_ERR, "</div>"))
    session$setInputs(region_end_text = "xyz")
    expect_identical(html_of(output$region_error_ui),
                     paste0('<div class="alert alert-warning">', PARSE_START_ERR, " | ", PARSE_END_ERR, "</div>"))
    session$setInputs(region_mode = "gene")
    expect_null(output$region_error_ui)
    session$setInputs(region_mode = "coords", region_start_text = "1kb", region_end_text = "4kb",
                      region_range = c(1000L, 4000L))
    expect_null(output$region_error_ui)
  })
})

# ---- Region: gene lookup ------------------------------------------------------------

test_that("a wrong-case gene symbol resolves, and the box is sent the canonical symbol", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    session$setInputs(gene_name = "elf5")
    expect_identical(state$gene_symbol, "Elf5")
    expect_identical(state$gene_region, list(chrom = "chr2", start = 103410000L, end = 103455000L))
    expect_null(state$gene_error)
    expect_null(state$gene_note)
    expect_identical(sent_values(rec, "gene_name"), list("Elf5"))
  })
})

test_that("an exactly typed symbol sends nothing back to the gene box", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    session$setInputs(gene_name = "Elf5")
    expect_identical(state$gene_symbol, "Elf5")
    expect_length(sent_values(rec, "gene_name"), 0)
  })
})

test_that("the echoed symbol is not looked up again: the typed ID's region and note stay", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    session$setInputs(gene_name = "ENSMUSG00000093931")
    by_id <- state$gene_region
    note <- state$gene_note
    expect_identical(state$gene_symbol, "Amy2a1")
    expect_identical(note, "ENSMUSG00000093931 is also named Amy2a3")
    expect_identical(sent_values(rec, "gene_name"), list("Amy2a1"))
    session$setInputs(gene_name = "Amy2a1")  # the browser's echo
    expect_identical(state$gene_region, by_id)
    expect_identical(state$gene_note, note)
    expect_length(sent_values(rec, "gene_name"), 1)
    # A flank change re-resolves the typed ID, not the echoed symbol: same centre.
    session$setInputs(s_gene_flank_pct = 50)
    wide <- state$gene_region
    expect_gt(wide$end - wide$start, by_id$end - by_id$start)
    expect_equal((wide$start + wide$end) / 2, (by_id$start + by_id$end) / 2, tolerance = 1)
    expect_identical(state$gene_note, note)
  })
})

test_that("an unknown gene sets the error and clears the region", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    session$setInputs(gene_name = "Elf5")
    session$setInputs(gene_name = "Nofoundgene")
    expect_identical(state$gene_error,
                     "'Nofoundgene' not found in mm10: enter a gene symbol (e.g. Pax5) or an Ensembl gene ID")
    expect_null(state$gene_region)
    expect_null(state$gene_note)
    expect_length(sent_values(rec, "gene_name"), 0)
    expect_identical(html_of(output$gene_result_ui),
                     paste0('<div class="alert alert-danger">', state$gene_error, "</div>"))
  })
})

test_that("a blank or whitespace gene box clears the whole lookup", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    session$setInputs(gene_name = "Elf5")
    expect_false(is.null(state$gene_region))
    session$setInputs(gene_name = "   ")
    expect_null(state$gene_region); expect_null(state$gene_symbol)
    expect_null(state$gene_error);  expect_null(state$gene_note)
    expect_null(output$gene_result_ui)
    session$setInputs(gene_name = "Nofoundgene")
    session$setInputs(gene_name = "")
    expect_null(state$gene_error)
    # Leading and trailing spaces are trimmed before the lookup.
    session$setInputs(gene_name = "  Elf5 ")
    expect_identical(state$gene_symbol, "Elf5")
  })
})

test_that("the gene flank comes from s_gene_flank_pct, else the saved setting, else the default", {
  need_fixture_data()
  cp <- char_composer(list(gene_flank_pct = 50))
  on.exit(unlink(cp$proj, recursive = TRUE), add = TRUE)
  rec <- composer_recorder()
  shiny::testServer(cp$app, {
    connect(session, rec)
    session$setInputs(gene_name = "Elf5")
    saved50 <- state$gene_region
    expect_identical(saved50, list(chrom = "chr2", start = 103400000L, end = 103460000L))
    session$setInputs(s_gene_flank_pct = 15)
    expect_identical(state$gene_region, list(chrom = "chr2", start = 103410000L, end = 103455000L))
    for (typed in list(NA, "abc")) {  # an empty box, or junk: the saved value
      session$setInputs(s_gene_flank_pct = typed)
      expect_identical(state$gene_region, saved50, info = deparse(typed))
    }
  })
})

test_that("switching region mode keeps the gene lookup's region, error and box", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    session$setInputs(region_mode = "gene", gene_name = "Elf5")
    found <- html_of(output$gene_result_ui)
    session$setInputs(region_mode = "coords")
    expect_false(is.null(state$gene_region))
    session$setInputs(region_mode = "gene")
    expect_identical(html_of(output$gene_result_ui), found)
    session$setInputs(gene_name = "Nofoundgene")
    err <- state$gene_error
    session$setInputs(region_mode = "coords")
    expect_identical(state$gene_error, err)
    session$setInputs(region_mode = "gene")
    expect_identical(state$gene_error, err)
  })
})

test_that("gene_result_ui shows the found region, its size and any note", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    expect_null(output$gene_result_ui)
    session$setInputs(gene_name = "Elf5")
    expect_identical(html_of(output$gene_result_ui), paste0(
      '<div class="alert alert-success">\n  <strong>Found: </strong>\n',
      "  chr2:103.41Mb \u2013 103.455Mb (45Kb)\n</div>"))
    session$setInputs(gene_name = "ENSMUSG00000093931")
    html <- html_of(output$gene_result_ui)
    expect_match(html, "chr3:113.311Mb \u2013 113.325Mb (14Kb)", fixed = TRUE)
    expect_match(html, "<br/>\n  ENSMUSG00000093931 is also named Amy2a3\n</div>", fixed = TRUE)
  })
})

# ---- Reference / Comparison ---------------------------------------------------------

test_that("picking the other box's condition clears the other box", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    session$setInputs(view_mode = "two", sel_reference = "Ba")
    expect_identical(last_sent(rec, "sel_reference"), "Ba")
    expect_identical(last_sent(rec, "sel_comparison"), "")
    expect_identical(c(sel_prev$control, sel_prev$treatment), c("Ba", ""))
    # Comparison picks Ba too: Reference is cleared.
    session$setInputs(sel_comparison = "Ba")
    expect_identical(last_sent(rec, "sel_reference"), "")
    expect_identical(last_sent(rec, "sel_comparison"), "Ba")
    expect_identical(c(sel_prev$control, sel_prev$treatment), c("", "Ba"))
    session$setInputs(sel_reference = "")  # the browser's echo
    # Reference picks Ba back: Comparison is cleared.
    session$setInputs(sel_reference = "Ba")
    expect_identical(last_sent(rec, "sel_reference"), "Ba")
    expect_identical(last_sent(rec, "sel_comparison"), "")
    expect_identical(c(sel_prev$control, sel_prev$treatment), c("Ba", ""))
  })
})

test_that("a swap in one step and distinct picks keep both boxes", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    session$setInputs(view_mode = "two", sel_reference = "Ba", sel_comparison = "LP")
    expect_identical(c(sel_prev$control, sel_prev$treatment), c("Ba", "LP"))
    session$setInputs(sel_reference = "LP", sel_comparison = "Ba")
    expect_identical(c(sel_prev$control, sel_prev$treatment), c("LP", "Ba"))
    expect_identical(last_sent(rec, "sel_reference"), "LP")
    expect_identical(last_sent(rec, "sel_comparison"), "Ba")
    # Both boxes are sent on every change, each time.
    expect_identical(length(sent_values(rec, "sel_reference")), length(sent_values(rec, "sel_comparison")))
  })
})

test_that("blank picks never collide", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    session$setInputs(sel_reference = "", sel_comparison = "")
    expect_identical(c(sel_prev$control, sel_prev$treatment), c("", ""))
    session$setInputs(sel_reference = "P6", sel_comparison = "P7")
    session$setInputs(sel_reference = "")
    expect_identical(c(sel_prev$control, sel_prev$treatment), c("", "P7"))
    expect_identical(last_sent(rec, "sel_comparison"), "P7")
  })
})

# ---- Conditions ---------------------------------------------------------------------

test_that("has_mode, selected_conditions and has_conditions_ready follow view_mode", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    expect_false(has_mode())
    expect_null(selected_conditions())
    expect_false(has_conditions_ready())
    session$setInputs(sel_condition = "P6", sel_reference = "Ba", sel_comparison = "LP")
    expect_identical(selected_conditions(), "P6")  # no mode reads the single box
    expect_false(has_conditions_ready())
    session$setInputs(view_mode = "single")
    expect_true(has_mode())
    expect_identical(selected_conditions(), "P6")
    expect_true(has_conditions_ready())
    session$setInputs(view_mode = "two")
    expect_identical(selected_conditions(), c("Ba", "LP"))
    expect_true(has_conditions_ready())
    session$setInputs(sel_comparison = "")
    expect_false(has_conditions_ready())
    expect_identical(selected_conditions(), c("Ba", ""))
    session$setInputs(view_mode = "compose")
    expect_true(has_conditions_ready())
    expect_identical(selected_conditions(), "P6")  # Multi still reads the stale single box
    session$setInputs(view_mode = "single", sel_condition = "")
    expect_false(has_conditions_ready())
    session$setInputs(view_mode = "bogus")
    expect_true(has_mode())
    expect_false(has_conditions_ready())
  })
})

test_that("tracks_hint says what to pick first, and goes once the conditions are ready", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    hint <- function(msg) sprintf('<style>#tracks-block { --pt-tracks-hint: "%s"; }</style>', msg)
    expect_identical(html_of(output$tracks_hint), hint("Pick the number of samples above to start."))
    session$setInputs(view_mode = "single")
    expect_identical(html_of(output$tracks_hint), hint("Pick a condition above to see its tracks."))
    session$setInputs(view_mode = "two")
    expect_identical(html_of(output$tracks_hint), hint("Pick a reference and a comparison above to see their tracks."))
    session$setInputs(view_mode = "bogus")
    expect_identical(html_of(output$tracks_hint), hint("Pick the number of samples above to start."))
    session$setInputs(view_mode = "compose")
    expect_null(output$tracks_hint)
    session$setInputs(view_mode = "single", sel_condition = "Ba")
    expect_null(output$tracks_hint)
  })
})

# ---- Track pool -----------------------------------------------------------------------

test_that("at startup both track lists are emptied and the pool key is not_ready", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    expect_identical(figure_pool_key(), "not_ready")
    expect_length(rec$ranks, 2)
    expect_identical(vapply(rec$ranks, `[[`, "", "id"), c("available_tracks", "figure_tracks"))
    expect_null(last_rank(rec, "available_tracks")$keys)
    expect_null(last_rank(rec, "figure_tracks")$keys)
    # A mode with no condition yet is the same key: nothing is sent again. The boxes are
    # "" as the browser sends them (a NULL box would drop out of the pair).
    session$setInputs(view_mode = "single", sel_condition = "", sel_reference = "", sel_comparison = "")
    session$setInputs(view_mode = "two", sel_reference = "Ba")
    expect_identical(figure_pool_key(), "not_ready")
    expect_length(rec$ranks, 2)
  })
})

test_that("One sample fills the pool with that sample's and the sample-free tracks", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    session$setInputs(view_mode = "single", sel_condition = "Ba")
    expect_identical(figure_pool_key(), "single::Ba")
    expect_identical(last_rank(rec, "available_tracks")$keys, SINGLE_BA_KEYS)
    expect_identical(last_rank(rec, "figure_tracks")$keys, "genes")
    expect_identical(last_rank(rec, "available_tracks")$labels[1], "<strong>atac</strong> \u00B7 Ba")
    expect_identical(last_rank(rec, "figure_tracks")$labels, "<strong>Gene annotation</strong>")
  })
})

test_that("Two samples add only that ordered pair's comparison tracks", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    session$setInputs(view_mode = "two", sel_reference = "Ba", sel_comparison = "P6")
    expect_identical(figure_pool_key(), "two::Ba|P6")
    avail <- last_rank(rec, "available_tracks")$keys
    expect_identical(grep("_vs_", avail, value = TRUE),
                     c("bw_diff_stacked_atac_Ba_vs_P6", "junctions_mirrored_diff_Junctions_Ba_vs_P6"))
    expect_true(all(c("bw_atac_Ba", "bw_atac_P6", "genes") %in% avail))
    expect_false(any(grepl("_(LP|P7)$", avail)))
    expect_identical(last_rank(rec, "figure_tracks")$keys, "genes")
    session$setInputs(sel_reference = "P6", sel_comparison = "Ba")
    expect_identical(figure_pool_key(), "two::P6|Ba")
    expect_identical(grep("_vs_", last_rank(rec, "available_tracks")$keys, value = TRUE),
                     c("bw_diff_stacked_atac_P6_vs_Ba", "junctions_mirrored_diff_Junctions_P6_vs_Ba"))
  })
})

test_that("Multi offers every non-comparison track and ignores the condition boxes", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    session$setInputs(view_mode = "compose")
    expect_identical(figure_pool_key(), "compose")
    expect_identical(last_rank(rec, "available_tracks")$keys, COMPOSE_KEYS)
    expect_identical(last_rank(rec, "figure_tracks")$keys, "genes")
    n <- length(rec$ranks)
    session$setInputs(sel_condition = "P6", sel_reference = "Ba", sel_comparison = "LP")
    expect_length(rec$ranks, n)
  })
})

test_that("a trigger that keeps the pool key leaves the lists alone; a new key refills them", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    session$setInputs(view_mode = "single", sel_condition = "Ba")
    n <- length(rec$ranks)
    session$setInputs(sel_reference = "LP", sel_comparison = "P6")  # other mode's boxes
    expect_length(rec$ranks, n)
    session$setInputs(sel_condition = "Ba")
    expect_length(rec$ranks, n)
    session$setInputs(view_mode = "compose")
    expect_length(rec$ranks, n + 2)
    session$setInputs(view_mode = "single")  # back: single::Ba again, refilled
    expect_length(rec$ranks, n + 4)
    expect_identical(last_rank(rec, "available_tracks")$keys, SINGLE_BA_KEYS)
    # Region and gene changes never refill the lists.
    session$setInputs(region_mode = "gene", gene_name = "Elf5")
    session$setInputs(region_mode = "coords", region_chr = "chr2", region_range = c(1000L, 5000L))
    expect_length(rec$ranks, n + 4)
    # Clearing the condition empties both lists.
    session$setInputs(sel_condition = "")
    expect_identical(figure_pool_key(), "not_ready")
    expect_length(last_rank(rec, "available_tracks")$keys, 0)
    expect_length(last_rank(rec, "figure_tracks")$keys, 0)
  })
})

test_that("a Multi start view fills Figure tracks with its kept tracks, in order", {
  need_fixture_data()
  cp <- char_composer()
  on.exit(unlink(cp$proj, recursive = TRUE), add = TRUE)
  .pt_write_start_view(cp$proj, "compose", character(0), "Elf5",
                       c("bw_atac_LP", "genes", "bw_atac_Ba", "genes", "nonexistent_key"))
  app <- suppressMessages(.pt_composer_app(cp$json, .pt_read_settings(cp$json, .pt_default_settings())))
  rec <- composer_recorder()
  shiny::testServer(app, {
    connect(session, rec)
    session$setInputs(view_mode = "compose")
    expect_identical(last_rank(rec, "figure_tracks")$keys, c("bw_atac_LP", "genes", "bw_atac_Ba", "genes"))
    expect_identical(last_rank(rec, "available_tracks")$keys, setdiff(COMPOSE_KEYS, c("bw_atac_LP", "bw_atac_Ba")))
  })
})

test_that("a One-sample start view applies to its own condition only", {
  need_fixture_data()
  cp <- char_composer()
  on.exit(unlink(cp$proj, recursive = TRUE), add = TRUE)
  .pt_write_start_view(cp$proj, "single", "P6", "Elf5", c("bw_atac_P6", "genes", "hic_hic_P6"))
  app <- suppressMessages(.pt_composer_app(cp$json, .pt_read_settings(cp$json, .pt_default_settings())))
  rec <- composer_recorder()
  shiny::testServer(app, {
    connect(session, rec)
    session$setInputs(view_mode = "single", sel_condition = "P6")
    expect_identical(last_rank(rec, "figure_tracks")$keys, c("bw_atac_P6", "genes", "hic_hic_P6"))
    expect_identical(last_rank(rec, "available_tracks")$keys,
                     c("protein_domains", "genes", "junctions_Junctions_P6", "transcripts_all",
                       "transcripts", "transcripts_annotated"))
    session$setInputs(sel_condition = "Ba")
    expect_identical(last_rank(rec, "figure_tracks")$keys, "genes")
    expect_identical(last_rank(rec, "available_tracks")$keys, SINGLE_BA_KEYS)
    session$setInputs(view_mode = "compose")
    expect_identical(last_rank(rec, "figure_tracks")$keys, "genes")
    expect_identical(last_rank(rec, "available_tracks")$keys, COMPOSE_KEYS)
  })
})

# ---- Region resolution and locks ------------------------------------------------------

test_that("resolved_region follows region_mode, never view_mode", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    expect_null(resolved_region())
    session$setInputs(region_chr = "chr2", region_range = c(1000L, 5000L))
    expect_null(resolved_region())  # no region mode: a gene is needed
    session$setInputs(gene_name = "Elf5")
    # current behaviour, possibly a bug: with no region mode chosen yet, a found gene
    # already resolves the region (only "coords" is checked for).
    expect_identical(resolved_region(),
                     list(rgn = list(chrom = "chr2", start = 103410000L, end = 103455000L), gene_name = "Elf5"))
    session$setInputs(region_mode = "coords")
    expect_identical(resolved_region(),
                     list(rgn = list(chrom = "chr2", start = 1000L, end = 5000L), gene_name = NA_character_))
    session$setInputs(region_range = c(5000L, 1000L))
    expect_null(resolved_region())
    session$setInputs(region_mode = "gene")
    expect_identical(resolved_region()$gene_name, "Elf5")
    for (mode in c("single", "two", "compose")) {
      session$setInputs(view_mode = mode)
      expect_identical(resolved_region()$rgn$start, 103410000L, info = mode)
    }
    session$setInputs(gene_name = "Nofoundgene")
    expect_null(resolved_region())
  })
})

test_that("condition_lock locks the region controls until a mode and its conditions are set", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    expect_identical(last_msg(rec, "condition_lock"), list(locked = TRUE))
    session$setInputs(view_mode = "single")
    expect_identical(last_msg(rec, "condition_lock"), list(locked = TRUE))
    session$setInputs(sel_condition = "Ba")
    expect_identical(last_msg(rec, "condition_lock"), list(locked = FALSE))
    session$setInputs(view_mode = "two", sel_reference = "Ba")
    expect_identical(last_msg(rec, "condition_lock"), list(locked = TRUE))
    session$setInputs(sel_comparison = "LP")
    expect_identical(last_msg(rec, "condition_lock"), list(locked = FALSE))
    session$setInputs(view_mode = "compose")
    expect_identical(last_msg(rec, "condition_lock"), list(locked = FALSE))
  })
})

test_that("region_lock locks the track picker until mode, conditions and a region are set", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    expect_identical(last_msg(rec, "region_lock"), list(locked = TRUE))
    session$setInputs(view_mode = "compose")
    expect_identical(last_msg(rec, "region_lock"), list(locked = TRUE))
    session$setInputs(region_mode = "coords", region_chr = "chr2", region_range = c(1000L, 5000L))
    expect_identical(last_msg(rec, "region_lock"), list(locked = FALSE))
    session$setInputs(region_range = c(5000L, 1000L))
    expect_identical(last_msg(rec, "region_lock"), list(locked = TRUE))
    session$setInputs(region_mode = "gene", gene_name = "Elf5")
    expect_identical(last_msg(rec, "region_lock"), list(locked = FALSE))
    session$setInputs(view_mode = "single")  # no condition picked
    expect_identical(last_msg(rec, "region_lock"), list(locked = TRUE))
    session$setInputs(sel_condition = "Ba")
    expect_identical(last_msg(rec, "region_lock"), list(locked = FALSE))
  })
})

# ---- Invalid tracks and the Generate button -------------------------------------------

test_that("invalid_figure_tracks flags gene-required tracks in Coordinates mode only", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    expect_identical(invalid_figure_tracks(), character(0))
    keys <- c("genes", "transcripts", "tpm_usage_Ba", "protein_domains", "transcripts_all",
              "transcripts_annotated", "bw_atac_Ba", "junctions_mirrored_diff_Junctions_Ba_vs_P6", "nonexistent")
    session$setInputs(figure_tracks = keys)
    expect_length(invalid_figure_tracks(), 0)  # no region mode
    session$setInputs(region_mode = "gene")
    expect_length(invalid_figure_tracks(), 0)
    session$setInputs(region_mode = "coords")
    expect_identical(unlist(invalid_figure_tracks()),
                     c("transcripts", "tpm_usage_Ba", "protein_domains", "transcripts_all", "transcripts_annotated"))
    session$setInputs(figure_tracks = c("genes", "bw_atac_Ba"))
    expect_length(invalid_figure_tracks(), 0)
  })
})

test_that("generate_warning_ui names the blocking tracks by label", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    expect_null(output$generate_warning_ui)
    session$setInputs(region_mode = "coords", figure_tracks = c("genes", "transcripts", "tpm_usage_Ba"))
    expect_identical(html_of(output$generate_warning_ui), paste0(
      '<div class="alert alert-warning" style="margin-bottom:0;">\n',
      "  <strong>Remove or fix before generating: </strong>\n",
      "  Transcripts (filtered), Transcript expression \u00B7 Ba\n</div>"))
    session$setInputs(region_mode = "gene")
    expect_null(output$generate_warning_ui)
  })
})

test_that("generate_ui is enabled only with a mode, its conditions, a region and valid tracks", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    enabled <- function() {
      html <- html_of(output$generate_ui)
      expect_match(html, "Generate figure \u2192", fixed = TRUE)
      expect_match(html, "generate-btn", fixed = TRUE)
      grepl('id="generate"', html, fixed = TRUE) && !grepl("disabled", html, fixed = TRUE)
    }
    expect_false(enabled())
    session$setInputs(view_mode = "compose")
    expect_false(enabled())
    session$setInputs(region_mode = "coords", region_chr = "chr2", region_range = c(1000L, 5000L))
    expect_false(enabled())  # no tracks
    session$setInputs(figure_tracks = "genes")
    expect_true(enabled())
    expect_identical(html_of(output$generate_ui), paste0(
      '<button class="btn btn-default action-button btn-success btn-lg generate-btn" id="generate" ',
      'type="button"><span class="action-label">Generate figure \u2192</span></button>'))
    session$setInputs(figure_tracks = c("genes", "transcripts"))
    expect_false(enabled())
    expect_identical(html_of(output$generate_ui),
                     '<button class="btn btn-success btn-lg generate-btn" disabled="disabled">Generate figure \u2192</button>')
    session$setInputs(figure_tracks = "genes", region_range = c(5000L, 1000L))
    expect_false(enabled())
    session$setInputs(region_range = c(1000L, 5000L), view_mode = "single")
    expect_false(enabled())
    session$setInputs(sel_condition = "Ba")
    expect_true(enabled())
    session$setInputs(view_mode = "two", sel_reference = "Ba", sel_comparison = "LP")
    expect_true(enabled())
    session$setInputs(region_mode = "gene")
    expect_false(enabled())  # no gene yet
    session$setInputs(gene_name = "Elf5")
    expect_true(enabled())
  })
})

# ---- Generate ---------------------------------------------------------------------------

test_that("Generate renders the picked tracks for the resolved region with the trimmed title", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    do.call(session$setInputs, utils::modifyList(COORDS, list(figure_tracks = c("genes", "bw_atac_Ba"))))
    session$setInputs(plot_title = "  My title  ", generate = 1)
    expect_length(rec$renders, 1)
    r <- last_render(rec)
    expect_identical(r$region, list(chrom = "chr2", start = 103300000L, end = 103500000L, gene = NA_character_))
    expect_identical(r$entries[[1]], list(viz = "genes", height = NULL, data = NULL))
    expect_identical(r$entries[[2]][c("viz", "sample", "colour", "name")],
                     list(viz = "bw", sample = "Ba", colour = "#E41A1C", name = "atac"))
    expect_identical(r$settings$title, "My title")
    expect_match(r$settings$output_path, "\\.svg$")
    expect_identical(r$settings$fontsize, 8L)
    expect_identical(state$svg_path, r$settings$output_path)
    expect_true(file.exists(state$svg_path))
    expect_identical(state$last_region, r$region)
    expect_identical(state$last_labels, c("genes", "bw_atac_Ba"))
    expect_identical(state$last_title, "My title")
    expect_identical(state$rendered_gene, NA_character_)
    expect_null(state$plot_error)
    expect_false(state$rendering)
    expect_match(state$rcode, "^# PlotTrellis figure script")
    expect_match(state$rcode, "# Region: chr2:103,300,000-103,500,000 (mm10)", fixed = TRUE)
    expect_identical(tail(msg_types(rec), 2), c("show_output", "hide_spinner"))
  })
})

test_that("Generate in Gene mode carries the gene to the render and the Gene info panel", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    session$setInputs(view_mode = "compose", region_mode = "gene", gene_name = "elf5", figure_tracks = "genes")
    session$setInputs(gene_name = "Elf5")  # the echo
    session$setInputs(generate = 1)
    expect_identical(last_render(rec)$region,
                     list(chrom = "chr2", start = 103410000L, end = 103455000L, gene = "Elf5"))
    expect_identical(state$rendered_gene, "Elf5")
    expect_identical(last_render(rec)$settings$title, "")  # no title box value
  })
})

test_that("Generate uses the live settings, comparison role colours included", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    session$setInputs(view_mode = "two", sel_reference = "Ba", sel_comparison = "LP", region_mode = "coords",
                      region_chr = "chr2", region_range = c(103300000L, 103500000L),
                      figure_tracks = c("bw_diff_stacked_atac_Ba_vs_LP", "genes"))
    session$setInputs(generate = 1)
    e <- last_render(rec)$entries[[1]]
    expect_identical(e[c("viz", "control_label", "treatment_label", "control_colour", "treatment_colour")],
                     list(viz = "bw_diff_stacked", control_label = "Ba", treatment_label = "LP",
                          control_colour = "#4D4D4D", treatment_colour = "#56B4E9"))
    session$setInputs(s_reference_colour = "#123456", s_fontsize = 14, generate = 2)
    expect_identical(last_render(rec)$entries[[1]]$control_colour, "#123456")
    expect_equal(last_render(rec)$settings$fontsize, 14)
  })
})

test_that("Generate does nothing without a mode, tracks, a region, or a complete pair", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    session$setInputs(region_mode = "coords", region_chr = "chr2", region_range = c(1000L, 5000L),
                      figure_tracks = "genes", generate = 1)
    expect_length(rec$renders, 0)  # no mode
    session$setInputs(view_mode = "compose", figure_tracks = character(0), generate = 2)
    expect_length(rec$renders, 0)  # no tracks
    session$setInputs(figure_tracks = "genes", region_range = c(5000L, 1000L), generate = 3)
    expect_length(rec$renders, 0)  # no valid region
    session$setInputs(region_range = c(1000L, 5000L), view_mode = "two", sel_reference = "Ba", generate = 4)
    expect_length(rec$renders, 0)  # no comparison
    session$setInputs(sel_comparison = "LP", generate = 5)
    expect_length(rec$renders, 1)
  })
})

test_that("Generate in One-sample mode renders even with no condition picked", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    session$setInputs(view_mode = "single", region_mode = "coords", region_chr = "chr2",
                      region_range = c(1000L, 5000L), figure_tracks = "genes")
    expect_match(html_of(output$generate_ui), "disabled", fixed = TRUE)
    # current behaviour, possibly a bug: the Generate observer checks the pair in Two-sample
    # mode only, so a One-sample generate event with no condition still renders (the button
    # itself is disabled, so only a stray event reaches this).
    session$setInputs(generate = 1)
    expect_length(rec$renders, 1)
  })
})

test_that("nothing renders while a render is in progress", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    do.call(session$setInputs, COORDS)
    session$setInputs(generate = 1)
    expect_length(rec$renders, 1)
    state$rendering <- TRUE
    session$setInputs(generate = 2, zoom_in_3x = 1, apply_settings = 1)
    expect_length(rec$renders, 1)
    state$rendering <- FALSE
    session$setInputs(generate = 3)
    expect_length(rec$renders, 2)
  })
})

test_that("a failed render shows its error and keeps the previous figure; the next success clears it", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    do.call(session$setInputs, COORDS)
    session$setInputs(plot_title = "First", generate = 1)
    first <- state$svg_path
    rec$fail <- "boom"
    session$setInputs(plot_title = "Second", generate = 2)
    expect_identical(state$plot_error, "boom")
    expect_false(file.exists(last_render(rec)$settings$output_path))
    expect_identical(state$svg_path, first)
    expect_true(file.exists(first))
    expect_identical(state$last_title, "First")
    expect_false(state$rendering)
    expect_identical(tail(msg_types(rec), 2), c("show_output", "hide_spinner"))
    expect_identical(html_of(output$render_error_ui), paste0(
      '<div class="alert alert-danger">\n  <strong>Error: </strong>\n  boom\n  <br/>\n',
      '  <pre style="font-size:11px;margin-top:8px;">boom</pre>\n</div>'))
    expect_match(html_of(output$svg_display_ui), "<svg ", fixed = TRUE)  # the first figure stays
    rec$fail <- NULL
    session$setInputs(generate = 3)
    expect_null(state$plot_error)
    expect_false(file.exists(first))  # the replaced figure's file is removed
    expect_identical(state$last_title, "Second")
    expect_error(output$render_error_ui, class = "shiny.silent.error")
  })
})

test_that("a failed first render leaves no figure, script or zoom", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    do.call(session$setInputs, COORDS)
    rec$fail <- "first failure"
    session$setInputs(generate = 1)
    expect_identical(state$plot_error, "first failure")
    expect_null(state$svg_path)
    expect_null(state$rcode)
    expect_null(state$last_region)
    expect_match(html_of(output$render_error_ui), "first failure", fixed = TRUE)
    expect_error(output$svg_display_ui, class = "shiny.silent.error")
    expect_error(output$zoom_ui, class = "shiny.silent.error")
    expect_identical(output$rcode_output, "")
  })
})

# ---- Zoom and Apply to figure ---------------------------------------------------------

test_that("each zoom button rescales the last region about its centre", {
  need_fixture_data()
  rec <- composer_recorder()
  want <- list(zoom_in_15x = c(103333333L, 103466667L), zoom_in_3x = c(103366667L, 103433333L),
               zoom_in_10x = c(103390000L, 103410000L), zoom_out_15x = c(103250000L, 103550000L),
               zoom_out_3x = c(103100000L, 103700000L), zoom_out_10x = c(102400000L, 104400000L))
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    do.call(session$setInputs, COORDS)
    i <- 0
    for (btn in names(want)) {
      i <- i + 1
      session$setInputs(generate = i)  # back to the same base region
      do.call(session$setInputs, setNames(list(i), btn))
      r <- last_render(rec)$region
      expect_identical(c(r$start, r$end), want[[btn]], info = btn)
      expect_identical(r$chrom, "chr2", info = btn)
    }
  })
})

test_that("a zoom moves the slider and both text boxes to the new region", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    do.call(session$setInputs, COORDS)
    session$setInputs(generate = 1)
    rec$inputs <- list()
    session$setInputs(zoom_in_3x = 1)
    expect_equal(as.numeric(last_sent(rec, "region_range")), c(103366667, 103433333))
    expect_identical(last_sent(rec, "region_start_text"), "103.367Mb")
    expect_identical(last_sent(rec, "region_end_text"), "103.433Mb")
  })
})

test_that("zoom out clamps at the chromosome's start and end", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    do.call(session$setInputs, utils::modifyList(COORDS, list(region_range = c(1L, 100001L))))
    session$setInputs(generate = 1, zoom_out_10x = 1)
    expect_identical(last_render(rec)$region[c("start", "end")], list(start = 1L, end = 550001L))
    session$setInputs(region_chr = "chrM", region_range = c(10001L, 16001L), generate = 2, zoom_out_10x = 2)
    expect_identical(last_render(rec)$region[c("chrom", "start", "end")],
                     list(chrom = "chrM", start = 1L, end = as.integer(CHRM_LEN)))
  })
})

test_that("zoom in keeps at least 1 kb in view", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    do.call(session$setInputs, utils::modifyList(COORDS, list(region_range = c(100001L, 101001L))))
    session$setInputs(generate = 1, zoom_in_10x = 1)
    expect_identical(last_render(rec)$region[c("start", "end")], list(start = 100001L, end = 101001L))
    session$setInputs(region_range = c(100001L, 100101L), generate = 2, zoom_in_15x = 1)
    expect_identical(last_render(rec)$region[c("start", "end")], list(start = 99551L, end = 100551L))
  })
})

test_that("zoom does nothing before a figure exists", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    do.call(session$setInputs, COORDS)
    rec$inputs <- list()
    session$setInputs(zoom_in_3x = 1, zoom_out_10x = 1, apply_settings = 1)
    expect_length(rec$renders, 0)
    expect_length(sent_values(rec, "region_range"), 0)
  })
})

test_that("zoom and Apply keep the last figure's tracks, gene and title, not the current inputs", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    session$setInputs(view_mode = "compose", region_mode = "gene", gene_name = "Elf5",
                      figure_tracks = c("genes", "bw_atac_Ba"), plot_title = "Kept")
    session$setInputs(generate = 1)
    session$setInputs(figure_tracks = "genes", plot_title = "Changed", gene_name = "Pax5")
    session$setInputs(zoom_in_15x = 1)
    z <- last_render(rec)
    expect_identical(vapply(z$entries, `[[`, "", "viz"), c("genes", "bw"))
    expect_identical(z$settings$title, "Kept")
    expect_identical(z$region$gene, "Elf5")
    expect_identical(state$rendered_gene, "Elf5")
    session$setInputs(s_fontsize = 11, apply_settings = 1)
    a <- last_render(rec)
    expect_identical(a$region, z$region)  # Apply is a zoom by 1
    expect_equal(a$settings$fontsize, 11)
    expect_identical(a$settings$title, "Kept")
  })
})

test_that("apply_settings_ui and zoom_ui appear once a figure exists", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    expect_error(output$apply_settings_ui, class = "shiny.silent.error")
    expect_error(output$zoom_ui, class = "shiny.silent.error")
    do.call(session$setInputs, COORDS)
    session$setInputs(generate = 1)
    expect_identical(html_of(output$apply_settings_ui), paste0(
      '<button class="btn btn-default action-button btn-primary generate-btn" id="apply_settings" ',
      'type="button"><span class="action-label">Apply to figure</span></button>'))
    zoom <- html_of(output$zoom_ui)
    ids <- regmatches(zoom, gregexpr('id="zoom_[a-z0-9_]+"', zoom))[[1]]
    expect_identical(ids, sprintf('id="%s"', c("zoom_in_15x", "zoom_in_3x", "zoom_in_10x",
                                               "zoom_out_15x", "zoom_out_3x", "zoom_out_10x")))
    expect_match(zoom, '<div class="zoom-group-label">Zoom in</div>', fixed = TRUE)
    expect_match(zoom, '<div class="zoom-group-label">Zoom out</div>', fixed = TRUE)
    expect_identical(lengths(regmatches(zoom, gregexpr("1.5\u00D7", zoom, fixed = TRUE))), 2L)
  })
})

# ---- Output page and downloads ----------------------------------------------------------

test_that("before any figure the output page is empty and the downloads refuse", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    expect_error(output$render_error_ui, class = "shiny.silent.error")
    expect_error(output$svg_display_ui, class = "shiny.silent.error")
    expect_identical(output$rcode_output, "")
    expect_error(output$dl_svg, class = "shiny.silent.error")
    expect_error(output$dl_rcode, class = "shiny.silent.error")
  })
})

test_that("a rendered figure is inlined, its script shown, and both download", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    do.call(session$setInputs, COORDS)
    session$setInputs(generate = 1)
    expect_identical(html_of(output$svg_display_ui), paste0(
      '<div style="width:100%; max-width:1400px; margin:0 auto; overflow-x:auto;">',
      "<svg style='width:100%;height:auto;' xml:space='preserve' width='10' height='10'><rect/></svg></div>"))
    expect_identical(output$rcode_output, state$rcode)
    svg <- output$dl_svg
    expect_match(basename(svg), "^plottrellis_[0-9]{8}_[0-9]{6}\\.svg$")
    expect_identical(readLines(svg), "<svg xml:space='preserve' width='10' height='10'><rect/></svg>")
    expect_true(file.exists(state$svg_path))  # a copy, not a move
    rcode <- output$dl_rcode
    expect_match(basename(rcode), "^plottrellis_[0-9]{8}_[0-9]{6}\\.R$")
    expect_identical(readLines(rcode), strsplit(state$rcode, "\n", fixed = TRUE)[[1]])
  })
})

test_that("the downloaded SVG has svglite's text width lock removed; the shown one keeps it", {
  need_fixture_data()
  rec <- composer_recorder()
  local_mocked_bindings(render_all = function(entries, region_params, settings, device_fn) {
    writeLines("<text x='1' textLength='5px' lengthAdjust='spacingAndGlyphs'>a</text>", settings$output_path)
    invisible(NULL)
  })
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    do.call(session$setInputs, COORDS)
    session$setInputs(generate = 1)
    expect_match(readLines(state$svg_path), "textLength='5px'", fixed = TRUE)
    expect_identical(readLines(output$dl_svg), "<text x='1'>a</text>")
  })
})

test_that("Cancel in a session with no render process only hides the spinner", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    rec$msgs <- list()
    session$setInputs(cancel_render = 1)
    expect_identical(msg_types(rec), "hide_spinner")
    expect_length(rec$notes, 0)
    # The test session renders in-process, so there is nothing to cancel: the flag stays.
    state$rendering <- TRUE
    session$setInputs(cancel_render = 2)
    expect_true(state$rendering)
    expect_length(rec$notes, 0)
  })
})

# ---- Figure settings ------------------------------------------------------------------

test_that("display_keys are the non-core settings, and live_cfg starts as saved_cfg", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    expect_identical(display_keys, setdiff(names(.pt_default_settings()), .PT_CORE_SETTINGS))
    expect_identical(live_cfg()[display_keys], saved_cfg()[display_keys])
    expect_false(settings_dirty())
    expect_null(output$settings_dirty_ui)
    expect_identical(settings_ui_gen(), 0)
    expect_identical(saved_cfg()$cache_dir, file.path(shared_composer()$proj, "cache"))
  })
})

test_that("settings_dirty compares numbers by value, everything else exactly", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    session$setInputs(s_fontsize = 8, s_page_w = 9L)  # saved 8L and 9
    expect_false(settings_dirty())
    session$setInputs(s_page_w = 9 + 1e-10)
    expect_false(settings_dirty())
    session$setInputs(s_page_w = 9.01)
    expect_true(settings_dirty())
    expect_identical(html_of(output$settings_dirty_ui),
                     '<span class="text-warning">Unsaved changes (this session only).</span>')
    session$setInputs(s_page_w = 9, s_zrange_cap = NA, s_hic_height = NA)  # empty optional boxes are "auto"
    expect_false(settings_dirty())
    session$setInputs(s_font = " DejaVu Sans ")  # trimmed
    expect_false(settings_dirty())
    session$setInputs(s_fontsize = "8")  # a string is not the number
    expect_true(settings_dirty())
    session$setInputs(s_fontsize = 8, s_show_title = FALSE)
    expect_true(settings_dirty())
    session$setInputs(s_show_title = TRUE)
    expect_false(settings_dirty())
    expect_null(output$settings_dirty_ui)
  })
})

test_that("Save writes the settings, updates saved_cfg and confirms", {
  need_fixture_data()
  cp <- char_composer()
  on.exit(unlink(cp$proj, recursive = TRUE), add = TRUE)
  rec <- composer_recorder()
  before <- jsonlite::read_json(cp$json)
  shiny::testServer(cp$app, {
    connect(session, rec)
    session$setInputs(s_fontsize = 12)
    expect_true(settings_dirty())
    session$setInputs(save_settings = 1)
    expect_length(rec$notes, 1)
    expect_identical(rec$notes[[1]][c("ui", "type")], list(ui = "Saved as the project default.", type = "message"))
    expect_equal(saved_cfg()$fontsize, 12)
    expect_false(settings_dirty())
    expect_null(output$settings_dirty_ui)
  })
  after <- jsonlite::read_json(cp$json)
  expect_equal(after$fontsize, 12)
  for (k in .PT_CORE_SETTINGS) expect_identical(after[[k]], before[[k]], info = k)
  expect_true(file.exists(file.path(cp$proj, "project_samples.txt")))
})

test_that("Save with the settings file gone reports an error and writes nothing", {
  need_fixture_data()
  cp <- char_composer()
  on.exit(unlink(cp$proj, recursive = TRUE), add = TRUE)
  rec <- composer_recorder()
  shiny::testServer(cp$app, {
    connect(session, rec)
    unlink(cp$json)
    session$setInputs(s_fontsize = 14, save_settings = 1)
    expect_length(rec$notes, 1)
    n <- rec$notes[[1]]
    expect_identical(n$type, "error")
    expect_identical(n$ui, paste0("Couldn't save the settings: ", cp$json, " no longer exists; restart the app."))
    expect_true("duration" %in% names(n$args))
    expect_null(n$args$duration)  # stays until closed
    expect_equal(saved_cfg()$fontsize, 8L)
    expect_true(settings_dirty())
  })
  expect_false(file.exists(cp$json))
  expect_false(file.exists(file.path(cp$proj, "project_samples.txt")))
})

test_that("Revert re-renders the settings panel from saved_cfg", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    expect_match(html_of(output$display_settings_ui), '<input id="s_fontsize" type="number"[^>]* value="8"')
    session$setInputs(s_fontsize = 14)
    expect_true(settings_dirty())
    session$setInputs(revert_settings = 1)
    expect_identical(settings_ui_gen(), 1)
    expect_match(html_of(output$display_settings_ui), '<input id="s_fontsize" type="number"[^>]* value="8"')
    # The panel is rebuilt; the session's input keeps 14 until the browser reports back.
    expect_true(settings_dirty())
    session$setInputs(revert_settings = 2)
    expect_identical(settings_ui_gen(), 2)
  })
})

test_that("the settings panel holds every s_ input and the comparison group", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    html <- html_of(output$display_settings_ui)
    for (k in display_keys) expect_match(html, sprintf('id="s_%s"', k), fixed = TRUE, info = k)
    expect_match(html, "Comparison (diff) tracks", fixed = TRUE)
  })
})

# ---- Small server pieces ----------------------------------------------------------------

test_that("a pt-once click is released after the flush that handled it", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    session$setInputs(pt_once_click = "save_settings")
    expect_identical(last_msg(rec, "pt_once_ready"), "save_settings")
  })
})

test_that("Edit project settings and Close app warn about unsaved Figure settings only when dirty", {
  need_fixture_data()
  rec <- composer_recorder()
  shiny::testServer(shared_composer()$app, {
    connect(session, rec)
    session$setInputs(open_settings = 1)
    expect_false(grepl("not saved and will be lost", rec$modals[[1]], fixed = TRUE))
    session$setInputs(s_fontsize = 14, open_settings = 2)
    expect_match(rec$modals[[2]], "Your Figure settings changes are not saved and will be lost.", fixed = TRUE)
    session$setInputs(close_app_btn = 1)
    expect_match(rec$modals[[3]], "not saved and will be lost", fixed = TRUE)
  })
})
