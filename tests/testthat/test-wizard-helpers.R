# Unit tests, with edge cases, for the setup wizard's helper functions (R/wizard-checks.R,
# R/util-settings.R, R/wizard-ui.R, R/wizard-genome.R, R/wizard-dirty.R) and the
# settings readers and defaults it builds on (R/app-main.R). Cases other test files already cover
# are not repeated here. No fixture data needed.
#
# The helpers defined inside wizard_server() are reached through `h`, inside a
# testServer() block (whose evaluation frame holds the server's variables).
h <- function(name) get(name, envir = parent.frame())

hp_server <- function(proj) {
  json <- file.path(proj, "project_settings.json")
  defaults <- .pt_default_settings()
  wizard_server(json, .pt_read_settings(json, defaults), defaults)
}
hp_dir <- function(env = parent.frame()) {
  d <- tempfile("pt-wizh-")
  dir.create(d)
  withr::defer(unlink(d, recursive = TRUE), envir = env)
  d
}
hp_touch <- function(root, files) {
  for (f in files) {
    dir.create(dirname(file.path(root, f)), recursive = TRUE, showWarnings = FALSE)
    file.create(file.path(root, f))
  }
}
hp_html <- function(x) paste(as.character(x), collapse = "")
hp_count <- function(html, pattern) lengths(regmatches(html, gregexpr(pattern, html, fixed = TRUE)))
hp_gtf <- function(lines, ext = ".gtf", env = parent.frame()) {
  f <- tempfile("pt-wizh-gtf-", fileext = ext)
  con <- if (grepl("\\.gz$", ext, ignore.case = TRUE)) gzfile(f, "w") else file(f, "w")
  writeLines(lines, con)
  close(con)
  withr::defer(unlink(f), envir = env)
  f
}
hp_gene <- function(id, attr_sep = " ") sprintf('chr1\tsrc\tgene\t1\t100\t.\t+\t.\tgene_id%s"%s"; gene_name "G";', attr_sep, id)
MOUSE <- list(species = "mouse", orgdb = "org.Mm.eg.db")
HUMAN <- list(species = "human", orgdb = "org.Hs.eg.db")

# ---- wizard_scoped_choices() -------------------------------------------------------------

test_that("wizard_scoped_choices(): no folder, a blank one or a missing one lists nothing", {
  expect_identical(wizard_scoped_choices(NULL, "bw"), character(0))
  expect_identical(wizard_scoped_choices("", "bw"), character(0))
  expect_identical(wizard_scoped_choices(file.path(tempdir(), "no-such-folder"), "quantdir"), character(0))
})

test_that("wizard_scoped_choices(): quantdir lists the folders holding a quant.sf, nested ones too", {
  root <- hp_dir()
  hp_touch(root, c("s1/quant.sf", "deep/er/s2/quant.sf", "s3/QUANT.SF", "s4/quant.sf.gz", "s5/other.txt", "quant.sf"))
  # current behaviour, possibly a bug: a quant.sf in the data folder itself is listed as "",
  # the same value as the picker's blank choice.
  expect_setequal(wizard_scoped_choices(root, "quantdir"), c("", "deep/er/s2", "s1"))
  expect_length(wizard_scoped_choices(root, "quantdir"), 3)
})

test_that("wizard_scoped_choices(): each kind's extensions, ignoring case, in subfolders too", {
  root <- hp_dir()
  hp_touch(root, c("a.bw", "B.BW", "c.bigWig", "sub/d.bigwig", "e.bw.bak",
                   "f.hic", "g.cool", "h.mcool", "i.HIC",
                   "j.broadPeak", "k.gappedPeak", "l.BED", "m.narrowPeak",
                   "n.bedpe", "o.bed12", "p.BEDPE.GZ"))
  # (Sets, not sequences: list.files() sorts by the locale.)
  expect_setequal(wizard_scoped_choices(root, "bw"), c("a.bw", "B.BW", "c.bigWig", "sub/d.bigwig"))
  expect_setequal(wizard_scoped_choices(root, "hic"), c("f.hic", "g.cool", "h.mcool", "i.HIC"))
  expect_setequal(wizard_scoped_choices(root, "bed"), c("j.broadPeak", "k.gappedPeak", "l.BED", "m.narrowPeak"))
  expect_setequal(wizard_scoped_choices(root, "paired-end"), c("l.BED", "n.bedpe", "o.bed12", "p.BEDPE.GZ"))
})

test_that("wizard_scoped_choices(): a trailing slash on the folder gives the same relative paths", {
  root <- hp_dir()
  hp_touch(root, c("x/a.bw", "q/quant.sf"))
  expect_identical(wizard_scoped_choices(paste0(root, "/"), "bw"), "x/a.bw")
  expect_identical(wizard_scoped_choices(paste0(root, "/"), "quantdir"), "q")
})

test_that("wizard_scoped_choices(): a type that isn't a known kind lists every file", {
  root <- hp_dir()
  hp_touch(root, c("a.bw", "b.txt"))
  # current behaviour, possibly a bug: an unknown type has no pattern, so nothing is filtered.
  expect_identical(wizard_scoped_choices(root, "nope"), c("a.bw", "b.txt"))
})

# ---- wizard_name_has_bad_dot() and wiz_clean_optional_text() -----------------------------

test_that("wizard_name_has_bad_dot(): a middle dot is bad unless one has a space on each side", {
  expect_true(wizard_name_has_bad_dot("ATAC\u00b7rep1"))
  expect_true(wizard_name_has_bad_dot("ATAC \u00b7rep1"))
  expect_true(wizard_name_has_bad_dot("ATAC\u00b7 rep1"))
  expect_true(wizard_name_has_bad_dot("\u00b7"))
  expect_false(wizard_name_has_bad_dot("ATAC \u00b7 rep1"))
  expect_false(wizard_name_has_bad_dot("ATAC.rep1"))
  expect_false(wizard_name_has_bad_dot(""))
  # current behaviour, possibly a bug: one spaced dot lets a second, bare one through.
  expect_false(wizard_name_has_bad_dot("a \u00b7 b\u00b7c"))
  # Any white space counts, a tab too.
  expect_false(wizard_name_has_bad_dot("a\t\u00b7\tb"))
})

test_that("wiz_clean_optional_text(): blank is NULL, text is trimmed", {
  expect_null(wiz_clean_optional_text(NULL))
  expect_null(wiz_clean_optional_text(""))
  expect_null(wiz_clean_optional_text("   "))
  expect_null(wiz_clean_optional_text(" \t\n "))
  expect_null(wiz_clean_optional_text(NA))
  expect_identical(wiz_clean_optional_text("  org.Mm.eg.db "), "org.Mm.eg.db")
  expect_identical(wiz_clean_optional_text("NA"), "NA")
  expect_error(wiz_clean_optional_text(character(0)), "argument is of length zero")
})

test_that(".pt_require_genome_label(): stops on a blank label, naming the settings file", {
  expect_error(.pt_require_genome_label("", "/p/project_settings.json"),
               paste0(PT_GENOME_LABEL_REQUIRED, ' Set "genome_label" in /p/project_settings.json, ',
                      "or fill it in again with plottrellis_wizard()."), fixed = TRUE)
  expect_error(.pt_require_genome_label(NULL, "x.json"), "Genome label is required", fixed = TRUE)
  expect_error(.pt_require_genome_label("  ", "x.json"), "Genome label is required", fixed = TRUE)
  expect_null(.pt_require_genome_label("m10", "x.json"))  # an unknown one is not checked here
  expect_invisible(.pt_require_genome_label("mm10", "x.json"))
})

# ---- wiz_save_problems() -----------------------------------------------------------------

test_that("wiz_save_problems(): spaces-only fields count as blank; a file is not a data folder", {
  one <- list(list(condition = "A"))
  expect_identical(wiz_save_problems("  ", one, "g.gtf", "org.Mm.eg.db", "mm10"), "Data root is required.")
  f <- tempfile(); file.create(f); on.exit(unlink(f), add = TRUE)
  expect_identical(wiz_save_problems(f, one, "g.gtf", "org.Mm.eg.db", "mm10"), "Data root directory does not exist.")
  expect_identical(wiz_save_problems(tempdir(), one, " \t", "org.Mm.eg.db", "mm10"), "GTF path is required.")
  expect_identical(wiz_save_problems(tempdir(), one, "g.gtf", "   ", "mm10"), "OrgDb package is required.")
  expect_identical(wiz_save_problems(tempdir(), NULL, "g.gtf", "org.Mm.eg.db", "mm10"), "Add at least one sample before saving.")
})

test_that("wiz_save_problems(): the GTF path is only required, not checked; labels are trimmed; NCBI names pass", {
  one <- list(list(condition = "A"))
  expect_identical(wiz_save_problems(tempdir(), one, "/no/such.gtf", "org.Mm.eg.db", " mm10 "), character(0))
  expect_identical(wiz_save_problems(tempdir(), one, "g.gtf", "org.Mm.eg.db", "GRCm39"), character(0))
  expect_identical(wiz_save_problems(tempdir(), one, "g.gtf", "org.Mm.eg.db", "MM10"),
                   "Genome label 'MM10' isn't a known genome build id (e.g. mm10, mm39, hg38).")
})

test_that("wiz_save_problems(): every kind of problem at once, in the order Save lists them", {
  expect_identical(
    wiz_save_problems(file.path(tempdir(), "no-such-dir"), list(), "", "org.Xx.eg.db", "m10"),
    c("Data root directory does not exist.", "Add at least one sample before saving.", "GTF path is required.",
      "OrgDb package 'org.Xx.eg.db' is not installed. Install it first, e.g. BiocManager::install('org.Xx.eg.db').",
      "Genome label 'm10' isn't a known genome build id (e.g. mm10, mm39, hg38)."))
})

# ---- wiz_gtf_check() ---------------------------------------------------------------------

test_that("wiz_gtf_check(): a gzipped GTF and an upper-case extension are read", {
  gz <- hp_gtf(hp_gene("ENSMUSG00000000001.4"), ext = ".gtf.gz")
  chk <- wiz_gtf_check(gz)
  expect_identical(chk, list(state = "ok", text = paste0("Looks like a GTF with Ensembl gene IDs (mouse). ",
                                                         "Make sure it is the release your reads were quantified with."),
                             species = MOUSE, build = NULL))
  up <- hp_gtf(hp_gene("ENSG00000141510"), ext = ".GTF.GZ")
  expect_identical(wiz_gtf_check(up)$species, HUMAN)
  # A path with spaces around it is trimmed first.
  expect_identical(wiz_gtf_check(paste0("  ", gz, " "))$state, "ok")
})

test_that("wiz_gtf_check(): a corrupt .gz reads as no lines, quietly", {
  bad <- tempfile(fileext = ".gtf.gz"); on.exit(unlink(bad), add = TRUE)
  writeBin(as.raw(c(0x1f, 0x8b, 0x08, 0x00, 0xff, 0xfe, 0x00, 0x01)), bad)
  expect_silent(chk <- wiz_gtf_check(bad))
  expect_identical(chk$state, "warn")
  expect_match(chk$text, "No gene_id attribute in the first 2000 lines", fixed = TRUE)
})

test_that("wiz_gtf_check(): header lines count towards n_lines, and the message names n_lines", {
  f <- hp_gtf(c(rep("#!genome-build GRCm39", 3), hp_gene("ENSMUSG00000000001")))
  expect_identical(wiz_gtf_check(f, n_lines = 3)$text, "No gene_id attribute in the first 3 lines. Is this a GTF file?")
  chk <- wiz_gtf_check(f, n_lines = 4)
  expect_identical(chk$state, "ok")
  expect_identical(chk$build, "GRCm39")
})

test_that("wiz_gtf_check(): gene_id after a tab or a semicolon counts, glued to another word does not", {
  tab <- hp_gtf('chr1\tsrc\tgene\t1\t100\t.\t+\t.\tgene_name "G";\tgene_id "ENSMUSG00000000001";')
  expect_identical(wiz_gtf_check(tab)$state, "ok")
  semi <- hp_gtf('chr1\tsrc\tgene\t1\t100\t.\t+\t.\tgene_name "G";gene_id "ENSMUSG00000000001";')
  expect_identical(wiz_gtf_check(semi)$state, "ok")
  # The attribute column after its tab is the first attribute: matched by the tab.
  first <- hp_gtf(hp_gene("ENSMUSG00000000001"))
  expect_identical(wiz_gtf_check(first)$state, "ok")
  glued <- hp_gtf('chr1\tsrc\tgene\t1\t100\t.\t+\t.\txgene_id "ENSMUSG00000000001";')
  expect_match(wiz_gtf_check(glued)$text, "No gene_id attribute", fixed = TRUE)
  # Two spaces before the quote are not the GTF form.
  two <- hp_gtf(hp_gene("ENSMUSG00000000001", attr_sep = "  "))
  expect_match(wiz_gtf_check(two)$text, "No gene_id attribute", fixed = TRUE)
})

test_that("wiz_gtf_check(): an empty gene_id is shown as the non-Ensembl example", {
  f <- hp_gtf(hp_gene(""))
  chk <- wiz_gtf_check(f)
  expect_identical(chk$state, "warn")
  expect_match(chk$text, "The gene IDs here look like '', not Ensembl IDs", fixed = TRUE)
})

test_that("wiz_gtf_check(): a GENCODE provider line with no description names the species only", {
  f <- hp_gtf(c("##provider: GENCODE", "##format: gtf", hp_gene("ENSMUSG00000000001")))
  chk <- wiz_gtf_check(f)
  expect_identical(chk$text, "Looks like a GENCODE GTF (mouse). Make sure it is the release your reads were quantified with.")
  expect_null(chk$build)
})

test_that("wiz_gtf_check(): a GENCODE description without 'evidence-based' and a human build", {
  f <- hp_gtf(c("##description: human genome (GRCh37), version 19 (Ensembl 74)", "##provider: GENCODE",
                hp_gene("ENSG00000223972.4")))
  chk <- wiz_gtf_check(f)
  expect_identical(chk$text, paste0("Looks like a GENCODE GTF: human genome (GRCh37), version 19 (Ensembl 74). ",
                                    "Make sure it is the release your reads were quantified with."))
  expect_identical(chk$build, "GRCh37")
  # A description with no GRC build: none found.
  g <- hp_gtf(c("##description: a subset", hp_gene("ENSG00000223972")))
  expect_null(wiz_gtf_check(g)$build)
})

test_that("wiz_gtf_check(): another species' warning still reports the build", {
  f <- hp_gtf(c("#!genome-build GRCh38", hp_gene("ENSDARG00000000001")))
  chk <- wiz_gtf_check(f)
  expect_identical(chk$state, "warn")
  expect_identical(chk$text, paste0("Looks like an Ensembl GTF, genome build GRCh38, but the gene IDs (ENSDARG...) ",
                                    "are not mouse or human. PlotTrellis supports mouse and human only."))
  expect_identical(chk$build, "GRCh38")
  expect_null(chk$species)
})

test_that("wiz_gtf_check(): blank results carry NULL text, species and build", {
  expect_identical(wiz_gtf_check("  "), list(state = "blank", text = NULL, species = NULL, build = NULL))
})

# ---- wiz_orgdb_check() and wiz_genome_label_check() --------------------------------------

test_that("wiz_orgdb_check(): trimmed, NULL is blank, the human package matches a human GTF", {
  expect_identical(wiz_orgdb_check(NULL), list(state = "blank", text = NULL))
  expect_identical(wiz_orgdb_check("   ", MOUSE), list(state = "blank", text = NULL))
  skip_if_not_installed("org.Mm.eg.db")
  expect_identical(wiz_orgdb_check(" org.Mm.eg.db "), list(state = "ok", text = "Installed."))
  skip_if_not_installed("org.Hs.eg.db")
  expect_identical(wiz_orgdb_check("org.Hs.eg.db", HUMAN), list(state = "ok", text = "Installed, matches the GTF (human)."))
  # Not installed wins over a species mismatch.
  expect_match(wiz_orgdb_check("org.Xx.eg.db", HUMAN)$text, "Not installed in this R.", fixed = TRUE)
})

test_that("wiz_genome_label_check(): each known build's label, trimmed", {
  expect_identical(wiz_genome_label_check("hg19", "GRCh37"), list(state = "ok", text = "Matches the GTF's build (GRCh37)."))
  expect_identical(wiz_genome_label_check("mm39", "GRCm39")$state, "ok")
  expect_identical(wiz_genome_label_check(" mm10 ", "GRCm38")$state, "ok")
  expect_identical(wiz_genome_label_check("hg38", "GRCh38")$state, "ok")
  expect_identical(wiz_genome_label_check("hg38", "GRCh37")$text, "The GTF is for GRCh37, whose UCSC name is hg19, not hg38.")
})

test_that("wiz_genome_label_check(): the build wins over the species; unknown labels pass", {
  # A known build: the species isn't looked at.
  expect_identical(wiz_genome_label_check("mm10", "GRCm38", HUMAN)$state, "ok")
  # An unknown build falls through to the species check.
  expect_identical(wiz_genome_label_check("hg38", "GRCz11", MOUSE)$text, "hg38 is a human build, but the GTF's gene IDs are mouse.")
  # A label that isn't one of the four: nothing to compare.
  expect_identical(wiz_genome_label_check("dm6", NULL, MOUSE), list(state = "blank", text = NULL))
  expect_identical(wiz_genome_label_check("MM10", NULL, HUMAN)$state, "blank")
})

# ---- wiz_error_list(), .pt_colour_input(), wiz_section(), wiz_gtf_guide() ----------------

test_that("wiz_error_list(): no heading, no messages, and any kind but warn is an error", {
  html <- hp_html(wiz_error_list(NULL, c("a", "b")))
  expect_no_match(html, "<b>", fixed = TRUE)
  expect_identical(hp_count(html, "<li>"), 2L)
  empty <- hp_html(wiz_error_list("Heading", character(0)))
  expect_match(empty, "<b>Heading</b>", fixed = TRUE)
  expect_identical(hp_count(empty, "<li>"), 0L)
  other <- hp_html(wiz_error_list("H", "x", kind = "info"))
  expect_match(other, '<div class="wiz-error wiz-error-box" role="alert">', fixed = TRUE)
  # Messages are escaped as text.
  expect_match(hp_html(wiz_error_list("H", "<i>&")), "<li>&lt;i&gt;&amp;</li>", fixed = TRUE)
})

test_that(".pt_colour_input(): a colour swatch reporting its value on input and on render", {
  html <- hp_html(.pt_colour_input("s_ref", "Reference colour:", "Tip", "#4D4D4D"))
  expect_match(html, '<label title="Tip">Reference colour:</label>', fixed = TRUE)
  expect_match(html, '<input type="color" id="s_ref" value="#4D4D4D"', fixed = TRUE)
  expect_match(html, "oninput=\"Shiny.setInputValue(&#39;s_ref&#39;, this.value, {priority: &#39;event&#39;});\"", fixed = TRUE)
  expect_match(html, "<script>Shiny.setInputValue('s_ref', \"#4D4D4D\");</script>", fixed = TRUE)
})

test_that("wiz_section(): numbered heading, optional intro, then the content", {
  html <- hp_html(wiz_section(2, "Data folder", "Intro text.", tags$p("body")))
  expect_match(html, '<span class="wiz-step">2</span>', fixed = TRUE)
  expect_match(html, "<h3>Data folder</h3>", fixed = TRUE)
  expect_match(html, '<p class="wiz-section-intro">Intro text.</p>', fixed = TRUE)
  expect_lt(regexpr("Intro text.", html, fixed = TRUE), regexpr("<p>body</p>", html, fixed = TRUE))
  expect_no_match(hp_html(wiz_section(1, "T", NULL)), "wiz-section-intro", fixed = TRUE)
})

test_that("wiz_gtf_guide(): a closed details block with the four genome rows and the download links", {
  html <- hp_html(wiz_gtf_guide())
  expect_match(html, '<details class="wiz-guide">', fixed = TRUE)
  expect_match(html, "<summary>Which GTF? Where do I get it?</summary>", fixed = TRUE)
  expect_no_match(html, "<details class=\"wiz-guide\" open", fixed = TRUE)
  expect_identical(hp_count(html, "<tr>"), 5L)  # a header row and four genomes
  for (g in c("Mouse, GRCm39", "Mouse, GRCm38", "Human, GRCh38", "Human, GRCh37"))
    expect_match(html, sprintf("<td>%s</td>", g), fixed = TRUE)
  for (l in c("mm39", "mm10", "hg38", "hg19")) expect_match(html, sprintf("<code>%s</code>", l), fixed = TRUE)
  expect_match(html, 'href="https://www.gencodegenes.org/mouse/" target="_blank"', fixed = TRUE)
})

# ---- Display settings: collector, defaults, reader -----------------------------------------

test_that(".pt_collect_display_settings(): the input beats the settings, which beat the default", {
  dflt <- .pt_default_settings()
  got <- .pt_collect_display_settings(list(s_font = " Arial ", s_show_title = FALSE), "s_",
                                      list(font = "Helvetica", page_w = 12, show_title = TRUE))
  expect_identical(got$font, "Arial")         # trimmed
  expect_identical(got$show_title, FALSE)     # FALSE is a value, not a blank
  expect_identical(got$page_w, 12)
  expect_identical(got$fontsize, dflt$fontsize)
  expect_identical(got$chr_fmt, "dot")
})

test_that(".pt_collect_display_settings(): a blank number falls back; a blank text does not", {
  got <- .pt_collect_display_settings(list(s_page_w = NA, s_font = "", s_hic_norm = "  "), "s_", list(page_w = 12))
  expect_identical(got$page_w, 12)
  # current behaviour, possibly a bug: the comment says a blank input falls back, but an
  # empty text box is kept as "" (only NULL and NA fall back).
  expect_identical(got$font, "")
  # W93: a blank Normalisation box falls back to the saved value, else KR.
  expect_identical(got$hic_norm, "KR")
})

test_that(".pt_collect_display_settings(): the optional numbers are blank for auto, the settings only when missing", {
  # A blank (NA) input is auto even when the settings have a value.
  expect_null(.pt_collect_display_settings(list(s_zrange_cap = NA), "s_", list(zrange_cap = 50))$zrange_cap)
  # A missing input takes the settings' value; a non-number or {} in the settings is auto.
  expect_identical(.pt_collect_display_settings(list(), "s_", list(zrange_cap = 50))$zrange_cap, 50)
  expect_null(.pt_collect_display_settings(list(), "s_", list(hic_height = list()))$hic_height)
  expect_null(.pt_collect_display_settings(list(), "s_", list(hic_resolution = Inf))$hic_resolution)
  expect_identical(.pt_collect_display_settings(list(s_hic_resolution = 5000), "s_", list())$hic_resolution, 5000)
  # The optional numbers never fall back to a default (theirs is NULL anyway).
  got <- .pt_collect_display_settings(list(), "", list())
  expect_true(all(c("zrange_cap", "hic_height", "hic_resolution") %in% names(got)))
  expect_null(got$zrange_cap)
})

test_that(".pt_collect_display_settings(): names in a fixed order, with no prefix as the wizard calls it", {
  got <- .pt_collect_display_settings(list(), "", .pt_default_settings())
  expect_identical(names(got), c(
    "font", "fontsize", "page_w", "chr_fmt", "show_track_axes", "show_title", "show_scale_bar", "hic_norm",
    "hic_style", "hic_colour_scale", "zrange_cap", "atac_height", "hic_height", "hic_res_mode", "hic_max_bins_per_view",
    "hic_resolution", "diff_colour_lo", "diff_colour_hi", "junction_opacity_cap", "min_junction_support",
    "reference_colour", "comparison_colour", "gene_flank_pct", "min_display_tpm", "uniprot_reviewed_only",
    "arc_clip_partial", "junction_height", "diff_shared_alpha"))
  d <- .pt_display_defaults()
  expect_identical(got[setdiff(names(got), c("zrange_cap", "hic_height", "hic_resolution"))],
                   d[setdiff(names(got), c("zrange_cap", "hic_height", "hic_resolution"))])
})

test_that(".pt_display_defaults() and .pt_default_settings(): keys, order, the NULL autos and the mm39 genome", {
  d <- .pt_display_defaults()
  expect_identical(names(d), c(
    "font", "fontsize", "page_w", "chr_fmt", "show_track_axes", "show_title", "show_scale_bar", "hic_norm",
    "hic_style", "hic_colour_scale", "zrange_cap", "atac_height", "hic_height", "hic_res_mode", "hic_max_bins_per_view",
    "hic_resolution", "gene_flank_pct", "min_junction_support", "min_display_tpm", "uniprot_reviewed_only",
    "arc_clip_partial", "junction_height", "diff_colour_lo", "diff_colour_hi", "junction_opacity_cap",
    "diff_shared_alpha", "reference_colour", "comparison_colour"))
  expect_identical(names(Filter(is.null, d)), c("zrange_cap", "hic_height", "hic_resolution"))
  expect_identical(d$fontsize, 8L)
  s <- .pt_default_settings()
  expect_identical(names(s), c(.PT_CORE_SETTINGS, names(d)))
  expect_identical(s[.PT_CORE_SETTINGS], list(gtf_path = "", orgdb = "org.Mm.eg.db", genome_label = "mm39",
                                              data_dir = "", samples = list()))
  expect_identical(s[names(d)], d)
})

test_that(".pt_read_settings(): no file gives the defaults; samples stay lists; unknown keys are kept", {
  defaults <- .pt_default_settings()
  expect_identical(.pt_read_settings(file.path(tempdir(), "no-such-settings.json"), defaults), defaults)
  json <- tempfile(fileext = ".json"); on.exit(unlink(json), add = TRUE)
  writeLines('{"samples": [{"condition": "A", "colour": "#000000", "data": [{"type": "bw", "name": "X", "path": "a.bw"}]},
                          {"condition": "B", "colour": "#111111", "data": []}],
              "extra_key": 3, "page_w": 7}', json)
  cfg <- .pt_read_settings(json, defaults)
  expect_identical(cfg$samples, list(
    list(condition = "A", colour = "#000000", data = list(list(type = "bw", name = "X", path = "a.bw"))),
    list(condition = "B", colour = "#111111", data = list())))
  expect_identical(cfg$extra_key, 3L)
  expect_identical(cfg$page_w, 7L)
  expect_identical(utils::tail(names(cfg), 1), "extra_key")
  # No samples key, or a null one: none.
  writeLines('{"page_w": 7}', json)
  expect_identical(.pt_read_settings(json, defaults)$samples, list())
  writeLines('{"samples": null}', json)
  expect_identical(.pt_read_settings(json, defaults)$samples, list())
})

test_that(".pt_clean_display_settings(): one bad key warns in the singular; core and unknown keys untouched", {
  cfg <- c(.pt_default_settings(), list(extra = "x"))
  cfg$page_w <- "wide"; cfg$gtf_path <- 5
  expect_warning(got <- .pt_clean_display_settings(cfg, "S.json"),
                 "^Ignoring unusable value in S.json for page_w; using the default\\.$")
  expect_identical(got$page_w, 9)
  expect_identical(got$gtf_path, 5)
  expect_identical(got$extra, "x")
})

test_that(".pt_clean_display_settings(): which values count as unusable, and which as empty", {
  d <- .pt_display_defaults()
  # Empty boxes: the default, silently, a missing key included.
  cfg <- list(font = NA, hic_norm = "NA", page_w = NULL, zrange_cap = NA)
  expect_no_warning(got <- .pt_clean_display_settings(cfg, "S.json"))
  expect_identical(got$font, d$font)
  expect_identical(got$hic_norm, d$hic_norm)
  expect_identical(got$page_w, d$page_w)
  expect_true("zrange_cap" %in% names(got))
  expect_null(got$zrange_cap)
  expect_identical(names(got)[1:5], c("font", "hic_norm", "page_w", "zrange_cap", "fontsize"))
  # Wrong types: a list for a value, text for a logical, a number for text, Inf for a number.
  bad <- list(fontsize = list(1), show_title = "TRUE", chr_fmt = 1, atac_height = Inf, hic_height = "2",
              hic_resolution = list())
  expect_warning(got <- .pt_clean_display_settings(bad, "S.json"),
                 "for fontsize, chr_fmt, show_title, atac_height, hic_height; using the defaults.", fixed = TRUE)
  expect_identical(got$fontsize, 8L)
  expect_identical(got$show_title, TRUE)
  expect_identical(got$chr_fmt, "dot")
  expect_identical(got$atac_height, 0.85)
  expect_null(got$hic_height)
  expect_identical(got$hic_resolution, list())  # {} for an auto number is kept as it is
  # An integer for a double default is fine, as is a small positive page width.
  expect_no_warning(got <- .pt_clean_display_settings(list(page_w = 1L, gene_flank_pct = 0L), "S.json"))
  expect_identical(got$page_w, 1L)
})

# ---- Writing the project files ------------------------------------------------------------

test_that(".pt_write_project_files(): the JSON keeps nulls as null and is written in place, no temp file left", {
  proj <- hp_dir()
  hp_touch(proj, "data/a.bw")
  json <- file.path(proj, "project_settings.json")
  settings <- list(data_dir = file.path(proj, "data"), samples = list(), gtf_path = "", zrange_cap = NULL,
                   page_w = 9, show_title = TRUE)
  expect_true(.pt_write_project_files(settings, json))
  raw <- readLines(json)
  expect_true(any(grepl('"zrange_cap": null', raw, fixed = TRUE)))
  expect_true(any(grepl('"samples": []', raw, fixed = TRUE)))
  expect_identical(sort(list.files(proj)), c("data", "project_samples.txt", "project_settings.json"))
  # An existing file is replaced.
  expect_true(.pt_write_project_files(utils::modifyList(settings, list(page_w = 11)), json))
  expect_identical(jsonlite::read_json(json)$page_w, 11L)
})

test_that(".pt_write_project_files(): a folder that isn't there stops with the reason", {
  json <- file.path(tempdir(), "no-such-folder-x", "project_settings.json")
  expect_error(.pt_write_project_files(list(page_w = 9), json),
               "^cannot open file '.*project_settings-.*\\.json': No such file or directory; cannot open the connection$")
})

test_that(".pt_write_project_files(): the samples copy failing gives FALSE, the settings still written", {
  proj <- hp_dir()
  json <- file.path(proj, "project_settings.json")
  dir.create(file.path(proj, "project_samples.txt"))
  expect_false(.pt_write_project_files(list(data_dir = proj, samples = list(), page_w = 9), json))
  expect_identical(jsonlite::read_json(json)$page_w, 9L)
})

test_that(".pt_write_samples_copy(): absolute folder and GTF, the genome always, other settings as deltas", {
  proj <- hp_dir()
  hp_touch(proj, c("data/atac/a.bw", "ref/g.gtf"))
  json <- file.path(proj, "project_settings.json")
  withr::with_dir(proj, {
    settings <- list(data_dir = "data", gtf_path = "ref/g.gtf", orgdb = "org.Mm.eg.db", genome_label = "mm10",
                     samples = list(list(condition = "Ba", colour = "#1B9E77",
                                         data = list(list(type = "bw", name = "ATAC", path = "atac/a.bw")))),
                     page_w = 9, fontsize = 11, zrange_cap = NULL, not_a_setting = "x")
    expect_true(.pt_write_samples_copy(settings, json))
  })
  expect_identical(readLines(file.path(proj, "project_samples.txt")), c(
    paste0("# data_root: ", normalizePath(file.path(proj, "data"))),
    paste0("# gtf_path: ", normalizePath(file.path(proj, "ref", "g.gtf"))),
    "# orgdb: org.Mm.eg.db", "# genome_label: mm10", "# fontsize: 11",
    "condition\tcolour\tbw:ATAC", "Ba\t#1B9E77\tatac/a.bw"))
  # current behaviour, possibly a bug: a relative GTF path to no file stays relative, and the
  # importer then reads it against the copy's own folder, not R's working folder.
  withr::with_dir(proj, expect_true(.pt_write_samples_copy(list(data_dir = "data", gtf_path = "ref/none.gtf", samples = list()), json)))
  expect_identical(readLines(file.path(proj, "project_samples.txt"))[2], "# gtf_path: ref/none.gtf")
})

test_that(".pt_write_samples_copy(): a blank GTF is left blank; a start view goes after the # lines", {
  proj <- hp_dir()
  hp_touch(proj, "data/atac/a.bw")
  json <- file.path(proj, "project_settings.json")
  settings <- list(data_dir = file.path(proj, "data"), gtf_path = "", orgdb = "org.Mm.eg.db",
                   samples = list(list(condition = "Lp", colour = "#1B9E77",
                                       data = list(list(type = "bw", name = "ATAC", path = "atac/a.bw")))),
                   start_view = .pt_start_view_entry("single", "Lp", "Elf5", c("genes", "bw_ATAC_Lp")))
  expect_true(.pt_write_samples_copy(settings, json))
  expect_identical(readLines(file.path(proj, "project_samples.txt")), c(
    paste0("# data_root: ", normalizePath(file.path(proj, "data"))), "# gtf_path: ", "# orgdb: org.Mm.eg.db",
    "# start_gene: Elf5", "# start_tracks: genes, ATAC", "# start_samples: NA, Lp",
    "condition\tcolour\tbw:ATAC", "Lp\t#1B9E77\tatac/a.bw"))
  # A start view that isn't a list is ignored; one whose lines fail is left out quietly.
  settings$start_view <- "Elf5"
  expect_true(.pt_write_samples_copy(settings, json))
  expect_false(any(grepl("^# start_", readLines(file.path(proj, "project_samples.txt")))))
  settings$start_view <- list(view_mode = "single", conditions = list("Gone"), gene = "Elf5", tracks = list("bw_ATAC_Gone"))
  expect_true(.pt_write_samples_copy(settings, json))
  expect_false(any(grepl("^# start_", readLines(file.path(proj, "project_samples.txt")))))
})

test_that(".pt_write_samples_copy(): a settings list it can't export gives FALSE and writes nothing", {
  proj <- hp_dir()
  json <- file.path(proj, "project_settings.json")
  expect_false(.pt_write_samples_copy(list(data_dir = NULL, samples = list()), json))
  expect_false(file.exists(file.path(proj, "project_samples.txt")))
})

# ---- .pt_kept_start_view() ---------------------------------------------------------------

test_that(".pt_kept_start_view(): the same GTF spelled another way keeps the view; a non-list view is none", {
  root <- hp_dir()
  hp_touch(root, c("data/atac/a.bw", "ref/g.gtf"))
  sv <- list(view_mode = "single", conditions = list("Ba"), gene = "Elf5", tracks = list("bw_ATAC_Ba", "genes"))
  smp <- list(list(condition = "Ba", colour = "#1B9E77", data = list(list(type = "bw", name = "ATAC", path = "atac/a.bw"))))
  old <- list(gtf_path = file.path(root, "ref", "g.gtf"), orgdb = "org.Mm.eg.db", start_view = sv)
  new <- list(gtf_path = file.path(root, "ref", "..", "ref", "g.gtf"), orgdb = "org.Mm.eg.db",
              data_dir = file.path(root, "data"), samples = smp)
  expect_identical(.pt_kept_start_view(old, new), sv)
  # Another sample added, the file path changed: kept.
  more <- new
  more$samples <- c(smp, list(list(condition = "Lp", colour = "#000000", data = list())))
  more$samples[[1]]$data[[1]]$path <- "elsewhere.bw"
  expect_identical(.pt_kept_start_view(old, more), sv)
  expect_null(.pt_kept_start_view(utils::modifyList(old, list(start_view = "Elf5")), new))
  # A blank GTF on both sides is the same GTF.
  expect_identical(.pt_kept_start_view(utils::modifyList(old, list(gtf_path = "")), utils::modifyList(new, list(gtf_path = ""))), sv)
  # The view's sample renamed: none.
  ren <- new; ren$samples[[1]]$condition <- "Bb"
  expect_null(.pt_kept_start_view(old, ren))
})

# ---- The form and the page ---------------------------------------------------------------

WIZ_FORM_IDS <- c("wiz_genome_preset", "wiz_gtf_path", "wiz_gtf_status", "wiz_orgdb", "wiz_orgdb_status",
                  "wiz_genome_label", "wiz_genome_label_status", "wiz_data_root", "wiz_data_root_status",
                  "wiz_add_sample", "wiz_samples_ui", "wiz_completion_error", "wiz_save", "wiz_export_samples")

test_that("wiz_shared_form(): every input and output id once; Close only when editing", {
  for (is_edit in c(FALSE, TRUE)) {
    html <- hp_html(wiz_shared_form(.pt_default_settings(), is_edit))
    for (id in WIZ_FORM_IDS) expect_identical(hp_count(html, sprintf('id="%s"', id)), 1L, info = paste(id, is_edit))
    expect_identical(hp_count(html, 'id="wiz_close"'), as.integer(is_edit))
    for (k in setdiff(names(.pt_default_settings()), .PT_CORE_SETTINGS))
      expect_identical(hp_count(html, sprintf('id="wiz_%s"', k)), 0L, info = k)
  }
})

test_that("wiz_shared_form(): values from the settings; the data folder only when editing", {
  cfg <- utils::modifyList(.pt_default_settings(), list(gtf_path = "/g.gtf", orgdb = "org.Hs.eg.db",
                                                        genome_label = "hg38", data_dir = "/d"))
  fresh <- hp_html(wiz_shared_form(cfg, FALSE))
  edit  <- hp_html(wiz_shared_form(cfg, TRUE))
  for (html in list(fresh, edit)) {
    expect_match(html, 'id="wiz_gtf_path" type="text" class="shiny-input-text form-control" value="/g.gtf"', fixed = TRUE)
    expect_match(html, 'id="wiz_orgdb" type="text" class="shiny-input-text form-control" value="org.Hs.eg.db"', fixed = TRUE)
    expect_match(html, 'id="wiz_genome_label" type="text" class="shiny-input-text form-control" value="hg38"', fixed = TRUE)
  }
  expect_match(fresh, 'id="wiz_data_root" type="text" class="shiny-input-text form-control" value=""', fixed = TRUE)
  expect_match(edit, 'id="wiz_data_root" type="text" class="shiny-input-text form-control" value="/d"', fixed = TRUE)
  expect_match(fresh, "Save and open the app", fixed = TRUE)
  expect_match(edit, "Save and return to the app", fixed = TRUE)
  # The presets, after a blank choice.
  expect_match(fresh, '<option value="" selected>Choose, or fill the fields below yourself</option>', fixed = TRUE)
  for (p in names(GENOME_PRESETS)) expect_match(fresh, sprintf('<option value="%s">%s</option>', p, GENOME_PRESETS[[p]]$label), fixed = TRUE)
  # Unset genome fields are blank.
  blank <- hp_html(wiz_shared_form(list(), FALSE))
  expect_match(blank, 'id="wiz_gtf_path" type="text" class="shiny-input-text form-control" value=""', fixed = TRUE)
})

test_that("wizard_ui(): its title and heading say setup or settings by whether the file exists", {
  json <- tempfile(fileext = ".json"); on.exit(unlink(json), add = TRUE)
  page <- function() { r <- htmltools::renderTags(wizard_ui(json)); paste(r$head, r$html) }
  fresh <- page()
  expect_match(fresh, "<title>PlotTrellis: set up your project</title>", fixed = TRUE)
  expect_match(fresh, "<h1>Welcome to PlotTrellis</h1>", fixed = TRUE)
  expect_identical(hp_count(fresh, 'id="wiz_body"'), 1L)
  expect_match(fresh, "pt_once_ready", fixed = TRUE)
  # The form itself is rendered by the server, not in the page.
  expect_no_match(fresh, 'id="wiz_gtf_path"', fixed = TRUE)
  writeLines("{}", json)
  edit <- page()
  expect_match(edit, "<title>PlotTrellis: project settings</title>", fixed = TRUE)
  expect_match(edit, "<h1>Project settings</h1>", fixed = TRUE)
})

# ---- The helpers inside wizard_server() --------------------------------------------------

hp_mocks <- function(env = parent.frame()) {
  rec <- new.env()
  rec$modals <- character(0); rec$later <- list()
  local_mocked_bindings(showModal = function(ui, ...) rec$modals <- c(rec$modals, hp_html(ui)), .env = env)
  local_mocked_bindings(later = function(func, delay = 0, ...) {
    rec$later[[length(rec$later) + 1]] <- list(func = func, delay = delay)
    invisible(NULL)
  }, .package = "later", .env = env)
  rec
}
hp_project <- function(env = parent.frame()) {
  proj <- hp_dir(env)
  hp_touch(proj, c("data/atac/a.bw", "data/quant/S1/quant.sf"))
  proj
}

test_that(".pt_wizard_status_line(): nothing for blank, a marked line per state, red for any other state", {
  proj <- hp_project()
  shiny::testServer(hp_server(proj), {
    sl <- .pt_wizard_status_line
    expect_null(sl(list(state = "blank", text = "x")))
    expect_identical(hp_html(sl(list(state = "ok", text = "Fine."))), '<div class="wiz-status wiz-ok">\u2713 Fine.</div>')
    expect_identical(hp_html(sl(list(state = "warn", text = "Hm."))), '<div class="wiz-status wiz-warn">\u26a0 Hm.</div>')
    expect_identical(hp_html(sl(list(state = "missing", text = "No."))), '<div class="wiz-status wiz-error">\u2717 No.</div>')
    expect_identical(hp_html(sl(list(state = "other", text = NULL))), '<div class="wiz-status wiz-error">\u2717 </div>')
  })
})

test_that("add_sample_modal(): a fresh dialog from the defaults; an error one from the inputs", {
  proj <- hp_project()
  shiny::testServer(hp_server(proj), {
    session$setInputs(wiz_data_root = file.path(proj, "data"), wiz_sample_name = "Typed", wiz_quantdir_path = "typed/q")
    fresh <- hp_html(h("add_sample_modal")())
    expect_match(fresh, '<h4 class="modal-title">Add a sample</h4>', fixed = TRUE)
    expect_match(fresh, 'id="wiz_sample_name" type="text" class="shiny-input-text form-control" value=""', fixed = TRUE)
    expect_match(fresh, '<option value="" selected></option>\n<option value="quant/S1">quant/S1</option></select>', fixed = TRUE)
    expect_no_match(fresh, "wiz-error-box", fixed = TRUE)
    err <- hp_html(h("add_sample_modal")(error = c("One.", "Two.")))
    expect_match(err, "<li>One.</li>", fixed = TRUE)
    expect_match(err, 'value="Typed"', fixed = TRUE)
    # A typed salmon path not under the data folder is added to the choices, selected.
    expect_match(err, '<option value="typed/q" selected>typed/q</option>', fixed = TRUE)
    # An edit of a name that is no longer there reads as an Edit with no sample's values.
    modal_edit_name("Gone")
    gone <- hp_html(h("add_sample_modal")())
    expect_match(gone, '<h4 class="modal-title">Edit sample</h4>', fixed = TRUE)
    expect_match(gone, 'id="wiz_sample_name" type="text" class="shiny-input-text form-control" value=""', fixed = TRUE)
  })
})

test_that("wiz_prefill_rows(): first-seen type:name pairs, no salmon rows, paths kept only on request", {
  proj <- hp_project()
  shiny::testServer(hp_server(proj), {
    samples <- list(
      list(condition = "A", data = list(list(type = "bw", name = "ATAC", path = "a.bw"),
                                        list(type = "quantdir", name = "Q", path = "q"),
                                        list(type = "bed", name = "ATAC", path = "p.bed"))),
      list(condition = "B", data = list(list(type = "bw", name = "ATAC", path = "b.bw"),
                                        list(type = "hic", name = "HiC", path = "m.hic"))))
    expect_identical(modal_row_id(), 0)
    rows <- h("wiz_prefill_rows")(samples, keep_path = FALSE)
    expect_identical(rows, list(list(id = 1, type = "bw", name = "ATAC", path = ""),
                                list(id = 2, type = "bed", name = "ATAC", path = ""),
                                list(id = 3, type = "hic", name = "HiC", path = "")))
    expect_identical(modal_row_id(), 3)
    kept <- h("wiz_prefill_rows")(samples[2], keep_path = TRUE)
    expect_identical(kept, list(list(id = 4, type = "bw", name = "ATAC", path = "b.bw"),
                                list(id = 5, type = "hic", name = "HiC", path = "m.hic")))
    expect_identical(h("wiz_prefill_rows")(list(), keep_path = TRUE), list())
    expect_identical(modal_row_id(), 5)
  })
})

test_that("add_modal_row(): appends a bare row with the next id", {
  proj <- hp_project()
  shiny::testServer(hp_server(proj), {
    modal_instances(list(list(id = 7, type = "bw", name = "X", path = "")))
    modal_row_id(7)
    h("add_modal_row")("paired-end")
    expect_identical(modal_instances(), list(list(id = 7, type = "bw", name = "X", path = ""),
                                             list(id = 8, type = "paired-end")))
    expect_identical(modal_row_id(), 8)
  })
})

test_that("wiz_display_settings() and wiz_collect_flat_settings(): the saved settings, imported ones on top, genome from the inputs", {
  proj <- hp_project()
  jsonlite::write_json(list(data_dir = file.path(proj, "data"), page_w = 12, font = "Arial", samples = list()),
                       file.path(proj, "project_settings.json"), auto_unbox = TRUE)
  shiny::testServer(hp_server(proj), {
    ds <- h("wiz_display_settings")()
    expect_identical(ds$page_w, 12L)
    expect_identical(ds$font, "Arial")
    import_visual(list(page_w = 5, zrange_cap = 40))
    ds <- h("wiz_display_settings")()
    expect_identical(ds$page_w, 5)
    expect_identical(ds$zrange_cap, 40)
    expect_identical(ds$font, "Arial")
    session$setInputs(wiz_gtf_path = " /g.gtf ", wiz_orgdb = "", wiz_genome_label = "mm10")
    flat <- h("wiz_collect_flat_settings")()
    expect_identical(names(flat), c("gtf_path", "orgdb", "genome_label", names(ds)))
    expect_identical(flat$gtf_path, "/g.gtf")
    expect_true("orgdb" %in% names(flat))
    expect_null(flat$orgdb)
    expect_identical(flat$genome_label, "mm10")
    expect_identical(flat[names(ds)], ds)
  })
})

test_that(".pt_wizard_do_close(): the closing note, a message, and a stop with restart scheduled", {
  proj <- hp_project()
  rec <- hp_mocks()
  shiny::testServer(hp_server(proj), {
    expect_message(.pt_wizard_do_close(), "^PlotTrellis: settings unchanged, returning to the app.\n$")
    expect_match(rec$modals[[1]], '<h4 class="modal-title">Closing</h4>', fixed = TRUE)
    expect_match(rec$modals[[1]], "<strong>You can close this tab now.</strong>", fixed = TRUE)
    expect_length(rec$later, 1)
    expect_identical(rec$later[[1]]$delay, 0.2)
    expect_identical(deparse(body(rec$later[[1]]$func)), 'stopApp(returnValue = "restart")')
  })
})

test_that("modal_edit_index(): NULL when adding, the row when editing, NA once removed", {
  proj <- hp_project()
  shiny::testServer(hp_server(proj), {
    wiz_samples(list(list(condition = "A", data = list()), list(condition = "B", data = list())))
    expect_null(h("modal_edit_index")())
    modal_edit_name("B")
    expect_identical(h("modal_edit_index")(), 2L)
    wiz_samples(wiz_samples()[2])
    expect_identical(h("modal_edit_index")(), 1L)
    modal_edit_name("A")
    expect_identical(h("modal_edit_index")(), NA_integer_)
  })
})
