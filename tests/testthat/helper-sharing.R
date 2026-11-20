# Helpers for test-sharing.R: one person (the sender) shares a project with another (the
# receiver), who must get the same figure. Temp folders stand in for their home folders.

# A copy of a bundled example's files in `dest` (created; may hold spaces), its samples.txt
# kept with its relative `# data_root: .` and `# gtf_path:` lines but without start lines.
# Returns the copy's samples.txt.
share_copy_example <- function(dest, dataset = "human") {
  src <- dirname(if (identical(dataset, "default")) pt_example() else pt_example(dataset))
  dir.create(dest, recursive = TRUE)
  stopifnot(all(file.copy(list.files(src, full.names = TRUE), dest, recursive = TRUE,
                          copy.mode = FALSE)))
  st <- file.path(dest, "samples.txt")
  writeLines(grep("^# start_", readLines(st), invert = TRUE, value = TRUE), st)
  st
}

# A copied example is outside the package, so its chromosome sizes would be fetched from UCSC:
# serves the ones the bundled examples ship instead, whichever genome builds they are.
local_example_chrom_sizes <- function(env = parent.frame()) {
  dirs <- dirname(c(pt_example(), pt_example("human")))
  label <- vapply(dirs, function(d)
    sub("^#\\s*genome_label:\\s*", "",
        grep("^#\\s*genome_label:", readLines(file.path(d, "samples.txt")), value = TRUE)[1]), "")
  shipped <- setNames(vapply(dirs, function(d)
    list.files(d, "chrom_sizes\\.rds$", full.names = TRUE)[1], ""), label)
  local_mocked_bindings(build_chrom_sizes = function(genome_label, out_path) {
    stopifnot(genome_label %in% names(shipped))
    file.copy(shipped[[genome_label]], out_path)
  }, .env = env)
}

share_project <- function(...) suppressWarnings(suppressMessages(pt_project(...)))

# The figure every sharing test compares: `gene`, every track type of the project, its second
# sample against its first, the reference (so the comparison tracks are drawn too). Left out: protein
# domains (fetched from UniProt, not read from the shared files: a network hiccup between two
# renders would change one) and every transcript of the region (slow, and drawn from the GTF
# as the transcripts track is).
share_figure <- function(proj, gene = "PSME2") {
  tracks <- setdiff(pt_tracks(proj)$track, c("protein_domains", "transcripts_all"))
  pt_figure(proj, gene, tracks, comparison = names(proj$colours)[2:1])
}

# `fig` saved as SVG; returns the file's md5. A failed track would put "track failed" in the
# SVG, so a figure compared here must draw every track.
share_svg <- function(fig, file = tempfile(fileext = ".svg")) {
  pt_save(fig, file)
  expect_false(any(grepl("track failed", readLines(file, warn = FALSE), fixed = TRUE)))
  unname(tools::md5sum(file))
}

# The first gene of an example's GTF the project's OrgDb knows (no gene name hard-coded, so
# a test survives a change of bundled example).
share_first_gene <- function(proj) {
  l <- readLines(gzfile(proj$cfg$gtf_path))
  l <- l[!startsWith(l, "#") & grepl("\tgene\t", l, fixed = TRUE)]
  intersect(sub('.*gene_name "([^"]+)".*', "\\1", l), proj$genes)[1]
}

# Every string held in projects and figures (paths included); a project's TxDb object left
# out.
share_strings <- function(...) {
  x <- lapply(list(...), function(o) unclass(o)[setdiff(names(o), "txdb")])
  unlist(rapply(x, identity, classes = "character", how = "unlist"), use.names = FALSE)
}

# Moves `path` out of the way (renamed) until the calling test ends, so nothing can still read
# it; put back afterwards. Open TxDbs are closed first (a file in use can't be renamed on
# Windows).
local_moved_away <- function(path, env = parent.frame()) {
  .pt_close_txdbs()
  away <- paste0(path, "-moved-away")
  stopifnot(file.rename(path, away))
  withr::defer(file.rename(away, path), envir = env)
  invisible(away)
}

# The sender: a copy of the human example in "sender home/my data" (a folder name with a
# space), its project and figure; made once per test run, on the first call.
.share <- new.env(parent = emptyenv())
share_sender <- function() {
  if (is.null(.share$sender)) {
    home <- file.path(tempfile("pt-share-"), "sender home")
    st <- normalizePath(share_copy_example(file.path(home, "my data")))
    home <- normalizePath(home)
    proj <- share_project(file.path(home, "PlotTrellis project"), st)
    fig <- share_figure(proj)
    .share$sender <- list(home = home, data = dirname(st), samples_txt = st, dir = proj$dir,
                          fig = fig, md5 = share_svg(fig))
  }
  .share$sender
}
