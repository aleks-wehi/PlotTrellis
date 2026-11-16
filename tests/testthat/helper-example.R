# A copy of a bundled example in a new temporary folder, for tests that change its data
# (pt_example() and plottrellis_example() give the installed files, which must not change).
# Returns the copy's samples.txt.
example_copy <- function(dataset = c("mammary", "human")) {
  src <- dirname(pt_example(match.arg(dataset)))
  dest <- tempfile("pt-example-copy-")
  dir.create(dest)
  stopifnot(all(file.copy(list.files(src, full.names = TRUE), dest, recursive = TRUE, copy.mode = FALSE)))
  file.path(dest, "samples.txt")
}

# A copied example is outside the package, so its chromosome sizes are not the shipped ones
# and ensure_genome_ready() would fetch them from UCSC: serves the shipped ones instead.
local_offline_chrom_sizes <- function(env = parent.frame()) {
  shipped <- c(mm10 = file.path(dirname(pt_example("mammary")), "gencode.vM25.example_chrom_sizes.rds"),
               hg38 = file.path(dirname(pt_example("human")), "annotation_chrom_sizes.rds"))
  local_mocked_bindings(build_chrom_sizes = function(genome_label, out_path) {
    stopifnot(genome_label %in% names(shipped))
    file.copy(shipped[[genome_label]], out_path)
  }, .env = env)
}
