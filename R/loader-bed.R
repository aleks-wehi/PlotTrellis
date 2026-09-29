# Reads only column 1 (cheap even on a large file) purely to get the observed
# chrom vocabulary -- the actual data read/draw happens in render_bed's own
# plotRanges() call, which is tolerant of non-standard extra BED columns in a
# way rtracklayer::import() is not.
load_bed <- function(entry, region_params) {
  raw <- tryCatch(data.table::fread(entry$path, header = FALSE, select = 1L), error = function(e) NULL)
  if (is.null(raw)) {
    entry$data <- list(failed = TRUE, no_chrom = FALSE)
    return(entry)
  }
  vocab <- unique(as.character(raw[[1]]))
  chrom <- resolve_chrom(region_params$chrom, vocab)
  if (is.null(chrom)) {
    entry$data <- list(failed = TRUE, no_chrom = TRUE)
    return(entry)
  }
  entry$data <- list(failed = FALSE, chrom = chrom)
  entry
}
