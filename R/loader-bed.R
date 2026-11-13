# Reads only column 1 (cheap even on a large file) purely to get the observed
# chrom vocabulary -- the actual data read/draw happens in render_bed's own
# plotRanges() call, which is tolerant of non-standard extra BED columns in a
# way rtracklayer::import() is not. A file that is missing or can't be read stops;
# an empty file, or one with only header and comment lines, is a track with nothing in
# view (`empty`).
#
# A .gz file is read by fread() (here, and by plotRanges() when drawing) through the R.utils package.
load_bed <- function(entry, region_params) {
  .pt_check_file(entry$path)
  if (file.size(entry$path) == 0) {
    entry$data <- list(no_chrom = FALSE, empty = TRUE)
    return(entry)
  }
  raw <- tryCatch(data.table::fread(entry$path, header = FALSE, select = 1L), error = function(e) .pt_read_error(entry$path, e))
  vocab <- unique(as.character(raw[[1]]))
  vocab <- vocab[!grepl("^(#|track$|browser$)", vocab)]  # header and comment lines
  if (length(vocab) == 0) {
    entry$data <- list(no_chrom = FALSE, empty = TRUE)
    return(entry)
  }
  chrom <- resolve_chrom(region_params$chrom, vocab)
  if (is.null(chrom)) {
    entry$data <- list(no_chrom = TRUE)
    return(entry)
  }
  entry$data <- list(no_chrom = FALSE, chrom = chrom)
  entry
}
