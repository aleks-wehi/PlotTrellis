load_bw_file <- function(path, region_params) {
  .pt_check_file(path)
  # rtracklayer gives the reason a file can't be read as a warning ("<path> is not a big wig
  # file") before failing with "UCSC library operation failed": the warning's is the one shown,
  # in the track's failure message only. A file that reads fine passes its warnings on.
  warns <- character(0)
  vocab <- tryCatch(withCallingHandlers({
    info <- GenomeInfoDb::seqinfo(rtracklayer::BigWigFile(path))
    as.character(GenomeInfoDb::seqnames(info))
  }, warning = function(w) {
    warns <<- c(warns, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) .pt_read_error(path, if (length(warns)) simpleError(warns[length(warns)]) else e))
  for (w in warns) warning(w, call. = FALSE)
  chrom <- resolve_chrom(region_params$chrom, vocab)
  if (is.null(chrom)) {
    return(list(gr = GRanges(), no_chrom = TRUE, chrom = NA_character_))
  }
  which_gr <- GRanges(chrom, IRanges(region_params$start, region_params$end))
  x <- tryCatch(rtracklayer::import(path, which = which_gr), error = function(e) .pt_read_error(path, e))
  list(gr = x, no_chrom = FALSE, chrom = chrom)
}
