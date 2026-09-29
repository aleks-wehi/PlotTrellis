load_bw_file <- function(path, region_params) {
  vocab <- tryCatch(as.character(GenomeInfoDb::seqnames(GenomeInfoDb::seqinfo(rtracklayer::BigWigFile(path)))), error = function(e) NULL)
  chrom <- resolve_chrom(region_params$chrom, vocab)
  if (is.null(chrom)) {
    return(list(gr = GRanges(), failed = TRUE, no_chrom = TRUE, chrom = NA_character_))
  }
  which_gr <- GRanges(chrom, IRanges(region_params$start, region_params$end))
  x <- tryCatch(rtracklayer::import(path, which = which_gr), error = function(e) NULL)
  list(gr = x %||% GRanges(), failed = is.null(x), no_chrom = FALSE, chrom = chrom)
}
