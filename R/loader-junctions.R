# Always pools UNFILTERED (min_junction_support = 0) -- filtering is display
# logic, applied downstream in each settings_<viz>/render_<viz>. A file that is
# missing or can't be read stops; an empty file has no junctions anywhere.
load_junctions <- function(path, region_params) {
  .pt_check_file(path)
  if (file.size(path) == 0) return(data.frame(intron_start = integer(0), intron_end = integer(0), nSupport = integer(0)))
  tryCatch(
    pool_junctions(list(list(junction_path = path)), region_params$chrom, region_params$start, region_params$end, min_junction_support = 0),
    error = function(e) .pt_read_error(path, e)
  )
}
