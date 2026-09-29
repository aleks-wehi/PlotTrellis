# Always pools UNFILTERED (min_junction_support = 0) -- filtering is display
# logic, applied downstream in each settings_<viz>/render_<viz>.
load_junctions <- function(path, region_params) {
  pool_junctions(list(list(junction_path = path)), region_params$chrom, region_params$start, region_params$end, min_junction_support = 0)
}
