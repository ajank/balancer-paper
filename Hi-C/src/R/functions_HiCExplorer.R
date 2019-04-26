require(data.table)

fread_map <- function(fname)
{
  genome <- strsplit(tail(strsplit(fname, "/")[[1]], 1), "_")[[1]][1]
  bin_size <- as.integer(tail(strsplit(sub("_corrected$", "", sub(".txt.gz$", "", fname)), "_")[[1]], 1))

  message("reading ", fname, " (genome ", genome, ", bin_size ", bin_size, ")")
  t <- fread(paste("pigz -dc -p 8", fname), header = T, showProgress = F, skip = 1, drop = 1)
  cat("finished reading map\n")
  t <- as.matrix(t)
  stopifnot(nrow(t) == ncol(t))
  gc()

  # note that in case of HiCExplorer contact maps, entries at diagonal are counted only once
  # in case of a not normalized contact map, (sum(t) + sum(diag(t))) / 2 is the total number of contacts

  bins <- as.data.table(do.call(rbind, strsplit(colnames(t), "\\|")))
  setnames(bins, c("V1", "V2", "V3"), c("bin_id", "separator", "coords"))
  bins[, bin_id := as.integer(bin_id) + 1L]
  stopifnot(bins$bin_id == seq_len(nrow(bins)))

  stopifnot(bins$separator == "--")
  bins[, separator := NULL]

  bins <- cbind(bins, as.data.table(do.call(rbind, strsplit(bins$coords, "[:-]"))))
  bins[, coords := NULL]
  setnames(bins, c("V1", "V2", "V3"), c("chrom", "start", "end"))
  bins[, start := as.integer(start) + 1L]
  bins[, end := as.integer(end)]

  # check if bins are sorted (--chromosomeOrder argument of hicExport)
  setkey(bins, chrom, start)
  if (!all(bins$bin_id == seq_len(nrow(bins))))
  {
    message("bins in .txt.gz file are not in order, reordering them")
    t <- t[bins$bin_id, ]
    t <- t[, bins$bin_id]
  }

  colnames(t) <- NULL
  bins[, bin_id := NULL]

  chrom_map <- fread(paste0("data/hic_maps/chrom_map_", genome, "_", bin_size, ".txt"), header = T, sep = "\t")
  stopifnot(unique(bins$chrom) == chrom_map$chrom)
  stopifnot(match(unique(bins$chrom), bins$chrom) == chrom_map$map_start)
  stopifnot(nrow(bins) == tail(chrom_map$map_end, 1))

  return(list(params = list(genome = genome, bin_size = bin_size, bins = nrow(bins)), chrom_map = chrom_map, bins = bins, map = t))
}
