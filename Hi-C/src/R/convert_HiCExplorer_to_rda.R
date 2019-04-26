args <- commandArgs(trailingOnly = TRUE)
fin <- args[1]
fout <- args[2]
stopifnot(!is.na(fout))

source("src/R/functions_HiCExplorer.R")

# note that in case of HiCExplorer contact maps, entries at diagonal are counted only once
# in case of a not normalized contact map, (sum(t) + sum(diag(t))) / 2 is the total number of contacts

# Normalization works as follows:
# require(RcppCNPy)
# f <- npyLoad("data/hicexplorer/export/correction_factors.npy.gz")
# t(t(traw / f) / f) - t \approx 0

#
#  read HiCExplorer contact map, and adjust it to previous conventions
#

# determine which bins are masked
fm <- fread_map(fin)
bins_masked <- apply(fm$map, 2, sum) == 0

rm(fm)
gc()

# calculate more statistics from the raw contact map
fin_raw <- sub("_corrected.txt.gz$", ".txt.gz", fin)
fm_raw <- fread_map(fin_raw)
count <- sum(fm_raw$map)

# masking the raw contact map
fm_raw$map[bins_masked, ] <- 0L
fm_raw$map[, bins_masked] <- 0L
count_after_masking <- sum(fm_raw$map)

rm(fm_raw)
gc()

# read the normalized contact map again
fm <- fread_map(fin)
fm$params$count <- count
fm$params$count_after_masking <- count_after_masking
fm$bins$masked <- bins_masked
fm$params$bins_after_masking <- sum(!fm$bins$masked)
message("out of ", fm$params$bins, " bins, removing ", sum(fm$bins$masked), " and keeping ", fm$params$bins_after_masking)

# mask the normalized contact map (convert all-0 rows and columns to all-NA)
fm$map[fm$bins$masked, ] <- NA
fm$map[, fm$bins$masked] <- NA

# scale the normalized contact map to previous conventions
fm$params$map_divisor <- sum(fm$map, na.rm = T) / fm$params$bins_after_masking
fm$map <- fm$map / fm$params$map_divisor

#
#  reformat to previous .Rdata format
#

params <- fm$params
chrom_map <- fm$chrom_map
bins <- fm$bins
map <- list(genome = list(observed = fm$map))

bin_size <- params$bin_size
chrom_lst <- as.character(chrom_map$chrom)

for (chrom in chrom_lst)
{
  m <- match(chrom, chrom_map$chrom)
  int <- chrom_map$map_start[m]:chrom_map$map_end[m]
  map[[chrom]] <- list(observed = map$genome$observed[int, int])
}

# for smaller bin sizes, keep only "cis" contacts
if (bin_size < 2000L)
  map$genome <- list()

#
#  save the .rda file
#

cat("saving Rdata file without fit\n")
save(params, chrom_map, bins, map, file = fout)
cat("saved\n")
