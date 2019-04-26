args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3)
{
  message("Usage:  Rscript process_balancer_log2FoldChange.R <dataset> {dm6|dm6bal3} <bin size>")
  q()
}

dataset <- args[1]
genome <- args[2]
bin_size <- as.integer(args[3])

# dataset <- "HiC_DB_6-8h"
# genome <- "dm6"
# bin_size <- as.integer(5e3)

source("src/R/functions_HiCExplorer.R")

sum_maps <- function(rda_dataset, datasets, genome, bin_size)
{
  sm <- NULL

  for (d in datasets)
  {
    fname <- paste0("data/hicexplorer/txt/", genome, "_", d, "_filtered_", bin_size, ".txt.gz")
    map <- fread_map(fname)$map

    if (is.null(sm))
      sm <- map
    else
      sm <- sm + map
  }

  # transform to coverage-normalized map
  sm <- sm / sum(sm, na.rm = T)

  return(sm)
}

process_log2FoldChange <- function(dataset, rda_dataset1, datasets1, rda_dataset2, datasets2, genome, bin_size)
{
  map1 <- sum_maps(rda_dataset1, datasets1, genome, bin_size)
  map2 <- sum_maps(rda_dataset2, datasets2, genome, bin_size)

  fin <- paste0("data/hicexplorer/txt/", genome, "_", rda_dataset1, "_filtered_", bin_size, "_corrected.txt.gz")
  fm <- fread_map(fin)
  chrom_map <- fm$chrom_map
  bins <- fm$bins
  params <- list(genome = genome, bin_size = bin_size, bins = nrow(bins))
  rm(fm)

  message("calculating ratio")
  map <- list(genome = list(log2FoldChange = log2(map1 / map2)))

  for (i in seq_len(nrow(chrom_map)))
  {
    chrom <- as.character(chrom_map$chrom[i])
    int <- chrom_map$map_start[i]:chrom_map$map_end[i]
    map[[chrom]] <- list(log2FoldChange = map$genome$log2FoldChange[int, int])
  }

  fout <- paste0("scratch/data/hicexplorer/rda/", genome, "_", dataset, "_filtered_", bin_size, "_corrected.rda")
  message("Saving ", fout)
  save(params, chrom_map, bins, map, file = fout)
  message("saved\n")
}

process_log2FoldChange(paste0(dataset, "_combined_BAL_vs_VRG_naive"),
  paste0(dataset, "_combined_BAL"), paste0(dataset, c("_R1_BAL", "_R2_BAL")),
  paste0(dataset, "_combined_VRG"), paste0(dataset, c("_R1_VRG", "_R2_VRG")),
  genome, bin_size)
# process_log2FoldChange("HiC_DB_6-8h_BAL_vs_VRG_naive", "HiC_DB_6-8h_BAL", "HiC_DB_6-8h_VRG", genome, bin_size)
# process_log2FoldChange("HiC_DB_6-8h_BAL_vs_VRG_oneexact_naive", "HiC_DB_6-8h_BAL_oneexact", "HiC_DB_6-8h_VRG_oneexact", genome, bin_size)
# process_log2FoldChange("HiC_DB_6-8h_BAL_vs_VRG_twoexact_naive", "HiC_DB_6-8h_BAL_twoexact", "HiC_DB_6-8h_VRG_twoexact", genome, bin_size)
# process_log2FoldChange("HiC_DB_6-8h_BAL_vs_VRG_5prim_naive", "HiC_DB_6-8h_BAL_5prim", "HiC_DB_6-8h_VRG_5prim", genome, bin_size)
# process_log2FoldChange("HiC_DB_6-8h_BAL_vs_VRG_oneexact_5prim_naive", "HiC_DB_6-8h_BAL_oneexact_5prim", "HiC_DB_6-8h_VRG_oneexact_5prim", genome, bin_size)
# process_log2FoldChange("HiC_DB_6-8h_BAL_vs_VRG_twoexact_5prim_naive", "HiC_DB_6-8h_BAL_twoexact_5prim", "HiC_DB_6-8h_VRG_twoexact_5prim", genome, bin_size)
