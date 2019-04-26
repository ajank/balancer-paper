options <- commandArgs(trailingOnly = TRUE)

if (length(options) < 1)
  stop("Usage:  Rscript CHiCAGO_makeDesignFiles_balancer_cap2.R designDir")

designDir <- options[1]

genome <- "dm6"
settings <- list()

if (grepl("_all", designDir))
{
  message("Using option: all")
  settings <- modifyList(settings, list(minFragLen = 0L, maxFragLen = 2000000000L, minNPerBait = 0L, removeAdjacent = FALSE, tlb.filterTopPercent = 0))
}

if (grepl("_min40bp", designDir))
{
  message("Using option: min40bp")
  settings <- modifyList(settings, list(minFragLen = 44L)) # including 2 bp of each 4 bp cut site, so effectively 40 bp between cut sites
}

require(data.table)
require(GenomicRanges)

frag <- fread(paste0("analysis/digest_DpnII_", genome, "_balancer_allele_specific_variation.tab"), sep = "\t", header = T)

if (grepl("_chr23", designDir))
{
  message("Using option: chr23 (taking only chr2+3 DpnII fragments)")
  frag <- frag[grepl("^chr[23][LR]", chrom), ]
}

if (grepl("_noARS", designDir))
{
  message("Using option: noARS (excluding DpnII fragments with allele-specific Restriction Sites)")
  frag <- frag[affected_RS == F, ]
}

if (grepl("_noASV", designDir))
{
  message("Using option: noASV (excluding DpnII fragments with allele-specific variation)")
  frag <- frag[affected_CNV == F & affected_RS == F, ]
}

if (grepl("_incr_res", designDir))
{
  message("Using option: incr_res (increase resolution by merging neighboring DpnII fragments at larger distances from the viewpoint)")
  load("analysis/balancer_cap2/contacts_all/design/dm6_dm6bal3_DpnII_fragment_distance.Rdata") # diagonal_dist, dist_fit

  init_dist <- NA
  step_dist <- NA
  limit <- 4
  if (grepl("_init10kb", designDir))
    init_dist <- 10e3
  if (grepl("_step5kb", designDir))
    step_dist <- 5e3
  if (grepl("_step10kb", designDir))
    step_dist <- 10e3
  if (grepl("_step15kb", designDir))
    step_dist <- 15e3

  message("  init_dist: ", init_dist)
  message("  step_dist: ", step_dist)
  message("  limit: ", limit)
  stopifnot(1 + log2(pmax(diagonal_dist - init_dist, 0) / step_dist) > limit)

  dist_fit[, incr_res_level := pmax(floor(1 + log2(pmax(pmin(distance, distance_other) - init_dist, 0) / step_dist)), 0)]
  dist_frag <- dist_fit[, list(incr_res_level = min(incr_res_level)), by = "otherEndID"]

  frag <- merge(frag, dist_frag, all.x = T, by.x = "frag_id", by.y = "otherEndID")
  frag$incr_res_level[is.na(frag$incr_res_level) | frag$incr_res_level > limit] <- limit
  frag[, incr_res_level := as.integer(incr_res_level)]

  chrom <- ""
  first <- -1L
  count <- 1L
  min_incr_res_level <- 0L
  frag_id <- frag$frag_id # important optimization for time
  
  for (i in seq_len(nrow(frag)))
  {
    min_incr_res_level <- min(frag$incr_res_level[i], min_incr_res_level)
    if (first + count == frag_id[i] && count < 2^min_incr_res_level && chrom == frag$chrom[i])
    {
      frag_id[i] <- first
      count <- count + 1L
    }
    else
    {
      chrom <- frag$chrom[i]
      first <- frag_id[i]
      count <- 1L
      min_incr_res_level <- frag$incr_res_level[i]
    }
  }

  frag$frag_id <- frag_id

  frag <- frag[, list(start = min(start), end = max(end)), by = c("chrom", "frag_id")]
}

load("analysis/balancer_cap2/viewpoints.Rdata")

if (grepl("_away_from_breakpoint", designDir))
{
  message("Using option: away_from_breakpoint (taking only baits > 1.5 Mb from a breakpoint)")
  lib <- subset(lib, breakpoint_distance > 1.5e6)
}


suppressWarnings(dir.create(designDir))

rmap <- frag[, c("chrom", "start", "end", "frag_id"), with = F]
write.table(rmap, file = paste0(designDir, "/", genome, "_DpnII.rmap"), col.names = F, row.names = F, sep = "\t", quote = F)

baitmap <- rmap[frag_id %in% lib$frag_id, ]
baitmap$name <- names(lib)[match(baitmap$frag_id, lib$frag_id)]
write.table(baitmap, file = paste0(designDir, "/", genome, "_DpnII.baitmap"), col.names = F, row.names = F, sep = "\t", quote = F)

write.table(cbind(names(settings), settings), file = paste0(designDir, "/", genome, "_DpnII.settingsFile"), col.names = F, row.names = F, sep = "\t", quote = F)
unlink(paste0(designDir, "/", genome, "_DpnII.baitmap_4col.txt"))

pyparams <- c("minFragLen", "maxFragLen", "maxLBrownEst", "binsize", "removeb2b", "removeAdjacent")
pysettings <- settings[names(settings) %in% pyparams]
pysettingscmd <- paste0(sapply(seq_along(pysettings), function(i) paste0(" --", names(pysettings)[i], "=", pysettings[i])), collapse = "")
system(paste0("python2.7 /g/furlong/jankowsk/chicago/chicagoTools/makeDesignFiles.py", pysettingscmd, " --designDir=", designDir))
