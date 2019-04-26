options <- commandArgs(trailingOnly = TRUE)

if (length(options) < 2)
  stop("Usage:  Rscript process_DESeq2_balancer_cap2.R designDir outputDir")

designDir <- options[1]
outputDir <- options[2]

message("designDir: ", designDir)
message("outputDir: ", outputDir)

require(Chicago)
require(GenomicRanges)
options(warn = 1)

#
#  read CHiCAGO data
#

cd.VRG <- readRDS(paste0(outputDir, "/dm6_DBal_BAL_N1_4-8h_CAP2_Rep1Rep2_VRG.Rds"))
cd.BAL <- readRDS(paste0(outputDir, "/dm6_DBal_BAL_N1_4-8h_CAP2_Rep1Rep2_BAL.Rds"))
cd.All <- readRDS(paste0(outputDir, "/dm6_DBal_BAL_N1_4-8h_CAP2_Rep1Rep2_All.Rds"))

#
#  prepare count data matrix for DESeq2
#

cols <- c("baitID", "otherEndID", "N.1", "N.2", "s_j", "s_i")
pm_all <- cd.All@x[, c(cols, "distSign", "log.p", "score"), with = F]

pm_all <- merge(pm_all, cd.VRG@x[, cols, with = F], by = c("baitID", "otherEndID"), all.x = T, suffixes = c("", ".VRG"))
pm_all$N.1.VRG[is.na(pm_all$N.1.VRG)] <- 0L
pm_all$N.2.VRG[is.na(pm_all$N.2.VRG)] <- 0L
pm_all$s_j.VRG[is.na(pm_all$s_j.VRG)] <- 1
pm_all$s_i.VRG[is.na(pm_all$s_i.VRG)] <- 1

pm_all <- merge(pm_all, cd.BAL@x[, cols, with = F], by = c("baitID", "otherEndID"), all.x = T, suffixes = c("", ".BAL"))
pm_all$N.1.BAL[is.na(pm_all$N.1.BAL)] <- 0L
pm_all$N.2.BAL[is.na(pm_all$N.2.BAL)] <- 0L
pm_all$s_j.BAL[is.na(pm_all$s_j.BAL)] <- 1
pm_all$s_i.BAL[is.na(pm_all$s_i.BAL)] <- 1

#
#  if needed, downsample VRG data to match BAL
#

if (grepl("VRGdownBAL", designDir))
{
  set.seed(42)

  ind <- rep(1:nrow(pm_all), pm_all$N.1.VRG)
  ds <- rle(sort(sample(ind, sum(pm_all$N.1.BAL))))
  pm_all$N.1.VRG <- 0L
  pm_all$N.1.VRG[ds$values] <- ds$lengths
  stopifnot(sum(pm_all$N.1.VRG) == sum(pm_all$N.1.BAL))

  ind <- rep(1:nrow(pm_all), pm_all$N.2.VRG)
  ds <- rle(sort(sample(ind, sum(pm_all$N.2.BAL))))
  pm_all$N.2.VRG <- 0L
  pm_all$N.2.VRG[ds$values] <- ds$lengths
  stopifnot(sum(pm_all$N.2.VRG) == sum(pm_all$N.2.BAL))
}

#
#  add the coordinates of bait and other end
#

baitmap <- fread(paste0(designDir, "/dm6_DpnII.baitmap"))
names(baitmap) <- c("chr", "start", "end", "ID", "name")

pm_all <- merge(pm_all, baitmap, by.x = "baitID", by.y = "ID", all.x = T)
setnames(pm_all, "chr", "baitChr")
setnames(pm_all, "start", "baitStart")
setnames(pm_all, "end", "baitEnd")
setnames(pm_all, "name", "baitName")

oe <- fread(paste0(designDir, "/dm6_DpnII.rmap"))
names(oe) <- c("chr", "start", "end", "ID")

pm_all <- merge(pm_all, oe, by.x = "otherEndID", by.y = "ID", all.x = T)
setnames(pm_all, "chr", "otherEndChr")
setnames(pm_all, "start", "otherEndStart")
setnames(pm_all, "end", "otherEndEnd")

#
#  save all interactions
#

save(pm_all, file = paste0(outputDir, "/all_interactions.Rdata"))

#
#  limit to the pre-selected set of restriction fragment pairs (according to distance in reference and balancer assembly)
#

load(paste0(designDir, "/dm6_dm6bal3_DpnII_fragment_distance.Rdata")) # dist_fit

pm <- merge(pm_all, dist_fit, by = c("baitID", "otherEndID"))

message(length(unique(pm_all$baitID)), " baits within all interactions")
message(length(unique(pm$baitID)), " baits within the interactions tested by DESeq2")
rm(pm_all)

#
#  now do the real DESeq2 job
#

require(DESeq2)

library("BiocParallel")
register(MulticoreParam(16))

colData <- data.frame(allele = factor(NULL, c("VRG", "BAL")), replicate = factor())
colData <- rbind(colData, data.frame(allele = "VRG", replicate = c("R1", "R2")))
colData <- rbind(colData, data.frame(allele = "BAL", replicate = c("R1", "R2")))

# order of the columns must match the order of rows in colData
countData <- pm[, c("N.1.VRG", "N.2.VRG", "N.1.BAL", "N.2.BAL"), with = F]
countData <- as.matrix(countData)

dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ allele)
dds <- DESeq(dds, parallel = T)
res <- results(dds, alpha = 0.05, parallel = T)
print(summary(res))

pdf(paste0(outputDir, "/DESeq2_diagnostics.pdf"))
plotDispEsts(dds)
plotMA(dds)
dev.off()

pm$int.log2FoldChange <- res$log2FoldChange
pm$int.padj <- res$padj

#
#  aggregate the DESeq2 results with captured fragment annotations
#

load("analysis/balancer_cap2/viewpoints.Rdata") # lib
baits <- as.data.table(lib)
baits$name <- names(lib)
setnames(baits, "seqnames", "chrom")
setnames(baits, "frag_id", "baitID")

# take only the baits which were used in the design files
stopifnot(baitmap$ID %in% baits$baitID)
baits <- baits[baitID %in% baitmap$ID, ]
baits_gr <- GRanges(baits$chrom, IRanges(baits$start, baits$end + 4L))

cols <- c("baitID", "is.control", "breakpoint_distance", "embryo.log2FoldChange", "embryo.padj")
pm <- merge(pm, baits[, cols, with = F], by = "baitID", all.x = T)

#
#  check if the interactions are in the same breakpoint interval
#

pm1_gr <- GRanges(pm$baitChr, IRanges(pm$baitStart, pm$baitEnd + 4L))
pm2_gr <- GRanges(pm$otherEndChr, IRanges(pm$otherEndStart, pm$otherEndEnd + 4L))
ov_filler <- data.table(queryHits = 1:nrow(pm), subjectHits = 0L)

balancer.dir <- "/g/furlong/project/39_Balancer_Hi-C"
load(paste0(balancer.dir, "/analysis/breakpoints_dm6_int_gr.Rdata")) # int_gr

ov1 <- findOverlaps(pm1_gr, int_gr)
ov2 <- findOverlaps(pm2_gr, int_gr)
ov_dt <- merge(merge(ov_filler, ov1, all = T, by = "queryHits", suffixes = c("", ".1")), ov2, all = T, by = "queryHits", suffixes = c("", ".2"))
same_int <- ov_dt$queryHits[ov_dt$subjectHits.1 == ov_dt$subjectHits.2]
pm$across_breakpoint <- !(seq_len(nrow(pm)) %in% same_int)

#
#  read all deletions and duplications
#

sv_gr <- GRanges()

read_bed <- function(filename)
{
  bed <- fread(filename, header = F, sep = "\t")
  bed <- bed[, 1:3, with = F]
  names(bed) <- c("chrom", "start", "end")
  bed_gr <- GRanges(bed$chrom, IRanges(bed$start + 1L, bed$end))
  return(bed_gr)
}

bed <- read_bed("analysis/balancer_cap2/annotations/DEL.final.BAL.bed")
bed$BAL <- 0L
bed$VRG <- 1L
sv_gr <- c(sv_gr, bed)

bed <- read_bed("analysis/balancer_cap2/annotations/DEL.final.VRG.bed")
bed$BAL <- 1L
bed$VRG <- 0L
sv_gr <- c(sv_gr, bed)

bed <- read_bed("analysis/balancer_cap2/annotations/DEL.final.common.bed")
bed$BAL <- 0L
bed$VRG <- 0L
sv_gr <- c(sv_gr, bed)

bed <- read_bed("analysis/balancer_cap2/annotations/DUP.final.BAL.bed")
bed$BAL <- 2L
bed$VRG <- 1L
sv_gr <- c(sv_gr, bed)

bed <- read_bed("analysis/balancer_cap2/annotations/DUP.final.VRG.bed")
bed$BAL <- 1L
bed$VRG <- 2L
sv_gr <- c(sv_gr, bed)

#
#  annotate the baits using deletions and duplications
#

# summarize the contributions of all SVs
SV_contribution <- function(gr)
{
  ov <- findOverlaps(gr, sv_gr)
  ov_dt <- as.data.table(ov)
  ov_dt$overlapWidth <- width(pintersect(gr[queryHits(ov)], sv_gr[subjectHits(ov)]))
  ov_dt$BAL <- sv_gr$BAL[subjectHits(ov)]
  ov_dt$VRG <- sv_gr$VRG[subjectHits(ov)]

  aggr <- ov_dt[, list(BAL = sum(overlapWidth * (BAL - 1L)), VRG = sum(overlapWidth * (VRG - 1L))), by = "queryHits"]
  return(aggr)
}

svc <- SV_contribution(baits_gr)
baits[, widthBAL := width(baits_gr)]
baits$widthBAL[svc$queryHits] <- baits$widthBAL[svc$queryHits] + svc$BAL
baits[, widthVRG := width(baits_gr)]
baits$widthVRG[svc$queryHits] <- baits$widthVRG[svc$queryHits] + svc$VRG

#
#  annotate the interactions using deletions and duplications
#

svc <- SV_contribution(pm1_gr)
pm[, baitWidthBAL := width(pm1_gr)]
pm$baitWidthBAL[svc$queryHits] <- pm$baitWidthBAL[svc$queryHits] + svc$BAL
pm[, baitWidthVRG := width(pm1_gr)]
pm$baitWidthVRG[svc$queryHits] <- pm$baitWidthVRG[svc$queryHits] + svc$VRG

svc <- SV_contribution(pm2_gr)
pm[, oeWidthBAL := width(pm2_gr)]
pm$oeWidthBAL[svc$queryHits] <- pm$oeWidthBAL[svc$queryHits] + svc$BAL
pm[, oeWidthVRG := width(pm2_gr)]
pm$oeWidthVRG[svc$queryHits] <- pm$oeWidthVRG[svc$queryHits] + svc$VRG

#
#  annotate whether restriction fragments are affected by DEL/DUP
#

annotated_DpnII <- fread("analysis/digest_DpnII_dm6_balancer_allele_specific_variation.tab", header = T)
baits$affected_RS <- annotated_DpnII$affected_RS[baits$baitID]
baits$affected_CNV <- annotated_DpnII$affected_CNV[baits$baitID]
pm$bait_affected_RS <- annotated_DpnII$affected_RS[pm$baitID]
pm$bait_affected_CNV <- annotated_DpnII$affected_CNV[pm$baitID]
pm$oe_affected_RS <- annotated_DpnII$affected_RS[pm$otherEndID]
pm$oe_affected_CNV <- annotated_DpnII$affected_CNV[pm$otherEndID]

#
#  overlap other ends with mesoderm CRMs, CAD4 and DHSes
#

meso_gr <- rtracklayer::import("analysis/balancer_cap2/annotations/MesoCRM_dm6_Nature_Zinzen2009.gff")

ov2 <- findOverlaps(pm2_gr, meso_gr)
pm$oeMesoCRM <- seq_len(nrow(pm)) %in% queryHits(ov2)

cad_positives <- fread("analysis/balancer_cap2/annotations/CAD4_plus_vienna_minus_inactive_corrected_names_dm6.bed", header = F, sep = "\t")
cad_positives <- cad_positives[, 1:3, with = F]
names(cad_positives) <- c("chrom", "start", "end")
cad_positives_gr <- GRanges(cad_positives$chrom, IRanges(cad_positives$start + 1L, cad_positives$end))

ov2 <- findOverlaps(pm2_gr, cad_positives_gr)
pm$oeCAD4 <- seq_len(nrow(pm)) %in% queryHits(ov2)

dhs <- fread("analysis/balancer_cap2/annotations/DNase_HS_sites_stages9-11_HotSpot_peaks_FDR-1pc_liftedToDm6.bed", header = F, sep = "\t")
dhs <- dhs[, 1:3, with = F]
names(dhs) <- c("chrom", "start", "end")
dhs_gr <- GRanges(dhs$chrom, IRanges(dhs$start + 1L, dhs$end))

ov2 <- findOverlaps(pm2_gr, dhs_gr)
pm$oeDHS <- seq_len(nrow(pm)) %in% queryHits(ov2)

#
#  save the differential interactions
#

pmd <- pm[int.padj < 0.05 & abs(int.log2FoldChange) > log2(1.5), ]
print(table(sign(pmd$int.log2FoldChange)))
save(baits, pm, pmd, file = paste0(outputDir, "/DESeq2_interactions.Rdata"))
