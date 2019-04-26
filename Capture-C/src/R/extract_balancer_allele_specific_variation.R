library(data.table)
library(GenomicRanges)

Sascha.dir <- "/g/korbel/shared/projects/drosophila_balancer"

## Read allele-specific deletions and duplications

DEL <- fread(paste0(Sascha.dir, "/analyses/SV_filter/FINAL/DEL.final.bed"))
DUP <- fread(paste0(Sascha.dir, "/analyses/SV_filter/FINAL/DUP.final.bed"))
colnames(DUP) <- colnames(DEL) <- c("seqnames", "start", "end", "type")
DEL$type <- paste0("DEL_", gsub("0/1_1/1.*", "VRG-specific", gsub("0/1_0/0.*", "BAL-specific", gsub("1/1_1/1.*", "common", DEL$type))))
DUP$type <- paste0("DUP_", gsub("0/1_0/1.*", "VRG-specific", gsub("0/1_0/0.*", "BAL-specific", DUP$type)))
CNV <- rbind(DEL, DUP)
# exclude common deletions
CNV <- CNV[!grepl("common", type), ]
CNV_gr <- makeGRangesFromDataFrame(CNV)

## Overlap DpnII restriction sites with DEL/DUP

annotated_DpnII <- fread("analysis/digest_DpnII_dm6_balancer.tab", header = T)

# Integrate annotations of disrupted sequence of restriction sites (by SNPs)
annotated_DpnII[, affected_RS_by_SNP :=
  BAL_affected_RS1_by_SNP != 0 |
  BAL_affected_RS2_by_SNP != 0 |
  VRG_affected_RS1_by_SNP != 0 |
  VRG_affected_RS2_by_SNP != 0 |
  BAL_internal_RS != 0 |
  VRG_internal_RS != 0]
stopifnot(with(annotated_DpnII, !(affected_RS_by_SNP & avoided_count == 0)))

# Find restriction sites afected by DEL/DUP
aRS1 <- with(annotated_DpnII, GRanges(chrom, IRanges(start + 0, start + 3)))
ovl <- as.data.table(findOverlaps(aRS1, CNV_gr))
annotated_DpnII$affected_RS1_by_CNV <- F
annotated_DpnII$affected_RS1_by_CNV[unique(ovl$queryHits)] <- T

aRS2 <- with(annotated_DpnII, GRanges(chrom, IRanges(end + 1, end + 4)))
ovl <- as.data.table(findOverlaps(aRS2, CNV_gr))
annotated_DpnII$affected_RS2_by_CNV <- F
annotated_DpnII$affected_RS2_by_CNV[unique(ovl$queryHits)] <- T

# Aggregate all allotations of affected restriction sites
annotated_DpnII[, affected_RS_by_CNV := affected_RS1_by_CNV | affected_RS2_by_CNV]
annotated_DpnII[, affected_RS := affected_RS_by_SNP | affected_RS_by_CNV]

# Get DpnII fragments that are affected by DEL/DUP
afrag <- with(annotated_DpnII, GRanges(chrom, IRanges(start + 0, end + 4)))
ovl <- as.data.table(findOverlaps(afrag, CNV_gr))
annotated_DpnII$affected_CNV <- F
annotated_DpnII$affected_CNV[unique(ovl$queryHits)] <- T

# filtering by annotateDpnII_balancer should be stricter, since it also considers common short indels
with(annotated_DpnII, stopifnot(!affected_CNV | affected_RS2_by_CNV | excluded_count > 0))


print_stats <- function(dt)
{
  with(dt, print(table(affected_CNV, affected_RS)))
  with(dt, message("Fraction of DpnII sites with allele-specific CNV: ", sum(affected_CNV) / length(affected_CNV)))
  with(dt, message("Fraction of DpnII sites with allele-specific Restriction Site: ", sum(affected_RS) / length(affected_RS)))
  with(dt, message("Either of the above: ", sum(affected_CNV | affected_RS) / length(affected_CNV)))
}
cat("\nStatistics for all DpnII fragments:\n")
print_stats(annotated_DpnII)

cat("\nStatistics for DpnII fragments on chr[23][LR]:\n")
print_stats(annotated_DpnII[grepl("^chr[23][LR]", chrom), ])

write.table(annotated_DpnII, file = "analysis/digest_DpnII_dm6_balancer_allele_specific_variation.tab", quote = F, row.names = F, sep = "\t")

#
#  save coordinates of restriction fragments affected by CNV/with affected RS
#

DpnII_affected_RS <- annotated_DpnII[affected_RS == T, c("chrom", "start", "end"), with = F]
DpnII_affected_RS[, start := start - 1L]
write.table(DpnII_affected_RS, file = "analysis/balancer_cap2/annotations/DpnII_dm6_balancer_affected_RS.bed",
  quote = F, row.names = F, col.names = F, sep = "\t")

DpnII_affected_CNV <- annotated_DpnII[affected_CNV == T, c("chrom", "start", "end"), with = F]
DpnII_affected_CNV[, start := start - 1L]
write.table(DpnII_affected_CNV, file = "analysis/balancer_cap2/annotations/DpnII_dm6_balancer_affected_CNV.bed",
  quote = F, row.names = F, col.names = F, sep = "\t")
