require(data.table)
require(DESeq2)
require(GenomicRanges)

if (!exists("data.loaded"))
  stop('source("src/R/functions_balancer_genes.R") first.')

if (!exists("SV.loaded"))
{
  DHS_dt <- fread("/g/furlong/project/37_Capture-C/analysis/balancer_cap2/annotations/DNase_HS_sites_stages9-11_HotSpot_peaks_FDR-1pc_liftedToDm6.bed", header = F, sep = "\t", col.names = c("chrom", "start", "end"))
  DHS_dt <- DHS_dt[grepl("^chr[23]", chrom), ]
  DHS_dt[, start := start + 1L]
  DHS_gr <- GRanges(DHS_dt)
  DHS_gr$dhs_id <- seq_along(DHS_gr)

  DHS_James_dt <- fread("/g/furlong/project/53_DHS_variation/data/Charles_IDR/DHS-clusters_smooth80_refined50_300bp_min2peaks_expanded150_annotated.txt", header = T, sep = "\t")
  setnames(DHS_James_dt, "chr", "chrom")
  DHS_James_dt <- DHS_James_dt[grepl("^chr[23]", chrom), ]
  DHS_James_dt[, start := start + 1L]
  DHS_James_gr <- GRanges(DHS_James_dt)

  MEI_dt <- fread(paste0(Sascha.dir, "/analyses/MEI_detection/own2/final.predicted_MEI.bed"), col.names = c("chrom", "start", "end", "type"))
  MEI_dt <- MEI_dt[grepl("^chr[23]", chrom), ]
  MEI_dt[, start := start + 1L]
  MEI_dt[, type := ifelse(type == "wild type", "vrg", ifelse(type == "balancer", "bal", type))]
  MEI_bal_gr <- GRanges(MEI_dt[type == "bal", ])
  MEI_vrg_gr <- GRanges(MEI_dt[type == "vrg", ])

  MEI_gr <- reduce(c(MEI_bal_gr, MEI_vrg_gr))
  MEI_gr$type <- factor(sub("^ ", "", paste0(
    ifelse(overlapsAny(MEI_gr, MEI_bal_gr), " bal", ""),
    ifelse(overlapsAny(MEI_gr, MEI_vrg_gr), " vrg", "")
  )))

  SNV_dt <- fread("analysis/balancer/SNVs.tab", header = T, sep = "\t")
  SNV_dt <- SNV_dt[grepl("^chr[23]", chrom), ]
  SNV_dt[, start := start + 1L]
  SNV_dt[, type := factor(ifelse(VRG == BAL, "common", ifelse(REF == VRG, "bal", ifelse(REF == BAL, "vrg", "unknown"))))]
  SNV_gr <- GRanges(SNV_dt)

  SNV_common_dt <- fread("analysis/balancer/SNVs_common.tab", header = T, sep = "\t")
  SNV_common_dt <- SNV_common_dt[grepl("^chr[23]", chrom), ]
  SNV_common_dt[, start := start + 1L]
  SNV_common_dt[, type := factor(ifelse(VRG == BAL, "common", ifelse(REF == VRG, "bal", ifelse(REF == BAL, "vrg", "unknown"))))]
  SNV_common_gr <- GRanges(SNV_common_dt)

  DEL_dt <- fread(evalq(FN.bed.sv.del, env.dm6), col.names = c("chrom", "start", "end", "name"))
  DEL_dt <- DEL_dt[grepl("^chr[23]", chrom), ]
  DEL_dt[, start := start + 1L]
  DEL_dt[, genotype := sapply(strsplit(name, "_DEL"), "[", 1)]
  DEL_dt[, type := ifelse(genotype == "0/1_0/0", "bal", ifelse(genotype == "0/1_1/1", "vrg", ifelse(genotype == "1/1_1/1", "common", "unknown")))]
  DEL_bal_gr <- GRanges(DEL_dt[type == "bal", ])
  DEL_vrg_gr <- GRanges(DEL_dt[type == "vrg", ])

  DUP_dt <- fread(evalq(FN.bed.sv.dup, env.dm6), col.names = c("chrom", "start", "end", "name"))
  DUP_dt <- DUP_dt[grepl("^chr[23]", chrom), ]
  DUP_dt[, start := start + 1L]
  DUP_dt[, genotype := sapply(strsplit(name, "_DUP"), "[", 1)]
  DUP_dt[, type := ifelse(genotype == "0/1_0/0", "bal", ifelse(genotype == "0/1_0/1", "vrg", "unknown"))]
  DUP_bal_gr <- GRanges(DUP_dt[type == "bal", ])
  DUP_vrg_gr <- GRanges(DUP_dt[type == "vrg", ])

  DUP_manual_dt <- fread("/g/korbel/shared/projects/drosophila_balancer/analyses/SV_filter/MANUAL/manual_calls.txt", col.names = c("chrom", "start", "end", "name", "comment"))
  DUP_manual_dt[, comment := NULL]
  DUP_manual_dt <- DUP_manual_dt[grepl("^chr[23]", chrom), ]
  DUP_manual_dt[, start := start + 1L]
  DUP_manual_dt[, genotype := sapply(strsplit(name, "_DUP"), "[", 1)]
  DUP_manual_dt[, type := ifelse(genotype == "0/1_0/0", "DUP_manual_bal", ifelse(genotype == "0/1_0/1", "DUP_manual_vrg", "DUP_manual_unknown"))]
  stopifnot(DUP_manual_dt$type == "DUP_manual_bal")
  DUP_manual_bal_gr <- GRanges(DUP_manual_dt[type == "DUP_manual_bal", ])

  c1 <- c(DEL_bal_gr, DEL_vrg_gr)
  c1$type <- paste0("DEL_", c1$type)
  c2 <- c(DUP_bal_gr, DUP_vrg_gr)
  c2$type <- paste0("DUP_", c2$type)
  CNV_unreduced_gr <- c(c1, c2)
  rm(c1, c2)

  # non-common deletions and duplications, reduced so that they do not overlap
  CNV_gr <- reduce(c(DEL_bal_gr, DEL_vrg_gr, DUP_bal_gr, DUP_vrg_gr))
  CNV_gr$type <- factor(sub("^ ", "", paste0(
    ifelse(overlapsAny(CNV_gr, DEL_bal_gr), " DEL_bal", ""),
    ifelse(overlapsAny(CNV_gr, DEL_vrg_gr), " DEL_vrg", ""),
    ifelse(overlapsAny(CNV_gr, DUP_bal_gr), " DUP_bal", ""),
    ifelse(overlapsAny(CNV_gr, DUP_vrg_gr), " DUP_vrg", "")
  )))

  # CAD4_dt <- fread("/g/korbel/shared/projects/drosophila_balancer/analyses/tracks/enhancer/CAD4_plus_vienna_dm6.core80percent.onlyChr2and3.bed",
  #   header = F, sep = "\t", col.names = c("chrom", "start", "end"))
  # CAD4_dt <- CAD4_dt[grepl("^chr[23]", chrom), ]
  # CAD4_dt[, start := start + 1L]
  # CAD4_gr <- GRanges(CAD4_dt)

  CAD4_dt <- fread("/g/furlong/project/CAD/CAD4_corrected_names/CAD4_plus_vienna_minus_inactive_corrected_names_dm6.bed",
    header = F, sep = "\t", col.names = c("chrom", "start", "end", "name", "score", "strand"))
  CAD4_dt <- CAD4_dt[grepl("^chr[23]", chrom), ]
  CAD4_dt[, start := start + 1L]
  CAD4_gr <- GRanges(CAD4_dt)

  # 8008 CRMs, Zinzen et al. 2009
  mesoCRM_gr <- rtracklayer::import("/g/furlong/project/37_Capture-C/analysis/balancer_cap2/annotations/MesoCRM_dm6_Nature_Zinzen2009.gff")

  SV.loaded <- T
}

if (!exists("diffHiC.loaded"))
{
  prefix <- "dm6_HiC_DB_6-8h_combined_BAL_vs_VRG_not_across_breakpoint_dist_100kb_filtered_5000_corrected"
  thresh <- 0.1

  bins <- NULL
  load(paste0("data/hicexplorer/rda/", prefix, ".rda")) # bins
  dds <- NULL
  res <- NULL
  load(paste0("data/DESeq2/", prefix, ".rda")) # dds, res

  diffHiC_dt1 <- as.data.table(res)
  diffHiC_dt1[, bin := elementMetadata(dds)$bin1]
  diffHiC_dt1[, bait_bin := elementMetadata(dds)$bin2]
  diffHiC_dt1 <- diffHiC_dt1[!is.na(padj) & padj < thresh]
  diffHiC_dt1[, diffc_id := seq_len(nrow(diffHiC_dt1))]

  diffHiC_dt2 <- as.data.table(res)
  diffHiC_dt2[, bin := elementMetadata(dds)$bin2]
  diffHiC_dt2[, bait_bin := elementMetadata(dds)$bin1]
  diffHiC_dt2 <- diffHiC_dt2[!is.na(padj) & padj < thresh]
  diffHiC_dt2[, diffc_id := seq_len(nrow(diffHiC_dt2))]
  diffHiC_dt2 <- diffHiC_dt2[bin != bait_bin, ] # these were included already in diffHiC_dt1

  diffHiC_dt <- rbind(diffHiC_dt1, diffHiC_dt2)
  rm(diffHiC_dt1, diffHiC_dt2)
  diffHiC_dt <- cbind(diffHiC_dt, bins[diffHiC_dt$bin])
  diffHiC_dt[, baitChr := bins$chrom[bait_bin]]
  diffHiC_dt[, baitStart := bins$start[bait_bin]]
  diffHiC_dt[, baitEnd := bins$end[bait_bin]]
  setkey(diffHiC_dt, chrom, start)

  diffHiC_dt[, DHS_overlap := countOverlaps(with(diffHiC_dt, GRanges(chrom, IRanges(start, end))), DHS_gr) > 0]
  diffHiC_dt[, DHS_James_overlap := countOverlaps(with(diffHiC_dt, GRanges(chrom, IRanges(start, end))), DHS_James_gr) > 0]
  diffHiC_dt[, bait_DHS_overlap := countOverlaps(with(diffHiC_dt, GRanges(baitChr, IRanges(baitStart, baitEnd))), DHS_gr) > 0]

  diffHiC_dt[, CAD4_overlap := countOverlaps(with(diffHiC_dt, GRanges(chrom, IRanges(start, end))), CAD4_gr) > 0]
  diffHiC_dt[, bait_CAD4_overlap := countOverlaps(with(diffHiC_dt, GRanges(baitChr, IRanges(baitStart, baitEnd))), CAD4_gr) > 0]

  annotate_overlap <- function(gr, anno_gr)
  {
    anno_gr <- reduce(anno_gr)
    ov <- findOverlaps(gr, anno_gr)
    ov_dt <- data.table(query = queryHits(ov), overlap_width = width(pintersect(gr[queryHits(ov)], anno_gr[subjectHits(ov)])))
    ov_dt <- ov_dt[, list(overlap_width = sum(overlap_width)), by = "query"]
    res <- rep(0L, length(gr))
    res[ov_dt$query] <- ov_dt$overlap_width
    return(res)
  }

  # annotate CNV overlap of the "other end" bin
  diffHiC_dt[, CNV_overlap := annotate_overlap(with(diffHiC_dt, GRanges(chrom, IRanges(start, end))), CNV_gr)]
  # annotate CNV overlap of the "bait" bin
  diffHiC_dt[, bait_CNV_overlap := annotate_overlap(with(diffHiC_dt, GRanges(baitChr, IRanges(baitStart, baitEnd))), CNV_gr)]

  diffHiC_gr <- with(diffHiC_dt, GRanges(chrom, IRanges(start, end)))
  elementMetadata(diffHiC_gr) <- diffHiC_dt[, c("diffc_id", "log2FoldChange", "padj", "DHS_overlap", "DHS_James_overlap", "CNV_overlap", "CAD4_overlap", "baitChr", "baitStart", "baitEnd", "bait_DHS_overlap", "bait_CNV_overlap", "bait_CAD4_overlap"), with = F]
  diffHiC_gr <- diffHiC_gr[diffHiC_gr$CNV_overlap < 1e3 & diffHiC_gr$bait_CNV_overlap < 1e3]

  diffHiC.loaded <- T
}

if (!exists("diffCaptureC.loaded"))
{
  pm <- NULL
  pmd <- NULL
  load("../37_Capture-C/analysis/balancer_cap2/contacts_noARS_all/DESeq2_interactions.Rdata") # baits, pmd
  load("../37_Capture-C/analysis/balancer_cap2/viewpoints.Rdata") # lib_to_gene
  diffCaptureC_pm <- pm

  # exclude the other ends affected by CNVs
  diffCaptureC_dt <- pmd[oe_affected_CNV == F, ]
  stopifnot(!diffCaptureC_dt$bait_affected_RS)
  stopifnot(!diffCaptureC_dt$bait_affected_CNV)
  stopifnot(!diffCaptureC_dt$oe_affected_RS)

  # basic annotations
  diffCaptureC_dt[, diffc_id := seq_len(nrow(diffCaptureC_dt))]
  setnames(diffCaptureC_dt, "int.log2FoldChange", "log2FoldChange")
  setnames(diffCaptureC_dt, "int.padj", "padj")
  stopifnot(with(diffCaptureC_dt, ceiling((baitStart + baitEnd) / 2) + distSign == ceiling((otherEndStart + otherEndEnd) / 2)))

  # note that 13 baits are assigned to more than one gene
  # option 1: consider the differential contacts of these baits multiple times
  diffCaptureC_dt <- merge(diffCaptureC_dt, with(lib_to_gene, data.table(baitID = frag_id, gene_id)), by = "baitID", all.x = T)
  genes_diffCaptureC <- genes[gene_id %in% lib_to_gene$gene_id, ]
  genes_diffCaptureC <- merge(genes_diffCaptureC, as.data.table(lib_to_gene)[, c("gene_id", "frag_id")], by = "gene_id", all.x = T)
  # option 2: consider them once, and do not associate them with any gene
  # diffCaptureC_dt$gene_id <- baits$gene_id[match(diffCaptureC_dt$baitID, baits$baitID)]
  # genes_diffCaptureC <- genes[gene_id %in% baits$gene_id, ]

  # exclude the genes associated with 7 baits that overlap SVs
  genes_diffCaptureC <- genes_diffCaptureC[frag_id %in% baits$baitID]

  # even more annotations
  diffCaptureC_dt[, DHS_overlap := countOverlaps(with(diffCaptureC_dt, GRanges(otherEndChr, IRanges(otherEndStart, otherEndEnd))), DHS_gr) > 0]
  diffCaptureC_dt[, DHS_James_overlap := countOverlaps(with(diffCaptureC_dt, GRanges(otherEndChr, IRanges(otherEndStart, otherEndEnd))), DHS_James_gr) > 0]
  diffCaptureC_dt[, bait_DHS_overlap := countOverlaps(with(diffCaptureC_dt, GRanges(baitChr, IRanges(baitStart, baitEnd))), DHS_gr) > 0]

  diffCaptureC_dt[, CAD4_overlap := countOverlaps(with(diffCaptureC_dt, GRanges(otherEndChr, IRanges(otherEndStart, otherEndEnd))), CAD4_gr) > 0]
  diffCaptureC_dt[, bait_CAD4_overlap := countOverlaps(with(diffCaptureC_dt, GRanges(baitChr, IRanges(baitStart, baitEnd))), CAD4_gr) > 0]

  # convert to GenomicRanges
  diffCaptureC_gr <- with(diffCaptureC_dt, GRanges(otherEndChr, IRanges(otherEndStart, otherEndEnd)))
  elementMetadata(diffCaptureC_gr) <- diffCaptureC_dt[, c("diffc_id", "baitID", "gene_id", "log2FoldChange", "padj", "DHS_overlap", "DHS_James_overlap", "CAD4_overlap", "baitChr", "baitStart", "baitEnd", "baitName", "distSign", "bait_DHS_overlap", "bait_CAD4_overlap"), with = F]

  diffCaptureC.loaded <- T
}

if (!exists("diffIS.loaded"))
{
  source("src/R/functions_IS.R")
  load("analysis/balancer/differential_IS.Rdata") # lm.fit, b

  IS_bal_dt <- read_IS_all_columns("HiC_DB_6-8h_combined_BAL", source = "filtered_5000_delta_0.1", genome = "dm6bal3")
  IS_bal_dt <- rbind(IS_bal_dt, IS_bal_dt[, list(value = mean(value), group = NA), by = .(chrom, start, end, midpoint)])
  IS_bal_dt[, allele := "b"]

  IS_vrg_dt <- read_IS_all_columns("HiC_DB_6-8h_combined_VRG", source = "filtered_5000_delta_0.1")
  IS_vrg_dt <- rbind(IS_vrg_dt, IS_vrg_dt[, list(value = mean(value), group = NA), by = .(chrom, start, end, midpoint)])
  IS_vrg_dt[, allele := "w"]
  IS_vrg_dt$value <- predict(lm.fit, data.table(IS.VRG = IS_vrg_dt$value))

  IS_dt <- rbind(IS_bal_dt, IS_vrg_dt)

  diffIS_dt <- b
  dt <- diffIS_dt[diff.IS.class %in% c("lower_5percent", "upper_5percent"), ]
  diffIS_gr <- with(dt, GRanges(chrom, IRanges(start + 1L, end)))
  diffIS_gr$log2FoldChange <- dt$diff.IS
  diffIS_bal_gr <- diffIS_gr[dt$diff.IS.class == "lower_5percent"] # more insulated in balancer
  diffIS_vrg_gr <- diffIS_gr[dt$diff.IS.class == "upper_5percent"] # less insulated in balancer

  diffIS.loaded <- T
}

diffTADs_dt <- fread("analysis/balancer/differential_TADs.tab", header = T)
diffTADs_common_gr <- GRanges(diffTADs_dt[class %in% c("matched", "shifted"), ])
diffTADs_bal_gr <- GRanges(diffTADs_dt[class == "balancer-specific", ])
diffTADs_vrg_gr <- GRanges(diffTADs_dt[class == "wild-type-specific", ])

# DHS peaks >=5% deleted
subset_overlap_percent <- function(query, subject, percentage)
{
  subject <- reduce(subject)
  ov <- findOverlaps(query, subject)
  ov_dt <- data.table(query = queryHits(ov), query_width = width(query[queryHits(ov)]),
    overlap_width = width(pintersect(query[queryHits(ov)], subject[subjectHits(ov)])))
  ov_dt <- ov_dt[, list(overlap_width = sum(overlap_width)), by = c("query", "query_width")]
  sel <- ov_dt[overlap_width / query_width >= percentage]$query
  return(query[sel])
}
DHS_deleted_vrg_gr <- subset_overlap_percent(DHS_gr, DEL_vrg_gr, 0.05)
DHS_deleted_bal_gr <- subset_overlap_percent(DHS_gr, DEL_bal_gr, 0.05)
DHS_deleted_dt <- rbind(data.table(as.data.table(DHS_deleted_vrg_gr), type = "vrg"), data.table(as.data.table(DHS_deleted_bal_gr), type = "bal"))
setnames(DHS_deleted_dt, "seqnames", "chrom")

if (!exists("CORE.loaded"))
{
  # CORE_dt <- fread("data/CORE/CORE_RNAseq_Nechaev_dm6_aliases_mapped.tsv", header = T, sep = "\t")
  CORE_dt <- fread("data/CORE/CORE_CAGE_Hoskins_dm6_aliases_mapped.tsv", header = T, sep = "\t")
  CORE_gr <- GRanges(CORE_dt)
  CORE.loaded <- T
}
