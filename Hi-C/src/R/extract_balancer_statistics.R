options(warn = 1)

source("src/R/functions_balancer_genes.R")
source("src/R/functions_balancer_annotations.R")
source("src/R/logFC_vs_dist.R")

#
#  Differential Hi-C and Capture-C contacts
#

{
  prefix <- "dm6_HiC_DB_6-8h_combined_BAL_vs_VRG_not_across_breakpoint_dist_100kb_filtered_5000_corrected"
  thresh <- 0.1
  dds <- NULL
  res <- NULL
  load(paste0("data/DESeq2/", prefix, ".rda")) # dds, res

  message(nrow(res), " Hi-C contacts within at most 100 kb")
  message(sum(is.finite(res$padj)), " Hi-C contacts tested (contacts with low counts weren't tested)")
  message(length(unique(diffHiC_dt$diffc_id)), " differential Hi-C contacts")
  message(length(unique(diffHiC_gr$diffc_id)), " differential Hi-C contacts after filtering out the 5 kb bins with >= 1 kb overlap of CNVs")
  message(length(unique(diffHiC_gr$diffc_id[abs(diffHiC_gr$log2FoldChange) > log2(1.5)])), " of them have |fold change|>1.5")
  message("")
}

{
  pm <- NULL
  pmd <- NULL
  load("../37_Capture-C/analysis/balancer_cap2/contacts_noARS_all/DESeq2_interactions.Rdata") # baits, pm, pmd
  load("../37_Capture-C/analysis/balancer_cap2/viewpoints.Rdata") # lib_to_gene

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

  message(length(unique(lib_to_gene$frag_id)), " Capture-C baits designed")

  # note that 13 baits are assigned to more than one gene
  # option 1: consider the differential contacts of these baits multiple times
  diffCaptureC_dt <- merge(diffCaptureC_dt, with(lib_to_gene, data.table(baitID = frag_id, gene_id)), by = "baitID", all.x = T)
  genes_diffCaptureC <- genes[gene_id %in% lib_to_gene$gene_id, ]
  genes_diffCaptureC <- merge(genes_diffCaptureC, as.data.table(lib_to_gene)[, c("gene_id", "frag_id")], by = "gene_id", all.x = T)
  # option 2: consider them once, and do not associate them with any gene
  # diffCaptureC_dt$gene_id <- baits$gene_id[match(diffCaptureC_dt$baitID, baits$baitID)]
  # genes_diffCaptureC <- genes[gene_id %in% baits$gene_id, ]

  stopifnot(!anyDuplicated(genes_diffCaptureC$gene_id))
  cat(length(unique(genes_diffCaptureC$frag_id)), " baits assigned to ", nrow(genes_diffCaptureC), " genes, namely:", sep = "")
  print(table(genes_diffCaptureC$signf))

  # exclude the genes associated with 7 baits that overlap SVs
  genes_diffCaptureC <- genes_diffCaptureC[frag_id %in% baits$baitID]

  message(nrow(baits), " Capture-C baits after excluding DpnII fragments with allele-specific variation")
  cat(length(unique(genes_diffCaptureC$frag_id)), " baits assigned to ", nrow(genes_diffCaptureC), " genes, namely:", sep = "")
  print(table(genes_diffCaptureC$signf))
}

{
  message("\n", nrow(pm), " Capture-C contacts within at most 100 kb")
  message(sum(is.finite(pm$int.padj)), " Capture-C contacts tested (contacts with low counts weren't tested)")
  message(nrow(pmd), " differential Capture-C contacts")
  message(length(unique(diffCaptureC_dt$diffc_id)), " differential Capture-C contacts after filtering out the other ends affected by CNVs")
  message(length(unique(diffCaptureC_dt$diffc_id[abs(diffCaptureC_dt$log2FoldChange) > log2(1.5)])), " of them have |fold change|>1.5")

  gr <- diffCaptureC_gr + 1e3
  gr_reduced_dt <- as.data.table(gr)[, as.data.table(reduce(GRanges(.SD))), by = "baitID"]
  message(nrow(gr_reduced_dt), " reduced (extended other ends by +/- 1 kb and merged overlapping ones) differential Capture-C contacts after filtering")
}
