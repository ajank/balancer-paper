options(warn = 1)

require(cowplot)
require(ggplot2)
require(rcartocolor)
require(scales) # comma

source("src/R/functions_balancer_genes.R")
source("src/R/functions_balancer_annotations.R")
source("src/R/logFC_vs_dist.R")

theme_set(theme_cowplot(font_size = 11)) # reduce default font size
ts <- theme_get()$plot.subtitle
ts$hjust <- 0.5
theme_update(plot.subtitle = ts) # , legend.title = theme_get()$legend.text
theme_update(strip.text.x = element_text(size = 11), strip.text.y = element_text(size = 11))

options(mc.cores = 16)


diffHiC_BAL_gr <- diffHiC_gr[diffHiC_gr$log2FoldChange > 0]
diffHiC_VRG_gr <- diffHiC_gr[diffHiC_gr$log2FoldChange < 0]

genes_to_diffHiC_dt <- logFC_vs_dist(genes, diffHiC_gr, bp_gr, mode = "all", extra_filtering = "bait", max_dist = 99e6)
genes_to_diffHiC_dt$diffc_id <- diffHiC_gr$diffc_id[genes_to_diffHiC_dt$center_id]
genes_to_diffHiC_dt$HiClog2FoldChange <- diffHiC_gr$log2FoldChange[genes_to_diffHiC_dt$center_id]
genes_to_diffHiC_dt$DHS_overlap <- diffHiC_gr$DHS_overlap[genes_to_diffHiC_dt$center_id]
genes_to_diffHiC_dt$DHS_James_overlap <- diffHiC_gr$DHS_James_overlap[genes_to_diffHiC_dt$center_id]
genes_to_diffHiC_dt$CNV_overlap <- diffHiC_gr$CNV_overlap[genes_to_diffHiC_dt$center_id]
genes_to_diffHiC_dt$CAD4_overlap <- diffHiC_gr$CAD4_overlap[genes_to_diffHiC_dt$center_id]
genes_to_diffHiC_dt$bait_DHS_overlap <- diffHiC_gr$bait_DHS_overlap[genes_to_diffHiC_dt$center_id]
genes_to_diffHiC_dt$bait_CNV_overlap <- diffHiC_gr$bait_CNV_overlap[genes_to_diffHiC_dt$center_id]
genes_to_diffHiC_dt$bait_CAD4_overlap <- diffHiC_gr$bait_CAD4_overlap[genes_to_diffHiC_dt$center_id]

TSS_to_diffHiC_dt <- logFC_vs_dist(TSS, diffHiC_gr, bp_gr, mode = "all", extra_filtering = "bait", max_dist = 99e6)
TSS_to_diffHiC_dt$diffc_id <- diffHiC_gr$diffc_id[TSS_to_diffHiC_dt$center_id]
TSS_to_diffHiC_dt$HiClog2FoldChange <- diffHiC_gr$log2FoldChange[TSS_to_diffHiC_dt$center_id]
TSS_to_diffHiC_dt$DHS_overlap <- diffHiC_gr$DHS_overlap[TSS_to_diffHiC_dt$center_id]
TSS_to_diffHiC_dt$DHS_James_overlap <- diffHiC_gr$DHS_James_overlap[TSS_to_diffHiC_dt$center_id]
TSS_to_diffHiC_dt$CNV_overlap <- diffHiC_gr$CNV_overlap[TSS_to_diffHiC_dt$center_id]
TSS_to_diffHiC_dt$CAD4_overlap <- diffHiC_gr$CAD4_overlap[TSS_to_diffHiC_dt$center_id]
TSS_to_diffHiC_dt$bait_DHS_overlap <- diffHiC_gr$bait_DHS_overlap[TSS_to_diffHiC_dt$center_id]
TSS_to_diffHiC_dt$bait_CNV_overlap <- diffHiC_gr$bait_CNV_overlap[TSS_to_diffHiC_dt$center_id]
TSS_to_diffHiC_dt$bait_CAD4_overlap <- diffHiC_gr$bait_CAD4_overlap[TSS_to_diffHiC_dt$center_id]

with(TSS_to_diffHiC_dt, message("diffHiC:\n",
  length(unique(diffHiC_dt$diffc_id)), " differential contacts\n",
  length(unique(diffHiC_gr$diffc_id)), " differential contacts considered after CNV-based filtering\n",
    "  of them, ", length(unique(diffHiC_gr$diffc_id[diffHiC_gr$DHS_overlap])), " overlap DHSes\n",
    "  of them, ", length(unique(diffHiC_gr$diffc_id[diffHiC_gr$CNV_overlap])), " overlap CNVs\n",
    "  of them, ", length(unique(diffHiC_gr$diffc_id[diffHiC_gr$CAD4_overlap])), " overlap CAD4 entries\n",
  length(unique(diffc_id)), " differential contacts assigned to ", nrow(genes), " genes\n",
    "  of them, ", length(unique(diffc_id[DHS_overlap])), " overlap DHSes\n",
    "  of them, ", length(unique(diffc_id[CNV_overlap])), " overlap CNVs\n",
    "  of them, ", length(unique(diffc_id[CAD4_overlap])), " overlap CAD4 entries\n",
  length(unique(diffc_id[!is.na(signf) & signf != "n"])), " differential contacts assigned to ", sum(genes$signf != "n"), " testable genes\n",
    "  of them, ", length(unique(diffc_id[DHS_overlap & !is.na(signf) & signf != "n"])), " overlap DHSes\n",
    "  of them, ", length(unique(diffc_id[CNV_overlap & !is.na(signf) & signf != "n"])), " overlap CNVs\n",
    "  of them, ", length(unique(diffc_id[CAD4_overlap & !is.na(signf) & signf != "n"])), " overlap CAD4 entries\n",
  length(unique(diffc_id[!is.na(signf) & signf == "s"])), " differential contacts assigned to ", sum(genes$signf == "s"), " DE genes\n",
    "  of them, ", length(unique(diffc_id[DHS_overlap & !is.na(signf) & signf == "s"])), " overlap DHSes\n",
    "  of them, ", length(unique(diffc_id[CNV_overlap & !is.na(signf) & signf == "s"])), " overlap CNVs\n",
    "  of them, ", length(unique(diffc_id[CAD4_overlap & !is.na(signf) & signf == "s"])), " overlap CAD4 entries\n",
  length(unique(diffc_id[!is.na(signf) & signf == "s" & abs(log2FoldChange) > log2(1.5)])), " differential contacts assigned to ", sum(genes$signf == "s" & abs(genes$log2FoldChange) > log2(1.5)), " strongly DE genes\n",
    "  of them, ", length(unique(diffc_id[DHS_overlap & !is.na(signf) & signf == "s" & abs(log2FoldChange) > log2(1.5)])), " overlap DHSes\n",
    "  of them, ", length(unique(diffc_id[CNV_overlap & !is.na(signf) & signf == "s" & abs(log2FoldChange) > log2(1.5)])), " overlap CNVs\n",
    "  of them, ", length(unique(diffc_id[CAD4_overlap & !is.na(signf) & signf == "s" & abs(log2FoldChange) > log2(1.5)])), " overlap CAD4 entries\n"
))


genes_to_diffCaptureC_dt <- logFC_vs_dist(genes_diffCaptureC, diffCaptureC_gr, bp_gr, mode = "all", extra_filtering = "gene_id", max_dist = 99e6)
genes_to_diffCaptureC_dt$CaptureClog2FoldChange <- diffCaptureC_gr$log2FoldChange[genes_to_diffCaptureC_dt$center_id]
genes_to_diffCaptureC_dt$DHS_overlap <- diffCaptureC_gr$DHS_overlap[genes_to_diffCaptureC_dt$center_id]
genes_to_diffCaptureC_dt$DHS_James_overlap <- diffCaptureC_gr$DHS_James_overlap[genes_to_diffCaptureC_dt$center_id]
genes_to_diffCaptureC_dt$CAD4_overlap <- diffCaptureC_gr$CAD4_overlap[genes_to_diffCaptureC_dt$center_id]

TSS_to_diffCaptureC_dt <- logFC_vs_dist(genes_diffCaptureC, diffCaptureC_gr, bp_gr, mode = "all", extra_filtering = "gene_id", max_dist = 99e6)
TSS_to_diffCaptureC_dt$CaptureClog2FoldChange <- diffCaptureC_gr$log2FoldChange[TSS_to_diffCaptureC_dt$center_id]
TSS_to_diffCaptureC_dt$DHS_overlap <- diffCaptureC_gr$DHS_overlap[TSS_to_diffCaptureC_dt$center_id]
TSS_to_diffCaptureC_dt$DHS_James_overlap <- diffCaptureC_gr$DHS_James_overlap[TSS_to_diffCaptureC_dt$center_id]
TSS_to_diffCaptureC_dt$CAD4_overlap <- diffCaptureC_gr$CAD4_overlap[TSS_to_diffCaptureC_dt$center_id]

stopifnot(nrow(diffCaptureC_dt[oe_affected_CNV == F]) == length(diffCaptureC_gr))
with(TSS_to_diffCaptureC_dt, message("diffCaptureC:\n",
  nrow(diffCaptureC_dt), " differential contacts\n",
  length(diffCaptureC_gr), " differential contacts considered after CNV-based filtering\n",
    "  of them, ", sum(diffCaptureC_gr$DHS_overlap), " differential contacts overlap DHSes\n",
    "  of them, ", sum(diffCaptureC_gr$CAD4_overlap), " differential contacts overlap CAD4 entries\n",
  length(signf), " differential contacts assigned to ", nrow(genes_diffCaptureC), " genes\n",
    "  of them, ", sum(DHS_overlap), " differential contacts overlap DHSes\n",
    "  of them, ", sum(CAD4_overlap), " differential contacts overlap CAD4 entries\n",
  sum(signf != "n"), " differential contacts assigned to ", sum(genes_diffCaptureC$signf != "n"), " testable genes\n",
    "  of them, ", sum(DHS_overlap & signf != "n"), " differential contacts overlap DHSes\n",
    "  of them, ", sum(CAD4_overlap & signf != "n"), " differential contacts overlap CAD4 entries\n",
  sum(signf == "s"), " differential contacts assigned to ", sum(genes_diffCaptureC$signf == "s"), " DE genes\n",
    "  of them, ", sum(DHS_overlap & signf == "s"), " differential contacts overlap DHSes\n",
    "  of them, ", sum(CAD4_overlap & signf == "s"), " differential contacts overlap CAD4 entries\n",
  sum(signf == "s" & abs(log2FoldChange) > log2(1.5)), " differential contacts assigned to ", sum(genes_diffCaptureC$signf == "s" & abs(genes_diffCaptureC$log2FoldChange) > log2(1.5)), " strongly DE genes\n",
    "  of them, ", sum(DHS_overlap & signf == "s" & abs(log2FoldChange) > log2(1.5)), " differential contacts overlap DHSes\n",
    "  of them, ", sum(CAD4_overlap & signf == "s" & abs(log2FoldChange) > log2(1.5)), " differential contacts overlap CAD4 entries\n"
))


#
#  Permutation test: are differential Capture-C contacts enriched in DHS overlap?
#

shuffle_gr <- function(gr)
{
  if (!("distSign" %in% names(elementMetadata(gr))))
    gr$distSign <- ceiling((start(gr) + end(gr)) / 2) - ceiling((gr$baitStart + gr$baitEnd) / 2)
  # else
  #   stopifnot(ceiling((gr$baitStart + gr$baitEnd) / 2) + gr$distSign == ceiling((start(gr) + end(gr)) / 2))

  gr$shift_start <- start(gr) - ceiling((start(gr) + end(gr)) / 2)
  gr$shift_end <- end(gr) - ceiling((start(gr) + end(gr)) / 2)

  gr$distSign <- sample(gr$distSign) * (2L * rbinom(length(gr), 1, 0.5) - 1L)

  start(gr) <- 0L
  end(gr) <- end(gr) + 2 * 99e6
  start(gr) <- ceiling((gr$baitStart + gr$baitEnd) / 2) + gr$distSign + gr$shift_start
  end(gr) <- ceiling((gr$baitStart + gr$baitEnd) / 2) + gr$distSign + gr$shift_end

  # remove the elements that are now outside chromosome
  gr <- gr[start(gr) >= 1L]
  len <- chrom_map$length[match(as.character(seqnames(gr)), chrom_map$chrom)]
  gr <- gr[end(gr) <= len]
  return(gr)
}

calculate_anno_overlap <- function(gr, anno_gr, shuffle = F, remove_self = T)
{
  if (shuffle)
    gr <- shuffle_gr(gr)

  if (remove_self)
  {
    ov1 <- findOverlaps(GRanges(gr$baitChr, IRanges(gr$baitStart, gr$baitEnd)), anno_gr)
    ov1_dt <- data.table(gr_id = queryHits(ov1), anno_gr_id = subjectHits(ov1), ov1 = rep(1L, length(ov1)))

    ov2 <- findOverlaps(gr, anno_gr)
    ov2_dt <- data.table(gr_id = queryHits(ov2), anno_gr_id = subjectHits(ov2), ov1 = rep(0L, length(ov2)))

    # extract the rows present in ov2_dt but not ov1_dt
    ovd_dt <- rbind(ov1_dt, ov2_dt)[, list(ov1 = sum(ov1)), by = c("gr_id", "anno_gr_id")]
    return(length(unique(ovd_dt[ov1 == 0]$gr_id)) / length(gr))
  }
  else
    return(mean(countOverlaps(gr, anno_gr) > 0))
}

read_DHS_gr <- function(stage = "stages9-11")
{
  DHS_dt <- fread(paste0("/g/furlong/project/37_Capture-C/analysis/balancer_cap2/annotations/DNase_HS_sites_", stage, "_HotSpot_peaks_FDR-1pc_liftedToDm6.bed"), header = F, sep = "\t")
  DHS_dt <- DHS_dt[, 1:3, with = F]
  colnames(DHS_dt) <- c("chrom", "start", "end")
  DHS_dt <- DHS_dt[grepl("^chr[23]", chrom), ]
  DHS_dt[, start := start + 1L]
  DHS_gr <- GRanges(DHS_dt)
  return(DHS_gr)
}

res <- NULL

calculate_pvalue <- function(gr, DHS_gr, label, DHS_label, num_shuffles = 1000L, plot = F, plot_xlab = NULL)
{
  observed <- calculate_anno_overlap(gr, DHS_gr, shuffle = F)
  # shuffled <- sapply(1:num_shuffles, function(i) calculate_anno_overlap(gr, DHS_gr, shuffle = T))
  shuffled <- unlist(mclapply(1:num_shuffles, function(i) calculate_anno_overlap(gr, DHS_gr, shuffle = T), mc.preschedule = T))
  expected <- mean(shuffled)
  pval <- 2 * min(sum(shuffled <= observed), sum(shuffled >= observed)) / num_shuffles
  res <<- rbind(res, data.table(dataset = label, ncontacts = length(gr), DHS = DHS_label, observed_overlaps = observed, expected_overlaps = expected, pval = pval))

  if (plot)
  {
    p <- ggplot(data.table(shuffled = shuffled), aes(shuffled)) +
      labs(x = plot_xlab, y = "Count") +
      geom_histogram(binwidth = 0.001, boundary = 0, color = NA, fill = "#c0c0c0") +
      geom_vline(data = data.table(expected), aes(xintercept = expected), lty = 2) +
      geom_vline(data = data.table(observed), aes(xintercept = observed), color = "darkorange") +
      geom_text(data = data.table(observed),
        aes(label = ifelse(pval < 1 / num_shuffles, sprintf("  p < %0.3f", 1 / num_shuffles), sprintf("  p = %0.3f", pval)),
        x = observed, y = Inf), hjust = 0, vjust = 1, color = "darkorange")
    print(p)
  }

  invisible(NULL)
}


gr <- diffCaptureC_gr

pdf("analysis/balancer/diffCaptureC_vs_DHS.pdf", width = 5.1, height = 2.5)
calculate_pvalue(gr, DHS_gr, "Capture-C differential contacts", "DHS stages 9-11", plot = T, plot_xlab = "Fraction of Capture\u00adC differential contacts that overlaps DHSes")
dev.off()
pdf("analysis/balancer/diffCaptureC_vs_DHS_James.pdf", width = 5.1, height = 2.5)
calculate_pvalue(gr, DHS_James_gr, "Capture-C differential contacts", "DHS_James", plot = T, plot_xlab = "Fraction of Capture\u00adC differential contacts that overlaps merged DHSes (James)")
dev.off()
# calculate_pvalue(gr, read_DHS_gr("stage5"), "Capture-C differential contacts", "DHS stage 5")
# calculate_pvalue(gr, read_DHS_gr("stage9"), "Capture-C differential contacts", "DHS stage 9")
# calculate_pvalue(gr, read_DHS_gr("stage10"), "Capture-C differential contacts", "DHS stage 10")
# calculate_pvalue(gr, read_DHS_gr("stage11"), "Capture-C differential contacts", "DHS stage 11")
# calculate_pvalue(gr, read_DHS_gr("stage14"), "Capture-C differential contacts", "DHS stage 14")
calculate_pvalue(gr, CAD4_gr, "Capture-C differential contacts", "CAD4")

pdf("analysis/balancer/diffCaptureC_vs_DHS_plus_1kb.pdf", width = 5.1, height = 2.5)
calculate_pvalue(gr, DHS_gr + 1e3, "Capture-C differential contacts", "DHS stages 9-11, with 1 kb margin", plot = T, plot_xlab = "Fraction of Capture\u00adC differential contacts that is < 1 kb from a DHS")
dev.off()
pdf("analysis/balancer/diffCaptureC_vs_DHS_James_plus_1kb.pdf", width = 5.1, height = 2.5)
calculate_pvalue(gr, DHS_James_gr + 1e3, "Capture-C differential contacts", "DHS_James, with 1 kb margin", plot = T, plot_xlab = "Fraction of Capture\u00adC differential contacts that is < 1 kb from a merged DHS (James)")
dev.off()
calculate_pvalue(gr, CAD4_gr + 1e3, "Capture-C differential contacts", "CAD4, with 1 kb margin")

gr <- diffCaptureC_gr[unique(genes_to_diffCaptureC_dt$center_id)]

calculate_pvalue(gr, DHS_gr, "Capture-C differential contacts assigned to genes", "DHS stages 9-11")
calculate_pvalue(gr, CAD4_gr, "Capture-C differential contacts assigned to genes", "CAD4")

calculate_pvalue(gr, DHS_gr + 1e3, "Capture-C differential contacts assigned to genes", "DHS stages 9-11, with 1 kb margin")
calculate_pvalue(gr, CAD4_gr + 1e3, "Capture-C differential contacts assigned to genes", "CAD4, with 1 kb margin")

gr <- diffCaptureC_gr[unique(genes_to_diffCaptureC_dt$center_id[genes_to_diffCaptureC_dt$signf != "n"])]

calculate_pvalue(gr, DHS_gr, "Capture-C differential contacts assigned to testable genes", "DHS stages 9-11")
calculate_pvalue(gr, CAD4_gr, "Capture-C differential contacts assigned to testable genes", "CAD4")

calculate_pvalue(gr, DHS_gr + 1e3, "Capture-C differential contacts assigned to testable genes", "DHS stages 9-11, with 1 kb margin")
calculate_pvalue(gr, CAD4_gr + 1e3, "Capture-C differential contacts assigned to testable genes", "CAD4, with 1 kb margin")

gr <- diffCaptureC_gr[unique(genes_to_diffCaptureC_dt$center_id[genes_to_diffCaptureC_dt$signf == "s"])]

calculate_pvalue(gr, DHS_gr, "Capture-C differential contacts assigned to DE genes", "DHS stages 9-11")
calculate_pvalue(gr, DHS_James_gr, "Capture-C differential contacts assigned to DE genes", "DHS James")
calculate_pvalue(gr, CAD4_gr, "Capture-C differential contacts assigned to DE genes", "CAD4")

calculate_pvalue(gr, DHS_gr + 1e3, "Capture-C differential contacts assigned to DE genes", "DHS stages 9-11, with 1 kb margin")
calculate_pvalue(gr, DHS_James_gr + 1e3, "Capture-C differential contacts assigned to DE genes", "DHS James, with 1 kb margin")
calculate_pvalue(gr, CAD4_gr + 1e3, "Capture-C differential contacts assigned to DE genes", "CAD4, with 1 kb margin")

gr <- diffCaptureC_gr[unique(genes_to_diffCaptureC_dt$center_id[genes_to_diffCaptureC_dt$signf == "s" & abs(genes_to_diffCaptureC_dt$log2FoldChange) > log2(1.5)])]

calculate_pvalue(gr, DHS_gr, "Capture-C differential contacts assigned to strongly DE genes", "DHS stages 9-11")
calculate_pvalue(gr, CAD4_gr, "Capture-C differential contacts assigned to strongly DE genes", "CAD4")

calculate_pvalue(gr, DHS_gr + 1e3, "Capture-C differential contacts assigned to strongly DE genes", "DHS stages 9-11, with 1 kb margin")
calculate_pvalue(gr, CAD4_gr + 1e3, "Capture-C differential contacts assigned to strongly DE genes", "CAD4, with 1 kb margin")

print(res[pval < 0.01])

#
#  Orthogonal approach: distance to the closest feature (e.g. closest DHS)
#

calculate_DHS_distances <- function(gr, DHS_gr, shuffle = F)
{
  if (shuffle)
    gr <- shuffle_gr(gr)

  h <- distanceToNearest(gr, DHS_gr)
  gh <<- h
  ggr <<- gr
  stopifnot(length(h) == length(gr))
  stopifnot(queryHits(h) == seq_along(gr))

  return(elementMetadata(h)$distance)
}

res <- NULL

calculate_pvalue <- function(gr, DHS_gr, label, DHS_label, num_shuffles = 1000L, plot = F, plot_xlab = NULL)
{
  observed <- calculate_DHS_distances(gr, DHS_gr, shuffle = F)
  # shuffled <- unlist(lapply(1:num_shuffles, function(i) calculate_DHS_distances(gr, DHS_gr, shuffle = T)))
  shuffled <- unlist(mclapply(1:num_shuffles, function(i) calculate_DHS_distances(gr, DHS_gr, shuffle = T), mc.preschedule = T))
  pval <- ks.test(observed, shuffled)$p.value
  res <<- rbind(res, data.table(dataset = label, ncontacts = length(gr), DHS = DHS_label, x0_observed = mean(observed == 0), x0_shuffled = mean(shuffled == 0), mean_observed_dist = mean(observed), mean_shuffled_dist = mean(shuffled), pval = pval))

  if (plot)
  {
    # dt <- rbind(data.table(subset = "observed", value = observed), data.table(subset = "shuffled", value = shuffled))
    # p <- ggplot(dt, aes(value, color = subset, fill = subset)) +
    #   xlim(1, 20e3) +
    #   geom_density(alpha = 0.1)
    # print(p)

    x <- 0:100 * 200
    dt <- rbind(data.table(subset = "observed", x = x, y = ecdf(observed)(x)), data.table(subset = "shuffled", x = x, y = ecdf(shuffled)(x)))
    p <- ggplot(dt, aes(x = x / 1e3, y = y, color = subset)) +
      geom_line() +
      scale_color_carto_d(name = NULL, palette = "Bold") +
      theme(legend.justification = c(1, 0), legend.position = c(1, 0)) +
      geom_text(data = data.table(x = 0 , y = 1), aes(x = x / 1e3, y = y, hjust = 0, vjust = 1,
        label = ifelse(pval < 0.001, sprintf("p = %0.1e", pval), sprintf("p = %0.3f", pval))), inherit.aes = F) +
      labs(title = label, x = plot_xlab, y = "Fraction of differential contacts")
    if (max(width(DHS_gr)) > 1)
      p <- p + geom_point(data = rbind(data.table(subset = "observed", x = 0, y = ecdf(observed)(0)), data.table(subset = "shuffled", x = 0, y = ecdf(shuffled)(0))))
    print(p)
  }

  invisible(NULL)
}

#
#  ...for differential Hi-C contacts
#

pdf("analysis/balancer/diffHiC_dist_to_DHS.pdf", width = 5, height = 2.8)

gr <- diffHiC_gr
label <- paste0("Hi\u00adC differential contacts (", length(unique(gr$diffc_id)), ")")
calculate_pvalue(gr, DHS_gr, label, "DHS stages 9-11", plot = T, plot_xlab = "Distance to DHS (kb)")
# calculate_pvalue(gr, CAD4_gr, label, "CAD4", plot = T, plot_xlab = "Distance to CAD4 entry (kb)")

gr <- diffHiC_gr[with(genes_to_diffHiC_dt, unique(center_id))]
label <- paste0("Hi\u00adC differential contacts assigned to a gene (", length(unique(gr$diffc_id)), ")")
calculate_pvalue(gr, DHS_gr, label, "DHS stages 9-11", plot = T, plot_xlab = "Distance to DHS (kb)")
# calculate_pvalue(gr, CAD4_gr, label, "CAD4", plot = T, plot_xlab = "Distance to CAD4 entry (kb)")

gr <- diffHiC_gr[with(genes_to_diffHiC_dt, unique(center_id[signf != "n"]))]
label <- paste0("Hi\u00adC differential contacts assigned to a testable gene (", length(unique(gr$diffc_id)), ")")
calculate_pvalue(gr, DHS_gr, label, "DHS stages 9-11", plot = T, plot_xlab = "Distance to DHS (kb)")
# calculate_pvalue(gr, CAD4_gr, label, "CAD4", plot = T, plot_xlab = "Distance to CAD4 entry (kb)")

gr <- diffHiC_gr[with(genes_to_diffHiC_dt, unique(center_id[signf == "s"]))]
label <- paste0("Hi\u00adC differential contacts assigned to a DE gene (", length(unique(gr$diffc_id)), ")")
calculate_pvalue(gr, DHS_gr, label, "DHS stages 9-11", plot = T, plot_xlab = "Distance to DHS (kb)")
# calculate_pvalue(gr, CAD4_gr, label, "CAD4", plot = T, plot_xlab = "Distance to CAD4 entry (kb)")

gr <- diffHiC_gr[with(genes_to_diffHiC_dt, unique(center_id[signf == "s" & abs(log2FoldChange) > log2(1.5)]))]
label <- paste0("Hi\u00adC differential contacts assigned to a strongly DE gene (", length(unique(gr$diffc_id)), ")")
calculate_pvalue(gr, DHS_gr, label, "DHS stages 9-11", plot = T, plot_xlab = "Distance to DHS (kb)")
# calculate_pvalue(gr, CAD4_gr, label, "CAD4", plot = T, plot_xlab = "Distance to CAD4 entry (kb)")

gr <- diffHiC_gr[with(genes_to_diffHiC_dt, unique(center_id[signf == "s" & log2FoldChange > 0]))]
label <- paste0("Hi\u00adC differential contacts assigned to an upregulated DE gene (", length(unique(gr$diffc_id)), ")")
calculate_pvalue(gr, DHS_gr, label, "DHS stages 9-11", plot = T, plot_xlab = "Distance to DHS (kb)")
# calculate_pvalue(gr, CAD4_gr, label, "CAD4", plot = T, plot_xlab = "Distance to CAD4 entry (kb)")

gr <- diffHiC_gr[with(genes_to_diffHiC_dt, unique(center_id[signf == "s" & log2FoldChange < 0]))]
label <- paste0("Hi\u00adC differential contacts assigned to a downregulated DE gene (", length(unique(gr$diffc_id)), ")")
calculate_pvalue(gr, DHS_gr, label, "DHS stages 9-11", plot = T, plot_xlab = "Distance to DHS (kb)")
# calculate_pvalue(gr, CAD4_gr, label, "CAD4", plot = T, plot_xlab = "Distance to CAD4 entry (kb)")

gr <- diffHiC_gr[with(genes_to_diffHiC_dt, unique(center_id[signf == "i"]))]
label <- paste0("Hi\u00adC differential contacts assigned to a non\u00adDE gene (", length(unique(gr$diffc_id)), ")")
calculate_pvalue(gr, DHS_gr, label, "DHS stages 9-11", plot = T, plot_xlab = "Distance to DHS (kb)")
# calculate_pvalue(gr, CAD4_gr, label, "CAD4", plot = T, plot_xlab = "Distance to CAD4 entry (kb)")

dev.off()


pdf("analysis/balancer/diffHiC_dist_to_DHS_James.pdf", width = 5, height = 2.8)

gr <- diffHiC_gr
label <- paste0("Hi\u00adC differential contacts (", length(unique(gr$diffc_id)), ")")
calculate_pvalue(gr, DHS_James_gr, label, "DHS_James", plot = T, plot_xlab = "Distance to James' DHS (kb)")

gr <- diffHiC_gr[with(genes_to_diffHiC_dt, unique(center_id))]
label <- paste0("Hi\u00adC differential contacts assigned to a gene (", length(unique(gr$diffc_id)), ")")
calculate_pvalue(gr, DHS_James_gr, label, "DHS_James", plot = T, plot_xlab = "Distance to James' DHS (kb)")

gr <- diffHiC_gr[with(genes_to_diffHiC_dt, unique(center_id[signf != "n"]))]
label <- paste0("Hi\u00adC differential contacts assigned to a testable gene (", length(unique(gr$diffc_id)), ")")
calculate_pvalue(gr, DHS_James_gr, label, "DHS_James", plot = T, plot_xlab = "Distance to James' DHS (kb)")

gr <- diffHiC_gr[with(genes_to_diffHiC_dt, unique(center_id[signf == "s"]))]
label <- paste0("Hi\u00adC differential contacts assigned to a DE gene (", length(unique(gr$diffc_id)), ")")
calculate_pvalue(gr, DHS_James_gr, label, "DHS_James", plot = T, plot_xlab = "Distance to James' DHS (kb)")

gr <- diffHiC_gr[with(genes_to_diffHiC_dt, unique(center_id[signf == "s" & abs(log2FoldChange) > log2(1.5)]))]
label <- paste0("Hi\u00adC differential contacts assigned to a strongly DE gene (", length(unique(gr$diffc_id)), ")")
calculate_pvalue(gr, DHS_James_gr, label, "DHS_James", plot = T, plot_xlab = "Distance to James' DHS (kb)")

gr <- diffHiC_gr[with(genes_to_diffHiC_dt, unique(center_id[signf == "s" & log2FoldChange > 0]))]
label <- paste0("Hi\u00adC differential contacts assigned to an upregulated DE gene (", length(unique(gr$diffc_id)), ")")
calculate_pvalue(gr, DHS_James_gr, label, "DHS_James", plot = T, plot_xlab = "Distance to James' DHS (kb)")

gr <- diffHiC_gr[with(genes_to_diffHiC_dt, unique(center_id[signf == "s" & log2FoldChange < 0]))]
label <- paste0("Hi\u00adC differential contacts assigned to a downregulated DE gene (", length(unique(gr$diffc_id)), ")")
calculate_pvalue(gr, DHS_James_gr, label, "DHS_James", plot = T, plot_xlab = "Distance to James' DHS (kb)")

gr <- diffHiC_gr[with(genes_to_diffHiC_dt, unique(center_id[signf == "i"]))]
label <- paste0("Hi\u00adC differential contacts assigned to a non\u00adDE gene (", length(unique(gr$diffc_id)), ")")
calculate_pvalue(gr, DHS_James_gr, label, "DHS_James", plot = T, plot_xlab = "Distance to James' DHS (kb)")

dev.off()


pdf("analysis/balancer/diffHiC_dist_to_CNV.pdf", width = 5, height = 2.8)

gr <- diffHiC_gr
label <- paste0("Hi\u00adC differential contacts (", length(unique(gr$diffc_id)), ")")
calculate_pvalue(gr, CNV_gr, label, "CNV", plot = T, plot_xlab = "Distance to CNV (kb)")

gr <- diffHiC_gr[with(genes_to_diffHiC_dt, unique(center_id))]
label <- paste0("Hi\u00adC differential contacts assigned to a gene (", length(unique(gr$diffc_id)), ")")
calculate_pvalue(gr, CNV_gr, label, "CNV", plot = T, plot_xlab = "Distance to CNV (kb)")

gr <- diffHiC_gr[with(genes_to_diffHiC_dt, unique(center_id[signf != "n"]))]
label <- paste0("Hi\u00adC differential contacts assigned to a testable gene (", length(unique(gr$diffc_id)), ")")
calculate_pvalue(gr, CNV_gr, label, "CNV", plot = T, plot_xlab = "Distance to CNV (kb)")

gr <- diffHiC_gr[with(genes_to_diffHiC_dt, unique(center_id[signf == "s"]))]
label <- paste0("Hi\u00adC differential contacts assigned to a DE gene (", length(unique(gr$diffc_id)), ")")
calculate_pvalue(gr, CNV_gr, label, "CNV", plot = T, plot_xlab = "Distance to CNV (kb)")

gr <- diffHiC_gr[with(genes_to_diffHiC_dt, unique(center_id[signf == "s" & abs(log2FoldChange) > log2(1.5)]))]
label <- paste0("Hi\u00adC differential contacts assigned to a strongly DE gene (", length(unique(gr$diffc_id)), ")")
calculate_pvalue(gr, CNV_gr, label, "CNV", plot = T, plot_xlab = "Distance to CNV (kb)")

gr <- diffHiC_gr[with(genes_to_diffHiC_dt, unique(center_id[signf == "s" & log2FoldChange > 0]))]
label <- paste0("Hi\u00adC differential contacts assigned to an upregulated DE gene (", length(unique(gr$diffc_id)), ")")
calculate_pvalue(gr, CNV_gr, label, "CNV", plot = T, plot_xlab = "Distance to CNV (kb)")

gr <- diffHiC_gr[with(genes_to_diffHiC_dt, unique(center_id[signf == "s" & log2FoldChange > 0]))]
label <- paste0("Hi\u00adC differential contacts assigned to a downregulated DE gene (", length(unique(gr$diffc_id)), ")")
calculate_pvalue(gr, CNV_gr, label, "CNV", plot = T, plot_xlab = "Distance to CNV (kb)")

gr <- diffHiC_gr[with(genes_to_diffHiC_dt, unique(center_id[signf == "i"]))]
label <- paste0("Hi\u00adC differential contacts assigned to a non\u00adDE gene (", length(unique(gr$diffc_id)), ")")
calculate_pvalue(gr, CNV_gr, label, "CNV", plot = T, plot_xlab = "Distance to CNV (kb)")

dev.off()


#
#  More tricky approach: genes stratified by fold change
#

genes[, level := factor(
  paste0(signf,
   ifelse(signf == "s", ifelse(abs(log2FoldChange) <= log2(1.5), "m", ifelse(abs(log2FoldChange) <= log2(3), "h", "v")), "")
  ), c("n", "i", "sm", "sh", "sv"))]
genes_to_diffHiC_dt[, level := genes$level[match(gene_id, genes$gene_id)]]

level.labels <- c("n" = "genes not tested",
  "i" = "non\u00adDE genes",
  "iu" = "upregulated non\u00adDE genes", "id" = "downregulated non\u00adDE genes",
  "s" = "DE genes",
  "sm" = "DE genes, moderate fold change", "sh" = "DE genes, high fold change", "sv" = "DE genes, very high fold change",
  "su" = "upregulated DE genes", "sd" = "downregulated DE genes")

# bre <- as.numeric(round(quantile(abs(unique(genes_to_diffHiC_dt[signf != "n"][, c("gene_id", "log2FoldChange")])$log2FoldChange), c(1, 3, 5) / 7), 2))
# bre <- c(-Inf, -rev(bre), bre, Inf)
# genes_to_diffHiC_dt$level <- cut(genes_to_diffHiC_dt$log2FoldChange, breaks = bre)

set.seed(4242)
pdf(paste0("analysis/balancer/genes_by_fold_change_vs_genomic_features.pdf"), width = 5, height = 2, onefile = T)

for (l in rev(c("n", "i", "s", "sm", "sh", "sv")))
{
  message("level: ", l)
  glm_pvalues <<- NULL
  stopifnot(grepl("^chr[23]", genes$chrom))

  genes_subset <- genes[grepl(paste0("^", l), level), ]
  # genes_subset_diffCaptureC <- genes_diffCaptureC[grepl(paste0("^", l), level), ]

  genes_subset <- sort_GRanges_by_log2FoldChange(genes_subset)
  # genes_subset_diffCaptureC <- sort_GRanges_by_log2FoldChange(genes_subset_diffCaptureC)

  # print(make_composite_plot_genes(genes_subset, diffHiC_gr, bp_gr, ggtitle = "Differential Hi\u00adC contacts from the TSS",
  #   extra_filtering = "bait", ribbon_label = "Differential\ncontact density", legend_name = "Hi\u00adC contact\nlog2 fold change"))
  ggtitle <- paste0("Differential Hi\u00adC contacts from the TSS\n", level.labels[l])
  gmap <- extract_map_genes(genes_subset, diffHiC_gr, bp_gr, ggtitle = ggtitle,
    mode = "all", fc_method = "most_significant", extra_filtering = "bait")
  print(make_ribbon_plot_genes(genes_subset, diffHiC_gr, gmap, ymax = 0.15, ggtitle = ggtitle))

  # print(make_composite_plot_genes(genes_subset_diffCaptureC, diffCaptureC_gr, bp_gr, ggtitle = "Differential Capture\u00adC contacts from the TSS",
  #   extra_filtering = "gene_id", sample_size = Inf, ribbon_label = "Differential\ncontact density", legend_name = "Capture\u00adC contact\nlog2 fold change", ylabeller = function(x) paste0("  ", format(x, digits = 1))))
}

dev.off()

#
#  one would expect that for each gene, 1 interaction with a DHS is lost while another one is gained. Is it the case?
#

dt <- NULL

TSS_subset <- TSS[!is.na(signf), ]
TSS_subset_to_diffHiC_dt <- logFC_vs_dist(TSS_subset, diffHiC_gr, bp_gr, mode = "all", extra_filtering = "bait", max_dist = 99e6)

for (l in c("s", "i", "n"))
{
  gr <- diffHiC_gr[with(TSS_subset_to_diffHiC_dt, unique(center_id[signf == l]))]
  # gr <- diffHiC_gr[with(TSS_subset_to_diffHiC_dt, unique(center_id[abs(log2FoldChange) > log2(1.5) & signf == l]))]
  d <- as.data.table(gr)[, list(b = sum(log2FoldChange > 0), w = sum(log2FoldChange < 0)), by = c("baitChr", "baitStart", "baitEnd")]
  dt <- rbind(dt, with(d, data.table(
    dataset = "Hi-C",
    signf = l,
    signf.label = level.labels[l],
    n = nrow(d),
    n_multiple_DC = sum(w + b >= 2),
    n_multiple_DC_different_direction = sum(w + b >= 2 & b * w != 0),
    fraction_multiple_DC = sum(w + b >= 2) / nrow(d),
    fraction_multiple_DC_different_direction = sum(w + b >= 2 & b * w != 0) / nrow(d)
  )))
}

for (l in c("s", "i", "n"))
{
  gr <- diffCaptureC_gr[diffCaptureC_gr$gene_id %in% with(genes, gene_id[signf == l])]
  # gr <- diffCaptureC_gr[diffCaptureC_gr$gene_id %in% with(genes, gene_id[abs(log2FoldChange) > log2(1.5) & signf == l])]
  d <- as.data.table(gr)[, list(b = sum(log2FoldChange > 0), w = sum(log2FoldChange < 0)), by = c("baitChr", "baitStart", "baitEnd")]
  dt <- rbind(dt, with(d, data.table(
    dataset = "Capture-C",
    signf = l,
    signf.label = level.labels[l],
    n = nrow(d),
    n_multiple_DC = sum(w + b >= 2),
    n_multiple_DC_different_direction = sum(w + b >= 2 & b * w != 0),
    fraction_multiple_DC = sum(w + b >= 2) / nrow(d),
    fraction_multiple_DC_different_direction = sum(w + b >= 2 & b * w != 0) / nrow(d)
  )))
}

print(dt)

#
#  ...for differential Capture-C contacts
#

# diffCaptureC_nonMEI_gr <- diffCaptureC_gr[diffCaptureC_gr$gene_id %in% genes[MEI == F]$gene_id]

gr <- diffCaptureC_gr
tt <- table(gr$gene_id)
tt <- names(tt[tt > 30])
gr2 <- gr[!(gr$gene_id %in% tt)]
gr3 <- gr[gr$gene_id %in% genes[signf == "s"]$gene_id]
gr4 <- gr2[gr2$gene_id %in% genes[signf == "s"]$gene_id]


pdf("analysis/balancer/diffCaptureC_dist_to_DHS.pdf", width = 5, height = 2.8)

gr <- diffCaptureC_gr
label <- paste0("Capture\u00adC differential contacts assigned to a gene (", length(gr), ")")
calculate_pvalue(gr, DHS_gr, label, "DHS stages 9-11", plot = T, plot_xlab = "Distance to DHS (kb)")
# calculate_pvalue(gr, CAD4_gr, label, "CAD4", plot = T, plot_xlab = "Distance to CAD4 entry (kb)")

gr <- diffCaptureC_gr[diffCaptureC_gr$gene_id %in% genes[signf != "n"]$gene_id]
label <- paste0("Capture\u00adC differential contacts assigned to a testable gene (", length(gr), ")")
calculate_pvalue(gr, DHS_gr, label, "DHS stages 9-11", plot = T, plot_xlab = "Distance to DHS (kb)")
# calculate_pvalue(gr, CAD4_gr, label, "CAD4", plot = T, plot_xlab = "Distance to CAD4 entry (kb)")

gr <- diffCaptureC_gr[diffCaptureC_gr$gene_id %in% genes[signf == "s"]$gene_id]
label <- paste0("Capture\u00adC differential contacts assigned to a DE gene (", length(gr), ")")
calculate_pvalue(gr, DHS_gr, label, "DHS stages 9-11", plot = T, plot_xlab = "Distance to DHS (kb)")
# calculate_pvalue(gr, CAD4_gr, label, "CAD4", plot = T, plot_xlab = "Distance to CAD4 entry (kb)")

gr <- diffCaptureC_gr[diffCaptureC_gr$gene_id %in% genes[signf == "s" & abs(log2FoldChange) > log2(1.5)]$gene_id]
label <- paste0("Capture\u00adC differential contacts assigned to a strongly DE gene (", length(gr), ")")
calculate_pvalue(gr, DHS_gr, label, "DHS stages 9-11", plot = T, plot_xlab = "Distance to DHS (kb)")
# calculate_pvalue(gr, CAD4_gr, label, "CAD4", plot = T, plot_xlab = "Distance to CAD4 entry (kb)")

gr <- diffCaptureC_gr[diffCaptureC_gr$gene_id %in% genes[signf == "s" & log2FoldChange > 0]$gene_id]
label <- paste0("Capture\u00adC differential contacts assigned to an upregulated DE gene (", length(gr), ")")
calculate_pvalue(gr, DHS_gr, label, "DHS stages 9-11", plot = T, plot_xlab = "Distance to DHS (kb)")
# calculate_pvalue(gr, CAD4_gr, label, "CAD4", plot = T, plot_xlab = "Distance to CAD4 entry (kb)")

gr <- diffCaptureC_gr[diffCaptureC_gr$gene_id %in% genes[signf == "s" & log2FoldChange < 0]$gene_id]
label <- paste0("Capture\u00adC differential contacts assigned to a downregulated DE gene (", length(gr), ")")
calculate_pvalue(gr, DHS_gr, label, "DHS stages 9-11", plot = T, plot_xlab = "Distance to DHS (kb)")
# calculate_pvalue(gr, CAD4_gr, label, "CAD4", plot = T, plot_xlab = "Distance to CAD4 entry (kb)")

gr <- diffCaptureC_gr[diffCaptureC_gr$gene_id %in% genes[signf == "i"]$gene_id]
label <- paste0("Capture\u00adC differential contacts assigned to a non\u00adDE gene (", length(gr), ")")
calculate_pvalue(gr, DHS_gr, label, "DHS stages 9-11", plot = T, plot_xlab = "Distance to DHS (kb)")
# calculate_pvalue(gr, CAD4_gr, label, "CAD4", plot = T, plot_xlab = "Distance to CAD4 entry (kb)")

dev.off()


pdf("analysis/balancer/diffCaptureC_dist_to_DHS_James.pdf", width = 5, height = 2.8)

gr <- diffCaptureC_gr
label <- paste0("Capture\u00adC differential contacts assigned to a gene (", length(gr), ")")
calculate_pvalue(gr, DHS_James_gr, label, "DHS_James", plot = T, plot_xlab = "Distance to James' DHS (kb)")

gr <- diffCaptureC_gr[diffCaptureC_gr$gene_id %in% genes[signf != "n"]$gene_id]
label <- paste0("Capture\u00adC differential contacts assigned to a testable gene (", length(gr), ")")
calculate_pvalue(gr, DHS_James_gr, label, "DHS_James", plot = T, plot_xlab = "Distance to James' DHS (kb)")

gr <- diffCaptureC_gr[diffCaptureC_gr$gene_id %in% genes[signf == "s"]$gene_id]
label <- paste0("Capture\u00adC differential contacts assigned to a DE gene (", length(gr), ")")
calculate_pvalue(gr, DHS_James_gr, label, "DHS_James", plot = T, plot_xlab = "Distance to James' DHS (kb)")

gr <- diffCaptureC_gr[diffCaptureC_gr$gene_id %in% genes[signf == "s" & abs(log2FoldChange) > log2(1.5)]$gene_id]
label <- paste0("Capture\u00adC differential contacts assigned to a strongly DE gene (", length(gr), ")")
calculate_pvalue(gr, DHS_James_gr, label, "DHS_James", plot = T, plot_xlab = "Distance to James' DHS (kb)")

gr <- diffCaptureC_gr[diffCaptureC_gr$gene_id %in% genes[signf == "s" & log2FoldChange > 0]$gene_id]
label <- paste0("Capture\u00adC differential contacts assigned to an upregulated DE gene (", length(gr), ")")
calculate_pvalue(gr, DHS_James_gr, label, "DHS_James", plot = T, plot_xlab = "Distance to James' DHS (kb)")

gr <- diffCaptureC_gr[diffCaptureC_gr$gene_id %in% genes[signf == "s" & log2FoldChange < 0]$gene_id]
label <- paste0("Capture\u00adC differential contacts assigned to a downregulated DE gene (", length(gr), ")")
calculate_pvalue(gr, DHS_James_gr, label, "DHS_James", plot = T, plot_xlab = "Distance to James' DHS (kb)")

gr <- diffCaptureC_gr[diffCaptureC_gr$gene_id %in% genes[signf == "i"]$gene_id]
label <- paste0("Capture\u00adC differential contacts assigned to a non\u00adDE gene (", length(gr), ")")
calculate_pvalue(gr, DHS_James_gr, label, "DHS_James", plot = T, plot_xlab = "Distance to James' DHS (kb)")

dev.off()


pdf("analysis/balancer/diffCaptureC_dist_to_DHS_control.pdf", width = 5, height = 2.8)
calculate_pvalue(gr2, DHS_gr, "Capture-C differential contacts", "DHS stages 9-11", plot = T, plot_xlab = "Distance to DHS (kb)")
calculate_pvalue(gr3, DHS_gr, "Capture-C differential contacts", "DHS stages 9-11", plot = T, plot_xlab = "Distance to DHS (kb)")
calculate_pvalue(gr4, DHS_gr, "Capture-C differential contacts", "DHS stages 9-11", plot = T, plot_xlab = "Distance to DHS (kb)")
calculate_pvalue(gr2, DHS_James_gr, "Capture-C differential contacts", "DHS_James", plot = T, plot_xlab = "Distance to James' DHS (kb)")
calculate_pvalue(gr3, DHS_James_gr, "Capture-C differential contacts", "DHS_James", plot = T, plot_xlab = "Distance to James' DHS (kb)")
calculate_pvalue(gr4, DHS_James_gr, "Capture-C differential contacts", "DHS_James", plot = T, plot_xlab = "Distance to James' DHS (kb)")
calculate_pvalue(gr2, CAD4_gr, "Capture-C differential contacts", "CAD4", plot = T, plot_xlab = "Distance to CAD4 entry (kb)")
calculate_pvalue(gr3, CAD4_gr, "Capture-C differential contacts", "CAD4", plot = T, plot_xlab = "Distance to CAD4 entry (kb)")
calculate_pvalue(gr4, CAD4_gr, "Capture-C differential contacts", "CAD4", plot = T, plot_xlab = "Distance to CAD4 entry (kb)")
dev.off()


pdf("analysis/balancer/diffCaptureC_dist_to_CNV.pdf", width = 5, height = 2.8)

gr <- diffCaptureC_gr
label <- paste0("Capture\u00adC differential contacts assigned to a gene (", length(gr), ")")
calculate_pvalue(gr, CNV_gr, label, "CNV", plot = T, plot_xlab = "Distance to CNV (kb)")

gr <- diffCaptureC_gr[diffCaptureC_gr$gene_id %in% genes[signf != "n"]$gene_id]
label <- paste0("Capture\u00adC differential contacts assigned to a testable gene (", length(gr), ")")
calculate_pvalue(gr, CNV_gr, label, "CNV", plot = T, plot_xlab = "Distance to CNV (kb)")

gr <- diffCaptureC_gr[diffCaptureC_gr$gene_id %in% genes[signf == "s"]$gene_id]
label <- paste0("Capture\u00adC differential contacts assigned to a DE gene (", length(gr), ")")
calculate_pvalue(gr, CNV_gr, label, "CNV", plot = T, plot_xlab = "Distance to CNV (kb)")

gr <- diffCaptureC_gr[diffCaptureC_gr$gene_id %in% genes[signf == "s" & abs(log2FoldChange) > log2(1.5)]$gene_id]
label <- paste0("Capture\u00adC differential contacts assigned to a strongly DE gene (", length(gr), ")")
calculate_pvalue(gr, CNV_gr, label, "CNV", plot = T, plot_xlab = "Distance to CNV (kb)")

gr <- diffCaptureC_gr[diffCaptureC_gr$gene_id %in% genes[signf == "s" & log2FoldChange > 0]$gene_id]
label <- paste0("Capture\u00adC differential contacts assigned to an upregulated DE gene (", length(gr), ")")
calculate_pvalue(gr, CNV_gr, label, "CNV", plot = T, plot_xlab = "Distance to CNV (kb)")

gr <- diffCaptureC_gr[diffCaptureC_gr$gene_id %in% genes[signf == "s" & log2FoldChange < 0]$gene_id]
label <- paste0("Capture\u00adC differential contacts assigned to a downregulated DE gene (", length(gr), ")")
calculate_pvalue(gr, CNV_gr, label, "CNV", plot = T, plot_xlab = "Distance to CNV (kb)")

gr <- diffCaptureC_gr[diffCaptureC_gr$gene_id %in% genes[signf == "i"]$gene_id]
label <- paste0("Capture\u00adC differential contacts assigned to a non\u00adDE gene (", length(gr), ")")
calculate_pvalue(gr, CNV_gr, label, "CNV", plot = T, plot_xlab = "Distance to CNV (kb)")

dev.off()


pdf("analysis/balancer/diffCaptureC_dist_to_CNV_control.pdf", width = 5, height = 2.8)
calculate_pvalue(gr2, CNV_gr, "Capture-C differential contacts", "CNV", plot = T, plot_xlab = "Distance to CNV (kb)")
calculate_pvalue(gr3, CNV_gr, "Capture-C differential contacts", "CNV", plot = T, plot_xlab = "Distance to CNV (kb)")
calculate_pvalue(gr4, CNV_gr, "Capture-C differential contacts", "CNV", plot = T, plot_xlab = "Distance to CNV (kb)")
dev.off()


pdf("analysis/balancer/diffCaptureC_dist_to_promoters.pdf", width = 5, height = 2.8)

gr <- diffCaptureC_gr
label <- paste0("Capture\u00adC differential contacts assigned to a gene (", length(gr), ")")
calculate_pvalue(gr, TSS_no_bait_gr, label, "promoter", plot = T, plot_xlab = "Distance to promoter (kb)")
calculate_pvalue(gr, TSS_no_bait_gr[TSS_no_bait_gr$signf != "n"], label, "promoter expressed", plot = T, plot_xlab = "Distance to promoter of a testable gene (kb)")
calculate_pvalue(gr, TSS_no_bait_gr[TSS_no_bait_gr$signf == "s"], label, "promoter DE", plot = T, plot_xlab = "Distance to promoter of a DE gene (kb)")
calculate_pvalue(gr, TSS_no_bait_gr[TSS_no_bait_gr$signf == "i"], label, "promoter non-DE", plot = T, plot_xlab = "Distance to promoter of a non\u00adDE gene (kb)")

gr <- diffCaptureC_gr[diffCaptureC_gr$gene_id %in% genes[signf != "n"]$gene_id]
label <- paste0("Capture\u00adC differential contacts assigned to a testable gene (", length(gr), ")")
calculate_pvalue(gr, TSS_no_bait_gr, label, "promoter", plot = T, plot_xlab = "Distance to promoter (kb)")
calculate_pvalue(gr, TSS_no_bait_gr[TSS_no_bait_gr$signf != "n"], label, "promoter expressed", plot = T, plot_xlab = "Distance to promoter of a testable gene (kb)")
calculate_pvalue(gr, TSS_no_bait_gr[TSS_no_bait_gr$signf == "s"], label, "promoter DE", plot = T, plot_xlab = "Distance to promoter of a DE gene (kb)")
calculate_pvalue(gr, TSS_no_bait_gr[TSS_no_bait_gr$signf == "i"], label, "promoter non-DE", plot = T, plot_xlab = "Distance to promoter of a non\u00adDE gene (kb)")

gr <- diffCaptureC_gr[diffCaptureC_gr$gene_id %in% genes[signf == "s"]$gene_id]
label <- paste0("Capture\u00adC differential contacts assigned to a DE gene (", length(gr), ")")
calculate_pvalue(gr, TSS_no_bait_gr, label, "promoter", plot = T, plot_xlab = "Distance to promoter (kb)")
calculate_pvalue(gr, TSS_no_bait_gr[TSS_no_bait_gr$signf != "n"], label, "promoter expressed", plot = T, plot_xlab = "Distance to promoter of a testable gene (kb)")
calculate_pvalue(gr, TSS_no_bait_gr[TSS_no_bait_gr$signf == "s"], label, "promoter DE", plot = T, plot_xlab = "Distance to promoter of a DE gene (kb)")
calculate_pvalue(gr, TSS_no_bait_gr[TSS_no_bait_gr$signf == "i"], label, "promoter non-DE", plot = T, plot_xlab = "Distance to promoter of a non\u00adDE gene (kb)")

gr <- diffCaptureC_gr[diffCaptureC_gr$gene_id %in% genes[signf == "i"]$gene_id]
label <- paste0("Capture\u00adC differential contacts assigned to a non\u00adDE gene (", length(gr), ")")
calculate_pvalue(gr, TSS_no_bait_gr, label, "promoter", plot = T, plot_xlab = "Distance to promoter (kb)")
calculate_pvalue(gr, TSS_no_bait_gr[TSS_no_bait_gr$signf != "n"], label, "promoter expressed", plot = T, plot_xlab = "Distance to promoter of a testable gene (kb)")
calculate_pvalue(gr, TSS_no_bait_gr[TSS_no_bait_gr$signf == "s"], label, "promoter DE", plot = T, plot_xlab = "Distance to promoter of a DE gene (kb)")
calculate_pvalue(gr, TSS_no_bait_gr[TSS_no_bait_gr$signf == "i"], label, "promoter non-DE", plot = T, plot_xlab = "Distance to promoter of a non\u00adDE gene (kb)")

dev.off()


pdf("analysis/balancer/diffCaptureC_dist_to_promoters_control.pdf", width = 5, height = 2.8)
calculate_pvalue(gr2, TSS_no_bait_gr, "Capture-C differential contacts", "promoter", plot = T, plot_xlab = "Distance to promoter (kb)")
calculate_pvalue(gr3, TSS_no_bait_gr, "Capture-C differential contacts", "promoter", plot = T, plot_xlab = "Distance to promoter (kb)")
calculate_pvalue(gr4, TSS_no_bait_gr, "Capture-C differential contacts", "promoter", plot = T, plot_xlab = "Distance to promoter (kb)")
calculate_pvalue(gr2, TSS_no_bait_gr[TSS_no_bait_gr$signf != "n"], "Capture-C differential contacts", "promoter expressed", plot = T, plot_xlab = "Distance to promoter of a testable gene (kb)")
calculate_pvalue(gr3, TSS_no_bait_gr[TSS_no_bait_gr$signf != "n"], "Capture-C differential contacts", "promoter expressed", plot = T, plot_xlab = "Distance to promoter of a testable gene (kb)")
calculate_pvalue(gr4, TSS_no_bait_gr[TSS_no_bait_gr$signf != "n"], "Capture-C differential contacts", "promoter expressed", plot = T, plot_xlab = "Distance to promoter of a testable gene (kb)")
calculate_pvalue(gr2, TSS_no_bait_gr[TSS_no_bait_gr$signf == "s"], "Capture-C differential contacts", "promoter DE", plot = T, plot_xlab = "Distance to promoter of a DE gene (kb)")
calculate_pvalue(gr3, TSS_no_bait_gr[TSS_no_bait_gr$signf == "s"], "Capture-C differential contacts", "promoter DE", plot = T, plot_xlab = "Distance to promoter of a DE gene (kb)")
calculate_pvalue(gr4, TSS_no_bait_gr[TSS_no_bait_gr$signf == "s"], "Capture-C differential contacts", "promoter DE", plot = T, plot_xlab = "Distance to promoter of a DE gene (kb)")
calculate_pvalue(gr2, TSS_no_bait_gr[TSS_no_bait_gr$signf == "i"], "Capture-C differential contacts", "promoter DE", plot = T, plot_xlab = "Distance to promoter of a non\u00adDE gene (kb)")
calculate_pvalue(gr3, TSS_no_bait_gr[TSS_no_bait_gr$signf == "i"], "Capture-C differential contacts", "promoter DE", plot = T, plot_xlab = "Distance to promoter of a non\u00adDE gene (kb)")
calculate_pvalue(gr4, TSS_no_bait_gr[TSS_no_bait_gr$signf == "i"], "Capture-C differential contacts", "promoter DE", plot = T, plot_xlab = "Distance to promoter of a non\u00adDE gene (kb)")
dev.off()

#
#  Subset of differential Hi-C contacts that could be detected by Capture-C baits
#

ov <- findOverlaps(GRanges(diffHiC_gr$baitChr, IRanges(diffHiC_gr$baitStart, diffHiC_gr$baitEnd)),
  with(baits, GRanges(chrom, IRanges(start, end))))
diffHiC_CaptureCsubset_gr <- diffHiC_gr[unique(queryHits(ov))]

pdf("analysis/balancer/diffHiC_CaptureCsubset_dist_to_DHS_James.pdf", width = 5, height = 2.8)

gr <- diffHiC_CaptureCsubset_gr
label <- paste0("Hi\u00adC differential contacts, subsetted\nto Capture\u00adC viewpoints (", length(unique(gr$diffc_id)), ")")
calculate_pvalue(gr, DHS_James_gr, label, "DHS_James", plot = T, plot_xlab = "Distance to James' DHS (kb)")

gr <- diffHiC_CaptureCsubset_gr[diffHiC_CaptureCsubset_gr$diffc_id %in% with(genes_to_diffHiC_dt, unique(diffc_id))]
label <- paste0("Hi\u00adC differential contacts, subsetted\nto Capture\u00adC viewpoints, assigned to a gene (", length(unique(gr$diffc_id)), ")")
calculate_pvalue(gr, DHS_James_gr, label, "DHS_James", plot = T, plot_xlab = "Distance to James' DHS (kb)")

gr <- diffHiC_CaptureCsubset_gr[diffHiC_CaptureCsubset_gr$diffc_id %in% with(genes_to_diffHiC_dt, unique(diffc_id[signf != "n"]))]
label <- paste0("Hi\u00adC differential contacts, subsetted\nto Capture\u00adC viewpoints, assigned to a testable gene (", length(unique(gr$diffc_id)), ")")
calculate_pvalue(gr, DHS_James_gr, label, "DHS_James", plot = T, plot_xlab = "Distance to James' DHS (kb)")

gr <- diffHiC_CaptureCsubset_gr[diffHiC_CaptureCsubset_gr$diffc_id %in% with(genes_to_diffHiC_dt, unique(diffc_id[signf == "s"]))]
label <- paste0("Hi\u00adC differential contacts, subsetted\nto Capture\u00adC viewpoints, assigned to a DE gene (", length(unique(gr$diffc_id)), ")")
calculate_pvalue(gr, DHS_James_gr, label, "DHS_James", plot = T, plot_xlab = "Distance to James' DHS (kb)")

gr <- diffHiC_CaptureCsubset_gr[diffHiC_CaptureCsubset_gr$diffc_id %in% with(genes_to_diffHiC_dt, unique(diffc_id[signf == "s" & abs(log2FoldChange) > log2(1.5)]))]
label <- paste0("Hi\u00adC differential contacts, subsetted\nto Capture\u00adC viewpoints, assigned to a strongly DE gene (", length(unique(gr$diffc_id)), ")")
calculate_pvalue(gr, DHS_James_gr, label, "DHS_James", plot = T, plot_xlab = "Distance to James' DHS (kb)")

gr <- diffHiC_CaptureCsubset_gr[diffHiC_CaptureCsubset_gr$diffc_id %in% with(genes_to_diffHiC_dt, unique(diffc_id[signf == "s" & log2FoldChange > 0]))]
label <- paste0("Hi\u00adC differential contacts, subsetted\nto Capture\u00adC viewpoints, assigned to an upregulated DE gene (", length(unique(gr$diffc_id)), ")")
calculate_pvalue(gr, DHS_James_gr, label, "DHS_James", plot = T, plot_xlab = "Distance to James' DHS (kb)")

gr <- diffHiC_CaptureCsubset_gr[diffHiC_CaptureCsubset_gr$diffc_id %in% with(genes_to_diffHiC_dt, unique(diffc_id[signf == "s" & log2FoldChange < 0]))]
label <- paste0("Hi\u00adC differential contacts, subsetted\nto Capture\u00adC viewpoints, assigned to a downregulated DE gene (", length(unique(gr$diffc_id)), ")")
calculate_pvalue(gr, DHS_James_gr, label, "DHS_James", plot = T, plot_xlab = "Distance to James' DHS (kb)")

gr <- diffHiC_CaptureCsubset_gr[diffHiC_CaptureCsubset_gr$diffc_id %in% with(genes_to_diffHiC_dt, unique(diffc_id[signf == "i"]))]
label <- paste0("Hi\u00adC differential contacts, subsetted\nto Capture\u00adC viewpoints, assigned to a non\u00adDE gene (", length(unique(gr$diffc_id)), ")")
calculate_pvalue(gr, DHS_James_gr, label, "DHS_James", plot = T, plot_xlab = "Distance to James' DHS (kb)")

dev.off()

#
#  Taking only short/long stretches of Capture-C contacts
#

# extend by +/- 1 kb, and reduce
gr <- diffCaptureC_gr + 1e3
gr_reduced <- GRanges(as.data.table(gr)[, as.data.table(reduce(GRanges(.SD))),
  by = c("baitID", "gene_id", "baitChr", "baitStart", "baitEnd", "baitName")])

pdf("analysis/balancer/diffCaptureC_clustered_dist_to_DHS_James.pdf", width = 5, height = 2.8)

gr <- gr_reduced
label <- paste0("Capture\u00adC differential contacts assigned to a gene (", length(gr), ")")
calculate_pvalue(gr, DHS_gr, label, "DHS stages 9-11", plot = T, plot_xlab = "Distance to DHS (kb)")
# calculate_pvalue(gr, CAD4_gr, label, "CAD4", plot = T, plot_xlab = "Distance to CAD4 entry (kb)")

gr <- gr_reduced[gr_reduced$gene_id %in% genes[signf != "n"]$gene_id]
label <- paste0("Capture\u00adC differential contacts assigned to a testable gene (", length(gr), ")")
calculate_pvalue(gr, DHS_gr, label, "DHS stages 9-11", plot = T, plot_xlab = "Distance to DHS (kb)")
# calculate_pvalue(gr, CAD4_gr, label, "CAD4", plot = T, plot_xlab = "Distance to CAD4 entry (kb)")

gr <- gr_reduced[gr_reduced$gene_id %in% genes[signf == "s"]$gene_id]
label <- paste0("Capture\u00adC differential contacts assigned to a DE gene (", length(gr), ")")
calculate_pvalue(gr, DHS_gr, label, "DHS stages 9-11", plot = T, plot_xlab = "Distance to DHS (kb)")
# calculate_pvalue(gr, CAD4_gr, label, "CAD4", plot = T, plot_xlab = "Distance to CAD4 entry (kb)")

gr <- gr_reduced[gr_reduced$gene_id %in% genes[signf == "s" & abs(log2FoldChange) > log2(1.5)]$gene_id]
label <- paste0("Capture\u00adC differential contacts assigned to a strongly DE gene (", length(gr), ")")
calculate_pvalue(gr, DHS_gr, label, "DHS stages 9-11", plot = T, plot_xlab = "Distance to DHS (kb)")
# calculate_pvalue(gr, CAD4_gr, label, "CAD4", plot = T, plot_xlab = "Distance to CAD4 entry (kb)")

gr <- gr_reduced[gr_reduced$gene_id %in% genes[signf == "s" & log2FoldChange > 0]$gene_id]
label <- paste0("Capture\u00adC differential contacts assigned to an upregulated DE gene (", length(gr), ")")
calculate_pvalue(gr, DHS_gr, label, "DHS stages 9-11", plot = T, plot_xlab = "Distance to DHS (kb)")
# calculate_pvalue(gr, CAD4_gr, label, "CAD4", plot = T, plot_xlab = "Distance to CAD4 entry (kb)")

gr <- gr_reduced[gr_reduced$gene_id %in% genes[signf == "s" & log2FoldChange < 0]$gene_id]
label <- paste0("Capture\u00adC differential contacts assigned to a downregulated DE gene (", length(gr), ")")
calculate_pvalue(gr, DHS_gr, label, "DHS stages 9-11", plot = T, plot_xlab = "Distance to DHS (kb)")
# calculate_pvalue(gr, CAD4_gr, label, "CAD4", plot = T, plot_xlab = "Distance to CAD4 entry (kb)")

gr <- gr_reduced[gr_reduced$gene_id %in% genes[signf == "i"]$gene_id]
label <- paste0("Capture\u00adC differential contacts assigned to a non\u00adDE gene (", length(gr), ")")
calculate_pvalue(gr, DHS_gr, label, "DHS stages 9-11", plot = T, plot_xlab = "Distance to DHS (kb)")
# calculate_pvalue(gr, CAD4_gr, label, "CAD4", plot = T, plot_xlab = "Distance to CAD4 entry (kb)")

dev.off()


## GRanges object with 931 ranges and 11 metadata columns:
#         seqnames               ranges strand |     gene_id log2FoldChange         padj DHS_overlap CAD4_overlap     baitChr baitStart   baitEnd  distSign bait_DHS_overlap bait_CAD4_overlap


#
#  Legacy plots
#

pdf("analysis/balancer/diffHiC_vs_DE.pdf", width = 8, height = 5)

make_feature_plot(genes, diffHiC_gr, bp_gr, "Differential Hi-C contacts (p<0.05)", mode = "closest", abs_dist = T, abs_fc = T)
make_feature_plot(genes, diffHiC_gr, bp_gr, "Differential Hi-C contacts (p<0.05)", mode = "closest", abs_dist = T)
make_feature_plot(genes, diffHiC_gr, bp_gr, "Differential Hi-C contacts (p<0.05)", mode = "closest", abs_fc = T)
make_feature_plot(genes, diffHiC_gr, bp_gr, "Differential Hi-C contacts (p<0.05)", mode = "closest")

make_feature_plot(genes, diffHiC_BAL_gr, bp_gr, "Hi-C contacts stronger in balancer (p<0.05)", mode = "closest", abs_dist = T, abs_fc = T)
make_feature_plot(genes, diffHiC_BAL_gr, bp_gr, "Hi-C contacts stronger in balancer (p<0.05)", mode = "closest", abs_dist = T)
make_feature_plot(genes, diffHiC_BAL_gr, bp_gr, "Hi-C contacts stronger in balancer (p<0.05)", mode = "closest", abs_fc = T)
make_feature_plot(genes, diffHiC_BAL_gr, bp_gr, "Hi-C contacts stronger in balancer (p<0.05)", mode = "closest")

make_feature_plot(genes, diffHiC_VRG_gr, bp_gr, "Hi-C contacts stronger in wild\u00adtype (p<0.05)", mode = "closest", abs_dist = T, abs_fc = T)
make_feature_plot(genes, diffHiC_VRG_gr, bp_gr, "Hi-C contacts stronger in wild\u00adtype (p<0.05)", mode = "closest", abs_dist = T)
make_feature_plot(genes, diffHiC_VRG_gr, bp_gr, "Hi-C contacts stronger in wild\u00adtype (p<0.05)", mode = "closest", abs_fc = T)
make_feature_plot(genes, diffHiC_VRG_gr, bp_gr, "Hi-C contacts stronger in wild\u00adtype (p<0.05)", mode = "closest")

dev.off()


dt_diffHiC <- unique(as.data.table(TSS_to_diffHiC_dt)[, c("gene_id", "log2FoldChange", "signf", "diffc_id", "HiClog2FoldChange")])
dt_diffCaptureC <- unique(as.data.table(TSS_to_diffCaptureC_dt)[, c("gene_id", "log2FoldChange", "signf", "frag_id", "CaptureClog2FoldChange")])

xmax <- max(abs(c(TSS_to_diffHiC_dt$HiClog2FoldChange, TSS_to_diffCaptureC_dt$CaptureClog2FoldChange)))
ymax <- max(abs(c(TSS_to_diffHiC_dt$log2FoldChange, TSS_to_diffCaptureC_dt$log2FoldChange)))

pdf("analysis/balancer/diffHiC_vs_DE_xyplot.pdf", width = 3.5, height = 4)

p <- ggplot(dt_diffHiC[signf != "n"], aes(HiClog2FoldChange, log2FoldChange, color = signf)) +
  geom_vline(aes(xintercept = 0), lty = 2) +
  geom_hline(aes(yintercept = 0), lty = 2) +
  geom_point(alpha = 0.6) +
  xlab(expression(paste("Hi\u00adC differential contact ", log[2] * " fold change"))) +
  ylab(expression(paste("Gene expression ", log[2] * " fold change"))) +
  # scale_y_continuous(limits = c(-5, 5), expand = c(0.005, 0), oob = scales::squish) +
  scale_x_continuous(limits = c(-1, 1) * xmax) +
  scale_y_continuous(limits = c(-1, 1) * ymax) +
  coord_fixed() +
  scale_color_manual(values = gene.colors, labels = gene.labels, name = NULL) +
  background_grid(major = "xy", minor = "xy") +
  theme(legend.justification = c(0.5, 1), legend.position = c(0.5, 1)) +
  theme(legend.position = "bottom")
print(p)

dev.off()


pdf("analysis/balancer/diffCaptureC_vs_DE_xyplot.pdf", width = 3.5, height = 4)

p <- ggplot(dt_diffCaptureC[signf != "n"], aes(CaptureClog2FoldChange, log2FoldChange, color = signf)) +
  geom_vline(aes(xintercept = 0), lty = 2) +
  geom_hline(aes(yintercept = 0), lty = 2) +
  geom_point(alpha = 0.6) +
  xlab(expression(paste("Capture\u00adC differential contact ", log[2] * " fold change"))) +
  ylab(expression(paste("Gene expression ", log[2] * " fold change"))) +
  # scale_y_continuous(limits = c(-5, 5), expand = c(0.005, 0), oob = scales::squish) +
  scale_x_continuous(limits = c(-1, 1) * xmax) +
  scale_y_continuous(limits = c(-1, 1) * ymax) +
  coord_fixed() +
  scale_color_manual(values = gene.colors, labels = gene.labels, name = NULL) +
  background_grid(major = "xy", minor = "xy") +
  theme(legend.justification = c(0.5, 1), legend.position = c(0.5, 1)) +
  theme(legend.position = "bottom")
print(p)

dev.off()

message("diffHiC: ", nrow(dt_diffHiC[signf != "n"]), " gene-differential contact pairs, Hi-C vs. expression log2FoldChange: Pearson r = ",
  with(dt_diffHiC[signf != "n"], cor(HiClog2FoldChange, log2FoldChange)))
message("diffCaptureC: ", nrow(dt_diffCaptureC[signf != "n"]), " gene-differential contact pairs, Hi-C vs. expression log2FoldChange: Pearson r = ",
  with(dt_diffCaptureC[signf != "n"], cor(CaptureClog2FoldChange, log2FoldChange)))
