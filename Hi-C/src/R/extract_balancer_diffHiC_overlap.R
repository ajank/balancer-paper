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

#
#  Classification of genomic regions: active TSSes, inactive TSSes, insulators, enhancers (Thomas et al. 2011), other
#

# check the FPKM distribution
FPKM_thr <- 10.5
pdf("analysis/balancer/genes_FPKM_distribution.pdf", width = 5, height = 2.5)
p <- ggplot(genes, aes(log10(FPKM))) +
  geom_density(fill = "gray80") +
  xlab(expression(paste(log[10] * " FPKM"))) +
  ylab("Relative frequency") +
  geom_vline(aes(xintercept = log10(FPKM_thr)), color = RColorBrewer::brewer.pal(9, "Set1")[1], lty = 2) +
  NULL
print(p)
dev.off()

# take the genes with FPKM calculated and FPKM >= FPKM_thr as active
# reasons for missing FPKM values:
#   1. no reads (NA in genes$FPKM)
#   2. not considered by DESeq2 due to overlapping same-strand transcripts (not in "genes")
m <- match(promoter_gr$gene_id, genes$gene_id)
sel <- !is.na(genes$FPKM[m]) & genes$FPKM[m] >= FPKM_thr
TSS_active_gr <- TSS_gr[sel]
TSS_inactive_gr <- TSS_gr[!sel]
promoter_active_gr <- promoter_gr[sel]
promoter_inactive_gr <- promoter_gr[!sel]

# combine insulator binding into clusters
insulator.dir <- "/g/korbel/shared/projects/drosophila_balancer/analyses/tracks/TADs/chipData"
insulator.files <- c(
  "CTCF (C)" = "ChipChip.CTCF.embryo_0-12h.White_K_769.bed",
  "CTCF (N)" = "ChipChip.CTCF.embryo_0-12h.White_K_770.bed",
  "su(Hw)" = "ChipChip.SuHw.embryo_0-12h.White_K_27.bed",
  "Zw5" = "ChipChip.ZW5.embryo_2-4h.Karpen_G_5265.bed",
  "BEAF-32" = "ChipChip.BEAF-32.embryo_0-12h.White_K_21.bed",
  "Cp190" = "ChipChip.CP190.embryo_0-12h.White_K_22.bed"
)

insulator_dt <- NULL
for (n in names(insulator.files))
{
  dt <- fread(paste0(insulator.dir, "/", insulator.files[n]), header = F, sep = "\t")
  dt <- dt[, 1:3]
  names(dt) <- c("chrom", "start", "end")
  dt[, start := start + 1L]
  dt <- dt[grepl("^chr[23]", chrom), ]
  dt$source <- n
  insulator_dt <- rbind(insulator_dt, dt)
}

insulator_gr <- GRanges(insulator_dt)

unmetadata <- function(gr)
{
  elementMetadata(gr) <- NULL
  return(gr)
}

enhancer_gr <- c(unmetadata(CAD4_gr), unmetadata(mesoCRM_gr))

#
#  objects to store all the statistics
#

if (!exists("diffHiC_ov"))
  diffHiC_ov <- list()

#
#  generate all the statistics
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

overlap_to_shuffled <- function(gr_list, diffHiC_subset_list, anno_list, num_shuffles = 1000L)
{
  dt_aggr <- NULL
  dt_shuffled <- NULL

  for (a in names(anno_list))
  {
    anno_gr <- anno_list[[a]]
    for (l in names(gr_list))
    {
      gr <- gr_list[[l]]
      observed <- calculate_anno_overlap(gr, anno_gr, shuffle = F)
      # shuffled <- sapply(1:num_shuffles, function(i) calculate_anno_overlap(gr, anno_gr, shuffle = T))
      shuffled <- unlist(mclapply(1:num_shuffles, function(i) calculate_anno_overlap(gr, anno_gr, shuffle = T), mc.preschedule = T))
      pval <- 2 * min(sum(shuffled < observed), sum(shuffled > observed)) / num_shuffles
      # pval <- binom.test(observed * length(gr), length(gr), mean(shuffled))$p.value

      dt_aggr <- rbind(dt_aggr, data.table(
        anno = a,
        level = l,
        observed = observed,
        mean_shuffled = mean(shuffled),
        pval = pval,
        n_DC_gr = length(gr),
        n_DC = length(unique(gr$diffc_id)),
        n_bait_TSS = length(unique(diffHiC_subset_list[[l]]$ID)),
        n_bait_genes = length(unique(diffHiC_subset_list[[l]]$gene_id)),
        n_bait_dhs = length(unique(diffHiC_subset_list[[l]]$dhs_id)),
        n_anno_gr = length(anno_gr)
      ))
      dt_shuffled <- rbind(dt_shuffled, data.table(anno = a, level = l, shuffled = shuffled))
    }
  }

  dt_aggr[, anno := factor(anno, names(anno_list))]
  dt_shuffled[, anno := factor(anno, names(anno_list))]
  dt_aggr[, level := factor(level, names(gr_list))]
  dt_shuffled[, level := factor(level, names(gr_list))]

  # dt_aggr[, pval_to_next := NA_real_]
  # for (i in seq_len(nrow(dt_aggr) - 1))
  #   if (dt_aggr$anno[i] == dt_aggr$anno[i + 1])
  #     dt_aggr$pval_to_next[i] <- with(dt_aggr, fisher.test(cbind(
  #       c(observed[i], 1 - observed[i]) * n_DC_gr[i],
  #       c(observed[i + 1], 1 - observed[i + 1]) * n_DC_gr[i + 1]
  #     ))$p.value)

  return(list(aggr = dt_aggr, shuffled = dt_shuffled))
}


H3K27ac_dt <- fread("/g/furlong/DataWarehouse/Data/datawarehouse/furlong_data/dmelanogaster/sequencing/ChIP-seq/mesoderm/K27ac/6-8h/peaks/dm6/K27ac_merged_6-8h_BiTS-ChIP_200bp-bw_IDR-5pc_optimal_peaks_vs_H3.bed", header = F)
names(H3K27ac_dt)[1:6] <- c("chrom", "start", "end", "name", "score", "strand")
H3K27ac_dt[, start := start + 1L]
H3K27ac_dt <- H3K27ac_dt[grepl("^chr[23]", chrom), ]
H3K27ac_gr <- GRanges(H3K27ac_dt)

eRNA_dt <- fread("/g/furlong/project/33_Hi-C/data/eRNA/Supp_Data_File_S5_DHS_eRNA_dm6.txt", header = T)
eRNA_dt[, start := start + 1L]
eRNA_dt <- eRNA_dt[grepl("^chr[23]", seqnames), ]
eRNA_active_dt <- eRNA_dt[procap_68h_pos + procap_68h_neg >= 60, ]
eRNA_active_gr <- GRanges(eRNA_active_dt)


# id_list <- c("DE_fc1.5_CaptureC_fc1.5", "DE_fc1.5_CaptureC_reduced_fc1.5", "DE_fc1.5_HiC_fc1.5", "DE_fc1.5_HiC_fc1.5_disruptedTADs", "DE_all_CaptureC_all", "DE_all_CaptureC_reduced_all", "DE_all_HiC_all", "DE_all_HiC_all_disruptedTADs")
id_list <- c("DE_all_HiC_all_activeenh", "DE_all_HiC_all_activeenh_disruptedTADs")
for (id in id_list)
{
  message("\n", id)
  # if (is.list(diffHiC_ov[[id]]))
  #   next
  lev <- c("i", "s") # c("n", "i", "s", "sm", "sh", "sv")

  TSS_subset <- TSS
  TSS_subset[, level := factor(
    paste0(signf,
     ifelse(signf == "s", ifelse(abs(log2FoldChange) <= log2(1.5), "m", ifelse(abs(log2FoldChange) <= log2(3), "h", "v")), "")
    ), c("n", "i", "sm", "sh", "sv"))]

  if (grepl("HiC_all", id))
    diffHiC_gr_subset <- diffHiC_gr
  else if (grepl("HiC_pval0.05_fc1.5", id))
    diffHiC_gr_subset <- diffHiC_gr[diffHiC_gr$padj < 0.05 & abs(diffHiC_gr$log2FoldChange) > log2(1.5)]
  else if (grepl("HiC_pval0.05", id))
    diffHiC_gr_subset <- diffHiC_gr[diffHiC_gr$padj < 0.05]
  else if (grepl("HiC_fc1.5", id))
    diffHiC_gr_subset <- diffHiC_gr[abs(diffHiC_gr$log2FoldChange) > log2(1.5)]
  else if (grepl("HiC_fc2", id))
    diffHiC_gr_subset <- diffHiC_gr[abs(diffHiC_gr$log2FoldChange) > log2(2)]
  else if (grepl("CaptureC_all", id))
    diffHiC_gr_subset <- diffCaptureC_gr
  else if (grepl("CaptureC_fc1.5", id))
    diffHiC_gr_subset <- diffCaptureC_gr[abs(diffCaptureC_gr$log2FoldChange) > log2(1.5)]
  else if (grepl("CaptureC_fc2", id))
    diffHiC_gr_subset <- diffCaptureC_gr[abs(diffCaptureC_gr$log2FoldChange) > log2(2)]
  else if (grepl("CaptureC_reduced", id))
  {
    # extend by +/- 1 kb
    if (grepl("CaptureC_reduced_all", id))
      gr <- diffCaptureC_gr + 1e3
    else if (grepl("CaptureC_reduced_fc1.5", id))
      gr <- diffCaptureC_gr[abs(diffCaptureC_gr$log2FoldChange) > log2(1.5)] + 1e3
    else if (grepl("CaptureC_reduced_fc2", id))
      gr <- diffCaptureC_gr[abs(diffCaptureC_gr$log2FoldChange) > log2(2)] + 1e3
    else
      stop("unknown diffHiC_gr_subset; id: ", id)

    # temporary data table to assign new diffc_ids to reduced other ends
    gr_reduced_dt <- as.data.table(gr)[, as.data.table(reduce(GRanges(.SD))),
      by = c("baitID", "baitChr", "baitStart", "baitEnd", "baitName")]
    gr_reduced_dt[, diffc_id := seq_len(nrow(gr_reduced_dt))]

    # final GenomicRanges object
    dt <- as.data.table(gr)[, as.data.table(reduce(GRanges(.SD))),
      by = c("baitID", "gene_id", "baitChr", "baitStart", "baitEnd", "baitName")]
    dtm <- merge(dt, gr_reduced_dt,
      by = c("baitID", "baitChr", "baitStart", "baitEnd", "baitName", "seqnames", "start", "end", "width", "strand"), all.x = T)
    stopifnot(nrow(dt) == nrow(dtm))
    diffHiC_gr_subset <- GRanges(dtm)
  }
  else
    stop("unknown diffHiC_gr_subset; id: ", id)

  if (grepl("disruptedTADs", id))
  {
    halves <- fread("analysis/balancer/diffTAD_thresh12.5kb_breakpoints_halves.tab", header = T)
    halves_gr <- with(halves, GRanges(w_chrom, IRanges(w_start, w_end)))
    ov <- findOverlaps(diffHiC_gr_subset, halves_gr)
    diffHiC_gr_subset <- diffHiC_gr_subset[unique(queryHits(ov))]
  }

  TSS_subset_to_diffHiC_dt <- logFC_vs_dist(TSS_subset, diffHiC_gr_subset, bp_gr, mode = "all", extra_filtering = "bait", max_dist = 99e6)
  TSS_subset_to_diffHiC_dt[, level := TSS_subset$level[match(gene_id, TSS_subset$gene_id)]]
  print(summary(TSS_subset_to_diffHiC_dt$level))

  diffHiC_ov[[id]]$TSS_subset <- TSS_subset
  diffHiC_ov[[id]]$diffHiC_gr_subset <- diffHiC_gr_subset
  diffHiC_ov[[id]]$TSS_subset_to_diffHiC_dt <- TSS_subset_to_diffHiC_dt

  fc_to_diffHiC_subset <- list()
  fc_gr_list <- list()
  fc_to_overlap <- list()

  fc_gr_list[["all"]] <- diffHiC_gr_subset

  if (grepl("DE_all", id))
  {
    TSS_s <- GRanges(TSS_subset[!is.na(signf) & signf == "s", ])
  }
  else if (grepl("DE_fc1.5", id))
  {
    TSS_s <- GRanges(TSS_subset[!is.na(signf) & signf == "s" & abs(log2FoldChange) > log2(1.5), ])
  }
  else
    stop("unknown TSS_subset; id: ", id)

  if (grepl("activeenh", id))
  {
    anno_to_overlap <- list(
      TSS = TSS_gr,
      # TSS_active = TSS_active_gr,
      TSS_testable = GRanges(TSS_subset[!is.na(signf) & signf %in% c("s", "i"), ]),
      TSS_s = TSS_s,
      TSS_i = GRanges(TSS_subset[!is.na(signf) & signf == "i", ]),
      # TSS_inactive = TSS_inactive_gr,
      H3K27ac = H3K27ac_gr,
      eRNA_active = eRNA_active_gr,
      DHS_distal = DHS_gr[!DHS_gr %over% promoter_gr],
      H3K27ac_distal = H3K27ac_gr[!H3K27ac_gr %over% promoter_gr],
      eRNA_active_distal = eRNA_active_gr[!eRNA_active_gr %over% promoter_gr],
      DHS_active = DHS_gr[(DHS_gr %over% H3K27ac_gr) | (DHS_gr %over% eRNA_active_gr)],
      DHS_active_distal = DHS_gr[((DHS_gr %over% H3K27ac_gr) | (DHS_gr %over% eRNA_active_gr)) & !(DHS_gr %over% promoter_gr)]
    )
  }
  else
  {
    anno_to_overlap <- list(
      TSS = TSS_gr,
      # TSS_active = TSS_active_gr,
      TSS_testable = GRanges(TSS_subset[!is.na(signf) & signf %in% c("s", "i"), ]),
      TSS_s = TSS_s,
      TSS_i = GRanges(TSS_subset[!is.na(signf) & signf == "i", ]),
      # TSS_inactive = TSS_inactive_gr,
      DHS_distal = DHS_gr[!DHS_gr %over% promoter_gr],
      DHS_enhancer_distal = DHS_gr[(DHS_gr %over% enhancer_gr) & !(DHS_gr %over% promoter_gr)],
      DHS_insulator_distal = DHS_gr[(DHS_gr %over% insulator_gr) & !(DHS_gr %over% promoter_gr)]
      # promoter_active = promoter_active_gr,
      # promoter_inactive = promoter_inactive_gr,
      # DHS_James_distal = DHS_James_gr[!DHS_James_gr %over% promoter_gr],
      # insulator = insulator_gr,
      # CAD4 = CAD4_gr,
      # CAD4_distal = CAD4_gr[!CAD4_gr %over% promoter_gr],
      # mesoCRM = mesoCRM_gr,
      # mesoCRM_distal = mesoCRM_gr[!mesoCRM_gr %over% promoter_gr],
      # enhancer = enhancer_gr,
    )
  }

  for (a in names(anno_to_overlap))
  {
    if ("gene_id" %in% names(elementMetadata(diffHiC_gr_subset)))
    {
      # if we are using Capture-C data, rely on gene assignments there
      fc_to_diffHiC_subset[[a]] <- fc_gr_list[[a]] <- diffHiC_gr_subset[diffHiC_gr_subset$gene_id %in% anno_to_overlap[[a]]$gene_id]
    }
    else
    {
      fc_to_diffHiC_subset[[a]] <- feature_vs_dist(anno_to_overlap[[a]], diffHiC_gr_subset, bp_gr,
        mode = "all", extra_filtering = "bait", max_dist = 99e6)
      fc_gr_list[[a]] <- diffHiC_gr_subset[with(fc_to_diffHiC_subset[[a]], unique(center_id))]
    }
    fc_to_overlap[[a]] <- anno_to_overlap[[a]]
  }

  # TSS_to_a <- feature_vs_dist(TSS_gr, diffHiC_gr_subset, bp_gr, mode = "all", extra_filtering = "bait", max_dist = 99e6)
  # fc_gr_list[["all_distal"]] <- diffHiC_gr_subset[-with(TSS_to_a, unique(center_id))]

  diffHiC_ov[[id]]$fc_to_diffHiC_subset <- fc_to_diffHiC_subset
  diffHiC_ov[[id]]$fc_gr_list <- fc_gr_list
  diffHiC_ov[[id]]$fc_to_overlap <- fc_to_overlap
  # diffHiC_ov[[id]]$anno_by_fc <- overlap_to_shuffled(fc_gr_list, fc_to_diffHiC_subset, anno_to_overlap)
  diffHiC_ov[[id]]$fc_by_fc <- overlap_to_shuffled(fc_gr_list, fc_to_diffHiC_subset, fc_to_overlap)

  dir_to_diffHiC_subset <- list()
  dir_gr_list <- list()
  dir_to_overlap <- list()

  # dir_gr_list[["n"]] <- diffHiC_gr_subset[with(TSS_subset_to_diffHiC_dt, unique(center_id[!is.na(signf) & signf == "n"]))]
  # dir_to_overlap[["n"]] <- TSS_gr[TSS_gr$gene_id %in% with(TSS_subset, gene_id[!is.na(signf) & signf == "n"])]

  dir_to_diffHiC_subset[["iu"]] <- TSS_subset_to_diffHiC_dt[!is.na(signf) & signf == "i" & log2FoldChange > 0, ]
  dir_gr_list[["iu"]] <- diffHiC_gr_subset[unique(dir_to_diffHiC_subset[["iu"]]$center_id)]
  dir_to_overlap[["iu"]] <- TSS_gr[TSS_gr$gene_id %in% with(TSS_subset, gene_id[!is.na(signf) & signf == "i" & log2FoldChange > 0])]

  dir_to_diffHiC_subset[["id"]] <- TSS_subset_to_diffHiC_dt[!is.na(signf) & signf == "i" & log2FoldChange < 0, ]
  dir_gr_list[["id"]] <- diffHiC_gr_subset[unique(dir_to_diffHiC_subset[["id"]]$center_id)]
  dir_to_overlap[["id"]] <- TSS_gr[TSS_gr$gene_id %in% with(TSS_subset, gene_id[!is.na(signf) & signf == "i" & log2FoldChange < 0])]

  dir_to_diffHiC_subset[["su"]] <- TSS_subset_to_diffHiC_dt[!is.na(signf) & signf == "s" & log2FoldChange > 0, ]
  dir_gr_list[["su"]] <- diffHiC_gr_subset[unique(dir_to_diffHiC_subset[["su"]]$center_id)]
  dir_to_overlap[["su"]] <- TSS_gr[TSS_gr$gene_id %in% with(TSS_subset, gene_id[!is.na(signf) & signf == "s" & log2FoldChange > 0])]

  dir_to_diffHiC_subset[["sd"]] <- TSS_subset_to_diffHiC_dt[!is.na(signf) & signf == "s" & log2FoldChange < 0, ]
  dir_gr_list[["sd"]] <- diffHiC_gr_subset[unique(dir_to_diffHiC_subset[["sd"]]$center_id)]
  dir_to_overlap[["sd"]] <- TSS_gr[TSS_gr$gene_id %in% with(TSS_subset, gene_id[!is.na(signf) & signf == "s" & log2FoldChange < 0])]

  diffHiC_ov[[id]]$dir_by_dir <- overlap_to_shuffled(dir_gr_list, dir_to_diffHiC_subset, dir_to_overlap)

  diffHiC_ov[[id]]$gmap <- list()
  for (l in rev(lev))
  {
    if (grepl("HiC_", id))
      g <- TSS_subset[grepl(paste0("^", l), level), ]
    else if (grepl("CaptureC_", id))
      g <- TSS_subset[grepl(paste0("^", l), level) & gene_id %in% genes_diffCaptureC$gene_id, ]
    else
      stop("unknown diffHiC_gr_subset; id: ", id)

    diffHiC_ov[[id]]$gmap[[l]] <- extract_map_genes(g, diffHiC_gr_subset, bp_gr,
      mode = "all", fc_method = "most_significant", extra_filtering = "bait")
  }

  # # extend to matching using MatchIt
  # m.out <- matchit(signf == "s" ~ baseMean + FPKM, data = genes_tested[, c("gene_id", "baseMean", "signf", "FPKM")])
  # summary(m.out)
  # head(match.data(m.out))

  ov_list = diffHiC_ov[[id]]
  save(ov_list, file = paste0("analysis/balancer/diffHiC_overlap/", id, "_ov_list.Rdata"))
}
