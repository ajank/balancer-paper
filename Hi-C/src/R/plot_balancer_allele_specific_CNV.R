library(cowplot)
library(dplyr)
library(data.table)
library(GenomicRanges)
library(scales)

source("src/R/functions_balancer_genes.R")
source("src/R/functions_balancer_annotations.R")

theme_set(theme_cowplot(font_size = 11)) # reduce default font size
ts <- theme_get()$plot.subtitle
ts$hjust <- 0.5
theme_update(plot.subtitle = ts) # , legend.title = theme_get()$legend.text
theme_update(strip.text.x = element_text(size = 11), strip.text.y = element_text(size = 11))


# read exon overlaps from file
ovl <- fread("analysis/balancer/genes_allele_specific_CNV.tab", header = T, sep = "\t")

## Add information on completely deleted/duplicated genes
ovl_full_with_subjectHits <- as.data.table(findOverlaps(GRanges(genes), CNV_unreduced_gr, type = "within"))
ovl_full_with_subjectHits[, gene_id := genes$gene_id[queryHits]]
ovl_full_with_subjectHits[, CNV_type := CNV_unreduced_gr$type[subjectHits]]
ovl_full <- ovl_full_with_subjectHits[, .(gene_id, CNV_type)]
# Remove double entries
ovl_full <- unique(ovl_full)
## move the full overlaps out of ovl
ovl <- ovl[!paste(gene_id, CNV_type) %in% with(ovl_full, paste(gene_id, CNV_type))]

## Add information on genes with deleted/duplicated TSS
ovl_tss <- as.data.table(findOverlaps(TSS_gr, CNV_unreduced_gr, type = "within"))
ovl_tss[, gene_id := TSS_gr$gene_id[queryHits]]
ovl_tss[, CNV_type := CNV_unreduced_gr$type[subjectHits]]
ovl_tss <- ovl_tss[, .(gene_id, CNV_type)]
# Remove double entries
ovl_tss <- unique(ovl_tss)
## move the TSS overlaps out of ovl
ovl <- ovl[!paste(gene_id, CNV_type) %in% with(ovl_tss, paste(gene_id, CNV_type))]

## Add information of genes split by breakpoint
genes_split_by_breakpoint <- fread("analysis/balancer/genes_split_by_breakpoint.tab", header = F)$V1
ovl_split <- data.table(gene_id = genes_split_by_breakpoint, CNV_type = "split_by_breakpoint")

## Rename CNV type for visualization
ovl[, CNV_type_label := paste0(CNV_type, ifelse(overlap_width >= 100, "_100+bp", "_1-99bp"))]
ovl_full[, CNV_type_label := paste0(CNV_type, "_full")]
ovl_tss[, CNV_type_label := paste0(CNV_type, " (TSS)")]
ovl_split[, CNV_type_label := CNV_type]


## Merge with ASE information:

go <- merge(genes, ovl, by="gene_id", all.x=T)
# go$sign = factor(go$padj < 0.05, c(F,T), c("insignificant", "significant"))

pdf("analysis/balancer/genes_vs_CNV.pdf", width = 4.5, height = 3)

p <- ggplot(go[!is.na(CNV_type) & signf != "n",]) + 
  aes(x = CNV_type, y = log2FoldChange, fill = CNV_type, color = "dummy") + 
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~ signf, labeller = as_labeller(gene.labels)) + 
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  geom_jitter(alpha = 0.3) +
  scale_fill_manual(name = "Deletions", values = RColorBrewer::brewer.pal(4, "Paired")[c(1, 3, 2, 4)],
    labels = c(DEL_bal = "balancer", DEL_vrg = "wild\u00adtype", DUP_bal = "balancer", DUP_vrg = "wild\u00adtype")) +
  scale_color_manual(name = "Duplications", values = "black") +
  # background_grid(major = "y", minor = "y") +
  geom_hline(yintercept=c(-1, 1), linetype="dotted") +
  xlab(NULL) +
  ylim(c(-10, 10)) +
  ylab(expression(paste("Gene expression ", log[2] * " fold change"))) +
  guides(fill = guide_legend(order = 1), color = guide_legend(order = 2))
print(p)

dev.off()


message("\nFull overlap:")
go <- merge(genes, ovl_full_with_subjectHits, by = "gene_id")
for (ct in unique(go$CNV_type))
  with(go[signf != "n" & CNV_type == ct, ], message(ct, ": ", length(unique(gene_id)), " testable genes, ", length(unique(subjectHits)), " CNVs"))

go <- merge(genes, ovl_split, by="gene_id")
message("\nSplit by breakpoint: ")
ov <- findOverlaps(GRanges(go), bp_gr)
message(length(unique(go$gene_id)), " genes, ", length(unique(subjectHits(ov))), " breakpoints")
ov <- findOverlaps(GRanges(go[signf != "n"]), bp_gr)
message(length(unique(go[signf != "n"]$gene_id)), " testable genes, ", length(unique(subjectHits(ov))), " breakpoints")


## Main paper figure (full duplications and split by breakpoint):

pdf("analysis/balancer/genes_vs_DUP_and_split.pdf", width = 3.5, height = 2.75)

go <- merge(genes, ovl_full, by="gene_id")
# go[, CNV_general_type := grepl("split", CNV_type)]

p <- ggplot(go[signf != "n",]) + 
  aes(x = CNV_type_label, y = log2FoldChange, fill = CNV_type) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  # facet_wrap(~ CNV_general_type, scales = "free") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_jitter(alpha = 0.3) +
  scale_fill_manual(values = c(DUP_vrg = "#4daf4a", DUP_bal = "#377eb8", split_by_breakpoint = "#984ea3")) +
  scale_x_discrete(labels = c(DUP_vrg_full = "001", DUP_bal_full = "002", split_by_breakpoint = "003")) +
  scale_y_continuous(breaks = pretty_breaks(7)) +
  background_grid(major = "y", minor = "y") +
  # geom_segment(data = data.table(x = c(1, 2), yi = c(1, -1)), aes(x = x - 0.5, xend = x + 0.5, y = yi, yend = yi, fill = NULL), linetype="dashed") +
  geom_hline(yintercept=c(-1, 1), linetype="dashed") +
  xlab(NULL) +
  coord_fixed(ratio = 0.75) +
  # ylim(c(-10, 10)) +
  ylab(expression(paste("Gene expression ", log[2] * " fold change"))) +
  # guides(fill = guide_legend(order = 1), color = guide_legend(order = 2))
  guides(fill = F)
print(p)

r <- range(go$log2FoldChange)

go <- merge(genes, ovl_split, by="gene_id")
# go[, CNV_general_type := grepl("split", CNV_type)]

p <- ggplot(go[signf != "n",]) + 
  aes(x = CNV_type_label, y = log2FoldChange, fill = CNV_type) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  # facet_wrap(~ CNV_general_type, scales = "free") + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_jitter(alpha = 0.3) +
  scale_fill_manual(values = c(DUP_vrg = "#4daf4a", DUP_bal = "#377eb8", split_by_breakpoint = "#984ea3")) +
  scale_x_discrete(labels = c(DUP_vrg_full = "001", DUP_bal_full = "002", split_by_breakpoint = "003")) +
  scale_y_continuous(breaks = pretty_breaks(7)) +
  background_grid(major = "y", minor = "y") +
  geom_hline(yintercept=c(-1, 1), linetype="dashed") +
  xlab(NULL) +
  coord_fixed(ratio = 0.75 * (max(r) - min(r)) / (max(go$log2FoldChange) - min(go$log2FoldChange))) +
  ylab(expression(paste("Gene expression ", log[2] * " fold change"))) +
  # guides(fill = guide_legend(order = 1), color = guide_legend(order = 2))
  guides(fill = F)
print(p)

dev.off()


## Supplementary paper figure (deletions):

ovl_any <- as.data.table(findOverlaps(GRanges(genes), CNV_unreduced_gr))
ovl_any[, gene_id := genes$gene_id[queryHits]]
ovl_any[, CNV_type := CNV_unreduced_gr$type[subjectHits]]
ovl_width <- width(pintersect(GRanges(genes)[ovl_any$queryHits], CNV_unreduced_gr[ovl_any$subjectHits]))
ovl_any[, overlap_width := ovl_width]
ovl_any <- ovl_any[, .(gene_id, CNV_type, overlap_width)]
# Remove double entries
ovl_any <- unique(ovl_any)

ovl_with_width <- merge(ovl_any, with(genes, data.table(gene_id, gene_width = end - start + 1L)), all.x = T)
# === >50% of the gene ===
ovl_most <- ovl_with_width[overlap_width > gene_width / 2, ]
ovl_most[, CNV_type_label := paste0(CNV_type, " (most of the gene)")]

# === The 5’ most 50% of the gene  (i.e. 50% from the TSS) ===
gr <- GRanges(genes)
start(gr) <- ifelse(strand(gr) == "+", start(gr), (start(gr) + end(gr)) / 2)
end(gr) <- ifelse(strand(gr) == "+", (start(gr) + end(gr)) / 2, end(gr))
ovl_most_5prim_half <- as.data.table(findOverlaps(gr, CNV_unreduced_gr, type = "within"))
ovl_most_5prim_half[, gene_id := genes$gene_id[queryHits]]
ovl_most_5prim_half[, CNV_type := CNV_unreduced_gr$type[subjectHits]]
ovl_most_5prim_half <- ovl_most_5prim_half[, .(gene_id, CNV_type)]
# Remove double entries
ovl_most_5prim_half <- unique(ovl_most_5prim_half)
ovl_most_5prim_half[, CNV_type_label := paste0(CNV_type, " (5' half of the gene)")]

# >1kb from the 5’ end (TSS)
gr <- promoters(TSS_gr, upstream = 0, downstream = 1000)
ovl_tss_1kb <- as.data.table(findOverlaps(gr, CNV_unreduced_gr, type = "within"))
ovl_tss_1kb[, gene_id := genes$gene_id[queryHits]]
ovl_tss_1kb[, CNV_type := CNV_unreduced_gr$type[subjectHits]]
ovl_tss_1kb <- ovl_tss_1kb[, .(gene_id, CNV_type)]
# Remove double entries
ovl_tss_1kb <- unique(ovl_tss_1kb)
ovl_tss_1kb[, CNV_type_label := paste0(CNV_type, " (at least 1 kb from TSS)")]


pdf("analysis/balancer/genes_vs_DEL.pdf", width = 3, height = 4)

go <- merge(genes, ovl, by = "gene_id")
go <- go[signf != "n" & grepl("DEL", CNV_type), ]
go <- merge(go, go[, list(label = paste0(CNV_type_label, " (", format(.N, big.mark = ","), ")")), by = "CNV_type_label"], by = "CNV_type_label")

p <- ggplot(go) + 
  aes(x = label, y = log2FoldChange, fill = CNV_type) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_jitter(alpha = 0.3) +
  scale_fill_manual(values = c(DEL_vrg = "#4daf4a", DEL_bal = "#377eb8", split_by_breakpoint = "#984ea3")) +
  scale_y_continuous(breaks = pretty_breaks(7)) +
  background_grid(major = "y", minor = "y") +
  geom_hline(yintercept=c(-1, 1), linetype="dashed") +
  xlab(NULL) +
  ylab(expression(paste("Gene expression ", log[2] * " fold change"))) +
  guides(fill = F)
print(p)

dev.off()


pdf("analysis/balancer/genes_vs_DUP_partial.pdf", width = 3, height = 4)

go <- merge(genes, ovl, by = "gene_id")
go <- go[signf != "n" & grepl("DUP", CNV_type), ]
go <- merge(go, go[, list(label = paste0(CNV_type_label, " (", format(.N, big.mark = ","), ")")), by = "CNV_type_label"], by = "CNV_type_label")

p <- ggplot(go) + 
  aes(x = label, y = log2FoldChange, fill = CNV_type) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_jitter(alpha = 0.3) +
  scale_fill_manual(values = c(DUP_vrg = "#4daf4a", DUP_bal = "#377eb8", split_by_breakpoint = "#984ea3")) +
  scale_y_continuous(breaks = pretty_breaks(7)) +
  background_grid(major = "y", minor = "y") +
  geom_hline(yintercept=c(-1, 1), linetype="dashed") +
  xlab(NULL) +
  ylab(expression(paste("Gene expression ", log[2] * " fold change"))) +
  guides(fill = F)
print(p)

dev.off()


pdf("analysis/balancer/genes_vs_DEL_and_DUP_partial.pdf", width = 4, height = 4)

go <- merge(genes, ovl, by = "gene_id")
go <- go[signf != "n" & grepl("DEL|DUP", CNV_type), ]
go[, CNV_type := sub("(DUP_.*)_.*", "\\1", CNV_type)]
go <- merge(go, go[, list(label = paste0(CNV_type_label, " (", format(.N, big.mark = ","), ")")), by = "CNV_type_label"], by = "CNV_type_label")

p <- ggplot(go) + 
  aes(x = label, y = log2FoldChange, fill = CNV_type) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  geom_jitter(alpha = 0.3) +
  scale_fill_manual(values = c(DEL_vrg = "#4daf4a", DEL_bal = "#377eb8", DUP_vrg = "#4daf4a", DUP_bal = "#377eb8", split_by_breakpoint = "#984ea3")) +
  scale_y_continuous(breaks = pretty_breaks(7)) +
  background_grid(major = "y", minor = "y") +
  geom_hline(yintercept=c(-1, 1), linetype="dashed") +
  xlab(NULL) +
  ylab(expression(paste("Gene expression ", log[2] * " fold change"))) +
  guides(fill = F)
print(p)

dev.off()
