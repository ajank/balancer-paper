options <- commandArgs(trailingOnly = TRUE)

if (length(options) < 1)
  stop("Usage:  Rscript plot_DESeq2_stats.R outputDir")

outputDir <- options[1]

# outputDir <- "analysis/balancer_cap2/contacts_noARS_all_lowseqdepth"

require(Chicago)
require(GenomicRanges)
options(warn = 1)

#
#  read CHiCAGO data
#

load(paste0(outputDir, "/DESeq2_interactions.Rdata")) # baits, pm, pmd

# FIXME this should be done earlier
pm <- pm[grepl("^chr[23][LR]", otherEndChr), ]
pmd <- pmd[grepl("^chr[23][LR]", otherEndChr), ]

message(length(unique(pm$baitName)), " baits within the interactions called by CHiCAGO")
message(nrow(baits), " baits in Capture-C experiment")
message(nrow(pmd), " differential interactions in ", outputDir)

# Create a custom color scale
library(RColorBrewer)
mycol <- brewer.pal(9, "Paired")[c(2, 4, 6)]

#
#  check what fraction of restriction fragments is affected by allele-specific genomic variation
#

require(data.table)
annotated_DpnII <- fread("analysis/digest_DpnII_dm6_balancer_allele_specific_variation.tab", header = T)

stats <- NULL
add_stats <- function(prefix, sel)
{
  chr23 <- grepl("^chr[23][LR]", annotated_DpnII$chrom)
  pmbait <- pm$baitID %in% annotated_DpnII$frag_id[sel]
  pmdbait <- pmd$baitID %in% annotated_DpnII$frag_id[sel]
  pmoe <- pm$otherEndID %in% annotated_DpnII$frag_id[sel]
  pmdoe <- pmd$otherEndID %in% annotated_DpnII$frag_id[sel]
  col <- data.frame(c(
    sum(chr23 & sel) / sum(chr23),
    sum(pmbait) / length(pmbait),
    sum(pmdbait) / length(pmdbait),
    sum(pmoe) / length(pmoe),
    sum(pmdoe) / length(pmdoe),
    sum(pmbait | pmoe) / length(pmbait),
    sum(pmdbait | pmdoe) / length(pmdbait)
  ))
  rownames(col) <- c(
    "Among all DpnII sites (chr2+3 only)",
    "Among baits of tested interactions",
    "Among baits of differential interactions",
    "Among other ends of tested interactions",
    "Among other ends of differential interactions",
    "Among either ends of tested interactions",
    "Among either ends of differential interactions"
  )
  colnames(col) <- prefix

  if (is.null(stats))
    stats <<- col
  else
    stats <<- cbind(stats, col)
}
add_stats("Allele-specific CNV", annotated_DpnII$affected_CNV)
add_stats("Allele-specific RS", annotated_DpnII$affected_RS)
add_stats("Allele-specific CNV or RS", annotated_DpnII$affected_CNV | annotated_DpnII$affected_RS)
print(stats)

#
#  carry over bait annotations
#

baits[, gene_class_with_DE := ifelse(gene_class == "one gene",
  ifelse(DE_class == "one DE gene", "DE gene", "non-DE gene"),
  "multiple genes")]
baits[, gene_class_with_DE := factor(gene_class_with_DE, c("DE gene", "non-DE gene", "multiple genes"))]

baits[, gene_class_with_DE2 := ifelse(gene_class == "one gene",
  ifelse(DE_class == "one DE gene", ifelse(embryo.log2FoldChange > 0, "gene up in balancer", "gene up in wild-type"), "non-DE gene"),
  "multiple genes")]
baits[, gene_class_with_DE2 := factor(gene_class_with_DE2, c("gene up in wild-type", "gene up in balancer", "non-DE gene", "multiple genes"))]

pm[, DE_class := baits$DE_class[match(baitID, baits$baitID)]]
pmd[, DE_class := baits$DE_class[match(baitID, baits$baitID)]]

pm[, gene_class_with_DE := baits$gene_class_with_DE[match(baitID, baits$baitID)]]
pmd[, gene_class_with_DE := baits$gene_class_with_DE[match(baitID, baits$baitID)]]

pm[, gene_class_with_DE2 := baits$gene_class_with_DE2[match(baitID, baits$baitID)]]
pmd[, gene_class_with_DE2 := baits$gene_class_with_DE2[match(baitID, baits$baitID)]]

#
#  check all the baits, calculate how many interactions they have
#

sum.All <- pm[, list(.N), by = "baitID"]
baits$sum.All <- 0L
baits$sum.All[match(sum.All$baitID, baits$baitID)] <- sum.All$N

sum.diff <- pmd[, list(oeCAD4 = sum(oeCAD4), absLog2FoldChange = sum(abs(int.log2FoldChange)), .N), by = "baitID"]
baits$sum.diff <- 0L
baits$sum.diff[match(sum.diff$baitID, baits$baitID)] <- sum.diff$N
baits$sum.diff.oeCAD4 <- 0L
baits$sum.diff.oeCAD4[match(sum.diff$baitID, baits$baitID)] <- sum.diff$oeCAD4
baits$sum.diff.absLog2FoldChange <- 0
baits$sum.diff.absLog2FoldChange[match(sum.diff$baitID, baits$baitID)] <- sum.diff$absLog2FoldChange

sum.diff.BAL <- pmd[int.log2FoldChange > 0, list(.N), by = "baitID"]
baits$sum.diff.BAL <- 0L
baits$sum.diff.BAL[match(sum.diff.BAL$baitID, baits$baitID)] <- sum.diff.BAL$N

extr.diff <- pmd[, list(dist = abs(distSign)[which.min(int.padj)], score = score[which.min(int.padj)], int.log2FoldChange = int.log2FoldChange[which.min(int.padj)], int.padj = min(int.padj)), by = "baitID"]
baits$extr.dist <- NA
baits$extr.score <- NA
baits$extr.int.log2FoldChange <- NA
baits$extr.int.padj <- NA
baits$extr.dist[match(extr.diff$baitID, baits$baitID)] <- extr.diff$dist
baits$extr.score[match(extr.diff$baitID, baits$baitID)] <- extr.diff$score
baits$extr.int.log2FoldChange[match(extr.diff$baitID, baits$baitID)] <- extr.diff$int.log2FoldChange
baits$extr.int.padj[match(extr.diff$baitID, baits$baitID)] <- extr.diff$int.padj


sum.nab.All <- pm[across_breakpoint == F, list(.N), by = "baitID"]
baits$sum.nab.All <- 0L
baits$sum.nab.All[match(sum.nab.All$baitID, baits$baitID)] <- sum.nab.All$N

sum.nab.diff <- pmd[across_breakpoint == F, list(oeCAD4 = sum(oeCAD4), absLog2FoldChange = sum(abs(int.log2FoldChange)), .N), by = "baitID"]
baits$sum.nab.diff <- 0L
baits$sum.nab.diff[match(sum.nab.diff$baitID, baits$baitID)] <- sum.nab.diff$N
baits$sum.nab.diff.oeCAD4 <- 0L
baits$sum.nab.diff.oeCAD4[match(sum.nab.diff$baitID, baits$baitID)] <- sum.nab.diff$oeCAD4
baits$sum.nab.diff.absLog2FoldChange <- 0
baits$sum.nab.diff.absLog2FoldChange[match(sum.nab.diff$baitID, baits$baitID)] <- sum.nab.diff$absLog2FoldChange

sum.nab.diff.BAL <- pmd[across_breakpoint == F & int.log2FoldChange > 0, list(.N), by = "baitID"]
baits$sum.nab.diff.BAL <- 0L
baits$sum.nab.diff.BAL[match(sum.nab.diff.BAL$baitID, baits$baitID)] <- sum.nab.diff.BAL$N

extr.nab.diff <- pmd[across_breakpoint == F, list(dist = abs(distSign)[which.min(int.padj)], score = score[which.min(int.padj)], int.log2FoldChange = int.log2FoldChange[which.min(int.padj)], int.padj = min(int.padj)), by = "baitID"]
baits$extr.nab.dist <- NA
baits$extr.nab.score <- NA
baits$extr.nab.int.log2FoldChange <- NA
baits$extr.nab.int.padj <- NA
baits$extr.nab.dist[match(extr.nab.diff$baitID, baits$baitID)] <- extr.nab.diff$dist
baits$extr.nab.score[match(extr.nab.diff$baitID, baits$baitID)] <- extr.nab.diff$score
baits$extr.nab.int.log2FoldChange[match(extr.nab.diff$baitID, baits$baitID)] <- extr.nab.diff$int.log2FoldChange
baits$extr.nab.int.padj[match(extr.nab.diff$baitID, baits$baitID)] <- extr.nab.diff$int.padj


#
#  now do the plotting
#

require(GenomicRanges)
require(ggplot2)


pdf(paste0(outputDir, "/plot_bait_stats.pdf"), width = 8, height = 3.75)

mycol <- brewer.pal(9, "Set1")[c(4, 5, 9)]

# p <- ggplot(baits, aes(sum.All, fill = DE_class)) +
#   labs(x = "Number of tested interactions", y = "Viewpoint count") +
#   scale_fill_brewer(palette = "Paired") +
#   geom_histogram(binwidth = 200, boundary = 0)
# print(p)

p <- ggplot(baits, aes(sum.diff, fill = gene_class_with_DE)) +
  labs(x = "Number of differential interactions", y = "Viewpoint count") +
  # scale_x_continuous(limits = c(-2.5, 52.5), oob = scales::squish, expand = c(0, 0)) +
  scale_fill_manual(values = mycol, name = "Viewpoint class") +
  geom_histogram(binwidth = 1, boundary = 0)
print(p)

# p <- ggplot(baits, aes(extr.dist / 1e3, fill = DE_class)) +
#   labs(x = "Distance to the most signif. differential interaction (kb)", y = "Viewpoint count") +
#   scale_x_continuous(limits = c(-0.05, 1.05) * 200, oob = scales::squish, expand = c(0, 0)) +
#   scale_fill_brewer(palette = "Paired") +
#   geom_histogram(binwidth = 5, boundary = 0)
# print(p)

dev.off()


pdf(paste0(outputDir, "/plot_interaction_stats.pdf"), width = 8, height = 3.75)

# p <- ggplot(pm, aes(abs(distSign) / 1e3, fill = DE_class)) +
#   labs(x = "Distance to tested interactions (kb)", y = "Interaction count") +
#   scale_x_continuous(limits = c(-0.05, 1.05) * 0.5e3, expand = c(0, 0)) +
#   scale_fill_brewer(palette = "Paired") +
#   geom_histogram(binwidth = 10, boundary = 0)
# print(p)

p <- ggplot(pmd, aes(abs(distSign) / 1e3, fill = gene_class_with_DE)) +
  labs(x = "Distance to differential interactions (kb)", y = "Differential interaction count") +
  scale_x_continuous(limits = c(-0.05, 1.05) * 0.5e3, expand = c(0, 0)) +
  scale_fill_manual(values = mycol, name = "Viewpoint class") +
  geom_histogram(binwidth = 10, boundary = 0)
print(p)

p <- ggplot(pmd, aes(pmin(abs(distSign) / 1e3, 100 + 1), fill = gene_class_with_DE)) +
  labs(x = "Distance to differential interactions (kb)", y = "Differential interaction count") +
  scale_x_continuous(limits = c(-0.05, 1.05) * 100, expand = c(0, 0)) +
  scale_fill_manual(values = mycol, name = "Viewpoint class") +
  geom_histogram(binwidth = 2, boundary = 0)
print(p)

dev.off()


pdf(paste0(outputDir, "/plot_bait_xyplots.pdf"), width = 8, height = 3.75)

g <- ggplot(baits, aes(embryo.log2FoldChange, sum.All)) +
  geom_point(aes(color = DE_class), alpha = 0.5) +
  scale_y_continuous(limits = c(-0.05, 1.05) * 2e3, oob = scales::squish, expand = c(0, 0)) +
  labs(x = "log2FoldChange (gene expression, embryo)", y = "Number of tested interactions") +
  scale_color_manual(values = mycol)
print(g)

g <- ggplot(baits, aes(-log10(embryo.padj), sum.All)) +
  geom_point(aes(color = DE_class), alpha = 0.5) +
  scale_y_continuous(limits = c(-0.05, 1.05) * 2e3, oob = scales::squish, expand = c(0, 0)) +
  labs(x = "-log10(padj) (gene expression, embryo)", y = "Number of tested interactions") +
  scale_color_manual(values = mycol)
print(g)

g <- ggplot(baits, aes(embryo.log2FoldChange, sum.diff)) +
  geom_point(aes(color = DE_class), alpha = 0.5) +
  scale_y_continuous(limits = c(-0.05, 1.05) * 10, oob = scales::squish, expand = c(0, 0)) +
  labs(x = "log2FoldChange (gene expression, embryo)", y = "Number of differential interactions") +
  scale_color_manual(values = mycol)
print(g)

g <- ggplot(baits, aes(-log10(embryo.padj), sum.diff)) +
  geom_point(aes(color = DE_class), alpha = 0.5) +
  scale_y_continuous(limits = c(-0.05, 1.05) * 10, oob = scales::squish, expand = c(0, 0)) +
  labs(x = "-log10(padj) (gene expression, embryo)", y = "Number of differential interactions") +
  scale_color_manual(values = mycol)
print(g)



g <- ggplot(baits, aes(embryo.log2FoldChange, extr.dist / 1e3)) +
  geom_point(aes(color = DE_class), alpha = 0.5) +
  scale_y_continuous(limits = c(-12.5, 262.5), oob = scales::squish, expand = c(0, 0)) +
  labs(x = "log2FoldChange (gene expression, embryo)", y = "Distance to the most signif. differential interaction") +
  scale_color_manual(values = mycol)
print(g)

g <- ggplot(baits, aes(-log10(embryo.padj), extr.dist / 1e3)) +
  geom_point(aes(color = DE_class), alpha = 0.5) +
  scale_y_continuous(limits = c(-12.5, 262.5), oob = scales::squish, expand = c(0, 0)) +
  labs(x = "-log10(padj) (gene expression, embryo)", y = "Distance to the most signif. differential interaction") +
  scale_color_manual(values = mycol)
print(g)

g <- ggplot(baits, aes(embryo.log2FoldChange, extr.int.log2FoldChange)) +
  geom_point(aes(color = DE_class), alpha = 0.5) +
  # scale_y_continuous(limits = c(-0.05, 1.05) * 100 / 15, oob = scales::squish, expand = c(0, 0)) +
  labs(x = "log2FoldChange (gene expression, embryo)", y = "log2FoldChange of most signif. differential interaction") +
  scale_color_manual(values = mycol)
print(g)

g <- ggplot(baits, aes(-log10(embryo.padj), extr.int.log2FoldChange)) +
  geom_point(aes(color = DE_class), alpha = 0.5) +
  # scale_y_continuous(limits = c(-0.05, 1.05) * 100 / 15, oob = scales::squish, expand = c(0, 0)) +
  labs(x = "-log10(padj) (gene expression, embryo)", y = "log2FoldChange of most signif. differential interaction") +
  scale_color_manual(values = mycol)
print(g)

g <- ggplot(baits, aes(embryo.log2FoldChange, -log10(extr.int.padj))) +
  geom_point(aes(color = DE_class), alpha = 0.5) +
  scale_y_continuous(limits = c(-1, 21), oob = scales::squish, expand = c(0, 0)) +
  labs(x = "log2FoldChange (gene expression, embryo)", y = "-log10(padj) of most signif. differential interaction") +
  scale_color_manual(values = mycol)
print(g)

g <- ggplot(baits, aes(-log10(embryo.padj), -log10(extr.int.padj))) +
  geom_point(aes(color = DE_class), alpha = 0.5) +
  scale_y_continuous(limits = c(-1, 21), oob = scales::squish, expand = c(0, 0)) +
  labs(x = "-log10(padj) (gene expression, embryo)", y = "-log10(padj) of most signif. differential interactions") +
  scale_color_manual(values = mycol)
print(g)



g <- ggplot(baits, aes(embryo.log2FoldChange, sum.nab.All)) +
  geom_point(aes(color = DE_class), alpha = 0.5) +
  scale_y_continuous(limits = c(-0.05, 1.05) * 2e3, oob = scales::squish, expand = c(0, 0)) +
  labs(x = "log2FoldChange (gene expression, embryo)", y = "Number of tested interactions\nnot crossing a breakpoint") +
  scale_color_manual(values = mycol)
print(g)

g <- ggplot(baits, aes(embryo.log2FoldChange, sum.nab.diff)) +
  geom_point(aes(color = DE_class), alpha = 0.5) +
  scale_y_continuous(limits = c(-0.05, 1.05) * 10, oob = scales::squish, expand = c(0, 0)) +
  labs(x = "log2FoldChange (gene expression, embryo)", y = "Number of differential interactions\nnot crossing a breakpoint") +
  scale_color_manual(values = mycol)
print(g)

g <- ggplot(baits, aes(embryo.log2FoldChange, sum.diff.BAL / sum.diff)) +
  geom_point(aes(color = DE_class), alpha = 0.5) +
  scale_y_continuous(limits = c(-0.05, 1.05), oob = scales::squish, expand = c(0, 0)) +
  labs(x = "log2FoldChange (gene expression, embryo)", y = "BAL-enriched/all differential interactions") +
  scale_color_manual(values = mycol)
print(g)

g <- ggplot(baits, aes(-log10(embryo.padj), sum.diff.BAL / sum.diff)) +
  geom_point(aes(color = DE_class), alpha = 0.5) +
  scale_y_continuous(limits = c(-0.05, 1.05), oob = scales::squish, expand = c(0, 0)) +
  labs(x = "-log10(padj) (gene expression, embryo)", y = "BAL-enriched/all differential interactions") +
  scale_color_manual(values = mycol)
print(g)


g <- ggplot(baits, aes(embryo.log2FoldChange, sum.diff.absLog2FoldChange)) +
  geom_point(aes(color = DE_class), alpha = 0.5) +
  scale_y_continuous(limits = c(-0.05, 1.05) * 20, oob = scales::squish, expand = c(0, 0)) +
  labs(x = "log2FoldChange (gene expression, embryo)", y = "Sum of absolute interaction log2FoldChange") +
  scale_color_manual(values = mycol)
print(g)

g <- ggplot(baits, aes(-log10(embryo.padj), sum.diff.absLog2FoldChange)) +
  geom_point(aes(color = DE_class), alpha = 0.5) +
  scale_y_continuous(limits = c(-0.05, 1.05) * 20, oob = scales::squish, expand = c(0, 0)) +
  labs(x = "-log10(padj) (gene expression, embryo)", y = "Sum of absolute interaction log2FoldChange") +
  scale_color_manual(values = mycol)
print(g)


g <- ggplot(baits, aes(embryo.log2FoldChange, sum.nab.diff.BAL / sum.nab.diff)) +
  geom_point(aes(color = DE_class), alpha = 0.5) +
  scale_y_continuous(limits = c(-0.05, 1.05), oob = scales::squish, expand = c(0, 0)) +
  labs(x = "log2FoldChange (gene expression, embryo)", y = "BAL-enriched/all differential interactions\nnot crossing a breakpoint") +
  scale_color_manual(values = mycol)
print(g)

g <- ggplot(baits, aes(-log10(embryo.padj), sum.nab.diff.BAL / sum.nab.diff)) +
  geom_point(aes(color = DE_class), alpha = 0.5) +
  scale_y_continuous(limits = c(-0.05, 1.05), oob = scales::squish, expand = c(0, 0)) +
  labs(x = "-log10(padj) (gene expression, embryo)", y = "BAL-enriched/all differential interactions\nnot crossing a breakpoint") +
  scale_color_manual(values = mycol)
print(g)


g <- ggplot(baits, aes(embryo.log2FoldChange, extr.nab.dist / 1e3)) +
  geom_point(aes(color = DE_class), alpha = 0.5) +
  scale_y_continuous(limits = c(-12.5, 262.5), oob = scales::squish, expand = c(0, 0)) +
  labs(x = "log2FoldChange (gene expression, embryo)", y = "Distance to the most signif. differential interaction\nnot crossing a breakpoint") +
  scale_color_manual(values = mycol)
print(g)

g <- ggplot(baits, aes(-log10(embryo.padj), extr.nab.dist / 1e3)) +
  geom_point(aes(color = DE_class), alpha = 0.5) +
  scale_y_continuous(limits = c(-12.5, 262.5), oob = scales::squish, expand = c(0, 0)) +
  labs(x = "-log10(padj) (gene expression, embryo)", y = "Distance to the most signif. differential interaction\nnot crossing a breakpoint") +
  scale_color_manual(values = mycol)
print(g)

g <- ggplot(baits, aes(embryo.log2FoldChange, extr.nab.int.log2FoldChange)) +
  geom_point(aes(color = DE_class), alpha = 0.5) +
  # scale_y_continuous(limits = c(-0.05, 1.05) * 100 / 15, oob = scales::squish, expand = c(0, 0)) +
  labs(x = "log2FoldChange (gene expression, embryo)", y = "log2FoldChange of most signif. differential interaction\nnot crossing a breakpoint") +
  scale_color_manual(values = mycol)
print(g)

g <- ggplot(baits, aes(-log10(embryo.padj), extr.nab.int.log2FoldChange)) +
  geom_point(aes(color = DE_class), alpha = 0.5) +
  # scale_y_continuous(limits = c(-0.05, 1.05) * 100 / 15, oob = scales::squish, expand = c(0, 0)) +
  labs(x = "-log10(padj) (gene expression, embryo)", y = "log2FoldChange of most signif. differential interaction\nnot crossing a breakpoint") +
  scale_color_manual(values = mycol)
print(g)

g <- ggplot(baits, aes(embryo.log2FoldChange, -log10(extr.nab.int.padj))) +
  geom_point(aes(color = DE_class), alpha = 0.5) +
  scale_y_continuous(limits = c(-1, 21), oob = scales::squish, expand = c(0, 0)) +
  labs(x = "log2FoldChange (gene expression, embryo)", y = "-log10(padj) of most signif. differential interaction\nnot crossing a breakpoint") +
  scale_color_manual(values = mycol)
print(g)

g <- ggplot(baits, aes(-log10(embryo.padj), -log10(extr.nab.int.padj))) +
  geom_point(aes(color = DE_class), alpha = 0.5) +
  scale_y_continuous(limits = c(-1, 21), oob = scales::squish, expand = c(0, 0)) +
  labs(x = "-log10(padj) (gene expression, embryo)", y = "-log10(padj) of most signif. differential interaction\nnot crossing a breakpoint") +
  scale_color_manual(values = mycol)
print(g)


g <- ggplot(baits, aes(embryo.log2FoldChange, sum.nab.diff.absLog2FoldChange)) +
  geom_point(aes(color = DE_class), alpha = 0.5) +
  scale_y_continuous(limits = c(-0.05, 1.05) * 20, oob = scales::squish, expand = c(0, 0)) +
  labs(x = "log2FoldChange (gene expression, embryo)", y = "Sum of absolute interaction log2FoldChange\nof interactions not crossing a breakpoint") +
  scale_color_manual(values = mycol)
print(g)

g <- ggplot(baits, aes(-log10(embryo.padj), sum.nab.diff.absLog2FoldChange)) +
  geom_point(aes(color = DE_class), alpha = 0.5) +
  scale_y_continuous(limits = c(-0.05, 1.05) * 20, oob = scales::squish, expand = c(0, 0)) +
  labs(x = "-log10(padj) (gene expression, embryo)", y = "Sum of absolute interaction log2FoldChange\nof interactions not crossing a breakpoint") +
  scale_color_manual(values = mycol)
print(g)


dev.off()


# pdf(paste0(outputDir, "/plot_interaction_stats.pdf"), width = 6, height = 4.5)

# p <- ggplot(pmd, aes(int.log2FoldChange, embryo.log2FoldChange)) +
#   geom_point(aes(color = across_breakpoint), alpha = 0.3) +
#   coord_fixed() +
#   # geom_vline(xintercept = 5, lty = 2) + 
#   # geom_hline(yintercept = 5, lty = 2) + 
#   # ggtitle(paste("Randomly sampled 1,000 out of", prettyNum(nrow(pm), big.mark = ","), "interactions")) +
#   labs(x = "log2FoldChange (differential Capture-C interaction)", y = "log2FoldChange (gene expression, embryo)") +
#   scale_color_brewer(palette = "Set1", direction = -1) +
#   theme(legend.position = "bottom")
# print(p)

# p <- ggplot(pmd[across_breakpoint == F, ], aes(int.log2FoldChange, embryo.log2FoldChange)) +
#   geom_point(aes(color = DE_class), alpha = 0.5) +
#   coord_fixed() +
#   # geom_vline(xintercept = 5, lty = 2) + 
#   # geom_hline(yintercept = 5, lty = 2) + 
#   # ggtitle(paste("Randomly sampled 1,000 out of", prettyNum(nrow(pm), big.mark = ","), "interactions")) +
#   labs(x = "log2FoldChange (differential Capture-C interaction)", y = "log2FoldChange (gene expression, embryo)") +
#   scale_color_brewer(palette = "Paired") +
#   theme(legend.position = "bottom")
# print(p)

# p <- ggplot(pmd[across_breakpoint == F, ], aes(score, embryo.log2FoldChange)) +
#   geom_point(aes(color = DE_class), alpha = 0.5) +
#   scale_x_continuous(limits = c(0, 20), oob = scales::squish, expand = c(0, 0)) +
#   geom_vline(xintercept = 3, lty = 2) + 
#   # geom_hline(yintercept = 5, lty = 2) + 
#   # ggtitle(paste("Randomly sampled 1,000 out of", prettyNum(nrow(pm), big.mark = ","), "interactions")) +
#   labs(x = "tested interaction score", y = "log2FoldChange (gene expression, embryo)") +
#   scale_color_brewer(palette = "Paired") +
#   theme(legend.position = "bottom")
# print(p)

# dev.off()

write.table(baits, file = paste0(outputDir, "/bait_stats.tsv"), sep = "\t", row.names = F, quote = F)

#
#  various statistics
#

print(summary(baits))

# are the interactions of the same viewpoint changing in the same direction?
m.log2FoldChange <- outer(pmd$int.log2FoldChange, pmd$int.log2FoldChange) >= 0
m.sameBaitID <- outer(pmd$baitID, pmd$baitID, "==")
print(wilcox.test(table(as.vector(m.log2FoldChange), as.vector(m.sameBaitID))))

# # are the interactions of the same viewpoint changing in the same direction?
# m.log2FoldChange <- outer(pmd[across_breakpoint == F, ]$int.log2FoldChange, pmd[across_breakpoint == F, ]$int.log2FoldChange) >= 0
# m.sameBaitID <- outer(pmd[across_breakpoint == F, ]$baitID, pmd[across_breakpoint == F, ]$baitID, "==")
# print(wilcox.test(table(as.vector(m.log2FoldChange), as.vector(m.sameBaitID))))

fraction_BAL <- function(pmd, shuffle = F)
{
  if (shuffle)
    pmd$baitID <- sample(pmd$baitID)
  dt <- pmd[, list(sum(int.log2FoldChange > 0) / .N), by = "baitID"]
  return(dt$V1)
}

fraction_BAL.pmd <- fraction_BAL(pmd)
fraction_BAL.control <- unlist(lapply(1:1000, function(i) fraction_BAL(pmd, shuffle = T)))
print(wilcox.test(fraction_BAL.pmd, fraction_BAL.control))

# fraction_BAL.pmd <- fraction_BAL(pmd[across_breakpoint == F, ])
# fraction_BAL.control <- unlist(lapply(1:1000, function(i) fraction_BAL(pmd[across_breakpoint == F, ], shuffle = T)))
# print(wilcox.test(fraction_BAL.pmd, fraction_BAL.control))
