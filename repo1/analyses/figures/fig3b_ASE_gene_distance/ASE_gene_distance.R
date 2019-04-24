library(ggplot2)
library(data.table)
library(scales)
library(GenomicRanges)
source("scripts/colors.R")



ASE.embryo <- fread(snakemake@input[["ase"]]) # ASE.embryo <- fread("/g/korbel/shared/projects/drosophila_balancer/analyses/ase/deseq/DESeq.N1_6-8h.standardFormat.txt")
gff        <- fread(snakemake@input[["gff"]]) # gff        <- fread("/g/korbel/shared/projects/drosophila_balancer/analyses/SNV_annotation/reformatted.gff")


genes      <- gff[, .(chrom = chrom[1], start = min(start), end = max(end), strand = strand[1]), by = gene_id]
gene_pos     = genes[,chrom := substr(chrom,1,4)][order(chrom,start),]
gene_pos.ase = gene_pos[gene_id %in% ASE.embryo[padj<0.05,]$gene_id,]
gene_dist.ase = gene_pos.ase[, .(distance = diff(start)), by = chrom]

gene_pos.bg  = gene_pos[gene_id %in% ASE.embryo[padj>=0.05,]$gene_id,]
gene_dist.bg = NULL
for (i in 1:500) {
  rnd_select = sort(sample(1:nrow(gene_pos.bg), size = nrow(gene_pos.ase)))
  gene_dist.bg = rbind(gene_dist.bg,
                       gene_pos.bg[rnd_select,][, .(distance = diff(start), run = i), by = chrom])
}

#p1 <- ggplot(gene_dist.ase) +
#  geom_histogram(data = gene_dist.bg,  aes(x = distance, y = ..count../100, fill="non-ASE genes"), binwidth=0.2) +
#  geom_histogram(data = gene_dist.ase, aes(x = distance, y = ..count.., fill="ASE genes"), binwidth=0.2, alpha = 0.4) +
#  scale_x_log10(breaks = c(1e2,1e3,1e4,1e5,1e6), label=comma, limits = c(30,2e6)) +
#  ylab("Count") + xlab("Distance between genes (bp)") +
#  my_theme +
#  theme(legend.position = c(0.17, 0.83), legend.title = element_blank(),
#        legend.background = element_rect(fill="white", size=0.2)) +
#  scale_fill_manual(values = c(`ASE genes` = "darkorange", `non-ASE genes` = "#bbbbbb"))
# ggsave(snakemake@output[[1]], plot = p1, width = 6, height = 4)



conf.interval_bins <- ggplot(gene_dist.bg) + aes(x = distance, y = ..count..) + geom_histogram(binwidth = 0.2) + scale_x_log10() + facet_wrap(~run)
conf.interval_data = as.data.table(ggplot_build(conf.interval_bins)$data[[1]])
conf.interval_data = conf.interval_data[, .(mean = mean(count), q10 = quantile(count,0.05), q90 = quantile(count,0.95)), by = x]

p1 <- ggplot(conf.interval_data) +
  geom_bar(aes(10**x, mean, fill = "control genes"), stat = "identity", col = "#999999") +
  geom_errorbar(aes(10**x, ymin = q10, ymax = q90), col = "#999999", width = 0.1) +
  scale_x_log10(breaks = c(1e2,1e3,1e4,1e5,1e6), label=comma, limits = c(30,2e6)) +
  ylab("Count") + xlab("Distance between genes (bp)") +
  my_theme +
  theme(legend.position = c(0.17, 0.83), legend.title = element_blank(),
        legend.background = element_rect(fill="white", size=0.2)) +
  scale_fill_manual(values = c(`ASE genes` = "darkorange", `control genes` = "#bbbbbb")) +
  geom_histogram(data = gene_dist.ase, aes(distance, fill = "ASE genes"), size= 0, binwidth = 0.2, alpha = 0.4)
ggsave(snakemake@output[[1]], plot = p1, width = 6, height = 4)
