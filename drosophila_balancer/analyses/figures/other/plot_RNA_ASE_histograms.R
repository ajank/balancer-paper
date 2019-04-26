suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(scales))
suppressMessages(library(grid))
source("../scripts/colors.R")

# Min number of reads per gene
minCount = 20

gff <- unique(fread("../../SNV_annotation/reformatted.gff")[, .(gene_id, chrom)])
gff$chrom = substr(gff$chrom, 1,4)

# Only look at chrom 2 and 3
gff <- gff[chrom %in% paste0("chr", c("2","3")),]


# Read data from three samples (2x 6-8h, 1x 4-8h)
bal = NULL
vrg = NULL
for ( sample in c("N1_pool_6-8h","N1sex_pool_6-8h","N1_pool_4-8h")) {
  bal = rbind(bal, rbind(
      fread(paste0("../../readSeparation/counts/",sample,"_rep1.htseq-count-rev.alt.txt"))[!grepl('^__',V1), .(gene_id = V1, sample = sample, `balancer` = V2, replicate = "Replicate 1")],
      fread(paste0("../../readSeparation/counts/",sample,"_rep2.htseq-count-rev.alt.txt"))[!grepl('^__',V1), .(gene_id = V1, sample = sample, `balancer` = V2, replicate = "Replicate 2")]))
  vrg = rbind(vrg, rbind(
      fread(paste0("../../readSeparation/counts/",sample,"_rep1.htseq-count-rev.ref.txt"))[!grepl('^__',V1), .(gene_id = V1, sample = sample, `wild type` = V2, replicate = "Replicate 1")],
      fread(paste0("../../readSeparation/counts/",sample,"_rep2.htseq-count-rev.ref.txt"))[!grepl('^__',V1), .(gene_id = V1, sample = sample, `wild type` = V2, replicate = "Replicate 2")]))
}

m <- merge(bal, vrg, by=c("gene_id", "replicate", "sample"))[balancer + `wild type` > minCount,]

# This merge also removes all genes that are not on 2, 3, or X
m <- merge(m, gff, by = "gene_id")


#tr1 = grobTree(textGrob(paste0("Genes with at least ", minCount, " fragments (N =", nrow(m), ")"), x=0.98,  y=0.88, hjust=1, vjust=1))

p = ggplot(m) + aes(balancer/(balancer+`wild type`), y = ..count.., fill = sample) + 
      #geom_density(fill="black", alpha=0.5, adjust = 0.5) + 
      geom_histogram(binwidth = 0.01) +
      scale_x_continuous(labels = percent) +
      xlab(paste0("Balancer fraction (genes with >=", minCount, " fragments)")) + 
      my_theme +
  geom_vline(xintercept = 0.25, linetype = "dotted", color = "black")

p1 = p + facet_grid(replicate ~ sample) +
     geom_label(aes(x = 0.99, y = Inf, label = paste("N =",N)),
                data = m[, .N, by = .(replicate, sample)], fill = "white",
                hjust = 1, vjust =1, position = position_nudge(y = -10, x = 0)) +
    theme(legend.position = "None")
ggsave(p1, filename = "ASE_histograms.pdf", width=8, height=5)
