suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(scales))
suppressMessages(library(grid))
source("scripts/colors.R")

# Min number of reads per gene
minCount = 20

gff <- unique(fread(snakemake@input[["gff"]])[, .(gene_id, chrom)])
gff$chrom = substr(gff$chrom, 1,4)
gff <- gff[chrom %in% paste0("chr", c("2","3","X")),]


bal1 = fread(snakemake@input[["bal1"]])[!grepl('^__',V1), .(gene_id = V1, `balancer` = V2, replicate = "Replicate 1")]
bal2 = fread(snakemake@input[["bal2"]])[!grepl('^__',V1), .(gene_id = V1, `balancer` = V2, replicate = "Replicate 2")]
bal  = rbind(bal1, bal2)
vrg1 = fread(snakemake@input[["vrg1"]])[!grepl('^__',V1), .(gene_id = V1, `wild type` = V2, replicate = "Replicate 1")]
vrg2 = fread(snakemake@input[["vrg2"]])[!grepl('^__',V1), .(gene_id = V1, `wild type` = V2, replicate = "Replicate 2")]
vrg  = rbind(vrg1, vrg2)



m <- merge(bal, vrg, by=c("gene_id", "replicate"))[balancer + `wild type` > minCount,]
# This merge also removes all genes that are not on 2, 3, or X
m <- merge(m, gff, by = "gene_id")


#tr1 = grobTree(textGrob(paste0("Genes with at least ", minCount, " fragments (N =", nrow(m), ")"), x=0.98,  y=0.88, hjust=1, vjust=1))

p = ggplot(m) + aes(balancer/(balancer+`wild type`), y = ..count..) + 
      #geom_density(fill="black", alpha=0.5, adjust = 0.5) + 
      geom_histogram(binwidth = 0.01, fill = "black", alpha = 0.5) +
      scale_x_continuous(labels = percent) +
      ggtitle(paste0("RNA-seq of ", snakemake@params[["title"]])) +
      xlab(paste0("Balancer fraction (genes with >=", minCount, " fragments)")) + 
      my_theme 

p1 = p + facet_grid(. ~ replicate) +
     geom_label(aes(x = 0.99, y = Inf, label = paste("N =",N)),
                data = m[, .N, by = replicate], 
                hjust = 1, vjust =1, position = position_nudge(y = -10, x = 0))
ggsave(p1, filename = snakemake@output[["hist"]], width=6, height=4)

p2 = p + facet_grid(chrom ~ replicate, scales = "free_y") +
     geom_label(aes(x = 0.99, y = Inf, label = paste("N =",N)),
                data = m[, .N, by = .(replicate, chrom)], 
                hjust = 1, vjust =1, position = position_nudge(y = -0.1, x = 0))
ggsave(p2, filename = snakemake@output[["chroms"]], width=6, height=7)


