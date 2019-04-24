library(data.table)
library(scales)
library(ggplot2)
source("scripts/colors.R")

gff = fread(snakemake@input[["gff"]])
gff <- gff[, .(chrom = substr(chrom[1],1,4)), by = gene_id]

d = fread(snakemake@input[["f1"]])
d <- merge(d, gff, by = "gene_id")
d <- merge(d, d[,.N,by=chrom], by = "chrom")
d$significant = factor(d$padj < 0.05, levels=c(T,F), labels = c("significant ASE", "no ASE"))

x <- d[,.(K = .N), by = .(chrom,significant,N)]

x <- x[, .(fraction = K/N, 
      lower_conf = binom.test(K,N,K/N)$conf.int[1],
      upper_conf = binom.test(K,N,K/N)$conf.int[2]),
  by = .(chrom,significant,K,N)]

x <- x[, chrom := paste0(chrom, " (", K, "/", N, ")")][]

p <- ggplot(x[significant == "significant ASE"]) + 
  aes(chrom, y = fraction) + 
  geom_bar(stat = "identity", fill = "grey", col = "black") +
  geom_errorbar(aes(ymin = lower_conf, ymax = upper_conf), width = 0.5) +
  ylab("ASE fraction") +
  xlab(NULL) + 
  coord_flip() +
  scale_y_continuous(labels = percent) +
  my_theme
ggsave(snakemake@output[["numbers"]], plot = p, width = 4, height = 2.5)


y1 = quantile(d[chrom != "chrX",]$log2FoldChange, probs = seq(0.01,0.99,0.01))
y2 = quantile(d[chrom == "chrX",]$log2FoldChange, probs = seq(0.01,0.99,0.01))
q <- ggplot(data.table(y1,y2)) + aes(y1,y2) + geom_line() + geom_point() +
  xlab("Log fold change quantiles (chr 2 & 3)") + ylab("Log fold change quantiles (chr X)") +
  my_theme

ggsave(snakemake@output[["qqplot"]], plot = q, width = 4, height = 4)