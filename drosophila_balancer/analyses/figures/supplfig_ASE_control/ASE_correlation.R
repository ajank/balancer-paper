library(data.table)
library(ggplot2)
library(cowplot)
source("scripts/DESeq_plots.R")

# Load ASE analyses to get the log2FoldChange (bal vs. vrg) of different data sets
single = rbind(cbind(fread(snakemake@input[["cyo_chr2"]]), chrom = "chr2"),
               cbind(fread(snakemake@input[["tm3_chr3"]]), chrom = "chr3"))
f1old  = rbind(cbind(fread(snakemake@input[["old_chr2"]]), chrom = "chr2"),
               cbind(fread(snakemake@input[["old_chr3"]]), chrom = "chr3"))
f1new  = rbind(cbind(fread(snakemake@input[["new_chr2"]]), chrom = "chr2"),
               cbind(fread(snakemake@input[["new_chr3"]]), chrom = "chr3"))

# Now also load interaction term analyses, to see which genes are significantly influenced
new_vs_old = rbind(cbind(fread(snakemake@input[["intf1_chr2"]]), chrom = "chr2"),
                   cbind(fread(snakemake@input[["intf1_chr3"]]), chrom = "chr3"))
si_vs_doub = rbind(cbind(fread(snakemake@input[["intn1_chr2"]]), chrom = "chr2"),
                   cbind(fread(snakemake@input[["intn1_chr3"]]), chrom = "chr3"))


# Plot preparations
ALPHA = 0.1
fac_labels = c(paste0("significant (FDR ",round(ALPHA*100),"%)"), "not significant")
fac_colors = c("dodgerblue", "#cecece")
names(fac_colors) = fac_labels

# Plot F1 new vs. F1 old
new_vs_old = new_vs_old[, .(gene_id, significant = factor(padj < ALPHA, levels = c(T,F), labels = fac_labels), int.log2fc = log2FoldChange)]
a = merge(f1old, f1new, by = "gene_id", suffixes = c(".old",".new"))
a = merge(a, new_vs_old, by = "gene_id")
message("       f1old: ", nrow(f1old))
message("       f1new: ", nrow(f1new))
message("  new_vs_old: ", nrow(new_vs_old))
message("intersection: ", nrow(a))

rsqu_a = round(cor(a$log2FoldChange.old, a$log2FoldChange.new),3)
p1 <- ggplot(a) +
  aes(log2FoldChange.old, log2FoldChange.new, color = significant, shape = chrom.old) +
  geom_point(alpha = 0.5) +
  xlab(expression(paste("Balancer : wild type ratio (log2) in ", F[1]^CyOTM3))) +
  ylab(expression(paste("Balancer : wild type ratio (log2) in ", F[1]^female))) +
  scale_color_manual(values = fac_colors, name = "ASE ratio difference") +
  my_theme + theme(legend.position = 'bottom') +
  ggtitle(paste("ASE ratios across two different samples")) +
  coord_cartesian(ylim = c(-10,10), xlim = c(-10,10)) +
  annotate("text", x=-9,  y=9.5, hjust=0, vjust=1, label = paste("italic(R) ^ 2 :", rsqu_a), parse = T, size = 3) +
  annotate("text", x=-9,  y=8,   hjust=0, vjust=1, label = paste0("#sign.: ", nrow(a[significant == fac_labels[1],]), "/", nrow(a)), size = 3)


legend <- g_legend(p1)
p1 <- p1 + theme(legend.position = "none")
p1



# Plot double balancer vs. single - balanced
si_vs_doub = si_vs_doub[, .(gene_id, significant = factor(padj < ALPHA, levels = c(T,F), labels = fac_labels), int.log2fc = log2FoldChange)]
b = merge(f1old, single, by = "gene_id", suffixes = c(".double",".single"))
b = merge(b, si_vs_doub, by = "gene_id")
message("  si_vs_doub: ", nrow(si_vs_doub))
message("      single: ", nrow(single))
message("      double: ", nrow(f1old))
message("intersection: ", nrow(b))

rsqu_b = round(cor(b$log2FoldChange.double, b$log2FoldChange.single),3)
p2 <- ggplot(b) +
  aes(log2FoldChange.double, log2FoldChange.single, color = significant, shape = chrom.double) +
  geom_point(alpha = 0.5) +
  xlab(expression(paste("Balancer : wild type ratio (log2) in ", F[1]^CyOTM3))) +
  ylab(expression(paste("Balancer : wild type ratio (log2) in ", F[1]^CyO, " resp. ", F[1]^TM3))) +
  scale_color_manual(values = fac_colors, name = "ASE ratio difference") +
  my_theme + theme(legend.position = 'none') +
  ggtitle(paste("ASE ratios in single- and double-balanced cross")) +
  coord_cartesian(ylim = c(-10,10), xlim = c(-10,10)) +
  annotate("text", x=-9,  y=9.5, hjust=0, vjust=1, label = paste("italic(R) ^ 2 :", rsqu_b), parse = T, size = 3) +
  annotate("text", x=-9,  y=8,   hjust=0, vjust=1, label = paste0("#sign.: ", nrow(b[significant == fac_labels[1],]), "/", nrow(b)), size = 3)

xx <- plot_grid(a = p1, b = p2, align = 'h')
final <- plot_grid(xx, legend, ncol = 1, rel_heights = c(8,1))
ggsave(snakemake@output[[1]], plot = final, width = 10, height = 5.5)


write.table(b, file = paste0(snakemake@output[[1]], ".data.txt"), quote=F, row.names=F, col.names=T, sep="\t")
