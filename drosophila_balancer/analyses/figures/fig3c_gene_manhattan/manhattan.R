library(data.table)
library(ggplot2)
library(scales)
source("scripts/colors.R")
format_Mb   <- function(x) {paste(comma(x/1e6), "Mb")}


ase = fread(snakemake@input[["ase"]])
setkey(ase, gene_id)

genes = fread(snakemake@input[["exons"]])
colnames(genes) = c("chrom","start","end","gene_id")
setkey(genes, gene_id)
genes = genes[grepl('^FBgn0*', gene_id), .(chrom = chrom[1], pos = (min(start) + max(end))/2), by = gene_id]

ase = merge(ase, genes)
ase[, sign := ifelse(padj<0.05, ifelse(log2FoldChange > 0, "bal", "vrg"), "not")]



minuslog10 <- function(l) {
  l = -l
  l <- paste0("10^",l)
  # turn the 'e+' into plotmath format
  #l <- gsub("e", "%*%10^", l)
  # return this as an expression
  parse(text=l)
}


plt = ggplot(ase) + 
  aes(x = pos, y = - log10(padj), shape = sign, col = sign, alpha = sign) + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", col = "darkorange") +
  geom_segment(aes(y = 0, yend = - log10(padj), x = pos, xend = pos, col = sign), size = 0.2) + # Comment out to remove vertical lines!
  geom_point(size = 1.1) + 
  scale_alpha_manual(values = c(not = 0.5, bal = 1, vrg = 1)) +
  scale_shape_manual(values = c(not = 16, bal = 15, vrg = 17)) +
  scale_color_manual(values = c(not = "grey", 
                                bal = as.character(my_colors["balancer"]),
                                vrg = as.character(my_colors["wildtype"]) )) + 
  facet_wrap(~ chrom, nrow = 1, scales = "free_x") + 
  scale_x_continuous(breaks = pretty_breaks(3), labels = format_Mb, expand = c(0,5e5)) +
  scale_y_continuous(breaks = c(0,5,10,15), labels = minuslog10) +
  theme_minimal() + 
  theme(legend.position = "none", 
        panel.spacing = unit(0.1, "cm"), 
        panel.background = element_rect(fill = "white", colour = "grey60")) +
  xlab("Genomic position") +
  ylab("Adjusted p-value")
ggsave(snakemake@output[[1]], plot = plt, width = 14, height = 3)



plt = ggplot(ase) + 
  geom_segment(data = ase[log2FoldChange>=0], aes(y = 0, yend = - log10(padj), x = pos, xend = pos, col = sign), size = 0.2) + # Comment out to remove vertical lines!
  geom_segment(data = ase[log2FoldChange<0],  aes(y = 0, yend = + log10(padj), x = pos, xend = pos, col = sign), size = 0.2) + # Comment out to remove vertical lines!
  geom_point(data = ase[log2FoldChange>=0], aes(x = pos, y = - log10(padj), shape = sign, col = sign, alpha = sign), size = 1.1) + 
  geom_point(data = ase[log2FoldChange<0],  aes(x = pos, y = + log10(padj), shape = sign, col = sign, alpha = sign), size = 1.1) + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", col = "darkorange") +
  geom_hline(yintercept = +log10(0.05), linetype = "dashed", col = "darkorange") +
  scale_alpha_manual(values = c(not = 0.5, bal = 1, vrg = 1)) +
  scale_shape_manual(values = c(not = 16, bal = 15, vrg = 17)) +
  scale_color_manual(values = c(not = "grey", 
                                bal = as.character(my_colors["balancer"]),
                                vrg = as.character(my_colors["wildtype"]) )) + 
  facet_wrap(~ chrom, nrow = 1, scales = "free_x") + 
  scale_x_continuous(breaks = pretty_breaks(3), labels = format_Mb, expand = c(0,5e5)) +
  scale_y_continuous(breaks = c(-15,-10,-5,0,5,10,15), labels = function(x) {minuslog10(abs(x))}) +
  theme_minimal() + 
  theme(legend.position = "none", 
        panel.spacing = unit(0.1, "cm"), 
        panel.background = element_rect(fill = "white", colour = "grey60")) +
  xlab("Genomic position") +
  ylab("Adjusted p-value\nSeparated for up-/down-regulated genes")
ggsave(paste0(snakemake@output[[1]],".v2.pdf"), plot = plt, width = 14, height = 3)


plt = ggplot(ase) + 
  aes(x = pos, y = ..count.., fill = factor(sign, levels = c("not", "bal","vrg"), ordered = T),
      col = factor(sign, levels = c("not", "bal","vrg"), ordered = T)) + 
  geom_density(adjust = 0.3333, alpha = 1, position = "stack") +
  facet_wrap(~ chrom, nrow = 1, scales = "free_x") + 
  scale_x_continuous(breaks = pretty_breaks(3), labels = format_Mb, expand = c(0,1e6)) +
  scale_fill_manual(values = c(not = "grey", 
                                bal = as.character(my_colors["balancer"]),
                                vrg = as.character(my_colors["wildtype"]) )) +
  scale_color_manual(values = c(not = "grey", 
                                bal = as.character(my_colors["balancer"]),
                                vrg = as.character(my_colors["wildtype"]) )) + 
  theme_minimal() + 
  theme(legend.position = "none", 
        panel.spacing = unit(0.1, "cm"), 
        panel.background = element_rect(fill = "white", colour = "grey60"),
        axis.text.y = element_blank()) +
  xlab("Genomic position") +
  ylab("Relative density")
ggsave(paste0(snakemake@output[[1]],".hist.pdf"), plot = plt, width = 14, height = 2.5)
