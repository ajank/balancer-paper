library(data.table)
library(ggplot2)
library(cowplot)

theme_set(theme_cowplot(font_size = 11)) # reduce default font size
ts <- theme_get()$plot.subtitle
ts$hjust <- 0.5
theme_update(plot.subtitle = ts) # , legend.title = theme_get()$legend.text

g_legend <- function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    legend
}

my_colors = c(del       = "deepskyblue2",
              dup       = "#e41a1c",        # 228,26,28
              wildtype  = "#4daf4a",
              balancer  = "#377eb8",
              common    = "#999999",
              hetero    = "darkorange",
              error     = "firebrick1")
my_theme =  theme_minimal() + theme(strip.background = element_rect(), plot.title = element_text(hjust = 0.5,face = "bold"))

path_prefix <- "/g/korbel/shared/projects/drosophila_balancer"
input <- c(
  cyo_chr2   = paste0(path_prefix, "/analyses/ase/deseq/DE_control.cyo.chrom2.standardFormat.txt"),
  tm3_chr3   = paste0(path_prefix, "/analyses/ase/deseq/DE_control.tm3.chrom3.standardFormat.txt"),
  old_chr2   = paste0(path_prefix, "/analyses/ase/deseq/DE_control.oldf1.chrom2.standardFormat.txt"),
  old_chr3   = paste0(path_prefix, "/analyses/ase/deseq/DE_control.oldf1.chrom3.standardFormat.txt"),
  new_chr2   = paste0(path_prefix, "/analyses/ase/deseq/DE_control.newf1.chrom2.standardFormat.txt"),
  new_chr3   = paste0(path_prefix, "/analyses/ase/deseq/DE_control.newf1.chrom3.standardFormat.txt"),
  intf1_chr2 = paste0(path_prefix, "/analyses/ase/deseq/DE_control.old_vs_new.chrom2.standardFormat.txt"),
  intf1_chr3 = paste0(path_prefix, "/analyses/ase/deseq/DE_control.old_vs_new.chrom3.standardFormat.txt"),
  intn1_chr2 = paste0(path_prefix, "/analyses/ase/deseq/DE_control.cyo_vs_double.chrom2.standardFormat.txt"),
  intn1_chr3 = paste0(path_prefix, "/analyses/ase/deseq/DE_control.tm3_vs_double.chrom3.standardFormat.txt")
)

# Load ASE analyses to get the log2FoldChange (bal vs. vrg) of different data sets
single = rbind(cbind(fread(input[["cyo_chr2"]]), chrom = "chr2"),
               cbind(fread(input[["tm3_chr3"]]), chrom = "chr3"))
f1old  = rbind(cbind(fread(input[["old_chr2"]]), chrom = "chr2"),
               cbind(fread(input[["old_chr3"]]), chrom = "chr3"))
f1new  = rbind(cbind(fread(input[["new_chr2"]]), chrom = "chr2"),
               cbind(fread(input[["new_chr3"]]), chrom = "chr3"))

# Now also load interaction term analyses, to see which genes are significantly influenced
new_vs_old = rbind(cbind(fread(input[["intf1_chr2"]]), chrom = "chr2"),
                   cbind(fread(input[["intf1_chr3"]]), chrom = "chr3"))
si_vs_doub = rbind(cbind(fread(input[["intn1_chr2"]]), chrom = "chr2"),
                   cbind(fread(input[["intn1_chr3"]]), chrom = "chr3"))


# Plot preparations
ALPHA = 0.05
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

r_a = round(cor(a$log2FoldChange.old, a$log2FoldChange.new),3)
p1 <- ggplot(a) +
  aes(log2FoldChange.old, log2FoldChange.new, color = significant, shape = chrom.old) +
  geom_point(alpha = 0.5) +
  xlab(expression(paste("Gene expression\n", log[2] * " fold change in ", F[1]^"CyO,TM3"))) +
  ylab(expression(paste("Gene expression\n", log[2] * " fold change in ", F[1]^"CyO,TM3", " female"))) +
  scale_shape(name = "Chromosome") +
  scale_color_manual(values = fac_colors, name = "Effect of sample on DE") +
  my_theme + theme(legend.position = 'bottom') +
  # my_theme + theme(legend.position = 'bottom', legend.box = "vertical") +
  # guides(color = guide_legend(ncol = 1)) +
  ggtitle(paste("DE fold changes across two different samples")) +
  coord_fixed(ylim = c(-10,10), xlim = c(-10,10)) +
  annotate("text", x=-9,  y=9.5, hjust=0, vjust=1, label = paste("italic(r) ==", r_a), parse = T, size = 4) +
  annotate("text", x=-9,  y=8,   hjust=0, vjust=1, label = paste0("significant: ", nrow(a[significant == fac_labels[1],]), "/", format(nrow(a), big.mark = ","), " genes"), size = 4)

pdf("analysis/balancer/plot_DE_correlation.pdf", width = 10, height = 4.5)
# legend <- g_legend(p1)
# p1 <- p1 + theme(legend.position = "none")
print(p1)



# Plot double balancer vs. single - balanced
si_vs_doub = si_vs_doub[, .(gene_id, significant = factor(padj < ALPHA, levels = c(T,F), labels = fac_labels), int.log2fc = log2FoldChange)]
b = merge(f1old, single, by = "gene_id", suffixes = c(".double",".single"))
b = merge(b, si_vs_doub, by = "gene_id")
message("  si_vs_doub: ", nrow(si_vs_doub))
message("      single: ", nrow(single))
message("      double: ", nrow(f1old))
message("intersection: ", nrow(b))

r_b = round(cor(b$log2FoldChange.double, b$log2FoldChange.single),3)
p2 <- ggplot(b) +
  aes(log2FoldChange.double, log2FoldChange.single, color = significant, shape = chrom.double) +
  geom_point(alpha = 0.5) +
  xlab(expression(paste("Gene expression\n", log[2] * " fold change in ", F[1]^"CyO,TM3"))) +
  ylab(expression(paste("Gene expression\n", log[2] * " fold change in ", F[1]^CyO, " or ", F[1]^TM3, ", accordingly"))) +
  scale_shape(name = "Chromosome") +
  scale_color_manual(values = fac_colors, name = "Effect of genotype on DE") +
  my_theme + theme(legend.position = 'bottom') +
  # my_theme + theme(legend.position = 'bottom', legend.box = "vertical") +
  # guides(color = guide_legend(ncol = 1)) +
  ggtitle(paste("DE fold changes in single\u00ad and double\u00adbalanced cross")) +
  coord_fixed(ylim = c(-10,10), xlim = c(-10,10)) +
  annotate("text", x=-9,  y=9.5, hjust=0, vjust=1, label = paste("italic(r) ==", r_b), parse = T, size = 4) +
  annotate("text", x=-9,  y=8,   hjust=0, vjust=1, label = paste0("significant: ", nrow(b[significant == fac_labels[1],]), "/", format(nrow(b), big.mark = ","), " genes"), size = 4)
print(p2)

# xx <- plot_grid(a = p1, b = p2, align = 'h')
# final <- plot_grid(xx, legend, ncol = 1, rel_heights = c(8,1))
# ggsave("analysis/balancer/plot_DE_correlation_merged.pdf", plot = final, width = 10, height = 5.5)

# write.table(b, file = paste0(snakemake@output[[1]], ".data.txt"), quote=F, row.names=F, col.names=T, sep="\t")


r1 <- fread("analysis/balancer/DESeq.between_haplotypes.N1mat_Rep1_6-8h.standardFormat.txt", header = T)
r2 <- fread("analysis/balancer/DESeq.between_haplotypes.N1mat_Rep2_6-8h.standardFormat.txt", header = T)
b = merge(r1, r2, by = "gene_id", suffixes = c(".r1",".r2"))
b[, significant := "not significant"]

r_b = round(cor(b$log2FoldChange.r1, b$log2FoldChange.r2),3)
s <- format(nrow(b), big.mark = ",")
p2 <- ggplot(b) +
  aes(log2FoldChange.r1, log2FoldChange.r2, color = significant) +
  geom_point(alpha = 0.5) +
  xlab(expression(paste("Gene expression\n", log[2] * " fold change in replicate 1"))) +
  ylab(expression(paste("Gene expression\n", log[2] * " fold change in replicate 2"))) +
  scale_color_manual(values = fac_colors, name = NULL) +
  my_theme + theme(legend.position = 'bottom') +
  # my_theme + theme(legend.position = 'bottom', legend.box = "vertical") +
  # guides(color = guide_legend(ncol = 1)) +
  ggtitle(paste("DE fold changes across replicates")) +
  coord_fixed(ylim = c(-10,10), xlim = c(-10,10)) +
  annotate("text", x=-9,  y=9.5, hjust=0, vjust=1, label = paste("italic(r) == ' ' *", r_b), parse = T, size = 4) +
  annotate("text", x=-9,  y=8,   hjust=0, vjust=1, label = paste0("italic(n) == ' ", format(nrow(b), big.mark = ","), " ' * genes"), parse = T, size = 4) +
  NULL
print(p2)


r1 <- fread("analysis/balancer/DESeq.between_replicates.N1mat_balancer_6-8h.standardFormat.txt", header = T)
r2 <- fread("analysis/balancer/DESeq.between_replicates.N1mat_virginizer_6-8h.standardFormat.txt", header = T)
b = merge(r1, r2, by = "gene_id", suffixes = c(".r1",".r2"))
b[, significant := "not significant"]

r_b = round(cor(b$log2FoldChange.r1, b$log2FoldChange.r2),3)
p2 <- ggplot(b) +
  aes(log2FoldChange.r1, log2FoldChange.r2, color = significant) +
  geom_point(alpha = 0.5) +
  xlab(expression(paste("Gene expression\n", log[2] * " fold change in balancer"))) +
  ylab(expression(paste("Gene expression\n", log[2] * " fold change in wild\u00adtype"))) +
  scale_color_manual(values = fac_colors, name = NULL) +
  my_theme + theme(legend.position = 'bottom') +
  # my_theme + theme(legend.position = 'bottom', legend.box = "vertical") +
  # guides(color = guide_legend(ncol = 1)) +
  ggtitle(paste("Replicate 2/replicate 1 fold changes across haplotypes")) +
  coord_fixed(ylim = c(-10,10), xlim = c(-10,10)) +
  annotate("text", x=-9,  y=9.5, hjust=0, vjust=1, label = paste("italic(r) == ' ' *", r_b), parse = T, size = 4) +
  annotate("text", x=-9,  y=8,   hjust=0, vjust=1, label = paste0("italic(n) == ' ", format(nrow(b), big.mark = ","), " ' * genes"), parse = T, size = 4) +
  NULL
print(p2)

dev.off()
