library(data.table)
library(ggplot2)
library(cowplot)
source("scripts/colors.R")

#load("/g/furlong/project/37_Capture-C/analysis/balancer_cap2/contacts_noARS_all_lowseqdepth/DESeq2_interactions.Rdata", verbose =T)
#load("/g/furlong/project/37_Capture-C/analysis/balancer_cap2/contacts_noARS_all/DESeq2_interactions.Rdata", verbose =T)
load(snakemake@input[[1]])

# 1) remove all non-unique view points
pmd <- pmd[!is.na(embryo.log2FoldChange), ]



interactions <- pmd[order(embryo.log2FoldChange),]
interactions <- interactions[,
              c("baitName", "ase", "contact") := 
              .(factor(interactions$baitName, levels = unique(interactions$baitName), ordered = T),
                ifelse(embryo.padj < 0.05, ifelse(embryo.log2FoldChange < 0, "w", "b") , "n"),
                factor(int.log2FoldChange < 0, levels = c(T,F), labels = c("w", "b"))   )][]
interactions <- interactions[,
              c("ase_N_genes") :=
              .(length(unique(baitName))), by = ase][]
int.ase_labels = c(w = "higher expression in wild type",
                   b = "higher expression in balancer",
                   n = "no significant ASE")
int.contact_labels = c(w = "stronger in wild type",
                       b = "stronger in balancer")
int.ase_colors = c(w = as.character(my_colors["wildtype"]), 
                   b = as.character(my_colors["balancer"]),
                   n = "darkgrey")


# 1) box plot - actually no longer a boxplot
color_zero_point = interactions[, (0 - min(embryo.log2FoldChange)) / (max(embryo.log2FoldChange) - min(embryo.log2FoldChange))]
p1 <- ggplot(interactions) + 
  aes(x = baitName, y = int.log2FoldChange, col = embryo.log2FoldChange, shape = embryo.padj < 0.05) + 
  geom_line(aes(group = baitName), size = 0.3) +
  geom_point(size = 1.2) + 
  scale_color_gradientn(name = "DE strength (log2 fold change)", 
                        colors = c("red", "darkseagreen3", "dodgerblue"), 
                        values = c(0,color_zero_point - 0.11, color_zero_point, color_zero_point + 0.08, 1) ) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_shape_manual(values = c(1,19), name = "Significant ASE gene") + 
  guides(shape = guide_legend(override.aes = list(size=3))) +
  scale_y_continuous(breaks = seq(-10,10,2)) +
  my_theme + 
  xlab(NULL) +
  ylab("Strength of diff. interactions (log2 fold change)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
        legend.position = "bottom", 
        legend.box.just = "left",
        panel.grid.major.x = element_blank()) +
  ggtitle("Strength of significant differential interactions (ordered by ASE)")
ggsave(p1, filename = snakemake@output[["boxplot"]], width = 13, height = 5.5)


# 2) scatter plot of all diff. interactions
p2 <- ggplot(interactions) + 
  aes(embryo.log2FoldChange, int.log2FoldChange, alpha = embryo.padj < 0.05, col = as.character(as.integer(baitName)%%8)) + 
  geom_point(size = 1.5) +
  geom_line(aes(group = baitName), size = 0.2) +
  guides(color = F) + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  scale_color_brewer(type = "qual", palette = 6) +
  xlab("Diff. gene expression (log2 fold change)") +
  ylab("Strength of diff .interaction (log2 fold change)") +
  my_theme + 
  theme(legend.position = "bottom", legend.box.just = "left") +
  ggtitle("Significant differential interactions") +
  scale_alpha_manual(values = c(T = 0.9, F = 0.2), name = NULL)
ggsave(snakemake@output[["scatter"]], plot = p2, width = 8, height = 4)


# 3) Histogram
p3 <- ggplot(interactions[, .(.N, ase = ase[1]), by = baitName]) + 
  geom_histogram(aes(N, fill = ase), binwidth = 1) + 
  my_theme + 
  xlab("Number of diff. interactions per gene") + 
  ylab("Count") + 
  scale_fill_manual(values = int.ase_colors,
                    labels = int.ase_labels,
                    name = "ASE") +
  theme(legend.position = c(0.95, 0.95), legend.justification = c(1,1), legend.background = element_rect(size = 0.2))
ggsave(snakemake@output[["histogram"]], plot = p3, width = 4, height = 3)



# 4) Violin plots
num_int_per_gene <- interactions[,. (.N, ase = ase[1]), by = baitName]
p4 <- ggplot(num_int_per_gene) + 
  aes(ase, N, col = ase)  + 
  geom_violin() + 
  geom_jitter(alpha = 0.6) + 
  coord_flip() +
  guides(color = F) + 
  scale_x_discrete(labels = int.ase_labels) +
  scale_color_manual(values = int.ase_colors) +
  ylab("Diff. interactions per gene") +
  xlab(NULL)

wilcox.test(num_int_per_gene[ase=="w", N],
            num_int_per_gene[ase=="b", N])
wilcox.test(num_int_per_gene[ase=="n", N],
            num_int_per_gene[ase=="b", N])
wilcox.test(num_int_per_gene[ase=="n", N],
            num_int_per_gene[ase=="w", N])

ggsave(snakemake@output[["num_int_violin"]], plot = p4, width = 9, height = 3.5)


# 5) statistics per interaction 1

# add number of genes to labels
label.dt = merge(data.table(label = int.ase_labels, ase = names(int.ase_labels)), interactions[, .N, by = ase], by = "ase")
new_labels        = str_wrap(label.dt[, paste0(label, " (", N, " genes)")], width = 20)
names(new_labels) = label.dt[,ase]


p5 <- ggplot(interactions) +
  aes(ase, col = ase, fill = ase, alpha = contact) + 
  geom_bar(position = position_dodge(), size = 1) + 
  coord_flip() + 
  ylab("Count") +
  ggtitle("Number of diff. interactions stratified by genes") + 
  xlab(NULL) + 
  scale_x_discrete(labels   = new_labels) + 
  scale_color_manual(values = int.ase_colors) +
  scale_fill_manual( values = int.ase_colors) +
  scale_alpha_manual(values = c(w = 0.4, b = 1),
                     labels = int.contact_labels,
                     name = "Sign. diff. interactions") +
  my_theme +
  guides(color = F, fill = F, alpha = guide_legend(override.aes=list(fill="#444444",colour="#444444"))) + 
  theme(legend.position = "bottom") +
  geom_text(data = interactions[, .N, by = .(ase,contact)], 
            position = position_dodge(width = 1), 
            aes(x = ase, 
                y = ifelse(N < 5, N+1, ifelse(N<30, N-1,N-3)), 
                label = paste("N=",N), 
                group = contact, 
                hjust = ifelse(N < 5, 0, 1)),
            color = "black",
            alpha = 1,
            size = 3)
p5
ggsave(snakemake@output[["num_diff_int"]], plot = p5, width = 7, height = 3.5)



# 6) statistics per interaction 2
p6 <- ggplot(interactions) +
  aes(ase, abs(int.log2FoldChange), col = ase, fill = ase, alpha = contact) + 
  geom_violin(position = position_dodge(), size = 1) + 
  coord_flip() + 
  ylab("Diff. interaction ratio (log2)") +
  ggtitle("Number of diff. interactions stratified by genes") + 
  xlab(NULL) + 
  scale_x_discrete(labels   = new_labels) + 
  scale_color_manual(values = int.ase_colors) +
  scale_fill_manual( values = int.ase_colors) +
  scale_alpha_manual(values = c(w = 0.4, b = 1),
                     labels = int.contact_labels,
                     name = "Sign. diff. interactions") +
  my_theme +
  guides(color = F, fill = F, alpha = guide_legend(override.aes=list(fill="#444444",colour="#444444"))) + 
  theme(legend.position = "bottom")
p6
ggsave(snakemake@output[["strength_diff_int"]], plot = p6, width = 7, height = 3.5)




# 7) statistics per interaction 3
p7 <- ggplot(interactions) +
  aes(ase, col = ase, fill = ase, alpha = contact) + 
  geom_bar(position = "fill", size = 1) + 
  coord_flip() + 
  ylab("Count") +
  ggtitle("Number of diff. interactions stratified by genes") + 
  xlab(NULL) + 
  scale_x_discrete(labels   = new_labels) + 
  scale_color_manual(values = int.ase_colors) +
  scale_fill_manual( values = int.ase_colors) +
  scale_alpha_manual(values = c(w = 0.4, b = 1),
                     labels = int.contact_labels,
                     name = "Sign. diff. interactions") +
  my_theme +
  guides(color = F, fill = F, alpha = guide_legend(override.aes=list(fill="#444444",colour="#444444"))) + 
  theme(legend.position = "bottom") +
  geom_text(data = interactions[, .N, by = .(ase,contact)], 
            aes(x = ase, 
                y = ifelse(contact == "w", 0.02, 0.98), 
                label = paste("N=",N), 
                group = contact, 
                hjust = ifelse(contact == "w", 0.02, 0.98)),
            color = "black",
            alpha = 1,
            size = 3)
p7
ggsave(snakemake@output[["num_diff_relative"]], plot = p7, width = 7, height = 3.5)




