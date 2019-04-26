options(warn = 1)

require(cowplot)
require(ggplot2)
require(rcartocolor)
require(scales) # comma

source("src/R/functions_balancer_genes.R")
source("src/R/functions_balancer_annotations.R")

theme_set(theme_cowplot(font_size = 11)) # reduce default font size
ts <- theme_get()$plot.subtitle
ts$hjust <- 0.5
theme_update(plot.subtitle = ts) # , legend.title = theme_get()$legend.text
theme_update(strip.text.x = element_text(size = 11), strip.text.y = element_text(size = 11))


# 1) remove all non-unique view points
interactions <- diffCaptureC_dt[!is.na(embryo.log2FoldChange), ]
setkey(interactions, embryo.log2FoldChange)
interactions[, baitName := factor(baitName, unique(baitName))]

pdf("analysis/balancer/diffCaptureC_vs_DE.pdf", width = 15, height = 5.5)

color_zero_point = interactions[, (0 - min(embryo.log2FoldChange)) / (max(embryo.log2FoldChange) - min(embryo.log2FoldChange))]
p1 <- ggplot(interactions) + 
  aes(x = baitName, y = log2FoldChange, col = embryo.log2FoldChange, shape = embryo.padj < 0.05) + 
  geom_line(aes(group = baitName), size = 0.3) +
  geom_point(size = 1.2) + 
  scale_color_gradientn(name = "Gene expression log2 fold change", 
                        colors = c("red", "darkseagreen3", "dodgerblue"), 
                        values = c(0,color_zero_point - 0.11, color_zero_point, color_zero_point + 0.08, 1) ) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  scale_shape_manual(values = c(1,19), name = NULL, labels = c("DE gene", "non\u00adDE gene")) + 
  guides(shape = guide_legend(override.aes = list(size=3))) +
  scale_y_continuous(breaks = seq(-10,10,2)) +
  my_theme + 
  xlab(NULL) +
  ylab("Differential Capture\u00adC contact log2 fold change") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
        legend.position = "bottom", 
        legend.box.just = "left",
        panel.grid.major.x = element_blank()) +
  ggtitle("Strength of significant differential interactions (ordered by DE)")

print(p1)
dev.off()
