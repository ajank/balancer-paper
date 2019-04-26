source("scripts/DESeq_plots.R")
library(cowplot)

d = fread(snakemake@input[[1]])
parts <- plot_MA(d)

title <- ggdraw() + draw_label("Allele-specific gene expression in 6-8h embryos (chr 2 & 3)", fontface='bold')
upper_row <- plot_grid(parts[["maplot"]], parts[["histogram"]], align = "h", rel_widths = c(5,2))

final <- plot_grid(title, upper_row, parts$legend, ncol = 1, rel_heights = c(1,8,1))

ggsave(final, filename = snakemake@output[[1]], width = 7.5, height = 5)
