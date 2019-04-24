source("scripts/DESeq_plots.R")
library(cowplot)


plot_MA <- function(d, min_lfc = log2(1.5), alpha = 0.05) {
  
  assert_that(is.data.table(d))
  assert_that("log2FoldChange" %in% colnames(d))
  assert_that("padj" %in% colnames(d))
  
  d$sign = ifelse(d$log2FoldChange < 0, 'increased in wild type', 'increased in balancer')
  d[padj >= alpha,]$sign = 'not significant'
  
  y_limits = c(min(d$log2FoldChange) * 1.05,
               max(d$log2FoldChange) * 1.05)
  
  # Main MA plot
  maplot <- ggplot(d) +
    aes(baseMean, log2FoldChange, col = padj<0.05) +
    geom_point(alpha = 0.75, size=0.5) +
    scale_color_manual(values = c(`TRUE` = "orangered", `FALSE` = "grey"), 
                       labels =c(`TRUE`  = "significant (FDR 5%)",
                                 `FALSE` = "not significant")) +
    scale_x_log10(breaks=c(100,1000,10000,100000), label=comma) + 
    scale_y_continuous(breaks=seq(-20,20,2)) +
    coord_cartesian(ylim = y_limits, expand = F) +
    theme_minimal() +
    xlab ("Mean read count") + 
    ylab ("Fold Change (log2)") +
    theme(legend.position = "bottom") +
    guides(colour = guide_legend(title = "Significant genes (5% FDR)", 
                                 override.aes = list(size=3)))

    maplot
}

d = fread(snakemake@input[[1]])
ggsave(plot_MA(d), filename = snakemake@output[[1]], width = 7, height = 4)
