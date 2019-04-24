library(data.table)
library(ggplot2)
library(scales)
library(assertthat)
library(cowplot)
source("scripts/colors.R")

g_legend <- function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    legend
}

## Returns a list with 3 object
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
    aes(baseMean, log2FoldChange, col = sign) +
    geom_point(alpha = 0.75, size=0.75) +
    scale_color_manual(values = c(`increased in wild type`      = as.character(my_colors["wildtype"]),
                                  `increased in balancer`       = as.character(my_colors["balancer"]),
                                  `not significant`             = "#dddddd")) +
    #geom_hline(yintercept = c(- min_lfc, min_lfc), linetype = "dashed", col = "darkorange", size = 0.6) +
    scale_x_log10(breaks=c(100,1000,10000,100000), label=comma) + 
    scale_y_continuous(breaks=seq(-20,20,2)) +
    coord_cartesian(ylim = y_limits, expand = F) +
    theme_minimal() +
    xlab ("Mean read count") + 
    ylab ("Fold Change (log2)") +
    theme(legend.position = "bottom") +
    guides(colour = guide_legend(title = "Significant genes (5% FDR)", 
                                 override.aes = list(size=3)))

  # Histogram
  histogram <- ggplot(d[padj <= alpha, ]) + 
    aes(log2FoldChange, fill = sign) + 
    geom_histogram(binwidth=0.1) +
    theme_minimal() + 
    scale_x_continuous(breaks=seq(-20,20,2)) +
    coord_flip(xlim = y_limits, expand = F) +
    #geom_vline(xintercept = c(- min_lfc, min_lfc), linetype = "dashed", col = "darkorange", size = 0.6) +
    ylab ("Count") + 
    guides(fill = FALSE) +
    theme(axis.title.y = element_blank()) +
    scale_fill_manual(values = c(`increased in wild type`      = as.character(my_colors["wildtype"]),
                                 `increased in balancer`       = as.character(my_colors["balancer"]),
                                  `not significant`            = "#dddddd"))

  # legend separately
  legend <- g_legend(maplot)
  maplot <- maplot + theme(legend.position = 'none')

  # Arrange plots together
  return (list(maplot = maplot,
              legend = legend,
              histogram = histogram) )
}