library(data.table)
library(ggplot2)
library(scales)
require(cowplot)
require(rcartocolor)
require(viridis)

theme_set(theme_cowplot(font_size = 11)) # reduce default font size
ts <- theme_get()$plot.subtitle
ts$hjust <- 0.5
theme_update(plot.subtitle = ts) # , legend.title = theme_get()$legend.text

source("src/R/functions_balancer_genes.R")
source("src/R/functions_balancer_annotations.R")


#  "mean of normalized counts" vs. "log_2 fold change"

g_legend <- function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    legend
}


plot_MA_simple <- function(d, min_lfc = log2(1.5), alpha = 0.05)
{
  assert_that(is.data.table(d))
  assert_that("log2FoldChange" %in% colnames(d))
  assert_that("padj" %in% colnames(d))

  d$sign = ifelse(d$log2FoldChange < 0, 'increased in wild type', 'increased in balancer')
  d[padj >= alpha,]$sign = 'not significant'
  
  y_limits = c(min(d$log2FoldChange), max(d$log2FoldChange))

  # Main MA plot
  maplot <- ggplot(d) +
    aes(baseMean, log2FoldChange, col = signf) +
    geom_point(size=0.75, alpha = 0.6) +
    scale_color_manual(values = gene.colors, labels = c(s = paste0(gene.labels["s"], " (", scales::comma(sum(genes_tested$signf == "s")), ")"),
      i = paste0(gene.labels["i"], " (", scales::comma(sum(genes_tested$signf == "i")), ")")), name = NULL) +
    # scale_alpha_manual(values = gene.colors, labels = gene.labels, name = NULL) +
    #geom_hline(yintercept = c(- min_lfc, min_lfc), linetype = "dashed", col = "darkorange", size = 0.6) +
    scale_x_log10(breaks=c(100,1000,10000,100000), label=comma) + 
    scale_y_continuous(breaks=seq(-20,20,2)) +
    coord_cartesian(ylim = y_limits) +
    background_grid(major = "xy", minor = "x") +
    # theme_minimal() +
    xlab("Mean of normalized read counts") + 
    ylab(expression(paste("Gene expression ", log[2] * " fold change"))) +
    theme(legend.position = "bottom", legend.justification = c(0.5, 1)) +
    guides(alpha = guide_legend(title.position="top"),
      color = guide_legend(title.position="top", override.aes = list(size = 1.5)))
  return(maplot)
}

pdf("analysis/balancer/plot_DE_genes.pdf", width = 3, height = 3.5)
print(plot_MA_simple(genes_tested))
dev.off()


## Returns a list with 3 object
plot_MA <- function(d, min_lfc = log2(1.5), alpha = 0.05) {

  assert_that(is.data.table(d))
  assert_that("log2FoldChange" %in% colnames(d))
  assert_that("padj" %in% colnames(d))

  d$sign = ifelse(d$log2FoldChange < 0, 'increased in wild type', 'increased in balancer')
  d[padj >= alpha,]$sign = 'not significant'
  
  y_limits = c(min(d$log2FoldChange), max(d$log2FoldChange))

  # Main MA plot
  ncol <- ifelse(length(unique(d$trivial_DE_class)) >= 3, 3, 1)
  maplot <- ggplot(d) +
    aes(baseMean, log2FoldChange, alpha = trivial_DE_class, col = trivial_DE_class) +
    geom_point(size=0.75) +
    DE_alpha_counted + DE_color_counted +
    #geom_hline(yintercept = c(- min_lfc, min_lfc), linetype = "dashed", col = "darkorange", size = 0.6) +
    scale_x_log10(breaks=c(100,1000,10000,100000), label=comma) + 
    scale_y_continuous(breaks=seq(-20,20,2)) +
    coord_cartesian(ylim = y_limits) +
    background_grid(major = "xy", minor = "x") +
    # theme_minimal() +
    xlab("Mean of normalized read counts") + 
    ylab(expression(paste("Gene expression ", log[2] * " fold change"))) +
    theme(legend.position = "bottom") +
    # theme(legend.position = c(1, 0.2)) +
    guides(alpha = guide_legend(title.position="top", ncol = ncol),
      color = guide_legend(title.position="top", ncol = ncol, override.aes = list(size = 1.5)))

  # Histogram
  histogram <- ggplot(d[padj <= alpha, ]) + 
    aes(log2FoldChange, fill = trivial_DE_class) + 
    geom_histogram(binwidth=0.25, boundary=0) +
    background_grid(major = "xy", minor = "") +
    # theme_minimal() + 
    scale_x_continuous(breaks=seq(-20,20,2)) +
    coord_flip(xlim = y_limits) +
    #geom_vline(xintercept = c(- min_lfc, min_lfc), linetype = "dashed", col = "darkorange", size = 0.6) +
    ylab ("Count of DE genes") + 
    guides(fill = FALSE) +
    theme(axis.title.y = element_blank()) +
    DE_fill

  # legend separately
  # legend <- g_legend(maplot)
  # maplot <- maplot + theme(legend.position = 'none')

  # Arrange plots together
  return (list(maplot = maplot,
              # legend = legend,
              histogram = histogram) )
}


parts <- plot_MA(genes_tested)
upper_row <- plot_grid(ggplot(), parts[["maplot"]], parts[["histogram"]], align = "h", axis = "tb", rel_widths = c(2,5,2), nrow = 1)
ggsave(upper_row, filename = "analysis/balancer/plot_DE_genes_by_affected.pdf", width = 6.5, height = 5)

dt <- copy(genes_tested)
dt[, trivial_DE_class := ifelse(trivial_DE_class == "ZZZnon-DE", "ZZZnon-DE", "none")]
print(table(dt$signf))
parts <- plot_MA(dt)
upper_row <- plot_grid(ggplot(), parts[["maplot"]], parts[["histogram"]], align = "h", axis = "tb", rel_widths = c(2,5,2), nrow = 1)
ggsave(upper_row, filename = "analysis/balancer/plot_DE_genes_not_by_affected.pdf", width = 6.5, height = 5)


## Returns a list with 3 object
plot_MA_by_diff_contacts <- function(d, genes_to_diff, what = "DHS_James_overlap", min_lfc = log2(1.5), alpha = 0.05) {
  d <- copy(d)

  assert_that(is.data.table(d))
  assert_that("log2FoldChange" %in% colnames(d))
  assert_that("padj" %in% colnames(d))

  if ("HiClog2FoldChange" %in% names(genes_to_diff))
    genes_to_diff[, int.log2FoldChange := HiClog2FoldChange]
  else if ("CaptureClog2FoldChange" %in% names(genes_to_diff))
    genes_to_diff[, int.log2FoldChange := CaptureClog2FoldChange]

  dt <- genes_to_diff[, list(
    int.N = .N,
    DHS_James_overlap = any(DHS_James_overlap),
    fraction_DHS_James_overlap = mean(DHS_James_overlap),
    int.log2FoldChange = int.log2FoldChange[which.max(abs(int.log2FoldChange))]
  ), by = "gene_id"]
  d <- merge(d, dt, by = "gene_id", all.x = T)

  d$sign = ifelse(d$log2FoldChange < 0, 'increased in wild type', 'increased in balancer')
  d[padj >= alpha,]$sign = 'not significant'
  
  y_limits = c(min(d$log2FoldChange), max(d$log2FoldChange))

  # Main MA plot
  if (what == "DHS_James_overlap")
    maplot <- ggplot(d) +
      aes(baseMean, log2FoldChange, color = DHS_James_overlap, shape = signf) +
      geom_point(alpha=0.6) +
      scale_color_viridis(discrete = T, end = 2/3, na.value = "#808080") +
      # scale_color_manual(values = carto_pal(name = "Vivid")[c(2, 1)], na.value = "#808080") +
      scale_x_log10(breaks=c(100,1000,10000,100000), label=comma) + 
      xlab("Mean of normalized read counts") + 
      guides(color = guide_legend(title.position="top"))
  else if (what == "fraction_DHS_James_overlap")
    maplot <- ggplot(d) +
      aes(baseMean, log2FoldChange, color = fraction_DHS_James_overlap, shape = signf) +
      geom_point(alpha=0.6) +
      scale_color_viridis(na.value = "#808080") +
      scale_x_log10(breaks=c(100,1000,10000,100000), label=comma) + 
      xlab("Mean of normalized read counts") + 
      guides(color = guide_colorbar(title.position="top", barwidth = 10))
  else if (what == "fraction_DHS_James_overlap_x")
    maplot <- ggplot(d) +
      aes(fraction_DHS_James_overlap, log2FoldChange, color = int.N, shape = signf) +
      geom_point(alpha = 0.6, position = position_jitter(height = 0)) +
      scale_color_viridis(name = "Number of differential contacts", limits = c(0, NA)) +
      guides(color = guide_colorbar(title.position="top", barwidth = 10))
  else if (what == "int.N_x")
    maplot <- ggplot(d) +
      geom_point(alpha = 0.6, position = position_jitter(height = 0)) +
      aes(int.N, log2FoldChange, color = fraction_DHS_James_overlap, shape = signf) +
      scale_color_viridis(na.value = "#808080") +
      xlab("Number of differential contacts") + 
      guides(color = guide_colorbar(title.position="top", barwidth = 10))
  else if (what == "int.log2FoldChange")
    maplot <- ggplot(d) +
      aes(baseMean, log2FoldChange, color = int.log2FoldChange, shape = signf) +
      geom_point(alpha=0.6) +
      scale_color_gradient2(low = "#3182bd", high = "#de2d26", na.value = "#808080", limits = c(-2, 2), oob = scales::squish,
        name = expression("Interaction " * log[2] * " fold change")) +
      scale_x_log10(breaks=c(100,1000,10000,100000), label=comma) + 
      xlab("Mean of normalized read counts") + 
      guides(color = guide_colorbar(title.position="top"))
  else if (what == "int.N")
    maplot <- ggplot(d) +
      aes(baseMean, log2FoldChange, color = int.N, shape = signf) +
      geom_point(alpha=0.6) +
      scale_color_viridis(limits = c(0, NA)) +
      scale_x_log10(breaks=c(100,1000,10000,100000), label=comma) + 
      xlab("Mean of normalized read counts") + 
      guides(color = guide_colorbar(title.position="top", barwidth = 10))

  maplot <- maplot +
    # geom_point(size=0.75) +
    scale_shape_manual(name = "Differential expression", labels = gene.labels, values = c("s" = 19, "i" = 1, "n" = 3)) +
    # DE_alpha_counted + DE_color_counted +
    #geom_hline(yintercept = c(- min_lfc, min_lfc), linetype = "dashed", col = "darkorange", size = 0.6) +
    scale_y_continuous(breaks=seq(-20,20,2)) +
    coord_cartesian(ylim = y_limits) +
    background_grid(major = "xy", minor = "x") +
    # theme_minimal() +
    ylab(expression(paste("Gene expression ", log[2] * " fold change"))) +
    theme(legend.position = "bottom") +
    # theme(legend.position = c(1, 0.2)) +
    guides(shape = guide_legend(title.position="top"))

  # Histogram
  histogram <- ggplot(d[padj <= alpha, ]) + 
    aes(log2FoldChange, fill = DHS_James_overlap, shape = signf) + 
    geom_histogram(binwidth=0.25, boundary=0) +
    background_grid(major = "xy", minor = "") +
    # theme_minimal() + 
    scale_x_continuous(breaks=seq(-20,20,2)) +
    coord_flip(xlim = y_limits) +
    #geom_vline(xintercept = c(- min_lfc, min_lfc), linetype = "dashed", col = "darkorange", size = 0.6) +
    ylab ("Count of DE genes") + 
    guides(fill = FALSE) +
    theme(axis.title.y = element_blank()) +
    DE_fill

  # legend separately
  # legend <- g_legend(maplot)
  # maplot <- maplot + theme(legend.position = 'none')

  # Arrange plots together
  return (list(maplot = maplot,
              # legend = legend,
              histogram = histogram) )
}

pdf("analysis/balancer/plot_DE_genes_by_diffHiC.pdf", width = 5, height = 4.5)
# parts <- plot_MA_by_diff_contacts(genes_tested[gene_id %in% genes_to_diffHiC_dt$gene_id], genes_to_diffHiC_dt, what = "DHS_James_overlap")
# print(parts[["maplot"]])
# parts <- plot_MA_by_diff_contacts(genes_tested[gene_id %in% genes_to_diffHiC_dt$gene_id], genes_to_diffHiC_dt, what = "fraction_DHS_James_overlap")
# print(parts[["maplot"]])
parts <- plot_MA_by_diff_contacts(genes_tested[gene_id %in% genes_to_diffHiC_dt$gene_id], genes_to_diffHiC_dt, what = "fraction_DHS_James_overlap_x")
print(parts[["maplot"]])
parts <- plot_MA_by_diff_contacts(genes_tested[gene_id %in% genes_to_diffHiC_dt$gene_id], genes_to_diffHiC_dt, what = "int.N_x")
print(parts[["maplot"]])
# parts <- plot_MA_by_diff_contacts(genes_tested[gene_id %in% genes_to_diffHiC_dt$gene_id], genes_to_diffHiC_dt, what = "int.log2FoldChange")
# print(parts[["maplot"]])
# parts <- plot_MA_by_diff_contacts(genes_tested[gene_id %in% genes_to_diffHiC_dt$gene_id], genes_to_diffHiC_dt, what = "int.N")
# print(parts[["maplot"]])
dev.off()

pdf("analysis/balancer/plot_DE_genes_by_diffHiC_all_genes.pdf", width = 5, height = 4.5)
# parts <- plot_MA_by_diff_contacts(genes_tested, genes_to_diffHiC_dt, what = "DHS_James_overlap")
# print(parts[["maplot"]])
# parts <- plot_MA_by_diff_contacts(genes_tested, genes_to_diffHiC_dt, what = "fraction_DHS_James_overlap")
# print(parts[["maplot"]])
parts <- plot_MA_by_diff_contacts(genes_tested, genes_to_diffHiC_dt, what = "fraction_DHS_James_overlap_x")
print(parts[["maplot"]])
parts <- plot_MA_by_diff_contacts(genes_tested, genes_to_diffHiC_dt, what = "int.N_x")
print(parts[["maplot"]])
# parts <- plot_MA_by_diff_contacts(genes_tested, genes_to_diffHiC_dt, what = "int.log2FoldChange")
# print(parts[["maplot"]])
# parts <- plot_MA_by_diff_contacts(genes_tested, genes_to_diffHiC_dt, what = "int.N")
# print(parts[["maplot"]])
dev.off()

pdf("analysis/balancer/plot_DE_genes_by_diffCaptureC.pdf", width = 5, height = 4.5)
# parts <- plot_MA_by_diff_contacts(genes_tested[gene_id %in% genes_diffCaptureC$gene_id], genes_to_diffCaptureC_dt, what = "DHS_James_overlap")
# print(parts[["maplot"]])
# parts <- plot_MA_by_diff_contacts(genes_tested[gene_id %in% genes_diffCaptureC$gene_id], genes_to_diffCaptureC_dt, what = "fraction_DHS_James_overlap")
# print(parts[["maplot"]])
parts <- plot_MA_by_diff_contacts(genes_tested[gene_id %in% genes_diffCaptureC$gene_id], genes_to_diffCaptureC_dt, what = "fraction_DHS_James_overlap_x")
print(parts[["maplot"]])
parts <- plot_MA_by_diff_contacts(genes_tested[gene_id %in% genes_diffCaptureC$gene_id], genes_to_diffCaptureC_dt, what = "int.N_x")
print(parts[["maplot"]])
# parts <- plot_MA_by_diff_contacts(genes_tested[gene_id %in% genes_diffCaptureC$gene_id], genes_to_diffCaptureC_dt, what = "int.log2FoldChange")
# print(parts[["maplot"]])
# parts <- plot_MA_by_diff_contacts(genes_tested[gene_id %in% genes_diffCaptureC$gene_id], genes_to_diffCaptureC_dt, what = "int.N")
# print(parts[["maplot"]])
dev.off()
