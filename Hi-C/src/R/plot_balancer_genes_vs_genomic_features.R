options(warn = 1)

source("src/R/functions_balancer_genes.R")
source("src/R/functions_balancer_annotations.R")
source("src/R/logFC_vs_dist.R")

theme_set(theme_cowplot(font_size = 11)) # reduce default font size
ts <- theme_get()$plot.subtitle
ts$hjust <- 0.5
theme_update(plot.subtitle = ts) # , legend.title = theme_get()$legend.text


for (subset in c("all", "fc1.5"))
# for (subset in c("all"))
{
  stopifnot(grepl("^chr[23]", genes$chrom))

  if (subset == "all")
  {
    genes_subset <- genes
    genes_subset_diffCaptureC <- genes_diffCaptureC
  }
  else if (subset == "non-trivial")
  {
    genes_subset <- genes[trivial_DE == F, ]
    genes_subset_diffCaptureC <- genes_diffCaptureC[trivial_DE == F, ]
  }
  else if (subset == "fc1.5")
  {
    genes_subset <- genes[signf != "s" | abs(log2FoldChange) > log2(1.5), ]
    genes_subset_diffCaptureC <- genes_diffCaptureC[signf != "s" | abs(log2FoldChange) > log2(1.5), ]
  }
  else
    stop("unknown subset: ", subset)

  genes_subset <- sort_GRanges_by_log2FoldChange(genes_subset[signf != "n"])
  genes_subset_diffCaptureC <- sort_GRanges_by_log2FoldChange(genes_subset_diffCaptureC[signf != "n"])

  cat("\nsubset: ", subset)
  print(table(genes_subset$signf))
  cat("\nsubset_diffCaptureC: ", subset)
  print(table(genes_subset_diffCaptureC$signf))

  set.seed(4242)
  pdf(paste0("analysis/balancer/plot_DE_", subset, "_vs_genomic_features.pdf"), width = 5, height = 5, onefile = T)

  glm_pvalues <<- NULL

  # print(make_composite_plot_genes(genes_subset, DHS_gr, bp_gr, ggtitle = "DHSes in balancer"))

  # print(make_composite_plot_genes(genes_subset, MEI_bal_gr, bp_gr, ggtitle = "MEIs in balancer"))
  # print(make_composite_plot_genes(genes_subset, MEI_vrg_gr, bp_gr, ggtitle = "MEIs in wild\u00adtype"))

  # print(make_composite_plot_genes(genes_subset, DEL_bal_gr, bp_gr, ggtitle = "Deletions in balancer"))
  # print(make_composite_plot_genes(genes_subset, DEL_vrg_gr, bp_gr, ggtitle = "Deletions in wild\u00adtype"))
  print(make_composite_plot_genes(genes_subset, c(DEL_bal_gr, DEL_vrg_gr), bp_gr, ggtitle = "Allele\u00adspecific deletions", ylabeller = function(x) paste0("  ", x)))

  # print(make_composite_plot_genes(genes_subset, DUP_bal_gr, bp_gr, ggtitle = "Duplications in balancer"))
  # print(make_composite_plot_genes(genes_subset, DUP_vrg_gr, bp_gr, ggtitle = "Duplications in wild\u00adtype"))

  print(make_composite_plot_genes(genes_subset, SNV_gr, bp_gr, ggtitle = "Allele\u00adspecific SNVs", zmax = 100, ylabeller = function(x) paste0("   ", x)))

  # print(make_composite_plot_genes(genes_subset, SNV_common_gr, bp_gr, ggtitle = "Common SNVs", zmax = 100, ylabeller = function(x) paste0("   ", x)))

  # print(make_composite_plot_genes(genes_subset, bp_gr, GRanges(), ggtitle = "Balancer breakpoints"))

  # gr <- diffIS_bal_gr
  # gr$log2FoldChange <- NULL
  # print(make_composite_plot_genes(genes_subset, gr, bp_gr, ggtitle = "Regions more insulated in balancer"))
  # gr <- diffIS_vrg_gr
  # gr$log2FoldChange <- NULL
  # print(make_composite_plot_genes(genes_subset, gr, bp_gr, ggtitle = "Regions more insulated in wild\u00adtype"))

  # print(make_composite_plot_genes(genes_subset, diffTADs_common_gr, bp_gr, ggtitle = "Common TAD boundaries"))
  # print(make_composite_plot_genes(genes_subset, diffTADs_bal_gr, bp_gr, ggtitle = "Balancer\u00adspecific TAD boundaries"))
  # print(make_composite_plot_genes(genes_subset, diffTADs_vrg_gr, bp_gr, ggtitle = "Wild\u00adtype\u00adspecific TAD boundaries"))

  # print(make_composite_plot_genes(genes_subset, diffHiC_gr, bp_gr, ggtitle = "Differential Hi\u00adC contacts",
  #   legend_name = "Hi\u00adC contact\nlog2 fold change"))
  # print(make_composite_plot_genes(genes_subset, unique(diffHiC_gr), bp_gr, ggtitle = "Differential Hi\u00adC contacts (unique)",
  #   legend_name = "Hi\u00adC contact\nlog2 fold change"))

  # print(make_composite_plot_genes(genes_subset_diffCaptureC, diffCaptureC_gr, bp_gr, ggtitle = "Differential Capture\u00adC contacts",
  #   legend_name = "Capture\u00adC contact\nlog2 fold change"))
  # print(make_composite_plot_genes(genes_subset_diffCaptureC, unique(diffCaptureC_gr), bp_gr, ggtitle = "Differential Capture\u00adC contacts (unique)",
  #   legend_name = "Capture\u00adC contact\nlog2 fold change"))

  print(make_composite_plot_genes(genes_subset, diffHiC_gr, bp_gr, ggtitle = "Differential Hi\u00adC contacts from the TSS",
    extra_filtering = "bait", ribbon_label = "Differential\ncontact density", legend_name = "Differential\nHi\u00adC contact\nlog2 fold change"))
  print(make_composite_plot_genes(genes_subset_diffCaptureC, diffCaptureC_gr, bp_gr, ggtitle = "Differential Capture\u00adC contacts from the TSS",
    extra_filtering = "gene_id", sample_size = Inf, ribbon_label = "Differential\ncontact density", legend_name = "Differential\nCapture\u00adC contact\nlog2 fold change", ylabeller = function(x) paste0("  ", format(x, digits = 1))))

  # write.table(glm_pvalues, file = paste0("analysis/balancer/glm_pvalues_", subset, "_genes_vs_genomic_features.tsv"), sep = "\t", quote = F, row.names = F, col.names = T)

  dev.off()


  # zoomed-in plots for the supplement
  set.seed(4242)
  pdf(paste0("analysis/balancer/plot_DE_", subset, "_vs_genomic_features_zoom.pdf"), width = 5, height = 5, onefile = T)

  print(make_composite_plot_genes(genes_subset, diffHiC_gr, bp_gr, max_dist = 50e3, ggtitle = "Differential Hi\u00adC contacts from the TSS",
    extra_filtering = "bait", ribbon_label = "Differential\ncontact density", legend_name = "Differential\nHi\u00adC contact\nlog2 fold change"))
  print(make_composite_plot_genes(genes_subset_diffCaptureC, diffCaptureC_gr, bp_gr, max_dist = 50e3, ggtitle = "Differential Capture\u00adC contacts from the TSS",
    extra_filtering = "gene_id", sample_size = Inf, ribbon_label = "Differential\ncontact density", legend_name = "Differential\nCapture\u00adC contact\nlog2 fold change", ylabeller = function(x) paste0("  ", format(x, digits = 1))))

  dev.off()
}
