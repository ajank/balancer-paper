options(warn = 1)

require(cowplot)
require(ggplot2)
require(rcartocolor)
require(scales) # comma

source("src/R/functions_balancer_genes.R")
source("src/R/functions_balancer_annotations.R")
source("src/R/logFC_vs_dist.R")

theme_set(theme_cowplot(font_size = 11)) # reduce default font size
ts <- theme_get()$plot.subtitle
ts$hjust <- 0.5
theme_update(plot.subtitle = ts) # , legend.title = theme_get()$legend.text
theme_update(
  strip.text.x = element_text(margin = margin(t = 11 / 2, b = 11 / 2), size = 11),
  strip.text.y = element_text(margin = margin(t = 11 / 2, b = 11 / 2), size = 11)
)

mycol <- c(carto_pal(name = "Bold")[1], as.character(gene.colors[c("i")]))

#
#  load the data
#

diffHiC_ov <- list()

# id_list <- c("DE_all_CaptureC_all", "DE_all_CaptureC_reduced_all", "DE_all_HiC_all_disruptedTADs", "DE_all_HiC_all",
#   "DE_fc1.5_CaptureC_fc1.5", "DE_fc1.5_CaptureC_reduced_fc1.5", "DE_fc1.5_HiC_fc1.5_disruptedTADs", "DE_fc1.5_HiC_fc1.5")
id_list <- c("DE_all_CaptureC_reduced_all", "DE_all_HiC_all")

for (id in id_list)
{
  load(paste0("analysis/balancer/diffHiC_overlap/", id, "_ov_list.Rdata")) # ov_list
  diffHiC_ov[[id]] <- ov_list
}

level.labels <- c(n = "genes not tested",
  i = "non\u00adDE genes",
  iu = "upregulated\nnon\u00adDE genes", id = "downregulated\nnon\u00adDE genes",
  s = "DE genes",
  sm = "DE genes, moderate fold change", sh = "DE genes, high fold change", sv = "DE genes, very high fold change",
  su = "upregulated\nDE genes", sd = "downregulated\nDE genes",
  TSS = "    Promoters",
  TSS_active = "        Promoters of active genes",
  TSS_testable = "        Promoters of testable genes",
  TSS_s = "        Promoters of DE genes",
  TSS_i = "        Promoters of non\u00adDE genes",
  TSS_inactive = "        Promoters of inactive genes",
  DHS_distal = "    Distal DHSes",
  DHS_James_distal = "        Distal James' DHSes",
  DHS_insulator_distal = "        Distal insulators",
  DHS_enhancer_distal = "        Distal enhancers",
  insulator = "    Insulators",
  insulator_distal = "    Distal insulators",
  enhancer = "    Enhancers",
  enhancer_distal = "    Distal enhancers",
  all = "All differential contacts")

anno.labels <- c(n = "genes not tested",
  TSS = "Promoters",
  TSS_active = "Promoters\nof active\ngenes",
  TSS_s = "Promoters\nof DE\ngenes",
  TSS_i = "Promoters\nof non\u00adDE\ngenes",
  TSS_inactive = "Promoters\nof inactive\ngenes",
  DHS_distal = "Distal\nDHSes",
  DHS_insulator_distal = "Distal\ninsulators",
  DHS_enhancer_distal = "Distal\nenhancers",
  insulator_distal = "Distal\ninsulators",
  enhancer_distal = "Distal\nenhancers",
  all = "All")

#
#  do all the plotting
#

level_to_fill <- function(l)
  ifelse(as.character(l) %in% names(gene.colors), as.character(l), "n")

for (id in id_list)
{
  pdf(paste0("analysis/balancer/diffHiC_overlap/", id, "_statistics.pdf"), width = 6, height = 2)

  # p <- ggplot(diffHiC_ov[[id]]$fc_by_fc$aggr[anno == anno[1]], aes(level, n_genes, fill = level_to_fill(level))) +
  #   geom_bar(stat = "identity") +
  #   labs(x = "Originating from ...", y = "Number of genes") +
  #   geom_text(aes(y = n_genes + 200, label = comma(n_genes)), hjust = 0) +
  #   scale_x_discrete(labels = level.labels) +
  #   scale_y_continuous(labels = comma, limits = c(0, max(diffHiC_ov[[id]]$fc_by_fc$aggr$n_genes) + 1000)) +
  #   scale_fill_manual(values = gene.colors, labels = gene.labels, name = NULL, guide = F) +
  #   coord_flip()
  # print(p)

  # p <- ggplot(diffHiC_ov[[id]]$fc_by_fc$aggr[anno == anno[1]], aes(level, n_DC_genes, fill = level_to_fill(level))) +
  #   geom_bar(stat = "identity") +
  #   labs(x = "Originating from ...", y = "Number of genes assigned to differential Hi\u00adC contacts") +
  #   geom_text(aes(y = n_DC_genes + 20, label = comma(n_DC_genes)), hjust = 0) +
  #   scale_x_discrete(labels = level.labels) +
  #   scale_y_continuous(labels = comma, limits = c(0, max(diffHiC_ov[[id]]$fc_by_fc$aggr$n_DC_genes) + 150)) +
  #   scale_fill_manual(values = gene.colors, labels = gene.labels, name = NULL, guide = F) +
  #   coord_flip()
  # print(p)

  p <- ggplot(diffHiC_ov[[id]]$fc_by_fc$aggr[anno == anno[1]], aes(level, n_DC, fill = level_to_fill(level))) +
    geom_bar(stat = "identity") +
    labs(x = "Originating from ...", y = "Number of differential Hi\u00adC contacts") +
    geom_text(aes(y = n_DC + 20, label = comma(n_DC)), hjust = 0) +
    scale_x_discrete(labels = level.labels) +
    scale_y_continuous(labels = comma, limits = c(0, max(diffHiC_ov[[id]]$fc_by_fc$aggr$n_DC) + 300)) +
    scale_fill_manual(values = gene.colors, labels = gene.labels, name = NULL, guide = F) +
    coord_flip()
  print(p)

  p <- ggplot(diffHiC_ov[[id]]$fc_by_fc$aggr[anno == anno[1]], aes(level, n_DC_gr, fill = level_to_fill(level))) +
    geom_bar(stat = "identity") +
    labs(x = "Originating from ...", y = "Number of assignments of differential Hi\u00adC contacts to features") +
    geom_text(aes(y = n_DC_gr + 40, label = comma(n_DC_gr)), hjust = 0) +
    scale_x_discrete(labels = level.labels) +
    scale_y_continuous(labels = comma, limits = c(0, max(diffHiC_ov[[id]]$fc_by_fc$aggr$n_DC_gr) + 600)) +
    scale_fill_manual(values = gene.colors, labels = gene.labels, name = NULL, guide = F) +
    coord_flip()
  print(p)

  dev.off()
}


for (id in id_list)
{
  pdf(paste0("analysis/balancer/diffHiC_overlap/", id, "_vs_anno_by_signf.pdf"), width = 6, height = 2)

  for (a in unique(diffHiC_ov[[id]]$anno_by_fc$shuffled$anno))
  {
    p <- ggplot(diffHiC_ov[[id]]$anno_by_fc$shuffled[anno == a, ], aes(shuffled)) +
      # facet_wrap(~ level, ncol = 1, as.table = F, labeller = as_labeller(level.labels)) +
      facet_wrap(~ level, scales = "free_y", ncol = 1, as.table = F, labeller = as_labeller(level.labels), strip.position="left") +
      theme(strip.text.y = element_text(angle = 180, hjust = 0)) +
      xlim(c(0, NA)) +
      scale_y_continuous(expand = c(0.1, 0)) +
      labs(x = paste0("Fraction of Hi\u00adC differential contacts\nwith ", level.labels[a], " at the other end"), y = "Originating from ...") +
      # geom_jitter(aes(y = 0), color = mycol[2], alpha = 0.1) +
      geom_density(color = NA, fill = mycol[2], alpha = 0.4) +
      # geom_histogram(binwidth = 0.005, boundary = 0, color = NA, fill = mycol[2], alpha = 0.6) +
      # geom_vline(data = diffHiC_ov[[id]]$anno_by_fc$aggr[anno == a, ], aes(xintercept = mean_shuffled), color = mycol[2]) +
      # geom_vline(data = diffHiC_ov[[id]]$anno_by_fc$aggr[anno == a, ], aes(xintercept = observed), color = mycol[1]) +
      geom_point(data = diffHiC_ov[[id]]$anno_by_fc$aggr[anno == a, ], aes(x = mean_shuffled, y = 0), color = mycol[2]) +
      geom_point(data = diffHiC_ov[[id]]$anno_by_fc$aggr[anno == a, ], aes(x = observed, y = 0), color = mycol[1]) +
      geom_text(data = diffHiC_ov[[id]]$anno_by_fc$aggr[anno == a, ],
        aes(label = ifelse(pval < 0.05, ifelse(pval < 0.001, sprintf("p < 0.001\n\n", pval), sprintf("p = %0.3f\n\n", pval)), "ns\n\n"),
        # x = observed, y = 0), hjust = 0, vjust = 0.5, color = mycol[1]) +
        x = (observed + mean_shuffled) / 2, y = 0), hjust = 0.5, vjust = 0.5, size = 3) +
      # geom_text(data = diffHiC_ov[[id]]$anno_by_fc$aggr[anno == a, ],
      #   aes(label = ifelse(pval_to_next < 0.05, ifelse(pval_to_next < 0.001, sprintf("p = %0.1e", pval_to_next), sprintf("p = %0.3f", pval_to_next)), "ns"),
      #   x = Inf, y = Inf), hjust = 1, vjust = 1, size = 3) +
      theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y = element_blank()) +
      NULL
    print(p)
  }

  dev.off()
}

#
#  do differential contacts of DE genes/distal enhancers/... preferentially contact DE genes/distal enhancers/...?
#

for (id in id_list)
{
  pdf(paste0("analysis/balancer/diffHiC_overlap/", id, "_overlap.pdf"), width = 6, height = 4)

  for (a in unique(diffHiC_ov[[id]]$fc_by_fc$shuffled$anno))
  {
    p <- ggplot(diffHiC_ov[[id]]$fc_by_fc$shuffled[anno == a, ], aes(shuffled)) +
      facet_wrap(~ level, scales = "free_y", ncol = 1, as.table = F, labeller = as_labeller(level.labels), strip.position="left") +
      theme(strip.text.y = element_text(angle = 180, hjust = 0)) +
      xlim(c(0, NA)) +
      scale_y_continuous(expand = c(0.1, 0)) +
      labs(x = paste0("Fraction of Hi\u00adC differential contacts\nwith ", level.labels[a], " at the other end"), y = "Originating from ...") +
      geom_density(color = NA, fill = mycol[2], alpha = 0.4) +
      geom_point(data = diffHiC_ov[[id]]$fc_by_fc$aggr[anno == a, ], aes(x = mean_shuffled, y = 0), color = mycol[2]) +
      geom_point(data = diffHiC_ov[[id]]$fc_by_fc$aggr[anno == a, ], aes(x = observed, y = 0), color = mycol[1]) +
      geom_text(data = diffHiC_ov[[id]]$fc_by_fc$aggr[anno == a, ],
        aes(label = ifelse(pval < 0.05,
          ifelse(pval < 0.001, sprintf("p < 0.001\n\n", pval), sprintf("p = %0.3f\n\n", pval)), "ns\n\n"),
        x = (observed + mean_shuffled) / 2, y = 0), hjust = 0.5, vjust = 0.5, size = 3) +
      # geom_text(data = diffHiC_ov[[id]]$fc_by_fc$aggr[anno == a, ],
      #   aes(label = ifelse(pval_to_next < 0.05,
      #     ifelse(pval_to_next < 0.001, sprintf("p = %0.1e", pval_to_next), sprintf("p = %0.3f", pval_to_next)), "ns"),
      #   x = Inf, y = Inf), hjust = 1, vjust = 1, size = 3) +
      theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y = element_blank()) +
      NULL
    print(p)
  }

  dev.off()
}


for (id in id_list)
{
  level_anno_labeller <- function(dt)
  {
    if ("level" %in% names(dt))
    {
      ll <- level.labels
      for (l in names(ll))
      {
        dts <- diffHiC_ov[[id]]$fc_by_fc$aggr[level == l]
        stopifnot(length(l) == 1)
        tr <- gsub("^ *", "", ll[l])
        tr <- paste0(tolower(substring(tr, 1, 1)), substring(tr, 2))
        sp <- gsub("[^ ].*", "", ll[l])
        ll[l] <- paste0(sp, "Differential ", ifelse(grepl("HiC", id), "Hi\u00adC", "Capture\u00adC"), " contacts (", comma(dts$n_DC[1]), ")")
        if (l != "all")
          ll[l] <- paste0(
            ll[l], "\n", sp, "from ", tr,
            ifelse(dts$n_bait_genes[1] > 0, paste0(" (", comma(dts$n_bait_genes[1]), ")"), ""),
            ifelse(dts$n_bait_dhs[1] > 0, paste0(" (", comma(dts$n_bait_dhs[1]), ")"), "")
          )
      }
    }
    else
    {
      ll <- anno.labels
    }

    if (grepl("DE_fc1.5", id))
      for (l in names(ll))
        {
          ll[l] <- gsub("of DE genes", "of strongly DE genes", ll[l])
          ll[l] <- gsub("of DE\ngenes", "of strongly\nDE genes", ll[l])
        }

    return(as_labeller(ll)(dt))
  }

  dt <- diffHiC_ov[[id]]$fc_by_fc$shuffled
  dt[, type := "shuffled"]
  dti <- with(diffHiC_ov[[id]]$fc_by_fc$aggr, rbind(
    data.table(level = level, anno = anno, type = "shuffled", value = mean_shuffled),
    data.table(level = level, anno = anno, type = "observed", value = observed)
  ))

  # pdf(paste0("analysis/balancer/diffHiC_overlap/", id, "_overlap_oneplot.pdf"), width = 10, height = 5.5)

  # p <- ggplot(dt, aes(shuffled, fill = type)) +
  #   facet_grid(level ~ anno, labeller = level_anno_labeller, scales = "free", switch = "y") +
  #   theme(strip.text.y = element_text(angle = 180, hjust = 0)) +
  #   # xlim(c(0, NA)) +
  #   expand_limits(x = c(0, 0.23)) +
  #   scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
  #   scale_y_continuous(expand = c(0.1, 0)) +
  #   labs(x = paste0("Fraction of differential contacts having their other end in ..."),
  #     y = NULL) + # paste0("Differential ", ifelse(grepl("HiC", id), "Hi\u00adC", "Capture\u00adC"), " contacts originating from ...")) +
  #   geom_density(color = NA, alpha = 0.4) +
  #   geom_point(data = dti, aes(x = value, y = 0, color = type)) +
  #   scale_color_manual(values = mycol, labels = c(observed = "observed", shuffled = "randomly shuffled"),
  #     name = "Differential contacts") +
  #   scale_fill_manual(values = mycol, labels = c(observed = "observed", shuffled = "randomly shuffled"),
  #     name = "Differential contacts") +
  #   geom_text(data = diffHiC_ov[[id]]$fc_by_fc$aggr,
  #     aes(label = ifelse(pval < 0.05,
  #       ifelse(pval < 0.001, sprintf("p < 0.001\n\n", pval), sprintf("p = %0.3f\n\n", pval)), ""),
  #     # x = pmax((observed + mean_shuffled) / 2, 0.08), y = 0), hjust = 0.5, vjust = 0.5, size = 3) +
  #     x = 0, y = 0), inherit.aes = F, hjust = 0, vjust = 0.5, size = 3) +
  #   theme(strip.text.x = element_text(margin = margin(5,0,5,0, "pt")), strip.text.y = element_text(margin = margin(0,5,0,5, "pt"))) +
  #   theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y = element_blank()) +
  #   theme(legend.justification = c(0.5, 1), legend.position = "bottom") +
  #   NULL
  # print(p)

  # dev.off()


  dt <- diffHiC_ov[[id]]$fc_by_fc$shuffled[, list(
    lwr_shuffled = quantile(shuffled, 0.025, na.rm = T),
    upr_shuffled = quantile(shuffled, 0.975, na.rm = T)
  ), by = c("level", "anno")]
  dt <- merge(diffHiC_ov[[id]]$fc_by_fc$aggr, dt, by = c("level", "anno"), all = T)

  dt <- dt[level %in% c("all", "TSS_testable", "TSS_s", "TSS_i") & anno %in% c("TSS_s", "TSS_i", "DHS_enhancer_distal", "DHS_distal"), ]
  dt <- dt[, anno := factor(anno, c("TSS_s", "TSS_i", "DHS_enhancer_distal", "DHS_distal")), ]
  # y: diff, diff from promoters, diff strongly DE, diff non-DE
  # x: strongly DE, non-DE, distal enh, distal DHSes

  dti <- with(dt, rbind(
    data.table(level = level, anno = anno, type = "shuffled", value = mean_shuffled, lwr = lwr_shuffled, upr = upr_shuffled),
    data.table(level = level, anno = anno, type = "observed", value = observed, lwr = NA, upr = NA)
  ))

  # pdf(paste0("analysis/balancer/diffHiC_overlap/", id, "_barplot.pdf"), width = 7.5, height = 4.5)

  # p <- ggplot(dti, aes(x = anno, y = value, fill = type)) +
  #   facet_grid(level ~ anno, labeller = level_anno_labeller, scales = "free_x", switch = "y") +
  #   theme(strip.text.y = element_text(angle = 180, hjust = 0)) +
  #   geom_bar(stat = "identity", position = position_dodge(), width = 0.4) +
  #   geom_errorbar(aes(ymin = lwr, ymax = upr), width = .2, position = position_dodge(0.4)) +
  #   geom_text(data = dt, inherit.aes = F,
  #     aes(x = anno, y = max(observed, upr_shuffled) + if (grepl("_HiC", id)) 0.05 else 0, label = ifelse(pval < 0.05,
  #       ifelse(pval < 0.001, sprintf("p < 0.001", pval), sprintf("p = %0.3f", pval)), ""))
  #   ) +
  #   labs(y = NULL, x = "Fraction of differential contacts having their other end in ...") +
  #   scale_x_discrete(position = "top") +
  #   scale_y_continuous(position = "right", breaks = scales::pretty_breaks(n = 3)) +
  #   scale_fill_manual(values = mycol, labels = c(observed = "observed", shuffled = "randomly shuffled"),
  #     name = "Differential contacts") +
  #   theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) +
  #   theme(strip.text.x = element_text(margin = margin(5,0,5,0, "pt")), strip.text.y = element_text(margin = margin(0,5,0,5, "pt"))) +
  #   theme(legend.justification = c(0.5, 1), legend.position = "bottom") +
  #   NULL
  # print(p)

  # dev.off()


  dti1 <- copy(dti)
  dti1[, value := ifelse(type == "shuffled", NA, value)]

  pdf(paste0("analysis/balancer/diffHiC_overlap/", id, "_violinplot.pdf"), width = 7.5, height = 4.5)

  p <- ggplot(dti1, aes(x = anno, y = value, fill = type)) +
    facet_grid(level ~ anno, labeller = level_anno_labeller, scales = "free_x", switch = "y") +
    theme(strip.text.y = element_text(angle = 180, hjust = 0)) +
    geom_errorbar(aes(ymin = value, ymax = value, color = type), width = 0.6, position = position_dodge()) +
    geom_violin(data = dt_shuffled, aes(x = 1.3, y = shuffled, fill = "shuffled"), width = 0.6) +
    # geom_errorbar(aes(ymin = lwr, ymax = upr), width = .2, position = position_dodge(1.2)) +
    geom_boxplot(data = dt_shuffled, aes(x = 1.3, y = shuffled, fill = "shuffled"), width = 0.3, outlier.shape = NA) +
    geom_text(data = dt, inherit.aes = F,
      aes(x = anno, y = max(observed, upr_shuffled) + if (grepl("_HiC", id)) 0.05 else 0, label = ifelse(pval < 0.05,
        ifelse(pval < 0.001, sprintf("p < 0.001", pval), sprintf("p = %0.3f", pval)), ""))
    ) +
    labs(y = NULL, x = "Fraction of differential contacts having their other end in ...") +
    coord_cartesian(xlim = c(1.05, 1.25)) +
    scale_x_discrete(position = "top") +
    scale_y_continuous(position = "right", limits = c(0, NA), expand = c(0, 0)) +
    scale_color_manual(values = mycol, labels = c(observed = "observed", shuffled = "randomly shuffled"),
      name = "Differential contacts") +
    scale_fill_manual(values = mycol, labels = c(observed = "observed", shuffled = "randomly shuffled"),
      name = "Differential contacts") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank()) +
    theme(strip.text.x = element_text(margin = margin(5,0,5,0, "pt")), strip.text.y = element_text(margin = margin(0,5,0,5, "pt"))) +
    theme(legend.justification = c(0.5, 1), legend.position = "bottom") +
    NULL
  print(p)

  dev.off()


  # dt <- diffHiC_ov[[id]]$fc_by_fc$shuffled[, list(
  #   lwr_shuffled = quantile(shuffled, 0.025, na.rm = T),
  #   upr_shuffled = quantile(shuffled, 0.975, na.rm = T)
  # ), by = c("level", "anno")]
  # dt <- merge(diffHiC_ov[[id]]$fc_by_fc$aggr, dt, by = c("level", "anno"), all = T)

  # dti <- with(dt, rbind(
  #   data.table(level = level, anno = anno, type = "shuffled", value = mean_shuffled, lwr = lwr_shuffled, upr = upr_shuffled),
  #   data.table(level = level, anno = anno, type = "observed", value = observed, lwr = NA, upr = NA)
  # ))

  # pdf(paste0("analysis/balancer/diffHiC_overlap/", id, "_barplot_DE_genes.pdf"), width = 6, height = 3.5)

  # for (l in "TSS_s")
  # {
  #   p <- ggplot(dti[level == l], aes(x = anno, y = value, fill = type)) +
  #     geom_bar(stat = "identity", position = position_dodge(), width = 0.4) +
  #     geom_errorbar(aes(ymin = lwr, ymax = upr), width = .2, position = position_dodge(0.4)) +
  #     geom_text(data = dt[level == l], inherit.aes = F,
  #       aes(x = anno, y = max(observed, upr_shuffled) + 0.05, label = ifelse(pval < 0.05, paste0("p = ", pval), ""))
  #     ) +
  #     labs(x = NULL, y = "Fraction of differential contacts\nhaving their other end in ...",
  #       title = paste0("Differential ", ifelse(grepl("HiC", id), "Hi\u00adC", "Capture\u00adC"), " contacts (", comma(dt[level == l]$n_DC[1]),
  #       ") from the TSSes of DE genes (", comma(dt[level == l]$n_bait_genes[1]), ")")) +
  #     scale_x_discrete(labels = anno.labels) +
  #     scale_fill_manual(values = mycol, labels = c(observed = "observed", shuffled = "randomly shuffled"),
  #       name = "Differential contacts") +
  #     theme(legend.justification = c(0.5, 1), legend.position = "bottom") +
  #     theme(axis.text.x = element_text(size = 11)) +
  #     NULL
  #   print(p)
  # }

  # dev.off()
}

#
#  directionality of differential contacts: do differential contacts of DE genes preferentially contact DE genes?
#

for (id in id_list)
{
  pdf(paste0("analysis/balancer/diffHiC_overlap/", id, "_vs_dir_by_dir.pdf"), width = 6, height = 3)

  for (a in unique(diffHiC_ov[[id]]$dir_by_dir$shuffled$anno))
  {
    p <- ggplot(diffHiC_ov[[id]]$dir_by_dir$shuffled[anno == a, ], aes(shuffled)) +
      facet_wrap(~ level, scales = "free_y", ncol = 1, as.table = F, labeller = as_labeller(level.labels), strip.position="left") +
      theme(strip.text.y = element_text(angle = 180, hjust = 0)) +
      xlim(c(0, NA)) +
      scale_y_continuous(expand = c(0.1, 0)) +
      labs(x = paste0("Fraction of Hi\u00adC differential contacts\nwith ", level.labels[a], " at the other end"), y = "Originating from ...") +
      geom_density(color = NA, fill = mycol[2], alpha = 0.4) +
      geom_point(data = diffHiC_ov[[id]]$dir_by_dir$aggr[anno == a, ], aes(x = mean_shuffled, y = 0), color = mycol[2]) +
      geom_point(data = diffHiC_ov[[id]]$dir_by_dir$aggr[anno == a, ], aes(x = observed, y = 0), color = mycol[1]) +
      geom_text(data = diffHiC_ov[[id]]$dir_by_dir$aggr[anno == a, ],
        aes(label = ifelse(pval < 0.05,
          ifelse(pval < 0.001, sprintf("p < 0.001\n\n", pval), sprintf("p = %0.3f\n\n", pval)), "ns\n\n"),
        x = (observed + mean_shuffled) / 2, y = 0), hjust = 0.5, vjust = 0.5, size = 3) +
      # geom_text(data = diffHiC_ov[[id]]$dir_by_dir$aggr[anno == a, ],
      #   aes(label = ifelse(pval_to_next < 0.05,
      #     ifelse(pval_to_next < 0.001, sprintf("p = %0.1e", pval_to_next), sprintf("p = %0.3f", pval_to_next)), "ns"),
      #   x = Inf, y = Inf), hjust = 1, vjust = 1, size = 3) +
      theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y = element_blank()) +
      NULL
    print(p)
  }

  dev.off()
}

for (id in id_list)
{
  pdf(paste0("analysis/balancer/diffHiC_overlap/", id, "_dir_oneplot.pdf"), width = 8, height = 4)

  p <- ggplot(diffHiC_ov[[id]]$dir_by_dir$shuffled, aes(shuffled)) +
    facet_grid(level ~ anno, labeller = as_labeller(level.labels), switch = "y") +
    theme(strip.text.y = element_text(angle = 180, hjust = 0)) +
    xlim(c(0, NA)) +
    scale_y_continuous(expand = c(0.1, 0)) +
    labs(x = paste0("Fraction of differential contacts having their other end in ..."), y = "Differential Hi\u00adC contacts originating from ...") +
    geom_density(color = NA, fill = mycol[2], alpha = 0.4) +
    geom_point(data = diffHiC_ov[[id]]$dir_by_dir$aggr, aes(x = mean_shuffled, y = 0), color = mycol[2]) +
    geom_point(data = diffHiC_ov[[id]]$dir_by_dir$aggr, aes(x = observed, y = 0), color = mycol[1]) +
    geom_text(data = diffHiC_ov[[id]]$dir_by_dir$aggr,
      aes(label = ifelse(pval < 0.05,
        ifelse(pval < 0.001, sprintf("p < 0.001\n\n", pval), sprintf("p = %0.3f\n\n", pval)), "ns\n\n"),
      x = (observed + mean_shuffled) / 2, y = 0), hjust = 0.5, vjust = 0.5, size = 3) +
    # geom_text(data = diffHiC_ov[[id]]$dir_by_dir$aggr,
    #   aes(label = ifelse(pval_to_next < 0.05,
    #     ifelse(pval_to_next < 0.001, sprintf("p = %0.1e", pval_to_next), sprintf("p = %0.3f", pval_to_next)), "ns"),
    #   x = Inf, y = Inf), hjust = 1, vjust = 1, size = 3) +
    theme(strip.text.x = element_text(margin = margin(5,0,5,0, "pt")), strip.text.y = element_text(margin = margin(0,5,0,5, "pt"))) +
    theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.line.y = element_blank()) +
    NULL
  print(p)

  dev.off()
}


pdf("analysis/balancer/diffHiC_overlap/genes_vs_diffHiC.pdf", width = 5, height = 2)
for (id in id_list)
{
  set.seed(4242)
  print(make_ribbon_plot_genes(diffHiC_ov[[id]]$TSS_subset, diffHiC_ov[[id]]$diffHiC_gr_subset, do.call(rbind, diffHiC_ov[[id]]$gmap), ymax = 0.15, ggtitle = id))

  # print(make_composite_plot_genes(diffHiC_ov[[id]]$TSS_subset_diffCaptureC, diffHiC_ov[[id]]$diffHiC_gr_subset, bp_gr, ggtitle = "Differential Capture\u00adC contacts from the TSS",
  #   extra_filtering = "gene_id", sample_size = Inf, ribbon_label = "Differential\ncontact density", legend_name = "Capture\u00adC contact\nlog2 fold change", ylabeller = function(x) paste0("  ", format(x, digits = 1))))
}
dev.off()

#
#  xy-plot of log2FoldChange of genes having a differential Hi-C contact between them
#

for (id in id_list)
{
  # overlap "other ends"
  ov <- findOverlaps(diffHiC_ov[[id]]$diffHiC_gr_subset, GRanges(diffHiC_ov[[id]]$TSS_subset))
  ov_dt <- unique(data.table(id = queryHits(ov), gene_id = diffHiC_ov[[id]]$TSS_subset$gene_id[subjectHits(ov)]))
  gr <- diffHiC_ov[[id]]$diffHiC_gr_subset[ov_dt$id]
  gr$other_end_gene_id <- ov_dt$gene_id

  # overlap baits
  ov <- findOverlaps(GRanges(gr$baitChr, IRanges(gr$baitStart, gr$baitEnd)), GRanges(diffHiC_ov[[id]]$TSS_subset))
  ov_dt <- unique(data.table(id = queryHits(ov), bait_gene_id = diffHiC_ov[[id]]$TSS_subset$gene_id[subjectHits(ov)]))
  both_genes_to_diffHiC_gr <- gr[ov_dt$id]
  both_genes_to_diffHiC_gr$bait_gene_id <- ov_dt$bait_gene_id

  # annotation
  both_genes_to_diffHiC_gr$bait_signf <- genes$signf[match(both_genes_to_diffHiC_gr$bait_gene_id, genes$gene_id)]
  both_genes_to_diffHiC_gr$bait_log2FoldChange <- genes$log2FoldChange[match(both_genes_to_diffHiC_gr$bait_gene_id, genes$gene_id)]
  both_genes_to_diffHiC_gr$bait_FPKM <- genes$FPKM[match(both_genes_to_diffHiC_gr$bait_gene_id, genes$gene_id)]
  both_genes_to_diffHiC_gr$other_end_signf <- genes$signf[match(both_genes_to_diffHiC_gr$other_end_gene_id, genes$gene_id)]
  both_genes_to_diffHiC_gr$other_end_log2FoldChange <- genes$log2FoldChange[match(both_genes_to_diffHiC_gr$other_end_gene_id, genes$gene_id)]
  both_genes_to_diffHiC_gr$other_end_FPKM <- genes$FPKM[match(both_genes_to_diffHiC_gr$other_end_gene_id, genes$gene_id)]
  both_genes_to_diffHiC_gr$distSign <- ceiling((start(both_genes_to_diffHiC_gr) + end(both_genes_to_diffHiC_gr)) / 2) - ceiling((both_genes_to_diffHiC_gr$baitStart + both_genes_to_diffHiC_gr$baitEnd) / 2)
  both_genes_to_diffHiC_gr$paired_signf <- paste(both_genes_to_diffHiC_gr$bait_signf, both_genes_to_diffHiC_gr$other_end_signf)
  both_genes_to_diffHiC_gr$same_dir <- sign(both_genes_to_diffHiC_gr$bait_log2FoldChange) == sign(both_genes_to_diffHiC_gr$other_end_log2FoldChange)

  # filtering
  sel1 <- both_genes_to_diffHiC_gr$bait_gene_id != both_genes_to_diffHiC_gr$other_end_gene_id
  sel2 <- !is.na(both_genes_to_diffHiC_gr$bait_signf) & both_genes_to_diffHiC_gr$bait_signf != "n"
  sel3 <- !is.na(both_genes_to_diffHiC_gr$other_end_signf) & both_genes_to_diffHiC_gr$other_end_signf != "n"
  both_genes_to_diffHiC_gr <- both_genes_to_diffHiC_gr[sel1 & sel2 & sel3]
  diffHiC_ov[[id]]$both_genes_to_diffHiC_gr <- both_genes_to_diffHiC_gr

  if (grepl("DE_fc1.5", id))
    both_genes_to_diffHiC_gr <- both_genes_to_diffHiC_gr[
      with(both_genes_to_diffHiC_gr,
        (bait_signf != "s" | abs(bait_log2FoldChange) > log2(1.5)) &
        (other_end_signf != "s" | abs(other_end_log2FoldChange) > log2(1.5))
      )
    ]


  pdf(paste0("analysis/balancer/diffHiC_overlap/", id, "_both_genes_xyplot.pdf"), width = 6, height = 6)

  p <- ggplot(as.data.table(both_genes_to_diffHiC_gr), aes(bait_log2FoldChange, other_end_log2FoldChange, color = same_dir)) +
    facet_wrap(paired_signf ~ ., labeller = as_labeller(c("s s" = "Both genes DE", "s i" = "Gene 1 DE, gene 2 non\u00adDE", "i s" = "Gene 1 non\u00adDE, gene 2 DE", "i i" = "Both genes non\u00adDE"))) +
    geom_vline(aes(xintercept = 0), lty = 2) +
    geom_hline(aes(yintercept = 0), lty = 2) +
    geom_point(alpha = 0.6) +
    xlab(expression(paste("Gene 1 ", log[2] * " fold change"))) +
    ylab(expression(paste("Gene 2 ", log[2] * " fold change"))) +
    scale_x_continuous(limits = c(-5, 5), expand = c(0.005, 0), oob = scales::squish) +
    scale_y_continuous(limits = c(-5, 5), expand = c(0.005, 0), oob = scales::squish) +
    coord_fixed() +
    # scale_color_manual(values = gene.colors, labels = gene.labels, name = NULL) +
    # theme(legend.justification = c(0, 1), legend.position = c(0, 1)) +
    theme(legend.position = "none") +
    theme(strip.text.x = element_text(margin = margin(5,0,5,0, "pt"))) +
    NULL
  print(p)

  if ("log2FoldChange" %in% names(elementMetadata(both_genes_to_diffHiC_gr)))
  {
    p <- ggplot(as.data.table(both_genes_to_diffHiC_gr), aes(bait_log2FoldChange, other_end_log2FoldChange, color = log2FoldChange)) +
      facet_wrap(paired_signf ~ ., labeller = as_labeller(c("s s" = "Both genes DE", "s i" = "Gene 1 DE, gene 2 non\u00adDE", "i s" = "Gene 1 non\u00adDE, gene 2 DE", "i i" = "Both genes non\u00adDE"))) +
      geom_vline(aes(xintercept = 0), lty = 2) +
      geom_hline(aes(yintercept = 0), lty = 2) +
      geom_point(alpha = 0.6) +
      xlab(expression(paste("Gene 1 ", log[2] * " fold change"))) +
      ylab(expression(paste("Gene 2 ", log[2] * " fold change"))) +
      scale_x_continuous(limits = c(-5, 5), expand = c(0.005, 0), oob = scales::squish) +
      scale_y_continuous(limits = c(-5, 5), expand = c(0.005, 0), oob = scales::squish) +
      coord_fixed() +
      scale_color_gradient2("log2FoldChange of differential Hi\u00adC contacts", low = "#3182bd", high = "#de2d26", limits = c(-2, 2), oob = scales::squish) +
      theme(legend.position = "bottom") +
      # scale_color_manual(values = gene.colors, labels = gene.labels, name = NULL) +
      # theme(legend.justification = c(0, 1), legend.position = c(0, 1)) +
      theme(strip.text.x = element_text(margin = margin(5,0,5,0, "pt"))) +
      NULL
    print(p)

    dt <- unique(with(as.data.table(both_genes_to_diffHiC_gr), rbind(
      data.table(diffc_id = diffc_id, paired_signf = paired_signf, same_dir = same_dir,
        gene_id = bait_gene_id, gene_log2FoldChange = bait_log2FoldChange, log2FoldChange = log2FoldChange),
      data.table(diffc_id = diffc_id, paired_signf = paired_signf, same_dir = same_dir,
        gene_id = other_end_gene_id, gene_log2FoldChange = other_end_log2FoldChange, log2FoldChange = log2FoldChange)
    )))
    p <- ggplot(dt, aes(gene_log2FoldChange, log2FoldChange, color = same_dir)) +
      facet_wrap(paired_signf ~ ., labeller = as_labeller(c("s s" = "Both genes DE", "s i" = "Gene 1 DE, gene 2 non\u00adDE", "i s" = "Gene 1 non\u00adDE, gene 2 DE", "i i" = "Both genes non\u00adDE"))) +
      geom_vline(aes(xintercept = 0), lty = 2) +
      geom_hline(aes(yintercept = 0), lty = 2) +
      geom_point(alpha = 0.6) +
      xlab(expression(paste("Gene ", log[2] * " fold change"))) +
      ylab(expression(paste("Hi\u00adC differential contact ", log[2] * " fold change"))) +
      scale_x_continuous(limits = c(-5, 5), expand = c(0.005, 0), oob = scales::squish) +
      scale_y_continuous(limits = c(-5, 5), expand = c(0.005, 0), oob = scales::squish) +
      coord_fixed() +
      # scale_color_manual(values = gene.colors, labels = gene.labels, name = NULL) +
      # theme(legend.justification = c(0, 1), legend.position = c(0, 1)) +
      theme(legend.position = "none") +
      theme(strip.text.x = element_text(margin = margin(5,0,5,0, "pt"))) +
      NULL
    print(p)
  }

  dev.off()


  pdf(paste0("analysis/balancer/diffHiC_overlap/", id, "_both_DE_genes_xyplot.pdf"), width = 2, height = 4)

  udt <- unique(as.data.table(both_genes_to_diffHiC_gr)[, c("bait_gene_id", "bait_log2FoldChange", "other_end_gene_id", "other_end_log2FoldChange", "same_dir", "paired_signf")])
  p <- ggplot(udt[paired_signf == "s s"], aes(bait_log2FoldChange, other_end_log2FoldChange,
      color = ifelse(same_dir, "concordant", "discordant"))) +
    geom_vline(aes(xintercept = 0), lty = 2) +
    geom_hline(aes(yintercept = 0), lty = 2) +
    geom_point(alpha = 0.6) +
    xlab("Gene 1 expression\nlog2 fold change") +
    ylab("Gene 2 expression\nlog2 fold change") +
    # scale_x_continuous(limits = c(-5, 5), expand = c(0.005, 0), oob = scales::squish) +
    # scale_y_continuous(limits = c(-5, 5), expand = c(0.005, 0), oob = scales::squish) +
    coord_fixed() +
    scale_color_manual(values = c("#008837", "#7b3294"), name = NULL) +
    guides(color = guide_legend(ncol = 1)) +
    # theme(legend.justification = c(0, 1), legend.position = c(0, 1)) +
    theme(legend.position = "bottom") +
    theme(strip.text.x = element_text(margin = margin(5,0,5,0, "pt"))) +
    NULL
  print(p)

  dev.off()


  # taking all the DE genes linked by a differential contact, is there a statistically significant enrichment of same-sign DE?
  res <- NULL
  for (ps in unique(as.data.table(both_genes_to_diffHiC_gr)$paired_signf))
  {
    dt <- as.data.table(both_genes_to_diffHiC_gr)[paired_signf == ps, ]
    setnames(dt, "seqnames", "chrom")

    # take only oriented pairs
    sel <- with(dt, bait_signf != other_end_signf | distSign > 0 | (distSign == 0 & as.character(bait_gene_id) < as.character(other_end_gene_id)))
    stopifnot(grepl("CaptureC_", id) | all(sel) | mean(sel) == 0.5)
    dt <- dt[sel]
    udt <- unique(dt[, c("bait_gene_id", "bait_log2FoldChange", "other_end_gene_id", "other_end_log2FoldChange", "same_dir")])

    # what is the fraction of upregulated genes among all DE genes?
    gs <- unique(diffHiC_ov[[id]]$TSS_subset[signf == "s", c("gene_id", "log2FoldChange")])
    p <- mean(gs$log2FoldChange > 0)
    # calculate p-value for  enrichment of same-sign DE
    res <- rbind(res, data.table(
      paired_signf = ps,
      n_gene_DC_triplets = nrow(dt),
      n_gene_pairs = nrow(udt),
      n_concordant = sum(udt$same_dir),
      cor = cor(udt$bait_log2FoldChange, udt$other_end_log2FoldChange),
      pval = ifelse(nrow(udt) > 0, binom.test(sum(udt$same_dir), nrow(udt), p * p + (1 - p) * (1 - p))$p.value, NA)
    ))
  }

  if (id == "DE_all_HiC_all")
  {
    message(id)
    print(res)

    diffHiC_both_DE_dual <- as.data.table(both_genes_to_diffHiC_gr)[paired_signf == "s s", ]
    setnames(diffHiC_both_DE_dual, "seqnames", "chrom")

    diffHiC_both_DE_pairs <- diffHiC_both_DE_dual[distSign > 0 | (distSign == 0 & as.character(bait_gene_id) < as.character(other_end_gene_id)), ]
    diffHiC_both_DE_pairs[, gene_symbol := genes$gene_symbol[match(other_end_gene_id, genes$gene_id)]]
    diffHiC_both_DE_pairs[, bait_gene_symbol := genes$gene_symbol[match(bait_gene_id, genes$gene_id)]]
    dt <- diffHiC_both_DE_pairs[, c("chrom", "start", "end", "other_end_gene_id", "gene_symbol", "other_end_log2FoldChange", "baitChr", "baitStart", "baitEnd", "bait_gene_id", "bait_gene_symbol", "bait_log2FoldChange", "same_dir", "distSign", "log2FoldChange", "padj")]
    setnames(dt, "log2FoldChange", "HiClog2FoldChange")
    setnames(dt, "other_end_gene_id", "gene_id")
    setnames(dt, "other_end_log2FoldChange", "log2FoldChange")
    write.table(dt, file = "analysis/balancer/diffHiC_both_DE_pairs.tab", sep = "\t", quote = F, row.names = F)
    print(dt)
  }
}
