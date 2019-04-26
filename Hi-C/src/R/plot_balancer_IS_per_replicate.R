require(cowplot)
require(ggplot2)
require(GGally)
require(MASS)
source("src/R/functions_IS.R")

theme_set(theme_cowplot(font_size = 11)) # reduce default font size
ts <- theme_get()$plot.subtitle
ts$hjust <- 0.5
theme_update(plot.subtitle = ts) # , legend.title = theme_get()$legend.text
theme_update(strip.text.x = element_text(size = 11), strip.text.y = element_text(size = 11))
theme_update(
  strip.text.x = element_text(margin = margin(t = 11 / 2, b = 11 / 2), size = 11),
  strip.text.y = element_text(margin = margin(l = 11 / 2, r = 11 / 2), size = 11)
)


window_size <- 195e3

# read IS values and annotate "boundary_gr" accordingly
read_IS("HiC_DB_6-8h_combined_VRG", ".VRG")
read_IS("HiC_DB_6-8h_combined_VRGdownBAL", ".VRGdownBAL")
read_IS("HiC_DB_6-8h_combined_BAL", ".BAL", genome = "dm6bal3")

read_IS("HiC_DB_6-8h_R1_VRG", ".VRG_R1")
read_IS("HiC_DB_6-8h_R1_VRGdownBAL", ".VRGdownBAL_R1")
read_IS("HiC_DB_6-8h_R1_BAL", ".BAL_R1", genome = "dm6bal3")

read_IS("HiC_DB_6-8h_R2_VRG", ".VRG_R2")
read_IS("HiC_DB_6-8h_R2_VRGdownBAL", ".VRGdownBAL_R2")
read_IS("HiC_DB_6-8h_R2_BAL", ".BAL_R2", genome = "dm6bal3")

#
#  overlap with breakpoints
#

bp <- fread("analysis/breakpoints_dm6.tab", header = T, sep = "\t")
bp <- bp[!grepl("start$|end$", bp$id), ]
bp_gr <- GRanges(bp$chrom, IRanges(bp$breakpoint, bp$breakpoint - 1L))
bp_gr <- bp_gr + window_size
boundary_gr$affected_by_breakpoint <- overlapsAny(boundary_gr, bp_gr)

#
#  clean up data structures
#

b <- as.data.table(boundary_gr)
setnames(b, "seqnames", "chrom")
b <- b[grepl("^chr[23]", chrom), ]

b[, start := NULL]
b[, end := NULL]
b[, width := NULL]
b[, strand := NULL]
setnames(b, "interval_start", "start")
setnames(b, "interval_end", "end")
b[, start := start - 1L]
b[, midpoint := midpoint - 1L]

#
# calculate correlations
#

# bt <- b[!is.na(IS.VRG) & !is.na(IS.BAL), ]
# res <- data.frame(included = "all", count = nrow(bt),
#   cor.Pearson = cor(bt$IS.VRG, bt$IS.BAL), cor.Spearman = cor(bt$IS.VRG, bt$IS.BAL, method = "spearman"))

# bt <- b[!is.na(IS.VRG) & !is.na(IS.BAL) & !b$affected_by_breakpoint, ]
# res <- rbind(res, data.frame(included = "not affected by breakpoint", count = nrow(bt),
#   cor.Pearson = cor(bt$IS.VRG, bt$IS.BAL), cor.Spearman = cor(bt$IS.VRG, bt$IS.BAL, method = "spearman")))

# bt <- b[!is.na(IS.VRGdownBAL) & !is.na(IS.BAL), ]
# res <- rbind(res, data.frame(included = "wild-type downsampled, all", count = nrow(bt),
#   cor.Pearson = cor(bt$IS.VRGdownBAL, bt$IS.BAL), cor.Spearman = cor(bt$IS.VRGdownBAL, bt$IS.BAL, method = "spearman")))

# bt <- b[!is.na(IS.VRGdownBAL) & !is.na(IS.BAL) & !b$affected_by_breakpoint, ]
# res <- rbind(res, data.frame(included = "wild-type downsampled, not affected by breakpoint", count = nrow(bt),
#   cor.Pearson = cor(bt$IS.VRGdownBAL, bt$IS.BAL), cor.Spearman = cor(bt$IS.VRGdownBAL, bt$IS.BAL, method = "spearman")))

# bt <- b[!is.na(IS.VRGdownBAL) & !is.na(IS.VRG), ]
# res <- rbind(res, data.frame(included = "wild-type downsampled vs. full, all", count = nrow(bt),
#   cor.Pearson = cor(bt$IS.VRGdownBAL, bt$IS.VRG), cor.Spearman = cor(bt$IS.VRGdownBAL, bt$IS.VRG, method = "spearman")))

# bt <- b[!is.na(IS.VRGdownBAL) & !is.na(IS.VRG) & !b$affected_by_breakpoint, ]
# res <- rbind(res, data.frame(included = "wild-type downsampled vs. full, not affected by breakpoint", count = nrow(bt),
#   cor.Pearson = cor(bt$IS.VRGdownBAL, bt$IS.VRG), cor.Spearman = cor(bt$IS.VRGdownBAL, bt$IS.VRG, method = "spearman")))

# print(res)

#
#  plot everything
#

# b$class <- ifelse(b$affected_by_breakpoint, "affected by breakpoint", "not affected by breakpoint")
b$class <- factor(ifelse(b$affected_by_breakpoint, "+", "-"), c("+", "-"))
n <- ncol(b)

bp <- bp[, list(breakpoint = mean(breakpoint)), by = c("chrom", "id")]


# b[, VRG_R1 := IS.VRG_R1]
# b[, VRG_R2 := IS.VRG_R2]
b[, VRG_R1 := IS.VRGdownBAL_R1]
b[, VRG_R2 := IS.VRGdownBAL_R2]
b[, BAL_R1 := IS.BAL_R1]
b[, BAL_R2 := IS.BAL_R2]
b[, shape := "dot"]


res <- NULL
cv <- c("VRG_R1", "VRG_R2", "BAL_R1", "BAL_R2")
for (col1 in cv)
  for (col2 in cv)
  {
    bt <- b[!is.na(b[[col1]]) & !is.na(b[[col2]]), ]
    res <- rbind(res,
      data.frame(col1, col2, included = "all", count = nrow(bt),
        cor.Pearson = cor(bt[[col1]], bt[[col2]]), cor.Spearman = cor(bt[[col1]], bt[[col2]], method = "spearman"))
    )
  }
print(res)


pdf("analysis/balancer/plot_IS_Fig4_xyplot_per_replicate.pdf", width = 6, height = 6)

pm <- ggpairs(b,
  aes(color = class, pch = shape, alpha = 0.6),
  columns = n + 1:4,
  # lower = list(continuous = "smooth"),
  diag = list(continuous = "blank")
)

pm <- pm
for(i in 2:pm$nrow)
{
  for(j in 1:(i-1))
  {
    pm[i,j] <- pm[i,j] +
      # geom_abline(slope = 1, lty = 2) +
      # geom_smooth(method = "lm", show.legend = F) +
      scale_x_continuous(limits = c(-1.05, 1.05)) +
      scale_y_continuous(limits = c(-1.05, 1.05)) +
  # geom_point(aes(color = class), shape = ".", pch = ".") +
  # geom_abline(slope = 1, lty = 2) +
  geom_smooth(method = "lm", show.legend = F, color = "black", lwd = 0.5) +
  # coord_fixed() +
  # labs(x = "Insulation Score (wild type)", y = "Insulation Score (balancer)") +
  scale_shape_manual(values = c(dot = ".")) +
  scale_color_manual(name = "Genomic locus", values = c("#e41a1c", "#999999")) +
  # guides(colour = guide_legend(override.aes = list(shape=16))) +
  # theme(legend.position = "bottom") #+

  # geom_point(aes(color = class), shape = ".", alpha = 0.6) +
  # geom_smooth(method = "lm", show.legend = F, color = "black") +
  # coord_fixed() +
  # labs(x = "Insulation Score (wild type)", y = "Insulation Score (balancer)") +
  # guides(color = guide_legend(override.aes = list(pch = 16), title.position = "left", ncol = 1)) +
  # theme(legend.position = "bottom") #+
  NULL

  }
}

print(pm)
dev.off()


# pdf("analysis/balancer/plot_IS_Fig4_xyplot_combined_replicates.pdf", width = 10, height = 10)

# b[, VRG := IS.VRG]
# b[, VRGdownBAL := IS.VRGdownBAL]
# b[, BAL := IS.BAL]

# pm <- ggpairs(b,
#   aes(color = class),
#   columns = n + 7:9,
#   # lower = list(continuous = "smooth"),
#   diag = list(continuous = "blank")
# )

# pm <- pm
# for(i in 2:pm$nrow)
# {
#   for(j in 1:(i-1))
#   {
#     pm[i,j] <- pm[i,j] +
#       geom_abline(slope = 1, lty = 2) +
#       # geom_smooth(method = "lm", show.legend = F) +
#       scale_x_continuous(limits = c(-1, 1)) +
#       scale_y_continuous(limits = c(-1, 1))
#   }
# }

# print(pm)
# dev.off()



#   scale_color_brewer(palette = "Set1") +
  # geom_point(aes(color = class), shape = ".") +
  # coord_fixed() 

# p <- ggplot(b, aes(IS.VRG, IS.BAL)) +
#   labs(x = "Insulation Score (wild\u00adtype)", y = "Insulation Score (balancer)") +
#   theme(legend.position = "bottom") #+
#   facet_grid(~ class)
# print(p)
