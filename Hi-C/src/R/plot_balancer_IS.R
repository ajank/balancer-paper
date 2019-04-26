require(cowplot)
require(ggplot2)
require(MASS)
source("src/R/functions_IS.R")

theme_set(theme_cowplot(font_size = 11)) # reduce default font size
ts <- theme_get()$plot.subtitle
ts$hjust <- 0.5
theme_update(plot.subtitle = ts) # , legend.title = theme_get()$legend.text
theme_update(strip.text.x = element_text(size = 11), strip.text.y = element_text(size = 11))


window_size <- 195e3

# read IS values and annotate "boundary_gr" accordingly
read_IS("HiC_DB_6-8h_combined_VRG", source = "filtered_5000", ".VRG")
# read_IS("HiC_DB_6-8h_combined_VRGdownBAL", source = "filtered_5000", ".VRGdownBAL")
read_IS("HiC_DB_6-8h_combined_BAL", ".BAL", source = "filtered_5000", genome = "dm6bal3")

# read TAD boundaries
b.VRG <- read_boundaries("HiC_DB_6-8h_combined_VRG")
b.VRG <- b.VRG[grepl("^chr[23]", chrom), ]
b.BAL <- read_boundaries("HiC_DB_6-8h_combined_BAL", genome = "dm6bal3")
b.BAL <- b.BAL[grepl("^chr[23]", chrom), ]

# identify unique TAD boundaries
boundary_equiv_dist <- 5e3
b.VRG_gr <- with(b.VRG, GRanges(chrom, IRanges(midpoint - boundary_equiv_dist / 2 - 1, midpoint + boundary_equiv_dist / 2)))
b.BAL_gr <- with(b.BAL, GRanges(chrom, IRanges(midpoint - boundary_equiv_dist / 2 - 1, midpoint + boundary_equiv_dist / 2)))
ov <- findOverlaps(b.VRG_gr, b.BAL_gr, minoverlap = 0)

b.VRG[, unique := T]
b.VRG$unique[unique(queryHits(ov))] <- F
b.BAL[, unique := T]
b.BAL$unique[unique(subjectHits(ov))] <- F


#{"step":10000,"minDepth":50000,"maxDepth":200000,"binsize":5000}
# 7 columns with bin sizes: 50k 60k 80k 100k 130k 160k 195k
# > 50e3 + (0:6)^1.5 * 10e3
# INFO:hicFindTADs:computing spectrum for window sizes between 10 (50000 bp)and 40 (200000 bp) at the following window sizes 2 [50000, 60000, 78284, 101961, 130000, 161803, 196969]

#
#  overlap with breakpoints
#

bp <- fread("analysis/breakpoints_dm6.tab", header = T, sep = "\t")
bp <- bp[!grepl("start$|end$", bp$id), ]
bp_gr <- GRanges(bp$chrom, IRanges(bp$breakpoint, bp$breakpoint - 1L))
bp_gr <- bp_gr + window_size
boundary_gr$affected_by_breakpoint <- overlapsAny(boundary_gr, bp_gr)

#
#  identify differential IS regions
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

lm.fit <- lm(IS.BAL ~ IS.VRG, b)
b[, IS.BAL.predicted.from.VRG := predict(lm.fit, b)]
b[, diff.IS := IS.BAL - IS.BAL.predicted.from.VRG]

low.thr <- quantile(b$diff.IS, 0.05, na.rm = T)
upp.thr <- quantile(b$diff.IS, 0.95, na.rm = T)
mid.thr <- quantile(abs(b$diff.IS), 0.5, na.rm = T)

b[, diff.IS.class := "other"]
b[, diff.IS.class := ifelse(diff.IS < low.thr, "lower_5percent", diff.IS.class)]
b[, diff.IS.class := ifelse(-mid.thr < diff.IS & diff.IS < mid.thr, "middle_50percent", diff.IS.class)]
b[, diff.IS.class := ifelse(upp.thr < diff.IS, "upper_5percent", diff.IS.class)]
b[, diff.IS.class := factor(diff.IS.class)]

print(b)
write.table(b, file = "analysis/balancer/differential_IS.tab", sep = "\t", quote = F, row.names = F, col.names = T)
save(b, lm.fit, low.thr, upp.thr, mid.thr, file = "analysis/balancer/differential_IS.Rdata")

#
#  do the plotting
#

bl <- rbind(
  data.table(b, value = b$IS.VRG, allele = "w"),
  # data.table(b, value = b$IS.VRGdownBAL, allele = "wd"),
  data.table(b, value = b$IS.BAL, allele = "b")
)
bp <- bp[, list(breakpoint = mean(breakpoint)), by = c("chrom", "id")]

source("src/R/functions_DE_colors.R")
is.colors = c(
  # w = as.character(my_colors["wildtype"]), 
  # b = as.character(my_colors["balancer"]),
  w = "#66c2a5",
  wd = "#000000",
  b = "#fc8d62", 
  n = as.character(my_colors["common"]))
is.labels = c(w = "wild\u00adtype", wd = "wild\u00adtype downsampled", b = "balancer")


pdf("analysis/balancer/plot_IS_Fig4_all_chrom.pdf", width = 10, height = 5)

# bl[grepl("^chr[23]", chrom), ]
p <- ggplot(bl, aes(x = midpoint / 1e3, y = value, color = allele)) + 
  geom_vline(data = bp, aes(xintercept = breakpoint / 1e3), color = "#984ea3", alpha = 1, lty = "dashed") +
  geom_line() +
  geom_point(data = b.VRG, aes(y = -1.2), color = is.colors["w"], shape = 124, size = 2) +
  geom_point(data = b.BAL, aes(y = 1.2), color = is.colors["b"], shape = 124, size = 2) +
  xlab("Genomic position (kb)") +
  ylab("Insulation Score") +
  facet_grid(chrom ~ ., switch = "y") +
  background_grid(major = "xy", minor = "x") +
  coord_cartesian(ylim = c(-1.2, 1.2)) +
  scale_x_continuous(expand = c(0.01, 0)) +
  scale_y_continuous(breaks = -1:1) +
  scale_color_manual(values = is.colors, labels = is.labels, name = NULL) +
  theme(legend.position = c(0.85, 0.9), legend.background = element_rect(fill = "white"))
print(p)

dev.off()


pdf("analysis/balancer/plot_IS_Fig4_fragment_chr2L.pdf", width = 10, height = 1.5)

p <- ggplot(bl[chrom == "chr2L" & midpoint > 9e6 & midpoint <= 13e6],
  aes(x = midpoint / 1e3, y = value, color = allele)) + 
  geom_vline(data = bp[chrom == "chr2L" & breakpoint > 9e6 & breakpoint <= 13e6],
    aes(xintercept = breakpoint / 1e3), color = "#984ea3", alpha = 1, lty = "dashed") +
  geom_line() +
  geom_point(data = b.VRG[chrom == "chr2L" & midpoint > 9e6 & midpoint <= 13e6],
    aes(y = -1.2), color = is.colors["w"], shape = 124, size = 2) +
  geom_point(data = b.BAL[chrom == "chr2L" & midpoint > 9e6 & midpoint <= 13e6],
    aes(y = 1.2), color = is.colors["b"], shape = 124, size = 2) +
  geom_point(data = b.VRG[chrom == "chr2L" & midpoint > 9e6 & midpoint <= 13e6 & unique == T],
    aes(y = -1), color = is.colors["w"], shape = 124, size = 2) +
  geom_point(data = b.BAL[chrom == "chr2L" & midpoint > 9e6 & midpoint <= 13e6 & unique == T],
    aes(y = 1), color = is.colors["b"], shape = 124, size = 2) +
  xlab("Genomic position (kb)") +
  ylab("Insulation\nScore") +
  facet_grid(chrom ~ ., switch = "y") +
  background_grid(major = "xy", minor = "x") +
  coord_cartesian(ylim = c(-1.2, 1.2)) +
  scale_x_continuous(expand = c(0.01, 0)) +
  scale_y_continuous(breaks = -1:1) +
  scale_color_manual(values = is.colors, labels = is.labels, name = NULL)
print(p)

p <- ggplot(bl[chrom == "chr2L" & midpoint > 11e6 & midpoint <= 12e6],
  aes(x = midpoint / 1e3, y = value, color = allele)) + 
  geom_vline(data = bp[chrom == "chr2L" & breakpoint > 11e6 & breakpoint <= 12e6],
    aes(xintercept = breakpoint / 1e3), color = "#984ea3", alpha = 1, lty = "dashed") +
  geom_line() +
  geom_point(data = b.VRG[chrom == "chr2L" & midpoint > 11e6 & midpoint <= 12e6],
    aes(y = -1.2), color = is.colors["w"], shape = 124, size = 2) +
  geom_point(data = b.BAL[chrom == "chr2L" & midpoint > 11e6 & midpoint <= 12e6],
    aes(y = 1.2), color = is.colors["b"], shape = 124, size = 2) +
  geom_point(data = b.VRG[chrom == "chr2L" & midpoint > 11e6 & midpoint <= 12e6 & unique == T],
    aes(y = -1), color = is.colors["w"], shape = 124, size = 2) +
  geom_point(data = b.BAL[chrom == "chr2L" & midpoint > 11e6 & midpoint <= 12e6 & unique == T],
    aes(y = 1), color = is.colors["b"], shape = 124, size = 2) +
  xlab("Genomic position (kb)") +
  ylab("Insulation\nScore") +
  facet_grid(chrom ~ ., switch = "y") +
  background_grid(major = "xy", minor = "x") +
  coord_cartesian(ylim = c(-1.2, 1.2)) +
  scale_x_continuous(expand = c(0.01, 0)) +
  scale_y_continuous(breaks = -1:1) +
  scale_color_manual(values = is.colors, labels = is.labels, name = NULL)
print(p)

dev.off()

#
# calculate correlations
#

bt <- b[!is.na(IS.VRG) & !is.na(IS.BAL), ]
res <- data.frame(included = "all", count = nrow(bt),
  cor.Pearson = cor(bt$IS.VRG, bt$IS.BAL), cor.Spearman = cor(bt$IS.VRG, bt$IS.BAL, method = "spearman"))

bt <- b[!is.na(IS.VRG) & !is.na(IS.BAL) & !b$affected_by_breakpoint, ]
res <- rbind(res, data.frame(included = "far from breakpoints", count = nrow(bt),
  cor.Pearson = cor(bt$IS.VRG, bt$IS.BAL), cor.Spearman = cor(bt$IS.VRG, bt$IS.BAL, method = "spearman")))

# bt <- b[!is.na(IS.VRGdownBAL) & !is.na(IS.BAL), ]
# res <- rbind(res, data.frame(included = "wild-type downsampled, all", count = nrow(bt),
#   cor.Pearson = cor(bt$IS.VRGdownBAL, bt$IS.BAL), cor.Spearman = cor(bt$IS.VRGdownBAL, bt$IS.BAL, method = "spearman")))

# bt <- b[!is.na(IS.VRGdownBAL) & !is.na(IS.BAL) & !b$affected_by_breakpoint, ]
# res <- rbind(res, data.frame(included = "wild-type downsampled, far from breakpoints", count = nrow(bt),
#   cor.Pearson = cor(bt$IS.VRGdownBAL, bt$IS.BAL), cor.Spearman = cor(bt$IS.VRGdownBAL, bt$IS.BAL, method = "spearman")))

# bt <- b[!is.na(IS.VRGdownBAL) & !is.na(IS.VRG), ]
# res <- rbind(res, data.frame(included = "wild-type downsampled vs. full, all", count = nrow(bt),
#   cor.Pearson = cor(bt$IS.VRGdownBAL, bt$IS.VRG), cor.Spearman = cor(bt$IS.VRGdownBAL, bt$IS.VRG, method = "spearman")))

# bt <- b[!is.na(IS.VRGdownBAL) & !is.na(IS.VRG) & !b$affected_by_breakpoint, ]
# res <- rbind(res, data.frame(included = "wild-type downsampled vs. full, far from breakpoints", count = nrow(bt),
#   cor.Pearson = cor(bt$IS.VRGdownBAL, bt$IS.VRG), cor.Spearman = cor(bt$IS.VRGdownBAL, bt$IS.VRG, method = "spearman")))

print(res)

#
#  plot everything
#

b$class <- factor(ifelse(b$affected_by_breakpoint, "close to a breakpoint", "far from breakpoints"),
  c("close to a breakpoint", "far from breakpoints"))

pdf("analysis/balancer/plot_IS_Fig4_xyplot.pdf", width = 3, height = 3.5)

p <- ggplot(b, aes(IS.VRG, IS.BAL)) +
  geom_point(aes(color = class), shape = ".", alpha = 0.6) +
  # geom_abline(slope = 1, lty = 2) +
  geom_smooth(method = "lm", show.legend = F, color = "black") +
  coord_fixed() +
  labs(x = "Insulation Score (wild type)", y = "Insulation Score (balancer)") +
  scale_color_manual(name = "Genomic locus", values = c("#e41a1c", "#999999")) +
  guides(color = guide_legend(override.aes = list(pch = 16), title.position = "left", ncol = 1)) +
  theme(legend.position = "bottom") #+
  # facet_grid(~ class)
print(p)

dev.off()

pdf("analysis/balancer/plot_IS_Fig4_xyplot_only_not_affected.pdf", width = 3, height = 3.5)

p <- ggplot(b[class == "far from breakpoints"], aes(IS.VRG, IS.BAL)) +
  geom_point(aes(color = class), shape = ".", alpha = 0.6) +
  # geom_abline(slope = 1, lty = 2) +
  geom_smooth(method = "lm", show.legend = F, color = "black") +
  coord_fixed() +
  labs(x = "Insulation Score (wild type)", y = "Insulation Score (balancer)") +
  scale_color_manual(name = "Genomic locus", values = c("#999999")) +
  guides(color = guide_legend(override.aes = list(pch = 16), title.position = "left", ncol = 1)) +
  theme(legend.position = "bottom") #+
  # facet_grid(~ class)
print(p)

dev.off()

# p <- ggplot(b, aes(IS.VRGdownBAL, IS.BAL)) +
#   geom_point(aes(color = class), shape = ".", alpha = 0.6) +
#   geom_abline(slope = 1, lty = 2) +
#   geom_smooth(method = "lm", aes(color = class), show.legend = F) +
#   coord_fixed() +
#   labs(x = "Insulation Score (wild\u00adtype downsampled)", y = "Insulation Score (balancer)") +
#   scale_color_manual(name = "Genomic locus", values = c("#e41a1c", "#999999")) +
#   guides(color = guide_legend(override.aes = list(pch = 16), title.position = "left", ncol = 1)) +
#   theme(legend.position = "bottom") +
#   facet_grid(~ class)
# print(p)

# p <- ggplot(b, aes(IS.VRGdownBAL, IS.VRG)) +
#   geom_point(aes(color = class), shape = ".", alpha = 0.6) +
#   geom_abline(slope = 1, lty = 2) +
#   geom_smooth(method = "lm", aes(color = class), show.legend = F) +
#   coord_fixed() +
#   labs(x = "Insulation Score (wild\u00adtype downsampled)", y = "Insulation Score (wild type)") +
#   scale_color_manual(name = "Genomic locus", values = c("#e41a1c", "#999999")) +
#   guides(color = guide_legend(override.aes = list(pch = 16), title.position = "left", ncol = 1)) +
#   theme(legend.position = "bottom") +
#   facet_grid(~ class)
# print(p)


# norm.fit <- fitdistr(b$diff.IS[!is.na(b$diff.IS)], "normal")
# norm.para <- norm.fit$estimate

pdf("analysis/balancer/plot_differential_IS_Fig4_xyplot.pdf", width = 3, height = 3.5)

p <- ggplot(b, aes(IS.BAL.predicted.from.VRG, IS.BAL)) +
  geom_point(aes(color = class), shape = ".", alpha = 0.6) +
  geom_abline(slope = 1, lty = 2) +
  # geom_smooth(method = "lm") +
  coord_fixed() +
  labs(x = "Insulation Score (wild type)\nscaled to match balancer", y = "Insulation Score (balancer)") +
  scale_color_manual(name = "Genomic locus", values = c("#e41a1c", "#999999")) +
  guides(color = guide_legend(override.aes = list(pch = 16), title.position = "left", ncol = 1))
print(p)

dev.off()


pdf("analysis/balancer/plot_differential_IS_Fig4_hist.pdf", width = 6, height = 4)

# norm.fit.dt <- data.table(diff.IS = seq(-1, 1, 0.02))
# norm.fit.dt[, y := dnorm(diff.IS, norm.para[1], norm.para[2]) * 0.02 * sum(!is.na(b$diff.IS))]

dt <- data.table(ymin = -Inf, ymax = Inf)
p <- ggplot(b, aes(diff.IS)) +
  geom_histogram(binwidth = 0.02, color = "#333333", fill = "white") +
  geom_rect(data = dt, aes(ymin = ymin, ymax = ymax), xmin = -Inf, xmax = low.thr, inherit.aes = F, fill = is.colors["w"], alpha = 0.2) +
  geom_rect(data = dt, aes(ymin = ymin, ymax = ymax), xmin = upp.thr, xmax = Inf, inherit.aes = F, fill = is.colors["b"], alpha = 0.2) +
  geom_rect(data = dt, aes(ymin = ymin, ymax = ymax), xmin = -mid.thr, xmax = mid.thr, inherit.aes = F, fill = is.colors["n"], alpha = 0.4) +
  # geom_line(aes(y = y), data = norm.fit.dt, color = "red") +
  xlim(c(-1, 1)) +
  xlab(expression(paste(Delta, "IS"))) +
  ylab("Count")
print(p)

dev.off()


pdf("analysis/balancer/plot_differential_IS_Fig4_all_chrom.pdf", width = 10, height = 5)

p <- ggplot(b, aes(x = midpoint / 1e3, y = diff.IS)) + 
  geom_vline(data = bp, aes(xintercept = breakpoint / 1e3), color = "#984ea3", alpha = 1, lty = "dashed") +
  geom_line(color = "#333333") +
  geom_point(data = b[diff.IS.class == "upper_5percent", ], aes(y = -1.2), color = is.colors["w"], shape = 124, size = 2) +
  geom_point(data = b[diff.IS.class == "lower_5percent", ], aes(y = -1.2), color = is.colors["b"], shape = 124, size = 2) +
  # geom_point(data = b.VRG[unique == T, ], aes(y = 1.2), color = is.colors["w"], shape = 124, size = 2) +
  # geom_point(data = b.BAL[unique == T, ], aes(y = 1.2), color = is.colors["b"], shape = 124, size = 2) +
  xlab("Genomic position (kb)") +
  ylab(expression(paste(Delta, "IS"))) +
  facet_grid(chrom ~ ., switch = "y") +
  background_grid(major = "xy", minor = "x") +
  coord_cartesian(ylim = c(-1.2, 1.2)) +
  scale_x_continuous(expand = c(0.01, 0)) +
  scale_y_continuous(breaks = -1:1)
print(p)

dev.off()
