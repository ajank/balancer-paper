options(warn = 1)

require(cowplot)
require(ggplot2)

source("src/R/functions_balancer_genes.R")
source("src/R/functions_balancer_annotations.R")

theme_set(theme_cowplot(font_size = 11)) # reduce default font size
ts <- theme_get()$plot.subtitle
ts$hjust <- 0.5
theme_update(plot.subtitle = ts) # , legend.title = theme_get()$legend.text


# add pseudo-boundaries at ends of chromosomes
chrom_sel <- unique(diffTADs_dt$chrom)
chrom_len <- chrom_map$length[match(chrom_sel, chrom_map$chrom)]
d1 <- data.table(chrom = chrom_sel, start = 1L, end = 0L, midpoint = 1L, class = "matched", delta = 1)
d2 <- data.table(chrom = chrom_sel, start = chrom_len + 1L, end = chrom_len, midpoint = chrom_len + 1L, class = "matched", delta = 1)
db <- rbind(
  diffTADs_dt,
  data.table(d1, allele = "w"), data.table(d1, allele = "b"),
  data.table(d2, allele = "w"), data.table(d2, allele = "b"),
  fill = T
)
setkey(db, allele, chrom, end)

# In chromosomes 2 and 3, we identified 771 TADs in the wild-type and 761 TADs in the balancer haplotypes, with a median size of 125 kb. 
message("In chromosomes 2 and 3, we identified ", nrow(db[allele == "w"]) - 4, " TADs in the wild-type and ", nrow(db[allele == "b"]) - 4, " TADs in the balancer haplotypes, with a median size of ", median(tail(db$end, -1) - head(db$end, -1)) / 1e3, " kb. ")


# look at the allele-specific boundaries from the perspective of encompassing "matched blocks" between matched boundaries
compare_as_boundaries <- function(this_allele = "w", other_allele = "b")
{
  message("comparing ", this_allele, "-specific boundaries to ", other_allele)
  db_matched <- db[allele == this_allele & class == "matched", ]
  db_matched_gr <- GRanges(db_matched[, list(start = head(start, -1), end = tail(end, -1)), by = "chrom"])

  db_as_this_gr <- with(db[allele == this_allele & class != "matched", ], GRanges(chrom, IRanges(end, end)))
  db_as_other_gr <- with(db[allele == other_allele & class != "matched", ], GRanges(chrom, IRanges(end, end)))

  db_as_hits <- rbind(
    # data.table(matched_block = seq_along(db_matched_gr)),
    data.table(matched_block = subjectHits(findOverlaps(db_as_this_gr, db_matched_gr)))[, list(as_this = .N), by = "matched_block"],
    data.table(matched_block = subjectHits(findOverlaps(db_as_other_gr, db_matched_gr)))[, list(as_other = .N), by = "matched_block"],
    fill = T
  )[, list(as_this = sum(as_this, na.rm = T), as_other = sum(as_other, na.rm = T)), by = "matched_block"]

  # matrix: number of allele-specific boundaries in the given matched block in this vs. other allele
  print(with(db_as_hits, table(as_this, as_other)))

  # how many allele-specific boundaries are lost for sure (not shifted)?
  message("summary: ", sum(with(db_as_hits, pmax(as_this - as_other, 0))), " are lost out of ", sum(with(db_as_hits, as_this)), "\n")
}

compare_as_boundaries("w", "b")
compare_as_boundaries("b", "w")


# identify "perturbed" TADs; are DE genes enriched there, comparing to the non-perturbed TADs?
compare_perturbed <- function(db, genes, delta_thr = 0, prominence_thr = 0, fc_thr = 0, diffTAD_margin = Inf, allele = NA, shuffle = F)
{
  this_allele <- allele
  db_subset <- db[class == "matched" | (delta > delta_thr & prominence > prominence_thr), ]
  if (!is.na(this_allele))
    db_subset <- db_subset[allele == this_allele, ]

  sel <- which(db_subset$class != "matched")
  n_perturbed_boundaries <- length(sel)
  if (shuffle)
    sel <- sample(which(db_subset$delta != 1), n_perturbed_boundaries) # db_subset$delta == 1 for chromosome start/end positions

  stopifnot(db_subset$chrom[sel] == db_subset$chrom[sel - 1])
  stopifnot(db_subset$chrom[sel] == db_subset$chrom[sel + 1])
  perturbed_gr <- reduce(GRanges(db_subset$chrom[sel], IRanges(
    pmax(db_subset$end[sel - 1], db_subset$end[sel] - diffTAD_margin),
    pmin(db_subset$end[sel + 1], db_subset$end[sel] + diffTAD_margin)
  )))

  sel <- which(db_subset$class == "matched" & db_subset$delta != 1) # db_subset$delta == 1 for chromosome start/end positions
  stopifnot(db_subset$chrom[sel] == db_subset$chrom[sel - 1])
  stopifnot(db_subset$chrom[sel] == db_subset$chrom[sel + 1])
  control_gr <- reduce(GRanges(db_subset$chrom[sel], IRanges(
    pmax(db_subset$end[sel - 1], db_subset$end[sel] - diffTAD_margin),
    pmin(db_subset$end[sel + 1], db_subset$end[sel] + diffTAD_margin)
  )))

  g <- genes
  g[, is_perturbed := countOverlaps(GRanges(g$chrom, IRanges(g$tss, g$tss)), perturbed_gr) > 0]
  g[, is_control := countOverlaps(GRanges(g$chrom, IRanges(g$tss, g$tss)), control_gr) > 0]
  g[, is_significant := signf == "s" & abs(log2FoldChange) > log2(fc_thr)]
  # print(wilcox.test(abs(g[is_perturbed == T, ]$log2FoldChange), abs(g[is_perturbed == F, ]$log2FoldChange)))

  perturbed_genes <- nrow(g[signf != "n", ][is_perturbed == T])
  perturbed_DE_genes <- nrow(g[signf != "n", ][is_significant == T & is_perturbed == T])
  control_DE_ratio <- nrow(g[signf != "n", ][is_significant == T & is_control == T & is_perturbed == F]) / nrow(g[signf != "n", ][is_control == T & is_perturbed == F])

  tbl <- with(g, table(is_significant, is_perturbed))
  dt <- data.table(
    allele = allele,
    delta_thr = delta_thr,
    prominence_thr = prominence_thr,
    fc_thr = fc_thr,
    diffTAD_margin = diffTAD_margin,
    perturbed_boundaries = n_perturbed_boundaries,
    perturbed_regions = length(perturbed_gr),
    perturbed_genes = perturbed_genes,
    perturbed_DE_genes = perturbed_DE_genes,
    perturbed_DE_ratio = perturbed_DE_genes / perturbed_genes,
    # control_DE_ratio = nrow(g[signf != "n", ][is_significant == T & is_perturbed == F]) / nrow(g[signf != "n", ][is_perturbed == F]),
    control_DE_ratio = control_DE_ratio,
    pval = binom.test(perturbed_DE_genes, perturbed_genes, p = control_DE_ratio)$p.value,
    q0 = qbinom(0.025, perturbed_genes, prob = control_DE_ratio) / perturbed_genes, # inner 95%
    q1 = qbinom(0.05, perturbed_genes, prob = control_DE_ratio) / perturbed_genes, # inner 90%
    q2 = qbinom(0.25, perturbed_genes, prob = control_DE_ratio) / perturbed_genes, # inner 50%
    q3 = qbinom(0.75, perturbed_genes, prob = control_DE_ratio) / perturbed_genes, # inner 50%
    q4 = qbinom(0.95, perturbed_genes, prob = control_DE_ratio) / perturbed_genes, # inner 90%
    q5 = qbinom(0.975, perturbed_genes, prob = control_DE_ratio) / perturbed_genes, # inner 95%
    # if (nrow(tbl) == 2 & ncol(tbl) == 2) fisher.test(tbl)$p.value else NA,
    shuffle = shuffle
  )

  return(list(dt = dt, g = g))
}

# res <- NULL
# for (allele in c(NA, "w", "b"))
#   for (prominence_thr in c(0, 0.1, 0.2))
#     for (delta_thr in c(0, 0.1, 0.2))
#       if (prominence_thr == 0 | delta_thr == 0)
#         for (fc_thr in c(0, 1.5))
#           for (diffTAD_margin in c(Inf, 50e3, 25e3, 10e3))
#             res <- rbind(res, compare_perturbed(db, genes, allele = allele, prominence_thr = prominence_thr, delta_thr = delta_thr, fc_thr = fc_thr, diffTAD_margin = diffTAD_margin)$dt)

res <- NULL
for (allele in c(NA))
  for (prominence_thr in c(0, 0.1))
    for (delta_thr in c(0, 0.1))
      if (prominence_thr == 0 | delta_thr == 0)
        for (fc_thr in c(0, 1.5))
          for (diffTAD_margin in c(10e3, 20e3, 50e3, 100e3, Inf))
            res <- rbind(res, compare_perturbed(db, genes, allele = allele, prominence_thr = prominence_thr, delta_thr = delta_thr, fc_thr = fc_thr, diffTAD_margin = diffTAD_margin)$dt)

print(res[delta_thr==0 & prominence_thr == 0])

write.table(res, file = "analysis/balancer/differential_TADs_vs_genes.tab", sep = "\t", quote = F, row.names = F, col.names = T)

cp <- compare_perturbed(db, genes, allele = NA, prominence_thr = 0, delta_thr = 0, fc_thr = 0, diffTAD_margin = 10e3)
if ("signfs" %in% names(genes))
  print(table(cp$g[is_perturbed == T, ]$signfs))

#
#  plot the fraction of DE genes
#

  # Rdata_file = paste0("analysis/balancer/diffTAD_perturbed_genes_fc_thr", fc_thr, ".Rdata")
  # if (file.exists(Rdata_file))
  # {
  #   message("Loading ", Rdata_file)
  #   load(Rdata_file, verbose = T)
  # }
  # else
  # {
  #   message("Calculating genes around random points and writing ", Rdata_file)

  #   dt <- NULL
  #   for (diffTAD_margin in 1:40 * 5e3)
  #   {
  #     print(diffTAD_margin)
  #     for (allele in c(NA, "w", "b"))
  #       dt <- rbind(dt, compare_perturbed(db, genes,
  #         allele = allele, prominence_thr = 0, delta_thr = 0, fc_thr = fc_thr, diffTAD_margin = diffTAD_margin)$dt)
  #     for (i in 1:10000)
  #       dt <- rbind(dt, compare_perturbed(db, genes,
  #         allele = NA, prominence_thr = 0, delta_thr = 0, fc_thr = fc_thr, diffTAD_margin = diffTAD_margin, shuffle = T)$dt)
  #   }

  #   save(dt, file = Rdata_file)
  # }

  # qdt <- dt[shuffle == T,
  #   list(
  #     q0 = quantile(perturbed_DE_ratio, 0.025), # inner 95%
  #     q1 = quantile(perturbed_DE_ratio, 0.05), # inner 90%
  #     q2 = quantile(perturbed_DE_ratio, 0.25), # inner 50%
  #     q3 = quantile(perturbed_DE_ratio, 0.75), # inner 50%
  #     q4 = quantile(perturbed_DE_ratio, 0.95), # inner 90%
  #     q5 = quantile(perturbed_DE_ratio, 0.975) # inner 95%
  #   ),
  #   by = c("allele", "delta_thr", "prominence_thr", "fc_thr", "diffTAD_margin", "shuffle")]


# main paper figure
for (fc_thr in c(0, 1.5))
{
  pdf(paste0("analysis/balancer/diffTAD_perturbed_genes_fc_thr", fc_thr, ".pdf"), width = 5.5, height = 4, onefile = T)

  for (allele in c(NA, "w", "b"))
  {
    qdt <- NULL
    for (diffTAD_margin in 1:40 * 5e3)
    {
      print(diffTAD_margin)
      qdt <- rbind(qdt, compare_perturbed(db, genes,
        allele = allele, prominence_thr = 0, delta_thr = 0, fc_thr = fc_thr, diffTAD_margin = diffTAD_margin)$dt)
    }
    stopifnot(length(unique(qdt$perturbed_boundaries)) == 1)

    p <- ggplot(qdt, aes(diffTAD_margin / 1e3)) +
      geom_ribbon(data = qdt, aes(ymin = q0, ymax = q5, fill = "95% quantile")) +
      geom_ribbon(data = qdt, aes(ymin = q1, ymax = q4, fill = "90% quantile")) +
      geom_ribbon(data = qdt, aes(ymin = q2, ymax = q3, fill = "50% quantile")) +
      scale_fill_manual(values = c(`95% quantile` = "#eeeeee",
                                   `90% quantile` = "#dddddd",
                                   `50% quantile` = "#bebebe"),
                        name = "Expectation based on\nmatched boundaries") +
      # geom_line(data = dt[shuffle == F, ], aes(y = perturbed_DE_ratio, color = allele), size=1) +
      geom_line(aes(y = perturbed_DE_ratio, color = ifelse(is.na(allele), "s", allele)), size=1) +
      scale_color_manual(values = gene.colors, labels = paste0(
        ifelse(is.na(allele), "allele", switch(allele, "b" = "balancer", "w" = "wild\u00adtype", "?")),
        "\u00adspecific boundaries (", qdt$perturbed_boundaries[1], ")"), name = NULL) +
      # scale_x_log10(breaks = c(10e3, 100e3, 1e6, 10e6), labels = c("10 kb", "100 kb","1 Mb","10 Mb")) +
      # scale_y_continuous(breaks = seq(0,1,0.1), labels = percent) +
      scale_y_continuous(limits = c(0,NA)) +
      xlab("Distance from TAD boundary (kb)") +
      ylab(paste0("Fraction of", ifelse(fc_thr == 1.5, " strongly", ""), " DE genes among testable genes")) +
      # theme_minimal() +
      background_grid(major = "xy", minor = "x") +
      guides(color = guide_legend(order = 1), fill = guide_legend(order = 2)) +
      theme(legend.position = c(1, 0), legend.justification = c(1, 0), legend.background = element_rect(fill = "white"), legend.box = "horizontal")
      NULL
    print(p)
  }

  dev.off()
}

# zoom in for the supplement
for (fc_thr in c(0, 1.5))
{
  pdf(paste0("analysis/balancer/diffTAD_perturbed_genes_zoom_fc_thr", fc_thr, ".pdf"), width = 3, height = 4, onefile = T)

  for (allele in c(NA, "w", "b"))
  {
    qdt <- NULL
    for (diffTAD_margin in 1:15 * 5e3)
    {
      print(diffTAD_margin)
      qdt <- rbind(qdt, compare_perturbed(db, genes,
        allele = allele, prominence_thr = 0, delta_thr = 0, fc_thr = fc_thr, diffTAD_margin = diffTAD_margin)$dt)
    }
    stopifnot(length(unique(qdt$perturbed_boundaries)) == 1)

    p <- ggplot(qdt, aes(diffTAD_margin / 1e3)) +
      geom_ribbon(data = qdt, aes(ymin = q0, ymax = q5, fill = "95% quantile")) +
      geom_ribbon(data = qdt, aes(ymin = q1, ymax = q4, fill = "90% quantile")) +
      geom_ribbon(data = qdt, aes(ymin = q2, ymax = q3, fill = "50% quantile")) +
      scale_fill_manual(values = c(`95% quantile` = "#eeeeee",
                                   `90% quantile` = "#dddddd",
                                   `50% quantile` = "#bebebe"),
                        name = "Expectation based on\nmatched boundaries") +
      # geom_line(data = dt[shuffle == F, ], aes(y = perturbed_DE_ratio, color = allele), size=1) +
      geom_line(aes(y = perturbed_DE_ratio, color = ifelse(is.na(allele), "s", allele)), size=1) +
      scale_color_manual(values = gene.colors, labels = paste0(
        ifelse(is.na(allele), "allele", switch(allele, "b" = "balancer", "w" = "wild\u00adtype", "?")),
        "\u00adspecific boundaries (", qdt$perturbed_boundaries[1], ")"), name = NULL) +
      # scale_x_log10(breaks = c(10e3, 100e3, 1e6, 10e6), labels = c("10 kb", "100 kb","1 Mb","10 Mb")) +
      # scale_y_continuous(breaks = seq(0,1,0.1), labels = percent) +
      scale_y_continuous(limits = c(0,NA)) +
      xlab("Distance from TAD boundary (kb)") +
      ylab(paste0("Fraction of", ifelse(fc_thr == 1.5, " strongly", ""), " DE genes\namong testable genes")) +
      # theme_minimal() +
      background_grid(major = "xy", minor = "x") +
      guides(color = guide_legend(order = 1), fill = guide_legend(order = 2)) +
      theme(legend.position = c(1, 0), legend.justification = c(1, 0), legend.background = element_rect(fill = "white"), legend.box = "horizontal")
      NULL
    print(p)
  }

  dev.off()
}

#
#  plot the numbers of DE and non-DE genes per TAD
#

pdf("analysis/balancer/diffTAD_vs_genes.pdf", width = 10, height = 3, onefile = T)

this_allele <- "w"
# # look at the allele-specific boundaries from the perspective of encompassing "matched blocks" between matched boundaries
# compare_as_boundaries <- function(this_allele = "w", other_allele = "b")

  TADs <- db[allele == this_allele, ]
  TADs_gr <- GRanges(TADs[, list(start = head(start, -1), end = tail(end, -1),
    is_perturbed = head(class, -1) != "matched" | tail(class, -1) != "matched"), by = "chrom"])

  g <- genes
  g_gr <- GRanges(g$chrom, IRanges(g$tss, g$tss))

  ov <- findOverlaps(g_gr, TADs_gr)
  dt <- data.table(TAD_id = subjectHits(ov), signf = g$signf[queryHits(ov)])
  TAD_hits <- rbind(
    expand.grid(TAD_id = seq_along(TADs_gr), signf = levels(g$signf)),
    dt[, list(count = .N), by = c("TAD_id", "signf")],
    fill = T
  )[, list(count = sum(count, na.rm = T)), by = c("TAD_id", "signf")]

  TAD_hits <- merge(TAD_hits, data.table(TAD_id = seq_along(TADs_gr), TAD_width = width(TADs_gr), is_perturbed = TADs_gr$is_perturbed))

  TAD_hits_wide <- TAD_hits[, list(TAD_width = TAD_width[1], is_perturbed = is_perturbed[1], count_s = count[signf == "s"], count_i = count[signf == "i"]), by = "TAD_id"]

  TAD_hits_wide[, ratio := ifelse(count_s + count_i == 0, Inf, count_s / (count_s + count_i))]
  TAD_hits_wide[, count_si := count_s + count_i]
  setkey(TAD_hits_wide, is_perturbed, ratio)
  TAD_hits_wide[, ratio_sort_id := seq_len(nrow(TAD_hits_wide))]
  setkey(TAD_hits_wide, is_perturbed, TAD_width)
  TAD_hits_wide[, width_sort_id := seq_len(nrow(TAD_hits_wide))]
  setkey(TAD_hits_wide, is_perturbed, count_si)
  TAD_hits_wide[, count_si_sort_id := seq_len(nrow(TAD_hits_wide))]
  TAD_hits <- merge(TAD_hits, TAD_hits_wide[, c("TAD_id", "ratio_sort_id", "width_sort_id", "count_si_sort_id")], by = "TAD_id")

boundary_color = scale_color_manual(values = c(`matched` = "grey", 
  `balancer-specific` = as.character(gene.colors["b"]),
  `wild-type-specific` = as.character(gene.colors["w"])))
boundary_fill = scale_fill_manual(values = c(`matched` = "grey", 
  `balancer-specific` = as.character(gene.colors["b"]),
  `wild-type-specific` = as.character(gene.colors["w"])))

# p <- ggplot(TAD_hits_wide, aes(count_i, count_s, color = is_perturbed)) +
#   # facet_grid(class ~ .) +
#   geom_point(size = 0.5)
#   # boundary_color
# print(p)

# p <- ggplot(TAD_hits[signf != "n"], aes(ratio_sort_id, count, fill = signf)) +
#   facet_grid(~ is_perturbed, scales = "free_x", space="free") +
#   geom_bar(stat = "identity")
#   # boundary_color
# print(p)

# p <- ggplot(TAD_hits[signf != "n"], aes(width_sort_id, count, fill = signf)) +
#   facet_grid(~ is_perturbed, scales = "free_x", space="free") +
#   geom_bar(stat = "identity")
#   # boundary_color
# print(p)

TAD_hits[, facet_label := ifelse(is_perturbed, "Perturbed TADs", "Non\u00adperturbed TADs")]

p <- ggplot(TAD_hits[signf != "n"], aes(count_si_sort_id, count, fill = signf)) +
  facet_grid(~ facet_label, scales = "free_x", space="free") +
  scale_x_continuous(expand = c(0, 0)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(NULL,
    values = c(`s` = as.character(gene.colors["s"]), `i` = as.character(gene.colors["i"]), `n` = as.character(gene.colors["n"])),
    labels = c(`s` = "DE genes", `i` = "non\u00adDE genes", `n` = "genes not tested")
  ) +
  ylab("Gene count") +
  theme(axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank()) +
  theme(legend.justification = c(0, 1), legend.position = c(0, 1))
  # boundary_color
print(p)


  g <- g15
  g_gr <- GRanges(g$chrom, IRanges(g$tss, g$tss))

  ov <- findOverlaps(g_gr, TADs_gr)
  dt <- data.table(TAD_id = subjectHits(ov), signf = g$signf[queryHits(ov)])
  TAD_hits <- rbind(
    expand.grid(TAD_id = seq_along(TADs_gr), signf = levels(g$signf)),
    dt[, list(count = .N), by = c("TAD_id", "signf")],
    fill = T
  )[, list(count = sum(count, na.rm = T)), by = c("TAD_id", "signf")]

  TAD_hits <- merge(TAD_hits, data.table(TAD_id = seq_along(TADs_gr), TAD_width = width(TADs_gr), is_perturbed = TADs_gr$is_perturbed))

  TAD_hits_wide <- TAD_hits[, list(TAD_width = TAD_width[1], is_perturbed = is_perturbed[1], count_s = count[signf == "s"], count_i = count[signf == "i"]), by = "TAD_id"]

  TAD_hits_wide[, ratio := ifelse(count_s + count_i == 0, Inf, count_s / (count_s + count_i))]
  TAD_hits_wide[, count_si := count_s + count_i]
  setkey(TAD_hits_wide, is_perturbed, ratio)
  TAD_hits_wide[, ratio_sort_id := seq_len(nrow(TAD_hits_wide))]
  setkey(TAD_hits_wide, is_perturbed, TAD_width)
  TAD_hits_wide[, width_sort_id := seq_len(nrow(TAD_hits_wide))]
  setkey(TAD_hits_wide, is_perturbed, count_si)
  TAD_hits_wide[, count_si_sort_id := seq_len(nrow(TAD_hits_wide))]
  TAD_hits <- merge(TAD_hits, TAD_hits_wide[, c("TAD_id", "ratio_sort_id", "width_sort_id", "count_si_sort_id")], by = "TAD_id")

TAD_hits[, facet_label := ifelse(is_perturbed, "Perturbed TADs", "Non\u00adperturbed TADs")]

p <- ggplot(TAD_hits[signf != "n"], aes(count_si_sort_id, count, fill = signf)) +
  facet_grid(~ facet_label, scales = "free_x", space="free") +
  scale_x_continuous(expand = c(0, 0)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(NULL,
    values = c(`s` = as.character(gene.colors["s"]), `i` = as.character(gene.colors["i"]), `n` = as.character(gene.colors["n"])),
    labels = c(`s` = "DE genes with absolute FC < 1.5", `i` = "other expressed genes", `n` = "genes not tested")
  ) +
  ylab("Gene count") +
  theme(axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank()) +
  theme(legend.justification = c(0, 1), legend.position = c(0, 1))
  # boundary_color
print(p)


dev.off()
