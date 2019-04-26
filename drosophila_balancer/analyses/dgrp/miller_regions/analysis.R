library(data.table)
library(assertthat)
library(magrittr)
library(ggplot2)
source("../../figures/scripts/colors.R")

# Load data
dgrp = fread("SNPs_dgrp.txt")
wgs  = fread("SNPs_wgs.txt")

# index
setkey(dgrp, chrom, pos)
setkey(wgs, chrom, pos)

# Filter chromosomes
dgrp = dgrp[grepl('chr[23X][LR]?',chrom)]
wgs  = wgs[grepl('chr[23X][LR]?',chrom)]

# Tidy up:
wgs[gt == "1/1_1/1"]$gt = "common"
wgs[gt == "0/1_0/0"]$gt = "balancer"
wgs[gt == "0/1_1/1"]$gt = "wildtype"


# Merge
all = merge(dgrp,
            wgs,
            by = c("chrom","pos"),
            suffixes = c(".dgrp", ".ours"),
            all = T)

# Check that REF is correct, and replace columns
assert_that(all (all[, is.na(ref.dgrp) | is.na(ref.ours) | ref.ours == ref.dgrp]))
all[, `:=`(ref = ifelse(!is.na(ref.dgrp), ref.dgrp, ref.ours),
           ref.ours = NULL,
           ref.dgrp = NULL)]

# Global statistics:
# How many of our SNPs are in DGRP?
stats1 <- copy(all)
stats1[, chrom := substr(chrom,1,4)]
stats1 <- stats1[
  !is.na(gt),
  .(snps_in_total = .N,
    snps_in_dgrp = sum(!is.na(alt.dgrp)),
    snps_in_dgrp_with_matching_alt = sum(!is.na(alt.dgrp) & alt.dgrp == alt.ours)),
  by = c("chrom","gt")]
stats1[, fraction := 1-snps_in_dgrp_with_matching_alt / snps_in_total]
stats1
stats1_mean = stats1[, sum(snps_in_dgrp)/sum(snps_in_total)]

ggplot(stats1) +
  aes(gt, fraction, fill = gt) +
  geom_bar(stat="identity") +
  facet_grid(chrom ~ .) +
  ylab("Fraction of SNVs not in DGRP")+
  xlab(NULL) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = my_colors) +
  my_theme +
  coord_flip() +
  theme(legend.position = "none") +
  geom_hline(yintercept = 1-stats1_mean, alpha = 0.6, linetype = "dashed")
ggsave("fraction_of_snvs_not_in_DGRP.pdf", width=4, height=3)

# Global stats 2
# How many DGRP SNVs are covered?
stats2 <- copy(all[!is.na(af)])
stats2[, chrom := substr(chrom,1,4)]
stats2[, af_level := ifelse(af <= 0.05, ifelse(af < 0.01, "<1%", "1-5%"), 
                                        ifelse(af > 0.2,  ">20%", "5-20%"))]
stats2[, af_level := factor(af_level, levels = c("<1%","1-5%","5-20%",">20%"), ordered = T)]
stats2 <- stats2[,
  .(snps_in_total = .N,
    snps_in_samples = sum(!is.na(alt.ours)),
    snps_in_samples_with_matching_alt = sum(!is.na(alt.ours) & alt.dgrp == alt.ours)),
  by = c("chrom","af_level")]
stats2[, fraction := snps_in_samples_with_matching_alt / snps_in_total]
stats2
stats2_mean = stats2[, sum(snps_in_samples_with_matching_alt)/sum(snps_in_total)]

ggplot(stats2) +
  aes(af_level, fraction, fill = af_level) +
  geom_bar(stat="identity") +
  facet_grid(chrom ~ .) +
  ylab("Fraction of DGRP SNVs found in F1")+
  xlab(NULL) +
  scale_y_continuous(labels = scales::percent) +
  my_theme +
  coord_flip() +
  scale_fill_discrete(name = "Allele frequency") +
  theme(legend.position = "bottom") +
  geom_hline(yintercept = stats2_mean, alpha = 0.6, linetype = "dashed")
ggsave("fraction_of_DGRP_SNVs_found_in_F1.pdf", width=4, height=3)


# Position-specific effects
# Trying to make a "track" along the chromosomes showing the number of SNVs that are in DGRP
# and similar measures

# helper functions
max_coord <- function(chrom_) {
  max(dgrp[chrom==chrom_, pos])
}
min_coord <- function(chrom_) {
  min(dgrp[chrom==chrom_, pos])
}
format_Mb <- function(x) {
  paste0(round(x/100000)/10, "Mb")
}

# Draw breakpoint lines:
bps <- fread("../../SV_validation/hiC_SVs/alek.hiC_breakpoints.final.txt")


# Compute density using kernel density estimate!
get_density <- function(data, chrom_, bw = 500*1000, n = 512*16, kernel = "gaussian") {
  D <- density(data,
               from   = min_coord(chrom_),
               to     = max_coord(chrom_),
               bw     = bw,
               n      = n,
               kernel = kernel)
  E <- data.table(pos = D$x, dens = D$y * length(data))
  E
}

# Get local densities of "our" SNVs, stratified by chromosome and GT
dens.ours <- all[!is.na(gt),
    get_density(pos, chrom),
    by = .(chrom, gt)]

# what's the mean density across all chromosomes (2 & 3)
mean_snv_density = all[!is.na(gt), .N, by = chrom]
mean_snv_density[, `:=`(length = max_coord(chrom) - min_coord(chrom)), by = chrom]

ggplot(dens.ours) + 
  aes(pos, dens, col = gt) +
  geom_line() +
  ylab("SNV density") +
  facet_grid(chrom ~ .) +
  my_theme +
  ggtitle("Local SNV density (500 Mb Gaussian kernel)") +
  geom_hline(yintercept = mean_snv_density[, sum(N)/sum(length)/3],
             alpha = 0.5,
             linetype = "dashed") +
  scale_color_discrete(name = NULL) +
  theme(legend.position = "bottom") +
  xlab("Genomic coordinates") +
  scale_x_continuous(breaks = scales::pretty_breaks(6), labels = format_Mb) +
  geom_vline(data = bps[breakpoint2!=1], aes(xintercept=breakpoint2), color="blue", size=0.33, alpha=0.6)
ggsave("local_snv_density.pdf", width=8, height=6)



# Also get the local densities of our SNVs that are also in DGRP
dens.ours_in_dgrp <- all[!is.na(gt) & !is.na(alt.dgrp) & alt.dgrp == alt.ours,
    get_density(pos, chrom),
    by = .(chrom, gt)]

# Calculate (density) fraction of SNVs that are not in DGRP
dens <- merge(dens.ours,
              dens.ours_in_dgrp,
              by = c("chrom","pos","gt"),
              suffixes = c(".ours", ".ours_in_dgrp"))
dens[, fraction_snvs_not_in_dgrp := 1-dens.ours_in_dgrp / dens.ours]


ggplot(dens) + 
  aes(pos, fraction_snvs_not_in_dgrp, col = gt) +
  geom_line() +
  ylab("Fraction") +
  facet_grid(chrom ~ .) +
  geom_hline(yintercept = 1-stats1_mean, alpha = 0.5, linetype = "dashed") +
  my_theme +
  scale_color_discrete(name = NULL) +
  coord_cartesian(ylim = c(0,0.4)) +
  ggtitle("Local fractions of SNVs not in DGRP (500 Mb Gaussian kernel)") +
  theme(legend.position = "bottom") +
  xlab("Genomic coordinates") +
  scale_x_continuous(breaks = scales::pretty_breaks(6), labels = format_Mb) +
  scale_y_continuous(breaks = scales::pretty_breaks(3), labels = function(x) {scales::percent(x, accuracy=1)}) +
  geom_vline(data = bps[breakpoint2!=1], aes(xintercept=breakpoint2), color="blue", size=0.33, alpha=0.6)
ggsave("local_fraction_snvs_not_in_dgrp.pdf", width=8, height=6)
