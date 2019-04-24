library(ggplot2)
library(data.table)
library(GenomicRanges)
source("../../figures/scripts/colors.R")


# Loading data
dgrp = fread("hglft_genome_52e6_1d0620.bed")
colnames(dgrp) = c("chrom","start","end")
dgrp <- dgrp[grepl('chr[23][LR]', chrom)]

ours = fread("../../SV_filter/FINAL/DUP.final.bed")
colnames(ours) = c("chrom","start","end","gt")
ours <- ours[grepl('chr[23][LR]', chrom)]

# Tidy up:
ours[, gt:= substring(gt,1,7)]
ours[gt == "1/1_1/1"]$gt = "common"
ours[gt == "0/1_0/0"]$gt = "balancer"
ours[gt == "0/1_1/1"]$gt = "wildtype"
ours[gt == "0/1_0/1"]$gt = "wildtype"

# Plot size distribution
d = rbind(cbind(dgrp, source="DGRP"),
          ours[, .(chrom, start, end, source = "F1")])
ggplot(d[end-start>=100]) + 
    aes(x=end-start, y=..density.., fill=`source`) + 
    geom_histogram(binwidth=0.1, position=position_dodge(), alpha=1) +
    theme_minimal() +
    scale_x_log10(breaks = c(10,100,1000,10000,100000), labels=c("10bp","100bp","1kb","10kb","100kb"),
                  minor_breaks = c(50, 500, 5000, 50000, 500000)) +
    #scale_y_log10(breaks = c(1,2,11,101,1001), labels=c("None",1,10,100,1000),
    #              minor_breaks = c(6,51,501)) +
    ggtitle("DUP size distribution (100+ bp)") +
    theme(legend.position = c(0.85,0.75), legend.background = element_rect(fill="white", size=0.2)) +
    xlab("SV size")
ggsave("DUP_sizes.pdf", width = 4, height = 3)



# Copied from `/g/furlong/project/33_Hi-C/src/R/functions_balancer_annotations.R` 12.11.2018
subset_overlap_percent <- function(query, subject, percentage, reciprocal=F)
{
  subject <- reduce(subject)
  ov <- findOverlaps(query, subject)
  ov_dt <- data.table(query = queryHits(ov),
                      query_width = width(query[queryHits(ov)]),
                      subject_width = width(subject[subjectHits(ov)]),
                      overlap_width = width(pintersect(query[queryHits(ov)], subject[subjectHits(ov)])))
  ov_dt <- ov_dt[, list(overlap_width = sum(overlap_width)), by = c("query", "query_width","subject_width")]
  if (!reciprocal) {
    sel <- ov_dt[overlap_width / query_width >= percentage]$query
  } else {
    sel <- ov_dt[overlap_width / query_width >= percentage & overlap_width / subject_width >= percentage]$query
  }
  return(query[sel])
}

# FIlter DGRP DELs of at least 10bp (50% of 19bp)
dgrp_gr = makeGRangesFromDataFrame(dgrp[end-start>208])
bal_gr  = makeGRangesFromDataFrame(ours, keep.extra.columns = T)


# overlap DGRP vs. balancer
rec_proc_ovl = subset_overlap_percent(bal_gr, dgrp_gr, 0.5, reciprocal = T)
message("Reciprocal overlaps found = ", length(rec_proc_ovl))

d = data.table(chrom = as.character(seqnames(rec_proc_ovl)),
               start = start(rec_proc_ovl),
               end = end(rec_proc_ovl),
               gt = rec_proc_ovl$gt)
message("Numbers rec. overlap:")
d[order(end-start), .(chrom = substr(chrom,1,4), gt, size = paste(round((end-start)/100)/10, "kb"))]

