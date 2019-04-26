suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(GenomicRanges))
source("scripts/colors.R")



### Input & parameters
f_out = snakemake@output[[1]]
f_del = snakemake@input[["dels"]] # "../SV_filter/FINAL/DEL.final.bed"
f_dhs = snakemake@input[["peaks"]] # "../tracks/enhancer/DNase_HS_sites_stages9-11_HotSpot_peaks_FDR-1pc_liftedToDm6.bed"
genome = GRanges(c("chr2L", "chr2R" ,"chr3L", "chr3R"), 
                 IRanges(c(0,0,0,0),
                         c(23513712, 25286936, 28110227, 32079331)))





# Copied from `/g/furlong/project/33_Hi-C/src/R/functions_balancer_annotations.R` 12.11.2018
subset_overlap_percent <- function(query, subject, percentage)
{
  subject <- reduce(subject)
  ov <- findOverlaps(query, subject)
  ov_dt <- data.table(query = queryHits(ov), query_width = width(query[queryHits(ov)]),
                      overlap_width = width(pintersect(query[queryHits(ov)], subject[subjectHits(ov)])))
  ov_dt <- ov_dt[, list(overlap_width = sum(overlap_width)), by = c("query", "query_width")]
  sel <- ov_dt[overlap_width / query_width >= percentage]$query
  return(query[sel])
}


# Source: https://github.com/leekgroup/enrichedRanges/blob/master/R/randomInterval.R
randomInterval <- function(I, n = 1, ms = 10000, strand.specific = FALSE, randomize.strand = FALSE) 
{
  stopifnot(is(I, "GRanges"))
  
  seqs <- width(I)
  sample_sequence <- function(m, seqs = seqs) {
    ps <- pmax(0, as.numeric(seqs - m + 1))
    psum <- sum(ps)
    if (psum == 0) {
      stop(paste("no sequence is long enough to sample an interval of length ", 
                 m, sep = ""))
    }
    prob <- ps / psum
    sample(seq_len(length(seqs)), size = 1, prob = prob)
  }
  ## vector ms
  xs <- sapply(cbind(1:n, ms)[, 2], sample_sequence, seqs = seqs)
  
  starts <- mapply(function(min, max) {
    as.integer(round(runif(n = 1, min = min, max = max)))
  }, start(I)[xs], end(I)[xs])
  
  if(strand.specific) {
    if(randomize.strand) {
      strand <- sample(c("+", "-"), size = length(starts), replace = TRUE)
    } else {
      strand <- strand(I)[xs]
    }
    
  } else {
    strand <- "*"
  }
  
  seqnames <- seqnames(I)[xs]
  gr <- GRanges(seqnames = seqnames, ranges = IRanges(starts, starts + 
                                                        ms - 1), strand = strand, seqlengths = seqlengths(I)[unique(seqnames)] )
  
  return(gr)
}



DEL_new = fread(f_del)
DEL_new = DEL_new[, .(chrom = V1, 
                      start = V2 + 1,
                      end = V3,
                      GT = substr(V4,1,7))]

# remove chr X
DEL_new <- DEL_new[grepl('^chr[23][LR]', chrom)]

DEL_all = with(DEL_new, GRanges(chrom, IRanges(start,end)))
DEL_BAL = with(DEL_new[GT=="0/1_0/0",], GRanges(chrom, IRanges(start,end)))
DEL_VRG = with(DEL_new[GT=="0/1_1/1",], GRanges(chrom, IRanges(start,end)))



DHS = with(fread(f_dhs),
           GRanges(V1, IRanges(V2,V3)))
DHS <- DHS[grepl('^chr[23][LR]', seqnames(DHS))]
start(DHS) <- start(DHS) + 1


### Starting analysis
overlap_thr = 0.05
iterations  = 500


# Compile data table with overlap counts
data = NULL
data = rbind(data, data.table(type = "real",
                              number = length(subset_overlap_percent(DHS, DEL_BAL, overlap_thr)),
                              Genotype = "Balancer"))
data = rbind(data, data.table(type = "real",
                              number = length(subset_overlap_percent(DHS, DEL_VRG, overlap_thr)),
                              Genotype = "Wild type"))


# Sample background distributions
cat("Sampling background for balancer ")
for (i in 1:iterations) {
  DEL_rand = randomInterval(genome, length(DEL_BAL), width(DEL_BAL))
  data = rbind(data, data.table(type="rand",
                                number = length(subset_overlap_percent(DHS, DEL_rand, overlap_thr)),
                                Genotype = "Balancer"))
  cat(".")
}
cat("\n")
cat("Sampling background for wild type ")
for (i in 1:iterations) {
  DEL_rand = randomInterval(genome, length(DEL_VRG), width(DEL_VRG))
  data = rbind(data, data.table(type="rand",
                                number = length(subset_overlap_percent(DHS, DEL_rand, overlap_thr)),
                                Genotype = "Wild type"))
  cat(".")
}
cat("\n")
head(data)

pvals = data[, .(pval = (frank(.SD, number)[type == "real"])/.N), by = Genotype]

plt <- ggplot(data[type == "rand",]) + 
  aes(number, fill = Genotype) + 
  geom_histogram(binwidth = 5) + 
  geom_vline(data = data[type=="real",], aes(xintercept = number), color="darkorange") +
  facet_wrap( ~ Genotype, ncol=1) +
  geom_label(data = pvals, aes(label = paste0("hat(p)==",round(pval,3))), 
             x=Inf, y=Inf, hjust=1, vjust=1, inherit.aes=F, parse=T) +
  theme_minimal() +
  xlab("Number of DHS peaks >=5% deleted") +
  scale_fill_manual(values = c(`Wild type` = as.character(my_colors["wildtype"]),
                               `Balancer`  = as.character(my_colors["balancer"]),
                               `common`    = as.character(my_colors["common"]))) +
  guides(fill=FALSE)
ggsave(plt, filename = f_out, width = 4, height = 3)
