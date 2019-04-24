library(data.table)
library(ggplot2)
library(scales)
library(assertthat)
library(GenomicRanges)


# parameters
max_dist = 50



# All reads mapped to dm6
d = NULL
for (file in Sys.glob("mapped_to_dm6/*.txt.gz")) {
    message("Reading ", file)
    dd = fread(paste("zcat", file))
    colnames(dd) = c("id","flag","chr","pos")
    dd$file_dm6 = sub('[.]txt[.]gz$', '', sub('^.*/','',file))
    dd$sample_dm6 = sub('^.*(vrg|cross).*$', '\\1', file)
    d = rbind(d, dd)
}
setkey(d, id)


# All reads mapped to TE
e = NULL
for (file in Sys.glob("mapped_to_TE/*.txt")) {
  ee = fread(file, sep="\t", col.names = "id")
  ee$file_te = file
  ee$sample_te = tolower(gsub('^.*(VRG|CROSS).*$','\\1',file))
  e = rbind(e,ee)
}
setkey(e, id)


# Keep only reads mapped to dm6 & TE and remove weird chromosomes
g = d[e, on="id", nomatch=0]
setkey(g, chr, pos)
assert_that(all(g$sample_te == g$sample_dm6))


# clean up
g = g[grepl('^chr[23X][LR]?$', chr), .(id, chr, pos, file = file_dm6, sample = sample_dm6)]

# plot as histogram along the chrom.
message("Plot clusters along the chromosome (20kb res; chrom_overview.pdf)")
plt_chrom <- ggplot(g) +
    aes(pos, fill = sample) +
    geom_histogram(binwidth=20e3) +
    facet_grid(chr ~ ., scales = "free_y") +
    scale_x_continuous(breaks = pretty_breaks(10), labels = comma) +
    scale_y_continuous(labels = comma)
ggsave("fig.1_chrom_overview.pdf", plot=plt_chrom, width=14, height=6)



# Cluster locally
message("Chaining reads if less than ", max_dist, "bp apart...")
g[, diff := c(0,diff(pos)), by = chr]
g[, together := cumsum(ifelse(diff < max_dist, 0, 1)), by = chr]
setkey(g, chr, together)

# Summarize by cluster
h = g[, .(start = min(pos),
          end   = max(pos),
          vrg   = sum(sample=="vrg"),
          cross = sum(sample=="cross"),
          total = .N),
      by = .(chr,together)]



# Remove singletons
message("most clsuters are singletons:")
table(h$total == 1)
h = h[total > 1]


# Some overview plots for the clusters
ggsave("fig.2a_cluster.numReads.pdf", plot = ggplot(h)+aes(total)+geom_histogram(binwidth=0.1)+scale_x_log10()+ggtitle("Reads per cluster. Singletons removed"), width=8,height=4)
ggsave("fig.2b_cluster.size.pdf", plot = ggplot(h)+aes(end-start, fill=total<10)+geom_histogram(binwidth=0.1)+scale_x_log10()+ggtitle("Cluster size. Singletons removed"), width=8,height=4)



# Overlap with TE annotation
t = fread("flybase/te.bed", sep="\t")
colnames(t) = c("chr","start","end","strand","name")

h_ = makeGRangesFromDataFrame(h[, .(seqnames=chr, start=start, end=end, total)], keep.extra.columns=T)
t_ = makeGRangesFromDataFrame(t, keep.extra.columns=T)
ovl = findOverlaps(h_,t_)

h$ovl = "none"
h$ovl[queryHits(ovl)] = t$name[subjectHits(ovl)]


# Plot overlap with TE dataset
min_N_dt = data.table(min_N = c(1:10, 20, 30, 40, 50, 70, 90, 110, 150, 200, 250, 300, 400, 500, 750, 1000))
min_N_dt = min_N_dt[, h[total >= min_N, .(te = sum(ovl!="none"), total = .N)], by = min_N]
min_N_dt[, min_N := factor(min_N)]
plt = ggplot(min_N_dt) +
    geom_bar(aes(min_N,total), stat="identity",fill="grey") +
    geom_bar(aes(min_N,te),stat="identity", fill="dodgerblue") +
    geom_text(aes(x=min_N,y=te*1.15,label=te),size=3,angle = 90,hjust=0,col="dodgerblue3") +
    geom_text(aes(x=min_N,y=49.9e3,label=total),size=3,angle = 90,hjust=1,col="black") +
    coord_cartesian(ylim=c(0,50e3))
ggsave("fig.3a_te_overlap_vs_min_N.pdf",plot=plt,width=10,height=6)


max_size = data.table(max_size = c(20e3,10e3,5e3,4e3,3e3,2e3,1e3,800,600,400,200,100,90,80,70,60,50,40,30,25,20,15,10,5,4,3,2,1))
max_size = max_size[, h[end-start <= max_size, .(te = sum(ovl!="none"), total = .N)], by = max_size]
max_size[, max_size := factor(max_size)]
plt = ggplot(max_size) +
    geom_bar(aes(max_size,total), stat="identity",fill="grey") +
    geom_bar(aes(max_size,te),stat="identity", fill="dodgerblue") +
    geom_text(aes(x=max_size,y=te*1.15,label=te),size=3,angle = 90,hjust=0,col="dodgerblue3") +
    geom_text(aes(x=max_size,y=49.9e3,label=total),size=3,angle = 90,hjust=1,col="black")
ggsave("fig.3b_te_overlap_vs_max_size.pdf",plot=plt,width=10,height=6)



# # Find difference using DESeq - doest work due to coverage bring too low
# library(DESeq2)
# de = merge(h, h[, .(locus = paste0(chr[1], ":", min(start), "-", max(end))), by = .(chr, together)])
# de[, c("chr","start","end","together") := NULL]
#
# # Preparing counts and colData
# counts = dcast(de, locus ~ factor(file, levels = unique(de$file)), value.var = "N", drop = c(TRUE,TRUE))
# counts = counts[complete.cases(counts)]
# assert_that(all(unique(de$file) == colnames(counts)[2:ncol(counts)]))
#
# count_matrix = as.matrix(counts[, 2:ncol(counts)])
# rownames(count_matrix) = counts$locus
# colnames(count_matrix) = unique(de$file)
# coldata = unique(de[, .(file, sample)])[order(file),]
#
# # Running analysis
# deseq = DESeqDataSetFromMatrix(countData = count_matrix, colData = coldata, design = ~ 1 + sample)
# deseq = DESeq(deseq)
# deseq_r = results(deseq, contrast = c("sample", "cross", "vrg"))
# deseq_r = deseq_r[!is.na(deseq_r$pvalue),]
#
# plt_pvalue = ggplot(as.data.table(deseq_r)) + aes(pvalue) + geom_histogram(binwidth=0.01)
# ggsave("pvalues.pdf", plot = plt_pvalue, width = 6, height = 3)






# Only keep clusters with N >= 50 and max 2kb length!
H = h[total >= 70 & end-start<2e3]
plt_chrom <- ggplot(H) +
    aes((end+start)/2) +
    geom_histogram(binwidth=20e3) +
    facet_grid(chr ~ ., scales = "free_y") +
    scale_x_continuous(breaks = pretty_breaks(10), labels = comma) +
    scale_y_continuous(labels = comma) +
    ggtitle(paste0("Cluter distribution (N = ", nrow(H), ")"))
ggsave("fig.4_chrom_overview_after_filter.pdf", plot=plt_chrom, width=14, height=6)





# binomial test
meanP = H[, median(vrg/total)]
meanP = 1/3
H[, p_binom := mapply(function(k,n){binom.test(k,n,meanP)$p.value}, vrg, total)]
H[, padj := p.adjust(p_binom)]





## Compare to few (33) known site
known = fread("4.list_of_known_MEIs.txt")
colnames(known) = c("chrom","start","end")
ovl = findOverlaps(makeGRangesFromDataFrame(H),
                   makeGRangesFromDataFrame(known))
H$ovl_true_MEI = F
H$ovl_true_MEI[queryHits(ovl)] = T
H[, significant := factor(padj < 0.01, levels=c(T,F), labels = c("sig", "not"))]


# Final set
FINAL = H[significant == "sig"]
FINAL[, genotype := ifelse(vrg/total < 0.33, "balancer", "wild type")]
message(nrow(I[ovl_true_MEI == T & significant=="sig"]), " of ", nrow(known), " known MEIs are present in the final set (N = ", nrow(FINAL),")")


plt_sample_support2 <- ggplot(H) +
    aes(log10(total), log2(vrg/cross), col = factor(padj < 0.01 , levels = c(F,T), labels = c("not","sign")) ) +
    geom_point(alpha = 0.6, size = 0.1) +
    geom_hline(yintercept = log2(1/2), linetype = "dashed", color = "darkblue") +
    scale_y_continuous(breaks = -10:10) +
    scale_x_continuous(breaks = c(1,2,3,4,5), labels = c("10","100","1,000","10,000","100,000")) +
    scale_color_manual(values = c(sign = "dodgerblue", not = "grey"),
        name = "Binomial test",
        labels = c(sign = "significant (5% FDR)", not = "not significant")) +
    theme(legend.position = "bottom") +
    xlab("total reads at locus") +
    guides(colour = guide_legend(override.aes = list(size=4))) +
    geom_point(data = I[ovl_true_MEI==T], aes(log10(total), log2(vrg/cross)), col = "darkorange", size = 1)
ggsave("fig.5a_final.vrg_vs_cross.pdf", plot = plt_sample_support2, width=6, height = 6)


# Some overview plots for the clusters
message("Plot clusters along the chromosome (20kb res; final.chrom_overview.pdf)")
plt_chrom <- ggplot(FINAL) +
    aes((start+end)/2, fill=genotype) +
    geom_histogram(binwidth=20e3) +
    facet_grid(chr ~ ., scales = "free_y") +
    scale_x_continuous(breaks = pretty_breaks(10), labels = comma) +
    scale_y_continuous(labels = comma) +
    ggtitle(paste0("Final cluter distribution (N = ", nrow(FINAL), ")"))
ggsave("fig.5b_chrom_overview_final.pdf", plot=plt_chrom, width=14, height=6)

write.table(H, file = "final.predicted_MEI.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(FINAL[, .(chr, start, end, genotype)], file = "final.predicted_MEI.bed", col.names = F, row.names=F, quote=F, sep="\t")
