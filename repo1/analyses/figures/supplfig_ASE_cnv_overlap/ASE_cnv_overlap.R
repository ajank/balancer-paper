library(data.table)
library(ggplot2)
library(assertthat)
suppressPackageStartupMessages(library(GenomicRanges))

# Read data
ase = fread(snakemake@input[["ase"]])
gff = fread(snakemake@input[["gff"]])
del = fread(snakemake@input[["dels"]])[, .(chrom = V1, start = V2, end = V3, GT = V4)]
dup = fread(snakemake@input[["dups"]])[, .(chrom = V1, start = V2, end = V3, GT = V4)]
assert_that("gene_id"  %in% colnames(ase),
            "baseMean" %in% colnames(ase),
            "log2FoldChange" %in% colnames(ase),
            "padj"     %in% colnames(ase))
assert_that("chrom"    %in% colnames(gff),
            "start"    %in% colnames(gff),
            "end"      %in% colnames(gff),
            "gene_id"  %in% colnames(gff),
            "feature"  %in% colnames(gff))

# Remove a few wrong entries from this list and then assert that all parts of a gene are
# exactly consecutive
gff <- gff[!(gene_id %in% c("FBgn0020309","FBgn0032378","FBgn0264815"))]
assert_that(all(gff[, max(end) - min(start) +1, by = gene_id] == gff[, sum(end - start+1), by = gene_id]))

# Use only exons in the gff file!
gff <- gff[feature %in% c('five_prime_UTR','CDS','three_prime_UTR','exon')]

# Subset everything to chr2 + chr3
gff <- gff[grepl('^chr[23][LR]$', chrom)]
ase <- merge(ase, gff[, .(chrom = chrom[1]), by = gene_id], by = "gene_id")
del <- del[grepl('^chr[23][LR]$', chrom)]
dup <- dup[grepl('^chr[23][LR]$', chrom)]
genomeInfo <- Seqinfo(c("chr2L", "chr2R" ,"chr3L", "chr3R"), genome="dm6", 
                      seqlengths = c(23513712, 25286936, 28110227, 32079331))



# Beautify del and dup
del[grepl('0/1_0/0', GT)]$GT = "bal-spec"
del[grepl('0/1_1/1', GT)]$GT = "wt-spec"
del[grepl('1/1_1/1', GT)]$GT = "common"
assert_that(all(sort(unique(del$GT)) == c("bal-spec","common","wt-spec")))

dup[grepl('0/1_0/0', GT)]$GT = "bal-spec"
dup[grepl('0/1_1/1|0/1_0/1', GT)]$GT = "wt-spec"
assert_that(all(sort(unique(dup$GT)) == c("bal-spec","wt-spec")))


# Merge del and dup data
cnv = rbind(cbind(del, sv = "del"), cbind(dup, sv = "dup"))
# Remove non-haplotype-specific calls
cnv = cnv[GT!="common",]


# Find overlap between genes and CNVs
gff.gr = makeGRangesFromDataFrame(gff, keep.extra.columns = T, seqinfo = genomeInfo)
cnv.gr = makeGRangesFromDataFrame(cnv, keep.extra.columns = T, seqinfo = genomeInfo)
ovl <- as.data.table(findOverlaps(cnv.gr,gff.gr))

# Annotate this overlap set
ovl$sv = cnv.gr[ovl$queryHits]$sv
ovl$gt = cnv.gr[ovl$queryHits]$GT
ovl$gene_id = gff.gr[ovl$subjectHits]$gene_id
ovl$start   = start(gff.gr[ovl$subjectHits])
ovl$end     = end(gff.gr[ovl$subjectHits])

# Merge by gene and type of SV
ovl <- ovl[, .(size = sum(end-start+1)), by = .(gene_id, sv, gt)]

# Clean up genes with more than 1 cnv overlapping
remove_genes <- ovl[duplicated(ovl$gene_id)]$gene_id
ovl[gene_id %in% remove_genes]$sv = "multiple"
ovl[gene_id %in% remove_genes]$gt = ""
ovl <- unique(ovl[, .(gene_id,sv,gt)])



# Some statistics:
data <- merge(ovl, ase, by = "gene_id")
n_tot = nrow(ase)
n_ase = nrow(ase[padj<0.05])
n_ovl = nrow(data)
n_ovl_ase = nrow(data[padj<0.05])

#message("Number of genes with CNV overlap:           ", nrow(ovl))
message("Number of expressed genes:                   ", n_tot)
message("Number of ASE genes:                         ", n_ase)
message("Number of expressed genes with CNV overlap:  ", n_ovl)
message("Number of expressed ASE genes with CNV ovlp: ", n_ovl_ase, "\t(", n_ovl_ase/n_ase, "% of ASE genes)")

# Testing this ratio         
x = matrix(c(n_tot - n_ase - (n_ovl - n_ovl_ase),    n_ovl - n_ovl_ase,
             n_ase - n_ovl_ase,                      n_ovl_ase         ), 
    byrow = T, nrow=2)
colnames(x) = c("no_CNV", "CNV_overlap")
rownames(x) = c("no_ASE", "ASE_genes")
message("Testing this matrix for randomness:")
print(x)
fisher.test(x)


# setup plot
title = paste(nrow(data[padj<0.05]), "ASE genes whith at least one exon affected by CNVs")
ggplot(data[padj < 0.05]) +
  aes(paste(sv,gt), log2FoldChange, fill = paste(sv,gt)) + 
  geom_boxplot(outlier.alpha = 0, alpha = 0.6) + 
  geom_jitter(size = 1, aes(color = paste(sv,gt))) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  guides(fill = FALSE, col = FALSE) +
  theme_bw() +
  ggtitle(title)
ggsave(snakemake@output[[1]], width=6, height = 3)

