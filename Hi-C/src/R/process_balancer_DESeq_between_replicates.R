# === common.R ===

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(fdrtool))
suppressPackageStartupMessages(library(genefilter)) # rowVars
suppressPackageStartupMessages(library(assertthat))
suppressPackageStartupMessages(library("pheatmap"))
suppressPackageStartupMessages(library("RColorBrewer"))

# Variables:
gff <- read.table("/g/korbel/shared/projects/drosophila_balancer/analyses/SNV_annotation/reformatted.gff", header=T) %>%
  group_by(gene_id) %>%
  summarize(chrom = first(chrom),
            from = min(start),
            to = max(end),
            strand = first(strand),
            length = sum(end - start +1))


# Functions:

### Take an DESeq assay, turn it into long format
### and annotate column information
assay_long <- function(DEassay, colInfo = colData(DE)) 
{
  colDF <- data.frame(colInfo, sample = rownames(colInfo))
  data.frame(DEassay, 
             gene_id = row.names(DEassay), 
             row.names = NULL) %>% 
    gather(sample, count, -gene_id) %>%
    merge(., colDF, by = "sample", all.x = T)
}



### Give a list of gene_ids for which we should plot 
### haplotype counts from the assay
geneCountData <- function(gene_ids, dataSet = DE, intgroup=c("Haplotype","GT", "Replicate")) 
{
  get <- function(g) {
    data.frame(gene = g, row.names = NULL,
               plotCounts(dataSet, gene=g, returnData=TRUE,
                          intgroup = intgroup))
  }
  Reduce( rbind, lapply(gene_ids, FUN = get), init=NULL)
}


### Read total gene expression data for annotation
read_htseqcount <- function(f, gene_size, dir="/g/korbel/shared/projects/drosophila_balancer/analyses/gene_expression/counts/")
{
  expr_lvl <- read.table(paste(dir,f,sep="")) %>%
    dplyr::rename(gene_id = V1, count = V2) %>%
    dplyr::filter(!grepl("^__", gene_id))
  assert_that(noNA(expr_lvl))
  total_read_count = sum(expr_lvl$count)
  expr_lvl <- suppressWarnings(dplyr::left_join(expr_lvl, gene_size, by="gene_id")) %>%
    mutate(FPKM = count*1000000/total_read_count*1000/length) %>%
    select(gene_id, count, FPKM) %>% 
    filter(!is.na(FPKM))
  assert_that(noNA(expr_lvl))
  expr_lvl
}


### Map each gene to a chromosome
CHROM <- fread("/g/korbel/shared/projects/drosophila_balancer/data/dm6_annotation/dm6-all-filtered-r6.05.UCSC_names.genes.gff3")[V3 == "exon" & grepl('chr[23X][LR]?$', V1) & grepl('^gene_id=', V9),]
CHROM <- CHROM[, gene_id := substr(V9, 9, 19)]
CHROM <- CHROM[, .(chrom = V1[1]), by = gene_id]


# === DESeq.Rmd ===

library(data.table)
library(cowplot)
library(scales)

COVPERSAMPLE = 50
ALPHA = 0.05
MINLFC = log2(1.5)


filter_genes <- function(DE, thr_per_sample = COVPERSAMPLE, chrom_regex = 'chr[32][LR]') {
  DE <- DE[ rownames(DE) %in% filter(gff, grepl(chrom_regex, chrom))$gene_id, ]
  DE <- DE[ rowSums(counts(DE)) >= ncol(DE)/2*thr_per_sample, ]
  return (DE)
}

fdr_correction <- function(R) {
  corr = fdrtool(R$stat, plot = F)
  p <- rbind(data.frame(p = R$padj,    t = "before correction"),
      data.frame(p = corr$pval, t = "after correction")) %>%
    ggplot() + aes(p) + theme_minimal() +
      geom_histogram(binwidth=0.1) + scale_x_log10() +
      facet_grid(.~t, scales="free")
  #print(p)
  R$pvalue = corr$pval
  R$padj = corr$qval
  return(R)
}

my_summary <- function(R, min_lfc = MINLFC, alpha = ALPHA) {
  dt <- as.data.table(R)
  n_sig = nrow(dt[padj<ALPHA,])
  n_sig_lfc = nrow(dt[padj<ALPHA & abs(log2FoldChange)>MINLFC,])
  summary(R, alpha = alpha)
  cat(paste0("Total DE genes:              ", n_sig,"\n"))
  cat(paste0("Total DE genes with FC>", round(2**(MINLFC),1), ": ", n_sig_lfc,"\n"))
  return(NULL)
}

write_deseq_table <- function(R, file) {
  dt <- cbind(gene_id = rownames(R), as.data.table(R)) 
  write.table(dt, file, quote = F, row.names = F, sep = "\t")
}


# Embryo data (4 replicates)

n=c("SampleName",              "File",                                         "Haplotype",  "Sample", "Replicate")
x=c("N1_6.8_father.vrg_1_BAL", "N1_pool_6-8h_rep1.htseq-count-rev.alt.txt",    "balancer",   "N1pat", "1",
    "N1_6.8_father.vrg_1_VRG", "N1_pool_6-8h_rep1.htseq-count-rev.ref.txt",    "virginizer", "N1pat", "1",
    "N1_6.8_father.vrg_2_BAL", "N1_pool_6-8h_rep2.htseq-count-rev.alt.txt",    "balancer",   "N1pat", "2",
    "N1_6.8_father.vrg_2_VRG", "N1_pool_6-8h_rep2.htseq-count-rev.ref.txt",    "virginizer", "N1pat", "2",
    "N1_6.8_mother.vrg_1_BAL", "N1sex_pool_6-8h_rep1.htseq-count-rev.alt.txt", "balancer",   "N1mat", "3",
    "N1_6.8_mother.vrg_1_VRG", "N1sex_pool_6-8h_rep1.htseq-count-rev.ref.txt", "virginizer", "N1mat", "3",
    "N1_6.8_mother.vrg_2_BAL", "N1sex_pool_6-8h_rep2.htseq-count-rev.alt.txt", "balancer",   "N1mat", "4",
    "N1_6.8_mother.vrg_2_VRG", "N1sex_pool_6-8h_rep2.htseq-count-rev.ref.txt", "virginizer", "N1mat", "4" )
D_e = data.frame(matrix(x, byrow=T,ncol=5))
colnames(D_e) = n

DE_e <- DESeqDataSetFromHTSeqCount(D_e, directory = "/g/korbel/shared/projects/drosophila_balancer/analyses/readSeparation/counts", design = ~ Replicate + Haplotype)
DE_e <- filter_genes(DE_e)
A_e = DESeq(DE_e)
R_e = results(A_e, contrast=c("Haplotype", "balancer", "virginizer"))

R_e = fdr_correction(R_e)
my_summary(R_e, alpha = ALPHA)

# plot_MA(R_e)

# write_deseq_table(R_e, file = "DESeq.N1_6-8h.standardFormat.txt")

# Let's check the chromosomal distribution (in %):

xx = merge(data.table(as.data.table(R_e), gene_id = rownames(R_e)), CHROM, by="gene_id")
round(table(xx[padj<0.05,]$chrom) / table(xx$chrom), 4)*100
xx = xx[, .(total=.N, ase=sum(padj<0.05),ratio=sum(padj<0.05)/.N), by=chrom]
xx


# Embryo data (comparing N1pat replicates in the same haplotype)

for (hapl in c("balancer", "virginizer"))
{
  Dbal_e <- D_e[D_e$Haplotype == hapl & D_e$Sample == "N1pat", ]
  Dbal_e$Haplotype <- NULL
  Dbal_e$Sample <- NULL
  DE_e <- DESeqDataSetFromHTSeqCount(Dbal_e, directory = "/g/korbel/shared/projects/drosophila_balancer/analyses/readSeparation/counts", design = ~ Replicate)
  DE_e <- filter_genes(DE_e)
  A_e = DESeq(DE_e)
  R_e = results(A_e, contrast=c("Replicate", "2", "1"))

  R_e = fdr_correction(R_e)
  my_summary(R_e, alpha = ALPHA)
  write_deseq_table(R_e, file = paste0("analysis/balancer/DESeq.between_replicates.N1mat_", hapl, "_6-8h.standardFormat.txt"))
}

# Embryo data (comparing N1pat haplotypes in the same replicate)

for (repl in c("1", "2"))
{
  Dbal_e <- D_e[D_e$Replicate == repl & D_e$Sample == "N1pat", ]
  Dbal_e$Replicate <- NULL
  Dbal_e$Sample <- NULL
  DE_e <- DESeqDataSetFromHTSeqCount(Dbal_e, directory = "/g/korbel/shared/projects/drosophila_balancer/analyses/readSeparation/counts", design = ~ Haplotype)
  DE_e <- filter_genes(DE_e)
  A_e = DESeq(DE_e)
  R_e = results(A_e, contrast=c("Haplotype", "balancer", "virginizer"))

  R_e = fdr_correction(R_e)
  my_summary(R_e, alpha = ALPHA)
  write_deseq_table(R_e, file = paste0("analysis/balancer/DESeq.between_haplotypes.N1mat_Rep", repl, "_6-8h.standardFormat.txt"))
}
