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
gff <- read.table("../../SNV_annotation/reformatted.gff", header=T) %>%
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


### plotPCA of assay
plotPCA <- function(object, intgroup, ntop=1000) 
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(assay(object)[select,]))
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop=FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply( intgroup.df, 1, paste, collapse=" : "))
  } else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], group=group, 
                  intgroup.df, name=colnames(object))
  ggplot(data=d, aes_string(x="PC1", y="PC2", color="group")) + 
    geom_point(size=2) + 
    xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
    ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance"))
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
read_htseqcount <- function(f, gene_size, dir="../../gene_expression/counts/")
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
CHROM <- fread("../../../data/dm6_annotation/dm6-all-filtered-r6.05.UCSC_names.genes.gff3")[V3 == "exon" & grepl('chr[23X][LR]?$', V1) & grepl('^gene_id=', V9),]
CHROM <- CHROM[, gene_id := substr(V9, 9, 19)]
CHROM <- CHROM[, .(chrom = V1[1]), by = gene_id]
