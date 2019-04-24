#library(knitr)
library(ggplot2)
library(data.table)
library(scales)
#library("DESeq2")
#library(tidyr)
#library(reshape2)
library(dplyr)
#library(fdrtool)
#library(genefilter) # rowVars
#library(assertthat)
#library("pheatmap")
#library("RColorBrewer")
library(GenomicRanges)


### Read DE(ASE) information on genes

ase_path="../../ase/deseq/"
ASE.adult = fread(paste0(ase_path, "DESeq.F1_heads.standardFormat.txt"))

ASE.embryo = fread(paste0(ase_path, "DESeq.N1_6-8h.standardFormat.txt"))



### Read gene span and exons

#Note: this will take a while


exon_file = "flattened_exons.txt"
genomeInfo <- Seqinfo(c("chr2L", "chr2R" ,"chr3L", "chr3R"), 
                      genome="dm6", 
                      seqlengths = c(23513712, 25286936, 28110227, 32079331))



if (!file.exists(exon_file)) {
  # if the flattened exons have not been computed yet, do so now (takes time)
  message("GFF file is being flattened ...")
  GFF <- fread("../../../data/dm6_annotation/dmel-all-filtered-r6.05.UCSC_names.genes.gff3")
  exons <- GFF[V3 == "exon" & grepl('chr[23][LR]', V1) & grepl('^gene_id=', V9),]
  exons <- exons[, gene_id := substr(V9, 9, 19)]
  exons <- exons[, c("V1","V4","V5","V7","gene_id"), with = F]
  EXONS <- GRanges(seqinfo=genomeInfo)
  for (gene in unique(exons$gene_id)) {
    sub <- exons[gene_id == gene,]
    gr <- reduce(GRanges(sub$V1, IRanges(sub$V4, sub$V5), 
                         strand = sub$V7[1], seqinfo=genomeInfo))
    gr$gene_id = gene
    EXONS <- c(EXONS, gr)
  }
  # write to disk
  write.table(as.data.table(EXONS), file=exon_file, row.names=F, col.names=T, sep="\t", quote=F)
  return (NULL)
}

# load from file
exons <- fread(exon_file)
EXONS <- makeGRangesFromDataFrame(exons, seqinfo = genomeInfo)
if (ncol(exons)>5)
  mcols(EXONS) = exons[,6:ncol(exons),with=F]

GENE_SPAN <- as.data.frame(EXONS) %>%
  group_by(gene_id) %>% 
  summarise(seqnames = first(seqnames), 
            start = min(start), 
            end = max(end), 
            strand = first(strand)) %>% 
  as.data.table
GENE_SPAN <- GRanges(GENE_SPAN$seqnames, IRanges(GENE_SPAN$start, GENE_SPAN$end), 
                     strand = GENE_SPAN$strand, gene_id = GENE_SPAN$gene_id)

GENE_TSS <- GENE_SPAN
end(GENE_TSS[strand(GENE_TSS) == "+",])   = start(GENE_TSS[strand(GENE_TSS) == "+",])
start(GENE_TSS[strand(GENE_TSS) == "-",]) = end(GENE_TSS[strand(GENE_TSS) == "-",])


### Functions to determine number of ASE genes in a region

get_genes_around <- function(gene_span, points, span) {
  suppressWarnings(x <- trim(points + span))
  return(get_genes_in_region(gene_span, x))
}
get_genes_in_region <- function(gene_span, intvl) {
  out <-NULL
  for(i in 1:length(intvl)) {
    genes <- unique(subsetByOverlaps(gene_span, intvl[i])$gene_id)
    adult.expr = ASE.adult[gene_id %in% genes,]
    embryo.expr = ASE.embryo[gene_id %in% genes,]
    n_tested.ad <- nrow(adult.expr)
    n_tested.em <- nrow(embryo.expr)
    n_ase.ad    <- nrow(adult.expr[padj < 0.05,])
    n_ase.em    <- nrow(embryo.expr[padj < 0.05,])
    out <- rbind(out, 
                 data.table(n = length(genes), 
                            #ids = paste(genes,collapse=","),
                            n_tested.ad = n_tested.ad,
                            n_tested.em = n_tested.em,
                            n_ase.ad = n_ase.ad,
                            n_ase.em = n_ase.em,
                            r_ase.ad = n_ase.ad/(n_tested.ad+0.01),
                            r_ase.em = n_ase.em/(n_tested.em+0.01) ))
  }
  mcols(intvl) <- cbind(mcols(intvl), out)
  intvl
}
get_random_breakpoints <- function(n, GR) {
  # Expects non-overlapping GRs
  interval = sample(1:length(GR), n, replace=T, prob = width(GR)/sum(width(GR)))
  chr <- as.character(seqnames(GR))[interval]
  pos <- as.integer(runif(n, 
                          min = start(GR)[interval], 
                          max = end(GR)[interval]) )
  return(GRanges(chr, IRanges(pos, pos), seqinfo = genomeInfo))
}



subsetByOverlapsDT <- function(dt1, dt2, f = 0, add_frac=F)
{
  gr1 <- GRanges(dt1[[1]], IRanges(dt1[[2]], dt1[[3]]), seqinfo = genomeInfo)
  gr2 <- GRanges(dt2[[1]], IRanges(dt2[[2]], dt2[[3]]), seqinfo = genomeInfo)
  
  ovl <- as.data.table(findOverlaps(gr1,gr2))
  d <- dt1[ovl$queryHits, ]
  e <- dt2[ovl$subjectHits, ]
  
  if (f > 0) {
    frac <- (pmin(d[[3]], e[[3]]) - pmax(d[[2]], e[[2]]) + 1) / (d[[3]] - d[[2]] + 1)
    if (add_frac) {
      return(cbind(d[frac >= f,], frac = frac[frac >= f]))
    } else {
      return(d[frac >= f,])
    }
  } else {
    return(d)
  }
}



## Theme

my_theme1 <- theme_minimal() +
  theme(legend.position="bottom",
        legend.title=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = "gray30", color = NA, size = 2),
        strip.text = element_text(colour = "white"))

my_theme2 <- theme_minimal() +
  theme(legend.position="bottom",
        legend.title=element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(fill = "gray30", color = NA, size = 2),
        strip.text = element_text(colour = "white"))
