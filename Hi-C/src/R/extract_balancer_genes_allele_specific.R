library(dplyr)
library(data.table)
library(GenomicRanges)

source("src/R/functions_balancer_genes.R")
source("src/R/functions_balancer_annotations.R")

Sascha.dir <- "/g/korbel/shared/projects/drosophila_balancer"


## Read (flattened) exons from GFF

exon_file = paste0("data/gff/dmel-all-filtered-r6.05.UCSC_names.flattened_exons_chr23.txt")

genome <- "dm6"
chrom_sizes <- fread("/g/furlong/genome/D.melanogaster/Dm6/chrsizes/dm6.ucsc.chrom.sizes", col.names = c("chrom", "size"))
chrom_sizes <- chrom_sizes[order(chrom), ]
genomeInfo <- Seqinfo(seqnames = chrom_sizes$chrom, seqlengths = chrom_sizes$size, isCircular = rep(F, nrow(chrom_sizes)), genome = genome)

if (!file.exists(exon_file)) {
  # if the flattened exons have not been computed yet, do so now (takes time)
  GFF <- fread("zcat /g/furlong/genome/D.melanogaster/Dm6/6.05/gff/dmel-all-filtered-r6.05.UCSC_names.genes.gff.gz")
  exons <- GFF[V3 == "exon" & grepl('chr[23][LR]', V1) & grepl('^gene_id=', V9),]
  exons <- exons[, gene_id := substr(V9, 9, 19)]
  exons <- exons[, c("V1","V4","V5","V7","gene_id"), with = F]
  exon_gr <- GRanges(exons$V1, IRanges(exons$V4, exons$V5), strand = exons$V7)
  
  EXONS <- GRanges(seqinfo=genomeInfo)
  for (gene in unique(exons$gene_id)) {
    sub <- exons[gene_id == gene,]
    gr <- reduce(GRanges(sub$V1, IRanges(sub$V4, sub$V5), seqinfo=genomeInfo))
    gr$gene_id = gene
    EXONS <- c(EXONS, gr)
  }
  # write to disk
  write.table(as.data.table(EXONS), file=exon_file, row.names=F, col.names=T, sep="\t", quote=F)
}

# load from file
exons <- fread(exon_file)
EXONS <- makeGRangesFromDataFrame(exons, seqinfo = genomeInfo)
if (ncol(exons)>5)
  mcols(EXONS) = exons[,6:ncol(exons),with=F]


## Overlap exons with CNVs
ovl <- as.data.table(findOverlaps(makeGRangesFromDataFrame(EXONS, seqinfo = genomeInfo), CNV_gr))
# Map gene id
ovl[, gene_id := EXONS$gene_id[queryHits]]
ovl <- ovl[, .(gene_id)]
# Remove double entries
ovl <- unique(ovl)

message("\nTotal #genes affected by CNVs on chr[23][LR]: ", length(unique(ovl$gene_id)))
write.table(ovl, file = "analysis/balancer/genes_allele_specific_CNV.tab", quote=F, row.names = F, sep="\t")


## Overlap exons with MEIs
ovl <- as.data.table(findOverlaps(makeGRangesFromDataFrame(EXONS, seqinfo = genomeInfo), MEI_gr))
# Map gene id
ovl[, gene_id := EXONS$gene_id[queryHits]]
ovl <- ovl[, .(gene_id)]
# Remove double entries
ovl <- unique(ovl)

message("\nTotal #genes affected by MEIs on chr[23][LR]: ", length(unique(ovl$gene_id)))
write.table(ovl, file = "analysis/balancer/genes_allele_specific_MEI.tab", quote=F, row.names = F, sep="\t")
