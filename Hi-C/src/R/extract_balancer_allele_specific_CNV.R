library(dplyr)
library(data.table)
library(GenomicRanges)

source("src/R/functions_balancer_genes.R")
source("src/R/functions_balancer_annotations.R")


## Read (flattened) exons from GFF

exon_file = paste0("data/gff/dmel-all-filtered-r6.05.UCSC_names.flattened_exons_chr23.txt")
first_exon_file = paste0("data/gff/dmel-all-filtered-r6.05.UCSC_names.first_exon_chr23.txt")

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

  fe <- exons[, list(
    V1 = first(V1),
    V4 = if (first(V7) == "+") V4[which.min(V4 + V5 / 1e9)] else V4[which.max(V5 + V4 / 1e9)],
    V5 = if (first(V7) == "+") V5[which.min(V4 + V5 / 1e9)] else V5[which.max(V5 + V4 / 1e9)],
    V7 = first(V7)
  ), by = "gene_id"]
  EXONS <- GRanges(fe$V1, IRanges(fe$V4, fe$V5), seqinfo=genomeInfo)
  EXONS$gene_id <- fe$gene_id
  # write to disk
  write.table(as.data.table(EXONS), file=first_exon_file, row.names=F, col.names=T, sep="\t", quote=F)
}

# load from file
exons <- fread(exon_file)
EXONS <- makeGRangesFromDataFrame(exons, seqinfo = genomeInfo)
if (ncol(exons)>5)
  mcols(EXONS) = exons[,6:ncol(exons),with=F]


## Overlap exons with CNVs
ovl <- as.data.table(findOverlaps(EXONS, CNV_unreduced_gr))
# Map gene id
ovl[, gene_id := EXONS$gene_id[queryHits]]
ovl[, CNV_type := CNV_unreduced_gr$type[subjectHits]]
ovl_width <- width(pintersect(EXONS[ovl$queryHits], CNV_unreduced_gr[ovl$subjectHits]))
ovl[, overlap_width := ovl_width]
ovl <- ovl[, .(gene_id, CNV_type, overlap_width)]
# Remove double entries
ovl <- ovl[, list(overlap_width = sum(overlap_width)), by = c("gene_id", "CNV_type")]

message("\nTotal #genes affected by CNVs on chr[23][LR]: ", length(unique(ovl$gene_id)))
write.table(ovl, file = "analysis/balancer/genes_allele_specific_CNV.tab", quote=F, row.names = F, sep="\t")


## Conclusions

# For now I will consider genes to be explained by a CNV if

#   * one of their exons overlaps a CNV
#   * they are significantly ASE
#   * the log2FoldChange goes into the **right direction** from 0

# If they are significant in both adult and embryo, the criteria has to be reached in both.

FINAL = NULL
for (g in unique(go$gene_id)) {
  gg <- go[gene_id == g,]
  if (all(gg$CNV_type %in% c("DUP_bal", "DUP_manual_bal", "DEL_vrg") & 
          gg$log2FoldChange > 0) || 
      all(gg$CNV_type %in% c("DUP_vrg", "DEL_bal") &
          gg$log2FoldChange < 0) )
      FINAL = rbind(FINAL, gg)
}
FINAL

message("Out of ", sum(genes$signf == "s"), " DE genes, ", length(unique(FINAL[signf == "s"]$gene_id)), " could be plausibly explained by a CNV.")
message("Out of ", sum(genes$signf == "s" & abs(genes$log2FoldChange) > log2(1.5)), " strongly DE genes, ",
  length(unique(FINAL[signf == "s" & abs(log2FoldChange) > log2(1.5)]$gene_id)), " could be plausibly explained by a CNV.")


## Aberrant transcription start sites

# from Excel spreadsheet
aberrant <- fread("analysis/balancer/DE_genes_0.05_1.5.tsv", header = T, sep = "\t")
aberrant <- aberrant[!is.na(`start of expression`), ]
aberrant <- aberrant[!is.na(`closest know TSS`), ]
aberrant_genes <- genes[gene_id %in% aberrant$gene_id & grepl('chr[23][LR]', chrom), ]

message("\nTotal # DE genes with aberrant transcription start sites: ", nrow(aberrant_genes))
message("Total # DE genes either plausibly explained by a CNV or with aberrant transcription start sites: ", length(unique(c(aberrant_genes$gene_id, FINAL[signf == "s"]$gene_id))))


## Strong ASE signal probably due to TE insertions

# https://git.embl.de/meiers/balancer/raw/master/analyses/ase/explain_weird_genes/README.md
weird_explain <- fread("analysis/balancer/explain_weird_genes.tsv", header = T, sep = "\t")
weird <- genes[gene_id %in% weird_explain$Gene & grepl('chr[23][LR]', chrom), ]

message("\nTotal # DE genes with manually annotated TE insertions: ", nrow(weird))
message("Total # DE genes either plausibly explained by a CNV or with manually annotated TE insertions: ", length(unique(c(weird$gene_id, FINAL[signf == "s"]$gene_id))))
