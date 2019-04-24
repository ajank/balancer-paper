library(data.table)
library(ggplot2)
library(scales)
library(assertthat)
library(GenomicRanges)


# 1) get MEI callsw
mei = fread("final.predicted_MEI.bed")
colnames(mei) = c("chrom", "start", "end", "genotype")


# 2) get genes
#library(BSgenome.Dmelanogaster.UCSC.dm6)
#library(GenomicFeatures)
#txdb <- makeTxDbFromUCSC("dm6", tablename="refGene") # takes a bit of time
#genes = genes(txdb)

genes = fread("../../correlation/gene_pos/flattened_exons.txt")
genes = genes[, .(seqnames = seqnames[1], start = min(start), end = max(end), strand = strand[1]), by = gene_id]
ase = fread("../../ase/deseq/DESeq.N1_6-8h.standardFormat.txt")
genes = merge(genes, ase, by = "gene_id", all.x = T)
genes[, is_expressed := factor(!is.na(baseMean), levels = c(T,F), labels = c("expressed", "not expressed"))]
genes[, ase := ifelse(is.na(baseMean), NA, ifelse(padj<0.05, "ASE", "not ASE"))]
GENE_SPAN <- makeGRangesFromDataFrame(genes)


# 3) overlap
ovl = findOverlaps(GENE_SPAN, makeGRangesFromDataFrame(mei))
genes$hit_by_MEI = "not hit"
genes[queryHits(ovl)]$hit_by_MEI = "hit by MEI"
write.table(genes[hit_by_MEI == "hit by MEI", .(gene_id, is_expressed, ase)], file = "final.genes.txt", quote=F, sep="\t", row.names=F, col.names=T)


# SUmmary

table(genes$is_expressed, genes$hit_by_MEI)
table(genes$ase, genes$hit_by_MEI)
fisher.test(table(genes$ase, genes$hit_by_MEI))


# 4) stratify hit genes
cmbd = cbind(genes[queryHits(ovl)], mei[subjectHits(ovl), .(mei_start = start, mei_end = end, mei_gt = genotype)])
ggsave("MEIs_hitting_genes1.pdf", plot = ggplot(cmbd) + aes(mei_end - mei_start, fill = is_expressed) + geom_histogram(binwidth = 500, position=position_dodge()), width=6,height=4)
ggsave("MEIs_hitting_genes1.pdf", plot = ggplot(cmbd[is_expressed = "expressed"]) + aes(mei_end - mei_start, fill = ase) + geom_histogram(binwidth = 500, position=position_dodge()), width=6,height=4)
