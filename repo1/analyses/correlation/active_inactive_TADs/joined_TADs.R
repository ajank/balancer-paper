library(dplyr)
library(data.table)
library(GenomicRanges)
library(assertthat)
library(ggplot2)


find_and_split_interval <- function(tads, chrom, pos) 
{
  assert_that(class(tads)=="GRanges")
  assert_that(is.integer(pos) && length(pos)==1)
  tad1 <- subsetByOverlaps(tads, GRanges(chrom, IRanges(pos,pos)))
  if (length(tad1) == 0) {
    x = GRanges(chrom, IRanges(pos,pos))
    return (list(l = x, r = x))
  } else {
    l <- tad1 
    r <- tad1
    end(l) = pos
    start(r) = pos
    return (list(l = l, r = r))
  }
}

### TADs + breakpoints
genomeInfo <- Seqinfo(c("chr2L", "chr2R" ,"chr3L", "chr3R"), 
                      genome="dm6", 
                      seqlengths = c(23513712, 25286936, 28110227, 32079331))

dm6_overview <- fread("/g/furlong/project/39_Balancer_Hi-C/analysis/assembly_dm6bal3.tab", header = T)
dm6_overview <- dm6_overview[, strand := ifelse(breakpoint1 < breakpoint2, "+", "-")][]
setnames(dm6_overview, c("breakpoint1", "breakpoint2"), c("start", "end"))
dm6_overview <- dm6_overview[, .(chrom, start, end, col, strand)]
tads = fread("/g/korbel/shared/projects/drosophila_balancer/analyses/TADs/data/TADs_Sexton_2012_embryo_16-18h_dm6.bed", sep="\t")[grepl('^chr[23][LR]$', V1)]
tads[grepl('HP1|PcG',V4),]$V4 = "Inactive"
tadgr <- GRanges(tads$V1, IRanges(tads$V2, tads$V3), type = tads$V4, seqinfo = genomeInfo)
###


### Gene data
ase_path="../../ase/deseq/"
ASE.embryo = fread(paste0(ase_path, "DESeq.N1_6-8h.standardFormat.txt"))
exon_file = "../gene_pos/flattened_exons.txt"
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
###



COUNTS = NULL
for (i in c(1,2,3, 5,6,7,    9, 11,12,13,14,15,16,17,18)) {
  
    one = dm6_overview[i,]
    two = dm6_overview[i+1,]
    message("Joining ", one$chrom, ":", round(one$end/1e6,1), "Mb(", one$strand, ") and ", two$chrom, ":", round(two$start/1e6,1), "Mb(", two$strand, ")")
    
    gene_ids_hit = c(subsetByOverlaps(GENE_SPAN, GRanges(one$chrom, IRanges(one$end, one$end)))$gene_id,
                     subsetByOverlaps(GENE_SPAN, GRanges(two$chrom, IRanges(two$end, two$end)))$gene_id)
    
    tad1 <- find_and_split_interval(tadgr, one$chrom, one$end)
    message("    Tad 1 is split into ", round(width(tad1$l)/1000), "kb + ", round(width(tad1$r)/1000), "kb.")
    tad2 <- find_and_split_interval(tadgr, two$chrom, two$start)
    message("    Tad 2 is split into ", round(width(tad2$l)/1000), "kb + ", round(width(tad2$r)/1000), "kb.")
    if(width(tad1$l)<=1) {
        message ("    ...skip (TAD 1)")
        next;
    }
    
    genes_origTad <- subsetByOverlaps(GENE_SPAN, c(tad1$l,tad1$r))$gene_id
    COUNTS = rbind(COUNTS, 
                   data.table(genome   = "dm6",
                              tad_size = sum(width(c(tad1$l,tad1$r))),
                              tad_id   = paste0( one$chrom, ":", round(one$end/1e6,1), "Mb"),
                              tad_type = tad1$l$type,
                              n_genes  = length(genes_origTad),
                              n_expr   = nrow(ASE.embryo[gene_id %in% genes_origTad,]),
                              n_ase    = nrow(ASE.embryo[gene_id %in% genes_origTad & padj < 0.05,]),
                              n_up     = nrow(ASE.embryo[gene_id %in% genes_origTad & padj < 0.05 & log2FoldChange > 0,]),
                              n_genes_hit = length(genes_origTad[! genes_origTad %in% gene_ids_hit]),
                              n_expr_hit  = nrow(ASE.embryo[gene_id %in% genes_origTad[! genes_origTad %in% gene_ids_hit],]),
                              n_ase_hit   = nrow(ASE.embryo[gene_id %in% genes_origTad[! genes_origTad %in% gene_ids_hit] & padj<0.05,]))  )
    
    if(width(tad2$l)<=1) {
      message ("    ...skip (TAD 2)")
      next;
    }
    new_tad <- c( if (one$strand == "+") tad1$l else tad1$r,
                  if (two$strand == "+") tad2$r else tad2$l  )
    message("    New TAD is ", round(sum(width(new_tad))/1000), "kb (", new_tad$type[1], " + ", new_tad$type[2],")")
    
    genes_newTad  <- subsetByOverlaps(GENE_SPAN, new_tad)$gene_id
    COUNTS = rbind(COUNTS, 
                   data.table(genome   = "balancer",
                              tad_size = sum(width(new_tad)),
                              tad_id   = paste0(one$chrom, ":", round(one$end/1e6,1), "Mb(", one$strand, ") & ", two$chrom, ":", round(two$start/1e6,1), "Mb(", two$strand, ")"),
                              tad_type = paste(sort(new_tad$type),collapse=","),
                              n_genes  = length(genes_newTad),
                              n_expr   = nrow(ASE.embryo[gene_id %in% genes_newTad,]),
                              n_ase    = nrow(ASE.embryo[gene_id %in% genes_newTad & padj < 0.05,]),
                              n_up     = nrow(ASE.embryo[gene_id %in% genes_newTad & padj < 0.05 & log2FoldChange > 0,]),
                              n_genes_hit = length(genes_newTad[! genes_newTad %in% gene_ids_hit]),
                              n_expr_hit  = nrow(ASE.embryo[gene_id %in% genes_newTad[! genes_newTad %in% gene_ids_hit],]),
                              n_ase_hit   = nrow(ASE.embryo[gene_id %in% genes_newTad[! genes_newTad %in% gene_ids_hit] & padj<0.05,]))   ) 
}


### Without removing genes hit by breakpoints:
ggplot(COUNTS[genome=="balancer",]) + aes(tad_type, n_ase / n_expr) + 
  geom_boxplot() + geom_point() + 
  geom_text(aes(label = paste0(n_ase,"/",n_expr)), hjust=0, vjust=0) +
  ggtitle("Without removing genes hit by bp.")
ggsave(file="plot1.pdf")

### After removing genes hit by breakpoints:
ggplot(COUNTS[genome=="balancer",]) + aes(tad_type, n_ase_hit/n_expr_hit) + 
  geom_boxplot() + geom_point() + 
  geom_text(aes(label = paste0(n_ase_hit,"/",n_expr_hit)), hjust=0, vjust=0) +
  ggtitle("After removing genes hit by bp.")
ggsave(file="plot2.pdf")


#ggplot(COUNTS[genome=="dm6"]) + aes(tad_size, n_genes, col = tad_type) + geom_point()
#ggplot(COUNTS[genome=="dm6",]) + aes(tad_type, ase_frac) + geom_boxplot() + geom_jitter()
COUNTS[genome == "balancer",]




message("Calculating average number of genes/ase genes per TAD...")
### Average number of genes/ase genes in active/inactive TADs:
tadgr$n_genes = -1;
tadgr$n_expr  = -1;
tadgr$n_ase   = -1;
for (i in 1:length(tadgr)) {
  gene_ids_in_TAD <- subsetByOverlaps(GENE_SPAN, tadgr[i])$gene_id
  tadgr$n_genes[i] = length(gene_ids_in_TAD)
  tadgr$n_expr[i]  = nrow(ASE.embryo[gene_id %in% gene_ids_in_TAD,])
  tadgr$n_ase[i]   = nrow(ASE.embryo[gene_id %in% gene_ids_in_TAD & padj < 0.05,])
}
tadgr$ase_frac = tadgr$n_ase / tadgr$n_expr

# Violin plot: ASE fraction in active/inactive TADs
ggplot(as.data.table(tadgr)) + aes(type, n_ase/n_expr) + geom_violin() +
  geom_boxplot(data = COUNTS[genome=="balancer",], aes(tad_type, n_ase/n_expr)) +
  geom_jitter(data = COUNTS[genome=="balancer",], aes(tad_type, n_ase/n_expr))
ggsave(file="plot3.pdf")

# scatter plot: ASE fraction vs. TAD size
ggplot(as.data.table(tadgr)) + aes(end-start, n_ase/n_expr, col = type) + 
  geom_point() + scale_x_log10(breaks=c(1e4,1e5,1e6),labels=c("10kb","100kb","1Mb")) +
  annotation_logticks(side="b")
ggsave(file="plot4.pdf")

# difference in ASE frac between active and inactive TADs?
message("Test whether ASE fraction between active and inactive TADs is different:")
print(wilcox.test(tadgr$ase_frac[tadgr$type == "Active"], 
            tadgr$ase_frac[tadgr$type == "Inactive"]))
