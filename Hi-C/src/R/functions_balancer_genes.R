library(assertthat, quietly = T)
library(DESeq2, quietly = T)
library(dplyr, quietly = T)
library(ggplot2, quietly = T)

if (!exists("data.loaded"))
{
  Sascha.dir <- "/g/korbel/shared/projects/drosophila_balancer"
  Alek.dir <- "/g/furlong/project/33_Hi-C"
  source(paste0(Alek.dir, "/src/R/plot_map_binned_functions.R"))
  env.libs <- environment()

  env.dm6 <- new.env(parent = env.libs)
  env.dm6$genome <- "dm6"
  source("Hi-C-ggbio/files.dm6.R", local=env.dm6)
  source("Hi-C-ggbio/genes.R", local=env.dm6)
  data.loaded <- T
}

# breakpoint coordinates (possibly multiple coordinates for the same breakpoint id)
bp <- fread("analysis/breakpoints_dm6.tab", header = T)
# intervals between breakpoints
bp_unique <- unique(bp[, c("id", "chrom", "breakpoint"), with = F])
bp_int_dt <- bp_unique[, list(start = head(breakpoint, -1), end = tail(breakpoint, -1)), by = "chrom"]
bp_int_gr <- GRanges(bp_int_dt$chrom, IRanges(bp_int_dt$start + 1L, bp_int_dt$end))

# now we can limit to non-terminal breakpoints
bp <- bp[!grepl("start$|end$", bp$id), ]
# bp_gr <- unique(GRanges(bp$chrom, IRanges(bp$breakpoint, bp$breakpoint - 1L)))
bp_gr <- GRanges(bp[, list(start = min(breakpoint) + 1L, end = max(breakpoint)), by = c("chrom", "id", "duplication_size")])

# mean coordinate for each breakpoint
bp_mean <- bp[, list(breakpoint = mean(breakpoint)), by = c("chrom", "id")]
bp_mean_gr <- GRanges(bp_mean$chrom, IRanges(bp_mean$breakpoint + 1L, bp_mean$breakpoint))

# chromosome sizes
dm6_chrom_sizes <- fread("/g/furlong/project/39_Balancer_Hi-C/liftOver/dm6.chrom.sizes", header = F, col.names = c("chrom", "chrom_size"))
dm6bal3_chrom_sizes <- fread("/g/furlong/project/39_Balancer_Hi-C/liftOver/dm6bal3.chrom.sizes", header = F, col.names = c("chrom", "chrom_size"))


extract_all_genes <- function(genes, expr, prefix)
{
  expr.aggr <- expr %>% group_by(gene_id) %>% summarize(FPKM = sum(FPKM), meanCount = mean(count))
  genes.aggr <- merge(genes, expr.aggr, by = "gene_id", all.x = T)
  message(prefix, ": ", nrow(genes.aggr), " all genes")

  return(genes.aggr)
}

annotate_genes <- function(genes)
{
  genes <- genes[grepl("^chr[23]", chrom), ]
  genes[, MEI := F]
  genes[, trivial_DE := F]
  genes[, trivial_class := ""]

  # genes_bp_ov <- findOverlaps(genes_gr, bp_gr)
  # genes_bp_ov <- findOverlaps(TSS_gr, evalq(bp_gr, env.dm6))
  # genes$breakpoint_class <- factor("away from breakpoint", c("away from breakpoint", "close to breakpoint"))
  # genes$breakpoint_class[unique(queryHits(genes_bp_ov))] <- "close to breakpoint"

  # genes$gene_id_by_fold_change <- with(genes, reorder(gene_id, log2FoldChange))
  # genes$gene_id_by_p_value <- with(genes, reorder(gene_id, padj))
  # genes$gene_id_by_FPKM <- with(genes, reorder(gene_id, FPKM))

  # exclude genes that are split, deleted or duplicated by breakpoint
  genes_bp_ov <- findOverlaps(GRanges(genes), bp_gr, minoverlap = 0L)
  ov <- cbind(genes[queryHits(genes_bp_ov)], as.data.table(bp_gr[subjectHits(genes_bp_ov)]))
  write.table(ov, "analysis/balancer/genes_split_by_breakpoint.tab",
    sep = "\t", quote = F, row.names = F)
  sel <- seq_len(nrow(genes)) %in% queryHits(genes_bp_ov)
  genes$trivial_DE[sel] <- T
  genes$trivial_class[sel] <- "breakpoint"

  anno <- fread("analysis/balancer/genes_allele_specific_CNV.tab", header = T)
  m <- match(anno$gene_id, genes$gene_id)
  stopifnot(!is.na(m))
  genes$trivial_DE[m] <- T
  genes$trivial_class[m] <- paste(genes$trivial_class[m], "CNV")

  anno <- fread("analysis/balancer/genes_allele_specific_MEI.tab", header = T)
  m <- match(anno$gene_id, genes$gene_id)
  stopifnot(!is.na(m))
  genes$MEI[m] <- T
  # genes$trivial_DE[m] <- T
  # genes$trivial_class[m] <- paste(genes$trivial_class[m], "MEI")

  genes[, trivial_class := ifelse(trivial_class == "", "none", sub("^ *", "" , trivial_class))]

  return(genes)
}

genes.all <- evalq(extract_all_genes(genes.embryo, expr.embryo, "genes.all"), env.dm6)
genes <- annotate_genes(genes.all)

levels(genes$signf)[levels(genes$signf) == "significant"] <- "s"
levels(genes$signf)[levels(genes$signf) == "insignificant"] <- "i"
levels(genes$signf)[levels(genes$signf) == "not tested"] <- "n"

genes[, signf := factor(signf, c("s", "i", "n"))]
genes[, signf2 := factor(ifelse(signf == "n", "n", ifelse(signf == "i", "i", ifelse(log2FoldChange > 0, "b", "w"))), c("b", "w", "i", "n"))]
genes[, pos := (start + end) / 2]
genes[, tss := ifelse(strand == "+", start, end)]
genes[, sign := ifelse(padj<0.05, ifelse(log2FoldChange > 0, "bal", "vrg"), "not")]
genes[, trivial_DE_class := factor(ifelse(signf == "s", trivial_class, "ZZZnon-DE"), c("breakpoint", "breakpoint MEI", "CNV", "CNV MEI", "MEI", "none", "ZZZnon-DE"))]

# Embryo data (4 replicates)

n=c(
  "SampleName",              "File",                                                                "Haplotype",   "Replicate"
)
x=c(
  "N1_6.8_father.vrg_1",     "gene_expression/counts/N1_pool_6-8h_rep1.htseq-count-rev.txt",        "unseparated", "1",
  "N1_6.8_father.vrg_1_BAL", "readSeparation/counts/N1_pool_6-8h_rep1.htseq-count-rev.alt.txt",     "balancer",    "1",
  "N1_6.8_father.vrg_1_VRG", "readSeparation/counts/N1_pool_6-8h_rep1.htseq-count-rev.ref.txt",     "virginizer",  "1",
  "N1_6.8_father.vrg_2",     "gene_expression/counts/N1_pool_6-8h_rep2.htseq-count-rev.txt",        "unseparated", "2",
  "N1_6.8_father.vrg_2_BAL", "readSeparation/counts/N1_pool_6-8h_rep2.htseq-count-rev.alt.txt",     "balancer",    "2",
  "N1_6.8_father.vrg_2_VRG", "readSeparation/counts/N1_pool_6-8h_rep2.htseq-count-rev.ref.txt",     "virginizer",  "2",
  "N1_6.8_mother.vrg_1",     "gene_expression/counts/N1sex_pool_6-8h_rep1.htseq-count-rev.txt",     "unseparated", "3",
  "N1_6.8_mother.vrg_1_BAL", "readSeparation/counts/N1sex_pool_6-8h_rep1.htseq-count-rev.alt.txt",  "balancer",    "3",
  "N1_6.8_mother.vrg_1_VRG", "readSeparation/counts/N1sex_pool_6-8h_rep1.htseq-count-rev.ref.txt",  "virginizer",  "3",
  "N1_6.8_mother.vrg_2",     "gene_expression/counts/N1sex_pool_6-8h_rep2.htseq-count-rev.txt",     "unseparated", "4",
  "N1_6.8_mother.vrg_2_BAL", "readSeparation/counts/N1sex_pool_6-8h_rep2.htseq-count-rev.alt.txt",  "balancer",    "4",
  "N1_6.8_mother.vrg_2_VRG", "readSeparation/counts/N1sex_pool_6-8h_rep2.htseq-count-rev.ref.txt",  "virginizer",  "4"
)
D_e = data.frame(matrix(x, byrow=T,ncol=4))
colnames(D_e) = n
DE_e <- DESeqDataSetFromHTSeqCount(D_e, design = ~ Replicate + Haplotype, directory = "/g/korbel/shared/projects/drosophila_balancer/analyses")

dt <- as.data.table(counts(DE_e))
gene_counts <- data.table(
  gene_id     = rownames(counts(DE_e)),
  unseparated = with(dt, `N1_6.8_father.vrg_1` + `N1_6.8_father.vrg_2` + `N1_6.8_mother.vrg_1` + `N1_6.8_mother.vrg_2`),
  bal         = with(dt, `N1_6.8_father.vrg_1_BAL` + `N1_6.8_father.vrg_2_BAL` + `N1_6.8_mother.vrg_1_BAL` + `N1_6.8_mother.vrg_2_BAL`),
  wt          = with(dt, `N1_6.8_father.vrg_1_VRG` + `N1_6.8_father.vrg_2_VRG` + `N1_6.8_mother.vrg_1_VRG` + `N1_6.8_mother.vrg_2_VRG`)
)
gene_counts[, allelespec := bal + wt]

# impute the counts to "genes"
m <- match(genes$gene_id, gene_counts$gene_id)
genes[, count.unseparated := ifelse(is.na(m), 0L, gene_counts$unseparated[m])]
genes[, count.allelespec := ifelse(is.na(m), 0L, gene_counts$allelespec[m])]
genes[, signfs := factor(ifelse(signf == "n", ifelse(count.allelespec > 0, "ns", "nu"), as.character(signf)), c("s", "i", "ns", "nu"))]

#
#  read and annotate TSSes accordingly
#

TSS_gr <- rtracklayer::import.gff3("/g/furlong/genome/D.melanogaster/Dm6/6.05/gff/dmel-all-r6.05.UCSC_names.genes.gff.gz")
# add missing strand information for one gene
strand(TSS_gr)[TSS_gr$gene_id == "FBgn0002781"] <- "-"
stopifnot(as.vector(strand(TSS_gr) %in% c("+", "-")))

# remove unused columns
elementMetadata(TSS_gr) <- elementMetadata(TSS_gr)[, c("source", "type", "gene_id", "Parent", "ID", "Name")]

# keep only transcripts
TSS_gr <- TSS_gr[grep("^FBtr", TSS_gr$ID)]
TSS_gr$type <- droplevels(TSS_gr$type)
stopifnot(!is.na(as.character(TSS_gr$Parent)))
# print(summary(width(TSS_gr)))

# convert CharacterList to a character vector
stopifnot(vapply(TSS_gr$Parent, length, FUN.VALUE = integer(1)) == 1)
TSS_gr$Parent <- unlist(TSS_gr$Parent)
stopifnot(length(setdiff(genes$gene_id, TSS_gr$Parent)) == 0)

# take only chr[23]
TSS_gr <- TSS_gr[grepl("^chr[23]", seqnames(TSS_gr))]

# take only TSS coordinate
TSS_gr <- resize(TSS_gr, 1L, fix = "start")

# copy the expression status from "genes"
for (col in c("log2FoldChange", "baseMean", "padj", "signf", "FPKM", "meanCount", "MEI", "trivial_DE", "trivial_class", "signf2", "pos", "trivial_DE_class", "sign"))
  elementMetadata(TSS_gr)[, col] <- genes[match(TSS_gr$gene_id, genes$gene_id), col, with = F]
TSS_gr$tss <- start(TSS_gr)

# factorize gene_id
TSS_gr$gene_id <- factor(TSS_gr$gene_id)
genes$gene_id <- factor(genes$gene_id, levels(TSS_gr$gene_id))

promoter_gr <- resize(TSS_gr, 0L, fix = "start") + 1e3
message("Considering ", length(promoter_gr), " transcripts, ", length(unique(promoter_gr)), " unique TSSes and ", length(unique(as.character(promoter_gr$Parent))), " genes")

# some useful data.tables
TSS <- as.data.table(TSS_gr)
setnames(TSS, "seqnames", "chrom")
TSS_tested <- TSS[signf != "n", ]

promoter <- as.data.table(promoter_gr)
setnames(promoter, "seqnames", "chrom")
promoter_tested <- promoter[signf != "n", ]

genes_tested <- genes[signf != "n", ]

#
#  colors, labels etc.
#

gene.colors = c(
  w = "#4daf4a",
  b = "#377eb8",
  s = "#ff7f00",
  i = "#999999",
  n = "#000000"
)
ggbio.gene.colors <- c(
  s = "darkorange",
  i = "dodgerblue4",
  n = "#999999",
  ns = "#6f94b9",
  nu = "#c0c0c0"
)
gene.alpha = c(
  w = 1,
  b = 1,
  s = 1,
  i = 0.6,
  n = 0.6
)
gene.labels = c(b = "up in balancer", w = "down in balancer", s = "DE genes", i = "non\u00adDE genes", n = "not tested",
  ns = "lowly expressed separable genes", nu = "unseparable genes")

format_Mb   <- function(x) {paste0(scales::comma(x/1e6))}

minuslog10 <- function(l) {
  l = -l
  l <- paste0("10^",l)
  # turn the 'e+' into plotmath format
  #l <- gsub("e", "%*%10^", l)
  # return this as an expression
  parse(text=l)
}

paste_names <- function(v, suffix)
{
  names(v) <- paste0(names(v), suffix)
  return(v)
}

DE_labels <- c(breakpoint = "balancer rearrangement", `breakpoint MEI` = "balancer rearrangement and MEI", CNV = "CNV", `CNV MEI` = "CNV and MEI", MEI = "MEI", none = "none of the above", `ZZZnon-DE` = "non\u00adDE genes")
DE_name <- "DE genes affected by"

DE_labels_counted <- DE_labels
for (l in levels(genes_tested$trivial_DE_class))
{
  sel <- names(DE_labels) == l
  DE_labels_counted[sel] <- paste0(DE_labels_counted[sel], " (", scales::comma(sum(genes_tested$trivial_DE_class == l)), ")")
}

DE_alpha_values <- c(breakpoint = 1, `breakpoint MEI` = 1, CNV = 1, `CNV MEI` = 1, MEI = 1, none = 1, `ZZZnon-DE` = 0.5)
DE_alpha <- scale_alpha_manual(
  values = c(DE_alpha_values, paste_names(DE_alpha_values, "_up"), paste_names(DE_alpha_values, "_down")),
  breaks = names(DE_labels), labels = DE_labels, name = DE_name
)
DE_alpha_counted <- scale_alpha_manual(
  values = c(DE_alpha_values, paste_names(DE_alpha_values, "_up"), paste_names(DE_alpha_values, "_down")),
  breaks = names(DE_labels), labels = DE_labels_counted, name = DE_name
)

DE_shape_values <- c(breakpoint = 15, `breakpoint MEI` = 15, CNV = 15, `CNV MEI` = 15, MEI = 15, none = 15, `ZZZnon-DE` = 16)
DE_shape <- scale_shape_manual(
  values = c(DE_shape_values, paste_names(DE_shape_values, "_up") * 0 + 2, paste_names(DE_shape_values, "_down") * 0 + 6),
  breaks = names(DE_labels), labels = DE_labels, name = DE_name
)

# DE_color <- scale_color_manual(values = c(not = "grey",  bal = as.character(gene.colors["b"]), vrg = as.character(gene.colors["w"])))
# DE_color_values <- c(breakpoint = "#0F8554", `breakpoint MEI` = "#38A6A5", CNV = "#CC503E", `CNV MEI` = "#5F4690", MEI = "#1D6996",
#     none = as.character(gene.colors["s"]), `ZZZnon-DE` = as.character(gene.colors["i"]))
DE_color_values <- c(breakpoint = "#5F4690", `breakpoint MEI` = "#CC503E", CNV = "#38A6A5", `CNV MEI` = "#0F8554", MEI = "#1D6996",
    none = as.character(gene.colors["s"]), `ZZZnon-DE` = as.character(gene.colors["i"]))
DE_color <- scale_color_manual(
  values = c(DE_color_values, paste_names(DE_color_values, "_up"), paste_names(DE_color_values, "_down")),
  breaks = names(DE_labels), labels = DE_labels, name = DE_name
)
DE_color_counted <- scale_color_manual(
  values = c(DE_color_values, paste_names(DE_color_values, "_up"), paste_names(DE_color_values, "_down")),
  breaks = names(DE_labels), labels = DE_labels_counted, name = DE_name
)
DE_fill <- scale_fill_manual(
  values = c(DE_color_values, paste_names(DE_color_values, "_up"), paste_names(DE_color_values, "_down")),
  breaks = names(DE_labels), labels = DE_labels, name = DE_name
)
# > library(rcartocolor); carto_pal(9, "Prism")
# [1] "#5F4690" "#1D6996" "#38A6A5" "#0F8554" "#73AF48" "#EDAD08" "#E17C05" "#CC503E" "#666666"
