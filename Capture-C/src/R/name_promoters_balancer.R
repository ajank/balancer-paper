library(GenomicRanges)
library(rtracklayer)
library(data.table)

#
#  DpnII fragments selected to target gene promoters
#

genome <- "dm6"
load("analysis/balancer/promoters_balancer.Rdata") # promoters, lib
promoters <- promoters[order(promoters$dpnp_id)]

#
#  Convert old dpnp_id to new frag_id
#

load(paste0("analysis/digest_DpnII_", genome, ".Rdata"))
stopifnot(frag$frag_id == 1:nrow(frag))

promoters$frag_id <- match(paste(seqnames(promoters), start(promoters)), paste(frag$chrom, frag$start + 4))
stopifnot(!is.na(promoters$frag_id))
promoters$dpnp_id <- NULL

lib$frag_id <- match(paste(seqnames(lib), start(lib)), paste(frag$chrom, frag$start + 4))
stopifnot(!is.na(lib$frag_id))
lib$dpnp_id <- NULL

#
#  FlyBase promoters for annotation/renaming
#

transcripts <- import.gff3(switch(genome, dm3 = "data/gff/dmel-transcript-r5.57.gff", dm6 = "data/gff/dmel-transcript-r6.05.gff"))
seqlevels(transcripts)[seqlevels(transcripts) == "dmel_mitochondrion_genome"] <- "M"
seqlevels(transcripts) <- paste0("chr", seqlevels(transcripts))

stopifnot(sapply(transcripts$Parent, length) <= 1)
transcripts$Parent <- sapply(transcripts$Parent, "[", 1)

#
#  Add gene names, extracted from transcript names
#

transcripts$GeneName <- NA
for (regexp in c("-R[A-Z]$", "-R[A-Z][A-Z]$"))
{
  sel <- grepl(regexp, transcripts$Name)
  transcripts$GeneName[sel] <- sub(regexp, "", transcripts$Name[sel])
}
stopifnot(!is.na(transcripts$GeneName))

tss <- resize(transcripts, 0, fix = "start")
sel <- !is.na(tss$Parent)
gene_name_dt <- data.table(id = tss$Parent[sel], name = tss$GeneName[sel], strand = as.factor(strand(tss))[sel])

# fix incomplete annotation
gene_name_dt$strand[gene_name_dt$id == "FBgn0002781" & gene_name_dt$name == "mod(mdg4)"] <- "-"
# update gene names in accordance with updated FlyBase IDs
gene_name_dt$name[gene_name_dt$id == "FBgn0283501" & gene_name_dt$name == "CR45992"] <- "CR46115"
gene_name_dt$name[gene_name_dt$id == "FBgn0283451" & gene_name_dt$name == "CG11509"] <- "br"

gene_name_dt <- unique(gene_name_dt)
stopifnot(length(unique(gene_name_dt$id)) == nrow(gene_name_dt))
stopifnot(length(unique(gene_name_dt$name)) == nrow(gene_name_dt))

#
#  Annotate the viewpoints by TSSes within +/- 1 kb
#  (in case of viewpoints associated with TSSes farther away -- take the distance to that gene as a threshold)
#

# calculate the distance between viewpoints and TSSes that are already associated to them
promoters$distance <- sapply(seq_along(promoters),
  function (i) min(distance(promoters[i], tss[!is.na(tss$Parent) & tss$Parent == promoters$Parent[i]])))
lib_distance <- as.data.table(promoters)[, list(distance = max(distance)), by = frag_id]
stopifnot(lib$frag_id == lib_distance$frag_id)
lib$distance <- lib_distance$distance

# extend the viewpoints by 1 kb (or by the distance calculated above, if greater), and overlap them with TSS annotations
lib_ext <- resize(lib, width(lib) + 2L * pmax(1000L, lib$distance), fix = "center")
ov <- findOverlaps(lib_ext, tss, minoverlap = 0L)

# build new TSS annotation
lib_to_transcript <- lib[queryHits(ov)]
elementMetadata(lib_to_transcript) <- elementMetadata(lib_to_transcript)[, "frag_id", drop = F]
lib_to_transcript$transcript_id <- tss$ID[subjectHits(ov)]
lib_to_transcript$transcript_name <- tss$Name[subjectHits(ov)]
lib_to_transcript$transcript_strand <- as.factor(strand(tss))[subjectHits(ov)]
lib_to_transcript$TSS <- start(tss)[subjectHits(ov)]
lib_to_transcript$TSS_distance <- distance(lib_to_transcript, tss[subjectHits(ov)])
lib_to_transcript$gene_id <- tss$Parent[subjectHits(ov)]
lib_to_transcript$gene_id[lib_to_transcript$transcript_id == "FBtr0309715"] <- "FBgn0263566" # update FlyBase ID
m <- match(lib_to_transcript$gene_id, gene_name_dt$id)
lib_to_transcript$gene_name <- gene_name_dt$name[m]
lib_to_transcript$gene_strand <- gene_name_dt$strand[m]

stopifnot(lib_to_transcript$transcript_strand == lib_to_transcript$gene_strand)
lib_to_transcript$transcript_strand <- NULL
# assert that all the previously associated TSSes are fetched
as_old <- with(promoters, paste(frag_id, Parent))
as_new <- with(lib_to_transcript, paste(frag_id, gene_id))
stopifnot(length(setdiff(as_old, as_new)) == 0)

lib_to_transcript$is.control <- lib_to_transcript$frag_id %in% lib$frag_id[grepl("^control.", lib$Name.tss)]

#
#  Add expression data
#

library(dplyr, quietly = T)
library(assertthat, quietly = T)

Sascha.dir <- "/g/korbel/shared/projects/drosophila_balancer"
Alek.dir <- "/g/furlong/project/33_Hi-C"
source(paste0(Alek.dir, "/src/R/plot_map_binned_functions.R"))
env.libs <- environment()

env.dm6 <- new.env(parent = env.libs)
env.dm6$genome <- "dm6"
source("../33_Hi-C/Hi-C-ggbio/files.dm6.R", local=env.dm6)
source("../33_Hi-C/Hi-C-ggbio/genes.R", local=env.dm6)

# env.dm6bal3 <- new.env(parent = env.libs)
# env.dm6bal3$genome <- "dm6bal3"
# source("../33_Hi-C/Hi-C-ggbio/files.dm6bal3.R", local=env.dm6bal3)
# source("../33_Hi-C/Hi-C-ggbio/genes.R", local=env.dm6bal3)

expr.embryo.repl <- evalq(expr.embryo, env = env.dm6)
expr.embryo <- as.data.table(expr.embryo.repl)[, list(FPKM = mean(FPKM)), by = "gene_id"]
lib_to_transcript$embryo.FPKM <- expr.embryo$FPKM[match(lib_to_transcript$gene_id, expr.embryo$gene_id)]

genes.embryo <- evalq(genes.embryo, env = env.dm6)
m <- match(lib_to_transcript$gene_id, genes.embryo$gene_id)
lib_to_transcript$embryo.log2FoldChange <- genes.embryo$log2FoldChange[m]
lib_to_transcript$embryo.padj <- genes.embryo$padj[m]
lib_to_transcript$embryo.signf <- genes.embryo$signf[m]

#
#  Check a few genes that have no expression quantified (mostly overlapping same-strand)
#

sel <- is.na(lib_to_transcript$embryo.FPKM)
message("The following genes were not quantified:")
cat(unique(lib_to_transcript$gene_id[sel]), "\n")

lib_to_transcript$embryo.log2FoldChange[sel] <- 0
lib_to_transcript$embryo.padj[sel] <- 0.99
lib_to_transcript$embryo.signf[sel] <- "not tested"

stopifnot(with(lib_to_transcript, (embryo.signf == "significant") == (embryo.padj < 0.05)))

# save differential expression according to the more stringent criterion
lib_to_transcript$is.DE <- lib_to_transcript$embryo.padj < 0.05 & abs(lib_to_transcript$embryo.log2FoldChange) > log2(1.5)

#
#  Add information about differential MEI
#

mei_genes <- fread("/g/korbel/shared/projects/drosophila_balancer/analyses/MEI_detection/own2/final.genes.txt", header = T)
mei_genes$differential_MEI <- T
m <- merge(as.data.table(lib_to_transcript), mei_genes, by = "gene_id", all.x = T, sort = F)
m$differential_MEI[is.na(m$differential_MEI)] <- F
lib_to_transcript$differential_MEI <- m$differential_MEI

#
#  Calculate distance to translocation breakpoints
#

br_dt <- fread("../39_Balancer_Hi-C/analysis/breakpoints_dm6.tab")
# br_dt <- br_dt[!grepl("start", id) & !grepl("end", id)]
br <- GRanges(br_dt$chrom, IRanges(br_dt$breakpoint + 1L, br_dt$breakpoint))
dist <- distanceToNearest(lib_to_transcript, br)
lib_to_transcript$breakpoint_distance <- elementMetadata(dist)$distance

#
#  Take only the genes considered in the original gene assignment OR differentially expressed ones
#

lt <- with(lib_to_transcript, paste(frag_id, gene_id))
lp <- with(promoters, paste(frag_id, Parent))
stopifnot(all(lp %in% lt))

lib_to_transcript$is.assigned <- lt %in% lp
lib_to_transcript <- lib_to_transcript[lib_to_transcript$is.assigned | lib_to_transcript$is.DE]

#
#  Reduce to genes
#

dt <- as.data.table(lib_to_transcript)

dt <- dt[,
  list(
    TSS_distance = min(TSS_distance)
  ),
  by = c("seqnames", "start", "end", "frag_id", "gene_id", "gene_name", "gene_strand",
    "is.control", "embryo.FPKM", "embryo.log2FoldChange", "embryo.padj", "embryo.signf",
    "is.DE", "differential_MEI", "breakpoint_distance", "is.assigned")]
lib_to_gene <- GRanges(dt)
stopifnot(anyDuplicated(lib_to_gene$gene_id) == 0)

#
#  Reduce to DpnII fragments
#

dt <- as.data.table(lib_to_gene)

dt <- dt[,
  list(
    gene_id = paste(gene_id, collapse = ","),
    gene_name = paste(gene_name, collapse = ","),
    gene_strand = paste(gene_strand, collapse = ","),
    embryo.FPKM = paste(embryo.FPKM, collapse = ","),
    embryo.log2FoldChange = paste(embryo.log2FoldChange, collapse = ","),
    embryo.padj = paste(embryo.padj, collapse = ","),
    embryo.signf = paste(embryo.signf, collapse = ","),
    is.DE = paste(is.DE, collapse = ","),
    TSS_distance = paste(TSS_distance, collapse = ","),
    differential_MEI = paste(differential_MEI, collapse = ",")
  ),
  by = c("seqnames", "start", "end", "frag_id", "is.control", "breakpoint_distance")]

lib_gff <- GRanges(dt)
stopifnot(anyDuplicated(lib_gff$frag_id) == 0)

names(lib_gff) <- lib_gff$gene_name
lib_gff$gene_id <- CharacterList(strsplit(lib_gff$gene_id, ","))
lib_gff$gene_name <- CharacterList(strsplit(lib_gff$gene_name, ","))
lib_gff$gene_strand <- CharacterList(strsplit(lib_gff$gene_strand, ","))
lib_gff$embryo.FPKM <- CharacterList(strsplit(lib_gff$embryo.FPKM, ","))
lib_gff$embryo.log2FoldChange <- CharacterList(strsplit(lib_gff$embryo.log2FoldChange, ","))
lib_gff$embryo.padj <- CharacterList(strsplit(lib_gff$embryo.padj, ","))
lib_gff$embryo.signf <- CharacterList(strsplit(lib_gff$embryo.signf, ","))
lib_gff$is.DE <- CharacterList(strsplit(lib_gff$is.DE, ","))
lib_gff$TSS_distance <- CharacterList(strsplit(lib_gff$TSS_distance, ","))
lib_gff$differential_MEI <- CharacterList(strsplit(lib_gff$differential_MEI, ","))

#
#  Reduce to DpnII fragments and annotate them with DE gene (if only one gene is associated)
#

dt <- as.data.table(lib_to_gene)

dt <- dt[,
  list(
    DE_class = if (sum(is.DE) == 1) "one DE gene" else if (sum(is.DE) == 0) "no DE genes" else "more than one DE gene",
    gene_class = if (.N == 1) "one gene" else "more than one gene",
    gene_id = if (.N == 1) gene_id else NA_character_,
    gene_name = if (.N == 1) gene_name else NA_character_,
    gene_strand = if (.N == 1) gene_strand else factor("*", levels(gene_strand)),
    embryo.FPKM = if (.N == 1) embryo.FPKM else NA_real_,
    embryo.log2FoldChange = if (.N == 1) embryo.log2FoldChange else NA_real_,
    embryo.padj = if (.N == 1) embryo.padj else NA_real_,
    TSS_distance = if (.N == 1) TSS_distance else NA_integer_,
    differential_MEI = if (.N == 1) differential_MEI else NA
  ),
  by = c("seqnames", "start", "end", "frag_id", "is.control", "breakpoint_distance")]

dt[, DE_class := factor(DE_class, c("no DE genes", "one DE gene", "more than one DE gene"))]
dt[, gene_class := factor(gene_class, c("one gene", "more than one gene"))]
lib <- GRanges(dt)
stopifnot(anyDuplicated(lib$frag_id) == 0)

stopifnot(lib$frag_id == lib_gff$frag_id)
names(lib) <- names(lib_gff)

#
#  Save to Rdata, GFF and BED
#

save(lib_to_transcript, lib_to_gene, lib_gff, lib, file = "analysis/balancer_cap2/viewpoints.Rdata")

export.gff3(lib_gff, "analysis/balancer_cap2/viewpoints.gff")

# a hack to squeeze the fragment id into the bed file
score(lib_gff) <- lib_gff$frag_id
export.bed(lib_gff, "analysis/balancer_cap2/viewpoints.bed")
