
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(grid))
#source("/Volumes/korbel/shared/projects/drosophila_balancer/analyses/common/R/themes.R")
source("/g/korbel/shared/projects/drosophila_balancer/analyses/common/R/themes.R")


# Min number of reads per kb of exon
minCount = 20
cat("Using", minCount, "as threshold for minimal #reads per kb exon\n")

# Get gene sizes to calculate RPKM
f_gff <- "/g/korbel/shared/projects/drosophila_balancer/analyses/SNV_annotation/reformatted.gff"
#f_gff <- "/Volumes/korbel/shared/projects/drosophila_balancer/analyses/SNV_annotation/reformatted.gff"
gene_size <- read.table(f_gff, sep="\t", header=T) %>%
  filter(feature %in% c("exon","CDS","five_prime_UTR","three_prime_UTR")) %>%
  group_by(gene_id) %>%
  dplyr::summarize(length = sum(end - start +1),
                   num_intvls = length(strand),
                   strand = first(strand))


cat("Reading balancer counts...\n")
bal = read.table(snakemake@input[["bal"]]) %>% filter(!grepl('^__',V1))
cat("Reading virginizer counts...\n")
vrg = read.table(snakemake@input[["vrg"]]) %>% filter(!grepl('^__',V1))

m <- merge(bal, vrg, by="V1") %>%
        mutate(balancer = V2.x, virginizer = V2.y) %>%
        select(gene_id = V1, balancer, virginizer)
m <- merge(m, gene_size, by="gene_id") %>%
        mutate(virginizer = virginizer/length*1000,
               balancer = balancer/length*1000) %>%
        filter(virginizer+balancer >= minCount)

tr1 = grobTree(textGrob(paste("min.",minCount,"reads/kb of exon"), x=0.98,  y=0.88, hjust=1, vjust=1))
tr2 = grobTree(textGrob(paste("#genes:",dim(m)[1]), x=0.98,  y=0.78, hjust=1, vjust=1))

ggplot(m) + aes(balancer+virginizer, balancer/(balancer+virginizer)) + geom_point(alpha=0.4) + scale_x_log10(breaks=c(1,10,100,1000,10000)) + xlab("Total number of (separable) reads per kb exon") + thm$commonTheme + thm$pdfTheme
ggplot(m) + aes(balancer/(balancer+virginizer)) + geom_density()

p = ggplot(m) + aes(balancer/(balancer+virginizer)) + 
      geom_density(fill="black", alpha=0.5) + 
      annotation_custom(tr1) + annotation_custom(tr2) +
      ggtitle("Gene-level allelic fraction") + xlab("balancer fraction") + thm$commonTheme + thm$pdfTheme

ggsave(snakemake@output[[1]], p, width=6, height=4)
