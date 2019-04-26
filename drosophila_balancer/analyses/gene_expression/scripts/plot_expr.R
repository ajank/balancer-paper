library(tidyr) # separate
library(dplyr)
library(ggplot2)
library(assertthat)
library(scales)
library(data.table)

myTheme = theme_minimal() +
  theme(text = element_text(size=12),
        title = element_text(size=12),
        axis.text = element_text(size=12),
        legend.text = element_text(size=14),
        strip.text = element_text(size=14))


# Get size of genes
f_gff     <- snakemake@config[["gff"]]
gene_size <- fread(f_gff) %>%
  filter(feature %in% c("exon","CDS","five_prime_UTR","three_prime_UTR")) %>%
  group_by(gene_id, gene_id) %>%
  dplyr::summarize(length = sum(end - start +1),
                   num_intvls = length(strand),
                   strand = first(strand)) %>%
  as.data.table


# Get HTSeq counts
sn <- snakemake@params[["sample_name"]]
rd <- snakemake@params[["read_len"]]

expr_lvl <- fread(snakemake@input[["count"]])
expr_lvl <- expr_lvl[!grepl("^__", V1), .(gene_id = V1, count = V2),]
assert_that(noNA(expr_lvl))



# Remove all genes that are not listed in "gene_size" and 
# calculate FPKM
total_read_count = sum(expr_lvl$count)
expr_lvl <- merge(expr_lvl, gene_size, by="gene_id") %>%
    mutate(FPKM = count*1000000/total_read_count*1000/length,
           mean_cov = rd * count / length)
assert_that(noNA(expr_lvl))


message("Processing ", sn, "...")
plt <- ggplot(expr_lvl) + 
  aes(FPKM) + 
  geom_histogram(binwidth=0.1) + 
  scale_x_log10(breaks=c(0.1,1,10,100,1000), label = comma) + 
  myTheme +
  geom_vline(aes(xintercept=median(FPKM)), col="red", linetype="dashed") +
  ggtitle(paste0("Gene expression ", sn, " (n=", nrow(expr_lvl), ")"))
ggsave(plt, filename = snakemake@output[[1]], width=6, height=4)
