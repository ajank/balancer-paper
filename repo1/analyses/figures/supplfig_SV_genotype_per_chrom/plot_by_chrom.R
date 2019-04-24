library(ggplot2)
library(data.table)
library(scales)
source("scripts/colors.R")


# Read SNVs
snvs = fread(snakemake@input[["snvs"]])
colnames(snvs) =  c("chrom", "Genotype", "count")

# Fix Genotype column
snvs$Genotype = sub('__','_', snvs$Genotype)
# Filter
snvs <- snvs[grepl('^chr[23X]$', chrom) & Genotype %in% c("0/1_0/0", "0/1_1/1"),]
# Rename Genotypes
snvs <- snvs[, Genotype := factor(Genotype, levels = c("0/1_0/0", "0/1_1/1", "1/1_1/1"), labels = c("balancer", "wild type", "common"))][]
# total counts
snvs <- merge(snvs, snvs[, .(total = sum(count)), by = .(chrom)], by = "chrom")
snvs[["Variant type"]] = "SNVs"

# Read SVs
dels = fread(snakemake@input[["dels"]])
dups = fread(snakemake@input[["dups"]])
colnames(dels) = colnames(dups) = c("chrom","start","end","Genotype")
d = rbind(cbind(dels, `Variant type` = "DEL"),
          cbind(dups, `Variant type` = "DUP"))

# Fix a few variables:
d <- d[, Genotype := substr(Genotype, 1, 7)][]
d <- d[, chrom := substr(chrom, 1, 4)][]
d <- d[grepl('^chr[23X]$', chrom),]

# fix Genotypes in DUPs
d <- d["0/1_0/1" == Genotype, Genotype := "0/1_1/1"][]

# Rename Genotypes
d <- d[, Genotype := factor(Genotype, levels = c("0/1_0/0", "0/1_1/1", "1/1_1/1"), labels = c("balancer", "wild type", "common"))][]

# Size classes:
d$addon = ""
d[`Variant type` == "DEL" & end-start < 50,]$addon = " (<50bp)"
d[`Variant type` == "DEL" & end-start >= 50 & end-start < 160,]$addon = " (50-159bp)"
d[`Variant type` == "DEL" & end-start >= 160,]$addon = " (>=160bp)"
d <- d[, `Variant type`:= paste0(`Variant type`, addon)]
d <- d[, `Variant type`:= factor(`Variant type`, levels = c("DUP", "DEL (>=160bp)", "DEL (50-159bp)", "DEL (<50bp)"), ordered = T)][]

# Summary table
d <- d[Genotype != "common", ]
e = merge(d[, .(count = nrow(.SD)), by = .(`Variant type`, Genotype, chrom)], 
          d[, .(total = nrow(.SD)), by = .(`Variant type`, chrom)], 
          by = c("Variant type", "chrom"))

e = rbind(e, snvs)



plt <- ggplot(e) + 
    aes(x = `Variant type`, y = count/total, fill=`Genotype`) + 
    geom_bar(stat="identity") +
    theme_minimal() +
    theme(legend.position = "bottom") +
    scale_fill_manual(values = c(`wild type` = as.character(my_colors["wildtype"]),
                                 `balancer`  = as.character(my_colors["balancer"]),
                                 `common`    = as.character(my_colors["common"]))) +
    coord_flip() +
    xlab(NULL) +
    scale_y_continuous(labels = percent) +
    ylab("Fraction (common SV excluded)") +
    ggtitle("Variant genotypes") +
    geom_text(data = e[, .(count = sum(count)), by = .(`Variant type`, chrom)], 
              aes(x = `Variant type`, label = paste("n =",count)), 
              y = 0.97, hjust = 1, col = "white", size = 3, inherit.aes = F) +
  facet_grid(chrom ~ .)

ggsave(plt, filename = snakemake@output[[1]], width = 4, height = 5)
