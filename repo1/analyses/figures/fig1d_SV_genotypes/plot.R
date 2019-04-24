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
snvs <- snvs[grepl('^chr[23]$', chrom) & Genotype %in% c("0/1_0/0", "0/1_1/1"),]
# Rename Genotypes
snvs <- snvs[, Genotype := factor(Genotype, levels = c("0/1_0/0", "0/1_1/1", "1/1_1/1"), labels = c("balancer", "wild type", "common"))][]
# Coalpse chromosomes
snvs <- snvs[, .(count = sum(count)), by = .(Genotype)]
snvs$total = snvs[, .(total = sum(count))]
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
d <- d[grepl('^chr[23]$', chrom),]

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
table(d$Genotype,d$`Variant type`)


d <- d[Genotype != "common", ]
e = merge(d[, .(count = nrow(.SD)), by = .(`Variant type`, Genotype)], 
          d[, .(total = nrow(.SD)), by = .(`Variant type`)], 
          by = c("Variant type"))

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
    ggtitle("Variants on chromosomes 2 and 3") +
    geom_text(data = e[, .(count = sum(count)), by = `Variant type`], 
              aes(x = `Variant type`, label = paste("n =",count)), 
              y = 0.97, hjust = 1, col = "white", size = 3, inherit.aes = F)
#+
#    scale_x_log10(breaks = c(10,100,1000,10000,100000), labels=c("10bp","100bp","1kb","10kb","100kb")) +
#    scale_y_log10(breaks = c(1,2,11,101,1001), labels=c("None",1,10,100,1000),
#                  minor_breaks = c(6,51,501)) +
#    scale_fill_manual(values=c(DEL="dodgerblue1", DUP="red1")) +
#    ggtitle("CNV size distribution") +
#    theme(legend.position = c(0.85,0.75), legend.background = element_rect(fill="white", size=0.2)) +
#    ylab("Count") + xlab("SV size")
ggsave(plt, filename = snakemake@output[[1]], width = 4, height = 3)
