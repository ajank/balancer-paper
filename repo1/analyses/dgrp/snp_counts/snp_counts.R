library(data.table)
library(scales)
library(ggplot2)
d = fread("snp_counts.txt")
pdf("snp_counts.pdf", width=4, height=3)
ggplot(d) +
    geom_boxplot(aes(x=1,y=V2), fill="dodgerblue") +
	geom_jitter(aes(x=1,y=V2), size=0.5, alpha=0.5) +
	scale_y_continuous(label=comma) +
	theme_classic() +
	xlab(NULL) +
	theme(axis.text.x=element_blank()) +
	ggtitle("SNP counts in DGRP lines") +
	ylab("count") 
dev.off()
