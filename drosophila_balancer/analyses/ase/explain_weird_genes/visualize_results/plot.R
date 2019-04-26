d=read.table("annotation.tsv", header=T, sep="\t")
library(ggplot2)


p <- ggplot(d) + aes(x = reorder(gene_id, order(cis.lfc)), 
                y=cis.lfc, 
                fill=grepl("MEI",Sascha), 
                col=NULL) + 
  geom_bar(stat="identity") + 
  scale_fill_manual(name="MEI?", values=c("grey","darkred")) + 
  theme_minimal() + 
  xlab("LFC of significant genes")
ggsave("annotation.pdf", p, width=12,height=3)
