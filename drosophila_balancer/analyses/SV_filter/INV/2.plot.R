library(ggplot2)

f = file("stdin")
e = read.table(f)
#e = read.table("test")
colnames(e) = c("chrom", "start", "end", "id", "precise", "cross", "vrg")
e$GT = interaction(e$cross, e$vrg)


ggplot(e) + aes(end - start) + 
    geom_histogram(binwidth=.1) + 
    scale_x_log10(breaks = c(1e2, 1e3, 1e4, 1e5, 1e6, 1e7), labels= c("100b", "1kb", "10kb", "100kb", "1Mb", "10Mb")) + 
    facet_grid(GT ~ .) +
    theme_minimal()
ggsave(filename = "delly0.7.5.INV.filter1.pdf", width=6, height=9)
