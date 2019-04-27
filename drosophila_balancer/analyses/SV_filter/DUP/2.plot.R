library(ggplot2)
library(data.table)

args = commandArgs(trailingOnly=T)
if (length(args)==1 && grepl('\\.pdf', args[1])) { 
	out = args[1];
} else { 
	out = "out.pdf";
}
   

f = file("stdin")
e = read.table(f)

ggplot(e) + aes(V3, V4/V5, col=V2) + 
    geom_point(alpha=0.4) + 
    scale_x_log10(breaks = c(1e2, 1e3, 1e4, 1e5, 1e6, 1e7), labels= c("100b", "1kb", "10kb", "100kb", "1Mb", "10Mb")) + 
    scale_y_log10(breaks = c(0.25, 0.5, 1, 2, 4)) + 
    geom_hline(yintercept = c(0.8, 1.25), linetype="dotted") + theme_bw()

ggsave(filename = out, width=9, height=6)

# h = fread("/Volumes/korbel/shared/projects/drosophila_balancer/data/variants/SVs3/dup.snp.ints.bed")
# ggplot(h[,.(V11,V6)]) + geom_histogram(aes(V6), binwidth=0.05) + facet_wrap(~V11, scale="free")
