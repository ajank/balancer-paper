library(ggplot2)
library(data.table)
source("scripts/colors.R")

dels = fread(snakemake@input[["dels"]])
dups = fread(snakemake@input[["dups"]])
colnames(dels) = colnames(dups) = c("chrom","start","end","GT")
d = rbind(cbind(dels, `SV type` = "DEL"),
          cbind(dups, `SV type` = "DUP"))
d$chrom = substr(d$chrom,1,4)

plt <- ggplot(d) + 
    aes(x=end-start, y = ..count.. + 1, fill=`SV type`) + 
    geom_histogram(binwidth=0.1, position=position_dodge(), alpha=1) +
    theme_minimal() +
    scale_x_log10(breaks = c(10,100,1000,10000,100000), labels=c("10bp","100bp","1kb","10kb","100kb"),
                  minor_breaks = c(50, 500, 5000, 50000, 500000)) +
    scale_y_log10(breaks = c(1,2,11,101,1001), labels=c("None",1,10,100,1000),
                  minor_breaks = c(6,51,501)) +
    scale_fill_manual(values=c(DEL = as.character(my_colors["del"]),
                               DUP = as.character(my_colors["dup"]))) +
    ggtitle("CNV size distribution") +
    theme(legend.position = c(0.85,0.75), legend.background = element_rect(fill="white", size=0.2)) +
    ylab("Count") + xlab("SV size") +
    facet_grid(chrom ~ .)
    #annotation_logticks(size = 0.2, short = unit(0.5,"mm"), mid = unit(1,"mm"), long = unit(2,"mm")) #+ 
    #coord_cartesian(ylim = c(1.55,3000), xlim = c(6, 290000), expand = F)
ggsave(plt, filename = snakemake@output[[1]], width = 4, height = 6)