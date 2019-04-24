library(data.table)
library(ggplot2)
source("scripts/colors.R")


f_in = snakemake@input[[1]]
f_out = snakemake@output[[1]]

d = fread(f_in)
colnames(d) = c("type","number","genotype")

pvals = d[, .(pval = (frank(.SD, number)[type == "real"])/.N), by = genotype]

plt <- ggplot(d[type == "rand",]) + 
    aes(number, fill = genotype) + 
    geom_histogram(binwidth = 5) + 
    geom_vline(data = d[type=="real",], aes(xintercept = number), color="darkorange") +
    facet_wrap( ~ genotype, ncol=1) +
    geom_label(data = pvals, aes(label = paste0("hat(p)==",round(pval,3))), 
               x=Inf, y=Inf, hjust=1, vjust=1, inherit.aes=F, parse=T) +
    theme_minimal() +
    xlab("Number of exons touched by deletions") +
    scale_fill_manual(values = c(`wild type` = as.character(my_colors["wildtype"]),
                                 `balancer`  = as.character(my_colors["balancer"]),
                                 `common`    = as.character(my_colors["common"]))) +
    guides(fill=FALSE)
ggsave(plt, filename = f_out, width = 4, height = 3)
