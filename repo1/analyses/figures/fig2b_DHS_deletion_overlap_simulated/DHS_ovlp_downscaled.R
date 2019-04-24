library(data.table)
library(ggplot2)
source("scripts/colors.R")

f_in = snakemake@input[[1]]
f_out = snakemake@output[[1]]

d = fread(f_in)
d$type[ d$type == "virginizer" ] = "wild type"

# p_val = wilcox.test(d$count[d$type=="wild type"], d$count[d$type=="balancer"])$p.value

plt <- ggplot(d) + 
    aes(count, fill = type) + 
    geom_histogram(binwidth = 1, alpha = 0.5, position = "identity") + 
    theme_minimal() +
    xlab("Number of DHS peaks >=5% deleted") +
    scale_fill_manual(values = c(`wild type`  = as.character(my_colors["wildtype"]),
                                 `balancer`   = as.character(my_colors["balancer"]),
                                 `background` = "#bbbbbb")) +
    geom_vline(data = d[, .(mean = mean(count)), by = type], aes(xintercept = mean, col = type)) +
    scale_color_manual(values = c(`wild type` = as.character(my_colors["wildtype"]),
                                 `balancer`   = as.character(my_colors["balancer"]),
                                 `background` = "#bbbbbb")) +
    guides(col = FALSE) + 
    theme(legend.position = "bottom")
ggsave(plt, filename = f_out, width = 4, height = 3)
