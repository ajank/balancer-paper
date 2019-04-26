library(ggplot2)
library(data.table)
#source("scripts/colors.R")

# Colors from http://www.colourlovers.com/palette/1473/Ocean_Five

input = snakemake@input[[1]]
d = fread(input)
colnames(d) = c("chrom","start","end","GT")
d[, GT := substr(GT,1,7)][, GT := sub('0/1_0/1', '0/1_1/1', GT)]


col = "grey"
if (grepl('DEL', input)) col = "firebrick1"
if (grepl('DUP', input)) col = "dodgerblue"

plt <- ggplot(d) +
    aes(x=end-start, y = ..count.. + 1, fill = GT) +
    geom_histogram(binwidth=0.1, alpha=1) +
    theme_minimal() +
    scale_x_log10(breaks = c(10,100,1000,10000,100000), labels=c("10bp","100bp","1kb","10kb","100kb"),
                  minor_breaks = c((2:9)*10, (2:9)*100, (2:9)*1000, (2:9)*1e4, (2:9)*1e5)) +
    #scale_y_log10(breaks = c(1,2,11,101,1001, 10001, 100001), labels=c("None",1,10,100,"1,000", "10000", "100,000"),
    #              minor_breaks = c(6,51,501)) +
    theme(legend.position = c(0.85,0.75), legend.background = element_rect(fill="white", size=0.2)) +
    ylab("Count") + xlab("SV size") +
    #coord_cartesian(xlim = c(10,5e5)) +
    scale_fill_manual(values = c(`0/1_0/0` = "#00A0B0",
                                 `0/1_1/1` = "#EDC951",
                                 `1/1_1/1` = "#6A4A3C"),
                      labels = c(`0/1_0/0` = "balancer",
                                 `0/1_1/1` = "wild type",
                                 `1/1_1/1` = "common"),
                      name = NULL)
ggsave(plt, filename = snakemake@output[[1]], width = 4.8, height = 3)
