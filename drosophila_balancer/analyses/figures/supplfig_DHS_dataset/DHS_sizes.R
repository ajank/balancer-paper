library(ggplot2)
library(data.table)

format_Mb <- function(x) {
    return( paste(round(x/1e3,0), "kb")  )
}

peaks = fread(snakemake@input[["peaks"]])
colnames(peaks) = c("chrom","start","end")

plt <- ggplot(peaks) + 
    aes(x = end-start) + 
    geom_histogram(binwidth=100, alpha=1) +
    theme_minimal() +
    scale_x_continuous(labels = format_Mb) +
    ggtitle("DHS Peaks from Thomas et al.") +
    ylab("Count") + xlab("Peak size")
ggsave(plt, filename = snakemake@output[[1]], width = 4, height = 3)

message("DHS number of peaks: ", nrow(peaks))
message("DHS median size:     ", median(peaks[, end - start]))
