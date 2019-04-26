library(assertthat)
args = commandArgs(trailingOnly = T)
stopifnot(length(args)>0)
stopifnot(file.exists(args[1]))
if (length(args)<2 || !grepl('.pdf$', args[2])) {
    args[2] = "out.pdf"
}
message("Input: ", args[1])
message("Output: ", args[2])


library(data.table)
library(ggplot2)
d <- fread(args[1])
dim(d)
plt <- ggplot(d[test=="shuffle",]) + 
    aes(overlaps) + 
    geom_histogram(binwidth=20, fill='lightgrey') + 
    facet_grid(TADcalls ~ ArcProt, scale='free') + 
    theme_bw() + 
    geom_vline(data = d[test=="actual",], aes(xintercept = overlaps), color = "red")
w = 2 + 2*length(unique(d$ArcProt))
h = 1 + 2*length(unique(d$TADcalls))
ggsave(filename = args[2], width = w, height = h)

