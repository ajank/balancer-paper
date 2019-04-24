suppressWarnings(suppressPackageStartupMessages(library(ggplot2)))

args = commandArgs(trailingOnly = T)
d = read.table(file('stdin'))

segments = data.frame(x   = c(0.1, 0.5, 1),
					  y   = c(d$V2[d$V1==0.1], d$V2[d$V1==0.5], d$V2[d$V1==1]),
					  col = c("red","blue","black"))

outfile <- "ovl.pdf"
if (length(args)>0) outfile <- args[1] 

p<-ggplot(d) + aes(V1,V2) + geom_line() +  
    ggtitle(paste("10% ovl:", d$V2[d$V1==0.1])) + 
    theme_minimal() + xlab("ovl. required") + ylab("count") +
	geom_segment(data = segments, aes(x=0, xend=x, y=y, yend=y, col=col), linetype="dashed") +
	geom_segment(data = segments, aes(x=x, xend=x, y=y, yend=0, col=col), linetype="dashed") +
	geom_text(data = segments, aes(x=x, y=y, label=y), hjust=c(0,0,1), vjust=c(0,0,1)) +
	guides(colour=FALSE)
ggsave(filename = outfile, width=4, height=3)

