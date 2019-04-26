library(data.table)
library(ggplot2)
library(scales)
library(GenomicRanges)
library(assertthat)

dm6_overview <- fread("/g/furlong/project/39_Balancer_Hi-C/analysis/assembly_dm6bal3.tab", header = T)
dm6_overview[, strand := ifelse(breakpoint1 < breakpoint2, "+", "-")]
setnames(dm6_overview, c("breakpoint1", "breakpoint2"), c("start", "end"))



hack_prepare_dm6bal3_overview <- function(dt, chromosome = 'chr2')
{
    x <- dt[grepl(chromosome, chrom),]
    x$s = 0
    x$e = 0
    x$c = ""
    temp = 0;
    chrom = paste0(chromosome, "L")
    for (i in 1:nrow(x)) 
    { 
        if (i>1 && x$col[i] == x$col[i-1]) {
            chrom = paste0(chromosome, "R")
            temp = 0
        }
        x$s[i] = temp +1;
        x$e[i] = temp + abs(x$end[i] - x$start[i]); 
        temp = temp + abs(x$end[i] - x$start[i]);
        x$c[i] = chrom
    }
    x$chrom = x$c
    x$start = x$s
    x$end   = x$e
    x[, c("s","e","c"):=NULL]
    return(x)
}


plot_overview <- function(wh, genome = "dm6")
{
    assert_that(is.character(genome) && length(genome)==1)
    assert_that(genome == "dm6" || genome == "dm6bal3")
    ch <- as.character(seqnames(wh)[1])
    d <- dm6_overview[grepl(substr(ch,1,4), chrom),]
    if (genome=="dm6") {
        d[, y:= 0]
    } else {
        d <- hack_prepare_dm6bal3_overview(d, substr(ch,1,4))
        d[, y:= ifelse(strand=="+", 0.5, -0.5)]
    }
    plt <- ggplot(d) + 
        aes(xmin=start, xmax=end, ymin=y-2, ymax=y+2, fill=col) +
        geom_rect(col="black", size=0.3) + 
        scale_fill_identity(guide = F) + 
        facet_grid(.~chrom, scales="free", space="free") + 
        theme_classic() + 
        scale_x_continuous(label=comma) +
        coord_cartesian(ylim=c(-4,4), expand=0) + 
        geom_rect(data = data.table(chrom=as.character(seqnames(wh)), start=start(wh), end=end(wh)), 
                  inherit.aes = F, alpha=0.3, fill="darkorange", col="darkorange",
                  aes(xmin=start, xmax=end, ymin=-4, ymax=4)) +
        xlab(NULL) + 
        theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
    return(plt)
}
