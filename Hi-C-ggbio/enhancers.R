# Bed files for SVs and enhancers:
assert_that(exists('FN.bed.sv.dup'))
assert_that(exists('FN.bed.sv.del'))
assert_that(exists('FN.bed.enh.macs'))
assert_that(exists('FN.bed.enh.cad4'))
assert_that(exists('FN.bedG.DHS'))


# Read all BED files (global SV.dt)
message("[global] loading SVs and enhancers")
SV.dt <- rbind(
    fread("/g/furlong/project/37_Capture-C/analysis/balancer_cap2/viewpoints.bed", select=1:3)[, V4:= "Capture-C viewpoints"],
    fread(FN.bed.sv.del)[, V4:= ifelse(substr(V4,1,7)=="0/1_1/1","Deletions (balancer)","Deletions (wild\u00adtype)")],
    fread(FN.bed.sv.dup, sep="\t")[, V4:= ifelse(substr(V4,1,7)=="0/1_1/1", "Duplications (balancer)", "Duplications (wild\u00adtype)")],
    fread(FN.bed.enh.cad4, select=1:3)[, V4:= "Enhancers (CAD4)"],
    fread(FN.bed.enh.macs, select=1:3)[, V4:= "Enhancers (DHS)"] )
colnames(SV.dt) = c("chrom", "start", "end", "name")
SV.dt <- SV.dt[order(chrom, start),]
SV.dt[, name := as.factor(name)]


# DHS data
message("[global] loading DHS track")
DHS.dt <- fread(paste("zcat", FN.bedG.DHS), stringsAsFactors=T, col.names=c("chrom","start","end","score"))
DHS.dt <- DHS.dt[grepl('chr[23][LR]', chrom),]
DHS.cutoff <- 3*quantile(DHS.dt$score, 0.9)
DHS.dt[score > DHS.cutoff,]$score = DHS.cutoff




# Read bed as GRanges to calculate overlaps (global: SV.ovlps)
annotate_overlaps_of_enhancer_and_DEL <- function()
{
    cad4 <- rtracklayer::import(FN.bed.enh.cad4)
    macs <- rtracklayer::import(FN.bed.enh.macs)
    mcols(cad4) <- data.frame(name = "Enhancers (CAD4)")
    mcols(macs) <- data.frame(name = "Enhancers (DHS)")
    enh <- c(macs, cad4)
    
    dt <- as.data.table(subsetByOverlaps(enh, rtracklayer::import(FN.bed.sv.del)))[, del:= "DEL"]
	dt[, chrom := seqnames]
    dt <- dt[order(seqnames, start), .(chrom, start, end, name, del)]
    dt
}
SV.ovlps <- annotate_overlaps_of_enhancer_and_DEL()

plt_SV_enhancers <- function(wh) 
{
    # Manual color assignment:
    bed.scale_cols <- c("Capture-C viewpoints" = "#e41a1c",
                        "Deletions (wild\u00adtype)" = "deepskyblue2",
                        "Duplications (wild\u00adtype)" = "firebrick1",
                        "Deletions (balancer)" = "firebrick2",
                        "Duplications (balancer)" = "deepskyblue3",
                        "Enhancers (CAD4)" = "#6f6f6f",
                        "Enhancers (DHS)" = "black")
    SV.dt.wh <- SV.dt[chrom == as.character(seqnames(wh))[1] & start >= start(wh)[1] & end <= end(wh)[1], ]
    SV.ovlps.wh <- SV.ovlps[chrom == as.character(seqnames(wh))[1] & start >= start(wh)[1] & end <= end(wh)[1], ]
    pSV <- ggplot(SV.dt.wh) + 
        aes(x = start, xend = end, y = name, yend = name, col = name) +
        geom_segment(size=3.5) +
        guides(color = "none") +
        ylab(NULL) + 
        scale_y_discrete(position = "right") +
        scale_color_manual(values = bed.scale_cols) + 
        geom_point(data = SV.ovlps.wh,
                   aes( x = (end+start)/2, y = name, col = del),
                   inherit.aes = F, size=2, col = "darkorange", alpha = 0.8) +
        coord_cartesian(expand = 0)
    return(pSV)
}


plt_DHS_track <- function(wh) 
{
    DHS.wh = DHS.dt[chrom == as.character(seqnames(wh)) & end >= start(wh)[1] & start <= end(wh)[1], ]

    plt <- ggplot(DHS.wh) +
        aes(xmin=start, xmax=end+1, ymax=score) + 
        geom_rect(ymin=0, color = "#377eb8") +
        ylab(NULL) +
        scale_y_continuous(limits = c(0, DHS.cutoff), breaks = scales::pretty_breaks(n = 3)) +
        coord_cartesian(expand = 0)
    return(plt)
}
