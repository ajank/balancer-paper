assert_that(exists("FN.bed.RNA.embryos.all.fwd"))
assert_that(exists("FN.bed.RNA.embryos.all.rev"))
assert_that(exists("FN.bed.RNA.embryos.bal.fwd"))
assert_that(exists("FN.bed.RNA.embryos.bal.rev"))
assert_that(exists("FN.bed.RNA.embryos.vrg.fwd"))
assert_that(exists("FN.bed.RNA.embryos.vrg.rev"))
assert_that(exists("FN.bed.RNA.adults.all.fwd"))
assert_that(exists("FN.bed.RNA.adults.all.rev"))
assert_that(exists("FN.bed.RNA.adults.bal.fwd"))
assert_that(exists("FN.bed.RNA.adults.bal.rev"))
assert_that(exists("FN.bed.RNA.adults.vrg.fwd"))
assert_that(exists("FN.bed.RNA.adults.vrg.rev"))

### RNA track
message("[global] Loading RNA tracks (embryo)...")
rna.fwd     <- rtracklayer::import(FN.bed.RNA.embryos.all.fwd)
rna.rev     <- rtracklayer::import(FN.bed.RNA.embryos.all.rev)
rna.bal.fwd <- rtracklayer::import(FN.bed.RNA.embryos.bal.fwd)
rna.bal.rev <- rtracklayer::import(FN.bed.RNA.embryos.bal.rev)
rna.vrg.fwd <- rtracklayer::import(FN.bed.RNA.embryos.vrg.fwd)
rna.vrg.rev <- rtracklayer::import(FN.bed.RNA.embryos.vrg.rev)

### RNA track
message("[global] Loading RNA tracks (adult)...")
rna2.fwd     <- rtracklayer::import(FN.bed.RNA.adults.all.fwd)
rna2.rev     <- rtracklayer::import(FN.bed.RNA.adults.all.rev)
rna2.bal.fwd <- rtracklayer::import(FN.bed.RNA.adults.bal.fwd)
rna2.bal.rev <- rtracklayer::import(FN.bed.RNA.adults.bal.rev)
rna2.vrg.fwd <- rtracklayer::import(FN.bed.RNA.adults.vrg.fwd)
rna2.vrg.rev <- rtracklayer::import(FN.bed.RNA.adults.vrg.rev)

#' Plot bedgraph file from RNAseq
plot_RNA_track <- function(wh, ds = "embryo", type = "all", ylim = 100)
{
    assert_that(is.character(ds) && length(ds)==1)
    assert_that(ds == "embryo" || ds == "adult")
    assert_that(is.character(type) && length(type)==1)
    assert_that(type == "all" || type == "bal" || type == "vrg")
    
    if (ds == "embryo") {
        if (type == "all") {
            fwd <- rna.fwd
            rev <- rna.rev
        } else if (type == "bal") {
            fwd <- rna.bal.fwd
            rev <- rna.bal.rev
        } else {
            fwd <- rna.vrg.fwd
            rev <- rna.vrg.rev
        }
    } else {
        if (type == "all") {
            fwd <- rna2.fwd
            rev <- rna2.rev
        } else if (type == "bal") {
            fwd <- rna2.bal.fwd
            rev <- rna2.bal.rev
        } else {
            fwd <- rna2.vrg.fwd
            rev <- rna2.vrg.rev
        }
    }
    
    f1 <- subsetByOverlaps(fwd, wh)
    f2 <- data.table(x = start(f1), x2=end(f1)+1, y=score(f1) )
    f2$y[f2$y > ylim] = ylim
    r1 <- subsetByOverlaps(rev, wh)
    r2 <- data.table(x = start(r1), x2=end(r1), y=score(r1) )
    r2$y[r2$y > ylim] = ylim
    
    plt <- ggplot(as.data.frame(f2)) + 
        aes(xmin=x, ymax = y, xmax = x2+1, ymin = 0) +
        geom_rect(color = "#6f6f6f") + 
        ylab(NULL) + 
        geom_rect(data = r2, color = "#6f6f6f", aes(ymin = -y, ymax = 0)) +
        geom_hline(yintercept=0, size=0.3) + 
        scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
        coord_cartesian(expand = 0)
    
    return(plt)
}
